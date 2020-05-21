PROGRAM main

  ! Stock-Watson UCSV model; estimated with Particle Learning

  USE embox, only : hrulefill, loadmatrix, loadarray1, savemat, savematlogical, savevec, storeEstimates, storeEstimatesTranspose, loft, timestampstr, es30d16, int2str
  USE blaspack, only : vech, ivech
  USE gibbsbox, only : drawNDXpdf, drawNDXsysresample

  USE vslbox
  USE timerbox
  USE omp_lib

  IMPLICIT NONE

  ! ----------------------------------------------------------------------------------

  logical, parameter :: doTimestamp = .false., doGains = .true.

  INTEGER, PARAMETER :: dof0 = 3

  INTEGER, PARAMETER :: Ny = 1, Nx = 2 , NsigmaX = Nx * (Nx + 1) / 2, Nw = 2, Nsv = 2
  INTEGER, PARAMETER :: ndxtrendRE = 1 , ndxgapRE = 2

  INTEGER :: Nparticles, Nmixturedraws
  INTEGER :: T,i,j,k,status 

  ! filter particles
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)     :: PARTICLEweights, DRAWllf
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: DRAWsvol, DRAWxhat, DRAWxsig, DRAWsqrtSigmaX
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DRAWxgain

  DOUBLE PRECISION, DIMENSION(Nx,Nx) :: sqrtVx0
  DOUBLE PRECISION, DIMENSION(Nx)    :: Ex0
  DOUBLE PRECISION, DIMENSION(Nsv)   :: SVar0, Eh0, Vh0

  ! DRAWS of various scale parameters
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWhInno

  ! priors for scale parameters
  DOUBLE PRECISION :: hvarT(Nsv)
  INTEGER :: hvarDof(Nsv)

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: y
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yNaN
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: loglike

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Xdraws

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: theta1
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: theta2
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ndx

  TYPE(progresstimer) :: timer
  CHARACTER (LEN=200) :: filename, datafile, nandatafile, fileXT, datalabel, modellabel, this

  ! VSL Random Stuff"
  type (vsl_stream_state) :: VSLstream
  integer :: seed
  ! integer :: brng
  integer :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  ! OPEN MP
  INTEGER :: NTHREADS !, TID

  ! ----------------------------------------------------------------------------------
  ! MODEL PARAMETERS
  ! ----------------------------------------------------------------------------------
  ! runtime parameters :start:
  ! first: set default values

  ! quick defaults
  Nparticles    = 10 ** 2
  Nmixturedraws = 10 ** 3

  datalabel        = 'cambridge2018GDPD'
  modellabel       = 'UCSV'



  ! ----------------------------------------------------------------------------------
  ! INIT
  ! ----------------------------------------------------------------------------------

  call hrulefill
  print *, 'Particle Filter estimation of ' // modellabel



  ! INIT OMP
  NTHREADS = 1
  !$OMP PARALLEL SHARED(NTHREADS)
  !$ NTHREADS = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  print *, "Number of Threads:", NTHREADS

  ! VSL
  ! brng    = vsl_brng_mt19937
  seed    = 0
  errcode = vslnewstream(VSLstream, vsl_brng_mt2203, seed)  

  WRITE(*,'(a25, i20, i20)') 'LAUNCHING VSLSTREAM ', VSLstream%descriptor1, VSLstream%descriptor2
  print *, 'vsl_brng', vsl_brng_mt2203

  ! runtime parameters :end: 

  call getarguments(datalabel, Nparticles)

  ! CONSTRUCT FILE EXTENTSIONS
  fileXT = '.particles.' // trim(datalabel) // '.' // trim(modellabel) // '.dat'
  if (doTimeStamp) filext = '.' // timestampstr() //  filext

  datafile    = trim(datalabel) // '.yData.txt'
  nandatafile = trim(datalabel) // '.yNaN.txt'

  ! read data
  T = loft(datafile) 
  IF (T < 10) THEN
     print *, 'Less than 10 observations in input file!', datafile
     STOP 1
  END IF

  ALLOCATE (y(Ny,T), yNaN(Ny,T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (Y)'
  END IF

  ! print *, 'trying to read', T, 'obs from', datafile
  CALL readdata(y,datafile,Ny,T)
  CALL readnandata(yNaN,nandatafile,Ny,T)

  ! validate yNaN and y
  DO k=1,T
     DO i = 1, Ny
        if (yNaN(i,k) .AND. y(i,k) /= 0.0d0 ) then
           write (*,*) 'YNAN PATTERN DOES NOT MATCH ZEROS IN Y'
        end if
     END DO
  END DO

  ! REPORT PARAMETERS TO SCREEN
  CALL HRULEFILL
  print *, 'data= ' // datalabel
  print *, 'model= ' // modellabel
  print *, 'Ny= ', Ny
  print *, 'T= ', T
  print *, 'Nparticles= ', Nparticles
  CALL HRULEFILL

  ! trivial since Nsv = 1 here
  Svar0    = (/ 0.6 / 3.0d0, 0.6 * 2.0d0 / 3.0d0 /)
  Vh0      = 10.0d0
  Eh0      = log(Svar0) - Vh0 * 0.5d0

  hvarDof = dof0
  hvarT   = (0.2d0 ** 2) * (dble(hvarDof) - 2.0d0)

  ! Linear prior
  Ex0       = 0.0d0
  Ex0(1)    = 2.0d0


  ! sqrtVx0, expressed as lower  triangular-choleski factor (sqrtVx0 * sqrtVx0')
  sqrtVx0      = 0.0d0
  sqrtVx0(1,1) = 100.0d0 ! sqrt(2.0d0)


  ! allocate memory for draws
  ALLOCATE (PARTICLEweights(Nparticles,0:T),DRAWllf(Nparticles,T), DRAWxhat(Nx,Nparticles,0:T), DRAWxsig(Nx,Nparticles,0:T), DRAWsvol(Nparticles,Nsv,0:T), DRAWxgain(Nx,Ny,Nparticles,T), DRAWsqrtSigmaX(NsigmaX,Nparticles,0:T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (draws)'
  END IF

  ! scale parameters
  ALLOCATE (DRAWhInno(Nparticles,Nsv,T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (scale parameter draws)'
  END IF

  PARTICLEweights = 1.0d0 / dble(Nparticles)
  DRAWllf         = 0.0d0
  DRAWxhat        = 0.0d0
  DRAWxsig        = 0.0d0
  DRAWsvol        = 0.0d0
  DRAWxgain       = 0.0d0

  DRAWhInno     = 0.0d0

  CALL particlefilter(T, Ny, y, yNaN, Nparticles, PARTICLEweights, DRAWllf, DRAWxhat, DRAWsqrtSigmaX, DRAWxsig, DRAWxgain, Nx, NsigmaX, Nw, Ex0, sqrtVx0, DRAWsvol, Nsv, Eh0, Vh0, DRAWhInno, hvarT, hvarDof, VSLstream,timer)


  CALL HRULEFILL
  WRITE (*,*) 'PARTICLE FILTER IS DONE!'
  CALL HRULEFILL

  ! WRITE SETTINGS
  CALL HRULEFILL
  filename = 'settings' // trim(adjustl(filext))
  OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
  WRITE(4,'(a20,a20)') 'TIME: ', timestampstr()
  WRITE(4,'(a20,a20)') 'Data: ', datalabel
  WRITE(4,'(a20,a20)') 'Model: ', modellabel
  WRITE(4,'(a40)') repeat('-',40)
  WRITE(4,'(a20,I20)') 'Nparticles: ', Nparticles
  CLOSE(UNIT=4)
  CALL HRULEFILL

  ! ----------------------------------------------------------------------------
  ! STORE
  ! ----------------------------------------------------------------------------

  CALL HRULEFILL
  WRITE (*,*) 'STARTING W/STORAGE !!!'
  CALL HRULEFILL

  ! STORE ESTIMATES
  ! Note: manual reshape avoids segmentation faults

  filename = 'YDATA' // filext
  call savemat(y, filename)

  filename = 'YNAN' // filext
  call savematlogical(yNaN, filename)


  ! LIKELIHOOD
  filename = 'LOGLIKE' // filext
  ALLOCATE (loglike(T), STAT=status)
  loglike  = log(sum(exp(DRAWllf),1) / Nparticles)
  call savevec(loglike, filename)
  call hrulefill
  WRITE (*,*) 'STORED LFF'
  WRITE (*,*) '... the loglikelihood is ', sum(loglike)
  WRITE (*,*) '... w/average contribution ', sum(loglike) / T
  call hrulefill
  DEALLOCATE (DRAWllf)

  ALLOCATE (theta1(T), STAT=status)

  ! ESS
  filename  = 'ESS' // filext
  theta1    = 1 / sum(PARTICLEweights(:,1:T) ** 2, 1) / Nparticles 
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ESS'


  ! store trends and gaps
  THIS      = 'TAUHATRE' 
  filename  = trim(this) // filext
  theta1 = sum(PARTICLEweights(:,1:T) * DRAWxhat(1,:,1:T), 1)
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ' // trim(this)

  THIS      = 'GAPHATRE' 
  filename  = trim(this) // filext
  theta1 = sum(PARTICLEweights(:,1:T) * DRAWxhat(2,:,1:T), 1)
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ' // trim(this)

  DEALLOCATE (theta1)

  ! draw distribution for linear states
  ALLOCATE (ndx(Nmixturedraws,T),Xdraws(Nx,Nmixturedraws,T))

  ! sample particle indices (note: these will also be used later for the SV particles)
  print *, 'Drawing particle indices  ...'
  DO j=1,T
     call drawNDXpdf(ndx(:,j), Nmixturedraws, PARTICLEweights(:,j), Nparticles, VSLstream)
  END DO
  print *, 'Done drawing particle indices.'

  print *, 'Drawing Xdraws normals ...'
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Nx * Nmixturedraws, Xdraws, 0.0d0, 1.0d0)
  print *, 'done drawing Xdraws normals.'

  FORALL (i=1:Nx,j=1:T,k=1:Nmixturedraws) Xdraws(i,k,j) = DRAWxhat(i,ndx(k,j),j) + DRAWxsig(i,ndx(k,j),j) * Xdraws(i,k,j) 



  ! store trends and gaps
  THIS      = 'TAURE' 
  filename  = trim(this) // filext
  CALL storeEstimatesTranspose(Xdraws(ndxtrendRE,:,:),T,Nmixturedraws,filename)
  WRITE (*,*) 'STORED ' // trim(this)

  THIS      = 'GAPRE' 
  filename  = trim(this) // filext
  CALL storeEstimatesTranspose(Xdraws(ndxgapRE,:,:),T,Nmixturedraws,filename)
  WRITE (*,*) 'STORED ' // trim(this)

  DEALLOCATE (Xdraws, DRAWxsig)



  ! 2D Gain Matrices
  if (doGains) then

     ! store analytical moments of gain
     ALLOCATE (theta2(T,Ny))
     ! trend gains 
     DO i=1,Nx
        filename  = 'GAIN' // trim(int2str(i)) // filext
        FORALL (j=1:Ny) theta2(:,j) = sum(PARTICLEweights(:,1:T) * DRAWxgain(i,j,:,:), 1)
        call savemat(theta2, filename)
        WRITE (*,*) 'STORED GAIN', i
     END DO
     DEALLOCATE (theta2)
  end if ! doGains
  DEALLOCATE(DRAWxgain)


  ! 2) Nparticle draws for the other particles
  ALLOCATE (theta1(T), STAT=status)
  DO i=1,Nsv
     filename  = 'SVHAT' // trim(int2str(i)) // filext
     theta1 = sum(PARTICLEweights(:,1:T) * DRAWsvol(:,i,1:T), 1)
     call savevec(theta1, filename)
     WRITE (*,*) 'STORED SVHAT', i
  END DO


  ! draw distribution
  ALLOCATE (theta2(Nmixturedraws,T))
  DO i=1,Nsv
     filename  = 'SV' // trim(int2str(i)) // filext
     FORALL (j=1:T,k=1:Nmixturedraws) theta2(k,j) = DRAWsvol(ndx(k,j),i,j) 
     CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED SV', i
  END DO


  ! SCALE PARAMETERS
  ! hInno
  DO i=1,Nsv
     filename  = 'HINNOHAT'  // trim(int2str(i))  // filext
     theta1 = sum(PARTICLEweights(:,1:T) * DRAWhInno(:,i,1:T), 1)
     call savevec(theta1, filename)
     WRITE (*,*) 'STORED HINNOHAT'  // trim(int2str(i))  
     ! draw distribution
     filename  = 'HINNO'  // trim(int2str(i))   // filext
     FORALL (j=1:T,k=1:Nmixturedraws) theta2(k,j) = DRAWhInno(ndx(k,j),i,j) 
     CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED HINNO'  // trim(int2str(i))  
  END DO

  DEALLOCATE(theta1)
  DEALLOCATE(theta2)
  DEALLOCATE(ndx)

  ! ----------------------------------------------------------------------------
  ! FINISHED: STORE FILTER
  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------
  ! CLEANUP FILTER
  ! ----------------------------------------------------------------------------

  call hrulefill
  WRITE (*,*) 'LFF:'
  WRITE (*,*) '... the loglikelihood is ', sum(loglike)
  WRITE (*,*) '... w/average contribution ', sum(loglike) / T
  call hrulefill
  DEALLOCATE (loglike)
  
  DEALLOCATE (PARTICLEweights)
  DEALLOCATE (DRAWsvol)
  DEALLOCATE (DRAWhInno)     

  DEALLOCATE (DRAWxhat, DRAWsqrtSigmaX)
  DEALLOCATE (y, yNaN)
  ! VSLstreams
  errcode = vsldeletestream(VSLstream)     


  call hrulefill
  WRITE(*,*) 'DONE. BYE, BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill

  STOP

CONTAINS

  ! -----------------------------------------------------------------
  SUBROUTINE getarguments(datalabel,Nparticles)

    INTENT(INOUT) Nparticles,datalabel

    INTEGER :: counter
    INTEGER :: Nparticles
    CHARACTER (LEN=100) :: datalabel
    CHARACTER(len=32) :: arg

    counter = 0

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, datalabel) 
    END IF

    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') Nparticles
    END IF


  END SUBROUTINE getarguments
  ! -----------------------------------------------------------------


  ! -----------------------------------------------------------------
  SUBROUTINE readdata(y,filename,Ny,T)
    IMPLICIT NONE

    INTENT(IN) :: filename,Ny,T
    INTENT(INOUT) :: y
    CHARACTER (LEN=200) :: filename
    CHARACTER (LEN=500) :: fmtstr

    DOUBLE PRECISION, DIMENSION(:,:) :: y
    INTEGER i, T, Ny

    fmtstr = es30d16(Ny)
    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')
    DO i=1,T
       READ(4,fmtstr) y(:,i)
    END DO

    CLOSE(UNIT=4)

  END SUBROUTINE readdata

  ! -----------------------------------------------------------------

  ! -----------------------------------------------------------------
  SUBROUTINE readnandata(nanny,filename,Ny,T)
    IMPLICIT NONE

    INTENT(IN) :: filename,T,Ny
    INTENT(INOUT) :: nanny
    CHARACTER (LEN=100) :: filename
    CHARACTER (LEN=500) :: fmtstr

    LOGICAL, DIMENSION(:,:) :: nanny
    INTEGER :: work(Ny)

    INTEGER i, j, T, Ny

    fmtstr = '(I2' // repeat(',I2', Ny-1) // ')'

    !Open File for reading
    OPEN (UNIT=4, FILE=filename, STATUS='OLD', ACTION='READ')

    DO i=1,T
       READ(4,fmtstr) (work(j), j=1,Ny)
       WHERE (work == 1) 
          nanny(:,i) = .TRUE.
       ELSEWHERE
          nanny(:,i) = .FALSE.
       END WHERE
    END DO

    CLOSE(UNIT=4)

  END SUBROUTINE readnandata
  ! -----------------------------------------------------------------


END PROGRAM main
! -----------------------------------------------------------------



! @\newpage\subsection{particlefilter}@
SUBROUTINE particlefilter(T, Ny, y, yNaN, Nparticles, PARTICLEweights, DRAWllf, DRAWxhat, DRAWsqrtSigmaX, DRAWxsig, DRAWxgain, Nx, NsigmaX, Nw, Ex0, sqrtVx00, DRAWsvol, Nsv, Eh0, Vh0, DRAWhInno, hvarT, hvarDof, VSLstream, timer)

  USE embox, only : savemat, savevec, int2str, mean, hrulefill

  use gibbsbox, only : drawNDXsysresample, igammaDraws
  use statespacebox, only : DLYAPsqrt
  use blaspack, only : pi, vechU,ivechU, vech, ivech, eye, symmetric, qrot, qrquery

  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE

  INTENT(INOUT) :: PARTICLEweights, DRAWllf, DRAWxhat, DRAWsqrtSigmaX, DRAWxsig, DRAWxgain, DRAWsvol, DRAWhInno, VSLstream, timer
  INTENT(IN)    :: T,Ny,y,yNaN, Nx, NsigmaX, Nw, Ex0, sqrtVx00, Nsv,Eh0,Vh0, Nparticles, hvarT, hvarDof

  INTEGER :: J, I, K, N, T, Nparticles, Nx, Ny, Nsv, NsigmaX, Nw

  ! OPEN MP
  INTEGER :: TID


  type(progresstimer) :: timer

  ! DOUBLE PRECISION, DIMENSION(Ny-1) :: horizons

  double precision, parameter :: minParticleWeight = 1.0d-12

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(Nparticles,T)   :: DRAWllf
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: DRAWxsig
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: DRAWxhat
  DOUBLE PRECISION, DIMENSION(NsigmaX,Nparticles, 0:T)  :: DRAWsqrtSigmaX
  DOUBLE PRECISION, DIMENSION(Nx,Ny,Nparticles,T) :: DRAWxgain
  DOUBLE PRECISION, DIMENSION(Nparticles,Nsv,0:T) :: DRAWsvol

  ! scale parameters
  DOUBLE PRECISION, DIMENSION(Nparticles,Nsv,T) :: DRAWhInno
  DOUBLE PRECISION, DIMENSION(Nsv,Nparticles)   :: hInno ! helper variable, note the transpose

  ! particles
  DOUBLE PRECISION :: xposterior(Nx,Nparticles), xsig(Nx,Nparticles), h(Nsv, Nparticles), SVol(Nsv,Nparticles), llf(Nparticles)
  INTEGER :: ndx(Nparticles)
  DOUBLE PRECISION :: shufflevec(Nparticles)

  ! state space objects
  DOUBLE PRECISION :: xprior(Nx), SigmaX(Nx,Nx), logdetSigmaY

  ! scale parameters
  DOUBLE PRECISION :: hvarT(Nsv)
  INTEGER :: hvarDof(Nsv)

  DOUBLE PRECISION :: PREVhvarT(Nparticles,Nsv)
  INTEGER :: PREVhvarDof(Nsv)
  DOUBLE PRECISION :: hDELTA(Nsv,Nparticles)

  ! SQRT objects
  DOUBLE PRECISION :: sqrtSigmaX(Nx,Nx), sqrtSigmaY(Ny,Ny), Kgain(Nx,Ny,Nparticles), qrR(Ny+Nx+Nw,Ny+Nx)
  DOUBLE PRECISION, DIMENSION(NsigmaX, Nparticles) :: vecSqrtSigmaX
  INTEGER :: qrLwork

  DOUBLE PRECISION :: Ex0(Nx), sqrtVx0(Nx,Nx), sqrtVx00(Nx,Nx), A(Nx,Nx), B(Nx,Nw), Bsv(Nx,Nw), C(Ny,Nx,T), sqrtR(Ny,Ny), ytilde(Ny)

  DOUBLE PRECISION, DIMENSION(Nsv) :: Eh0, Vh0
  DOUBLE PRECISION :: minSVh(Nsv)


  INTEGER :: Nynonan

  ! CHARACTER (LEN=200) :: filename

  ! VSL
  INTEGER :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
  type (vsl_stream_state) :: VSLstream
  double precision :: uniformdraws(2,T) ! two draws: one for each APF step

  ! index variables for state space
  INTEGER :: ndxTrendRE, ndxGapRE,  ndxGapREstart, ndxGapREstop, shockndxTrend, shockndxGap
  ! INTEGER :: ndxgap(1)  

  ! CHARACTER (LEN=100) :: filename

  minSVh   = log(0.001d0 ** 2)

  ! state-vector indices
  ndxTrendRE = 1
  ndxGapRE   = 2

  ndxGapREstart = ndxGapRE
  ndxGapREstop  = ndxGapRE

  ! ndxgap = (/ ndxgapREstart : ndxGapREstop /)

  shockndxTrend = 1
  shockndxGap   = 2

  ! init sufficient statistics of scale parameters
  forall(k=1:Nparticles,j=1:Nsv) PREVhvarT(k,j) = hvarT(j)
  forall(j=1:Nsv)                PREVhvarDof(j) = hvarDof(j)

  ! prepare state space
  ! A
  A = 0.0d0
  ! unit root in trend
  A(ndxTrendRE,ndxTrendRE) = 1.0d0

  ! B
  B                           = 0.0d0
  B(ndxTrendRE,shockndxTrend) = 1.0d0
  B(ndxGapRE,shockndxGap)     = 1.0d0

  ! C
  C         = 0.0d0
  ! inflation
  C(1,ndxTrendRE,:)    = 1.0d0
  C(1,ndxGapRE,:)      = 1.0d0

  ! Time 0 particles

  ! SV0
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, h, 0.0d0, 1.0d0)
  FORALL (j=1:Nsv,k=1:Nparticles) h(j,k) = Eh0(j) + sqrt(Vh0(j)) * h(j,k)  
  SVol = exp(h * 0.5d0)

  ! RB priors for linear states
  FORALL(k=1:Nparticles) xposterior(:,k) = Ex0

  ! prepare prior variance of linear states
  vecSqrtSigmaX   = 0.0d0
  sqrtVx0         = transpose(sqrtVx00)
  ! DO k=1,Nparticles
  !    vecSqrtSigmaX(:,k)     = vechU(sqrtVx0,Nx)
  ! END DO
  FORALL (k=1:Nparticles)  vecSqrtSigmaX(:,k)     = vechU(sqrtVx0,Nx)


  FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,0)  = vecSqrtSigmaX(i,k)
  FORALL(i=1:Nx,k=1:Nparticles)  DRAWxhat(i,k,0)            = xposterior(i,k) 

  FORALL (i=1:Nsv,k=1:Nparticles) DRAWsvol(k,i,0) = SVol(i,k) ! transpose
  PARTICLEweights = 1 / dble(Nparticles)
  DRAWllf         = 0.0d0

  ! uniform draws for systematic resampling
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


  ! workspace query for qr decomposition
  qrR = 0.0d0
  qrlwork = qrquery(qrR)

  CALL initprogressbar(timer, 15.0d0)
  DO j=1,T

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: APF RESAMPLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     Nynonan = count(.not. yNaN(:,j))

     !$OMP PARALLEL DO SHARED(xposterior, xsig, vecSqrtSigmaX, SVol, llf, Nparticles, Nsv, j, y, yNaN, Ny, Nx, Nw, Nynonan, ndxTrendRE, ndxGapRE, ndxGapREstart, ndxGapREstop, shockndxTrend, shockndxGap) FIRSTPRIVATE(A,B,C,qrLWORK) PRIVATE(ytilde, sqrtR, Bsv, logdetSigmaY, errcode, xprior, SigmaX, sqrtSigmaX, sqrtSigmaY, qrR, TID) DEFAULT(NONE) SCHEDULE(STATIC)


     DO k = 1,Nparticles

        TID = 0
        !$ TID = OMP_GET_THREAD_NUM()


        ! 2) Fill Particles into state space
        sqrtR = 0.0d0

        ! zero out missing obs
        DO i = 1, Ny
           if (yNaN(i,j)) then 
              C(i,:,j)     = 0.0d0 
              ! sqrtR(i,:)   = 0.0d0
           end if
        END DO


        ! Bsv
        FORALL (i=1:Nsv) Bsv(:,i)  = B(:,i) * SVol(i,k) 

        ! 3) Kalman Filter

        ! xprior = A * xposterior(-1)
     	xprior      = 0.0d0
        call DGEMV('n',Nx,Nx,1.0d0,A,Nx,xposterior(:,k),1,0.0d0,xprior,1)

        sqrtSigmaX = ivechU(vecSqrtSigmaX(:,k),Nx)

        ! fill directly into qrR
        qrR = 0.0d0
        qrR(1:Ny,1:Ny) = transpose(sqrtR)
        qrR(Ny+Nx+1:Ny+Nx+Nw,Ny+1:Ny+Nx) = transpose(Bsv)
        ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
        call DGEMM('n','t',Nx,Nx,Nx,1.0d0,sqrtSigmaX,Nx,A,Nx,0.0d0,qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx),Nx)
        ! qrR(Ny+1:Ny+Nx+Nw,1:Ny) = qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx) * C'
        call DGEMM('n','t',Nx+Nw,Ny,Nx,1.0d0,qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx),Nx+Nw,C(:,:,j),Ny,0.0d0,qrR(Ny+1:Ny+Nx+Nw,1:Ny),Nx+Nw)

        ! QR decomposition
        call qrot(qrR, qrLWORK)

        ! map qr into Kalman objects
        sqrtSigmaY   = qrR(1:Ny,1:Ny) ! upper triangular
        sqrtSigmaX   = qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) ! upper triangular
        ! Kgain(:,:,k) = transpose(qrR(1:Ny,Ny+1:Ny+Nx))

        ! ytilde = y - C * xprior
        ytilde = y(:,j)
        call DGEMV('n',Ny,Nx,-1.0d0,C(:,:,j),Ny,xprior,1,1.0d0,ytilde,1)

        ! ! compute resid
        ! do i=1,Ny
        !    if (.NOT. yNaN(i,j)) then
        !       resid(i,k) = ytilde(i) * sqrtR(i,i) / sqrtSigmaY(i,i)
        !    end if
        ! end do

        ! singularity fix: insert unit dummies for missing values
        do i=1,Ny
           if (yNaN(i,j)) sqrtSigmaY(i,i) = 1.0d0
           ! ytilde(i) = 0.0d0 -- not needed since y(i)=0 and C(i) * xprior also zero
        end do

        logdetSigmaY = 0.0d0
        DO i=1,Ny
           logdetSigmaY = logdetSigmaY + log(abs(sqrtSigmaY(i,i)))
        END DO
        logdetSigmaY = 2.0d0 * logdetSigmaY

        ! rotate ytilde (up to sign, consistent with rotation of K -- needed for llf computation)
        call dtrsv('U', 'T', 'N', Ny, sqrtSigmaY, Ny, ytilde, 1) ! recall: sqrtSigmaY is returned as upper triangular, right factor

        ! DO NOT STORE POSTERIORS FOR APF STEP
        ! ! xposterior = xprior + K * ytilde
        ! xposterior(:,k) = xprior
        ! call DGEMV('N',Nx,Ny,1.0d0,Kgain(:,:,k),Nx,ytilde,1,1.0d0,xposterior(:,k),1)


        ! rotate Kalman gain
        ! call dtrsm('R', 'U', 'T', 'N', Nx, Ny, 1.0d0, sqrtSigmaY, Ny, Kgain(:,:,k), Nx) ! recall: sqrtSigmaY is returned as upper triangular, right factor

        ! ! remove unit dummies
        ! do i=1,Ny
        !    if (yNaN(i,j)) sqrtSigmaY(i,i) = 0.0d0
        ! end do

        ! compute log-likelihood
        ! llf
        llf(k)       = -0.5d0 * (Nynonan * log(2.0d0 * pi) + logdetSigmaY + sum(ytilde ** 2))

        ! DO NOT STORE POSTERIORS FOR APF STEP
        ! vecSqrtSigmaX(:,k) = vechU(sqrtSigmaX,Nx)
        ! ! store state variances 
        ! call dsyrk('u', 't', Nx, Nx, 1.0d0, sqrtSigmaX, Nx, 0.0d0, SigmaX, Nx)
        ! forall (i=1:Nx) xsig(i,k) = sqrt(SigmaX(i,i))


     END DO ! k particles
     !$OMP END PARALLEL DO 


     if (Nynonan > 0) then

        ! Reweight particles for next round   
        PARTICLEweights(:,j) = exp(llf)
        PARTICLEweights(:,j) = PARTICLEweights(:,j) / sum(PARTICLEweights(:,j))

        where (PARTICLEweights(:,j) < minParticleWeight) PARTICLEweights(:,j) = minParticleWeight
        PARTICLEweights(:,j) = PARTICLEweights(:,j) / sum(PARTICLEweights(:,j))

        call drawNDXsysresample(ndx, Nparticles, PARTICLEweights(:,j), Nparticles, uniformdraws(1,j))

        FORALL(k=1:Nparticles) shufflevec(k) = PARTICLEweights(ndx(k),j)
        PARTICLEweights(:,j) = shufflevec

        DO i=1,Nx
           FORALL(k=1:Nparticles) shufflevec(k) = xposterior(i,ndx(k))
           xposterior(i,:) = shufflevec
        END DO

        DO i=1,Nsigmax
           FORALL(k=1:Nparticles) shufflevec(k) = vecSqrtSigmaX(i,ndx(k))
           vecSqrtSigmaX(i,:) = shufflevec
        END DO

        DO i=1,Nsv
           FORALL(k=1:Nparticles) shufflevec(k) = h(i,ndx(k))
           h(i,:) = shufflevec
        END DO

        DO i=1,Nsv
           FORALL(k=1:Nparticles) shufflevec(k) = PREVhvarT(ndx(k),i)
           PREVhvarT(:,i) = shufflevec
        END DO

     else ! i.e. Nynonan == 0
        ! do nothing since no information received
        ! PARTICLEweights(:,j) = 1 / dble(Nparticles)
     end if

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: APF STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! 0) draw scale parameters
     ! hInno
     DO i=1,Nsv 
        call igammaDraws(DRAWhinno(:,i,j), Nparticles, PREVhvarT(:,i), PREVhvarDof(i), VSLstream) 
     END DO
     DRAWhinno(:,:,j) = sqrt(DRAWhinno(:,:,j))
     forall (i=1:Nsv,k=1:Nparticles) hInno(i,k) = DRAWhInno(k,i,j) ! helper variable, provides better aligned access to the j data


     ! 1) Draw Particles 
     errcode     = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, hDELTA, 0.0d0, 1.0d0)
     hDELTA      = hDELTA * hInno
     h           = h + hDELTA
     SVol        = exp(h * 0.5d0)

     Nynonan = count(.not. yNaN(:,j))


     !$OMP PARALLEL DO SHARED(xposterior, xsig, vecSqrtSigmaX, SVol, Kgain, llf, Nparticles, Nsv, j, y, yNaN, Ny, Nx, Nw, Nynonan, ndxTrendRE, ndxGapRE, ndxGapREstart, ndxGapREstop, shockndxTrend, shockndxGap) FIRSTPRIVATE(A,B,C,qrLWORK) PRIVATE(ytilde, sqrtR, Bsv, logdetSigmaY, errcode, xprior, SigmaX, sqrtSigmaX, sqrtSigmaY, qrR, TID) DEFAULT(NONE) SCHEDULE(STATIC)


     DO k = 1,Nparticles

        TID = 0
        !$ TID = OMP_GET_THREAD_NUM()

        xprior      = 0.0d0
        ! SigmaX      = 0.0d0
        ! call ivech(SigmaX, vecSigmaX(:,k))

        ! 2) Fill Particles into state space
        sqrtR = 0.0d0
        ! forall (i=1:Ny) sqrtR(i,i) = noisevol(i)

        ! zero out missing obs
        DO i = 1, Ny
           if (yNaN(i,j)) then 
              C(i,:,j)     = 0.0d0 
              ! sqrtR(i,:)   = 0.0d0
           end if
        END DO

        ! Bsv
        FORALL (i=1:Nsv) Bsv(:,i)  = B(:,i) * SVol(i,k) 

        ! 3) Kalman Filter

        ! xprior = A * xposterior(-1)
        call DGEMV('n',Nx,Nx,1.0d0,A,Nx,xposterior(:,k),1,0.0d0,xprior,1)

        ! ------------------------------------------------------------------------
        ! SQRT KALMAN
        ! ------------------------------------------------------------------------
        sqrtSigmaX = ivechU(vecSqrtSigmaX(:,k),Nx)

        ! fill directly into qrR
        qrR = 0.0d0
        qrR(1:Ny,1:Ny) = transpose(sqrtR)
        qrR(Ny+Nx+1:Ny+Nx+Nw,Ny+1:Ny+Nx) = transpose(Bsv)
        ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
        call DGEMM('n','t',Nx,Nx,Nx,1.0d0,sqrtSigmaX,Nx,A,Nx,0.0d0,qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx),Nx)
        ! qrR(Ny+1:Ny+Nx+Nw,1:Ny) = qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx) * C'
        call DGEMM('n','t',Nx+Nw,Ny,Nx,1.0d0,qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx),Nx+Nw,C(:,:,j),Ny,0.0d0,qrR(Ny+1:Ny+Nx+Nw,1:Ny),Nx+Nw)

        ! QR decomposition
        call qrot(qrR, qrLWORK)

        ! map qr into Kalman objects
        sqrtSigmaY   = qrR(1:Ny,1:Ny) ! upper triangular
        sqrtSigmaX   = qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) ! upper triangular
        Kgain(:,:,k) = transpose(qrR(1:Ny,Ny+1:Ny+Nx))

        ! ytilde and logdetSigmaY
        ytilde = y(:,j)
        call DGEMV('n',Ny,Nx,-1.0d0,C(:,:,j),Ny,xprior,1,1.0d0,ytilde,1)

        ! singularity fix: insert unit dummies for missing values
        do i=1,Ny
           if (yNaN(i,j)) sqrtSigmaY(i,i) = 1.0d0
           ! ytilde(i) = 0.0d0 -- not needed since y(i)=0 and C(i) * xprior also zero
        end do

        logdetSigmaY = 0.0d0
        DO i=1,Ny
           logdetSigmaY = logdetSigmaY + log(abs(sqrtSigmaY(i,i)))
        END DO
        logdetSigmaY = 2.0d0 * logdetSigmaY
        ! rotate/normalize ytilde (up to sign, consistent with rotation of K)
        call dtrsv('U', 'T', 'N', Ny, sqrtSigmaY, Ny, ytilde, 1) ! recall: sqrtSigmaY is returned as upper triangular, right factor

        ! xposterior = xprior + K * ytilde
        xposterior(:,k) = xprior
        call DGEMV('N',Nx,Ny,1.0d0,Kgain(:,:,k),Nx,ytilde,1,1.0d0,xposterior(:,k),1)


        ! rotate Kalman gain into space of non-normalized ytilde
        call dtrsm('R', 'U', 'T', 'N', Nx, Ny, 1.0d0, sqrtSigmaY, Ny, Kgain(:,:,k), Nx) ! recall: sqrtSigmaY is returned as upper triangular, right factor

        ! remove unit dummies -- can be omitted since sqrtSigmaY not used any further
        ! do i=1,Ny
        !    if (yNaN(i,j)) sqrtSigmaY(i,i) = 0.0d0
        ! end do

        ! compute log-likelihood
        ! llf
        llf(k)       = -0.5d0 * (Nynonan * log(2.0d0 * pi) + logdetSigmaY + sum(ytilde ** 2))

        vecSqrtSigmaX(:,k) = vechU(sqrtSigmaX,Nx)
        ! store state variances 
        call dsyrk('u', 't', Nx, Nx, 1.0d0, sqrtSigmaX, Nx, 0.0d0, SigmaX, Nx)
        forall (i=1:Nx) xsig(i,k) = sqrt(SigmaX(i,i))


        ! ------------------------------------------------------------------------
        ! DONE: SQRT KALMAN
        ! ------------------------------------------------------------------------

     END DO ! k particles
     !$OMP END PARALLEL DO 


     ! Store NON-reweighted statistics
     DRAWllf(:,j)      = llf
     FORALL(i=1:Nsv,k=1:Nparticles) DRAWsvol(k,i,j)   = SVol(i,k) ! note the transpose

     FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,j)  = vecSqrtSigmaX(i,k)
     FORALL(i=1:Nx,k=1:Nparticles) DRAWxhat(i,k,j)             = xposterior(i,k)
     FORALL(i=1:Nx,k=1:Nparticles) DRAWxsig(i,k,j)             = xsig(i,k) 
     FORALL(i=1:Nx,n=1:Ny,k=1:Nparticles) DRAWxgain(i,n,k,j)   = Kgain(i,n,k) ! Kprime(n,i,k)

     if (Nynonan > 0) then

        ! Reweight particles for next round   
        PARTICLEweights(:,j) = exp(llf) / PARTICLEweights(:,j) 
        PARTICLEweights(:,j) = PARTICLEweights(:,j) / sum(PARTICLEweights(:,j))

        where (PARTICLEweights(:,j) < minParticleWeight) PARTICLEweights(:,j) = minParticleWeight
        PARTICLEweights(:,j) = PARTICLEweights(:,j) / sum(PARTICLEweights(:,j))

        call drawNDXsysresample(ndx, Nparticles, PARTICLEweights(:,j), Nparticles, uniformdraws(2,j))


        DO i=1,Nx
           FORALL(k=1:Nparticles) shufflevec(k) = xposterior(i,ndx(k))
           xposterior(i,:) = shufflevec
        END DO

        DO i=1,Nsigmax
           FORALL(k=1:Nparticles) shufflevec(k) = vecSqrtSigmaX(i,ndx(k))
           vecSqrtSigmaX(i,:) = shufflevec
        END DO

        DO i=1,Nsv
           FORALL(k=1:Nparticles) shufflevec(k) = h(i,ndx(k))
           h(i,:) = shufflevec
        END DO
        DO i=1,Nsv ! also reshuflle SVol to prep the APF resample step
           FORALL(k=1:Nparticles) shufflevec(k) = SVol(i,ndx(k))
           SVol(i,:) = shufflevec
        END DO

        ! reshuffle sufficient statistics for scale parameters
        DO i=1,Nsv
           FORALL(k=1:Nparticles) shufflevec(k) = PREVhvarT(ndx(k),i)
           PREVhvarT(:,i) = shufflevec
        END DO

        ! propagate sufficient statistics (including reshuffling of DELTAs)
        FORALL(k=1:Nparticles,i=1:Nsv) PREVhvarT(k,i)       = PREVhvarT(k,i)    + hDELTA(i,ndx(k)) ** 2

        PREVhvarDof      = PREVhvarDof + 1
     else ! i.e. Nynonan == 0
        !   PARTICLEweights(:,j) = 1 / dble(Nparticles)
     end if


     CALL progressbarcomment(dble(j) / dble(T), timer, 'Particle Step')

  END DO ! j=1,T

END SUBROUTINE particlefilter

