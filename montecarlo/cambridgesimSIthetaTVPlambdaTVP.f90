PROGRAM main

  ! time-varying AR1

  USE embox, only    : mean, median, hrulefill, loadmatrix, loadarray1, savemat, savematlogical, savevec, storeEstimates, storeEstimatesTranspose, loft, timestampstr, es30d16, int2str
  USE blaspack, only : vech, ivech, eye
  USE gibbsbox, only : drawNDXpdf, drawNDXsysresample

  use cambridgebox, only : ucsvSI, generateUCSVSIdata

  USE vslbox
  USE timerbox
  USE omp_lib

  IMPLICIT NONE

  ! ----------------------------------------------------------------------------------

  CHARACTER (LEN=200) :: modellabel ='SIthetaTVPlambdaTVP'

  logical, parameter :: doTimestamp = .false., doGains = .true.
  logical :: doInflationNoise = .true.
  logical :: doStoreDGP = .false.

  logical :: doSecondResamplingStep = .false.
  integer :: dofAPFfattail = 0 ! set to zero to use normal APF
  logical, parameter :: doPercentError = .false.
  INTEGER :: Ngrid = 10, gridcounter

  ! smoothing parameters
  LOGICAL :: doSmoother
  INTEGER :: Nsmoother, smootherNparticles

  INTEGER, PARAMETER :: dof0 = 3

  INTEGER, PARAMETER :: p = 1, Nsurveys = 5, Ny = Nsurveys + 1, Nx = 2 * (1 + p), NsigmaX = Nx * (Nx + 1) / 2, Nw = 2, Nsv = 2
  INTEGER, PARAMETER :: ndxtrendRE = 1 , ndxgapRE = 2 , ndxtrendSI = 1 + p + 1 , ndxgapSI = ndxtrendSI + 1

  double precision :: noisevol0(Ny)  = sqrt(0.1d0) 
  double precision :: noisevol00(Ny) = sqrt(0.1d0) 

  INTEGER :: Nparticles, Nmixturedraws
  INTEGER :: T,i,j,k,status 

  ! filter particles
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)     :: PARTICLEweights, DRAWlambda
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: logMDD
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: DRAWsvol, DRAWxhat, DRAWxsig, DRAWsqrtSigmaX
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DRAWxgain

  ! smoothed parameters
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: PARAMweights
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)     :: PARAMhInno, PARAMsigma
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: PARAMsiglambda
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: PARAMsiga


  ! smoother particles
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)     :: SMOOTHERlambda
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: SMOOTHERsvol
  ! DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: SMOOTHERx
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: SMOOTHERxhat
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: SMOOTHERsigma 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: SMOOTHERhinno 

  DOUBLE PRECISION, DIMENSION(Nx,Nx) :: sqrtVx0,sqrtVx00
  DOUBLE PRECISION, DIMENSION(Nx)    :: Ex0
  DOUBLE PRECISION, DIMENSION(Nsv)   :: SVar0, Eh0, Vh0

  DOUBLE PRECISION :: lambda0, lambda0V
  DOUBLE PRECISION :: a0, a0V 

  ! arrays for initial values
  double precision, allocatable, dimension(:)   :: lambda00, a00 
  double precision, allocatable, dimension(:,:) :: x00

  ! DRAWS of various scale parameters
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: DRAWsiglambda, DRAWsiga
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWsigma, DRAWhInno

  ! priors for scale parameters
  DOUBLE PRECISION :: sigaT, siglambdaT, hvarT(Nsv), sigmaT(Ny)
  INTEGER :: sigaDof, siglambdaDof, hvarDof(Nsv), sigmaDof(Ny)
  DOUBLE PRECISION :: siga0, siglambda0,  hinno0(Nsv)
  DOUBLE PRECISION :: siga00, siglambda00

  ! DRIFT AR 1 (specialized to p=1)
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DRAWa
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SMOOTHERa

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: y
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yNaN

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: thetavec

  ! define a derived type for cambridgeSI true model parameters and states
  type(ucsvSI(:,:,:,:)), allocatable :: DGP

  ! define bias arrays
  double precision, dimension(:,:,:), allocatable :: biasX, biasSV
  double precision, dimension(:,:), allocatable   :: biasHinno
  double precision, dimension(:,:), allocatable   :: biasLambda, biasTheta
  double precision, dimension(:), allocatable     :: biasSigTheta, biasSigLambda
  double precision, dimension(:,:), allocatable   :: biasNoisevol

  double precision, dimension(:,:,:), allocatable :: smootherBiasX, smootherBiasSV
  double precision, dimension(:,:), allocatable   :: smootherBiasLambda, smootherBiasTheta

  ! other ...
  INTEGER :: ompcount

  TYPE(progresstimer) :: timer
  CHARACTER (LEN=200) :: filename, fileXT, datalabel, this

  ! VSL Random Stuff
  type (vsl_stream_state) :: VSLstream
  integer :: seed
  type (vsl_stream_state), allocatable, dimension(:) :: VSLarray

  ! integer :: brng
  integer :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  ! OPEN MP
  INTEGER :: NTHREADS, TID

  ! ----------------------------------------------------------------------------------
  ! MODEL PARAMETERS
  ! ----------------------------------------------------------------------------------
  ! runtime parameters :start:
  ! first: set default values

  ! thorough
  Nparticles    = 10 ** 5
  Nmixturedraws = 10 ** 4
  Nsmoother     = 10 ** 4


  ! quick
  Nparticles    = 10 ** 4
  Nmixturedraws = 10 ** 3
  Nsmoother     = 10 ** 4
  smootherNparticles = 10 ** 4

  T = 100
  ! ----------------------------------------------------------------------------------
  ! INIT
  ! ----------------------------------------------------------------------------------

  ! INIT OMP
  NTHREADS = 1

  !$OMP PARALLEL SHARED(NTHREADS)
  !$ NTHREADS = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  print *, "Number of Threads:", NTHREADS
  print *, 'KMP stack size', kmp_get_stacksize_s() 

  ! VSL
  allocate(VSLarray(0:NTHREADS-1))
  seed = 0
  errcode = vslnewstream(VSLstream, vsl_brng_mt2203, seed)
  ! VSL array
  DO i=1,NTHREADS
     errcode = vslnewstream(VSLarray(i-1), vsl_brng_mt2203 + i, seed)
     if (errcode /= 0) then
        print *, 'error init VSLarray', i
        stop 1
     end if

     WRITE(*,'(a25, i4, i20, i20)') 'LAUNCHING VSLARRAY ', i-1, VSLarray(i-1)%descriptor1, VSLarray(i-1)%descriptor2

  END DO

  ! runtime parameters :end: 

  call getarguments(Ngrid, T, doInflationNoise, Nparticles, doStoreDGP, doSmoother, Nsmoother, smootherNparticles, doSecondResamplingstep, dofAPFfattail)
  call hrulefill
  print *, 'Particle Filter estimation of ' // modellabel

  ! CONSTRUCT FILE EXTENTSIONS
  fileXT = trim(datalabel) // '.' // trim(modellabel)
  if (.not. doInflationNoise)  fileXT = trim(filext) // '.nonoise'
  if (doSecondResamplingStep)  fileXT = trim(filext) // '.2ndAPFresample'
  if (dofAPFfattail .gt. 0)    fileXT = trim(filext) // '.dofAPFfattail' // trim(int2str(dofAPFFattail))
  fileXT = trim(filext) // '.dat'
  if (doTimeStamp) filext = '.' // timestampstr() //  filext


  datalabel        = 'simdataT' // trim(int2str(T)) 

  fileXT =  '.' // trim(datalabel) // '.' // trim(modellabel)
  if (.not. doInflationNoise)  fileXT = trim(filext) // '.nonoise'
  if (doSmoother) then
     fileXT =  trim(filext) // '.Nsmoother'// trim(int2str(Nsmoother))  // '.smootherNparticles'// trim(int2str(smootherNparticles)) 
  end if
  fileXT =  trim(filext) // '.Nparticles'// trim(int2str(Nparticles)) // '.Ngrid'// trim(int2str(Ngrid)) // '.dat'


  ! REPORT PARAMETERS TO SCREEN
  CALL HRULEFILL
  print *, 'data=       ' // datalabel
  print *, 'seed=       ', seed
  print *, 'model=      ' // modellabel
  print *, 'Ngrid=      ', Ngrid
  print *, 'Ny=         ', Ny
  print *, 'T=          ', T
  print *, 'Nparticles= ', Nparticles
  if (doSmoother) then
     print *, 'Nsmoother=          ', Nsmoother
     print *, 'smootherNparticles= ', smootherNparticles
  end if
  print *, 'p=          ', p
  if (doSecondResamplingStep) then
     print *, 'APF uses second resampling step'
  end if
  if (dofAPFfattail .gt. 0) then
     print *, 'APF-MV-T w/ ', dofAPFfattail, ' dof'
  end if
  if (doInflationNoise) then
     print *, 'Model variant with Noise in Inflation'
  else
     print *, 'Model variant WITHOUT Noise in Inflation'
  end if
  if (doSmoother) then
     print *, 'With smoother'
  end if
  CALL HRULEFILL

  ! trivial since Nsv = 1 here
  Svar0    = (/ 0.6 / 3.0d0, 0.6 * 2.0d0 / 3.0d0 /)
  Vh0      = 10.0d0
  Eh0      = log(Svar0) - Vh0 * 0.5d0

  hvarDof = dof0
  hInno0  = 0.2d0
  hvarT   = (hInno0 ** 2) * (dble(hvarDof) - 2.0d0)

  ! lambda
  lambda0      = 0.5d0
  lambda0V     = 1.0d0
  siglambdaDof = dof0
  siglambda0   = 0.1d0
  siglambdaT   = (siglambda0 ** 2) * (dble(siglambdaDof) - 2.0d0)

  ! a
  a0        = 0.0d0
  a0V       = 1.0d0
  sigaDof   = dof0
  siga0     = 0.1d0
  sigaT     = (siga0 ** 2) * (dble(sigaDof) - 2.0d0)

  ! Linear prior
  Ex0       = 0.0d0
  Ex0(1)    = 2.0d0
  Ex0(3)    = 2.0d0
  
  ! sqrtVx0, expressed as lower  triangular-choleski factor (sqrtVx0 * sqrtVx0')
  ! ORIGINAL PRIOR
  call eye(sqrtVx0, 10.0d0) ! non-zero prior gap-variance matters only for initial conditions for sampling lagged gap as used for estimating AR(1) coefficient 
  sqrtVx0(1,1) = 100.0d0 ! sqrt(2.0d0)
  sqrtVx0(3,3) = 100.0d0 ! sqrt(2.0d0) 

  ! REVISED PRIOR
  ! call eye(sqrtVx0, 10.0d0) 
  ! sqrtVx0(1,1)   = 100.0d0 
  ! sqrtVx0(3,1)   = sqrtVx0(1,1)
  ! sqrtVx0(4,2)   = sqrtVx0(2,2)

  sigmaDof    = 20 ! dof0
  sigmaT      = (noisevol0 ** 2) * (dble(sigmaDof) - 2.0d0)

  ! draw *true* initial values
  allocate(lambda00(Ngrid), a00(Ngrid))
  ! errcode = vdrnguniform(VSLmethodUniform, VSLstream, Ngrid, lambda00, 0.0d0, 1.0d0 )
  lambda00     = lambda0 ! RESET
  ! errcode = vdrnguniform(VSLmethodUniform, VSLstream, Ngrid, a00, -1.0d0, 1.0d0 )
  a00          = a0

  ! x00
  sqrtVx00       = sqrtVx0
  sqrtVx00(1,1)  = 0.0d0 ! to avoid drawing huge numbers for the trend component ...
  sqrtVx00(3,1)  = sqrtVx00(1,1) ! to avoid drawing huge numbers for the trend component ...
  ! sqrtVx00(3,3)  = 20.0d0
  ! sqrtVx00(4,2)  = sqrtVx00(2,2)

  filename = 'sqrtVx00' // fileXT
  call savemat(sqrtVx00, filename)
  filename = 'sqrtVx0' // fileXT
  call savemat(sqrtVx0, filename)

  allocate(x00(Nx,Ngrid))
  forall (gridcounter=1:Ngrid) x00(:,gridcounter) = Ex0 
  allocate(thetavec(Nx,Ngrid))
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx * Ngrid, thetavec, 0.0d0, 1.0d0)
  DO gridcounter=1,Ngrid
     call DGEMV('N',Nx,Nx,1.0d0,sqrtVx00,Nx,thetavec(:,gridcounter),1,1.0d0,x00(:,gridcounter),1)
  END DO
  deallocate(thetavec)

  ! apply maximum to x00
  where (x00(ndxTrendSI,:) >  50.0d0 + 2.0d0) x00(ndxTrendSI,:) =  50.0d0 + 2.0d0
  where (x00(ndxTrendSI,:) < -50.0d0 + 2.0d0) x00(ndxTrendSI,:) = -50.0d0 + 2.0d0
  filename = 'x00' // fileXT
  call savemat(x00, filename)

  siga00       = siga0
  siglambda00  = siglambda0 

  ! release default stream
  errcode = vsldeletestream(VSLstream)

  allocate (biasX(Nx,0:T,Ngrid),biasSV(Nsv,0:T,Ngrid))
  allocate (biasTheta(0:T,Ngrid),biasLambda(0:T,Ngrid))
  allocate (biasHinno(Nsv,Ngrid))
  allocate (biasSigTheta(Ngrid))
  allocate (biasSigLambda(Ngrid))
  allocate (biasNoisevol(Ny,Ngrid))

  if (doSmoother) then
     allocate (smootherBiasX(Nx,0:T,Ngrid),smootherBiasSV(Nsv,0:T,Ngrid))
     allocate (smootherBiasTheta(0:T,Ngrid),smootherBiasLambda(0:T,Ngrid))
  end if

  CALL initprogressbar(timer, 15.0d0)
  TID = 0
  OMPcount = 0
  !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(biasX,biasSV,biasHinno,biasTheta,biasLambda,biasSigTheta,biasSigLambda,biasNoisevol)  SHARED(smootherBiasX,smootherBiasSV,smootherBiasTheta,smootherBiasLambda) SHARED(VSLarray) SHARED(timer,ompcount) SHARED(lambda00,a00,x00) SCHEDULE(STATIC) 
  
  do gridcounter = 1,Ngrid

     !$OMP ATOMIC UPDATE
     ompCount = ompCount + 1
     !$OMP END ATOMIC 

     !$ TID = OMP_GET_THREAD_NUM()

     ! GENERATE SIMULATED DATA
     allocate (ucsvSI(Ny,Nx,Nsv,T)::dgp) 
     ALLOCATE (y(Ny,T), yNaN(Ny,T), STAT=status)
     IF (status /= 0) THEN
        WRITE (*,*) 'Allocation problem (Y)'
     END IF

     call generateUCSVSIdata(dgp, doInflationNoise, T, Ny, x00(:,gridcounter), Eh0, hinno0, a00(gridcounter), siga00, lambda00(gridcounter), siglambda00, noisevol00, VSLarray(TID))

     yNaN = .false.
     y    = dgp%y

     if (doStoreDGP) then
        ! store truth
        filename = 'trueSV' // '.grid' // trim(int2str(gridcounter)) // fileXT
        call savemat(transpose(dgp%SV(:,0:T)), filename)

        filename = 'trueSTATE' // '.grid' // trim(int2str(gridcounter)) // fileXT 
        call savemat(transpose(dgp%x(:,0:T)), filename)

        filename = 'trueTHETA' // '.grid' // trim(int2str(gridcounter)) // fileXT
        call savevec(dgp%theta(0:T), filename)

        filename = 'trueLAMBDA' // '.grid' // trim(int2str(gridcounter)) // fileXT
        call savevec(dgp%lambda(0:T), filename)

        filename = 'trueNOISEVOL' // '.grid' // trim(int2str(gridcounter)) // fileXT
        call savevec(dgp%noisevol, filename)

        filename = 'trueNOISEVOL' // '.grid' // trim(int2str(gridcounter)) // fileXT
        call savevec(dgp%hinno, filename)

        filename = 'trueSIGA' // '.grid' // trim(int2str(gridcounter)) // fileXT
        call savevec((/ dgp%sigtheta /), filename)

        filename = 'trueSIGLAMBDA' // '.grid' // trim(int2str(gridcounter)) // fileXT
        call savevec((/ dgp%siglambda /), filename)

     end if


     ! allocate memory for draws
     ALLOCATE (PARTICLEweights(Nparticles,0:T),logMDD(T), DRAWlambda(Nparticles,0:T), DRAWxhat(Nx,Nparticles,0:T), DRAWxsig(Nx,Nparticles,0:T), DRAWsvol(Nparticles,Nsv,0:T), DRAWxgain(Nx,Ny,Nparticles,T), DRAWsqrtSigmaX(NsigmaX,Nparticles,0:T), STAT=status)
     IF (status /= 0) THEN
        WRITE (*,*) 'Allocation problem (draws)'
     END IF

     ALLOCATE (DRAWa(Nparticles,0:T), STAT=status)
     IF (status /= 0) THEN
        WRITE (*,*) 'Allocation problem (DRAWa)'
     END IF

     ! scale parameters
     ALLOCATE (DRAWsiglambda(Nparticles,T), DRAWsigma(Nparticles,Ny,T), DRAWhInno(Nparticles,Nsv,T), STAT=status)
     IF (status /= 0) THEN
        WRITE (*,*) 'Allocation problem (scale parameter draws)'
     END IF
     ALLOCATE (DRAWsiga(Nparticles,T), STAT=status)
     IF (status /= 0) THEN
        WRITE (*,*) 'Allocation problem (siga parameter draws)'
     END IF


     PARTICLEweights = 1.0d0 / dble(Nparticles)
     logMDD          = 0.0d0
     DRAWlambda      = 0.0d0
     DRAWa           = 0.0d0
     DRAWxhat        = 0.0d0
     DRAWxsig        = 0.0d0
     DRAWsvol        = 0.0d0
     DRAWxgain       = 0.0d0

     DRAWsiga      = 0.0d0
     DRAWsiglambda = 0.0d0
     DRAWhInno     = 0.0d0
     DRAWsigma     = 0.0d0 

     CALL plefilter(doSecondResamplingStep, dofAPFfattail, doInflationNoise, T, Ny, y, yNaN, Nparticles, PARTICLEweights, logMDD, DRAWxhat, DRAWsqrtSigmaX, DRAWxsig, DRAWxgain, Nx, NsigmaX, Nw, Ex0, sqrtVx0, p, DRAWa, a0, a0V, DRAWsiga, sigaT, sigaDof, DRAWlambda, lambda0, lambda0V, DRAWsiglambda, siglambdaT, siglambdaDof, DRAWsvol, Nsv, Eh0, Vh0, DRAWhInno, hvarT, hvarDof, DRAWsigma, sigmaT, sigmaDof, VSLarray(TID))

     ! before demeaning: prepare smoother

     ALLOCATE (PARAMweights(Nparticles))
     ALLOCATE (PARAMhInno(Nparticles,Nsv), PARAMsigma(Nparticles,Ny))
     ALLOCATE (PARAMsiglambda(Nparticles))
     ALLOCATE (PARAMsiga(Nparticles))

     PARAMweights   = PARTICLEweights(:,T)
     PARAMhInno     = DRAWhInno(:,:,T)
     PARAMsigma     = DRAWsigma(:,:,T)
     PARAMsiglambda = DRAWsiglambda(:,T)
     PARAMsiga      = DRAWsiga(:,T)


     ! STORE FILTERED STATES AND PARAMETERS
     ! DEMEAN results
     forall(j=0:T,k=1:Nparticles,i=1:Nx) DRAWxhat(i,k,j)   = dgp%x(i,j) - DRAWxhat(i,k,j)
     forall(j=0:T,k=1:Nparticles,i=1:Nsv) DRAWsvol(k,i,j)  = dgp%SV(i,j) - DRAWsvol(k,i,j)
     forall(j=1:T,k=1:Nparticles,i=1:Nsv) DRAWhinno(k,i,j) = dgp%hinno(i) - DRAWhinno(k,i,j)
     forall(j=0:T,k=1:Nparticles) DRAWa(k,j)               = dgp%Theta(j) - DRAWa(k,j)
     forall(j=0:T,k=1:Nparticles) DRAWlambda(k,j)          = dgp%Lambda(j) - DRAWlambda(k,j)

     forall(j=1:T,k=1:Nparticles) DRAWsiga(k,j)            = dgp%sigTheta - DRAWsiga(k,j)
     forall(j=1:T,k=1:Nparticles) DRAWsiglambda(k,j)       = dgp%siglambda - DRAWsiglambda(k,j)
     forall(j=1:T,k=1:Nparticles,i=1:Ny) DRAWsigma(k,i,j)   = (dgp%noisevol(i) ** 2) - DRAWsigma(k,i,j)

     ! Collect bias in filtered estimates
     forall(j=0:T,i=1:Nx) biasX(i,j,gridcounter) = sum(PARTICLEweights(:,j) * DRAWxhat(i,:,j))
     forall(j=0:T,i=1:Nsv) biasSV(i,j,gridcounter) = sum(PARTICLEweights(:,j) * DRAWsvol(:,i,j))
     forall(i=1:Nsv) biasHinno(i,gridcounter) = sum(PARTICLEweights(:,T) * DRAWhInno(:,i,T))

     forall(j=0:T) biasTheta(j,gridcounter) = sum(PARTICLEweights(:,j) * DRAWa(:,j))
     forall(j=0:T) biasLambda(j,gridcounter) = sum(PARTICLEweights(:,j) * DRAWlambda(:,j))
     biasSigTheta(gridcounter)  = sum(PARTICLEweights(:,T) * DRAWsiga(:,T))
     biasSigLambda(gridcounter) = sum(PARTICLEweights(:,T) * DRAWsiglambda(:,T))
     forall(j=1:Ny) biasNoisevol(j,gridcounter) = sum(PARTICLEweights(:,T) * DRAWsigma(:,j,T))

     ! ----------------------------------------------------------------------------
     ! CLEANUP FILTER
     ! ----------------------------------------------------------------------------

     DEALLOCATE (PARTICLEweights, DRAWlambda)
     DEALLOCATE (logMDD)
     DEALLOCATE (DRAWxgain)
     DEALLOCATE (DRAWsvol)
     DEALLOCATE (DRAWsiglambda, DRAWsigma, DRAWhInno)     
     DEALLOCATE (DRAWa, DRAWsiga)
     DEALLOCATE (DRAWxhat, DRAWxsig, DRAWsqrtSigmaX)


     IF (doSmoother) THEN
        ! ! ----------------------------------------------------------------------------
        ! ! SMOOTHER
        ! ! ----------------------------------------------------------------------------
        ALLOCATE (SMOOTHERsvol(Nsv,Nsmoother,0:T), SMOOTHERlambda(Nsmoother,0:T), SMOOTHERa(Nsmoother,0:T))
        ALLOCATE (SMOOTHERsigma(Ny,Nsmoother))
        ALLOCATE (SMOOTHERhinno(Nsv,Nsmoother))

        call hrulefill
        print *, 'STARTING w/SMOOTHER ...', gridcounter
        CALL plesmoother(doInflationNoise, Nsmoother, smootherNparticles, T, p, SMOOTHERa, SMOOTHERlambda, Nsv, SMOOTHERsvol, SMOOTHERhinno, SMOOTHERsigma, Nparticles, PARAMweights, PARAMsiga, PARAMsiglambda, PARAMhInno, PARAMsigma, a0, a0V, lambda0, lambda0V, Eh0, Vh0, Ex0, sqrtVx0, Ny, y, yNaN, Nx, NsigmaX, Nw, VSLarray(TID))

        ! debug: store smoother lambda results
        filename = 'SMOOTHERlambda.grid' // trim(int2str(gridcounter)) // fileXT
        call savemat(SMOOTHERlambda, filename)
        filename = 'TRUElambda.grid' // trim(int2str(gridcounter)) // fileXT
        call savevec(dgp%lambda, filename)

        ! ----------------------------------------------------------------------------------------------------------------------------------------------------------
        ! X SMOOTHER
        ! ----------------------------------------------------------------------------------------------------------------------------------------------------------

        ALLOCATE (SMOOTHERxhat(Nx,0:T,Nsmoother))
        SMOOTHERxhat = 0.0d0
        call particleSmootherXhat(T, Ny, y, yNaN, Nsv, Nsmoother, SMOOTHERa, p, SMOOTHERlambda, SMOOTHERsvol, SMOOTHERxhat, Nx, Nw, Ex0, sqrtVx0, SMOOTHERsigma)
        
        ! DEMEAN results
        forall(j=0:T,k=1:Nsmoother,i=1:Nx) SMOOTHERxhat(i,j,k)  = dgp%x(i,j) - SMOOTHERxhat(i,j,k)
        forall(j=0:T,k=1:Nsmoother,i=1:Nsv) SMOOTHERsvol(i,k,j) = dgp%SV(i,j) - SMOOTHERsvol(i,k,j)
        forall(j=0:T,k=1:Nsmoother) SMOOTHERa(k,j) = dgp%Theta(j) - SMOOTHERa(k,j)
        forall(j=0:T,k=1:Nsmoother) SMOOTHERlambda(k,j) = dgp%Lambda(j) - SMOOTHERlambda(k,j)


        ! Collect SMOOTHER bias 
        forall(j=0:T,i=1:Nx) smootherBiasX(i,j,gridcounter)   = mean(SMOOTHERxhat(i,j,:)) 
        forall(j=0:T,i=1:Nsv) smootherBiasSV(i,j,gridcounter) = mean(SMOOTHERsvol(i,:,j)) 
        forall(j=0:T) smootherBiasTheta(j,gridcounter)        = mean(SMOOTHERa(:,j))
        forall(j=0:T) smootherBiasLambda(j,gridcounter)       = mean(smootherLambda(:,j))

        ! clean up smoother
        DEALLOCATE (SMOOTHERsvol, SMOOTHERlambda, SMOOTHERa)
        DEALLOCATE (SMOOTHERxhat)
        DEALLOCATE (SMOOTHERsigma, SMOOTHERhinno)

        print *, '... DONE w/SMOOTHER.', gridcounter
        call hrulefill

        ! ----------------------------------------------------------------------------
        ! FINAL CLEANUP
        ! ----------------------------------------------------------------------------
     END IF ! doSmoother

     DEALLOCATE (PARAMweights)
     DEALLOCATE (PARAMsigma, PARAMhInno)     
     DEALLOCATE (PARAMsiglambda)
     DEALLOCATE (PARAMsiga)

     DEALLOCATE (dgp)
     DEALLOCATE (y, yNaN)

     !$OMP CRITICAL
     CALL progressbarcomment(dble(ompcount) / dble(Ngrid), timer, 'MONTE CARLO SIMS')
     !$OMP END CRITICAL

  end do ! Ngrid
  !$OMP END PARALLEL DO 

  ! WRITE SETTINGS
  CALL HRULEFILL
  filename = 'settings' // trim(adjustl(filext))
  OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
  WRITE(4,'(a20,a40)') 'TIME: ', timestampstr()
  WRITE(4,'(a20,I40)') 'SEED: ', seed
  WRITE(4,'(a20,a40)') 'Data: ', datalabel
  WRITE(4,'(a20,a40)') 'Model: ', modellabel
  WRITE(4,'(a60)') repeat('-',60)
  WRITE(4,'(a20,I40)') 'Nparticles: ', Nparticles
  if (doSmoother) THEN
     WRITE(4,'(a20,I40)') 'Nsmoother:', Nsmoother
     WRITE(4,'(a20,I40)') 'smootherNparticles:', smootherNparticles
  END IF
  WRITE(4,'(a20,I40)') 'p: ', p
  WRITE(4,'(a60)') repeat('-',60)
  WRITE(4,'(a20,F6.2)') 'a0: ', a0
  WRITE(4,'(a20,F6.2)') 'a0V: ', a0V
  WRITE(4,'(a20,F6.2)') 'siga0: ', siga0
  WRITE(4,'(a20,F6.2)') 'lambda0: ', lambda0
  WRITE(4,'(a20,F6.2)') 'lambda0V: ', lambda0V
  WRITE(4,'(a20,F6.2)') 'siglambda0: ', siglambda0
  do i=1,Ny
     WRITE(4,'(a20,I2,F6.2,F6.2)') 'noisevol0(0): ', i, noisevol0(i), noisevol00(i)
  end do
  WRITE(4,'(a60)') repeat('-',60)
  if (doInflationNoise) THEN
     WRITE(4,'(a60)') 'With noise in inflation'
  ELSE
     WRITE(4,'(a60)') 'WITHOUT noise in inflation'
  END IF
  if (doSecondResamplingStep) THEN
     WRITE(4,'(a60)') 'APF uses second resampling step'
  END IF
  if (dofAPFfattail .gt. 0) THEN
     WRITE(4,'(a20,i40)') 'APF MV-T dof:', dofAPFfattail
  END IF
  CLOSE(UNIT=4)
  CALL HRULEFILL

  filename  = 'lambda00' // filext
  CALL savevec(lambda00,filename)
  filename  = 'a00' // filext
  CALL savevec(a00,filename)
  filename  = 'x00' // filext
  CALL savemat(x00,filename)
  deallocate(lambda00, a00, x00)


  ! ----------------------------------------------------------------------------
  ! BEGIN: STORE FILTER-BIAS
  ! ----------------------------------------------------------------------------

  allocate(thetavec(Ngrid,T))

  do i = 1,Nx
     THIS = 'filterBiasX' // trim(int2str(i))
     filename  = trim(this) // filext
     FORALL (j=1:Ngrid) thetavec(j,:) = biasX(i,1:T,j)
     ! CALL storeEstimates(thetavec,T,Ngrid,filename)
     CALL savemat(thetavec,filename)
     WRITE (*,*) 'STORED ' // trim(this)
  end do

  do i = 1,Nsv
     THIS = 'filterBiasSV' // trim(int2str(i))
     filename  = trim(this) // filext
     FORALL (j=1:Ngrid) thetavec(j,:) =  biasSV(i,1:T,j)
     ! CALL storeEstimates(thetavec,T,Ngrid,filename)
     CALL savemat(thetavec,filename)
     WRITE (*,*) 'STORED ' // trim(this)
  end do

  THIS = 'filterBiasTheta' 
  filename  = trim(this) // filext
  FORALL (j=1:Ngrid) thetavec(j,:) =  biasTheta(1:T,j)
  ! CALL storeEstimates(thetavec,T,Ngrid,filename)
  CALL savemat(thetavec,filename)
  WRITE (*,*) 'STORED ' // trim(this)

  THIS = 'filterBiasLambda' 
  filename  = trim(this) // filext
  FORALL (j=1:Ngrid) thetavec(j,:) =  biasLambda(1:T,j)
  ! CALL storeEstimates(thetavec,T,Ngrid,filename)
  CALL savemat(thetavec,filename)
  WRITE (*,*) 'STORED ' // trim(this)

  THIS = 'filterBiasHinno' 
  filename  = trim(this) // filext
  CALL savemat(transpose(biasHinno),filename)
  WRITE (*,*) 'STORED ' // trim(this)

  THIS = 'filterBiasSigTheta' 
  filename  = trim(this) // filext
  ! CALL storeEstimates(thetavec,T,Ngrid,filename)
  CALL savevec(biasSigTheta,filename)
  WRITE (*,*) 'STORED ' // trim(this)

  THIS = 'filterBiasSigLambda' 
  filename  = trim(this) // filext
  ! CALL storeEstimates(thetavec,T,Ngrid,filename)
  CALL savevec(biasSigLambda,filename)
  WRITE (*,*) 'STORED ' // trim(this)

  THIS = 'filterBiasNoisevol' 
  filename  = trim(this) // filext
  CALL savemat(transpose(biasNoisevol),filename)
  WRITE (*,*) 'STORED ' // trim(this)

  deallocate(thetavec)

  deallocate (biasX,biasSV)
  deallocate (biasTheta,biasLambda)
  deallocate (biasSigTheta)
  deallocate (biasSigLambda)
  deallocate (biasNoisevol)
  deallocate (biasHinno)

  ! ----------------------------------------------------------------------------
  ! FINISHED: STORE FILTER-BIAS
  ! ----------------------------------------------------------------------------


  ! ----------------------------------------------------------------------------
  ! START: STORE SMOOTHERBIAS
  ! ----------------------------------------------------------------------------


  if (doSmoother) then
     allocate(thetavec(Ngrid,T))
     do i = 1,Nx
        THIS = 'smootherBiasX' // trim(int2str(i))
        filename  = trim(this) // filext
        FORALL (j=1:Ngrid) thetavec(j,:) = smootherBiasX(i,1:T,j)
        ! CALL storeEstimates(thetavec,T,Ngrid,filename)
        CALL savemat(thetavec,filename)
        WRITE (*,*) 'STORED ' // trim(this)
     end do

     do i = 1,Nsv
        THIS = 'smootherBiasSV' // trim(int2str(i))
        filename  = trim(this) // filext
        FORALL (j=1:Ngrid) thetavec(j,:) =  smootherBiasSV(i,1:T,j)
        ! CALL storeEstimates(thetavec,T,Ngrid,filename)
        CALL savemat(thetavec,filename)
        WRITE (*,*) 'STORED ' // trim(this)
     end do

     THIS = 'smootherBiasTheta' 
     filename  = trim(this) // filext
     FORALL (j=1:Ngrid) thetavec(j,:) =  smootherBiasTheta(1:T,j)
     ! CALL storeEstimates(thetavec,T,Ngrid,filename)
     CALL savemat(thetavec,filename)
     WRITE (*,*) 'STORED ' // trim(this)

     THIS = 'smootherBiasLambda' 
     filename  = trim(this) // filext
     FORALL (j=1:Ngrid) thetavec(j,:) =  smootherBiasLambda(1:T,j)
     ! CALL storeEstimates(thetavec,T,Ngrid,filename)
     CALL savemat(thetavec,filename)
     WRITE (*,*) 'STORED ' // trim(this)

     deallocate (smootherBiasX,smootherBiasSV)
     deallocate (smootherBiasTheta,smootherBiasLambda)

     deallocate(thetavec)

  end if
  ! ----------------------------------------------------------------------------
  ! FINISHED: STORE SMOOTHER-BIAS
  ! ----------------------------------------------------------------------------

  ! VSLstreams
  do i=0,NTHREADS-1
     errcode = vsldeletestream(VSLarray(i))
  end do
  deallocate(VSLarray)

  call hrulefill
  WRITE(*,*) 'DONE. BYE, BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill

  STOP

CONTAINS

  ! -----------------------------------------------------------------
  SUBROUTINE getarguments(Ngrid, T, doInflationNoise, Nparticles, doStoreDGP, doSmoother, Nsmoother, smootherNparticles, doSecondResamplingStep, dofAPFfattail)

    INTENT(INOUT) Ngrid, T, doInflationNoise, Nparticles, doSecondResamplingStep, doStoreDGP

    INTEGER :: counter, dummy
    INTEGER :: Nparticles,Nsmoother, smootherNparticles
    INTEGER :: Ngrid, T
    LOGICAL :: doSmoother,doInflationNoise
    LOGICAL :: doSecondResamplingStep
    LOGICAL :: doStoreDGP
    INTEGER, INTENT(INOUT) :: dofAPFfattail
    CHARACTER(len=32) :: arg

    counter = 0

    ! Ngrid
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') Ngrid
    END IF

    ! T
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') T
    END IF

    ! doInflationNoise
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          doInflationNoise = .true.
       else
          doInflationNoise = .false.
       end if
    END IF

    ! Nparticles
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') Nparticles
    END IF

    ! doStoreDGP
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          doStoreDGP = .true.
       else
          doStoreDGP = .false.
       end if

       ! doSmoother
       counter = counter + 1
       IF (command_argument_count() >= counter) THEN
          CALL get_command_argument(counter, arg) 
          READ(arg, '(i20)') dummy
          if (dummy > 0) then
             doSmoother = .true.
          else
             doSmoother = .false.
          end if
       END IF

       ! Nsmoother
       counter = counter + 1
       IF (command_argument_count() >= counter) THEN
          CALL get_command_argument(counter, arg) 
          READ(arg, '(i20)') Nsmoother
       END IF

       ! smootherNparticles
       counter = counter + 1
       IF (command_argument_count() >= counter) THEN
          CALL get_command_argument(counter, arg) 
          READ(arg, '(i20)') smootherNparticles
       END IF

    END IF

    ! doSecondResamplingStep
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       if (dummy > 0) then
          doSecondResamplingStep = .true.
       else
          doSecondResamplingStep = .false.
       end if
    END IF

    ! dofAPFfattail
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dofAPFfattail
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



! @\newpage\subsection{plefilter}@
SUBROUTINE plefilter(doSecondResamplingStep, dofAPFfattail, doInflationNoise, T, Ny, y, yNaN, Nparticles, PARTICLEweights, logMDD, DRAWxhat, DRAWsqrtSigmaX, DRAWxsig, DRAWxgain, Nx, NsigmaX, Nw, Ex0, sqrtVx00, p, DRAWa, a0, a0V, DRAWsiga, sigaT, sigaDof, DRAWlambda, lambda0, lambda0V, DRAWsiglambda, siglambdaT, siglambdaDof, DRAWsvol, Nsv, Eh0, Vh0, DRAWhInno, hvarT, hvarDof, DRAWsigma, sigmaT, sigmaDof, VSLstream)

  USE embox, only : savemat, savevec, int2str, mean, hrulefill

  use gibbsbox, only : drawNDXsysresample, igammaDraws
  use statespacebox, only : DLYAPsqrt
  use blaspack, only : pi, vechU,ivechU, vech, ivech, eye, symmetric, qrot, qrquery

  use cambridgebox, only : drawdeltatruncnorm

  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE


  INTENT(INOUT) :: PARTICLEweights, logMDD, DRAWa, DRAWlambda, DRAWxhat, DRAWsqrtSigmaX, DRAWxsig, DRAWxgain, DRAWsvol, DRAWsiga, DRAWsiglambda, DRAWhInno, DRAWsigma, VSLstream
  INTENT(IN)    :: T,Ny,y,yNaN, Nx, NsigmaX, Nw, Ex0, sqrtVx00, Nsv,Eh0,Vh0, Nparticles, p, lambda0, lambda0V, a0, a0V, sigaT, sigaDof, siglambdaT, siglambdaDof, hvarT, hvarDof, sigmaT, sigmaDof

  INTEGER :: J, I, K, N, T, Nparticles, Nx, Ny, Nsv, NsigmaX, p, Nw, Nsurveys

  ! OPEN MP
  INTEGER :: TID


  type(progresstimer) :: timer

  logical, intent(in) :: doInflationNoise

  double precision, parameter :: minParticleWeight = 1.0d-12

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(T)              :: logMDD
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: DRAWlambda
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: DRAWa
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: DRAWxsig
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: DRAWxhat
  DOUBLE PRECISION, DIMENSION(NsigmaX,Nparticles, 0:T)  :: DRAWsqrtSigmaX
  DOUBLE PRECISION, DIMENSION(Nx,Ny,Nparticles,T) :: DRAWxgain
  DOUBLE PRECISION, DIMENSION(Nparticles,Nsv,0:T) :: DRAWsvol

  ! APF llf correction
  DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
  DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax


  ! scale parameters
  DOUBLE PRECISION, DIMENSION(Nparticles,T)     :: DRAWsiga
  DOUBLE PRECISION, DIMENSION(Nparticles,T)     :: DRAWsiglambda
  DOUBLE PRECISION, DIMENSION(Nparticles,Nsv,T) :: DRAWhInno
  DOUBLE PRECISION, DIMENSION(Nsv,Nparticles)   :: hInno ! helper variable, note the transpose
  DOUBLE PRECISION, DIMENSION(Nparticles,Ny,T)  :: DRAWsigma

  ! particles
  DOUBLE PRECISION :: xposterior(Nx,Nparticles), xsig(Nx,Nparticles), h(Nsv, Nparticles), SVol(Nsv,Nparticles), llf(Nparticles)
  DOUBLE PRECISION :: kernelweights(Nparticles) ! for unnormalized weights
  DOUBLE PRECISION  :: kernelsum, loglikemax
  INTEGER :: ndx(Nparticles)
  DOUBLE PRECISION :: shufflevec(Nparticles)

  ! state space objects
  DOUBLE PRECISION :: xprior(Nx), SigmaX(Nx,Nx), logdetSigmaY

  DOUBLE PRECISION :: lambda(Nparticles), lambdaPREV(Nparticles), lambda0, lambda0V
  DOUBLE PRECISION :: adrift(Nparticles), aPREV(Nparticles),  a0, a0V

  ! scale parameters
  DOUBLE PRECISION :: sigaT, siglambdaT, hvarT(Nsv), sigmaT(Ny)
  INTEGER :: sigaDof, siglambdaDof, hvarDof(Nsv), sigmaDof(Ny)

  DOUBLE PRECISION :: PREVsiglambdaT(Nparticles), PREVhvarT(Nparticles,Nsv), PREVsigmaT(Nparticles,Ny)
  INTEGER :: PREVsiglambdaDof, PREVhvarDof(Nsv), PREVsigmaDof(Ny)
  DOUBLE PRECISION :: lambdaDELTA(Nparticles), hDELTA(Nsv,Nparticles), resid(Ny,Nparticles)
  double precision :: xdraw(Nx,Nparticles) ! used for constructing resid

  ! AR1:
  DOUBLE PRECISION :: PREVsigaT(Nparticles)
  INTEGER :: PREVsigaDof
  DOUBLE PRECISION :: aDELTA(Nparticles)


  ! SQRT objects
  DOUBLE PRECISION :: sqrtSigmaX(Nx,Nx), sqrtSigmaY(Ny,Ny), Kgain(Nx,Ny,Nparticles), qrR(Ny+Nx+Nw,Ny+Nx)
  DOUBLE PRECISION, DIMENSION(NsigmaX, Nparticles) :: vecSqrtSigmaX
  INTEGER :: qrLwork

  DOUBLE PRECISION :: Ex0(Nx), sqrtVx0(Nx,Nx), sqrtVx00(Nx,Nx), A(Nx,Nx), B(Nx,Nw), Bsv(Nx,Nw), C(Ny,Nx,T), sqrtR(Ny,Ny), ytilde(Ny)
  ! DOUBLE PRECISION :: infgapcompanion(p,p), infgapcompanionh(p,p)
  DOUBLE PRECISION :: ygap0sqrtvariance(2*p,2*p), gapshock0loadings(2*p,Nw), gaptransition(2*p,2*p)

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
  INTEGER :: ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap
  INTEGER :: ndxgap(2*p)

  logical, intent(in) :: doSecondResamplingStep
  integer, intent(in) :: dofAPFfattail

  Nsurveys = Ny - 1
  minSVh   = log(0.001d0 ** 2)

  ! state-vector indices
  ndxTrendRE = 1
  ndxGapRE   = 2
  ndxTrendSI = 2+p
  ndxGapSI   = ndxTrendSI + 1

  ndxGapREstart = ndxGapRE
  ndxGapREstop  = ndxGapREstart + p - 1
  ndxGapSIstart = ndxGapSI
  ndxGapSIstop  = ndxGapSIstart + p - 1

  ndxgap = (/ ndxgapREstart : ndxGapREstop, ndxgapSIstart : ndxGapSIstop /)

  shockndxTrend = 1
  shockndxGap   = 2

  ! init sufficient statistics of scale parameters
  PREVsigaT   = sigaT
  PREVsigaDof = sigaDof

  PREVsiglambdaT   = siglambdaT
  PREVsiglambdaDof = siglambdaDof

  forall(k=1:Nparticles,j=1:Nsv) PREVhvarT(k,j) = hvarT(j)
  forall(j=1:Nsv)                PREVhvarDof(j) = hvarDof(j)

  forall(k=1:Nparticles,j=1:Ny) PREVsigmaT(k,j) = sigmaT(j)
  forall(j=1:Ny)                PREVsigmaDof(j) = sigmaDof(j)


  ! prepare state space
  ! A
  A = 0.0d0
  ! unit root in trend
  A(ndxTrendRE,ndxTrendRE) = 1.0d0
  ! kompanion for gap
  IF (p > 1) THEN
     FORALL(j=1:(p-1)) A(ndxGapRE+j,ndxGapRE-1+j) = 1.0d0
  END IF

  ! B
  B                           = 0.0d0
  B(ndxTrendRE,shockndxTrend) = 1.0d0
  B(ndxGapRE,shockndxGap)     = 1.0d0

  ! C
  C         = 0.0d0
  ! inflation
  C(1,ndxTrendRE,:)    = 1.0d0
  C(1,ndxGapRE,:)      = 1.0d0
  ! surveys
  C(2:Ny,ndxTrendSI,:) = 1.0d0


  ! Time 0 particles

  ! SV0
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, h, 0.0d0, 1.0d0)
  FORALL (j=1:Nsv,k=1:Nparticles) h(j,k) = Eh0(j) + sqrt(Vh0(j)) * h(j,k)  
  SVol = exp(h * 0.5d0)

  ! Lambda(0)
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, lambdaPREV, lambda0, sqrt(lambda0V))
  lambda = lambdaPREV
  do while (ANY(lambda < 0.0d0) .OR. ANY(lambda > 1.0d0))
     errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, lambdaPREV, lambda0, sqrt(lambda0V)) ! a bit brute force to *always* redraw N random variables ..
     where (lambda < 0.0d0)   lambda = lambdaPREV
     where (lambda > 1.0d0)   lambda = lambdaPREV
  end do
  
  TID = 0
  !$ TID = OMP_GET_THREAD_NUM()
  
  ! a(0)
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, aPREV, a0, sqrt(a0V))
  adrift  = aPREV
  do while (ANY(adrift < -1.0d0) .OR. ANY(adrift > 1.0d0))
     errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, aPREV, a0, sqrt(a0V)) ! a bit brute force to *always* redraw N random variables ..
     where (adrift < -1.0d0)  adrift = aPREV
     where (adrift > 1.0d0)   adrift = aPREV
  end do

  ! RB priors for linear states
  FORALL(k=1:Nparticles) xposterior(:,k) = Ex0

  ! prepare prior variance of linear states
  vecSqrtSigmaX       = 0.0d0
  ! x !$OMP PARALLEL DO DEFAULT(NONE), &
  ! x !$OMP& SHARED(vecSqrtSigmaX, p, SVol, lambda, adrift, Ny, Nx, Nw, Nsv, Nparticles, sqrtVx00, ndxgap, shockndxGap), &
  ! x !$OMP& PRIVATE(gaptransition,gapshock0loadings, sqrtVx0, ygap0sqrtvariance, errcode) 
  DO k=1,Nparticles

     sqrtVx0         = transpose(sqrtVx00)

     if (abs(adrift(k)) > 1.0d-4 .AND. lambda(k) > 1.0d-4) then 
        gaptransition                  = 0.0d0
        gaptransition(1:p,1:p)         = adrift(k)
        gaptransition(p+1:2*p,1:p)     = (1 - lambda(k)) * adrift(k)
        gaptransition(p+1:2*p,p+1:2*p) = lambda(k) * adrift(k)

        ! Fill in unconditional variance of stationary states
        ! allow for trendshockslopes, thought they are all zero here
        gapshock0loadings                              = 0.0d0
        gapshock0loadings(1,shockndxGap)               = SVol(shockndxGap,k)
        gapshock0loadings(p+1,shockndxGap)             = (1 - lambda(k)) * SVol(shockndxGap,k)

        CALL DLYAPsqrt(ygap0sqrtvariance, gaptransition, gapshock0loadings, 2 * p, Nw, errcode) 
        if (errcode /= 0) then
           write (*,*) 'DLYAP error (ygap0sqrtvariance -- init plefilter)', errcode
           call savemat(gaptransition, 'gaptransition.debug')
           call savemat(gapshock0loadings, 'gapshock0loadings.debug')
           call savemat(ygap0sqrtvariance, 'ygap0sqrtvariance.debug')
           stop 1
        end if
        sqrtVx0(ndxgap,ndxgap) = ygap0sqrtvariance
     end if
     vecSqrtSigmaX(:,k)     = vechU(sqrtVx0,Nx)

  END DO
  ! x !$OMP END PARALLEL DO 


  FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,0)  = vecSqrtSigmaX(i,k)
  FORALL(i=1:Nx,k=1:Nparticles)  DRAWxhat(i,k,0)            = xposterior(i,k) 


  FORALL (i=1:Nsv,k=1:Nparticles) DRAWsvol(k,i,0) = SVol(i,k) ! transpose
  DRAWlambda(:,0) = lambda
  DRAWa(:,0)      = adrift
  PARTICLEweights  = 1.0d0 / dble(Nparticles)
  logMDD         = 0.0d0

  ! uniform draws for systematic resampling
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


  ! workspace query for qr decomposition
  qrR     = 0.0d0
  qrlwork = qrquery(qrR)

  if (TID .eq. 0) CALL initprogressbar(timer, 15.0d0)
  DO j=1,T

     Nynonan = count(.not. yNaN(:,j))
     ! print *, 'plefilter for t=', j, ' with ', Nynonan, 'obs' 
     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: APF RESAMPLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! Sigma
     i=1
     if (doInflationNoise) then
        call igammaDraws(DRAWsigma(:,i,j), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream)
     else
        DRAWsigma(:,i,j) = 0.0d0
        PREVsigmaT(:,i)  = 0.0d0
     end if
     DO i=2,Ny 
        call igammaDraws(DRAWsigma(:,i,j), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream) 
     END DO

     ! x !$OMP PARALLEL DO DEFAULT(NONE), SCHEDULE(STATIC), &
     ! x !$OMP& SHARED(xposterior, xsig, vecSqrtSigmaX, lambda, adrift, SVol, DRAWsigma, llf, Nparticles, Nsv, j, y, yNaN, Ny, Nx, Nw, Nynonan, Nsurveys, p, ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap), &
     ! x !$OMP& SHARED(dofAPFfattail), &
     ! x !$OMP& FIRSTPRIVATE(A,B,C,qrLWORK), &
     ! x !$OMP& PRIVATE(ytilde, sqrtR, Bsv, logdetSigmaY, errcode, xprior, SigmaX, sqrtSigmaX, sqrtSigmaY, qrR, TID) 


     DO k = 1,Nparticles

        ! TID = 0
        ! x !$ TID = OMP_GET_THREAD_NUM()


        ! 2) Fill Particles into state space

        ! update A 

        A(ndxGapRE,ndxGapRE)           = adrift(k)
        A(ndxtrendSI,ndxtrendRE)       = 1 - lambda(k)
        A(ndxtrendSI,ndxtrendSI)       = lambda(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - lambda(k)) * adrift(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = lambda(k) * adrift(k)

        ! update sqrtR
        sqrtR = 0.0d0
        forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i,j)) 


        ! update C
        FORALL (i=1:Nsurveys) C(1+i,ndxGapSI,j) = adrift(k) ** i
        ! zero out missing obs
        DO i = 1, Ny
           if (yNaN(i,j)) then 
              C(i,:,j)     = 0.0d0 
              sqrtR(i,:)   = 0.0d0
           end if
        END DO


        ! update B
        B(ndxtrendSI,shockndxTrend)    = 1 - lambda(k)
        B(ndxgapSI,shockndxGap)        = 1 - lambda(k)

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
        if (dofAPFfattail .gt. 0) then
           ! use t distribution with dofAPFtail
           ! note: only kernel of MV-T needed
           llf(k)       = -0.5d0 * (logdetSigmaY + (dofAPFfattail + Nynonan) * log(1.0d0 + sum(ytilde ** 2) / dofAPFfattail))
        else 
           llf(k)       = -0.5d0 * (Nynonan * log(2.0d0 * pi) + logdetSigmaY + sum(ytilde ** 2))
        end if


        ! DO NOT STORE POSTERIORS FOR APF STEP
        ! vecSqrtSigmaX(:,k) = vechU(sqrtSigmaX,Nx)
        ! ! store state variances 
        ! call dsyrk('u', 't', Nx, Nx, 1.0d0, sqrtSigmaX, Nx, 0.0d0, SigmaX, Nx)
        ! forall (i=1:Nx) xsig(i,k) = sqrt(SigmaX(i,i))


     END DO ! k particles
     ! x !$OMP END PARALLEL DO 


     if (Nynonan > 0) then

        ! Reweight particles for next round   
        loglikeAPFmax      = maxval(llf)
        llf                = llf - loglikeAPFmax

        APFlike = exp(llf)
        if (doSecondResamplingStep) then
           APFkernelweights     = APFlike / dble(Nparticles) 
           APFkernelsum         = sum(APFkernelweights)
           PARTICLEweights(:,j) = APFkernelweights / APFkernelsum

        else
           APFkernelweights     = APFlike * PARTICLEweights(:,j-1)
           APFkernelsum         = sum(APFkernelweights)
           PARTICLEweights(:,j) = APFkernelweights / APFkernelsum
        end if

        ! where (PARTICLEweights(:,j) < minParticleWeight) PARTICLEweights(:,j) = minParticleWeight
        ! PARTICLEweights(:,j) = PARTICLEweights(:,j) / sum(PARTICLEweights(:,j))

        call drawNDXsysresample(ndx, Nparticles, PARTICLEweights(:,j), Nparticles, uniformdraws(1,j))

        FORALL(k=1:Nparticles) shufflevec(k) = APFlike(ndx(k))
        APFlike = shufflevec

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

        FORALL(k=1:Nparticles) shufflevec(k) = lambda(ndx(k))
        lambda = shufflevec

        FORALL(k=1:Nparticles) shufflevec(k) = adrift(ndx(k))
        adrift = shufflevec

        ! reshuffle sufficient statistics for scale parameters
        FORALL(k=1:Nparticles) shufflevec(k) = PREVsigaT(ndx(k))
        PREVsigaT = shufflevec

        FORALL(k=1:Nparticles) shufflevec(k) = PREVsiglambdaT(ndx(k))
        PREVsiglambdaT = shufflevec

        DO i=1,Nsv
           FORALL(k=1:Nparticles) shufflevec(k) = PREVhvarT(ndx(k),i)
           PREVhvarT(:,i) = shufflevec
        END DO

        DO i=1,Ny
           FORALL(k=1:Nparticles) shufflevec(k) = PREVsigmaT(ndx(k),i)
           PREVsigmaT(:,i) = shufflevec
        END DO

        DO i=1,Ny
           FORALL(k=1:Nparticles) shufflevec(k) = DRAWsigma(ndx(k),i,j)
           DRAWsigma(:,i,j) = shufflevec
        END DO

     end if

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: APF STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! 0) draw scale parameters
     ! siga
     call igammaDraws(DRAWsiga(:,j), Nparticles, PREVsigaT, PREVsigaDof, VSLstream)
     DRAWsiga(:,j) = sqrt(DRAWsiga(:,j))
     ! siglambda
     call igammaDraws(DRAWsiglambda(:,j), Nparticles, PREVsiglambdaT, PREVsiglambdaDof, VSLstream)
     DRAWsiglambda(:,j) = sqrt(DRAWsiglambda(:,j))
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


     lambdaPREV  = lambda
     lambdaDELTA = drawdeltatruncnorm(Nparticles, lambdaPREV, DRAWsiglambda(:,j), 0.0d0, 1.0d0, VSLstream)
     lambda      = lambdaPREV + lambdaDELTA

     aPREV      = adrift
     aDELTA     = drawdeltatruncnorm(Nparticles, aPREV, DRAWsiga(:,j), -1.0d0, 1.0d0, VSLstream)
     adrift     = aPREV + aDELTA

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! x !$OMP PARALLEL DO SHARED(xposterior, xsig, vecSqrtSigmaX, lambda, adrift, SVol, Kgain, DRAWsigma, llf, Nparticles, Nsv, j, y, yNaN, Ny, Nx, Nw, Nynonan, Nsurveys, p, ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap) FIRSTPRIVATE(A,B,C,qrLWORK) PRIVATE(ytilde, sqrtR, Bsv, logdetSigmaY, errcode, xprior, SigmaX, sqrtSigmaX, sqrtSigmaY, qrR, TID) DEFAULT(NONE) SCHEDULE(STATIC)


     DO k = 1,Nparticles

        ! TID = 0
        ! x !$ TID = OMP_GET_THREAD_NUM()

        xprior      = 0.0d0
        ! SigmaX      = 0.0d0
        ! call ivech(SigmaX, vecSigmaX(:,k))

        ! 2) Fill Particles into state space

        ! update A 

        A(ndxGapRE,ndxGapRE)           = adrift(k)
        A(ndxtrendSI,ndxtrendRE)       = 1 - lambda(k)
        A(ndxtrendSI,ndxtrendSI)       = lambda(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - lambda(k)) * adrift(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = lambda(k) * adrift(k)

        ! update sqrtR
        sqrtR = 0.0d0
        forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i,j)) 


        ! update C
        FORALL (i=1:Nsurveys) C(1+i,ndxGapSI,j) = adrift(k) ** i

        ! zero out missing obs
        DO i = 1, Ny
           if (yNaN(i,j)) then 
              C(i,:,j)     = 0.0d0 
              sqrtR(i,:)   = 0.0d0
           end if
        END DO


        ! update B
        B(ndxtrendSI,shockndxTrend)    = 1 - lambda(k)
        B(ndxgapSI,shockndxGap)        = 1 - lambda(k)

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
     ! x !$OMP END PARALLEL DO 

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! sample noise from Kalman-posteriors for updating SIGMA (further below)
     errcode  = vdrnggaussian(VSLmethodGaussian, VSLstream, Nx * Nparticles, xdraw, 0.0d0, 1.0d0)
     ! scale xdraw by vecSqrtSigmaX
     do k=1,Nparticles
        call DTPMV('u','t','n',Nx,vecSqrtSigmaX(:,k),xdraw(:,k),1)
     end do
     ! add mean to xdraw
     forall (k=1:Nparticles,i=1:Nx) xdraw(i,k) = xposterior(i,k) + xdraw(i,k)

     ! resid = y - C x
     forall (k=1:Nparticles,i=1:Ny) resid(i,k) = y(i,j)
     do k=1,Nparticles
        ! update C
        FORALL (i=1:Nsurveys) C(1+i,ndxGapSI,j) = adrift(k) ** i
        call dgemv('n',Ny,Nx,-1.0d0,C(:,:,j),Ny,xdraw(:,k),1,1.0d0,resid(:,k),1) 
        ! zero out when there were missing values
        where (yNaN(:,j)) resid(:,k) = 0.0d0
     end do

     ! done sampling noise


     ! MDD
     loglikemax        = maxval(llf)
     llf               = llf - loglikemax
     kernelweights     = exp(llf) / APFlike 
     kernelsum         = sum(kernelweights)
     logMDD(j)         = log(kernelsum / Nparticles) + log(APFkernelsum) + loglikemax ! note: loglikeAPFmax correction should be added to APFlike and subtracted from kernelweights (b/o APFlike in denominator), thus cancelling each other

     ! Store NON-reweighted statistics
     FORALL(i=1:Nsv,k=1:Nparticles) DRAWsvol(k,i,j)   = SVol(i,k) ! note the transpose
     DRAWlambda(:,j)   = lambda
     DRAWa(:,j)        = adrift


     FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,j)  = vecSqrtSigmaX(i,k)
     FORALL(i=1:Nx,k=1:Nparticles) DRAWxhat(i,k,j)             = xposterior(i,k)
     FORALL(i=1:Nx,k=1:Nparticles) DRAWxsig(i,k,j)             = xsig(i,k) 
     FORALL(i=1:Nx,n=1:Ny,k=1:Nparticles) DRAWxgain(i,n,k,j)   = Kgain(i,n,k) ! Kprime(n,i,k)


     ! propagate sufficient statistics
     if (Nynonan > 0) then ! nothing to propagate if there was no observed data
        FORALL(k=1:Nparticles)         PREVsigaT(k)         = PREVsigaT(k)      + aDELTA(k) ** 2
        FORALL(k=1:Nparticles)         PREVsiglambdaT(k)    = PREVsiglambdaT(k) + lambdaDELTA(k) ** 2
        FORALL(k=1:Nparticles,i=1:Nsv) PREVhvarT(k,i)       = PREVhvarT(k,i)    + hDELTA(i,k) ** 2
        if (doInflationNoise) then
           i = 1
           FORALL(k=1:Nparticles)  PREVsigmaT(k,i)  = PREVsigmaT(k,i)   + resid(i,k) ** 2
        end if
        FORALL(k=1:Nparticles,i=2:Ny)  PREVsigmaT(k,i)      = PREVsigmaT(k,i)   + resid(i,k) ** 2 ! note: missing obs handled by zero values of resid

        PREVsigaDof      = PREVsigaDof + 1
	PREVsiglambdaDof = PREVsiglambdaDof + 1
        PREVhvarDof      = PREVhvarDof + 1
        WHERE (.NOT. yNaN(:,j))
           PREVsigmaDof     = PREVsigmaDof + 1 
        END WHERE

        ! particles weights
        PARTICLEweights(:,j) = kernelweights / kernelsum
        ! where (PARTICLEweights(:,j) < minParticleWeight) PARTICLEweights(:,j) = minParticleWeight
        ! PARTICLEweights(:,j) = PARTICLEweights(:,j) / sum(PARTICLEweights(:,j))

        ! resample
        if (doSecondResamplingStep) then

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
           ! SVol needs to be reshuffled, to prep the APF resample step
           DO i=1,Nsv 
              FORALL(k=1:Nparticles) shufflevec(k) = SVol(i,ndx(k))
              SVol(i,:) = shufflevec
           END DO

           FORALL(k=1:Nparticles) shufflevec(k) = lambda(ndx(k))
           lambda = shufflevec

           FORALL(k=1:Nparticles) shufflevec(k) = adrift(ndx(k))
           adrift = shufflevec

           ! reshuffle sufficient statistics for scale parameters
           FORALL(k=1:Nparticles) shufflevec(k) = PREVsigaT(ndx(k))
           PREVsigaT = shufflevec

           FORALL(k=1:Nparticles) shufflevec(k) = PREVsiglambdaT(ndx(k))
           PREVsiglambdaT = shufflevec

           DO i=1,Nsv
              FORALL(k=1:Nparticles) shufflevec(k) = PREVhvarT(ndx(k),i)
              PREVhvarT(:,i) = shufflevec
           END DO

           DO i=1,Ny
              FORALL(k=1:Nparticles) shufflevec(k) = PREVsigmaT(ndx(k),i)
              PREVsigmaT(:,i) = shufflevec
           END DO

        end if ! doSecondResamplingStep

     end if


     if (TID .eq. 0) CALL progressbarcomment(dble(j) / dble(T), timer, 'Particle Step')

  END DO ! j=1,T

END SUBROUTINE plefilter

! @\newpage\subsection{particlefilter}@
SUBROUTINE particlefilter(doSecondResamplingStep, dofAPFfattail, T, Ny, y, yNaN, Nparticles, PARTICLEweights, logMDD, DRAWxhat, DRAWsqrtSigmaX, Nx, NsigmaX, Nw, Ex0, sqrtVx00, p, DRAWa, a0, a0V, siga, DRAWlambda, lambda0, lambda0V, siglambda, DRAWh, Nsv, Eh0, Vh0, hInno, sigma, VSLstream)

  USE embox, only : savemat, savevec, int2str, mean, hrulefill

  use gibbsbox, only : drawNDXsysresample
  use statespacebox, only : DLYAPsqrt
  use blaspack, only : pi, vechU,ivechU, vech, ivech, eye, symmetric, qrot, qrquery

  use cambridgebox, only : drawdeltatruncnorm1
  use vslbox
  ! use omp_lib
  ! use timerbox

  IMPLICIT NONE



  INTENT(INOUT) :: PARTICLEweights, logMDD, DRAWa, DRAWlambda, DRAWxhat, DRAWsqrtSigmaX, DRAWh, VSLstream !, timer
  INTENT(IN)    :: T,Ny,y,yNaN, Nx, NsigmaX, Nw, Ex0, sqrtVx00, Nsv,Eh0,Vh0, Nparticles, p, lambda0, lambda0V, a0, a0V
  ! declared below -- INTENT(IN)    :: siga, siglambda, hInno, sigma

  INTEGER :: J, I, K, T, Nparticles, Nx, Ny, Nsv, NsigmaX, p, Nw, Nsurveys

  ! ! OPEN MP
  ! INTEGER :: TID


  ! type(progresstimer) :: timer

  double precision, parameter :: minParticleWeight = 1.0d-12

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(T)              :: logMDD
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: DRAWlambda
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: DRAWa
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: DRAWxhat
  DOUBLE PRECISION, DIMENSION(NsigmaX,Nparticles, 0:T)  :: DRAWsqrtSigmaX
  DOUBLE PRECISION, DIMENSION(Nsv,Nparticles,0:T) :: DRAWh

  ! APF llf correction
  DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
  DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax


  ! scale parameters
  DOUBLE PRECISION, INTENT(IN)                 :: siga
  DOUBLE PRECISION, INTENT(IN)                 :: siglambda
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nsv) :: hInno 
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Ny)  :: sigma

  ! particles
  DOUBLE PRECISION :: xposterior(Nx,Nparticles), h(Nsv, Nparticles), SVol(Nsv,Nparticles), llf(Nparticles)
  DOUBLE PRECISION :: kernelweights(Nparticles) ! for unnormalized weights
  DOUBLE PRECISION  :: kernelsum, loglikemax
  INTEGER :: ndx(Nparticles)
  DOUBLE PRECISION :: shufflevec(Nparticles)

  ! state space objects
  DOUBLE PRECISION :: xprior(Nx), logdetSigmaY

  DOUBLE PRECISION :: lambda(Nparticles), lambdaPREV(Nparticles), lambda0, lambda0V
  DOUBLE PRECISION :: adrift(Nparticles), aPREV(Nparticles),  a0, a0V

  DOUBLE PRECISION :: lambdaDELTA(Nparticles), hDELTA(Nsv,Nparticles)

  ! AR1:
  DOUBLE PRECISION :: aDELTA(Nparticles)


  ! SQRT objects
  DOUBLE PRECISION :: sqrtSigmaX(Nx,Nx), sqrtSigmaY(Ny,Ny), Kgain(Nx,Ny,Nparticles), qrR(Ny+Nx+Nw,Ny+Nx)
  DOUBLE PRECISION, DIMENSION(NsigmaX, Nparticles) :: vecSqrtSigmaX
  INTEGER :: qrLwork

  DOUBLE PRECISION :: Ex0(Nx), sqrtVx0(Nx,Nx), sqrtVx00(Nx,Nx), A(Nx,Nx), B(Nx,Nw), Bsv(Nx,Nw), C(Ny,Nx,T), sqrtR(Ny,Ny), ytilde(Ny)
  ! DOUBLE PRECISION :: infgapcompanion(p,p), infgapcompanionh(p,p)
  DOUBLE PRECISION :: ygap0sqrtvariance(2*p,2*p), gapshock0loadings(2*p,Nw), gaptransition(2*p,2*p) 

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
  INTEGER :: ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap
  INTEGER :: ndxgap(2*p)

  logical, intent(in) :: doSecondResamplingStep
  integer, intent(in) :: dofAPFfattail

  ! CHARACTER (LEN=100) :: filename

  Nsurveys = Ny - 1
  minSVh   = log(0.001d0 ** 2)

  ! state-vector indices
  ndxTrendRE = 1
  ndxGapRE   = 2
  ndxTrendSI = 2+p
  ndxGapSI   = ndxTrendSI + 1

  ndxGapREstart = ndxGapRE
  ndxGapREstop  = ndxGapREstart + p - 1
  ndxGapSIstart = ndxGapSI
  ndxGapSIstop  = ndxGapSIstart + p - 1

  ndxgap = (/ ndxgapREstart : ndxGapREstop, ndxgapSIstart : ndxGapSIstop /)

  shockndxTrend = 1
  shockndxGap   = 2

  ! prepare state space
  ! A
  A = 0.0d0
  ! unit root in trend
  A(ndxTrendRE,ndxTrendRE) = 1.0d0
  ! kompanion for gap
  IF (p > 1) THEN
     FORALL(j=1:(p-1)) A(ndxGapRE+j,ndxGapRE-1+j) = 1.0d0
  END IF

  ! B
  B                           = 0.0d0
  B(ndxTrendRE,shockndxTrend) = 1.0d0
  B(ndxGapRE,shockndxGap)     = 1.0d0

  ! C
  C         = 0.0d0
  ! inflation
  C(1,ndxTrendRE,:)    = 1.0d0
  C(1,ndxGapRE,:)      = 1.0d0
  ! surveys
  C(2:Ny,ndxTrendSI,:) = 1.0d0


  ! Time 0 particles

  ! SV0
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, h, 0.0d0, 1.0d0)
  FORALL (j=1:Nsv,k=1:Nparticles) h(j,k) = Eh0(j) + sqrt(Vh0(j)) * h(j,k)  
  SVol = exp(h * 0.5d0)

  ! Lambda(0)
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, lambdaPREV, lambda0, sqrt(lambda0V))
  lambda = lambdaPREV
  do while (ANY(lambda < 0.0d0) .OR. ANY(lambda > 1.0d0))
     errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, lambdaPREV, lambda0, sqrt(lambda0V)) ! a bit brute force to *always* redraw N random variables ..
     where (lambda < 0.0d0)   lambda = lambdaPREV
     where (lambda > 1.0d0)   lambda = lambdaPREV
  end do

  ! a(0)
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, aPREV, a0, sqrt(a0V))
  adrift  = aPREV
  do while (ANY(adrift < -1.0d0) .OR. ANY(adrift > 1.0d0))
     errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, aPREV, a0, sqrt(a0V)) ! a bit brute force to *always* redraw N random variables ..
     where (adrift < -1.0d0)  adrift = aPREV
     where (adrift > 1.0d0)   adrift = aPREV
  end do

  ! RB priors for linear states
  FORALL(k=1:Nparticles) xposterior(:,k) = Ex0

  ! prepare prior variance of linear states
  vecSqrtSigmaX       = 0.0d0

  DO k=1,Nparticles

     sqrtVx0         = transpose(sqrtVx00)

     if (abs(adrift(k)) > 1.0d-4 .AND. lambda(k) > 1.0d-4) then 
        gaptransition                  = 0.0d0
        gaptransition(1:p,1:p)         = adrift(k)
        gaptransition(p+1:2*p,1:p)     = (1 - lambda(k)) * adrift(k)
        gaptransition(p+1:2*p,p+1:2*p) = lambda(k) * adrift(k)

        ! Fill in unconditional variance of stationary states
        ! allow for trendshockslopes, thought they are all zero here
        gapshock0loadings                              = 0.0d0
        gapshock0loadings(1,shockndxGap)               = SVol(shockndxGap,k)
        gapshock0loadings(p+1,shockndxGap)             = (1 - lambda(k)) * SVol(shockndxGap,k)

        CALL DLYAPsqrt(ygap0sqrtvariance, gaptransition, gapshock0loadings, 2 * p, Nw, errcode) 
        if (errcode /= 0) then
           write (*,*) 'DLYAP error (ygap0sqrtvariance -- init particlefilter)', errcode
           call savemat(gaptransition, 'gaptransition.debug')
           call savemat(gapshock0loadings, 'gapshock0loadings.debug')
           call savemat(ygap0sqrtvariance, 'ygap0sqrtvariance.debug')
           stop 1
        end if
        sqrtVx0(ndxgap,ndxgap) = ygap0sqrtvariance
     end if
     vecSqrtSigmaX(:,k)     = vechU(sqrtVx0,Nx)

  END DO


  FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,0)  = vecSqrtSigmaX(i,k)
  FORALL(i=1:Nx,k=1:Nparticles)  DRAWxhat(i,k,0)            = xposterior(i,k) 


  FORALL (i=1:Nsv,k=1:Nparticles) DRAWh(i,k,0) = h(i,k) 
  DRAWlambda(:,0)  = lambda
  DRAWa(:,0)       = adrift
  PARTICLEweights   = 1.0d0 / dble(Nparticles)
  logMDD          = 0.0d0

  ! uniform draws for systematic resampling
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


  ! workspace query for qr decomposition
  qrR     = 0.0d0
  qrlwork = qrquery(qrR)

  ! CALL initprogressbar(timer, 15.0d0)
  DO j=1,T

     Nynonan = count(.not. yNaN(:,j))

     ! update sqrtR -- note: redo it for every j, b/o yNaN check below
     sqrtR = 0.0d0
     forall (i=1:Ny) sqrtR(i,i) = sqrt(sigma(i)) 
     DO i = 1, Ny
        if (yNaN(i,j)) then 
           sqrtR(i,:)   = 0.0d0
        end if
     END DO


     ! print *, 'particlefilter for t=', j, ' with ', Nynonan, 'obs' 
     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: APF RESAMPLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     DO k = 1,Nparticles

        ! 2) Fill Particles into state space

        ! update A 

        A(ndxGapRE,ndxGapRE)           = adrift(k)
        A(ndxtrendSI,ndxtrendRE)       = 1 - lambda(k)
        A(ndxtrendSI,ndxtrendSI)       = lambda(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - lambda(k)) * adrift(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = lambda(k) * adrift(k)

        ! update C
        FORALL (i=1:Nsurveys) C(1+i,ndxGapSI,j) = adrift(k) ** i
        DO i = 1, Ny
           if (yNaN(i,j)) then 
              C(i,:,j)     = 0.0d0 
           end if
        END DO


        ! update B
        B(ndxtrendSI,shockndxTrend)    = 1 - lambda(k)
        B(ndxgapSI,shockndxGap)        = 1 - lambda(k)

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
        if (dofAPFfattail .gt. 0) then
           ! use t distribution with dofAPFtail
           ! note: only kernel of MV-T needed
           llf(k)       = -0.5d0 * (logdetSigmaY + (dofAPFfattail + Nynonan) * log(1.0d0 + sum(ytilde ** 2) / dofAPFfattail))
        else 
           llf(k)       = -0.5d0 * (Nynonan * log(2.0d0 * pi) + logdetSigmaY + sum(ytilde ** 2))
        end if


     END DO ! k particles


     if (Nynonan > 0) then

        ! Reweight particles for next round   
        loglikeAPFmax      = maxval(llf)
        llf                = llf - loglikeAPFmax

        APFlike = exp(llf)
        if (doSecondResamplingStep) then
           APFkernelweights     = APFlike / dble(Nparticles) 
           APFkernelsum         = sum(APFkernelweights)
           PARTICLEweights(:,j) = APFkernelweights / APFkernelsum

        else
           APFkernelweights     = APFlike * PARTICLEweights(:,j-1)
           APFkernelsum         = sum(APFkernelweights)
           PARTICLEweights(:,j) = APFkernelweights / APFkernelsum
        end if

        ! where (PARTICLEweights(:,j) < minParticleWeight) PARTICLEweights(:,j) = minParticleWeight
        ! PARTICLEweights(:,j) = PARTICLEweights(:,j) / sum(PARTICLEweights(:,j))

        call drawNDXsysresample(ndx, Nparticles, PARTICLEweights(:,j), Nparticles, uniformdraws(1,j))

        FORALL(k=1:Nparticles) shufflevec(k) = APFlike(ndx(k))
        APFlike = shufflevec

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

        FORALL(k=1:Nparticles) shufflevec(k) = lambda(ndx(k))
        lambda = shufflevec

        FORALL(k=1:Nparticles) shufflevec(k) = adrift(ndx(k))
        adrift = shufflevec

     end if

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: APF STEP
     ! ------------------------------------------------------------------------------------------------------------------------------


     ! 1) Draw Particles 
     errcode     = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, hDELTA, 0.0d0, 1.0d0)
     forall (k=1:Nparticles,i=1:Nsv) hDELTA(i,k) = hDELTA(i,k) * hInno(i) 
     h           = h + hDELTA
     SVol        = exp(h * 0.5d0)


     lambdaPREV  = lambda
     lambdaDELTA = drawdeltatruncnorm1(Nparticles, lambdaPREV, siglambda, 0.0d0, 1.0d0, VSLstream)
     lambda      = lambdaPREV + lambdaDELTA

     aPREV      = adrift
     aDELTA     = drawdeltatruncnorm1(Nparticles, aPREV, siga, -1.0d0, 1.0d0, VSLstream)
     adrift     = aPREV + aDELTA

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     DO k = 1,Nparticles

        xprior      = 0.0d0

        ! 2) Fill Particles into state space

        ! update A 

        A(ndxGapRE,ndxGapRE)           = adrift(k)
        A(ndxtrendSI,ndxtrendRE)       = 1 - lambda(k)
        A(ndxtrendSI,ndxtrendSI)       = lambda(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - lambda(k)) * adrift(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = lambda(k) * adrift(k)

        ! update C
        FORALL (i=1:Nsurveys) C(1+i,ndxGapSI,j) = adrift(k) ** i
        ! zero out missing obs
        DO i = 1, Ny
           if (yNaN(i,j)) then 
              C(i,:,j)     = 0.0d0 
           end if
        END DO


        ! update B
        B(ndxtrendSI,shockndxTrend)    = 1 - lambda(k)
        B(ndxgapSI,shockndxGap)        = 1 - lambda(k)

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

        ! ------------------------------------------------------------------------
        ! DONE: SQRT KALMAN
        ! ------------------------------------------------------------------------


     END DO ! k particles

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------



     ! MDD
     loglikemax        = maxval(llf)
     llf               = llf - loglikemax
     kernelweights     = exp(llf) / APFlike 
     kernelsum         = sum(kernelweights)
     logMDD(j)         = log(kernelsum / Nparticles) + log(APFkernelsum) + loglikemax ! note: loglikeAPFmax correction should be added to APFlike and subtracted from kernelweights (b/o APFlike in denominator), thus cancelling each other

     ! Store NON-reweighted statistics
     FORALL(i=1:Nsv,k=1:Nparticles) DRAWh(i,k,j)   = h(i,k)
     DRAWlambda(:,j)   = lambda
     DRAWa(:,j)        = adrift


     FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,j)  = vecSqrtSigmaX(i,k)
     FORALL(i=1:Nx,k=1:Nparticles) DRAWxhat(i,k,j)             = xposterior(i,k)


     if (Nynonan > 0) then ! nothing to propagate if there was no observed data

        ! particles weights
        PARTICLEweights(:,j) = kernelweights / kernelsum
        ! where (PARTICLEweights(:,j) < minParticleWeight) PARTICLEweights(:,j) = minParticleWeight
        ! PARTICLEweights(:,j) = PARTICLEweights(:,j) / sum(PARTICLEweights(:,j))

        ! resample
        if (doSecondResamplingStep) then

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
           ! SVol needs to be reshuffled, to prep the APF resample step
           DO i=1,Nsv 
              FORALL(k=1:Nparticles) shufflevec(k) = SVol(i,ndx(k))
              SVol(i,:) = shufflevec
           END DO

           FORALL(k=1:Nparticles) shufflevec(k) = lambda(ndx(k))
           lambda = shufflevec

           FORALL(k=1:Nparticles) shufflevec(k) = adrift(ndx(k))
           adrift = shufflevec

        end if ! doSecondResamplingStep
     end if

  END DO ! j=1,T

END SUBROUTINE particlefilter

! @\newpage\subsection{plesmoother}@
SUBROUTINE plesmoother(doInflationNoise, Nsmoother, Nparticles, T, p, SMOOTHERa, SMOOTHERlambda, Nsv, SMOOTHERh, SMOOTHERhinno, SMOOTHERsigma, Nparamparticles, PARAMweights, PARAMsiga, PARAMsiglambda, PARAMhInno, PARAMsigma, a0, a0V, lambda0, lambda0V, Eh0, Vh0, Ex0, sqrtVx0, Ny, y, yNaN, Nx, NsigmaX, Nw, VSLstream)

  ! use embox
  use embox, only : savevec, savemat, int2str
  use gibbsbox, only : drawNDXpdf ! sysresample
  use blaspack, only : qrquery, qrot, ivechU
  use vslbox
  ! use omp_lib
  use timerbox

  IMPLICIT NONE

  integer, intent(in) :: Nparticles
  integer, parameter  :: dofAPFfattail = 0
  logical, parameter  :: doSecondResamplingStep = .false.
  logical, intent(in) :: doInflationNoise
  integer, intent(in) :: Nsmoother, T, p, Nparamparticles, Nsv, Ny, Nw, Nx, NsigmaX
  type (vsl_stream_state), intent(inout) :: VSLstream 

  DOUBLE PRECISION, INTENT(OUT), DIMENSION(Nsmoother,0:T)      :: SMOOTHERlambda
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(Nsmoother,0:T)      :: SMOOTHERa
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(Nsv,Nsmoother,0:T)  :: SMOOTHERh

  DOUBLE PRECISION, INTENT(OUT) :: SMOOTHERsigma(Ny,Nsmoother), SMOOTHERhinno(Nsv,Nsmoother)
  DOUBLE PRECISION :: SMOOTHERsiglambda(Nsmoother), SMOOTHERsiga(Nsmoother) 

  ! note: no intent(out) in previous line since smoothed draws for the params won't be needed outside the plesmoother

  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles)         :: PARAMweights
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles,Nsv)     :: PARAMhinno
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles)         :: PARAMsiga
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles)         :: PARAMsiglambda
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles, Ny)     :: PARAMsigma

  DOUBLE PRECISION, INTENT(IN) :: a0, a0V, lambda0, lambda0V, Eh0(Nsv), Vh0(Nsv), Ex0(Nx), sqrtVx0(Nx,Nx)

  DOUBLE PRECISION, INTENT(IN), DIMENSION(Ny,T) :: y
  LOGICAL, INTENT(IN), DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Ny,T) :: ytilde

  INTEGER :: k,i,j,n,errcode

  ! objects for the forward filters, conditional on parameters
  INTEGER :: filterndx
  DOUBLE PRECISION :: ufilterdraw(Nsmoother)

  DOUBLE PRECISION, DIMENSION(Nparticles,0:T)     :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(Nsv,Nparticles,0:T) :: FILTERh
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T)     :: FILTERlambda
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T)     :: FILTERa
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: FILTERxhat
  DOUBLE PRECISION, DIMENSION(Nsigmax,Nparticles,0:T) :: FILTERsqrtSigmaX
  DOUBLE PRECISION, DIMENSION(T) :: logMDD

  double precision, dimension(Nsv) :: prevSV, prevH
  double precision :: prevLAMBDA
  double precision :: prevA

  ! filter particles and lindsten helpers
  DOUBLE PRECISION, DIMENSION(Nparticles)     :: w, SVpdf, LAMBDApdf, eta, logdetsqrtLAMBDA
  DOUBLE PRECISION, DIMENSION(Nparticles)     :: Apdf
  double precision :: sqrtOmegaHat(Nx,Nx), lambdaHat(Nx), lambda(Nx), sqrtLAMBDA(Nx,Nx), zhat(Nx)
  DOUBLE PRECISION, DIMENSION(Nsv,Nparticles) :: hshock
  DOUBLE PRECISION, DIMENSION(Nparticles)     :: z
  DOUBLE PRECISION, DIMENSION(Nparticles)     :: lambdashock, lambdaCDFlb, lambdaCDFub
  DOUBLE PRECISION, DIMENSION(Nparticles)     :: ashock, aCDFlb, aCDFub

  ! qrR helpers for lindsten et al smoother
  integer :: qr1LWORK, qr2LWORK, qr3LWORK !, qr4LWORK
  double precision :: qr1R(Nx+Ny,Nx), qr2R(Nw+Nx,Nw+Nx), qr3R(Nx+Nx,Nx) !, qr4R(Nx+Ny,Nx)

  ! other lindsten et al objects  
  double precision :: ddot
  double precision :: sqrtSigmaX(Nx,Nx)

  ! smoother objects
  DOUBLE PRECISION :: udraw(0:T-1,Nsmoother)
  INTEGER :: ndx(Nsmoother)
  DOUBLE PRECISION :: cdf

  INTEGER, PARAMETER :: VSLmethodUniform = 0

  ! state space variables
  INTEGER :: Nsurveys
  INTEGER :: ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap
  INTEGER :: ndxgap(2*p)
  DOUBLE PRECISION :: A(Nx,Nx), Bsv(Nx,Nw), C(Ny,Nx,T), invSqrtR(Ny,Ny,T), invR(Ny,Ny,T)

  ! helper for progressbar
  INTEGER :: loopCount
  type(progresstimer) :: timer
  INTEGER :: TID

  ! CHARACTER (LEN=200) :: filename

  TID = 0
  !$ TID = OMP_GET_THREAD_NUM()

  ! qr-workspace queries
  qr1R = 0.0d0
  qr2R = 0.0d0
  qr3R = 0.0d0
  ! qr4R = 0.0d0

  qr1LWORK = qrquery(qr1R)
  qr2LWORK = qrquery(qr2R)
  qr3LWORK = qrquery(qr3R)
  ! qr4LWORK = qrquery(qr4R)

  ! ----------------------------------------------------------
  ! SETUP STATE SPACE 
  ! ----------------------------------------------------------

  Nsurveys   = Ny - 1
  ndxTrendRE = 1
  ndxGapRE   = 2
  ndxTrendSI = 2+p
  ndxGapSI   = ndxTrendSI + 1

  ndxGapREstart = ndxGapRE
  ndxGapREstop  = ndxGapREstart + p - 1
  ndxGapSIstart = ndxGapSI
  ndxGapSIstop  = ndxGapSIstart + p - 1

  ndxgap = (/ ndxgapREstart : ndxGapREstop, ndxgapSIstart : ndxGapSIstop /)

  shockndxTrend = 1
  shockndxGap   = 2

  ! prepare state space
  ! A
  A = 0.0d0
  ! unit root in trend
  A(ndxTrendRE,ndxTrendRE) = 1.0d0
  ! kompanion for gap
  IF (p > 1) THEN
     FORALL(j=1:(p-1)) A(ndxGapRE+j,ndxGapRE-1+j) = 1.0d0
  END IF

  ! B
  ! B      = 0.0d0
  ! B(ndxTrendRE,shockndxTrend) = 1.0d0
  ! B(ndxGapRE,shockndxGap)     = 1.0d0

  ! C
  C         = 0.0d0
  ! inflation
  C(1,ndxTrendRE,:)    = 1.0d0
  C(1,ndxGapRE,:)      = 1.0d0
  ! surveys
  C(2:Ny,ndxTrendSI,:) = 1.0d0

  ! ----------------------------------------------------------
  ! DONE STATE SPACE MATRICES
  ! ----------------------------------------------------------

  ! Draw parameters at time T
  ! t=T
  call drawNDXpdf(ndx, Nsmoother, PARAMweights, Nparamparticles, VSLstream)

  forall (k=1:Nsmoother) SMOOTHERsigma(:,k)      = PARAMsigma(ndx(k),:) 
  forall (k=1:Nsmoother) SMOOTHERsiga(k)         = PARAMsiga(ndx(k)) 
  forall (k=1:Nsmoother) SMOOTHERsiglambda(k)    = PARAMsiglambda(ndx(k)) 
  forall (k=1:Nsmoother) SMOOTHERhinno(:,k)      = PARAMhinno(ndx(k),:) 

  ! draw uniforms: Nsmoother x T
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, Nsmoother * T, udraw, 0.0d0, 1.0d0)
  ! draw uniforms for selecting filtered particles
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, Nsmoother, ufilterdraw, 0.0d0, 1.0d0)

  loopCount = 0
  if (TID .eq. 0) CALL initprogressbar(timer, 15.0d0)

  DO k=1,Nsmoother

     loopCount = loopCount + 1
     ! call fwd filter, Note: logMDD remains unused
     call particlefilter(doSecondResamplingStep, dofAPFfattail, T, Ny, y, yNaN, Nparticles, PARTICLEweights, logMDD, FILTERxhat, FILTERsqrtSigmaX, Nx, NsigmaX, Nw, Ex0, sqrtVx0, p, FILTERa, a0, a0V, SMOOTHERsiga(k), FILTERlambda, lambda0, lambda0V, SMOOTHERsiglambda(k), FILTERh, Nsv, Eh0, Vh0, SMOOTHERhInno(:,k), SMOOTHERsigma(:,k), VSLstream)

     ! draw time T particle 
     ! a) compute filterndx from ufilterdraw
     filterndx = 0
     cdf = 0 
     DO WHILE (cdf < ufilterdraw(k))
        filterndx   = filterndx + 1
        cdf = cdf + PARTICLEweights(filterndx,T)
     END DO
     ! b) assign particles
     SMOOTHERh(:,k,T)    = FILTERh(:,filterndx,T)
     SMOOTHERlambda(k,T) = FILTERlambda(filterndx,T) 
     SMOOTHERa(k,T)      = FILTERa(filterndx,T) 

     ! R
     invSqrtR = 0.0d0
     forall (i=2:Ny,j=1:T) invSqrtR(i,i,j) = 1.0d0 / sqrt(SMOOTHERsigma(i,k))
     invR = 0.0d0
     forall (i=2:Ny,j=1:T) invR(i,i,j) = 1.0d0 / SMOOTHERsigma(i,k)
     if (doInflationNoise) then
        i = 1
        forall (j=1:T) invSqrtR(i,i,j) = 1.0d0 / sqrt(SMOOTHERsigma(i,k))
        forall (j=1:T) invR(i,i,j) = 1.0d0 / SMOOTHERsigma(i,k)
     end if
     ! prepare invR and invSqrtR for missing values
     DO j=1,T
        DO i = 1, Ny
           if (yNaN(i,j)) then 
              invSqrtR(:,i,j) = 0.0d0 ! assumed upper triangular inv(R^(-1/2))'; redundant as long as diagonal
              invR(:,i,j)     = 0.0d0 ! sufficient, assuming upper triangular storage 
              invR(i,:,j)     = 0.0d0 ! a bit redundant since diagonal
           end if
        END DO
     END DO

     ! prepare ytilde
     ytilde = y
     do j=1,T
        call DTRMV('u', 'n', 'n', Ny, invR(:,:,j), Ny, ytilde(:,j), 1) 
     end do


     ! lambdaHat = 0.0d0
     lambda       = 0.0d0
     sqrtOmegaHat = 0.0d0

     prevH      = FILTERh(:,filterndx,T)
     prevSV     = exp(0.5d0 * prevH) 
     prevLAMBDA = FILTERlambda(filterndx,T)
     prevA      = FILTERa(filterndx,T)

     DO j=T-1,0,-1

        ! setup C(t+1)
        FORALL (i=1:Nsurveys) C(1+i,ndxGapSI,j+1) = prevA ** i
        DO i = 1, Ny
           if (yNaN(i,j+1)) C(i,:,j+1) = 0.0d0
        END DO

        ! update lambdahat(t+1)
        lambdaHat = lambda
        call DGEMV('t', Ny, Nx, 1.0d0, C(:,:,j+1), Ny, ytilde(:,j+1), 1, 1.0d0, lambdaHat, 1)

        ! update sqrtOmegaHat(t+1)
        qr1R               = 0.0d0
        qr1R(1:Nx,:)       = sqrtOmegaHat
        qr1R(Nx+1:Nx+Ny,:) = C(:,:,j+1)
        call DTRMM('l', 'u', 't', 'n', Ny, Nx, 1.0d0, invSqrtR(:,:,j+1), Ny, qr1R(Nx+1:Nx+Ny,:), Ny)

        call qrot(qr1R, qr1LWORK)
        sqrtOmegaHat = qr1R(1:Nx,1:Nx)

        ! set up Bsv(t+1)
        Bsv                              = 0.0d0
        Bsv(ndxTrendRE,shockndxTrend)    = prevSV(shockndxTrend)
        Bsv(ndxTrendSI,shockndxTrend)    = (1 - prevLAMBDA) * prevSV(shockndxTrend)
        Bsv(ndxGapRE,shockndxGap)        = prevSV(shockndxGap)
        Bsv(ndxGapSI,shockndxGap)        = (1 - prevLAMBDA) * prevSV(shockndxGap)


        ! set up A(t+1)
        A = 0.0d0
        A(ndxTrendRE,ndxTrendRE)         = 1.0d0
        A(ndxGapRE,ndxGapRE)             = prevA
        A(ndxtrendSI,ndxtrendRE)         = 1 - prevLAMBDA
        A(ndxtrendSI,ndxtrendSI)         = prevLAMBDA
        A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 -  prevLAMBDA) * prevA
        A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) =  prevLAMBDA * prevA

        ! update sqrtOmega(t) and lambda(t)
        qr2R = 0.0d0
        forall (n=1:Nw) qr2R(n,n)   = 1.0d0
        qr2R(Nw+1:Nw+Nx,1:Nw)       = Bsv
        qr2R(Nw+1:Nw+Nx,Nw+1:Nw+Nx) = A

        call dtrmm('l', 'u', 'n', 'n', Nx, Nw+Nx, 1.0d0, sqrtOmegaHat, Nx, qr2R(Nw+1:Nw+Nx,:), Nx)

        call qrot(qr2R, qr2LWORK)
        sqrtOmegaHat = qr2R(Nw+1:Nw+Nx,Nw+1:Nw+Nx)



        ! lambda (overwriting Bsv, A)
        call dtrsm('r', 'u', 'n', 'n', Nx, Nw, 1.0d0, qr2R(1:Nw,1:Nw), Nw, Bsv, Nx)
        call dgemm('n', 'n', Nx, Nx, Nw, -1.0d0, Bsv, Nx, qr2R(1:Nw,Nw+1:Nw+Nx), Nw, 1.0d0, A, Nx)
        call dgemv('t', Nx, Nx, 1.0d0, A, Nx, lambdaHat, 1, 0.0d0, lambda,1)

        ! draw smoothed particle k
        ! a) compute kernel weights

        ! - backward smoother PDF
        eta              = 0.0d0
        logdetsqrtLAMBDA = 0.0d0
        sqrtLAMBDA       = 0.0d0
        do n=1,Nparticles

           sqrtSigmaX = ivechU(FILTERsqrtSigmaX(:,n,j), Nx)

           ! sqrtLAMBDA
           qr3R = 0.0d0
           forall(i=1:Nx) qr3R(Nx+i,i) = 1.0d0
           qr3R(1:Nx,1:Nx) = sqrtOmegaHat 
           call dtrmm('r', 'u', 't', 'n', Nx, Nx, 1.0d0, sqrtSigmaX, Nx, qr3R(1:Nx,1:Nx), Nx)

           call qrot(qr3R, qr3LWORK)
           sqrtLAMBDA = qr3R(1:Nx,1:Nx) ! save yourself sqrtLAMBDA


           logdetsqrtLAMBDA(n) = log(abs(sqrtLAMBDA(1,1)))
           do i = 2,Nx
              logdetsqrtLAMBDA(n) = logdetsqrtLAMBDA(n) + log(abs(sqrtLAMBDA(i,i)))
           end do

           ! eta
           zhat             = FILTERXhat(:,n,j)

           eta(n)           = ddot(Nx, zhat, 1, lambda, 1)

           call dtrmv('u', 'n', 'n', Nx, sqrtOmegaHat, Nx, zhat, 1)
           eta(n)           = sum(zhat ** 2) - 2 * eta(n)

           call dtrmv('u', 't', 'n', Nx, sqrtOmegaHat, Nx, zhat, 1)
           zhat = lambda - zhat

           call dtrmv('u', 'n', 'n', Nx, sqrtSigmaX, Nx, zhat, 1) 

           call dtrsv('u', 't', 'n', Nx, sqrtLAMBDA, Nx, zhat, 1) 

           eta(n) = eta(n) - sum(zhat ** 2)

        end do

        ! - fix big etas (note: pdf's get rescaled anyway)
        eta = eta - maxval(eta)

        ! - pdf of SV move
        FORALL (i=1:Nsv,n=1:Nparticles) hshock(i,n) = (prevH(i) - FILTERh(i,n,j)) / SMOOTHERhInno(i,k) 
        SVpdf = sum(hshock ** 2, 1) 

        ! - pdf of LAMBDA move
        z = (1.0d0 - FILTERlambda(:,j)) / SMOOTHERsiglambda(k)
        call vdcdfnorm(Nparticles, z, lambdaCDFub)
        z = (0.0d0 - FILTERlambda(:,j)) / SMOOTHERsiglambda(k)
        call vdcdfnorm(Nparticles, z, lambdaCDFlb)
        lambdashock = (prevLAMBDA - FILTERlambda(:,j)) / SMOOTHERsiglambda(k)
        LAMBDApdf   = lambdashock ** 2 ! LAMBDApdf   = exp(-.5d0 * (lambdashock ** 2)) / (lambdaCDFub - lambdaCDFlb)


        ! - pdf of A move
        z = (1.0d0 - FILTERa(:,j)) / SMOOTHERsiga(k)
        call vdcdfnorm(Nparticles, z, aCDFub)
        z = (-1.0d0 - FILTERa(:,j)) / SMOOTHERsiga(k)
        call vdcdfnorm(Nparticles, z, aCDFlb)
        ashock      = (prevA - FILTERa(:,j)) / SMOOTHERsiga(k)
        Apdf        = ashock ** 2



	! try computing w in logs as much as possible
        w = -0.5d0 * (SVpdf + LAMBDApdf  + Apdf + eta) - logdetsqrtLAMBDA - log(lambdaCDFub - lambdaCDFlb) - log(aCDFub - aCDFlb) + log(PARTICLEweights(:,j))

        w = w - maxval(w)
        w = exp(w)
        w = w / sum(w)

        ! b) draw using kernel weights
        i   = 0
        cdf = 0.0d0
        DO WHILE (cdf < udraw(j,k))
           i   = i + 1
           cdf = cdf + w(i)
        END DO
        SMOOTHERh(:,k,j)    = FILTERh(:,i,j)
        prevH               = FILTERh(:,i,j)
        prevSV              = exp(0.5d0 * prevH)

        SMOOTHERlambda(k,j) = FILTERlambda(i,j)
        prevLAMBDA          = FILTERlambda(i,j)

        SMOOTHERa(k,j)      = FILTERa(i,j)
        prevA               = FILTERa(i,j)


     END DO ! j=T-1,0,-1

     if (TID .eq. 0) CALL progressbarcomment(dble(loopCount) / dble(Nsmoother), timer, 'Particle Smoother Step')

  END DO ! k

  ! convert log-variances back into vols
  SMOOTHERh = exp(SMOOTHERh * 0.5d0)

END SUBROUTINE plesmoother

! @\newpage\subsection{particleSmootherXhat}@
SUBROUTINE particleSmootherXhat(T, Ny, y, yNaN, Nsv, Nsmoother, SMOOTHERa, p, SMOOTHERlambda, SMOOTHERsvol, SMOOTHERxhat, Nx, Nw, Ex0, sqrtVx00, SMOOTHERsigma)

  ! use embox
  ! use gibbsbox
  use statespacebox, only : DLYAP, disturbancesmootherA3B3C3nanscalar ! samplerA3B3C3noisenan

  use vslbox
  ! use timerbox
  ! use omp_lib

  IMPLICIT NONE

  INTENT(INOUT) :: SMOOTHERxhat
  INTENT(IN)    :: T, Ny, Ex0, sqrtVx00, y, Nsv, Nsmoother, SMOOTHERlambda, SMOOTHERa, SMOOTHERsvol, Nx, Nw, p, SMOOTHERsigma

  INTEGER :: T, Ny, Nsv, Nsmoother, Nx, Nw, Nsurveys, p

  INTEGER :: j,k,i

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Nsv,Nsmoother,0:T) :: SMOOTHERsvol
  DOUBLE PRECISION, DIMENSION(Nsmoother,0:T) :: SMOOTHERlambda
  DOUBLE PRECISION, DIMENSION(Nsmoother,0:T) :: SMOOTHERa
  DOUBLE PRECISION, DIMENSION(Nx,0:T,Nsmoother) :: SMOOTHERxhat

  ! particles
  DOUBLE PRECISION :: x(Nx,0:T), xshock(Nx,T), ynoise(Ny,T)

  ! state space objects
  DOUBLE PRECISION :: SMOOTHERsigma(Ny,Nsmoother)
  DOUBLE PRECISION :: Ex0(Nx), sqrtVx00(Nx,Nx), sqrtVx0(Nx,Nx), A(Nx,Nx,T), B(Nx,Nw,T), C(Ny,Nx,T), noisevol(Ny,T)
  DOUBLE PRECISION :: Vx0(Nx,Nx), dummy(Nx,Nx)
  DOUBLE PRECISION :: ygap0variance(2*p,2*p), gaptransition(2*p,2*p), gapshock0loadings(2*p,Nw)

  INTEGER :: errcode

  ! index variables for state space
  INTEGER ::   ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap

  ! CHARACTER (LEN=100) :: filename

  Nsurveys = Ny - 1

  ! state-vector indices
  ndxTrendRE = 1
  ndxGapRE   = 2
  ndxTrendSI = 2+p
  ndxGapSI   = ndxTrendSI + 1

  ndxGapREstart = ndxGapRE
  ndxGapREstop  = ndxGapREstart + p - 1
  ndxGapSIstart = ndxGapSI
  ndxGapSIstop  = ndxGapSIstart + p - 1

  shockndxTrend = 1
  shockndxGap   = 2

  ! prepare state space
  ! A
  A = 0.0d0
  ! unit root in trend
  A(ndxTrendRE,ndxTrendRE,:) = 1.0d0
  ! kompanion for gap
  IF (p > 1) THEN
     FORALL(j=1:(p-1)) A(ndxGapRE+j,ndxGapRE-1+j,:) = 1.0d0
  END IF

  ! B
  B                             = 0.0d0
  B(ndxTrendRE,shockndxTrend,:) = 1.0d0
  B(ndxGapRE,shockndxGap,:)     = 1.0d0

  ! C
  C         = 0.0d0
  ! inflation
  C(1,ndxTrendRE,:)    = 1.0d0
  C(1,ndxGapRE,:)      = 1.0d0
  ! surveys
  C(2:Ny,ndxTrendSI,:) = 1.0d0

  ! prepare C and noisevol for missing values
  ! DO k=1,T
  !    DO j = 1, Ny
  !       if (yNaN(j,k)) then
  !          C(j,:,k) = 0.0d0
  !          ! noisevol(j,k) = 0.0d0
  !       end if
  !       if (yNaN(j,k) .AND. y(j,k) /= 0.0d0 ) then
  !          write (*,*) 'YNAN PATTERN DOES NOT MATCH ZEROS IN Y'
  !       end if
  !    END DO
  ! END DO

  sqrtVx0 = sqrtVx00 


  DO k=1,Nsmoother


     ! update A 
     forall (j=1:T) A(ndxGapRE,ndxGapRE,j)      = SMOOTHERa(k,j)
     forall (j=1:T) A(ndxtrendSI,ndxtrendRE,j)  = 1 - SMOOTHERlambda(k,j)
     forall (j=1:T) A(ndxtrendSI,ndxtrendSI,j)  = SMOOTHERlambda(k,j)
     forall (j=1:T) A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop,j) = (1 - SMOOTHERlambda(k,j)) * SMOOTHERa(k,j)
     forall (j=1:T) A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop,j) = SMOOTHERlambda(k,j)* SMOOTHERa(k,j)

     ! update noisevol
     forall(j=1:T) noisevol(:,j) = sqrt(SMOOTHERsigma(:,k)) 


     ! update C
     FORALL (i=1:Nsurveys, j=1:T) C(1+i,ndxGapSI,j) = SMOOTHERa(k,j) ** i
     DO j=1,T
        DO i = 1, Ny
           if (yNaN(i,j)) then
              noisevol(i,j) = 0.0d0
              C(i,:,j) = 0.0d0
           end if
        END DO
     END DO


     ! update B
     forall (j=1:T) B(ndxtrendRE,shockndxTrend,j)    = SMOOTHERsvol(shockndxTrend,k,j)
     forall (j=1:T) B(ndxgapRE,shockndxGap,j)        = SMOOTHERsvol(shockndxGap,k,j)
     forall (j=1:T) B(ndxtrendSI,shockndxTrend,j)    = (1 - SMOOTHERlambda(k,j)) * SMOOTHERsvol(shockndxTrend,k,j)
     forall (j=1:T) B(ndxgapSI,shockndxGap,j)        = (1 - SMOOTHERlambda(k,j)) * SMOOTHERsvol(shockndxGap,k,j)


     ! prior variance for linear states
     gaptransition                  = 0.0d0
     gaptransition(1:p,1:p)         = SMOOTHERa(k,0)
     gaptransition(p+1:2*p,1:p)     = (1 - SMOOTHERlambda(k,0)) * SMOOTHERa(k,0)
     gaptransition(p+1:2*p,p+1:2*p) = SMOOTHERlambda(k,0) * SMOOTHERa(k,0)

     gapshock0loadings                   = 0.0d0
     gapshock0loadings(1,shockndxGap)    = SMOOTHERsVOL(shockndxGap,k,0)
     gapshock0loadings(p+1,shockndxGap)  = (1 - SMOOTHERlambda(k,0)) * SMOOTHERsVOL(shockndxGap,k,0)

     CALL DLYAP(ygap0variance, gaptransition, gapshock0loadings, 2 * p, Nw, errcode) 
     if (errcode /= 0) then
        write (*,*) 'DLYAP error (ygap0variance)', errcode
        stop 1
     end if

     CALL DPOTRF('L', 2 * p, ygap0variance, 2 * p, errcode) 
     if (errcode /= 0) then
        write (*,*) 'DPOTRF error (ygap0variance)', errcode
        stop 1
     end if
     ! zero out the upper triangular
     ! FORALL (i=1:2*p-1) ygap0variance(i,i+1:2*p) = 0.0d0
     sqrtVx0(ndxgapREstart:ndxgapREstop,ndxgapREstart:ndxgapREstop) = ygap0variance(1:p,1:p)
     sqrtVx0(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = ygap0variance(p+1:2*p,1:p)
     sqrtVx0(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = ygap0variance(p+1:2*p,p+1:2*p)

     Vx0 = 0.0d0
     call DSYRK('U','N',Nx,Nx,1.0d0,sqrtVx0,Nx,0.0d0,Vx0,Nx)


     ! 3) KALMAN SMOOTHER:
     ! note: xshock and ynoise will not be used
     ! call samplerA3B3C3noisenan(x,xshock,ynoise,y,yNaN,T,Ny,Nx,Nw,A,B,C,noisevol,Ex0,sqrtVx0,VSLstream,errcode)
     call disturbancesmootherA3B3C3nanscalar(x,xshock,dummy,ynoise,y,yNaN,T,Ny,Nx,Nw,A,B,C,noisevol,Ex0,Vx0,errcode)

     FORALL (i=1:Nx) SMOOTHERxhat(i,:,k)  = x(i,:)

  END DO ! k


END SUBROUTINE particleSmootherXhat


