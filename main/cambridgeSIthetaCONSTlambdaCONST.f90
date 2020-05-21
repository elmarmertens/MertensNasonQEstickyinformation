PROGRAM main

  ! constant-parameter AR1 and constant lambda

  USE embox, only : hrulefill, loadmatrix, loadarray1, savemat, savematlogical, savevec, storeEstimates, storeEstimatesTranspose, loft, timestampstr, es30d16, int2str
  USE blaspack, only : vech, ivech, eye
  USE gibbsbox, only : drawNDXpdf, drawNDXsysresample

  USE vslbox
  USE timerbox
  USE omp_lib

  IMPLICIT NONE

  ! ----------------------------------------------------------------------------------

  CHARACTER (LEN=200) :: modellabel ='SIthetaCONSTlambdaCONST'

  logical, parameter :: doTimestamp = .false., doGains = .true.
  logical :: doSmoother = .true.
  logical :: doInflationNoise = .true.
  logical :: doLinearUncertainty = .false.


  logical :: doSecondResamplingStep = .false.
  integer :: dofAPFfattail = 0  ! set to zero to use normal APF
  
  double precision :: priornoisevol = sqrt(0.1d0)


  INTEGER, PARAMETER :: dof0 = 3

  INTEGER, PARAMETER :: p = 1, Nsurveys = 5, Ny = Nsurveys + 1, Nx = 2 * (1 + p), NsigmaX = Nx * (Nx + 1) / 2, Nw = 2, Nsv = 2
  INTEGER, PARAMETER :: Nxx = Nx + 1, NsigmaXX = Nxx * (Nxx + 1) / 2 ! for state space with lagged gap
  INTEGER, PARAMETER :: ndxtrendRE = 1 , ndxgapRE = 2 , ndxtrendSI = 1 + p + 1 , ndxgapSI = ndxtrendSI + 1

  INTEGER :: Nparticles, Nsmoother, NsmootherX, Nmixturedraws
  INTEGER :: smootherNparticles
  INTEGER :: T,i,j,k,status 

  ! filter particles
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)     :: PARTICLEweights, DRAWlambda
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: logMDD
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: DRAWsvol, DRAWxhat, DRAWxsig, DRAWsqrtSigmaX
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: DRAWxgain

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: DRAWgapnoisehat ! contains RE and SI

  ! smoothed parameters
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: PARAMweights
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)     :: PARAMhInno, PARAMsigma
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: PARAMlambda
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: PARAMa
  ! uncertainty measures
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: Xsig 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: STATEsvol
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: STATElambda
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: STATEa
  DOUBLE PRECISION, DIMENSION(Ny) :: STATEsigmasqrt
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: horizons, measurementErrorVol
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: SIndx

  ! inflation forecast draws
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: DRAWpihat, DRAWpifcst ! Nsurveys, Nparticles, 0:T

  ! Forecast Persistence
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: FcstPersistence ! Nsurveys, 0:T
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWpipersistence ! Nsurveys, Nparticles, 0:T


  ! smoother particles
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: SMOOTHERlambda
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: SMOOTHERsvol
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: SMOOTHERx

  DOUBLE PRECISION, DIMENSION(Nx,Nx) :: sqrtVx0
  DOUBLE PRECISION, DIMENSION(Nx)    :: Ex0
  DOUBLE PRECISION, DIMENSION(Nsv)   :: SVar0, Eh0, Vh0

  DOUBLE PRECISION :: lambdaAlpha, lambdaBeta
  DOUBLE PRECISION :: lambdaN = 1.0d0 ! can but need not be an integer

  DOUBLE PRECISION :: a0, aV0 

  ! DRAWS of various scale parameters

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRAWsigma, DRAWhInno
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: SMOOTHERhinno ! legacy code
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: SMOOTHERsigma ! needed for particlesmootherX

  ! priors for scale parameters
  DOUBLE PRECISION ::  hvarT(Nsv), sigmaT(Ny)
  INTEGER :: hvarDof(Nsv), sigmaDof(Ny)

  ! DRIFT AR 1 (specialized to p=1)
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DRAWa
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SMOOTHERa

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: y
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yNaN

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Xdraws

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: theta1
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: theta2
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ndx

  TYPE(progresstimer) :: timer
  CHARACTER (LEN=200) :: filename, datafile, nandatafile, fileXT, datalabel, this

  ! VSL Random Stuff
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

  ! thorough
  Nparticles    = 10 ** 5
  Nmixturedraws = 10 ** 4
  Nsmoother     = 10 ** 4
  NsmootherX    = 100


  ! quick
  Nparticles    = 10 ** 4
  Nmixturedraws = 10 ** 3
  Nsmoother     = 10 ** 4
  NsmootherX    = 100


  ! quicker
  Nparticles    = 10 ** 2
  Nmixturedraws = 10 ** 3
  Nsmoother     = 10 ** 2
  NsmootherX    = 100


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
  ! brng    = vsl_brng_mt19937
  seed    = 0
  errcode = vslnewstream(VSLstream, vsl_brng_mt2203, seed)  

  WRITE(*,'(a25, i20, i20)') 'LAUNCHING VSLSTREAM ', VSLstream%descriptor1, VSLstream%descriptor2
  print *, 'vsl_brng', vsl_brng_mt2203

  datalabel        = 'cambridge2018GDPD'
  call getarguments(doInflationNoise, Nparticles, doSmoother, Nsmoother, smootherNparticles, datalabel, doSecondResamplingStep, dofAPFfattail, lambdaN)

  ! runtime parameters :end: 

  call hrulefill
  if (doSmoother) then
     print *, 'Particle Filter (and Smoother) estimation of ' // modellabel
  else
     print *, 'Particle Filter estimation of ' // modellabel
  end if

  ! CONSTRUCT FILE EXTENTSIONS
  fileXT = '.particles.' // trim(datalabel) // '.' // trim(modellabel)
  if (.not. doInflationNoise)  fileXT = trim(filext) // '.nonoise'
  if (doSecondResamplingStep)  fileXT = trim(filext) // '.2ndAPFresample'
  if (lambdaN > 1.0d0)  fileXT = trim(filext) // '.lambdaN' // trim(int2str(int(lambdaN)))
  if (dofAPFfattail .gt. 0)    fileXT = trim(filext) // '.dofAPFfattail' // trim(int2str(dofAPFFattail))
  fileXT = trim(filext) // '.dat'
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
  print *, 'data=       ' // datalabel
  print *, 'model=      ' // modellabel
  print *, 'Ny=         ', Ny
  print *, 'T=          ', T
  print *, 'Nparticles= ', Nparticles
  if (doSmoother) print *, 'Nsmoother=          ', Nsmoother
  if (doSmoother) print *, 'smootherNparticles= ', smootherNparticles
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
  CALL HRULEFILL

  ! trivial since Nsv = 1 here
  Svar0    = (/ 0.6 / 3.0d0, 0.6 * 2.0d0 / 3.0d0 /)
  Vh0      = 10.0d0
  Eh0      = log(Svar0) - Vh0 * 0.5d0

  hvarDof = dof0
  hvarT   = (0.2d0 ** 2) * (dble(hvarDof) - 2.0d0)

  ! lambda
  lambdaAlpha  = 1.0d0
  lambdaBeta   = 1.0d0

  ! a
  a0        = 0.0d0
  aV0       = 1.0d0

  ! Linear prior
  Ex0       = 0.0d0
  Ex0(1)    = 2.0d0
  Ex0(2+p)  = 2.0d0

  ! sqrtVx0, expressed as lower  triangular-choleski factor (sqrtVx0 * sqrtVx0')
  call eye(sqrtVx0, 10.0d0) ! non-zero prior gap-variance matters only for initial conditions for sampling lagged gap as used for estimating AR(1) coefficient 
  sqrtVx0(1,1) = 100.0d0 ! sqrt(2.0d0)
  sqrtVx0(2+p,2+p) = 100.0d0 ! sqrt(2.0d0) 

  sigmaDof    = 20 ! dof0
  sigmaT      = (priornoisevol ** 2) * (dble(sigmaDof) - 2.0d0)


  ! allocate memory for draws
  ALLOCATE (PARTICLEweights(Nparticles,0:T),logMDD(T), DRAWlambda(Nparticles,T), DRAWxhat(Nx,Nparticles,0:T), DRAWxsig(Nx,Nparticles,0:T), DRAWsvol(Nparticles,Nsv,0:T), DRAWxgain(Nx,Ny,Nparticles,T), DRAWsqrtSigmaX(NsigmaX,Nparticles,0:T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (draws)'
  END IF

  ALLOCATE (DRAWgapnoisehat(Nparticles,0:T,2), STAT=status)

  ALLOCATE (DRAWa(Nparticles,0:T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (DRAWa)'
  END IF

  ALLOCATE (DRAWpihat(Nsurveys,Nparticles,0:T), DRAWpifcst(Nsurveys,Nparticles,0:T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (DRAWpihat)'
  END IF

  ALLOCATE (DRAWpipersistence(Nsurveys,Nparticles,1:T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (DRAWpihat)'
  END IF

  ! scale parameters
  ALLOCATE (DRAWsigma(Nparticles,Ny,T), DRAWhInno(Nparticles,Nsv,T), STAT=status)
  IF (status /= 0) THEN
     WRITE (*,*) 'Allocation problem (scale parameter draws)'
  END IF

  PARTICLEweights = 1.0d0 / dble(Nparticles)
  logMDD          = 0.0d0
  DRAWlambda      = 0.0d0
  DRAWa           = 0.0d0
  DRAWxhat        = 0.0d0
  DRAWgapnoisehat = 0.0d0
  DRAWxsig        = 0.0d0
  DRAWsvol        = 0.0d0
  DRAWxgain       = 0.0d0

  DRAWpihat       = 0.0d0
  DRAWpifcst      = 0.0d0
  DRAWpipersistence = 0.0d0

  DRAWhInno     = 0.0d0
  DRAWsigma     = 0.0d0 

  CALL plefilter(doSecondResamplingStep, dofAPFfattail, doInflationNoise, T, Ny, y, yNaN, Nparticles, PARTICLEweights, logMDD, DRAWpihat, DRAWpifcst, DRAWpipersistence, DRAWxhat, DRAWsqrtSigmaX, DRAWxsig, DRAWxgain, Nxx, NsigmaXX, Nx, NsigmaX, Nw, Ex0, sqrtVx0, DRAWgapnoisehat, p, DRAWa, a0, aV0, DRAWlambda, lambdaAlpha, lambdaBeta, lambdaN, DRAWsvol, Nsv, Eh0, Vh0, DRAWhInno, hvarT, hvarDof, DRAWsigma, sigmaT, sigmaDof, VSLstream,timer)


  CALL HRULEFILL
  WRITE (*,*) 'PARTICLE FILTER IS DONE!'
  CALL HRULEFILL

  ! WRITE SETTINGS
  CALL HRULEFILL
  filename = 'settings' // trim(adjustl(filext))
  OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
  WRITE(4,'(a20,a40)') 'TIME: ', timestampstr()
  WRITE(4,'(a20,a40)') 'Data: ', datalabel
  WRITE(4,'(a20,a40)') 'Model: ', modellabel
  WRITE(4,'(a60)') repeat('-',60)
  WRITE(4,'(a20,I40)') 'Nparticles: ', Nparticles
  if (doSmoother) THEN
     WRITE(4,'(a20,I40)') 'Nsmoother:', Nsmoother
     WRITE(4,'(a20,I40)') 'smootherNparticles:', smootherNparticles
  END IF
  WRITE(4,'(a20,I40)') 'p: ', p
  WRITE(4,'(a20,10F8.2)') 'lambdaN0: ', lambdaN
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
  filename = 'LOGMDD' // filext
  call savevec(logMDD, filename)
  call hrulefill
  WRITE (*,*) 'STORED MDD'
  WRITE (*,*) '... the log-MDD is ', sum(logMDD)
  WRITE (*,*) '... w/average contribution ', sum(logMDD) / T
  call hrulefill

  ALLOCATE (theta1(T), STAT=status)

  ! ESS
  filename  = 'ESS' // filext
  theta1    = 1 / sum(PARTICLEweights(:,1:T) ** 2, 1) / Nparticles 
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ESS'

  ! PARTICLEWEIGHTS
  filename  = 'PARTICLEWEIGHTS' // filext
  call savemat(PARTICLEWEIGHTS(:,1:T), filename)
  WRITE (*,*) 'STORED PARTICLEWEIGHTS'


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

  THIS      = 'TAUHATSI' 
  filename  = trim(this) // filext
  theta1 = sum(PARTICLEweights(:,1:T) * DRAWxhat(3,:,1:T), 1)
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ' // trim(this)

  THIS      = 'GAPHATSI' 
  filename  = trim(this) // filext
  theta1 = sum(PARTICLEweights(:,1:T) * DRAWxhat(4,:,1:T), 1)
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ' // trim(this)
  DEALLOCATE (theta1)

  ! store GAPNOISE
  THIS      = 'GAPNOISEHATRE' 
  filename  = trim(this) // filext
  theta1 = sum(PARTICLEweights(:,1:T) * DRAWgapnoisehat(:,1:T,1), 1)
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ' // trim(this)
  THIS      = 'GAPNOISEHATSI' 
  filename  = trim(this) // filext
  theta1 = sum(PARTICLEweights(:,1:T) * DRAWgapnoisehat(:,1:T,2), 1)
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED ' // trim(this)

  ! store inflation forecasts
  ! RE
  ALLOCATE (theta2(T,Nsurveys))
  THIS      = 'PIHAT' 
  filename  = trim(this) // filext
  forall (j=1:T,i=1:Nsurveys) theta2(j,i) = sum(PARTICLEweights(:,j) * DRAWpihat(i,:,j))
  call savemat(theta2, filename)
  WRITE (*,*) 'STORED ' // trim(this)

  ! SI
  THIS      = 'PIFCST' 
  filename  = trim(this) // filext
  forall (j=1:T,i=1:Nsurveys) theta2(j,i) = sum(PARTICLEweights(:,j) * DRAWpifcst(i,:,j))
  call savemat(theta2, filename)
  WRITE (*,*) 'STORED ' // trim(this)

  DEALLOCATE (DRAWpihat, DRAWpifcst)

  THIS      = 'PIPERSISTENCE' 
  filename  = trim(this) // filext
  forall (j=1:T,i=1:Nsurveys) theta2(j,i) = sum(PARTICLEweights(:,j) * DRAWpipersistence(i,:,j))
  call savemat(theta2, filename)
  WRITE (*,*) 'STORED ' // trim(this)
  DEALLOCATE(DRAWpipersistence)

  DEALLOCATE (theta2)

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

  THIS      = 'TAUSI' 
  filename  = trim(this) // filext
  CALL storeEstimatesTranspose(Xdraws(ndxtrendSI,:,:),T,Nmixturedraws,filename)
  WRITE (*,*) 'STORED ' // trim(this)

  THIS      = 'GAPSI' 
  filename  = trim(this) // filext
  CALL storeEstimatesTranspose(Xdraws(ndxgapSI,:,:),T,Nmixturedraws,filename)
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


  ! LAMBDA
  filename  = 'LAMBDAHAT' // filext
  theta1 = sum(PARTICLEweights(:,1:T) * DRAWlambda(:,1:T), 1)
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED LAMBDAHAT'
  ! draw distribution
  filename  = 'LAMBDA' // filext
  FORALL (j=1:T,k=1:Nmixturedraws) theta2(k,j) = DRAWlambda(ndx(k,j),j) 
  CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
  WRITE (*,*) 'STORED LAMBDA'

  ! change in lambda over time
  ! filename  = 'DELTALAMBDA' // filext
  ! FORALL (j=1:T,k=1:Nmixturedraws) theta2(k,j) = DRAWlambda(ndx(k,j),j) - DRAWlambda(ndx(k,j),1) 
  ! CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
  ! WRITE (*,*) 'STORED DELTALAMBDA'

  ! a
  filename  = 'AHAT'  // filext
  theta1 = sum(PARTICLEweights(:,1:T) * DRAWa(:,1:T), 1)
  call savevec(theta1, filename)
  WRITE (*,*) 'STORED AHAT'
  ! draw distribution
  filename  = 'A'  // filext
  FORALL (j=1:T,k=1:Nmixturedraws) theta2(k,j) = DRAWa(ndx(k,j),j) 
  CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
  WRITE (*,*) 'STORED A'

  ! SCALE PARAMETERS

  ! hInno
  DO i=1,Nsv
     filename  = 'HINNOHAT'  // trim(int2str(i))  // filext
     theta1 = sum(PARTICLEweights(:,1:T) * (DRAWhInno(:,i,1:T) ** 2), 1)
     call savevec(theta1, filename)
     WRITE (*,*) 'STORED HINNOHAT'  // trim(int2str(i))  
     ! draw distribution
     filename  = 'HINNO'  // trim(int2str(i))   // filext
     FORALL (j=1:T,k=1:Nmixturedraws) theta2(k,j) = DRAWhInno(ndx(k,j),i,j) ** 2
     CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED HINNO'  // trim(int2str(i))  
  END DO

  ! Sigma
  DO i=1,Ny
     filename  = 'SIGMAHAT'  // trim(int2str(i))  // filext
     theta1 = sum(PARTICLEweights(:,1:T) * DRAWsigma(:,i,1:T), 1)
     call savevec(theta1, filename)
     WRITE (*,*) 'STORED SIGMAHAT'  // trim(int2str(i))  
     ! draw distribution
     filename  = 'SIGMA'  // trim(int2str(i))   // filext
     FORALL (j=1:T,k=1:Nmixturedraws) theta2(k,j) = DRAWsigma(ndx(k,j),i,j) 
     CALL storeEstimatesTranspose(theta2,T,Nmixturedraws,filename)
     WRITE (*,*) 'STORED SIGMA'  // trim(int2str(i))  
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
  WRITE (*,*) 'MDD:'
  WRITE (*,*) '... the log-MDD is ', sum(logMDD)
  WRITE (*,*) '... w/average contribution ', sum(logMDD) / T
  call hrulefill
  DEALLOCATE (logMDD)

  DEALLOCATE (DRAWsvol)
  DEALLOCATE (DRAWxhat, DRAWsqrtSigmaX)
  DEALLOCATE (DRAWgapnoisehat)

  ALLOCATE (PARAMweights(Nparticles))
  ALLOCATE (PARAMhInno(Nparticles,Nsv), PARAMsigma(Nparticles,Ny))
  ALLOCATE (PARAMlambda(Nparticles))
  ALLOCATE (PARAMa(Nparticles))

  PARAMweights =  PARTICLEweights(:,T)
  PARAMhInno   =  DRAWhInno(:,:,T)
  PARAMsigma   =  DRAWsigma(:,:,T)
  PARAMlambda  =  DRAWlambda(:,T)
  PARAMa       = DRAWa(:,T)

  DEALLOCATE (PARTICLEweights)
  DEALLOCATE (DRAWhInno,DRAWsigma)
  DEALLOCATE (DRAWlambda)
  DEALLOCATE (DRAWa)

  IF (doSmoother) THEN

     ! ! ----------------------------------------------------------------------------
     ! ! SMOOTHER
     ! ! ----------------------------------------------------------------------------
     ALLOCATE (SMOOTHERsvol(Nsv,Nsmoother,0:T), SMOOTHERlambda(Nsmoother), SMOOTHERa(Nsmoother))
     ALLOCATE (SMOOTHERsigma(Ny,Nsmoother))
     ALLOCATE (SMOOTHERhinno(Nsv,Nsmoother))

     call hrulefill
     print *, 'STARTING w/SMOOTHER ...'
     CALL plesmoother(doInflationNoise, Nsmoother, smootherNparticles, T, p, SMOOTHERa, SMOOTHERlambda, Nsv, SMOOTHERsvol, SMOOTHERhinno, SMOOTHERsigma, Nparticles, PARAMweights, PARAMa, PARAMlambda, PARAMhInno, PARAMsigma, Eh0, Vh0, Ex0, sqrtVx0, Ny, y, yNaN, Nx, NsigmaX, Nw, VSLstream)

     DEALLOCATE (PARAMweights)
     DEALLOCATE (PARAMsigma, PARAMhInno)     
     DEALLOCATE (PARAMlambda)
     DEALLOCATE (PARAMa)



     ! ----------------------------------------------------------------------------------------------------------------------------------------------------------
     ! X SMOOTHER
     ! ----------------------------------------------------------------------------------------------------------------------------------------------------------
     print *, '... Smoother X ...'
     ALLOCATE (SMOOTHERx(Nx,0:T,NsmootherX,Nsmoother))
     SMOOTHERx = 0.0d0
     call particleSmootherX(T, Ny, y, yNaN, Nsv, Nsmoother, SMOOTHERa, p, SMOOTHERlambda, SMOOTHERsvol, NsmootherX, SMOOTHERx, Nx, Nw, Ex0, sqrtVx0, SMOOTHERsigma)

     ! STORE SMOOTHER X
     ALLOCATE (theta2(T,Nsmoother * NsmootherX))

     ! store smoothed trends and gaps
     THIS      = 'smootherTAURE' 
     filename  = trim(this) // filext
     FORALL (j=1:NsmootherX,k=1:Nsmoother) theta2(:,(k-1) * NsmootherX + j) = SMOOTHERx(ndxtrendRE,1:T,j,k)
     CALL storeEstimates(theta2,T,Nsmoother * NsmootherX,filename)
     WRITE (*,*) 'STORED ' // trim(this)

     THIS      = 'smootherTAUREdelta' 
     filename  = trim(this) // filext
     FORALL (j=1:NsmootherX,k=1:Nsmoother,i=1:T) theta2(i,(k-1) * NsmootherX + j) = SMOOTHERx(ndxtrendRE,i,j,k) - SMOOTHERx(ndxtrendRE,1,j,k)
     CALL storeEstimates(theta2,T,Nsmoother * NsmootherX,filename)
     WRITE (*,*) 'STORED ' // trim(this)

     THIS      = 'smootherGAPRE' 
     filename  = trim(this) // filext
     FORALL (j=1:NsmootherX,k=1:Nsmoother) theta2(:,(k-1) * NsmootherX + j) = SMOOTHERx(ndxgapRE,1:T,j,k)
     CALL storeEstimates(theta2,T,Nsmoother * NsmootherX,filename)
     WRITE (*,*) 'STORED ' // trim(this)

     THIS      = 'smootherTAUSI' 
     filename  = trim(this) // filext
     FORALL (j=1:NsmootherX,k=1:Nsmoother) theta2(:,(k-1) * NsmootherX + j) = SMOOTHERx(ndxtrendSI,1:T,j,k)
     CALL storeEstimates(theta2,T,Nsmoother * NsmootherX,filename)
     WRITE (*,*) 'STORED ' // trim(this)

     THIS      = 'smootherGAPSI' 
     filename  = trim(this) // filext
     FORALL (j=1:NsmootherX,k=1:Nsmoother) theta2(:,(k-1) * NsmootherX + j) = SMOOTHERx(ndxgapSI,1:T,j,k)
     CALL storeEstimates(theta2,T,Nsmoother * NsmootherX,filename)
     WRITE (*,*) 'STORED ' // trim(this)

     ! ! store smoothed X
     ! DO i=1,Nx
     !    filename  = 'smootherX' // trim(int2str(i)) // filext
     !    FORALL (j=1:NsmootherX,k=1:Nsmoother) theta2(:,(k-1) * NsmootherX + j) = SMOOTHERx(i,1:T,j,k)
     !    CALL storeEstimates(theta2,T,Nsmoother * NsmootherX,filename)
     !    WRITE (*,*) 'SMOOTHER STORED X', i
     ! END DO
     DEALLOCATE (theta2)
     ! FINISHED STORING SMOOTHER X

     ! store rest of smoother
     ALLOCATE (theta2(Nsmoother,1:T))

     print *, 'STARTING w/SMOOTHER STORAGE ...'
     filename  = 'smootherLAMBDA' // filext
     CALL savevec(SMOOTHERlambda, filename)
     WRITE (*,*) 'SMOOTHER: STORED LAMBDA' 


     DO i=1,Nsv
        filename  = 'smootherSV' // trim(int2str(i)) // filext
        CALL storeEstimatesTranspose(SMOOTHERsvol(i,:,1:T),T,Nsmoother,filename)
        WRITE (*,*) 'SMOOTHER: STORED SV', i

        FORALL (j=1:T,k=1:Nsmoother) theta2(k,j) = SMOOTHERsvol(i,k,j) - SMOOTHERsvol(i,k,1) 
        filename  = 'smootherDELTASV' // trim(int2str(i)) // filext
        CALL storeEstimatesTranspose(theta2,T,Nsmoother,filename)
        WRITE (*,*) 'SMOOTHER: STORED DELTA SV', i

     END DO

     DEALLOCATE (theta2)

     ! ----------------------------------------------------------------------------------------------------------------------------------------------------------
     ! LINEAR STATE UNCERTAINTY 
     ! ----------------------------------------------------------------------------------------------------------------------------------------------------------
     ! OK to come after storing (and sorting) smoother draws, since this section operates only on posterior moments (not draw by draw)
     if (doLinearUncertainty) then
        print *, '... Linear-state Uncertainty ...'
        allocate (Xsig(Nx,0:T), STATEa(0:T), STATElambda(0:T), STATEsvol(Nsv,0:T))

        ! construct state trajectories (mean of filtered)
        ! forall (i=1:Ny) STATEsigmasqrt(i) = sum(PARTICLEweights(:,T) * sqrt(DRAWsigma(:,i,T)))
        ! STATEa       = sum(PARTICLEweights * DRAWa, 1)
        ! STATElambda  = sum(PARTICLEweights * DRAWlambda, 1)
        ! forall (i=1:Nsv,j=0:T) STATEsvol(i,j) = sum(PARTICLEweights(:,j) * DRAWsvol(:,i,j))

        STATEa       = sum(SMOOTHERa) / Nsmoother
        STATElambda  = sum(SMOOTHERlambda) / Nsmoother
        forall (i=1:Nsv,j=0:T) STATEsvol(i,j) = sum(SMOOTHERsvol(i,:,j)) / Nsmoother
        forall (i=1:Ny) STATEsigmasqrt(i)     = sum(sqrt(SMOOTHERsigma(i,:))) / Nsmoother

        ! debug 
        filename  = 'STATEa'  // filext
        call savevec(STATEa(1:T), filename)
        WRITE (*,*) 'STORED STATEa'
        filename  = 'STATElambda'  // filext
        call savevec(STATElambda(1:T), filename)
        WRITE (*,*) 'STORED STATElambda'
        filename  = 'STATEsvol'  // filext
        call savemat(STATEsvol(:,1:T), filename)
        WRITE (*,*) 'STORED STATEsvol'
        filename  = 'STATEsigmasqrt'  // filext
        call savevec(STATEsigmasqrt, filename)
        WRITE (*,*) 'STORED STATEsigmasqrt'

        ALLOCATE (FcstPersistence(Nsurveys,0:T), STAT=status)
        IF (status /= 0) THEN
           WRITE (*,*) 'Allocation problem (FcstPersistence)'
        END IF


        ! UNCERTAINTY: ALL

        allocate (horizons(Ny), SIndx(Ny), measurementErrorVol(Ny))
        forall (i=1:Ny) horizons(i) = i-1
        forall (i=1:Ny) measurementErrorVol(i) = STATEsigmasqrt(i)
        SIndx    = .True.
        SIndx(1) = .False.
        call linearstateVarianceFilter(horizons, SIndx, measurementErrorVol, size(horizons), Xsig, FcstPersistence, T, Nsurveys, Nx, Nw, sqrtVx0, p, STATEa, STATElambda, STATEsvol, Nsv)
        filename  = 'XSIGall'  // filext
        call savemat(transpose(Xsig(:,1:T)), filename)
        WRITE (*,*) 'STORED XSIGall'
        deallocate (horizons, SIndx, measurementErrorVol)

        ! UNCERTAINTY: SRV
        allocate (horizons(Nsurveys), SIndx(Nsurveys), measurementErrorVol(Nsurveys))
        forall (i=1:Nsurveys) horizons(i) = i
        forall (i=1:Nsurveys) measurementErrorVol(i) = STATEsigmasqrt(i+1)
        SIndx    = .True.
        call linearstateVarianceFilter(horizons, SIndx, measurementErrorVol, size(horizons), Xsig, FcstPersistence, T, Nsurveys, Nx, Nw, sqrtVx0, p, STATEa, STATElambda, STATEsvol, Nsv)
        filename  = 'XSIGsrv'  // filext
        call savemat(transpose(Xsig(:,1:T)), filename)
        WRITE (*,*) 'STORED XSIGsrv'
        deallocate (horizons, SIndx, measurementErrorVol)

        ! UNCERTAINTY: INF
        allocate (horizons(1), SIndx(1), measurementErrorVol(1))
        horizons = 0
        measurementErrorVol = STATEsigmasqrt(1)
        SIndx    = .False.
        call linearstateVarianceFilter(horizons, SIndx, measurementErrorVol, size(horizons), Xsig, FcstPersistence, T, Nsurveys, Nx, Nw, sqrtVx0, p, STATEa, STATElambda, STATEsvol, Nsv)
        filename  = 'XSIGinf'  // filext
        call savemat(transpose(Xsig(:,1:T)), filename)
        WRITE (*,*) 'STORED XSIGinf'
        filename  = 'INFLATIONpersistence'  // filext
        call savemat(transpose(FcstPersistence(:,1:T)), filename)
        WRITE (*,*) 'STORED INFLATIONpersistence'
        deallocate (horizons, SIndx, measurementErrorVol)

        ! UNCERTAINTY: INF and SPFQ4
        allocate (horizons(2), SIndx(2), measurementErrorVol(2))
        horizons = (/0, 5/)
        measurementErrorVol(1) = STATEsigmasqrt(1)
        measurementErrorVol(2) = STATEsigmasqrt(Ny)
        SIndx(1)    = .False.
        SIndx(2)    = .True.
        call linearstateVarianceFilter(horizons, SIndx, measurementErrorVol, size(horizons), Xsig, FcstPersistence, T, Nsurveys, Nx, Nw, sqrtVx0, p, STATEa, STATElambda, STATEsvol, Nsv)
        filename  = 'XSIGinfspflong'  // filext
        call savemat(transpose(Xsig(:,1:T)), filename)
        WRITE (*,*) 'STORED XSIGinfspflong'
        deallocate (horizons, SIndx, measurementErrorVol)

        ! UNCERTAINTY: SPFQ4
        allocate (horizons(1), SIndx(1), measurementErrorVol(1))
        horizons = 5
        measurementErrorVol(1) = STATEsigmasqrt(Ny)
        SIndx(1)    = .True.
        call linearstateVarianceFilter(horizons, SIndx, measurementErrorVol, size(horizons), Xsig, FcstPersistence, T, Nsurveys, Nx, Nw, sqrtVx0, p, STATEa, STATElambda, STATEsvol, Nsv)
        filename  = 'XSIGspflong'  // filext
        call savemat(transpose(Xsig(:,1:T)), filename)
        WRITE (*,*) 'STORED XSIGspflong'
        deallocate (horizons, SIndx, measurementErrorVol)

        deallocate (Xsig, STATElambda, STATEsvol)
        deallocate (FcstPersistence)

     end if
     ! ----------------------------------------------------------------------------------------------------------------------------------------------------------
     ! FINISH UNCERTAINTY CALCULATIONS
     ! ----------------------------------------------------------------------------------------------------------------------------------------------------------


     ! clean up smoother
     DEALLOCATE (SMOOTHERsvol, SMOOTHERlambda, SMOOTHERa, SMOOTHERx)
     DEALLOCATE (SMOOTHERsigma)
     DEALLOCATE (SMOOTHERhinno)

     print *, '... DONE w/SMOOTHER.'
     call hrulefill

     ! ----------------------------------------------------------------------------
     ! FINAL CLEANUP
     ! ----------------------------------------------------------------------------
  END IF ! doSmoother


  DEALLOCATE (y, yNaN)
  ! VSLstreams
  errcode = vsldeletestream(VSLstream)     


  call hrulefill
  WRITE(*,*) 'DONE. BYE, BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill

  STOP

CONTAINS

  ! -----------------------------------------------------------------
  SUBROUTINE getarguments(doInflationNoise, Nparticles, doSmoother, Nsmoother, smootherNparticles, datalabel, doSecondResamplingStep, dofAPFfattail,lambdaN)

    INTENT(INOUT) doInflationNoise, Nparticles, doSmoother, Nsmoother, smootherNparticles, datalabel, doSecondResamplingStep
    INTENT(INOUT) lambdaN

    INTEGER :: counter, dummy
    INTEGER :: Nparticles,Nsmoother, smootherNparticles
    DOUBLE PRECISION :: lambdaN
    LOGICAL :: doSmoother,doInflationNoise
    LOGICAL :: doSecondResamplingStep
    INTEGER, INTENT(INOUT) :: dofAPFfattail
    CHARACTER (LEN=100) :: datalabel
    CHARACTER(len=32) :: arg

    counter = 0

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

    ! datalabel
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, datalabel) 
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

    ! lambdaN
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') dummy
       print *, dummy
       lambdaN = dble(dummy)
       print *, lambdaN
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
SUBROUTINE plefilter(doSecondResamplingStep, dofAPFfattail, doInflationNoise, T, Ny, y, yNaN, Nparticles, PARTICLEweights, logMDD, DRAWpihat, DRAWpifcst, DRAWpipersistence, DRAWxhat, DRAWsqrtSigmaX, DRAWxsig, DRAWxgain, Nxx, NsigmaXX, Nx, NsigmaX, Nw, Ex0, sqrtVx00, DRAWgapnoisehat, p, DRAWa, a0, aV0, DRAWlambda, lambdaAlpha0, lambdaBeta0, lambdaN0, DRAWsvol, Nsv, Eh0, Vh0, DRAWhInno, hvarT, hvarDof, DRAWsigma, sigmaT, sigmaDof, VSLstream, timer)

  USE embox, only : savemat, savevec, int2str, mean, hrulefill

  use gibbsbox, only : drawNDXsysresample, igammaDraws
  use statespacebox, only : DLYAPsqrt
  use blaspack, only : pi, vechU,ivechU, vech, ivech, eye, symmetric, qrot, qrquery
  use cambridgebox, only : drawbetas

  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE

  INTENT(INOUT) :: PARTICLEweights, logMDD, DRAWa, DRAWlambda, DRAWxhat, DRAWsqrtSigmaX, DRAWxsig, DRAWxgain, DRAWsvol, DRAWhInno, DRAWsigma, VSLstream, timer
  INTENT(IN)    :: T,Ny,y,yNaN, Nx, NsigmaX, Nw, Ex0, sqrtVx00, Nsv,Eh0,Vh0, Nparticles, p, lambdaAlpha0, lambdaBeta0, a0, aV0,  hvarT, hvarDof, sigmaT, sigmaDof
  INTENT(IN)    :: lambdaN0
  INTENT(IN) :: Nxx, NsigmaXX
  INTENT(INOUT) :: DRAWpihat, DRAWpifcst
  INTENT(INOUT) :: DRAWpipersistence

  INTEGER :: J, I, K, N, T, Nparticles, Nx, Ny, Nsv, NsigmaX, p, Nw, Nsurveys
  INTEGER :: Nxx, NsigmaXX

  ! OPEN MP
  INTEGER :: TID


  type(progresstimer) :: timer

  logical, intent(in) :: doInflationNoise

  double precision, parameter :: minParticleWeight = 1.0d-12

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(T)              :: logMDD
  DOUBLE PRECISION, DIMENSION(Nparticles,T) :: DRAWa
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: DRAWxsig
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: DRAWxhat
  DOUBLE PRECISION, DIMENSION(NsigmaX,Nparticles, 0:T)  :: DRAWsqrtSigmaX
  DOUBLE PRECISION, DIMENSION(Nx,Ny,Nparticles,T) :: DRAWxgain
  DOUBLE PRECISION, DIMENSION(Nparticles,Nsv,0:T) :: DRAWsvol

  DOUBLE PRECISION, INTENT(INOUT), DIMENSION(Nparticles,0:T,2) :: DRAWgapnoisehat
  DOUBLE PRECISION, DIMENSION(Nparticles,1:T) :: DRAWnoise

  ! APF llf correction
  DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
  DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax


  DOUBLE PRECISION, DIMENSION(Ny-1,Nparticles,0:T)  :: DRAWpihat, DRAWpifcst
  DOUBLE PRECISION, DIMENSION(Ny-1,Nparticles,1:T)  :: DRAWpipersistence
  DOUBLE PRECISION, DIMENSION(Ny-1,Nparticles)      :: pihat, pifcst
  DOUBLE PRECISION, DIMENSION(Ny-1,Nparticles)      :: pipersistence
  DOUBLE PRECISION :: CholeskiGain(Nxx,Ny)

  ! scale parameters
  DOUBLE PRECISION, DIMENSION(Nparticles,T)     :: DRAWlambda
  DOUBLE PRECISION, DIMENSION(Nparticles,Nsv,T) :: DRAWhInno
  DOUBLE PRECISION, DIMENSION(Nsv,Nparticles)   :: hInno ! helper variable, note the transpose
  DOUBLE PRECISION, DIMENSION(Nparticles,Ny,T)  :: DRAWsigma

  ! particles
  DOUBLE PRECISION :: xxposterior(Nxx,Nparticles), xxsig(Nxx,Nparticles), h(Nsv, Nparticles), SVol(Nsv,Nparticles), llf(Nparticles)
  DOUBLE PRECISION :: kernelweights(Nparticles) ! for unnormalized weights
  DOUBLE PRECISION  :: kernelsum, loglikemax
  INTEGER :: ndx(Nparticles)
  DOUBLE PRECISION :: shufflevec(Nparticles)

  ! state space objects
  DOUBLE PRECISION :: xxprior(Nxx), SigmaXX(Nxx,Nxx), logdetSigmaY

  DOUBLE PRECISION :: lambdaAlpha0, lambdaBeta0
  DOUBLE PRECISION :: lambdaN0
  DOUBLE PRECISION, DIMENSION(Nparticles) :: lambdaAlpha, lambdaBeta
  DOUBLE PRECISION, DIMENSION(Nparticles) :: lambda

  DOUBLE PRECISION :: adrift(Nparticles), zdraw(Nparticles),  a0, aV0 ! note: no more "drift" in a but keeping name for sake of comparability with other code
  DOUBLE PRECISION :: zstats(Nparticles)
  DOUBLE PRECISION, PARAMETER :: rejectioncritvalue = 3, unitcircle = 1.0d0 - 1.0d-5
  INTEGER :: counter
  INTEGER, PARAMETER :: maxcount = 10000

  ! scale parameters
  DOUBLE PRECISION :: hvarT(Nsv), sigmaT(Ny)
  INTEGER :: hvarDof(Nsv), sigmaDof(Ny)

  DOUBLE PRECISION :: PREVhvarT(Nparticles,Nsv), PREVsigmaT(Nparticles,Ny)
  INTEGER :: PREVhvarDof(Nsv), PREVsigmaDof(Ny)
  DOUBLE PRECISION :: hDELTA(Nsv,Nparticles), resid(Ny,Nparticles)
  double precision :: xxdraw(Nxx,Nparticles), gapdraw(Nparticles), gaplagdraw(Nparticles) ! used for constructing resid and gap draws used in AR1 update

  ! AR1:
  double precision, dimension(Nparticles) :: PREVa, PREVaV, sqrtPREVaV

  ! SQRT objects
  DOUBLE PRECISION :: sqrtSigmaXX(Nxx,Nxx), sqrtSigmaY(Ny,Ny), Kgain(Nxx,Ny,Nparticles), qrR(Ny+Nxx+Nw,Ny+Nxx)
  DOUBLE PRECISION, DIMENSION(NsigmaXX, Nparticles) :: vecSqrtSigmaXX
  ! DOUBLE PRECISION, DIMENSION(NsigmaX, Nparticles)  :: vecSqrtSigmaX
  INTEGER :: qrLwork

  DOUBLE PRECISION :: sqrtVxx0(Nxx,Nxx), A(Nxx,Nxx), B(Nxx,Nw),  Bsv(Nxx,Nw), C(Ny,Nxx,T), sqrtR(Ny,Ny), ytilde(Ny)
  DOUBLE PRECISION :: Ex0(Nx), sqrtVx00(Nx,Nx)

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
  INTEGER :: ndxGapRElag
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

  ndxGapRElag   = ndxGapSIstop + 1

  ndxgap = (/ ndxgapREstart : ndxGapREstop, ndxgapSIstart : ndxGapSIstop /)

  shockndxTrend = 1
  shockndxGap   = 2

  ! init sufficient statistics of scale parameters
  PREVa  = a0
  PREVaV = aV0

  forall(k=1:Nparticles,j=1:Nsv) PREVhvarT(k,j) = hvarT(j)
  forall(j=1:Nsv)                PREVhvarDof(j) = hvarDof(j)

  forall(k=1:Nparticles,j=1:Ny) PREVsigmaT(k,j) = sigmaT(j)
  forall(j=1:Ny)                PREVsigmaDof(j) = sigmaDof(j)


  ! prepare state space
  ! A
  A = 0.0d0
  ! lagged gap
  A(ndxGapRElag,ndxGapRE) = 1.0d0
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
  lambdaAlpha = lambdaAlpha0
  lambdaBeta  = lambdaBeta0
  lambda      = drawBetas(lambdaAlpha, lambdaBeta, Nparticles, VSLstream) 



  ! a(0)
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, zdraw, 0.0d0, 1.0d0)
  adrift  =  a0 + sqrt(aV0) * zdraw
  do while (ANY(adrift < -1.0d0) .OR. ANY(adrift > 1.0d0))
     errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, zdraw, 0.0d0, 1.0d0) ! a bit brute force to *always* redraw N random variables ..
     where (adrift < -1.0d0)   adrift  =  a0 + sqrt(aV0) * zdraw
     where (adrift > 1.0d0)    adrift  =  a0 + sqrt(aV0) * zdraw
  end do

  ! RB priors for linear states
  FORALL(k=1:Nparticles) xxposterior(1:Nx,k) = Ex0

  ! prepare prior variance of linear states
  vecSqrtSigmaXX       = 0.0d0
  !$OMP PARALLEL DO SHARED(vecSqrtSigmaXX, p, SVol, lambda, adrift, Ny, Nx, Nxx, Nw, Nsv, Nparticles, sqrtVx00, ndxgap, shockndxGap) PRIVATE(gaptransition,gapshock0loadings, sqrtVxx0, ygap0sqrtvariance, errcode) DEFAULT(NONE)
  DO k=1,Nparticles

     sqrtVxx0              = 0.0d0
     sqrtVxx0(1:Nx,1:Nx)   = transpose(sqrtVx00)
     sqrtVxx0(Nxx,Nxx)     = 10.0d0 ! prior variance over gapRElag

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
        sqrtVxx0(ndxgap,ndxgap) = ygap0sqrtvariance
     end if
     vecSqrtSigmaXX(:,k)     = vechU(sqrtVxx0,Nxx)
     ! vecSqrtSigmaX(:,k)      = vechU(sqrtVxx0(1:Nx,1:Nx),Nx) ! equal to vecSqrtSigmaXX(1:NsigmaX,k)

  END DO
  !$OMP END PARALLEL DO 


  FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,0)  = vecSqrtSigmaXX(i,k)
  FORALL(i=1:Nx,k=1:Nparticles)  DRAWxhat(i,k,0)            = xxposterior(i,k) 

  FORALL(i=1:Nsurveys,k=1:Nparticles)  DRAWpihat(i,k,0)     = EX0(ndxTrendRE)
  FORALL(i=1:Nsurveys,k=1:Nparticles)  DRAWpifcst(i,k,0)    = EX0(ndxTrendSI)



  FORALL (i=1:Nsv,k=1:Nparticles) DRAWsvol(k,i,0) = SVol(i,k) ! transpose
  PARTICLEweights  = 1.0d0 / dble(Nparticles)
  logMDD           = 0.0d0

  ! uniform draws for systematic resampling
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


  ! workspace query for qr decomposition
  qrR     = 0.0d0
  qrlwork = qrquery(qrR)

  CALL initprogressbar(timer, 15.0d0)
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

     if (j > 1) then
        ! lambda parameter
        lambda = drawBetas(lambdaAlpha, lambdaBeta, Nparticles, VSLstream) 

	! a parameter
        sqrtPREVaV =  sqrt(PREVaV) 
        errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, zdraw, 0.0d0, 1.0d0)
        adrift  = PREVa + sqrtPREVaV * zdraw

        ! when truncated distribution collapses on boudnary set value to boundary
        ! check whether all mass below -1
        zstats = (-1.0d0 - PREVa) / sqrtPREVaV
        where (zstats .gt. rejectioncritvalue) adrift = -unitcircle ! high prob of being below -1
        zstats = (1.0d0 - PREVa) / sqrtPREVaV
        where (zstats .lt. -rejectioncritvalue) adrift = unitcircle ! low prob of being below 1

        counter = 0
        do while (ANY(abs(adrift) .ge. 1.0d0))
           counter = counter + 1
           if (counter .gt. maxcount) then
              call savevec(zstats, 'z.debug')
              call savevec(adrift, 'a.debug')
              call savevec(PREVa, 'mu.debug')
              call savevec(sqrtPREVaV, 'sig.debug')
              print *, 'aborting b/o too many adrift rejections, t=', j, ' (after ', counter, ' rejection draws)'
              stop 11
           end if
           errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, zdraw, 0.0d0, 1.0d0)
           where (abs(adrift) .ge. 1.0d0)  adrift  = PREVa + sqrtPREVaV * zdraw
           ! where (adrift > 1.0d0)   adrift  = PREVa + sqrtPREVaV * zdraw
        end do
     end if

     !$OMP PARALLEL DO SHARED(xxposterior, xxsig, vecSqrtSigmaXX, lambda, adrift, SVol, DRAWsigma, llf, Nparticles, Nsv, j, y, yNaN, Ny, Nxx, Nw, Nynonan, Nsurveys, p, ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap) SHARED(dofAPFfattail) FIRSTPRIVATE(A,B,C,qrLWORK) PRIVATE(ytilde, sqrtR, Bsv, logdetSigmaY, errcode, xxprior, SigmaXX, sqrtSigmaXX, sqrtSigmaY, qrR, TID) DEFAULT(NONE) SCHEDULE(STATIC)


     DO k = 1,Nparticles

        TID = 0
        !$ TID = OMP_GET_THREAD_NUM()


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
     	xxprior      = 0.0d0
        call DGEMV('n',Nxx,Nxx,1.0d0,A,Nxx,xxposterior(:,k),1,0.0d0,xxprior,1)

        sqrtSigmaXX = ivechU(vecSqrtSigmaXX(:,k),Nxx)

        ! fill directly into qrR
        qrR = 0.0d0
        qrR(1:Ny,1:Ny) = transpose(sqrtR)
        qrR(Ny+Nxx+1:Ny+Nxx+Nw,Ny+1:Ny+Nxx) = transpose(Bsv)
        ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
        call DGEMM('n','t',Nxx,Nxx,Nxx,1.0d0,sqrtSigmaXX,Nxx,A,Nxx,0.0d0,qrR(Ny+1:Ny+Nxx,Ny+1:Ny+Nxx),Nxx)
        ! qrR(Ny+1:Ny+Nx+Nw,1:Ny) = qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx) * C'
        call DGEMM('n','t',Nxx+Nw,Ny,Nxx,1.0d0,qrR(Ny+1:Ny+Nxx+Nw,Ny+1:Ny+Nxx),Nxx+Nw,C(:,:,j),Ny,0.0d0,qrR(Ny+1:Ny+Nxx+Nw,1:Ny),Nxx+Nw)

        ! QR decomposition
        call qrot(qrR, qrLWORK)

        ! map qr into Kalman objects
        sqrtSigmaY   = qrR(1:Ny,1:Ny) ! upper triangular
        sqrtSigmaXX  = qrR(Ny+1:Ny+Nxx,Ny+1:Ny+Nxx) ! upper triangular
        ! Kgain(:,:,k) = transpose(qrR(1:Ny,Ny+1:Ny+Nx))

        ! ytilde = y - C * xprior
        ytilde = y(:,j)
        call DGEMV('n',Ny,Nxx,-1.0d0,C(:,:,j),Ny,xxprior,1,1.0d0,ytilde,1)

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
     !$OMP END PARALLEL DO 


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
           FORALL(k=1:Nparticles) shufflevec(k) = xxposterior(i,ndx(k))
           xxposterior(i,:) = shufflevec
        END DO

        DO i=1,NsigmaXX
           FORALL(k=1:Nparticles) shufflevec(k) = vecSqrtSigmaXX(i,ndx(k))
           vecSqrtSigmaXX(i,:) = shufflevec
        END DO

        DO i=1,Nsv
           FORALL(k=1:Nparticles) shufflevec(k) = h(i,ndx(k))
           h(i,:) = shufflevec
        END DO

        FORALL(k=1:Nparticles) shufflevec(k) = lambda(ndx(k))
        lambda = shufflevec

        FORALL(k=1:Nparticles) shufflevec(k) = lambdaAlpha(ndx(k))
        lambdaAlpha = shufflevec

        FORALL(k=1:Nparticles) shufflevec(k) = lambdaBeta(ndx(k))
        lambdaBeta = shufflevec

        FORALL(k=1:Nparticles) shufflevec(k) = adrift(ndx(k))
        adrift = shufflevec
        FORALL(k=1:Nparticles) shufflevec(k) = PREVa(ndx(k))
        PREVa = shufflevec
        FORALL(k=1:Nparticles) shufflevec(k) = PREVaV(ndx(k))
        PREVaV = shufflevec

        ! reshuffle sufficient statistics for scale parameters

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



     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     !$OMP PARALLEL DO SHARED(xxposterior, xxsig, vecSqrtSigmaXX, lambda, adrift, SVol, Kgain, DRAWsigma, llf, Nparticles, Nsv, j, y, yNaN, Ny, Nx, Nxx, Nw, Nynonan, Nsurveys, p, ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap) SHARED(pihat,pifcst,pipersistence) PRIVATE(CholeskiGain) FIRSTPRIVATE(A,B,C,qrLWORK) PRIVATE(ytilde, sqrtR, Bsv, logdetSigmaY, errcode, xxprior, SigmaXX, sqrtSigmaXX, sqrtSigmaY, qrR, TID) DEFAULT(NONE) SCHEDULE(STATIC)

     DO k = 1,Nparticles

        TID = 0
        !$ TID = OMP_GET_THREAD_NUM()

        xxprior      = 0.0d0
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
        call DGEMV('n',Nxx,Nxx,1.0d0,A,Nxx,xxposterior(:,k),1,0.0d0,xxprior,1)

        ! ------------------------------------------------------------------------
        ! SQRT KALMAN
        ! ------------------------------------------------------------------------
        sqrtSigmaXX = ivechU(vecSqrtSigmaXX(:,k),Nxx)

        ! fill directly into qrR
        qrR = 0.0d0
        qrR(1:Ny,1:Ny) = transpose(sqrtR)
        qrR(Ny+Nxx+1:Ny+Nxx+Nw,Ny+1:Ny+Nxx) = transpose(Bsv)
        ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
        call DGEMM('n','t',Nxx,Nxx,Nxx,1.0d0,sqrtSigmaXX,Nxx,A,Nxx,0.0d0,qrR(Ny+1:Ny+Nxx,Ny+1:Ny+Nxx),Nxx)
        ! qrR(Ny+1:Ny+Nx+Nw,1:Ny) = qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx) * C'
        call DGEMM('n','t',Nxx+Nw,Ny,Nxx,1.0d0,qrR(Ny+1:Ny+Nxx+Nw,Ny+1:Ny+Nxx),Nxx+Nw,C(:,:,j),Ny,0.0d0,qrR(Ny+1:Ny+Nxx+Nw,1:Ny),Nxx+Nw)

        ! QR decomposition
        call qrot(qrR, qrLWORK)

        ! map qr into Kalman objects
        sqrtSigmaY   = qrR(1:Ny,1:Ny) ! upper triangular
        sqrtSigmaXX  = qrR(Ny+1:Ny+Nxx,Ny+1:Ny+Nxx) ! upper triangular
        Kgain(:,:,k) = transpose(qrR(1:Ny,Ny+1:Ny+Nxx))

        ! pull out gain for computation of persistence
        CholeskiGain = Kgain(:,:,k)
        do i = 1,Ny
           ! if (sqrtSigmaY(i,i) < 0.0d0) CholeskiGain(:,i) = - CholeskiGain(:,i)  ! normalize to unit standard deviation shocks
           CholeskiGain(:,i) = CholeskiGain(:,i) / sqrtSigmaY(i,i)  ! normalize to unit-sized shocks
        end do
        forall(i=1:Nsurveys) pipersistence(i,k) = CholeskiGain(ndxTrendRE,1) + (adrift(k) ** i) * CholeskiGain(ndxGapRE,1)


        ! ytilde and logdetSigmaY
        ytilde = y(:,j)
        call DGEMV('n',Ny,Nxx,-1.0d0,C(:,:,j),Ny,xxprior,1,1.0d0,ytilde,1)

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
        xxposterior(:,k) = xxprior
        call DGEMV('N',Nxx,Ny,1.0d0,Kgain(:,:,k),Nxx,ytilde,1,1.0d0,xxposterior(:,k),1)


        ! rotate Kalman gain into space of non-normalized ytilde
        call dtrsm('R', 'U', 'T', 'N', Nxx, Ny, 1.0d0, sqrtSigmaY, Ny, Kgain(:,:,k), Nxx) ! recall: sqrtSigmaY is returned as upper triangular, right factor

        ! remove unit dummies -- can be omitted since sqrtSigmaY not used any further
        ! do i=1,Ny
        !    if (yNaN(i,j)) sqrtSigmaY(i,i) = 0.0d0
        ! end do

        ! compute log-likelihood
        ! llf
        llf(k)       = -0.5d0 * (Nynonan * log(2.0d0 * pi) + logdetSigmaY + sum(ytilde ** 2))

        vecSqrtSigmaXX(:,k) = vechU(sqrtSigmaXX,Nxx)
        ! vecSqrtSigmaX(:,k)  = vechU(sqrtSigmaXX(1:Nx,1:Nx),Nx) ! equal to vecSqrtSigmaXX(1:NsigmaX,k)
        ! vecSqrtSigmaX(:,k)  = vecSqrtSigmaXX(1:NsigmaX,k) -- just using vecSqrtSigmaXX now ...

        ! store state variances 
        call dsyrk('u', 't', Nxx, Nxx, 1.0d0, sqrtSigmaXX, Nxx, 0.0d0, SigmaXX, Nxx)
        forall (i=1:Nxx) xxsig(i,k) = sqrt(SigmaXX(i,i))

        ! ------------------------------------------------------------------------
        ! DONE: SQRT KALMAN
        ! ------------------------------------------------------------------------

        ! ------------------------------------------------------------------------
        ! INFLATION FORECASTS
        ! ------------------------------------------------------------------------
        forall(i=1:Nsurveys) pihat(i,k) = xxposterior(ndxTrendRE,k) + (adrift(k) ** i) * xxposterior(ndxGapRE,k)
        forall(i=1:Nsurveys) pifcst(i,k) = xxposterior(ndxTrendSI,k) + (adrift(k) ** i) * xxposterior(ndxGapSI,k)


     END DO ! k particles
     !$OMP END PARALLEL DO 

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     ! sample noise from Kalman-posteriors for updating SIGMA and gap-AR1 (further below)
     errcode  = vdrnggaussian(VSLmethodGaussian, VSLstream, Nxx * Nparticles, xxdraw, 0.0d0, 1.0d0)
     ! scale xdraw by vecSqrtSigmaXX
     do k=1,Nparticles
        call DTPMV('u','t','n',Nxx,vecSqrtSigmaXX(:,k),xxdraw(:,k),1)
     end do
     ! add mean to xdraw
     forall (k=1:Nparticles,i=1:Nxx) xxdraw(i,k) = xxposterior(i,k) + xxdraw(i,k)

     ! resid = y - C x
     forall (k=1:Nparticles,i=1:Ny) resid(i,k) = y(i,j)
     do k=1,Nparticles
        ! update C
        FORALL (i=1:Nsurveys) C(1+i,ndxGapSI,j) = adrift(k) ** i
        call dgemv('n',Ny,Nxx,-1.0d0,C(:,:,j),Ny,xxdraw(:,k),1,1.0d0,resid(:,k),1) 
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

     FORALL(i=1:Nsigmax,k=1:Nparticles) DRAWsqrtSigmaX(i,k,j)  = vecSqrtSigmaXX(i,k)
     FORALL(i=1:Nx,k=1:Nparticles) DRAWxhat(i,k,j)             = xxposterior(i,k)
     FORALL(i=1:Nx,k=1:Nparticles) DRAWxsig(i,k,j)             = xxsig(i,k) 
     FORALL(i=1:Nx,n=1:Ny,k=1:Nparticles) DRAWxgain(i,n,k,j)   = Kgain(i,n,k) ! Kprime(n,i,k)

     if (.NOT. yNaN(1,j)) then
        ! RE noise
        FORALL(k=1:Nparticles) DRAWnoise(k,j) = y(1,j) - xxposterior(ndxTrendRE,k) - xxposterior(ndxGapRE,k) 
        ! RE: i=1
        FORALL(k=1:Nparticles) DRAWgapnoisehat(k,j,1) = y(1,j) - xxposterior(ndxTrendRE,k) 
        ! SI: i=2
        FORALL(k=1:Nparticles) DRAWgapnoisehat(k,j,2) = xxposterior(ndxGapSI,k) + (1.0d0 - lambda(k)) * DRAWnoise(k,j)
     end if

     FORALL(i=1:Nsurveys,k=1:Nparticles) DRAWpihat(i,k,j)      = pihat(i,k)
     FORALL(i=1:Nsurveys,k=1:Nparticles) DRAWpifcst(i,k,j)     = pifcst(i,k)
     FORALL(i=1:Nsurveys,k=1:Nparticles) DRAWpipersistence(i,k,j) = pipersistence(i,k)

     ! propagate sufficient statistics
     if (Nynonan > 0) then ! nothing to propagate if there was no observed data
        FORALL(k=1:Nparticles,i=1:Nsv) PREVhvarT(k,i)       = PREVhvarT(k,i)    + hDELTA(i,k) ** 2
        if (doInflationNoise) then
           i = 1
           FORALL(k=1:Nparticles)  PREVsigmaT(k,i)  = PREVsigmaT(k,i)   + resid(i,k) ** 2
        end if

        PREVhvarDof      = PREVhvarDof + 1
        WHERE (.NOT. yNaN(:,j))
           PREVsigmaDof     = PREVsigmaDof + 1 
        END WHERE


        ! PROPAGATE AND RESHUFFLE lambda statistics 
        ! propagate 
        lambdaAlpha = lambdaAlpha + lambda * lambdaN0
        lambdaBeta  = lambdaBeta  + (1.0d0 - lambda) * lambdaN0

        FORALL(k=1:Nparticles) shufflevec(k) = lambdaAlpha(k)
        lambdaAlpha = shufflevec

        FORALL(k=1:Nparticles) shufflevec(k) = lambdaBeta(k)
        lambdaBeta = shufflevec


        ! propagate sufficient statistics for AR1-a 
        FORALL(k=1:Nparticles)
           ! note: X and y in regression get scaled by residual Variance
           ! and SV has already been reshuffled
           gapdraw(k)    = xxdraw(ndxgapRE,k) /  SVol(2,k)
           gaplagdraw(k) = xxdraw(ndxgapRElag,k) /  SVol(2,k)

           ! direct computation of posterior mean and variance
           ! - formulas adapted to scalar case
           ! - order is key: posterior mean depends on prior variance which gets overwritten by posterior in next step
           PREVa(k)      = (PREVa(k) + PREVaV(k) * gapdraw(k) * gaplagdraw(k)) /  (1.0d0 + PREVaV(k) * (gaplagdraw(k) ** 2))
           PREVaV(k)     = PREVaV(k) / (1.0d0 + PREVaV(k) * gaplagdraw(k) ** 2)

        END FORALL


        ! particles weights
        PARTICLEweights(:,j) = kernelweights / kernelsum
        ! where (PARTICLEweights(:,j) < minParticleWeight) PARTICLEweights(:,j) = minParticleWeight
        ! PARTICLEweights(:,j) = PARTICLEweights(:,j) / sum(PARTICLEweights(:,j))

        ! resample
        if (doSecondResamplingStep) then

           call drawNDXsysresample(ndx, Nparticles, PARTICLEweights(:,j), Nparticles, uniformdraws(2,j))


           DO i=1,Nx
              FORALL(k=1:Nparticles) shufflevec(k) = xxposterior(i,ndx(k))
              xxposterior(i,:) = shufflevec
           END DO

           DO i=1,NsigmaXX
              FORALL(k=1:Nparticles) shufflevec(k) = vecSqrtSigmaXX(i,ndx(k))
              vecSqrtSigmaXX(i,:) = shufflevec
           END DO
           ! DO i=1,NsigmaX
           !    FORALL(k=1:Nparticles) shufflevec(k) = vecSqrtSigmaX(i,ndx(k))
           !    vecSqrtSigmaX(i,:) = shufflevec
           ! END DO

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

           DO i=1,Nsv
              FORALL(k=1:Nparticles) shufflevec(k) = PREVhvarT(ndx(k),i)
              PREVhvarT(:,i) = shufflevec
           END DO

           DO i=1,Ny
              FORALL(k=1:Nparticles) shufflevec(k) = PREVsigmaT(ndx(k),i)
              PREVsigmaT(:,i) = shufflevec
           END DO

           ! gap-AR1
           FORALL(k=1:Nparticles) shufflevec(k) = PREVa(ndx(k))
           PREVa = shufflevec
           FORALL(k=1:Nparticles) shufflevec(k) = PREVaV(ndx(k))
           PREVaV = shufflevec

        end if ! doSecondResamplingStep

     else ! i.e. Nynonan == 0
        !   PARTICLEweights(:,j) = 1 / dble(Nparticles)
     end if


     CALL progressbarcomment(dble(j) / dble(T), timer, 'Particle Step')

  END DO ! j=1,T

END SUBROUTINE plefilter

! @\newpage\subsection{particlefilter}@
SUBROUTINE particlefilter(doSecondResamplingStep, dofAPFfattail, T, Ny, y, yNaN, Nparticles, PARTICLEweights, logMDD, DRAWxhat, DRAWsqrtSigmaX, Nx, NsigmaX, Nw, Ex0, sqrtVx00, p, thisA, lambda, DRAWh, Nsv, Eh0, Vh0, hInno, sigma, VSLstream)

  USE embox, only : savemat, savevec, int2str, mean, hrulefill

  use gibbsbox, only : drawNDXsysresample
  use statespacebox, only : DLYAPsqrt
  use blaspack, only : pi, vechU,ivechU, vech, ivech, eye, symmetric, qrot, qrquery
  use vslbox
  ! use omp_lib
  ! use timerbox

  IMPLICIT NONE



  INTENT(INOUT) :: PARTICLEweights, logMDD, DRAWxhat, DRAWsqrtSigmaX, DRAWh, VSLstream !, timer
  INTENT(IN)    :: T,Ny,y,yNaN, Nx, NsigmaX, Nw, Ex0, sqrtVx00, Nsv,Eh0,Vh0, Nparticles, p
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
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: DRAWxhat
  DOUBLE PRECISION, DIMENSION(NsigmaX,Nparticles, 0:T)  :: DRAWsqrtSigmaX
  DOUBLE PRECISION, DIMENSION(Nsv,Nparticles,0:T) :: DRAWh

  ! APF llf correction
  DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
  DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax


  ! scale parameters
  DOUBLE PRECISION, INTENT(IN)                 :: thisA
  DOUBLE PRECISION, INTENT(IN)                 :: lambda
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

  DOUBLE PRECISION :: hDELTA(Nsv,Nparticles)


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
  ! A (all constant parameters)
  A = 0.0d0
  ! unit root in trend
  A(ndxTrendRE,ndxTrendRE) = 1.0d0
  A(ndxGapRE,ndxGapRE)     = thisA 
  A(ndxtrendSI,ndxtrendRE)       = 1 - lambda
  A(ndxtrendSI,ndxtrendSI)       = lambda
  A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - lambda) * thisA
  A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = lambda * thisA
  ! kompanion for gap
  IF (p > 1) THEN
     FORALL(j=1:(p-1)) A(ndxGapRE+j,ndxGapRE-1+j) = 1.0d0
  END IF

  ! B
  B                           = 0.0d0
  B(ndxTrendRE,shockndxTrend) = 1.0d0
  B(ndxGapRE,shockndxGap)     = 1.0d0
  B(ndxtrendSI,shockndxTrend) = 1 - lambda
  B(ndxgapSI,shockndxGap)     = 1 - lambda

  ! C
  C         = 0.0d0
  ! inflation
  C(1,ndxTrendRE,:)    = 1.0d0
  C(1,ndxGapRE,:)      = 1.0d0
  ! surveys
  C(2:Ny,ndxTrendSI,:) = 1.0d0
  FORALL (i=1:Nsurveys,j=1:T) C(1+i,ndxGapSI,j) = thisA ** i
  ! zero out missing obs
  do j = 1, T
     do i = 1, Ny
        if (yNaN(i,j)) then 
           C(i,:,j)     = 0.0d0 
        end if
     end do
  end do

  ! Time 0 particles

  ! SV0
  errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, h, 0.0d0, 1.0d0)
  FORALL (j=1:Nsv,k=1:Nparticles) h(j,k) = Eh0(j) + sqrt(Vh0(j)) * h(j,k)  
  SVol = exp(h * 0.5d0)



  ! RB priors for linear states
  FORALL(k=1:Nparticles) xposterior(:,k) = Ex0

  ! prepare prior variance of linear states
  vecSqrtSigmaX       = 0.0d0

  DO k=1,Nparticles ! note: still a TVP model b/o SV

     sqrtVx0         = transpose(sqrtVx00)

     if (abs(thisA) > 1.0d-4 .AND. lambda > 1.0d-4) then 
        gaptransition                  = 0.0d0
        gaptransition(1:p,1:p)         = thisA
        gaptransition(p+1:2*p,1:p)     = (1 - lambda) * thisA
        gaptransition(p+1:2*p,p+1:2*p) = lambda * thisA

        ! Fill in unconditional variance of stationary states
        ! allow for trendshockslopes, thought they are all zero here
        gapshock0loadings                              = 0.0d0
        gapshock0loadings(1,shockndxGap)               = SVol(shockndxGap,k)
        gapshock0loadings(p+1,shockndxGap)             = (1 - lambda) * SVol(shockndxGap,k)

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
  PARTICLEweights  = 1.0d0 / dble(Nparticles)
  logMDD           = 0.0d0

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



     end if

     ! ------------------------------------------------------------------------------------------------------------------------------
     ! END: APF STEP
     ! ------------------------------------------------------------------------------------------------------------------------------


     ! 1) Draw Particles 
     errcode     = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, hDELTA, 0.0d0, 1.0d0)
     forall (k=1:Nparticles,i=1:Nsv) hDELTA(i,k) = hDELTA(i,k) * hInno(i) 
     h           = h + hDELTA
     SVol        = exp(h * 0.5d0)



     ! ------------------------------------------------------------------------------------------------------------------------------
     ! BEGIN: MAIN PARTICLE STEP
     ! ------------------------------------------------------------------------------------------------------------------------------

     DO k = 1,Nparticles

        xprior      = 0.0d0

        ! 2) Fill Particles into state space


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

        end if ! doSecondResamplingStep

     end if

  END DO ! j=1,T

END SUBROUTINE particlefilter

! @\newpage\subsection{plesmoother}@
SUBROUTINE plesmoother(doInflationNoise, Nsmoother, Nparticles, T, p, SMOOTHERa, SMOOTHERlambda, Nsv, SMOOTHERh, SMOOTHERhinno, SMOOTHERsigma, Nparamparticles, PARAMweights, PARAMa, PARAMlambda, PARAMhInno, PARAMsigma, Eh0, Vh0, Ex0, sqrtVx0, Ny, y, yNaN, Nx, NsigmaX, Nw, VSLstream)

  ! use embox
  use embox, only : savevec, savemat, int2str
  use gibbsbox, only : drawNDXpdf ! sysresample
  use blaspack, only : qrquery, qrot, ivechU
  use vslbox
  use omp_lib
  use timerbox

  IMPLICIT NONE

  integer, intent(in) :: Nparticles
  integer, parameter  :: dofAPFfattail = 0
  logical, parameter  :: doSecondResamplingStep = .false.
  logical, intent(in) :: doInflationNoise
  integer, intent(in) :: Nsmoother, T, p, Nparamparticles, Nsv, Ny, Nw, Nx, NsigmaX
  type (vsl_stream_state), intent(inout) :: VSLstream 
  type (vsl_stream_state) :: VSLlocal
  integer :: TID

  DOUBLE PRECISION, INTENT(OUT), DIMENSION(Nsmoother)          :: SMOOTHERlambda
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(Nsmoother)          :: SMOOTHERa
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(Nsv,Nsmoother,0:T)  :: SMOOTHERh

  DOUBLE PRECISION, INTENT(OUT) :: SMOOTHERsigma(Ny,Nsmoother), SMOOTHERhinno(Nsv,Nsmoother)

  ! note: no intent(out) in previous line since smoothed draws for the params won't be needed outside the plesmoother

  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles)         :: PARAMweights
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles,Nsv)     :: PARAMhinno
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles)         :: PARAMa
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles)         :: PARAMlambda
  DOUBLE PRECISION, INTENT(IN), DIMENSION(Nparamparticles, Ny)     :: PARAMsigma

  DOUBLE PRECISION, INTENT(IN) :: Eh0(Nsv), Vh0(Nsv), Ex0(Nx), sqrtVx0(Nx,Nx)

  DOUBLE PRECISION, INTENT(IN), DIMENSION(Ny,T) :: y
  LOGICAL, INTENT(IN), DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Ny,T) :: ytilde

  INTEGER :: k,i,j,n,errcode

  ! objects for the forward filters, conditional on parameters
  INTEGER :: filterndx
  DOUBLE PRECISION :: ufilterdraw(Nsmoother)

  DOUBLE PRECISION, DIMENSION(Nparticles,0:T)     :: PARTICLEweights
  DOUBLE PRECISION, DIMENSION(Nsv,Nparticles,0:T) :: FILTERh
  DOUBLE PRECISION, DIMENSION(Nx,Nparticles,0:T)  :: FILTERxhat
  DOUBLE PRECISION, DIMENSION(Nsigmax,Nparticles,0:T) :: FILTERsqrtSigmaX
  DOUBLE PRECISION, DIMENSION(T) :: logMDD

  double precision, dimension(Nsv) :: prevSV, prevH

  ! filter particles and lindsten helpers
  DOUBLE PRECISION, DIMENSION(Nparticles)     :: w, SVpdf, eta, logdetsqrtLAMBDA
  double precision :: sqrtOmegaHat(Nx,Nx), lambdaHat(Nx), lambda(Nx), sqrtLAMBDA(Nx,Nx), zhat(Nx)
  DOUBLE PRECISION, DIMENSION(Nsv,Nparticles) :: hshock

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

  ! CHARACTER (LEN=200) :: filename

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
  forall (k=1:Nsmoother) SMOOTHERa(k)            = PARAMa(ndx(k)) 
  forall (k=1:Nsmoother) SMOOTHERlambda(k)       = PARAMlambda(ndx(k)) 
  forall (k=1:Nsmoother) SMOOTHERhinno(:,k)      = PARAMhinno(ndx(k),:) 

  ! draw uniforms: Nsmoother x T
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, Nsmoother * T, udraw, 0.0d0, 1.0d0)
  ! draw uniforms for selecting filtered particles
  errcode = vdrnguniform(VSLmethodUniform, VSLstream, Nsmoother, ufilterdraw, 0.0d0, 1.0d0)
  loopCount = 0

  CALL initprogressbar(timer, 15.0d0)

  !$OMP PARALLEL DEFAULT(NONE), &
  !$OMP& SHARED(loopCount,timer,doInflationNoise), &
  !$OMP& SHARED(y,yNaN,T, Ny, Nx, Nw, Nsv, Nsurveys,p,NsigmaX), &
  !$OMP& SHARED(Nsmoother,Nparticles), &
  !$OMP& SHARED(ndxtrendRE, ndxgapRE, ndxtrendSI,shockndxTrend,ndxgapSI,shockndxGap, ndxgapREstart, ndxgapREstop, ndxgapSIstart, ndxgapSIstop, qr1LWORK, qr2LWORK, qr3LWORK), &
  !$OMP& SHARED(udraw, ufilterdraw), &
  !$OMP& PRIVATE(filterndx), & 
  !$OMP& SHARED(SMOOTHERhinno,SMOOTHERsigma), &
  !$OMP& SHARED(SMOOTHERa, SMOOTHERlambda,SMOOTHERh), &
  !$OMP& PRIVATE(VSLlocal,TID), & 
  !$OMP& PRIVATE(j, i, errcode), & 
  !$OMP& FIRSTPRIVATE(C), &
  !$OMP& PRIVATE(invR, invSqrtR, ytilde, Bsv, w, cdf, hshock, SVpdf, eta, sqrtSigmaX, logdetsqrtLAMBDA, sqrtLAMBDA, sqrtOmegaHAt, lambdaHat, lambda, zhat, prevSV, prevH, A, qr1R, qr2R, qr3R), &
  !$OMP& PRIVATE(logMDD, PARTICLEweights, FILTERxhat, FILTERsqrtSigmaX, FILTERh), &
  !$OMP& SHARED(Eh0, Vh0, Ex0, sqrtVx0)

  TID = 0
  !$ TID = OMP_GET_THREAD_NUM()
  errcode = vslnewstream(VSLlocal, vsl_brng_mt2203 + TID, 0)
  if (errcode /= 0) then
     print *,'VSL new stream failed'
     stop 1
  end if


  !$OMP DO 
  DO k=1,Nsmoother

     !$OMP ATOMIC
     loopCount = loopCount + 1

     ! call fwd filter, Note: logMDD remains unused
     call particlefilter(doSecondResamplingStep, dofAPFfattail, T, Ny, y, yNaN, Nparticles, PARTICLEweights, logMDD, FILTERxhat, FILTERsqrtSigmaX, Nx, NsigmaX, Nw, Ex0, sqrtVx0, p, SMOOTHERa(k), SMOOTHERlambda(k), FILTERh, Nsv, Eh0, Vh0, SMOOTHERhInno(:,k), SMOOTHERsigma(:,k), VSLlocal)

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

     ! setup C
     FORALL (i=1:Nsurveys,j=1:T) C(1+i,ndxGapSI,j) = SMOOTHERa(k) ** i
     DO j=1,T
        DO i = 1, Ny
           if (yNaN(i,j)) C(i,:,j) = 0.0d0
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

     DO j=T-1,0,-1


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
        Bsv(ndxTrendSI,shockndxTrend)    = (1 - SMOOTHERlambda(k)) * prevSV(shockndxTrend)
        Bsv(ndxGapRE,shockndxGap)        = prevSV(shockndxGap)
        Bsv(ndxGapSI,shockndxGap)        = (1 - SMOOTHERlambda(k)) * prevSV(shockndxGap)


        ! set up A(t+1)
        A = 0.0d0
        A(ndxTrendRE,ndxTrendRE)         = 1.0d0
        A(ndxGapRE,ndxGapRE)             = SMOOTHERa(k)
        A(ndxtrendSI,ndxtrendRE)         = 1 - SMOOTHERlambda(k)
        A(ndxtrendSI,ndxtrendSI)         = SMOOTHERlambda(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 -  SMOOTHERlambda(k)) * SMOOTHERa(k)
        A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) =  SMOOTHERlambda(k) * SMOOTHERa(k)

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

	! try computing w in logs as much as possible
        w = -0.5d0 * (SVpdf + eta) - logdetsqrtLAMBDA  + log(PARTICLEweights(:,j))

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



     END DO ! j=T-1,0,-1

     CALL progressbarcomment(dble(loopCount) / dble(Nsmoother), timer, 'Particle Smoother Step')

  END DO ! k
  !$OMP END DO 
  errcode = vsldeletestream(VSLlocal)     
  !$OMP END PARALLEL

  ! convert log-variances back into vols
  SMOOTHERh = exp(SMOOTHERh * 0.5d0)

END SUBROUTINE plesmoother

! @\newpage\subsection{particleSmootherX}@
SUBROUTINE particleSmootherX(T, Ny, y, yNaN, Nsv, Nsmoother, SMOOTHERa, p, SMOOTHERlambda, SMOOTHERsvol, NsmootherX, SMOOTHERx, Nx, Nw, Ex0, sqrtVx00, SMOOTHERsigma)

  ! use embox
  ! use gibbsbox
  use statespacebox, only : DLYAP, samplerA3B3C3noisenan

  use vslbox
  use timerbox
  use omp_lib

  IMPLICIT NONE

  INTENT(INOUT) :: SMOOTHERx
  INTENT(IN)    :: T, Ny, Ex0, sqrtVx00, y, Nsv, Nsmoother, SMOOTHERlambda, SMOOTHERa, SMOOTHERsvol, NsmootherX, Nx, Nw, p, SMOOTHERsigma

  INTEGER :: T, Ny, Nsv, Nsmoother, NsmootherX, Nx, Nw, Nsurveys, p

  INTEGER :: j,k,i

  type(progresstimer) :: timer

  DOUBLE PRECISION, DIMENSION(Ny,T) :: y
  LOGICAL, DIMENSION(Ny,T) :: yNaN
  DOUBLE PRECISION, DIMENSION(Nsv,Nsmoother,0:T) :: SMOOTHERsvol
  DOUBLE PRECISION, DIMENSION(Nsmoother)     :: SMOOTHERlambda
  DOUBLE PRECISION, DIMENSION(Nsmoother)     :: SMOOTHERa
  DOUBLE PRECISION, DIMENSION(Nx,0:T,NsmootherX,Nsmoother) :: SMOOTHERx

  ! particles
  DOUBLE PRECISION :: x(Nx,0:T), xshock(Nx,T), ynoise(Ny,T)

  ! state space objects
  DOUBLE PRECISION :: SMOOTHERsigma(Ny,Nsmoother)
  DOUBLE PRECISION :: Ex0(Nx), sqrtVx00(Nx,Nx), sqrtVx0(Nx,Nx), A(Nx,Nx,T), B(Nx,Nw,T), C(Ny,Nx,T), noisevol(Ny,T)
  DOUBLE PRECISION :: ygap0variance(2*p,2*p), gaptransition(2*p,2*p), gapshock0loadings(2*p,Nw)

  ! VSL
  TYPE (VSL_STREAM_STATE) :: VSLstream
  INTEGER :: errcode
  INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0

  ! OMP
  INTEGER :: TID, NTHREADS

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

  !$OMP PARALLEL SHARED(y, yNaN, SMOOTHERx, Ex0, Ny, p, Nx, Nw, Nsurveys, T, SMOOTHERsigma, SMOOTHERlambda, SMOOTHERa, SMOOTHERsvol, Nsv, NsmootherX, Nsmoother,ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap) PRIVATE(x,gaptransition,gapshock0loadings,ygap0variance,xshock,errcode,TID,NTHREADS,timer,VSLstream,ynoise,noisevol,i,j) FIRSTPRIVATE(A,B,C,sqrtVx0) DEFAULT(NONE)


  NTHREADS = 1
  TID = 0
  !$ TID = OMP_GET_THREAD_NUM()
  !$ NTHREADS = OMP_GET_NUM_THREADS()
  errcode = vslnewstream(VSLstream, vsl_brng_mt2203 + TID + 1, 0)  
  if (errcode /= 0) then
     print *,'VSL new stream failed'
     stop 1
  end if

  if (TID == 0) call initprogressbar(timer, 15.0d0)

  !$OMP DO 
  DO k=1,Nsmoother

     ! if (TID == 0) print *, 'k', k

     ! update A 
     forall (j=1:T) A(ndxGapRE,ndxGapRE,j)      = SMOOTHERa(k)
     forall (j=1:T) A(ndxtrendSI,ndxtrendRE,j)  = 1 - SMOOTHERlambda(k)
     forall (j=1:T) A(ndxtrendSI,ndxtrendSI,j)  = SMOOTHERlambda(k)
     forall (j=1:T) A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop,j) = (1 - SMOOTHERlambda(k)) * SMOOTHERa(k)
     forall (j=1:T) A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop,j) = SMOOTHERlambda(k)* SMOOTHERa(k)

     ! update noisevol
     forall(j=1:T) noisevol(:,j) = sqrt(SMOOTHERsigma(:,k)) 


     ! update C
     FORALL (i=1:Nsurveys, j=1:T) C(1+i,ndxGapSI,j) = SMOOTHERa(k) ** i
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
     forall (j=1:T) B(ndxtrendSI,shockndxTrend,j)    = (1 - SMOOTHERlambda(k)) * SMOOTHERsvol(shockndxTrend,k,j)
     forall (j=1:T) B(ndxgapSI,shockndxGap,j)        = (1 - SMOOTHERlambda(k)) * SMOOTHERsvol(shockndxGap,k,j)


     ! prior variance for linear states
     gaptransition                  = 0.0d0
     gaptransition(1:p,1:p)         = SMOOTHERa(k)
     gaptransition(p+1:2*p,1:p)     = (1 - SMOOTHERlambda(k)) * SMOOTHERa(k)
     gaptransition(p+1:2*p,p+1:2*p) = SMOOTHERlambda(k) * SMOOTHERa(k)

     gapshock0loadings                   = 0.0d0
     gapshock0loadings(1,shockndxGap)    = SMOOTHERsVOL(shockndxGap,k,0)
     gapshock0loadings(p+1,shockndxGap)  = (1 - SMOOTHERlambda(k)) * SMOOTHERsVOL(shockndxGap,k,0)

     CALL DLYAP(ygap0variance, gaptransition, gapshock0loadings, 2 * p, Nw, errcode) 
     if (errcode /= 0) then
        write (*,*) 'DLYAP error (ygap0variance, smootherX)', errcode
        stop 1
     end if

     CALL DPOTRF('L', 2 * p, ygap0variance, 2 * p, errcode) 
     if (errcode /= 0) then
        write (*,*) 'DPOTRF error (ygap0variance, smootherX)', errcode
        stop 1
     end if
     ! zero out the upper triangular
     ! FORALL (i=1:2*p-1) ygap0variance(i,i+1:2*p) = 0.0d0
     sqrtVx0(ndxgapREstart:ndxgapREstop,ndxgapREstart:ndxgapREstop) = ygap0variance(1:p,1:p)
     sqrtVx0(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = ygap0variance(p+1:2*p,1:p)
     sqrtVx0(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = ygap0variance(p+1:2*p,p+1:2*p)


     ! 3) KALMAN SMOOTHER: NsmootherX draws
     DO j=1,NsmootherX
        ! note: xshock and ynoise will note be used
        call samplerA3B3C3noisenan(x,xshock,ynoise,y,yNaN,T,Ny,Nx,Nw,A,B,C,noisevol,Ex0,sqrtVx0,VSLstream,errcode)
        FORALL (i=1:Nx) SMOOTHERx(i,:,j,k)  = x(i,:)

     END DO ! j


     ! timer assumes that equal chunks are spread across workers
     if (TID == 0) call progressbarcomment(dble(k * NTHREADS) / dble(Nsmoother), timer, 'Smoothing Particles X')

  END DO ! k
  !$OMP END DO 


  errcode = vsldeletestream(VSLstream)     


  !$OMP END PARALLEL


END SUBROUTINE particleSmootherX


! @\newpage\subsection{linearstateVarianceFilter}@
SUBROUTINE linearstateVarianceFilter(horizons, SIndx, measurementErrorVol, Ny, Xsig, FcstPersistence, T, Nsurveys, Nx, Nw, sqrtVx0, p, STATEa, STATElambda, STATEsvol, Nsv)

  USE embox, only : savemat, savevec, int2str, mean, hrulefill

  use statespacebox, only : DLYAPsqrt
  use blaspack, only : eye, symmetric, qrot, qrquery

  use vslbox
  use omp_lib
  ! use timerbox

  IMPLICIT NONE

  INTENT(OUT) :: Xsig
  INTENT(IN)  :: STATEa, STATElambda, STATEsvol
  INTENT(IN)  :: T,Ny, Nx, Nw, sqrtVx0, Nsv, p
  INTENT(IN)  :: horizons, SIndx, measurementErrorVol

  INTENT(OUT) :: FcstPersistence
  INTENT(IN)  :: Nsurveys

  INTEGER :: J, I, T, Nx, Ny, Nsv, p, Nw
  INTEGER :: Nsurveys

  INTEGER :: errcode


  ! type(progresstimer) :: timer

  DOUBLE PRECISION, DIMENSION(Ny) :: horizons, measurementErrorVol
  LOGICAL, DIMENSION(Ny) :: SIndx

  DOUBLE PRECISION, DIMENSION(Nx, 0:T)  :: Xsig
  DOUBLE PRECISION, DIMENSION(Nsurveys, 0:T) :: FcstPersistence

  DOUBLE PRECISION, DIMENSION(0:T) :: STATElambda
  DOUBLE PRECISION, DIMENSION(0:T) :: STATEa
  DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: STATEsvol

  ! SQRT objects
  DOUBLE PRECISION :: sqrtSigmaX(Nx,Nx), SigmaX(Nx,Nx), qrR(Ny+Nx+Nw,Ny+Nx) ! , sqrtSigmaY(Ny,Ny)
  INTEGER :: qrLwork
  DOUBLE PRECISION :: Kgain(Nx,Ny), sqrtSigmaY(Ny,Ny)

  ! state space objects
  DOUBLE PRECISION :: sqrtVx0(Nx,Nx), A(Nx,Nx), B(Nx,Nw), Bsv(Nx,Nw), C(Ny,Nx), sqrtR(Ny,Ny)
  DOUBLE PRECISION :: ygap0sqrtvariance(2*p,2*p), gapshock0loadings(2*p,Nw), gaptransition(2*p,2*p)


  ! ! CHARACTER (LEN=200) :: filename

  ! index variables for state space
  INTEGER :: ndxTrendRE, ndxGapRE, ndxTrendSI, ndxGapSI, ndxGapREstart, ndxGapREstop, ndxGapSIstart, ndxGapSIstop, shockndxTrend, shockndxGap
  INTEGER :: ndxgap(2*p)

  ! CHARACTER (LEN=100) :: filename

  FcstPersistence = 0.0d0

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

  ! trend loadings
  where (SIndx)
     C(:,ndxTrendSI) = 1.0d0
  elsewhere
     C(:,ndxTrendRE)    = 1.0d0
  end where


  ! sqrtR
  sqrtR = 0.0d0
  forall (i=1:Ny) sqrtR(i,i) = measurementErrorVol(i)


  ! init time j= 0
  j = 0
  ! prepare prior variance of linear states
  sqrtSigmaX                    = transpose(sqrtVx0)

  if (abs(STATEa(j)) > 1.0d-8) then
     gaptransition                  = 0.0d0
     gaptransition(1:p,1:p)         = STATEa(j)
     gaptransition(p+1:2*p,1:p)     = (1 - STATElambda(j)) * STATEa(j)
     gaptransition(p+1:2*p,p+1:2*p) = STATElambda(j) * STATEa(j) 

     ! Fill in unconditional variance of stationary states
     ! allow for trendshockslopes, thought they are all zero here
     gapshock0loadings                              = 0.0d0
     gapshock0loadings(1,shockndxGap)               = STATEsvol(shockndxGap,j)
     gapshock0loadings(p+1,shockndxGap)             = (1 - STATElambda(j)) * STATEsvol(shockndxGap,j)

     CALL DLYAPsqrt(ygap0sqrtvariance, gaptransition, gapshock0loadings, 2 * p, Nw, errcode) 
     if (errcode /= 0) then
        write (*,*) 'DLYAP error (ygap0sqrtvariance -- linearstateVariance)', errcode
        stop 1
     end if
     sqrtSigmaX(ndxgap,ndxgap) = ygap0sqrtvariance
  end if

  ! store state variances 
  call dsyrk('u', 't', Nx, Nx, 1.0d0, sqrtSigmaX, Nx, 0.0d0, SigmaX, Nx)
  forall (i=1:Nx) xsig(i,j) = sqrt(SigmaX(i,i))


  ! workspace query for qr decomposition
  qrR = 0.0d0
  qrlwork = qrquery(qrR)

  ! CALL initprogressbar(timer, 15.0d0)
  DO j=1,T

     ! 2) Fill Particles into state space

     ! update A 

     A(ndxGapRE,ndxGapRE)           = STATEa(j)
     A(ndxtrendSI,ndxtrendRE)       = 1 - STATElambda(j)
     A(ndxtrendSI,ndxtrendSI)       = STATElambda(j)
     A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - STATElambda(j)) * STATEa(j)
     A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = STATElambda(j) * STATEa(j)

     ! update gap loadings in C
     where (SIndx)
        C(:,ndxGapSI) = STATEa(j) ** horizons
     elsewhere
        C(:,ndxGapRE) = STATEa(j) ** horizons
     end where

     ! update B
     B(ndxtrendSI,shockndxTrend)    = 1 - STATElambda(j)
     B(ndxgapSI,shockndxGap)        = 1 - STATElambda(j)

     ! Bsv
     FORALL (i=1:Nsv) Bsv(:,i)  = B(:,i) * STATEsvol(i,j) 

     ! 3) Kalman Filter

     ! ------------------------------------------------------------------------
     ! SQRT KALMAN
     ! ------------------------------------------------------------------------

     ! fill in qrR
     qrR = 0.0d0
     qrR(1:Ny,1:Ny) = transpose(sqrtR)
     qrR(Ny+Nx+1:Ny+Nx+Nw,Ny+1:Ny+Nx) = transpose(Bsv)
     ! qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) = sqrtSigmaX' * A' (sqrtSigmaX is already stored as transpose)
     call DGEMM('n','t',Nx,Nx,Nx,1.0d0,sqrtSigmaX,Nx,A,Nx,0.0d0,qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx),Nx)
     ! qrR(Ny+1:Ny+Nx+Nw,1:Ny) = qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx) * C'
     call DGEMM('n','t',Nx+Nw,Ny,Nx,1.0d0,qrR(Ny+1:Ny+Nx+Nw,Ny+1:Ny+Nx),Nx+Nw,C,Ny,0.0d0,qrR(Ny+1:Ny+Nx+Nw,1:Ny),Nx+Nw)

     ! QR decomposition
     call qrot(qrR, qrLWORK)

     ! map qr into Kalman objects
     sqrtSigmaY   = qrR(1:Ny,1:Ny) ! upper triangular
     sqrtSigmaX   = qrR(Ny+1:Ny+Nx,Ny+1:Ny+Nx) ! upper triangular
     Kgain(:,:) = transpose(qrR(1:Ny,Ny+1:Ny+Nx))
     ! sign Kgain
     do i=1,Ny
        ! if (sqrtSigmaY(i,i) < 0.0d0) then ! normalize to unit standard deviation shocks
        !    Kgain(:,i) = - Kgain(:,i)
        ! end if
        Kgain(:,i) = Kgain(:,i) / sqrtSigmaY(i,i) ! normalize to unit-sized shocks
     end do
     ! store FcstPersistence
     do i=1,Nsurveys
        FcstPersistence(i,j) = Kgain(ndxTrendRE,1) + (STATEa(j) ** i) * Kgain(ndxGapRE,1)
     end do

     ! store state variances 
     call dsyrk('u', 't', Nx, Nx, 1.0d0, sqrtSigmaX, Nx, 0.0d0, SigmaX, Nx)
     forall (i=1:Nx) xsig(i,j) = sqrt(SigmaX(i,i))

     ! ------------------------------------------------------------------------
     ! DONE: SQRT KALMAN
     ! ------------------------------------------------------------------------



     !    CALL progressbarcomment(dble(j) / dble(T), timer, 'Particle Step')

  END DO ! j=1,T

END SUBROUTINE linearstateVarianceFilter
