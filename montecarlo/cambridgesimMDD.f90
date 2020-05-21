PROGRAM main

  USE embox, only : hrulefill, loadmatrix, loadarray1, savemat, savematlogical, savevec, storeEstimates, loft, timestampstr, es30d16, int2str
  USE blaspack, only : eye

  USE cambridgebox, only : MDDthetaTVPlambdaTVP, MDDthetaCONSTlambdaTVP, MDDthetaTVPlambdaCONST, MDDthetaCONSTlambdaCONST, generateUCSVSIdata, ucsvSI
  USE vslbox
  USE timerbox
  USE omp_lib

  IMPLICIT NONE

  ! ----------------------------------------------------------------------------------

  CHARACTER (LEN=200) :: modellabel ='SIthetaTVPlambdaTVP'
  INTEGER  :: Ngrid = 500

  logical, parameter :: doTimestamp = .false., doSmoother = .false., doGains = .false.
  logical :: doInflationNoise = .true.
  logical :: doStoreDGP       = .false.

  INTEGER, PARAMETER :: p = 1, Nsurveys = 5, Ny = Nsurveys + 1, Nx = 2 * (1 + p), NsigmaX = Nx * (Nx + 1) / 2, Nw = 2, Nsv = 2
  INTEGER, PARAMETER :: Nxx = Nx + 1, NsigmaXX = Nxx * (Nxx + 1) / 2 ! for state space with lagged gap

  INTEGER, PARAMETER :: ndxtrendRE = 1 , ndxgapRE = 2 , ndxtrendSI = 1 + p + 1 , ndxgapSI = ndxtrendSI + 1
  ! DOUBLE PRECISION, DIMENSION(Nsurveys), PARAMETER :: horizons = (/1,2,3,4,5/)

  INTEGER, PARAMETER :: dof0 = 3
  double precision :: noisevol0(Ny), hinno0(Nsv), siga0, siglambda0

  INTEGER :: Nparticles ! Nsmoother, NsmootherX, Nmixturedraws
  INTEGER :: T,i,j,status 

  DOUBLE PRECISION, DIMENSION(Nx,Nx) :: sqrtVx0,sqrtVx00
  DOUBLE PRECISION, DIMENSION(Nx)    :: Ex0
  DOUBLE PRECISION, DIMENSION(Nsv)   :: SVar0, Eh0, Vh0

  DOUBLE PRECISION :: lambda0, lambda0V
  DOUBLE PRECISION :: a0, a0V 

  ! arrays for initial values
  double precision, allocatable, dimension(:)   :: lambda00, a00 
  double precision, allocatable, dimension(:,:) :: x00
  double precision, allocatable, dimension(:,:) :: draws2

  DOUBLE PRECISION, parameter :: lambdaAlpha = 1.0d0, lambdaBeta = 1.0d0

  ! priors for scale parameters
  DOUBLE PRECISION :: sigaT, siglambdaT, hvarT(Nsv), sigmaT(Ny)
  INTEGER :: sigaDof, siglambdaDof, hvarDof(Nsv), sigmaDof(Ny)


  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: y
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: yNaN

  ! grid stuff
  integer, parameter :: Nmodels = 3 ! 0:3
  integer :: thismodel
  ! define a derived type for cambridgeSI true model parameters and states
  type(ucsvSI(:,:,:,:)), allocatable :: DGP
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: logMDD

  integer :: MDDgrid = 20
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)       :: logMDDscore

  TYPE(progresstimer) :: itertimer
  INTEGER :: iterCount, iterTotal

  CHARACTER (LEN=200) :: filename, fileXT, datalabel

  ! VSL Random Stuff
  type (vsl_stream_state) :: VSLstream, defaultVSLstream
  integer :: seed
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

  MDDgrid       = 20
  Nparticles    = 10 ** 4
  Ngrid         = 16

  T = 200 
  call getarguments(Nparticles,Ngrid,MDDgrid,doInflationNoise,T)


  ! ----------------------------------------------------------------------------------
  ! INIT
  ! ----------------------------------------------------------------------------------

  ! INIT OMP
  NTHREADS = 1
  !$OMP PARALLEL SHARED(NTHREADS)
  !$ NTHREADS = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  print *, "Number of Threads:", NTHREADS

  errcode = vslnewstream(defaultVSLstream, vsl_brng_mt2203, 0)  
  if (errcode /= 0) then
     print *,'VSLstream failed to init'
     stop 1
  end if
  WRITE(*,'(a25, i20, i20)') 'LAUNCHING VSLSTREAM ', defaultVSLstream%descriptor1, defaultVSLstream%descriptor2


  ! runtime parameters :end: 





  datalabel        = 'simdataT' // trim(int2str(T)) 
  ! Key DGP parameters
  hInno0       = 0.2d0
  siglambda0   = 0.1d0
  siga0        = 0.1d0
  noisevol0    = sqrt(0.1d0) ! BASELINE
  ! noisevol0    = 0.01d0 ! RESET
 
  ! trivial since Nsv = 1 here
  Svar0    = (/ 0.6 / 3.0d0, 0.6 * 2.0d0 / 3.0d0 /)
  Vh0      = 10.0d0
  Eh0      = log(Svar0) - Vh0 * 0.5d0

  hvarDof = dof0
  hvarT   = (hInno0 ** 2) * (dble(hvarDof) - 2.0d0)

  ! lambda
  lambda0      = .5d0
  lambda0V     = 1.0d0
  siglambdaDof = dof0
  siglambdaT   = (siglambda0 ** 2) * (dble(siglambdaDof) - 2.0d0)

  ! a
  a0        = 0.0d0
  a0V       = 1.0d0
  sigaDof   = dof0
  sigaT     = (siga0 ** 2) * (dble(sigaDof) - 2.0d0)

  ! Linear prior
  Ex0       = 0.0d0
  Ex0(1)    = 2.0d0
  Ex0(3)    = 2.0d0
  
  ! sqrtVx0, expressed as lower  triangular-choleski factor (sqrtVx0 * sqrtVx0')
  call eye(sqrtVx0, 10.0d0) ! non-zero prior gap-variance matters only for initial conditions for sampling lagged gap as used for estimating AR(1) coefficient 
  sqrtVx0(1,1) = 100.0d0 ! sqrt(2.0d0)
  sqrtVx0(2+p,2+p) = 100.0d0 ! sqrt(2.0d0) 
  
  sigmaDof = 20 ! dof0
  sigmaT   = (noisevol0 ** 2) * (dble(sigmaDof) - 2.0d0)

  ! -------------------------------------------------------------
  ! START: draw initial values for theta and lambda
  ! -------------------------------------------------------------
  allocate(lambda00(Ngrid), a00(Ngrid))
  
  ! lambda00
  errcode = vdrnguniform(VSLmethodUniform, defaultVSLstream, Ngrid, lambda00, 0.0d0, 1.0d0 )

  ! a00
  ! errcode = vdrnguniform(VSLmethodUniform, defaultVSLstream, Ngrid, a00, -.9d0, 0.9d0 ) ! drawing too close to boundary leads to issues with simulations later
  a00 = a0

  ! x00
  sqrtVx00 = sqrtVx0
  sqrtVx00(1,1) = 0.0d0 ! to avoid drawing huge numbers for the trend component ...
  sqrtVx00(3,3) = 2.0d0
  ! call savemat(sqrtVx00, 'sqrtVx00.debug')

  allocate(x00(Nx,Ngrid))
  forall (i=1:Ngrid) x00(:,i) = Ex0 
  allocate(draws2(Nx,Ngrid))
  errcode = vdrnggaussian(VSLmethodGaussian, defaultVSLstream, Nx * Ngrid, draws2, 0.0d0, 1.0d0)
  DO i=1,Ngrid
     call DGEMV('N',Nx,Nx,1.0d0,sqrtVx00,Nx,draws2(:,i),1,1.0d0,x00(:,i),1)
  END DO
  deallocate(draws2)

  ! delete defaultVSLstream after last use
  errcode = vsldeletestream(defaultVSLstream)

  ! -------------------------------------------------------------
  ! STOP: draw initial values for theta and lambda
  ! -------------------------------------------------------------


  ! REPORT PARAMETERS TO SCREEN
  CALL HRULEFILL
  print *, 'data=       ' // datalabel
  print *, 'model=      ' // modellabel
  print *, 'Ngrid=      ', Ngrid
  print *, 'MDDgrid=    ', MDDgrid
  print *, 'Ny=         ', Ny
  print *, 'T=          ', T
  print *, 'Nparticles= ', Nparticles
  ! print *, 'Nsmoother= ', Nsmoother
  print *, 'p=          ', p
  if (doInflationNoise) then
     print *, 'Model variant with Noise in Inflation'
  else
     print *, 'Model variant WITHOUT Noise in Inflation'
  end if
  do i=1,Ny
     WRITE(*,'(a11,i1,a3, f6.4)') 'noisevol0 (', i, '): ', noisevol0(i)
  end do
  CALL HRULEFILL

  ALLOCATE (logMDD(T,0:Nmodels,Ngrid))

  CALL initprogressbar(itertimer, 15.0d0)
  iterCount = 0
  iterTotal = Ngrid * (Nmodels + 1) * MDDgrid

  !$OMP PARALLEL DO SHARED(doInflationNoise,logMDD,iterCount,itertimer,Nparticles,T,Ngrid,Ex0,sqrtVx0,A0,A0V,SigAT,SigAdof,LAMBDA0,LAMBDA0V,siglambdaT,siglambdaDof,Eh0,Vh0,hvarT,hvarDof,sigmaT,sigmaDof) SHARED(doStoreDGP) PRIVATE(TID,VSLstream,seed,errcode) PRIVATE(thismodel) SHARED(lambda00,a00,x00) FIRSTPRIVATE(y,yNaN,logMDDscore) DEFAULT(FIRSTPRIVATE) SCHEDULE(DYNAMIC)
  ! note: openMP does not like explicit clauses with DGP (since it is parameterized type), but defaulting to FIRSTPRIVATE works ...

  DO i=1,Ngrid

     ! VSL
     TID = 0
     !$ TID = OMP_GET_THREAD_NUM()
     seed = 1
     errcode = vslnewstream(VSLstream, vsl_brng_mt2203 + i, seed)  
     if (errcode /= 0) then
        print *,'VSLstream failed to init'
        stop 1
     end if
     ! WRITE(*,'(a25, i20, i20)') 'LAUNCHING VSLSTREAM ', VSLstream%descriptor1, VSLstream%descriptor2
     ! print *, 'vsl_brng', vsl_brng_mt2203


     ! generate data
     ALLOCATE (logMDDscore(T,MDDgrid)) 

     allocate (ucsvSI(Ny,Nx,Nsv,T)::DGP)
     ALLOCATE (y(Ny,T), yNaN(Ny,T), STAT=status)
     IF (status /= 0) THEN
        WRITE (*,*) 'Allocation problem (Y)'
     END IF


     ! generate data
     call generateUCSVSIdata(dgp, doInflationNoise, T, Ny, x00(:,i), Eh0, hinno0, a00(i), siga0, lambda00(i), siglambda0, noisevol0, VSLstream)

     ! store individual DGP's
     if (doStoreDGP) then
        filename = 'theta_grid' // trim(int2str(i)) // '.dgp.dat'
        call savevec(dgp%theta, filename)
        filename = 'lambda_grid' // trim(int2str(i)) // '.dgp.dat'
        call savevec(dgp%lambda, filename)
        filename = 'x_grid' // trim(int2str(i)) // '.dgp.dat'
        call savemat(dgp%x, filename)
        filename = 'y_grid' // trim(int2str(i)) // '.dgp.dat'
        call savemat(dgp%y, filename)
        filename = 'noise_grid' // trim(int2str(i)) // '.dgp.dat'
        call savemat(dgp%noise, filename)
     end if

     yNaN = .false.
     y = dgp%y

     do thismodel = 0,Nmodels


        do j=1,MDDgrid

           !$OMP ATOMIC
           iterCount = iterCount + 1

           logMDDscore = 0.0d0
           select case (thisModel)
           case (2)
              CALL MDDthetaTVPlambdaTVP(doInflationNoise,T, logMDDscore(:,j), Ny, y, yNaN, Nparticles, Nx, NsigmaX, Nw, Ex0, sqrtVx0, p, a0, a0V, sigaT, sigaDof, lambda0, lambda0V, siglambdaT, siglambdaDof, Nsv, Eh0, Vh0, hvarT, hvarDof, sigmaT, sigmaDof, VSLstream)
           case (0)
              CALL MDDthetaCONSTlambdaTVP(doInflationNoise,T, logMDDscore(:,j), Ny, y, yNaN, Nparticles, Nxx, NsigmaXX, Nx, Nw, Ex0, sqrtVx0, p, a0, a0V, lambda0, lambda0V, siglambdaT, siglambdaDof, Nsv, Eh0, Vh0, hvarT, hvarDof, sigmaT, sigmaDof, VSLstream)
           case (1)
              CALL MDDthetaCONSTlambdaCONST(doInflationNoise,T, logMDDscore(:,j), Ny, y, yNaN, Nparticles, Nxx, NsigmaXX, Nx, Nw, Ex0, sqrtVx0, p, a0, a0V, lambdaAlpha, lambdaBeta, Nsv, Eh0, Vh0, hvarT, hvarDof, sigmaT, sigmaDof, VSLstream)
           case (3)
              CALL MDDthetaTVPlambdaCONST(doInflationNoise,T, logMDDscore(:,j), Ny, y, yNaN, Nparticles, Nx, NsigmaX, Nw, Ex0, sqrtVx0, p, a0, a0V, sigaT, sigaDof, lambdaAlpha, lambdaBeta, Nsv, Eh0, Vh0, hvarT, hvarDof, sigmaT, sigmaDof, VSLstream)
           end select
           ! !$OMP CRITICAL
           ! CALL progressbarcomment( dble(iterCount) / dble(iterTotal) , itertimer, 'LLF grid')
           ! !$OMP END CRITICAL

        end do ! MDDgrid
        logMDD(:,thismodel,i) = sum(logMDDscore,2) / dble(MDDgrid)

        !$OMP CRITICAL
        CALL progressbarcomment( dble(iterCount) / dble(iterTotal) , itertimer, 'LLF grid')
        !$OMP END CRITICAL


     end do ! Nmodels


     ! VSLstreams
     errcode = vsldeletestream(VSLstream)
     ! print *,'releasing stream', i
     deallocate(dgp,y,yNaN)
     deallocate(logMDDscore)

     ! call hrulefill
     ! !$OMP CRITICAL
     ! CALL progressbarcomment(dble(gridCount) / dble(Ngrid), itertimer, 'LLF grid')
     ! !$OMP END CRITICAL

  end do ! Ngrid
  !$OMP END PARALLEL DO 

  ! CONSTRUCT FILE EXTENTSIONS
  fileXT = '.particles.MDDgrid' // trim(int2str(MDDgrid)) // '.' // trim(datalabel) // '.' // trim(modellabel)
  if (.not. doInflationNoise)  fileXT = trim(filext) // '.nonoise'
  fileXT = trim(filext) // '.Ngrid' // trim(int2str(Ngrid)) // '.Nparticles' // trim(int2str(Nparticles)) // '.dat'
  if (doTimeStamp) filext = '.' // timestampstr() //  filext

  ! WRITE SETTINGS
  CALL HRULEFILL
  filename = 'settings' // trim(adjustl(filext))
  OPEN (UNIT=4, FILE=filename, STATUS='REPLACE', ACTION='WRITE')
  WRITE(4,'(a20,a40)') 'TIME: ', timestampstr()
  WRITE(4,'(a20,a40)') 'Data: ', datalabel
  WRITE(4,'(a20,a40)') 'Model: ', modellabel
  WRITE(4,'(a60)') repeat('-',60)
  WRITE(4,'(a20,I40)') 'Nparticles: ', Nparticles
  WRITE(4,'(a20,I40)') 'Ngrid: ', Ngrid
  WRITE(4,'(a20,I40)') 'MDDgrid: ', MDDgrid
  WRITE(4,'(a20,I40)') 'p: ', p
  WRITE(4,'(a60)') repeat('-',60)
  WRITE(4,'(a20,F6.2)') 'a0: ', a0
  WRITE(4,'(a20,F6.2)') 'a0V: ', a0V
  WRITE(4,'(a20,F6.2)') 'siga0: ', siga0
  WRITE(4,'(a20,F6.2)') 'lambda0: ', lambda0
  WRITE(4,'(a20,F6.2)') 'lambda0V: ', lambda0V
  WRITE(4,'(a20,F6.2)') 'siglambda0: ', siglambda0
  do i=1,Ny
     WRITE(4,'(a20,I2,F6.2)') 'noisevol0(0): ', i, noisevol0(i)
  end do
  WRITE(4,'(a60)') repeat('-',60)
  if (doInflationNoise) THEN
     WRITE(4,'(a60)') 'With noise in inflation'
  ELSE
     WRITE(4,'(a60)') 'WITHOUT noise in inflation'
  END IF
  
  CLOSE(UNIT=4)
  CALL HRULEFILL

  do thismodel = 0,Nmodels

     filename = 'LOGMDD' // trim(int2str(thismodel)) // filext
     call savemat(logMDD(:,thismodel,:), filename)
     WRITE (*,*) 'STORED ' // filename

  end do

  DEALLOCATE (logMDD)
  DEALLOCATE (a00,lambda00,x00)

  call hrulefill
  WRITE(*,*) 'DONE. BYE, BYE. (' // trim(adjustl(filext)) // ')'
  call hrulefill

  STOP

CONTAINS

  ! -----------------------------------------------------------------
  SUBROUTINE getarguments(Nparticles,Ngrid,MDDgrid,doInflationNoise,T)

    IMPLICIT NONE 

    LOGICAL, INTENT(INOUT) :: doInflationNoise
    INTEGER, INTENT(INOUT)  :: Nparticles,Ngrid,MDDgrid,T


    INTEGER :: counter, dummy
    CHARACTER(len=32) :: arg

    counter = 0

    ! Nparticles
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') Nparticles
    END IF

    ! Ngrid
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') Ngrid
    END IF

    ! MDDgrid
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') MDDgrid
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

    ! T
    counter = counter + 1
    IF (command_argument_count() >= counter) THEN
       CALL get_command_argument(counter, arg) 
       READ(arg, '(i20)') T
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
