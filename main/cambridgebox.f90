MODULE cambridgebox

  use vslbox
  USE embox, only : savemat, savevec, int2str, mean, hrulefill
  use gibbsbox, only : drawNDXsysresample, igammaDraws
  use statespacebox, only : DLYAPsqrt
  use blaspack, only : pi, vechU,ivechU, vech, ivech, eye, symmetric, qrot, qrquery

  IMPLICIT NONE

  type :: ucsv(Ny,Nx,Nsv,T)
     integer, len :: Ny, Nx, Nsv, T
     double precision :: y(Ny,T), x(Nx,0:T), SV(Nsv,0:T), hinno(Nsv)
  end type ucsv

  type :: ucsvSI(Ny,Nx,Nsv,T)
     ! note: Nsurveys = Ny - 1
     integer, len :: Ny, Nx, Nsv, T
     double precision :: y(Ny,T), x(Nx,0:T), SV(Nsv,0:T), hinno(Nsv)
     double precision :: noise(Ny,T)
     double precision :: theta(0:T), lambda(0:T)
     double precision :: sigtheta, siglambda, noisevol(Ny)
  end type ucsvSI

contains

  ! @\newpage\subsection{generateUCSVdata}@
  SUBROUTINE generateUCSVdata(dgp, T, Ex0, Eh0, hinno, VSLstream)


    IMPLICIT NONE

    INTENT(INOUT) :: VSLstream
    INTENT(INOUT) :: dgp
    INTENT(IN)    :: T, Ex0, Eh0, hinno

    INTEGER :: T

    INTEGER :: ii,jj
    INTEGER, PARAMETER :: Nx = 2, Nsv = 2, Ny = 1

    INTEGER :: errcode

    type (ucsv(Ny,Nx,Nsv,T)) :: dgp

    double precision, dimension(Nsv,0:T) :: h,SV

    double precision :: Ex0(Nx), Eh0(Nsv), hinno(Nsv)

    DOUBLE PRECISION, DIMENSION(Ny,T) :: y
    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: z, x

    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    type (vsl_stream_state) :: VSLstream


    ! draw SV
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, (T + 1) * Nsv, h, 0.0d0, 1.0d0)
    ! note: h(:,0) draw currently unused
    h(:,0) = Eh0
    ! scale shocks
    ! forall(ii=1:Nsv,jj=1:T) h(ii,jj) = hinno(ii) * h(ii,jj)
    ! accumulate RW
    do jj=1,T
       h(:,jj) = h(:,jj-1) + hinno * h(:,jj)
    end do

    SV = exp(0.5d0 * h)


    ! draw x shocks
    ! note: z(:,0) currently unused
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, (T + 1) * Nx, z, 0.0d0, 1.0d0)
    ! scale shocks
    forall(ii=1:Nx,jj=1:T) z(ii,jj) = SV(ii,jj) * z(ii,jj)


    x(:,0) = Ex0

    ! accumulate RW
    ii = 1
    do jj=1,T
       x(ii,jj) = x(ii,jj-1) + z(ii,jj)
    end do


    ! gap
    x(2,:) = z(2,:)

    forall(jj=1:T) y(1,jj) = x(1,jj) + x(2,jj)


    ! copy into output structure
    dgp%y = y
    dgp%x = x
    dgp%SV = SV
    dgp%hinno = hinno


  END SUBROUTINE generateUCSVdata

  ! @\newpage\subsection{generateUCSVSIdata}@
  SUBROUTINE generateUCSVSIdata(dgp, doNoise, T, Ny, Ex0, Eh0, hinno, theta0, sigtheta, lambda0, siglambda, noisevol, VSLstream)


    IMPLICIT NONE

    INTENT(INOUT) :: VSLstream
    INTENT(INOUT) :: dgp
    INTENT(IN)    :: Ex0, Eh0, hinno

    INTEGER, PARAMETER :: Nx = 4, Nsv = 2
    INTEGER, INTENT(IN) :: T, Ny
    LOGICAL, INTENT(IN) :: doNoise

    double precision, intent(in), dimension(Ny) :: noisevol
    double precision, intent(in) :: theta0, sigtheta, lambda0, siglambda

    INTEGER :: ii,jj
    INTEGER :: Nsurveys

    INTEGER :: errcode

    type (ucsvSI(Ny,Nx,Nsv,T)) :: dgp

    double precision, dimension(Nsv,0:T) :: h,SV

    double precision :: Ex0(Nx), Eh0(Nsv), hinno(Nsv)

    DOUBLE PRECISION, DIMENSION(Ny,T) :: y, noise
    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: x
    DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: z

    DOUBLE PRECISION, DIMENSION(0:T) :: lambda, theta

    integer, parameter :: ndxtauRE = 1, ndxgapRE = 2, ndxtauSI = 3, ndxgapSI = 4

    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    type (vsl_stream_state) :: VSLstream

    logical :: OK
    double precision, parameter :: maxSVh = 7.0d0 ! maxSVh = 8 is like SV <  55 (maxSVh=10 is like SV < 100)

    ! default values
    if (sigtheta .gt. 0.0d0) then
       theta  = drawtruncRWvec(T, theta0, sigtheta, -1.0d0, 1.0d0, VSLstream)
    else
       theta = theta0
    end if
    if (siglambda .gt. 0.0d0) then
       lambda = drawtruncRWvec(T, lambda0, siglambda, 0.0d0, 1.0d0, VSLstream)
    else
       lambda = lambda0
    end if


    ! prelim
    Nsurveys = Ny - 1

    ! draw SV
    OK = .false.
    do while (.not. OK)
       errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, (T + 1) * Nsv, h, 0.0d0, 1.0d0)
       ! note: h(:,0) draw currently unused
       h(:,0) = Eh0
       ! accumulate RW
       do jj=1,T
          h(:,jj) = h(:,jj-1) + hinno * h(:,jj)
       end do
       OK = all(h < maxSVh)
    end do
    SV = exp(0.5d0 * h)


    ! draw x shocks
    ! note: z(:,0) currently unused
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, (T + 1) * Nsv, z, 0.0d0, 1.0d0)
    ! scale shocks
    forall(ii=1:Nsv,jj=1:T) z(ii,jj) = SV(ii,jj) * z(ii,jj)
    ! note: z(:,0) currently not used

    x(:,0) = Ex0

    ! accumulate trend and gap
    do jj=1,T
       x(ndxtauRE,jj) = x(ndxtauRE,jj-1) + z(ndxtauRE,jj)
       x(ndxgapRE,jj) = theta(jj) * x(ndxgapRE,jj-1) + z(ndxgapRE,jj)
    end do

    ! SI 
    do jj=1,T
       x(ndxtauSI,jj) = (1.0d0 - lambda(jj)) * x(ndxtauRE,jj) + lambda(jj) * x(ndxtauSI,jj-1)
       x(ndxgapSI,jj) = (1.0d0 - lambda(jj)) * x(ndxgapRE,jj) + lambda(jj) * theta(jj) * x(ndxgapSI,jj-1) 
    end do


    ! inflation
    forall(jj=1:T) y(1,jj) = x(ndxtauRE,jj) + x(ndxgapRE,jj)
    ! forecasts
    forall(ii=1:Nsurveys,jj=1:T) y(1+ii,jj) = x(ndxtauSI,jj) + (theta(jj) ** ii) * x(ndxgapSI,jj)


    ! add noise
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Ny, noise, 0.0d0, 1.0d0)
    forall(jj=1:T,ii=1:Ny) noise(ii,jj) = noisevol(ii) * noise(ii,jj)
    if (.not. donoise) noise(1,:)   = 0.0d0 
    y = y + noise

    ! copy into output structure
    dgp%y      = y
    dgp%x      = x
    dgp%noise  = noise
    dgp%SV     = SV
    dgp%hinno  = hinno
    dgp%lambda = lambda
    dgp%theta  = theta

    dgp%sigtheta   = sigtheta
    dgp%siglambda  = siglambda
    dgp%noisevol   = noisevol
    if (.not. donoise) dgp%noisevol(1)  = 0.0d0 ! note: cannot reset noisevol input parameter since intent(in)


  END SUBROUTINE generateUCSVSIdata

  ! @\newpage\subsection{generateUCSVSIdataGivenStates}@
  SUBROUTINE generateUCSVSIdataGivenStates(dgp, doNoise, T, Ny, Ex0, theta, lambda, SV, noisevol, VSLstream)


    IMPLICIT NONE

    INTEGER, PARAMETER :: Nx = 4, Nsv = 2
    INTEGER, INTENT(IN) :: T, Ny
    LOGICAL, INTENT(IN) :: doNoise

    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    type (vsl_stream_state), intent(inout) :: VSLstream

    type (ucsvSI(Ny,Nx,Nsv,T)), intent(inout) :: dgp
    double precision, intent(in) :: Ex0(Nx)

    double precision, intent(in), dimension(Nsv,0:T) :: SV
    double precision, intent(in), dimension(0:T) :: theta, lambda

    double precision, intent(in), dimension(Ny) :: noisevol

    INTEGER :: ii,jj
    INTEGER :: Nsurveys

    INTEGER :: errcode

    DOUBLE PRECISION, DIMENSION(Ny,T) :: y, noise
    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: x
    DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: z

    integer, parameter :: ndxtauRE = 1, ndxgapRE = 2, ndxtauSI = 3, ndxgapSI = 4


    ! prelim
    Nsurveys = Ny - 1

    ! draw x shocks
    ! note: z(:,0) currently unused
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, (T + 1) * Nsv, z, 0.0d0, 1.0d0)
    ! scale shocks
    forall(ii=1:Nsv,jj=1:T) z(ii,jj) = SV(ii,jj) * z(ii,jj)


    x(:,0) = Ex0

    ! accumulate trend and gap
    do jj=1,T
       x(ndxtauRE,jj) = x(ndxtauRE,jj-1) + z(ndxtauRE,jj)
       x(ndxgapRE,jj) = theta(jj) * x(ndxgapRE,jj-1) + z(ndxgapRE,jj)
    end do

    ! SI 
    do jj=1,T
       x(ndxtauSI,jj) = (1.0d0 - lambda(jj)) * x(ndxtauRE,jj) + lambda(jj) * x(ndxtauSI,jj-1)
       x(ndxgapSI,jj) = (1.0d0 - lambda(jj)) * x(ndxgapRE,jj) + lambda(jj) * theta(jj) * x(ndxgapSI,jj-1) 
    end do


    ! inflation
    forall(jj=1:T) y(1,jj) = x(ndxtauRE,jj) + x(ndxgapRE,jj)
    ! forecasts
    forall(ii=1:Nsurveys,jj=1:T) y(1+ii,jj) = x(ndxtauSI,jj) + (theta(jj) ** ii) * x(ndxgapSI,jj)


    ! add noise
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Ny, noise, 0.0d0, 1.0d0)
    forall(jj=1:T,ii=1:Ny) noise(ii,jj) = noisevol(ii) * noise(ii,jj)
    if (.not. donoise) noise(1,:)   = 0.0d0 
    y = y + noise

    ! copy into output structure
    dgp%y      = y
    dgp%x      = x
    dgp%noise  = noise
    dgp%SV         = SV
    dgp%lambda     = lambda
    dgp%theta      = theta

    dgp%hinno      = -999.0d0
    dgp%sigtheta   = -999.0d0
    dgp%siglambda  = -999.0d0

    dgp%noisevol   = noisevol
    if (.not. donoise) dgp%noisevol(1)  = 0.0d0 ! note: cannot reset noisevol since intent(in)


  END SUBROUTINE generateUCSVSIdataGivenStates

  ! @\newpage\subsection{generateCambridgeDataLambdaKink}@
  SUBROUTINE generateCambridgeDataLambdaKink(dgp, doNoise, T, Ny, Ex0, Eh0, hinno, theta0, sigtheta, lambdaStart, lambdaEnd, lambdaKinkDate, noisevol, VSLstream)

    ! implements only two values for lambda, switching at middle of sample

    IMPLICIT NONE

    INTENT(INOUT) :: VSLstream
    INTENT(INOUT) :: dgp
    INTENT(IN)    :: Ex0, Eh0, hinno

    INTEGER, PARAMETER :: Nx = 4, Nsv = 2
    INTEGER, INTENT(IN) :: T, Ny
    LOGICAL, INTENT(IN) :: doNoise

    double precision, intent(in), dimension(Ny) :: noisevol
    double precision, intent(in) :: theta0, sigtheta
    double precision, intent(in) :: lambdaStart, lambdaEnd
    integer, intent(in) :: lambdaKinkDate

    INTEGER :: ii,jj
    INTEGER :: Nsurveys

    INTEGER :: errcode

    type (ucsvSI(Ny,Nx,Nsv,T)) :: dgp

    double precision, dimension(Nsv,0:T) :: h,SV

    double precision :: Ex0(Nx), Eh0(Nsv), hinno(Nsv)

    DOUBLE PRECISION, DIMENSION(Ny,T) :: y, noise
    DOUBLE PRECISION, DIMENSION(Nx,0:T) :: x
    DOUBLE PRECISION, DIMENSION(Nsv,0:T) :: z

    DOUBLE PRECISION, DIMENSION(0:T) :: lambda, theta

    integer, parameter :: ndxtauRE = 1, ndxgapRE = 2, ndxtauSI = 3, ndxgapSI = 4

    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0
    type (vsl_stream_state) :: VSLstream

    logical :: OK
    double precision, parameter :: maxSVh = 7.0d0 ! maxSVh = 8 is like SV <  55 (maxSVh=10 is like SV < 100)

    ! default values
    if (sigtheta .gt. 0.0d0) then
       theta  = drawtruncRWvec(T, theta0, sigtheta, -1.0d0, 1.0d0, VSLstream)
    else
       theta = theta0
    end if

    lambda(0:lambdaKinkDate) = lambdaStart
    lambda(1+lambdaKinkDate:T) = lambdaEnd

    ! prelim
    Nsurveys = Ny - 1

    ! draw SV
    OK = .false.
    do while (.not. OK)
       errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, (T + 1) * Nsv, h, 0.0d0, 1.0d0)
       ! note: h(:,0) draw currently unused
       h(:,0) = Eh0
       ! accumulate RW
       do jj=1,T
          h(:,jj) = h(:,jj-1) + hinno * h(:,jj)
       end do
       OK = all(h < maxSVh)
    end do
    SV = exp(0.5d0 * h)


    ! draw x shocks
    ! note: z(:,0) currently unused
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, (T + 1) * Nsv, z, 0.0d0, 1.0d0)
    ! scale shocks
    forall(ii=1:Nsv,jj=1:T) z(ii,jj) = SV(ii,jj) * z(ii,jj)


    x(:,0) = Ex0

    ! accumulate trend and gap
    do jj=1,T
       x(ndxtauRE,jj) = x(ndxtauRE,jj-1) + z(ndxtauRE,jj)
       x(ndxgapRE,jj) = theta(jj) * x(ndxgapRE,jj-1) + z(ndxgapRE,jj)
    end do

    ! SI 
    do jj=1,T
       x(ndxtauSI,jj) = (1.0d0 - lambda(jj)) * x(ndxtauRE,jj) + lambda(jj) * x(ndxtauSI,jj-1)
       x(ndxgapSI,jj) = (1.0d0 - lambda(jj)) * x(ndxgapRE,jj) + lambda(jj) * theta(jj) * x(ndxgapSI,jj-1) 
    end do


    ! inflation
    forall(jj=1:T) y(1,jj) = x(ndxtauRE,jj) + x(ndxgapRE,jj)
    ! forecasts
    forall(ii=1:Nsurveys,jj=1:T) y(1+ii,jj) = x(ndxtauSI,jj) + (theta(jj) ** ii) * x(ndxgapSI,jj)


    ! add noise
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T * Ny, noise, 0.0d0, 1.0d0)
    forall(jj=1:T,ii=1:Ny) noise(ii,jj) = noisevol(ii) * noise(ii,jj)
    if (.not. donoise) noise(1,:)   = 0.0d0 
    y = y + noise

    ! copy into output structure
    dgp%y      = y
    dgp%x      = x
    dgp%noise  = noise
    dgp%SV     = SV
    dgp%hinno  = hinno
    dgp%lambda = lambda
    dgp%theta  = theta

    dgp%sigtheta   = sigtheta
    dgp%siglambda  = 0.0d0
    dgp%noisevol   = noisevol
    if (.not. donoise) dgp%noisevol(1)  = 0.0d0 ! note: cannot reset noisevol since intent(in)


  END SUBROUTINE generateCambridgeDataLambdaKink

  ! @\newpage\subsection{drawtruncRW}@
  FUNCTION drawtruncRW(N,T, x0, sig, lb, ub, VSLstream) RESULT(x)

    IMPLICIT NONE

    INTENT(IN) :: N, T, x0, sig, lb, ub
    INTENT(INOUT) :: VSLstream
    INTEGER :: N,T
    DOUBLE PRECISION :: sig, lb, ub
    DOUBLE PRECISION :: x0, x(N,0:T), dx(N)

    integer :: ii

    ! VSL
    TYPE (VSL_STREAM_STATE) :: VSLstream
    INTEGER :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    x(:,0) = x0

    do ii=1,T

       errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, N, dx, 0.0d0, sig)
       x(:,ii) = x(:,ii-1) + dx

       do while (ANY(x(:,ii) < lb) .OR. ANY(x(:,ii) > ub))
          errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, N, dx, 0.0d0, sig) ! a bit brute force to *always* redraw N random variables ..
          where (x(:,ii) < lb) x(:,ii) = x(:,ii-1) + dx
          where (x(:,ii) > ub) x(:,ii) = x(:,ii-1) + dx
       end do

    end do


  END FUNCTION drawtruncRW

  ! @\newpage\subsection{drawtruncRWvec}@
  FUNCTION drawtruncRWvec(T, x0, sig, lb, ub, VSLstream) RESULT(x)
    ! draw entire RW vector and then accept/reject trajectory

    IMPLICIT NONE

    INTENT(IN) :: T, x0, sig, lb, ub
    INTENT(INOUT) :: VSLstream
    INTEGER :: T
    DOUBLE PRECISION :: sig, lb, ub
    DOUBLE PRECISION :: x0, x(0:T), dx(T)

    integer :: ii
    logical :: OK

    ! VSL
    TYPE (VSL_STREAM_STATE) :: VSLstream
    INTEGER :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    x(0) = x0

    OK = .false.
    do while (.not. ok)   
       errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, T, dx, 0.0d0, sig)
       do ii=1,T
          x(ii) = x(ii-1) + dx(ii)
       end do

       OK = (ALL((x .ge. lb) .AND. (x .le. ub)))
    end do

  END FUNCTION drawtruncRWvec

  ! @\newpage\subsection{drawtruncnorm}@
  FUNCTION drawtruncnorm(N, x0, sig, lb, ub, VSLstream) RESULT(x)

    IMPLICIT NONE

    INTENT(IN) :: N, x0, sig, lb, ub
    INTENT(INOUT) :: VSLstream
    INTEGER :: N
    DOUBLE PRECISION :: sig, lb, ub
    DOUBLE PRECISION, DIMENSION(N) :: x0, x, dx


    ! VSL
    TYPE (VSL_STREAM_STATE) :: VSLstream
    INTEGER :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, N, dx, 0.0d0, sig)
    x = x0 + dx
    do while (ANY(x < lb) .OR. ANY(x > ub))
       errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, N, dx, 0.0d0, sig) ! a bit brute force to *always* redraw N random variables ..
       where (x < lb) 
          x = x0 + dx
       elsewhere (x > ub) 
          x = x0 + dx
       end where
    end do

  END FUNCTION drawtruncnorm

  ! @\newpage\subsection{drawdeltatruncnorm}@
  FUNCTION drawdeltatruncnorm(N, x0, sig, lb, ub, VSLstream) RESULT(delta)

    IMPLICIT NONE

    INTENT(IN) :: N, x0, sig, lb, ub
    INTENT(INOUT) :: VSLstream
    INTEGER :: N
    DOUBLE PRECISION :: lb, ub
    DOUBLE PRECISION, DIMENSION(N) :: x0, x, dx, sig, delta


    ! VSL
    TYPE (VSL_STREAM_STATE) :: VSLstream
    INTEGER :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, N, dx, 0.0d0, 1.0d0)
    dx = dx * sig
    x     = x0 + dx
    delta = dx
    do while (ANY(x < lb) .OR. ANY(x > ub))
       errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, N, dx, 0.0d0, 1.0d0) ! a bit brute force to *always* redraw N random variables ..
       dx = dx * sig
       where (x < lb) 
          x     = x0 + dx
          delta = dx
       elsewhere (x > ub) 
          x = x0 + dx
          delta = dx
       end where
    end do

  END FUNCTION drawdeltatruncnorm


  ! @\newpage\subsection{drawdeltatruncnorm1}@
  FUNCTION drawdeltatruncnorm1(N, x0, sig, lb, ub, VSLstream) RESULT(delta)

    IMPLICIT NONE

    INTENT(IN) :: N, x0, sig, lb, ub
    INTENT(INOUT) :: VSLstream
    INTEGER :: N
    DOUBLE PRECISION :: lb, ub
    DOUBLE PRECISION, DIMENSION(N) :: x0, x, dx, delta
    DOUBLE PRECISION :: sig


    ! VSL
    TYPE (VSL_STREAM_STATE) :: VSLstream
    INTEGER :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, N, dx, 0.0d0, 1.0d0)
    dx = dx * sig
    x     = x0 + dx
    delta = dx
    do while (ANY(x < lb) .OR. ANY(x > ub))
       errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, N, dx, 0.0d0, 1.0d0) ! a bit brute force to *always* redraw N random variables ..
       dx = dx * sig
       where (x < lb) 
          x     = x0 + dx
          delta = dx
       elsewhere (x > ub) 
          x = x0 + dx
          delta = dx
       end where
    end do

  END FUNCTION drawdeltatruncnorm1

  ! @\newpage\subsection{drawBeta}@
  FUNCTION drawBeta(alpha, beta, VSLstream) RESULT(draw) 

    ! single draw from beta distribution (simulated with normals)

    IMPLICIT NONE

    INTENT(IN) :: alpha, beta
    INTENT(INOUT) :: VSLstream

    INTEGER  :: alpha, beta

    DOUBLE PRECISION :: draw, achi2, bchi2

    DOUBLE PRECISION :: z(2 * (alpha + beta))


    ! VSL
    TYPE (VSL_STREAM_STATE) :: VSLstream
    INTEGER :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0


    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, 2 * (alpha + beta), z, 0.0d0, 1.0d0) 
    achi2 = sum(z(1:2*alpha) ** 2)
    bchi2 = sum(z(2*alpha+1:2*(alpha+beta)) ** 2)

    draw = achi2 / (achi2 + bchi2)

  END FUNCTION drawBeta

  ! @\newpage\subsection{drawBetas}@
  FUNCTION drawBetas(alpha, beta, Ndraws, VSLstream) RESULT(draws) ! single draw so far

    ! multiple draws from beta distribution (with separate parameters for each draw)
    ! jsut a wrapper to VSL

    IMPLICIT NONE


    INTEGER, INTENT(IN) :: Ndraws

    DOUBLE PRECISION, INTENT(IN), DIMENSION(Ndraws)  :: alpha, beta

    DOUBLE PRECISION, DIMENSION(Ndraws) :: draws

    ! VSL
    TYPE (VSL_STREAM_STATE), INTENT(INOUT) :: VSLstream
    INTEGER :: errcode
    INTEGER, PARAMETER :: VSLmethodGaussian = 0, VSLmethodUniform = 0, VSLmethodBeta = 0 ! VSL_RNG_METHOD_BETA_CJA is 0
    INTEGER :: n

    do n=1,Ndraws
       errcode = vdrngBeta(VSL_RNG_METHOD_BETA_CJA, VSLstream, 1, draws(n), alpha(n), beta(n), 0.0d0, 1.0d0) 
    end do

  END FUNCTION drawBetas

  ! @\newpage\subsection{simXhatSI}@
  FUNCTION simXhatSI(Nx, T, Ndraws, x0SI, xRE, lambda, p, infgapcompanion) RESULT(xSI)

    IMPLICIT NONE

    INTENT(IN) :: Nx, T, Ndraws, x0SI, xRE, lambda, p, infgapcompanion

    INTEGER :: Nx,T,Ndraws,p
    DOUBLE PRECISION :: x0SI(Nx), xRE(Nx,0:T), lambda(Ndraws,0:T), xSI(Nx,Ndraws,0:T), infgapcompanion(p,p)
    DOUBLE PRECISION, DIMENSION(Nx,Nx) :: A
    INTEGER :: ndxTrend = 1, ndxGapStart = 2, ndxGapStop
    INTEGER :: ii, hh

    ndxGapStop = ndxGapStart + p - 1


    A = 0.0d0
    A(ndxTrend,ndxTrend) = 1.0d0
    A(ndxGapStart:ndxGapStop,ndxGapStart:ndxGapStop) = infgapcompanion

    forall (ii=1:Ndraws) xSI(:,ii,0) = x0SI
    forall (ii=1:Ndraws,hh=1:T) xSI(:,ii,hh) = xRE(:,hh)
    do hh=1,T
       do ii=1,Ndraws
          call dgemv('n',Nx,Nx,lambda(ii,hh),A,Nx,xSI(:,ii,hh-1),1,(1.0d0 - lambda(ii,hh)),xSI(:,ii,hh),1)
       end do
    end do

    ! ! xhatSI trend component 
    ! forall (ii=1:Ndraws) xSI(ndxTrend,ii,0) = x0SI(ndxTrend)
    ! do hh=1,T
    !    forall (ii=1:Ndraws) xSI(ndxTrend,ii,hh) = lambda(ii,hh) * xSI(ndxTrend,ii,hh-1) + (1.0d0 - lambda(ii,hh)) * xRE(ndxTrend,hh)
    ! end do

    ! ! simulate xhatSI gap component
    ! forall (ii=1:Ndraws) xSI(ndxGapStart:ndxGapStop,ii,0)  = x0SI(ndxGapStart:ndxGapStop)
    ! if ((p == 1) .AND. (infgapcompanion(1,1) == 0)) then
    !    xSI(ndxGapStart:ndxGapStop,:,1:T)  = 0.0d0
    ! else
    !    do hh=1,T
    !       do ii=1,Ndraws
    !          call dgemv('n',p,p,lambda(ii,hh),infgapcompanion,p,xSI(ndxGapStart:ndxGapStop,ii,hh-1),1,(1.0d0 - lambda(ii,hh)),xSI(ndxGapStart:ndxGapStop,ii,hh),1)
    !       end do
    !    end do
    ! end if



  END FUNCTION simXhatSI


  ! @\newpage\subsection{MDDthetaTVPlambdaTVP}@
  SUBROUTINE MDDthetaTVPlambdaTVP(doInflationNoise,T, logMDD, Ny, y, yNaN, Nparticles, Nx, NsigmaX, Nw, Ex0, sqrtVx00, p, a0, a0V,  sigaT, sigaDof,  lambda0, lambda0V, siglambdaT, siglambdaDof, Nsv, Eh0, Vh0, hvarT, hvarDof, sigmaT, sigmaDof, VSLstream)


    IMPLICIT NONE

    INTENT(INOUT) :: VSLstream ! , timer
    INTENT(IN)    :: T,Ny,y,yNaN, Nx, NsigmaX, Nw, Ex0, sqrtVx00, Nsv,Eh0,Vh0, Nparticles, p, lambda0, lambda0V, a0, a0V, sigaT, sigaDof, siglambdaT, siglambdaDof, hvarT, hvarDof, sigmaT, sigmaDof
    INTENT(OUT)   :: logMDD

    INTEGER :: J, I, K, T, Nparticles, Nx, Ny, Nsv, NsigmaX, p, Nw, Nsurveys




    logical, intent(in) :: doInflationNoise

    ! double precision, parameter :: minParticleWeight = 1.0d-12

    DOUBLE PRECISION, DIMENSION(Ny,T) :: y
    LOGICAL, DIMENSION(Ny,T) :: yNaN
    DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: PARTICLEweights
    DOUBLE PRECISION, DIMENSION(T)   :: logMDD


    ! APF llf correction
    DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
    DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax
    ! ! scale parameters
    DOUBLE PRECISION, DIMENSION(Nparticles)     :: DRAWsiga
    DOUBLE PRECISION, DIMENSION(Nparticles)     :: DRAWsiglambda
    ! DOUBLE PRECISION, DIMENSION(Nparticles,Nsv) :: DRAWhInno
    DOUBLE PRECISION, DIMENSION(Nsv,Nparticles)   :: hInno ! helper variable, note the transpose
    DOUBLE PRECISION, DIMENSION(Nparticles,Ny)  :: DRAWsigma

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

    logical, parameter :: doSecondResamplingStep = .false.
    integer, parameter :: dofAPFfattail = 0.0d0

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
             write (*,*) 'DLYAP error (ygap0sqrtvariance -- init particlefilter MDDthetaTVPlambdaTVP)', errcode
             call savemat(gaptransition, 'gaptransition.debug')
             call savemat(gapshock0loadings, 'gapshock0loadings.debug')
             call savemat(ygap0sqrtvariance, 'ygap0sqrtvariance.debug')
             ! stop 1
          else
             sqrtVx0(ndxgap,ndxgap) = ygap0sqrtvariance
          end if
       end if
       vecSqrtSigmaX(:,k)     = vechU(sqrtVx0,Nx)

    END DO


    PARTICLEweights  = 1.0d0 / dble(Nparticles)
    logMDD         = 0.0d0

    ! uniform draws for systematic resampling
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


    ! workspace query for qr decomposition
    qrR     = 0.0d0
    qrlwork = qrquery(qrR)

    DO j=1,T

       Nynonan = count(.not. yNaN(:,j))
       ! print *, 'plefilter for t=', j, ' with ', Nynonan, 'obs' 
       ! ------------------------------------------------------------------------------------------------------------------------------
       ! BEGIN: APF RESAMPLE STEP
       ! ------------------------------------------------------------------------------------------------------------------------------

       ! Sigma
       i=1
       if (doInflationNoise) then
          call igammaDraws(DRAWsigma(:,i), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream)
       else
          DRAWsigma(:,i) = 0.0d0
          PREVsigmaT(:,i)  = 0.0d0
       end if
       DO i=2,Ny 
          call igammaDraws(DRAWsigma(:,i), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream) 
       END DO

       DO k = 1,Nparticles

          ! 2) Fill Particles into state space

          ! update A 

          A(ndxGapRE,ndxGapRE)           = adrift(k)
          A(ndxtrendSI,ndxtrendRE)       = 1 - lambda(k)
          A(ndxtrendSI,ndxtrendSI)       = lambda(k)
          A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - lambda(k)) * adrift(k)
          A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = lambda(k) * adrift(k)

          ! update sqrtR
          sqrtR = 0.0d0
          forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i)) 


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
             FORALL(k=1:Nparticles) shufflevec(k) = DRAWsigma(ndx(k),i)
             DRAWsigma(:,i) = shufflevec
          END DO

       end if

       ! ------------------------------------------------------------------------------------------------------------------------------
       ! END: APF STEP
       ! ------------------------------------------------------------------------------------------------------------------------------

       ! 0) draw scale parameters
       ! siga
       call igammaDraws(DRAWsiga, Nparticles, PREVsigaT, PREVsigaDof, VSLstream)
       DRAWsiga = sqrt(DRAWsiga)
       ! siglambda
       call igammaDraws(DRAWsiglambda, Nparticles, PREVsiglambdaT, PREVsiglambdaDof, VSLstream)
       DRAWsiglambda = sqrt(DRAWsiglambda)
       ! hInno
       DO i=1,Nsv 
          call igammaDraws(hinno(i,:), Nparticles, PREVhvarT(:,i), PREVhvarDof(i), VSLstream) 
       END DO
       hinno = sqrt(hinno)
       ! forall (i=1:Nsv,k=1:Nparticles) hInno(i,k) = DRAWhInno(k,i) ! helper variable, provides better aligned access to the j data


       ! 1) Draw Particles 
       errcode     = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, hDELTA, 0.0d0, 1.0d0)
       hDELTA      = hDELTA * hInno
       h           = h + hDELTA
       SVol        = exp(h * 0.5d0)


       lambdaPREV  = lambda
       lambdaDELTA = drawdeltatruncnorm(Nparticles, lambdaPREV, DRAWsiglambda, 0.0d0, 1.0d0, VSLstream)
       lambda      = lambdaPREV + lambdaDELTA

       aPREV      = adrift
       aDELTA     = drawdeltatruncnorm(Nparticles, aPREV, DRAWsiga, -1.0d0, 1.0d0, VSLstream)
       adrift     = aPREV + aDELTA

       ! ------------------------------------------------------------------------------------------------------------------------------
       ! BEGIN: MAIN PARTICLE STEP
       ! ------------------------------------------------------------------------------------------------------------------------------



       DO k = 1,Nparticles


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
          forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i)) 


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

          ! ------------------------------------------------------------------------
          ! DONE: SQRT KALMAN
          ! ------------------------------------------------------------------------



       END DO ! k particles

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



    END DO ! j=1,T

  END SUBROUTINE MDDthetaTVPlambdaTVP

  ! @\newpage\subsection{MDDthetaCONSTlambdaTVP}@
  SUBROUTINE MDDthetaCONSTlambdaTVP(doInflationNoise, T, logMDD, Ny, y, yNaN, Nparticles, Nxx, NsigmaXX, Nx, Nw, Ex0, sqrtVx00, p, a0, aV0, lambda0, lambda0V, siglambdaT, siglambdaDof, Nsv, Eh0, Vh0, hvarT, hvarDof, sigmaT, sigmaDof, VSLstream)


    IMPLICIT NONE

    INTENT(INOUT) :: VSLstream ! , timer
    INTENT(IN)    :: T,Ny,y,yNaN, Nx, Nw, Ex0, sqrtVx00, Nsv,Eh0,Vh0, Nparticles, p, lambda0, lambda0V, a0, aV0, siglambdaT, siglambdaDof, hvarT, hvarDof, sigmaT, sigmaDof
    INTENT(IN) :: Nxx, NsigmaXX
    INTENT(OUT) :: logMDD
    INTEGER :: Nxx, NsigmaXX

    INTEGER :: J, I, K, T, Nparticles, Nx, Ny, Nsv, p, Nw, Nsurveys




    logical, intent(in) :: doInflationNoise

    ! double precision, parameter :: minParticleWeight = 1.0d-12

    DOUBLE PRECISION, DIMENSION(Ny,T) :: y
    LOGICAL, DIMENSION(Ny,T) :: yNaN
    DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: PARTICLEweights
    DOUBLE PRECISION, DIMENSION(T)   :: logMDD


    ! APF llf correction
    DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
    DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax
    ! ! scale parameters
    DOUBLE PRECISION, DIMENSION(Nparticles)     :: DRAWsiglambda
    ! DOUBLE PRECISION, DIMENSION(Nparticles,Nsv) :: DRAWhInno
    DOUBLE PRECISION, DIMENSION(Nsv,Nparticles)   :: hInno ! helper variable, note the transpose
    DOUBLE PRECISION, DIMENSION(Nparticles,Ny)  :: DRAWsigma

    ! particles
    DOUBLE PRECISION :: xxposterior(Nxx,Nparticles),  h(Nsv, Nparticles), SVol(Nsv,Nparticles), llf(Nparticles)
    DOUBLE PRECISION :: kernelweights(Nparticles) ! for unnormalized weights
    DOUBLE PRECISION  :: kernelsum, loglikemax
    INTEGER :: ndx(Nparticles)
    DOUBLE PRECISION :: shufflevec(Nparticles)

    ! state space objects
    DOUBLE PRECISION :: xxprior(Nxx), logdetSigmaY
    DOUBLE PRECISION :: lambda(Nparticles), lambdaPREV(Nparticles), lambda0, lambda0V
    DOUBLE PRECISION :: adrift(Nparticles), zdraw(Nparticles),  a0, aV0 ! note: no more "drift" in a but keeping name for sake of comparability with other code
    DOUBLE PRECISION :: zstats(Nparticles)
    DOUBLE PRECISION, PARAMETER :: rejectioncritvalue = 3, unitcircle = 1.0d0 - 1.0d-5
    INTEGER :: counter
    INTEGER, PARAMETER :: maxcount = 10000

    ! scale parameters
    DOUBLE PRECISION :: siglambdaT, hvarT(Nsv), sigmaT(Ny)
    INTEGER :: siglambdaDof, hvarDof(Nsv), sigmaDof(Ny)

    DOUBLE PRECISION :: PREVsiglambdaT(Nparticles), PREVhvarT(Nparticles,Nsv), PREVsigmaT(Nparticles,Ny)
    INTEGER :: PREVsiglambdaDof, PREVhvarDof(Nsv), PREVsigmaDof(Ny)
    DOUBLE PRECISION :: lambdaDELTA(Nparticles), hDELTA(Nsv,Nparticles), resid(Ny,Nparticles)
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

    logical, parameter :: doSecondResamplingStep = .false.
    integer, parameter :: dofAPFfattail = 0.0d0

    ! CHARACTER (LEN=100) :: filename

    Nsurveys = Ny - 1
    minSVh   = log(0.001d0 ** 2)

    ! state-vector indices
    ndxTrendRE  = 1
    ndxGapRE    = 2
    ndxTrendSI  = 2+p
    ndxGapSI    = ndxTrendSI + 1

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

    PREVsiglambdaT   = siglambdaT
    PREVsiglambdaDof = siglambdaDof

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
    errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, lambdaPREV, lambda0, sqrt(lambda0V))
    lambda = lambdaPREV
    do while (ANY(lambda < 0.0d0) .OR. ANY(lambda > 1.0d0))
       errcode = vdrnggaussian(VSLmethodGaussian, VSLstream, Nparticles, lambdaPREV, lambda0, sqrt(lambda0V)) ! a bit brute force to *always* redraw N random variables ..
       where (lambda < 0.0d0)   lambda = lambdaPREV
       where (lambda > 1.0d0)   lambda = lambdaPREV
    end do

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
             write (*,*) 'DLYAP error (ygap0sqrtvariance -- init particlefilter MDDthetaCONSTlambdaTVP)', errcode
             call savemat(gaptransition, 'gaptransition.debug')
             call savemat(gapshock0loadings, 'gapshock0loadings.debug')
             call savemat(ygap0sqrtvariance, 'ygap0sqrtvariance.debug')
             ! stop 1
          else 
             sqrtVxx0(ndxgap,ndxgap) = ygap0sqrtvariance
          end if
       end if
       vecSqrtSigmaXX(:,k)     = vechU(sqrtVxx0,Nxx)

    END DO


    PARTICLEweights  = 1.0d0 / dble(Nparticles)
    logMDD           = 0.0d0

    ! uniform draws for systematic resampling
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


    ! workspace query for qr decomposition
    qrR     = 0.0d0
    qrlwork = qrquery(qrR)

    DO j=1,T

       Nynonan = count(.not. yNaN(:,j))
       ! print *, 'plefilter for t=', j, ' with ', Nynonan, 'obs' 
       ! ------------------------------------------------------------------------------------------------------------------------------
       ! BEGIN: APF RESAMPLE STEP
       ! ------------------------------------------------------------------------------------------------------------------------------

       ! Sigma
       i=1
       if (doInflationNoise) then
          call igammaDraws(DRAWsigma(:,i), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream)
       else
          DRAWsigma(:,i) = 0.0d0
          PREVsigmaT(:,i)  = 0.0d0
       end if
       DO i=2,Ny 
          call igammaDraws(DRAWsigma(:,i), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream) 
       END DO

       if (j > 1) then

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
          ! print *, 'adrift rejection draws, t=', j, counter

       end if


       DO k = 1,Nparticles

          ! 2) Fill Particles into state space

          ! update A 

          A(ndxGapRE,ndxGapRE)           = adrift(k)
          A(ndxtrendSI,ndxtrendRE)       = 1 - lambda(k)
          A(ndxtrendSI,ndxtrendSI)       = lambda(k)
          A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - lambda(k)) * adrift(k)
          A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = lambda(k) * adrift(k)

          ! update sqrtR
          sqrtR = 0.0d0
          forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i)) 


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

          FORALL(k=1:Nparticles) shufflevec(k) = adrift(ndx(k))
          adrift = shufflevec
          FORALL(k=1:Nparticles) shufflevec(k) = PREVa(ndx(k))
          PREVa = shufflevec
          FORALL(k=1:Nparticles) shufflevec(k) = PREVaV(ndx(k))
          PREVaV = shufflevec

          ! reshuffle sufficient statistics for scale parameters
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
             FORALL(k=1:Nparticles) shufflevec(k) = DRAWsigma(ndx(k),i)
             DRAWsigma(:,i) = shufflevec
          END DO

       end if

       ! ------------------------------------------------------------------------------------------------------------------------------
       ! END: APF STEP
       ! ------------------------------------------------------------------------------------------------------------------------------

       ! 0) draw scale parameters
       ! siglambda
       call igammaDraws(DRAWsiglambda, Nparticles, PREVsiglambdaT, PREVsiglambdaDof, VSLstream)
       DRAWsiglambda = sqrt(DRAWsiglambda)
       ! hInno
       DO i=1,Nsv 
          call igammaDraws(hinno(i,:), Nparticles, PREVhvarT(:,i), PREVhvarDof(i), VSLstream) 
       END DO
       hinno = sqrt(hinno)
       ! forall (i=1:Nsv,k=1:Nparticles) hInno(i,k) = DRAWhInno(k,i) ! helper variable, provides better aligned access to the j data


       ! 1) Draw Particles 
       errcode     = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, hDELTA, 0.0d0, 1.0d0)
       hDELTA      = hDELTA * hInno
       h           = h + hDELTA
       SVol        = exp(h * 0.5d0)


       lambdaPREV  = lambda
       lambdaDELTA = drawdeltatruncnorm(Nparticles, lambdaPREV, DRAWsiglambda, 0.0d0, 1.0d0, VSLstream)
       lambda      = lambdaPREV + lambdaDELTA

       ! ------------------------------------------------------------------------------------------------------------------------------
       ! BEGIN: MAIN PARTICLE STEP
       ! ------------------------------------------------------------------------------------------------------------------------------



       DO k = 1,Nparticles


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
          forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i)) 


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

          ! ------------------------------------------------------------------------
          ! DONE: SQRT KALMAN
          ! ------------------------------------------------------------------------



       END DO ! k particles

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


       ! propagate sufficient statistics
       if (Nynonan > 0) then ! nothing to propagate if there was no observed data
          FORALL(k=1:Nparticles)         PREVsiglambdaT(k)    = PREVsiglambdaT(k) + lambdaDELTA(k) ** 2
          FORALL(k=1:Nparticles,i=1:Nsv) PREVhvarT(k,i)       = PREVhvarT(k,i)    + hDELTA(i,k) ** 2
          if (doInflationNoise) then
             i = 1
             FORALL(k=1:Nparticles)  PREVsigmaT(k,i)  = PREVsigmaT(k,i)   + resid(i,k) ** 2
          end if
          FORALL(k=1:Nparticles,i=2:Ny)  PREVsigmaT(k,i)      = PREVsigmaT(k,i)   + resid(i,k) ** 2 ! note: missing obs handled by zero values of resid

          PREVsiglambdaDof = PREVsiglambdaDof + 1
          PREVhvarDof      = PREVhvarDof + 1
          WHERE (.NOT. yNaN(:,j))
             PREVsigmaDof     = PREVsigmaDof + 1 
          END WHERE

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

             ! gap-AR1
             FORALL(k=1:Nparticles) shufflevec(k) = PREVa(ndx(k))
             PREVa = shufflevec
             FORALL(k=1:Nparticles) shufflevec(k) = PREVaV(ndx(k))
             PREVaV = shufflevec
          end if ! doSecondResamplingStep


       end if



    END DO ! j=1,T

  END SUBROUTINE MDDthetaCONSTlambdaTVP


  ! @\newpage\subsection{MDDthetaTVPlambdaCONST}@
  SUBROUTINE MDDthetaTVPlambdaCONST(doInflationNoise,T, logMDD, Ny, y, yNaN, Nparticles, Nx, NsigmaX, Nw, Ex0, sqrtVx00, p, a0, a0V,  sigaT, sigaDof, lambdaAlpha0, lambdaBeta0, Nsv, Eh0, Vh0, hvarT, hvarDof, sigmaT, sigmaDof, VSLstream)

    IMPLICIT NONE

    INTENT(INOUT) :: VSLstream 
    INTENT(IN)    :: T,Ny,y,yNaN, Nx, NsigmaX, Nw, Ex0, sqrtVx00, Nsv,Eh0,Vh0, Nparticles, p, lambdaAlpha0, lambdaBeta0, a0, a0V, sigaT, sigaDof, hvarT, hvarDof, sigmaT, sigmaDof
    INTENT(OUT) :: logMDD

    INTEGER :: J, I, K, T, Nparticles, Nx, Ny, Nsv, NsigmaX, p, Nw, Nsurveys




    logical, intent(in) :: doInflationNoise

    ! double precision, parameter :: minParticleWeight = 1.0d-12

    DOUBLE PRECISION, DIMENSION(Ny,T) :: y
    LOGICAL, DIMENSION(Ny,T) :: yNaN
    DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: PARTICLEweights
    DOUBLE PRECISION, DIMENSION(T)   :: logMDD


    ! APF llf correction
    DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
    DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax
    ! ! scale parameters
    DOUBLE PRECISION, DIMENSION(Nparticles)     :: DRAWsiga
    ! DOUBLE PRECISION, DIMENSION(Nparticles)     :: DRAWlambda
    ! DOUBLE PRECISION, DIMENSION(Nparticles,Nsv) :: DRAWhInno
    DOUBLE PRECISION, DIMENSION(Nsv,Nparticles)   :: hInno ! helper variable, note the transpose
    DOUBLE PRECISION, DIMENSION(Nparticles,Ny)  :: DRAWsigma

    ! particles
    DOUBLE PRECISION :: xposterior(Nx,Nparticles), h(Nsv, Nparticles), SVol(Nsv,Nparticles), llf(Nparticles)
    DOUBLE PRECISION :: kernelweights(Nparticles) ! for unnormalized weights
    DOUBLE PRECISION  :: kernelsum, loglikemax
    INTEGER :: ndx(Nparticles)
    DOUBLE PRECISION :: shufflevec(Nparticles)

    ! state space objects
    DOUBLE PRECISION :: xprior(Nx), logdetSigmaY

    DOUBLE PRECISION :: lambdaAlpha0, lambdaBeta0
    DOUBLE PRECISION, PARAMETER :: lambdaN0 = 1.0d0
    DOUBLE PRECISION, DIMENSION(Nparticles) :: lambdaAlpha, lambdaBeta
    DOUBLE PRECISION, DIMENSION(Nparticles) :: lambda
    DOUBLE PRECISION :: adrift(Nparticles), aPREV(Nparticles),  a0, a0V

    ! scale parameters
    DOUBLE PRECISION :: sigaT, hvarT(Nsv), sigmaT(Ny)
    INTEGER :: sigaDof, hvarDof(Nsv), sigmaDof(Ny)

    DOUBLE PRECISION :: PREVhvarT(Nparticles,Nsv), PREVsigmaT(Nparticles,Ny)
    INTEGER :: PREVhvarDof(Nsv), PREVsigmaDof(Ny)
    DOUBLE PRECISION :: hDELTA(Nsv,Nparticles), resid(Ny,Nparticles)
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

    logical, parameter :: doSecondResamplingStep = .false.
    integer, parameter :: dofAPFfattail = 0.0d0

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

    ! init sufficient statistics of scale parameters
    PREVsigaT   = sigaT
    PREVsigaDof = sigaDof

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
    lambdaAlpha = lambdaAlpha0
    lambdaBeta  = lambdaBeta0
    lambda      = drawBetas(lambdaAlpha, lambdaBeta, Nparticles, VSLstream) 



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
             write (*,*) 'DLYAP error (ygap0sqrtvariance -- init particlefilter MDDthetaTVPlambdaTVP)', errcode
             call savemat(gaptransition, 'gaptransition.debug')
             call savemat(gapshock0loadings, 'gapshock0loadings.debug')
             call savemat(ygap0sqrtvariance, 'ygap0sqrtvariance.debug')
             ! stop 1
          else
             sqrtVx0(ndxgap,ndxgap) = ygap0sqrtvariance
          end if
       end if
       vecSqrtSigmaX(:,k)     = vechU(sqrtVx0,Nx)

    END DO


    PARTICLEweights  = 1.0d0 / dble(Nparticles)
    logMDD         = 0.0d0

    ! uniform draws for systematic resampling
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


    ! workspace query for qr decomposition
    qrR     = 0.0d0
    qrlwork = qrquery(qrR)

    DO j=1,T

       Nynonan = count(.not. yNaN(:,j))
       ! print *, 'plefilter for t=', j, ' with ', Nynonan, 'obs' 
       ! ------------------------------------------------------------------------------------------------------------------------------
       ! BEGIN: APF RESAMPLE STEP
       ! ------------------------------------------------------------------------------------------------------------------------------

       ! Sigma
       i=1
       if (doInflationNoise) then
          call igammaDraws(DRAWsigma(:,i), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream)
       else
          DRAWsigma(:,i) = 0.0d0
          PREVsigmaT(:,i)  = 0.0d0
       end if
       DO i=2,Ny 
          call igammaDraws(DRAWsigma(:,i), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream) 
       END DO

       if (j > 1) then
          ! lambda parameter
          lambda = drawBetas(lambdaAlpha, lambdaBeta, Nparticles, VSLstream) 
       end if

       DO k = 1,Nparticles

          ! 2) Fill Particles into state space

          ! update A 

          A(ndxGapRE,ndxGapRE)           = adrift(k)
          A(ndxtrendSI,ndxtrendRE)       = 1 - lambda(k)
          A(ndxtrendSI,ndxtrendSI)       = lambda(k)
          A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - lambda(k)) * adrift(k)
          A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = lambda(k) * adrift(k)

          ! update sqrtR
          sqrtR = 0.0d0
          forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i)) 


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

          FORALL(k=1:Nparticles) shufflevec(k) = lambdaAlpha(ndx(k))
          lambdaAlpha = shufflevec

          FORALL(k=1:Nparticles) shufflevec(k) = lambdaBeta(ndx(k))
          lambdaBeta = shufflevec

          FORALL(k=1:Nparticles) shufflevec(k) = adrift(ndx(k))
          adrift = shufflevec

          ! reshuffle sufficient statistics for scale parameters
          FORALL(k=1:Nparticles) shufflevec(k) = PREVsigaT(ndx(k))
          PREVsigaT = shufflevec

          DO i=1,Nsv
             FORALL(k=1:Nparticles) shufflevec(k) = PREVhvarT(ndx(k),i)
             PREVhvarT(:,i) = shufflevec
          END DO

          DO i=1,Ny
             FORALL(k=1:Nparticles) shufflevec(k) = PREVsigmaT(ndx(k),i)
             PREVsigmaT(:,i) = shufflevec
          END DO

          DO i=1,Ny
             FORALL(k=1:Nparticles) shufflevec(k) = DRAWsigma(ndx(k),i)
             DRAWsigma(:,i) = shufflevec
          END DO

       end if

       ! ------------------------------------------------------------------------------------------------------------------------------
       ! END: APF STEP
       ! ------------------------------------------------------------------------------------------------------------------------------

       ! 0) draw scale parameters
       ! siga
       call igammaDraws(DRAWsiga, Nparticles, PREVsigaT, PREVsigaDof, VSLstream)
       DRAWsiga = sqrt(DRAWsiga)

       ! hInno
       DO i=1,Nsv 
          call igammaDraws(hinno(i,:), Nparticles, PREVhvarT(:,i), PREVhvarDof(i), VSLstream) 
       END DO
       hinno = sqrt(hinno)
       ! forall (i=1:Nsv,k=1:Nparticles) hInno(i,k) = DRAWhInno(k,i) ! helper variable, provides better aligned access to the j data


       ! 1) Draw Particles 
       errcode     = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, hDELTA, 0.0d0, 1.0d0)
       hDELTA      = hDELTA * hInno
       h           = h + hDELTA
       SVol        = exp(h * 0.5d0)


       aPREV      = adrift
       aDELTA     = drawdeltatruncnorm(Nparticles, aPREV, DRAWsiga, -1.0d0, 1.0d0, VSLstream)
       adrift     = aPREV + aDELTA

       ! ------------------------------------------------------------------------------------------------------------------------------
       ! BEGIN: MAIN PARTICLE STEP
       ! ------------------------------------------------------------------------------------------------------------------------------



       DO k = 1,Nparticles


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
          forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i)) 


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

          ! ------------------------------------------------------------------------
          ! DONE: SQRT KALMAN
          ! ------------------------------------------------------------------------



       END DO ! k particles

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


       ! propagate sufficient statistics
       if (Nynonan > 0) then ! nothing to propagate if there was no observed data
          FORALL(k=1:Nparticles)         PREVsigaT(k)         = PREVsigaT(k)      + aDELTA(k) ** 2
          FORALL(k=1:Nparticles,i=1:Nsv) PREVhvarT(k,i)       = PREVhvarT(k,i)    + hDELTA(i,k) ** 2
          if (doInflationNoise) then
             i = 1
             FORALL(k=1:Nparticles)  PREVsigmaT(k,i)  = PREVsigmaT(k,i)   + resid(i,k) ** 2
          end if
          FORALL(k=1:Nparticles,i=2:Ny)  PREVsigmaT(k,i)      = PREVsigmaT(k,i)   + resid(i,k) ** 2 ! note: missing obs handled by zero values of resid

          PREVsigaDof      = PREVsigaDof + 1
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



    END DO ! j=1,T

  END SUBROUTINE MDDthetaTVPlambdaCONST

  ! @\newpage\subsection{MDDthetaCONSTlambdaCONST}@
  SUBROUTINE MDDthetaCONSTlambdaCONST(doInflationNoise, T, logMDD, Ny, y, yNaN, Nparticles, Nxx, NsigmaXX, Nx, Nw, Ex0, sqrtVx00, p, a0, aV0, lambdaAlpha0, lambdaBeta0, Nsv, Eh0, Vh0, hvarT, hvarDof, sigmaT, sigmaDof, VSLstream)


    IMPLICIT NONE

    INTENT(INOUT) :: VSLstream ! , timer
    INTENT(IN)    :: T,Ny,y,yNaN, Nx, Nw, Ex0, sqrtVx00, Nsv,Eh0,Vh0, Nparticles, p, lambdaAlpha0, lambdaBeta0,  a0, aV0, hvarT, hvarDof, sigmaT, sigmaDof
    INTENT(IN) :: Nxx, NsigmaXX
    INTENT(OUT) :: logMDD
    INTEGER :: Nxx, NsigmaXX

    INTEGER :: J, I, K, T, Nparticles, Nx, Ny, Nsv, p, Nw, Nsurveys




    logical, intent(in) :: doInflationNoise

    ! double precision, parameter :: minParticleWeight = 1.0d-12

    DOUBLE PRECISION, DIMENSION(Ny,T) :: y
    LOGICAL, DIMENSION(Ny,T) :: yNaN
    DOUBLE PRECISION, DIMENSION(Nparticles,0:T) :: PARTICLEweights
    DOUBLE PRECISION, DIMENSION(T)   :: logMDD


    ! APF llf correction
    DOUBLE PRECISION, DIMENSION(Nparticles)   :: APFkernelweights, APFlike
    DOUBLE PRECISION  :: APFkernelsum, loglikeAPFmax
    ! scale parameters
    ! DOUBLE PRECISION, DIMENSION(Nparticles,Nsv) :: DRAWhInno
    DOUBLE PRECISION, DIMENSION(Nsv,Nparticles)   :: hInno ! helper variable, note the transpose
    DOUBLE PRECISION, DIMENSION(Nparticles,Ny)  :: DRAWsigma

    ! particles
    DOUBLE PRECISION :: xxposterior(Nxx,Nparticles),  h(Nsv, Nparticles), SVol(Nsv,Nparticles), llf(Nparticles)
    DOUBLE PRECISION :: kernelweights(Nparticles) ! for unnormalized weights
    DOUBLE PRECISION  :: kernelsum, loglikemax
    INTEGER :: ndx(Nparticles)
    DOUBLE PRECISION :: shufflevec(Nparticles)

    ! state space objects
    DOUBLE PRECISION :: xxprior(Nxx), logdetSigmaY

    DOUBLE PRECISION :: lambdaAlpha0, lambdaBeta0
    DOUBLE PRECISION, PARAMETER :: lambdaN0 = 1.0d0
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

    logical, parameter :: doSecondResamplingStep = .false.
    integer, parameter :: dofAPFfattail = 0.0d0

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
             write (*,*) 'DLYAP error (ygap0sqrtvariance -- init particlefilter MDDthetaCONSTlambdaCONST)', errcode
             call savemat(gaptransition, 'gaptransition.debug')
             call savemat(gapshock0loadings, 'gapshock0loadings.debug')
             call savemat(ygap0sqrtvariance, 'ygap0sqrtvariance.debug')
             ! stop 1
          else
             sqrtVxx0(ndxgap,ndxgap) = ygap0sqrtvariance
          end if
       end if
       vecSqrtSigmaXX(:,k)     = vechU(sqrtVxx0,Nxx)

    END DO


    PARTICLEweights  = 1.0d0 / dble(Nparticles)
    logMDD           = 0.0d0

    ! uniform draws for systematic resampling
    errcode = vdrnguniform(VSLmethodUniform, VSLstream, 2*T, uniformdraws, 0.0d0, 1.0d0)


    ! workspace query for qr decomposition
    qrR     = 0.0d0
    qrlwork = qrquery(qrR)

    DO j=1,T

       Nynonan = count(.not. yNaN(:,j))
       ! print *, 'plefilter for t=', j, ' with ', Nynonan, 'obs' 
       ! ------------------------------------------------------------------------------------------------------------------------------
       ! BEGIN: APF RESAMPLE STEP
       ! ------------------------------------------------------------------------------------------------------------------------------

       ! Sigma
       i=1
       if (doInflationNoise) then
          call igammaDraws(DRAWsigma(:,i), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream)
       else
          DRAWsigma(:,i) = 0.0d0
          PREVsigmaT(:,i)  = 0.0d0
       end if
       DO i=2,Ny 
          call igammaDraws(DRAWsigma(:,i), Nparticles, PREVsigmaT(:,i), PREVsigmaDof(i), VSLstream) 
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


       DO k = 1,Nparticles

          ! 2) Fill Particles into state space

          ! update A 

          A(ndxGapRE,ndxGapRE)           = adrift(k)
          A(ndxtrendSI,ndxtrendRE)       = 1 - lambda(k)
          A(ndxtrendSI,ndxtrendSI)       = lambda(k)
          A(ndxgapSIstart:ndxgapSIstop,ndxgapREstart:ndxgapREstop) = (1 - lambda(k)) * adrift(k)
          A(ndxgapSIstart:ndxgapSIstop,ndxgapSIstart:ndxgapSIstop) = lambda(k) * adrift(k)

          ! update sqrtR
          sqrtR = 0.0d0
          forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i)) 


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
             FORALL(k=1:Nparticles) shufflevec(k) = DRAWsigma(ndx(k),i)
             DRAWsigma(:,i) = shufflevec
          END DO

       end if

       ! ------------------------------------------------------------------------------------------------------------------------------
       ! END: APF STEP
       ! ------------------------------------------------------------------------------------------------------------------------------

       ! 0) draw scale parameters
       ! hInno
       DO i=1,Nsv 
          call igammaDraws(hinno(i,:), Nparticles, PREVhvarT(:,i), PREVhvarDof(i), VSLstream) 
       END DO
       hinno = sqrt(hinno)
       ! forall (i=1:Nsv,k=1:Nparticles) hInno(i,k) = DRAWhInno(k,i) ! helper variable, provides better aligned access to the j data


       ! 1) Draw Particles 
       errcode     = vdrnggaussian(VSLmethodGaussian, VSLstream, Nsv * Nparticles, hDELTA, 0.0d0, 1.0d0)
       hDELTA      = hDELTA * hInno
       h           = h + hDELTA
       SVol        = exp(h * 0.5d0)



       ! ------------------------------------------------------------------------------------------------------------------------------
       ! BEGIN: MAIN PARTICLE STEP
       ! ------------------------------------------------------------------------------------------------------------------------------



       DO k = 1,Nparticles


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
          forall (i=1:Ny) sqrtR(i,i) = sqrt(DRAWsigma(k,i)) 


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

          ! ------------------------------------------------------------------------
          ! DONE: SQRT KALMAN
          ! ------------------------------------------------------------------------



       END DO ! k particles

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



    END DO ! j=1,T

  END SUBROUTINE MDDthetaCONSTlambdaCONST
END MODULE cambridgebox


