!By: Ana Lucia de Santos Medina

!Exam code:

!Program to integrate the linear advection eq using the algorithm defined in 1a in the periodic domain 0<x<750, dx= 0.15 m, and a
!constant advection velocity u= 0.17 m/s

!Initial shape:
!0.0 for x<65.0, 5.0 for 65<=x<85, 10.0 for 85<=x<=105, 0.0 for x>105

!Integrate forward and show solutions from t=0 to t=1500 and plot every 300s. Choose a dt to respect the CFL

!e) Apply a time filter.
!f) Add a diffusion term.

!Module with functions:

MODULE linear_advection
 IMPLICIT NONE
 CONTAINS
 
 !X vector constructor:
 FUNCTION xvec(x0,x1,dx) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x0, x1, dx
  INTEGER :: n, nx
  REAL :: x
  REAL, DIMENSION(:), ALLOCATABLE :: res
  
  nx = ANINT((x1-x0)/dx) + 1

  ALLOCATE(res(nx))

  x=x0
  DO n=1, nx
   res(n)=x
   x=x+dx
  END DO
  
 END FUNCTION xvec

 !Initial condition:
 FUNCTION phi0(x) RESULT(res)
  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: x
  INTEGER :: n
  REAL, DIMENSION(:), ALLOCATABLE :: res
   
  ALLOCATE(res(SIZE(x)))

  DO n= 1, SIZE(x)
   IF (x(n) < 65.0) THEN 
     res(n) = 0.0
   ELSE IF (65.0 <= x(n) .AND. x(n) < 85.0) THEN 
    res(n) = 5.0
   ELSE IF (85.0 <= x(n) .AND. x(n) <= 105.0) THEN
    res(n) = 10.0
   ELSE IF (x(n) > 105.0) THEN 
    res(n) = 0.0
   END IF 
  END DO

 END FUNCTION phi0

 !FTFS:
 FUNCTION ftfs(phi, c) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: c
  REAL, DIMENSION(:), INTENT(IN) :: phi
  REAL, DIMENSION(:), ALLOCATABLE :: res
  INTEGER :: n, sz

  sz=SIZE(phi)

  ALLOCATE(res(sz))

  DO n=1, sz-1
    res(n) = (1+c) * phi(n) - c * phi(n+1)
  END DO

  res(sz) = (1+c) * phi(sz) - c * phi(1)

 END FUNCTION ftfs

 !FTBS scheme:
 FUNCTION ftbs(phi, c) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: c
  REAL, DIMENSION(:), INTENT(IN) :: phi
  REAL, DIMENSION(:), ALLOCATABLE :: res
  INTEGER :: n, sz

  sz=SIZE(phi)

  ALLOCATE(res(sz))

  DO n=2, sz
    res(n) = (1-c) * phi(n) + c * phi(n-1)
  END DO
  
  res(1)= (1-c) * phi(1) + c * phi(sz)

 END FUNCTION ftbs

 !CTCS scheme:
 FUNCTION ctcs(phi_old, phi_now, c) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: c
  REAL, DIMENSION(:), INTENT(IN) :: phi_old, phi_now
  REAL, DIMENSION(:), ALLOCATABLE :: res
  INTEGER :: n, sz

  sz=SIZE(phi_old)

  ALLOCATE(res(sz))

  DO n=2, sz-1
    res(n) = phi_old(n) - c * (phi_now(n+1) - phi_now(n-1))
  END DO

  res(1) = phi_old(1) - c * (phi_now(2) - phi_now(sz))
  res(sz) = phi_old(sz) - c * (phi_now(1) - phi_now(sz-1))

 END FUNCTION ctcs

!FTCS:
 FUNCTION ftcs_diff(phi, c, dc) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: c, dc
  REAL, DIMENSION(:), INTENT(IN) :: phi
  REAL, DIMENSION(:), ALLOCATABLE :: res
  INTEGER :: n, sz

  sz=SIZE(phi)

  ALLOCATE(res(sz))

  DO n=2, sz-1
    res(n) = phi(n) - 0.5 * c * (phi(n+1) - phi(n-1)) + &
            dc * (phi(n+1) - 2.0*phi(n) + phi(n-1)) 
  END DO
  
  res(1) = phi(1) - 0.5 * c * (phi(2) - phi(sz)) + &
          dc * (phi(2) - 2.0 * phi(1) + phi(sz))
 
  res(sz) = phi(sz) - 0.5 * c * (phi(1) - phi(sz-1)) + &
          dc * (phi(1) - 2.0 * phi(sz) + phi(sz-1))

 END FUNCTION ftcs_diff

!CTCS scheme:
 FUNCTION ctcs_diff(phi_old, phi_now, c, dc) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: c, dc
  REAL, DIMENSION(:), INTENT(IN) :: phi_old, phi_now
  REAL, DIMENSION(:), ALLOCATABLE :: res
  INTEGER :: n, sz

  sz=SIZE(phi_old)

  ALLOCATE(res(sz))

  DO n=2, sz-1
    res(n) = phi_old(n) - c * (phi_now(n+1) - phi_now(n-1)) + &
            2.0 * dc * (phi_old(n+1) - 2.0 * phi_old(n) + phi_old(n-1))
  END DO

  res(1) = phi_old(1) - c * (phi_now(2) - phi_now(sz)) + &
          2.0 * dc * (phi_old(2) - 2.0 * phi_old(1) + phi_old(sz))
  res(sz) = phi_old(sz) - c * (phi_now(1) - phi_now(sz-1)) + &
          2.0 * dc * (phi_old(1) - 2.0 * phi_old(sz) + phi_old(sz-1))

 END FUNCTION ctcs_diff

END MODULE linear_advection

!Begin program:
PROGRAM exam
 USE linear_advection
 IMPLICIT NONE
 REAL :: x0, x1, dx, u, dt, c, t, t0, t1, tp, alpha, beta, k, dc
 REAL, DIMENSION(:), ALLOCATABLE :: x, phi_old, phi_now, phi00, phi_new, d
 REAL, DIMENSION(:,:), ALLOCATABLE :: phi_plot
 INTEGER :: diff, nx, n, ios, uni, nn, filter, np, ii
 uni=10

 diff = 0 !0 no diffusion, 1 with diffusion, No filter applied in diffusion case
 filter = 2  !0 no filter, 1 RA filter, 2 RAW filter 

 nn = 0
 u = 0.17
 x0 = 0.0  
 x1 = 750.0
 dx = 0.15 
 t0 = 0.0
 t1 = 1500.0
 tp = 300.0
 np = ANINT(t1/tp)
 k = 0.01 !diffusion coeff

 IF (filter == 1) THEN
  alpha = 0.1
 ELSEIF (filter == 2) THEN
  alpha = 0.05
 END IF

 beta = 0.53
 
 IF (diff == 0) THEN
  !No diffusion:
  c = 0.1
  dt = ABS(c * dx / u) 
 ELSE
  !Diffusion:
  dt = 0.01
  c = u * dt/dx
  dc = k * dt/(dx**2)
 END IF

 nx = ANINT((x1-x0)/dx) + 1
 
 ALLOCATE(x(nx), phi00(nx), phi_old(nx), phi_now(nx), phi_new(nx), phi_plot(np,nx), d(nx))

 x = xvec(x0,x1,dx)
 phi00 = phi0(x)
 phi_old = phi0(x)
 t =  t0

 IF (diff == 0) THEN
 !No diffusion:
  IF (u > 0) THEN
   phi_now = ftbs(phi_old,c)
  ELSE
   phi_now = ftfs(phi_old,c)
  END IF
 ELSE
!For diffusion:
  phi_now = ftcs_diff(phi_old, c, dc)
 END IF

 DO 
  IF (t>t1) EXIT

  IF (diff == 0) THEN
   !No diffusion:
   phi_new = ctcs(phi_old, phi_now, c)
  ELSE
  !With diffusion:
   phi_new = ctcs_diff(phi_old, phi_now, c, dc)
  END IF

  IF (diff == 0) THEN
   !Filters:
   IF (filter == 0) THEN
    phi_old = phi_now
    phi_now = phi_new
   !RA filter:
   ELSEIF (filter == 1) THEN
    d = alpha * (phi_old + phi_new - 2.0 * phi_now)
    phi_old = phi_now + d
    phi_now = phi_new
   !RAW filter:
   ELSEIF (filter == 2) THEN
    d = alpha * (phi_old + phi_new - 2.0 * phi_now)
    phi_old = phi_now + beta*d
    phi_now = phi_new + (1.0 - beta)*d
   END IF
  ELSE
   phi_old = phi_now
   phi_now = phi_new
  END IF

  t = t+dt
  
  IF (MOD(t,tp) < dt) THEN
   nn=nn+1

   IF (nn<=np) THEN
    phi_plot(nn,:) = phi_now
   END IF

  END IF

 END DO

 !Write data for plotting:

 !For diffusion: 
 IF (diff == 1) THEN
  OPEN(UNIT=uni, IOSTAT=ios, FILE='exam_diffusion.dat', STATUS='new', ACTION='write')
 ELSE
  IF (filter == 0 ) THEN
   OPEN(UNIT=uni, IOSTAT=ios, FILE='exam_nofilter.dat', STATUS='new', ACTION='write')
  ELSEIF (filter == 1) THEN
   OPEN(UNIT=uni, IOSTAT=ios, FILE='exam_RAfilter.dat', STATUS='new', ACTION='write')       
  ELSEIF (filter == 2) THEN
   OPEN(UNIT=uni, IOSTAT=ios, FILE='exam_RAWfilter.dat', STATUS='new', ACTION='write')       
  END IF
 END IF

 DO n=1, nx
  WRITE(uni, *) x(n), phi00(n), phi_plot(:,n)
 END DO 

 CLOSE(uni)

END PROGRAM exam

