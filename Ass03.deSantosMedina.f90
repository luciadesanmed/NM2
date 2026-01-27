!Assignment 3:
!Program to integrate the linear advection eq using the CTCS scheme in the domain 0<x<500 with advection velocity u=-0.31 m/s
!Let dx=0.1 m and assume periodic BC

!Initial shape:
!0.1 for x<200, 2.0 for 200<=x<250, 1.0 for 250<=x<=300, 0.1 for x>300

!Integrate forward and show solutions from t=0 to t=1000 every 200s and explain the characteristics of the solution

!Repeat the exercise implementing the RA filter with alfa=0.1
!Repeat the exercise implementing the RAW filter with alfa=0.05 and beta=0.53

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
   IF (x(n) < 200) THEN 
     res(n) = 0.1
   ELSE IF (200 <= x(n) .AND. x(n) < 250) THEN 
    res(n) = 2.0
   ELSE IF (250 <= x(n) .AND. x(n) <= 300) THEN
    res(n) = 1.0
   ELSE IF (x(n) > 300) THEN 
    res(n) = 0.1
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

END MODULE linear_advection

!Begin program:
PROGRAM ass3
 USE linear_advection
 IMPLICIT NONE
 REAL :: x0, x1, dx, u, dt, c, t, t0, t1, tp, alpha, beta
 REAL, DIMENSION(:), ALLOCATABLE :: x, phi_old, phi_now, phi00, phi_new, d
 REAL, DIMENSION(:,:), ALLOCATABLE :: phi_plot
 INTEGER :: nx, n, ios, uni, nn, ii, jj
 uni=10

 nn=0
 u = -0.31
 x0 = 0.0
 x1 = 500.0
 dx = 0.1
 t0 = 0.0
 t1 = 1000.0
 tp = 200.0
 alpha = 0.1 !RA FILTER -> 0.1, RAW -> 0.05
 beta = 0.53
 c = -0.1
 !c=1
 dt = ABS(c * dx / u) 
 !c = u * (dt/dx)
 nx = ANINT((x1-x0)/dx) + 1
 
 ALLOCATE(x(nx), phi00(nx), phi_old(nx), phi_now(nx), phi_new(nx), phi_plot(5,nx), d(nx))

 x = xvec(x0,x1,dx)
 phi00 = phi0(x)
 phi_old = phi0(x)
 t =  t0

  IF (u > 0) THEN
   phi_now = ftbs(phi_old,c)
  ELSE
   phi_now = ftfs(phi_old,c)
  END IF

 DO 
  IF (t>t1) EXIT
  phi_new = ctcs(phi_old, phi_now, c)
  
  !No filter:
  phi_old = phi_now
  phi_now = phi_new

  !RA filter:
  !d = alpha * (phi_old + phi_new - 2.0 * phi_now)
  !phi_old = phi_now + d
  !phi_now = phi_new

  !RAW filter:
  !d = alpha * (phi_old + phi_new - 2.0 * phi_now)
  !phi_old = phi_now + beta*d
  !phi_now = phi_new + (1.0 - beta)*d

  t = t+dt
  
  IF (MOD(t,tp) < dt) THEN
   nn=nn+1

   IF (nn<=5) THEN
    phi_plot(nn,:) = phi_now
   END IF

  END IF

 END DO

 !Write data for plotting:
 OPEN(UNIT=uni, IOSTAT=ios, FILE='data3_nofilter.dat', STATUS='new', ACTION='write')
 DO n=1, nx
  WRITE(uni, *) x(n), phi00(n), phi_plot(1,n), phi_plot(2,n), phi_plot(3,n), phi_plot(4,n), phi_plot(5,n)
 END DO 

 CLOSE(uni)
END PROGRAM ass3

