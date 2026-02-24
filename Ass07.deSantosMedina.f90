!Assignment 7:
!Program to integrate the advection-diffusion eq using the scheme in [12] in the domain 0<x<1000 with advection velocity u=0.95 m/s
!and diff coef k=0.029. Let dx=0.2 and assume periodic BC

!Initial shape:
!0.0 for x<400, 0.01*(x-400.0) for 400<=x<550, 2.0-0.01*(x-400.0) for 500<=x<=600, 0.0 for x>600

!Integrate forward and show solutions from t=0 to t=2000 every 500s. What jhappens if you increase the spatial resolution dx? Set
!dx=0.05

!Apply a RAW filter with alpa=0.1 and beta=0.53

!Module with functions:

MODULE advection_diffusion
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
   IF (x(n) < 400.0) THEN 
     res(n) = 0.0
   ELSE IF (400.0 <= x(n) .AND. x(n) < 500.0) THEN 
    res(n) = 0.01*(x(n)-400.0)
   ELSE IF (500.0 <= x(n) .AND. x(n) <= 600.0) THEN
    res(n) = 2.0 - 0.01*(x(n)-400.0)
   ELSE IF (x(n) > 600.0) THEN 
    res(n) = 0.0
   END IF 
  END DO

 END FUNCTION phi0

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

END MODULE advection_diffusion

!Begin program:
PROGRAM ass7
 USE advection_diffusion
 IMPLICIT NONE
 REAL :: x0, x1, dx, u, dt, c, t, t0, t1, tp, alpha, beta, k, dc
 REAL, DIMENSION(:), ALLOCATABLE :: x, phi_old, phi_now, phi00, phi_new, d
 REAL, DIMENSION(:,:), ALLOCATABLE :: phi_plot
 INTEGER :: nx, n, ios, uni, nn, ii, jj
 uni=10

 nn= 0
 u = 0.95
 x0 = 0.0
 x1 = 1000.0
 dx = 0.05
 t0 = 0.0
 t1 = 2000.0
 tp = 500.0
 alpha = 0.1
 beta = 0.53
 k = 0.029
 dt = 0.01
 c = u * dt/dx
 dc = k * dt/(dx**2)
 dt = ABS(c * dx / u) 
 nx = ANINT((x1-x0)/dx) + 1
 
 ALLOCATE(x(nx), phi00(nx), phi_old(nx), phi_now(nx), phi_new(nx), phi_plot(4,nx), d(nx))

 x = xvec(x0,x1,dx)
 phi00 = phi0(x)
 phi_old = phi0(x)
 t =  t0

 phi_now = ftcs_diff(phi_old, c, dc)

 DO 
  IF (t>t1) EXIT
  phi_new = ctcs_diff(phi_old, phi_now, c, dc)
  
  !RAW filter:
  d = alpha * (phi_old + phi_new - 2.0 * phi_now)
  phi_old = phi_now + beta*d
  phi_now = phi_new + (1.0 - beta)*d

  t = t+dt
  
  IF (MOD(t,tp) < dt) THEN
   nn=nn+1

   IF (nn<=4) THEN
    phi_plot(nn,:) = phi_now
    PRINT*, nn
    !PRINT*, phi_now
   END IF

  END IF

 END DO

 !Write data for plotting:
 OPEN(UNIT=uni, IOSTAT=ios, FILE='data7_dx05dt01.dat', STATUS='new', ACTION='write')
 DO n=1, nx
  WRITE(uni, *) x(n), phi00(n), phi_plot(1,n), phi_plot(2,n), phi_plot(3,n), phi_plot(4,n)
 END DO 

 CLOSE(uni)
END PROGRAM ass7

