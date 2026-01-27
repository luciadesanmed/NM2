!Assignment 2:
!Program to integrate the linear advection eq using an upwind scheme in the domain 0<x<100 with advection velocity
!u=0.087 m/s. dx=0.10 m, assuming periodic BC.
!Initial shape:
!0.0 for x<40, sin(pi*(x-40)/30)^2 for 40<x<70, 0.0 for x>70

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
  REAL :: pi
  INTEGER :: n
  REAL, DIMENSION(:), ALLOCATABLE :: res
   
  ALLOCATE(res(SIZE(x)))

  pi = 4.0 * ATAN(1.0)
  
  DO n= 1, SIZE(x)
   IF (x(n) < 40) THEN 
     res(n) = 0.0
   ELSE IF (40 < x(n) .AND. x(n) < 70) THEN 
    res(n) = (SIN(pi * (x(n)-40.0)/30.0))**2
   ELSE IF (x(n)>70) THEN 
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

  res(1) = (1-c) * phi(1) + c * phi(sz)

 END FUNCTION ftbs

END MODULE linear_advection

!Begin program:
PROGRAM ass2
 USE linear_advection
 IMPLICIT NONE
 REAL :: x0, x1, dx, u, dt, c, t, t0, t1, tp
 REAL, DIMENSION(:), ALLOCATABLE :: x, phi, phi00, phi_new
 REAL, DIMENSION(:,:), ALLOCATABLE :: phi_plot
 INTEGER :: nx, n, ios, uni, nn, ii, jj
 uni=10

 nn=0
 u = 0.087
 x0 = 0.0
 x1 = 100.0
 dx = 0.1
 t0 = 0.0
 t1 = 1000.0
 tp = 200.0
 dt = dx / u !With c=1
 c = u * (dt/dx) 
 nx = ANINT((x1-x0)/dx) + 1
 
 ALLOCATE(x(nx), phi00(nx), phi(nx), phi_new(nx), phi_plot(5,nx))

 x=xvec(x0,x1,dx)
 phi00=phi0(x)
 phi=phi0(x)
 t =  t0

 DO
  IF (t>t1) EXIT
  
  IF (u > 0) THEN
   phi_new = ftbs(phi,c)
  ELSE
   phi_new = ftfs(phi,c)
  END IF

  phi = phi_new
  t = t+dt
  
  IF (MOD(t,tp) < dt) THEN
   nn=nn+1

   IF (nn<=5) THEN
    phi_plot(nn,:) = phi
   END IF

  END IF

 END DO

 !Write data for plotting:
 OPEN(UNIT=uni, IOSTAT=ios, FILE='data2.dat', STATUS='new', ACTION='write')
 DO n=1, nx
  WRITE(uni, *) x(n), phi00(n), phi_plot(1,n), phi_plot(2,n), phi_plot(3,n), phi_plot(4,n), phi_plot(5,n)
 END DO 

 CLOSE(uni)
END PROGRAM ass2

