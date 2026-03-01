!Assignment 4:
!Program to integrate the linear advection eq using the SL scheme with linear interpolation in the domain 0<x<1000 with advection
!velocity u=0.75. Let dx=0.5 and assume periodic BC. Assume the initial shape to be:
!0.0 for x<400, 0.1(x-400.0) for 400<x<500, 20.0 -0.1(x-400.0) for 500<x<600, 0.0 for x>600

!Integrate forward and show solutions from t=0 ro t=2000 every 250s and explain the characteristics of the solution

!Repeat the exercise using cubic interpolation

!Do we need to worry about the CFL condition? Explain why or why not


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
   IF (x(n) < 400) THEN 
     res(n) = 0.0
   ELSE IF (400 <= x(n) .AND. x(n) < 500) THEN 
    res(n) = 0.1*(x(n)-400.0)
   ELSE IF (500 <= x(n) .AND. x(n) <= 600) THEN
    res(n) = 20.0-0.1*(x(n)-400.0)
   ELSE IF (x(n) > 600) THEN 
    res(n) = 0.0
   END IF 
  END DO

 END FUNCTION phi0

END MODULE linear_advection

!Begin program:
PROGRAM ass4
 USE linear_advection
 IMPLICIT NONE
 REAL :: x0, x1, dx, u, dt, t, t0, t1, tp, xdep, alpha
 REAL, DIMENSION(:), ALLOCATABLE :: x, phi_now, phi00, phi_new
 REAL, DIMENSION(:,:), ALLOCATABLE :: phi_plot
 INTEGER :: nt, nx, n, ios, uni, m, mp, mp2, mf, nn
 uni=10
 
 nn=0
 u = 0.75
 x0 = 0.0
 x1 = 1000.0
 dx = 0.5
 t0 = 0.0
 t1 = 2000.0
 tp = 250.0
 dt = 1.0
 nx = ANINT((x1-x0)/dx) + 1
 nt = ANINT(t1/tp)
 
 ALLOCATE(x(nx), phi00(nx), phi_now(nx), phi_new(nx), phi_plot(nt,nx))

 x = xvec(x0,x1,dx)
 phi00 = phi0(x)
 phi_now = phi0(x)
 t =  t0

 
 !With linear interpolation:
 DO 
  IF (t>t1) EXIT
  
  DO n=1, nx
   xdep = x0 + MOD(x(n) - u * dt, x1 - x0)
   IF (xdep < x0) THEN 
     xdep = x1 + xdep
   END IF
   m = FLOOR(xdep/dx)
   alpha = (xdep/dx)-m
   m = m+1
   mp = m+1
   IF (mp>nx) THEN 
    mp = mp-nx
   END IF
   phi_new(n)=(1-alpha)*phi_now(m)+alpha*phi_now(mp)
  END DO
  phi_now = phi_new
  t = t+dt
  
  IF (MOD(t,tp) < dt) THEN
   nn=nn+1

   IF (nn<=nt) THEN
    phi_plot(nn,:) = phi_now
   END IF

  END IF

 END DO

 !Write data for plotting:
 OPEN(UNIT=uni, IOSTAT=ios, FILE='data4_linear.dat', STATUS='new', ACTION='write')
 DO n=1, nx
  WRITE(uni, *) x(n), phi00(n), phi_plot(:,n)
 END DO 

 CLOSE(uni)

 t = t0
 phi_now = phi0(x)
 nn = 0

 !With cubic interpolation:
 DO 
  IF (t>t1) EXIT
  
  DO n=1, nx
   xdep = x0 + MOD(x(n) - u * dt, x1 - x0)
   IF (xdep < x0) THEN 
     xdep = x1 + xdep
   END IF
   m = FLOOR(xdep/dx)
   alpha = (xdep/dx)-m
   m = m+1
   mp = m+1
   mp2 = m+2
   mf = m-1

   IF (mp>nx) THEN 
    mp = mp - nx
   END IF
   IF (mp2>nx) THEN
    mp2 = mp2 - nx
   END IF
   IF (mf<1) THEN
    mf = nx
   END IF

   phi_new(n)=-(alpha*(1.0-alpha**2)/6.0)*phi_now(mp2) &
           + ((alpha*(1.0+alpha)*(2.0-alpha))/2.0)*phi_now(mp) &
           + (((1-alpha**2)*(2-alpha))/2.0)*phi_now(m) &
           - ((alpha*(1-alpha)*(2-alpha))/6.0)*phi_now(mf)
  END DO
  phi_now = phi_new
  t = t+dt
  
  IF (MOD(t,tp) < dt) THEN
   nn=nn+1

   IF (nn<=nt) THEN
    phi_plot(nn,:) = phi_now
   END IF

  END IF

 END DO

 !Write data for plotting:
 OPEN(UNIT=uni, IOSTAT=ios, FILE='data4_cubic.dat', STATUS='new', ACTION='write')
 DO n=1, nx
  WRITE(uni, *) x(n), phi00(n), phi_plot(:,n)
 END DO 

 CLOSE(uni)

END PROGRAM ass4

