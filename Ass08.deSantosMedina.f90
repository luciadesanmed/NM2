!Assignment 8:
!Program to integrate the gravity wave eqs in [1] using discretization in [5] and [7] in the domain 0<x<1000 with mean height such
!that gh=1 m^2/s^2
!Let dx=0.5 and assume periodic BC

!Initial shape:
!0.0 for x<400, (sin(pi*((x-400.0)/200.0)))^2 0.0 for x>600

!and initial velocity u(x,0)=0.0

!Integrate forward and show solutions from t=0 to t=2000s every 200s and explain the characteristics of the solution

!Attention: Select dt to have a stable scheme remembering the CFL condition

!Module with functions:

MODULE gravity_waves
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
  REAL :: pi
  REAL, DIMENSION(:), ALLOCATABLE :: res
   
  pi= 4.0*ATAN(1.0) 

  ALLOCATE(res(SIZE(x)))

  DO n= 1, SIZE(x)
   IF (x(n) < 400.0) THEN 
     res(n) = 0.0
   ELSE IF (400.0 <= x(n) .AND. x(n) <= 600.0) THEN 
    res(n) = (SIN(pi*((x(n)-400.0)/200.0)))**2
   ELSE IF (x(n) > 600.0) THEN 
    res(n) = 0.0
   END IF 
  END DO

 END FUNCTION phi0

 !U ft scheme 1, 2:
 FUNCTION ueq_ft(u, phi, dtdx) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: dtdx
  REAL, DIMENSION(:), INTENT(IN) :: phi, u
  REAL, DIMENSION(:), ALLOCATABLE :: res
  INTEGER :: n, sz

  sz=SIZE(phi)

  ALLOCATE(res(sz))

  DO n = 2, sz-1
   res(n) = u(n) - 0.5 * dtdx * (phi(n+1) - phi(n-1))
  END DO

  res(1) = u(1) - 0.5 * dtdx * (phi(2) - phi(sz))
 
  res(sz) = u(sz) - 0.5 * dtdx * (phi(1) - phi(sz-1))

 END FUNCTION ueq_ft

 !Phi ft scheme 1, 2:
 FUNCTION peq_ft(u, phi, c) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: c
  REAL, DIMENSION(:), INTENT(IN) :: phi, u
  REAL, DIMENSION(:), ALLOCATABLE :: res
  INTEGER :: n, sz

  sz = SIZE(phi)

  ALLOCATE(res(sz))

  DO n = 2, sz-1
   res(n) = phi(n) - 0.5 * c * (u(n+1) - u(n-1))
  END DO
  
  res(1) = phi(1) - 0.5 * c * (u(2) - u(sz))
  res(sz) = phi(sz) - 0.5 * c * (u(1) - u(sz-1))

  END FUNCTION peq_ft

 FUNCTION ueq(u, phi, dtdx) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: dtdx
  REAL, DIMENSION(:), INTENT(IN) :: phi, u
  REAL, DIMENSION(:), ALLOCATABLE :: res
  INTEGER :: n, sz

  sz = SIZE(phi)
  ALLOCATE(res(sz))

  DO n = 2, sz-1
   res(n) = u(n) - dtdx * (phi(n+1) - phi(n-1))
  END DO
 
  res(1) = u(1) - dtdx * (phi(2) - phi(sz))
  res(sz) = u(sz) - dtdx * (phi(1) - phi(sz-1))

 END FUNCTION ueq

 FUNCTION peq(u, phi, c) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: c
  REAL, DIMENSION(:), INTENT(IN) :: phi, u
  REAL, DIMENSION(:), ALLOCATABLE :: res
  INTEGER :: n, sz

  sz = SIZE(phi)
  ALLOCATE(res(sz))

  DO n = 2, sz-1
   res(n) = phi(n) - c * (u(n+1) - u(n-1))
  END DO

  res(1) = phi(1) - c * (u(2) - u(sz))
  res(sz) = phi(sz) - c * (u(1) - u(sz-1))

 END FUNCTION peq

END MODULE gravity_waves

!Begin program:
PROGRAM ass8
 USE gravity_waves
 IMPLICIT NONE
 REAL :: x0, x1, dx, dt, c, t, t0, t1, tp, dtdx, h
 REAL, DIMENSION(:), ALLOCATABLE :: x, p_old, p_now, p_new, u_old, u_now, u_new, phi00
 REAL, DIMENSION(:,:), ALLOCATABLE :: p_plot, u_plot
 INTEGER :: nt, nx, n, ios, uni, nn, ii, jj
 uni=10

 nn= 0
 h = 1.0
 x0 = 0.0
 x1 = 1000.0
 dx = 0.5
 t0 = 0.0
 t1 = 2000.0
 tp = 200.0
 c = 1
 dt = c * dx/ SQRT(h)
 dtdx = dt/dx
 nx = ANINT((x1-x0)/dx) + 1
 nt = ANINT(t1/tp)
 
 ALLOCATE(x(nx), phi00(nx), u_old(nx), u_now(nx), u_new(nx), p_old(nx), p_now(nx), p_new(nx), p_plot(nt,nx), u_plot(nt, nx))

 x = xvec(x0,x1,dx)
 phi00 = phi0(x)
 p_old = phi0(x)
 u_old = 0.0
 t =  t0

 p_now = peq_ft(u_old, p_old, c)
 u_now = ueq_ft(u_old, p_old, dtdx)

 DO 
  IF (t>t1) EXIT

  u_new = ueq(u_old, p_now, dtdx)
  p_new = peq(u_now, p_old, c)

  p_old = p_now
  p_now = p_new
  u_old = u_now
  u_now = u_new

  t = t+dt
  
  IF (MOD(t,tp) < dt) THEN
   nn=nn+1

   IF (nn<=nt) THEN
    p_plot(nn,:) = p_now
    u_plot(nn,:) = u_now
   END IF

  END IF

 END DO

 PRINT*, dt

 !Write data for plotting:
 OPEN(UNIT=uni, IOSTAT=ios, FILE='data8_scheme1.dat', STATUS='new', ACTION='write')
 DO n=1, nx
  WRITE(uni, *) x(n), phi00(n), p_plot(:,n), u_plot(:, n)
 END DO 

 CLOSE(uni)

 !Other scheme:
 c = 1 !Here I was testing for many dt
 dt = c * dx / SQRT(h)
 t = t0
 nn = 0
 
 u_now = 0.0
 u_new = 0.0
 p_new = 0.0
 p_plot = 0.0
 u_plot = 0.0
 p_now = phi0(x)
 
 DO
  IF (t>t1) EXIT
  t = t + dt

  u_new = ueq_ft(u_now, p_now, dtdx)
  p_new = peq_ft(u_new, p_now, c)
  
  p_now = p_new
  u_now = u_new

  IF (MOD(t,tp) < dt) THEN
   nn=nn+1
   IF (nn<=nt) THEN
    p_plot(nn, :) = p_now 
    u_plot(nn, :) = u_now    
   END IF
  END IF

 END DO
  
 OPEN(UNIT=uni, IOSTAT=ios, FILE='data8_scheme2.dat', STATUS='new', ACTION='write')
 DO n=1, nx
  WRITE(uni, *) x(n), phi00(n), p_plot(:,n), u_plot(:, n)
 END DO 

 CLOSE(uni)

 PRINT*, dt
 END PROGRAM ass8

