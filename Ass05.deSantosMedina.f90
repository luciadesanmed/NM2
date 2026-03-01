!Assignment 5:
!Use the forward difference scheme to solve the diffusion problem.
!Use a spatial resolution of dx=0.01 with diff coeff K=2.92e-5
!integrate for at least 6 h and show the sol every hour.
!IC: 273.15+20x+sin(50pix) for 0<x<0.5
!273.15+20-20x+sin(50pix) for 0.5<x<1
!BC: bc0=273.15, bc1=273.15

!Module with functions:

MODULE diffusion_mod
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

  x = x0
  DO n = 1, nx
   res(n) = x
   x = x + dx
  END DO
  
 END FUNCTION xvec

 !Initial condition:
 FUNCTION phi0(x) RESULT(res)
  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(IN) :: x
  INTEGER :: n
  REAL :: pi
  REAL, DIMENSION(:), ALLOCATABLE :: res
   
  pi = 4.0*ATAN(1.0)

  ALLOCATE(res(SIZE(x)))

  DO n= 1, SIZE(x)
   IF (0 <= x(n) .AND. x(n) <= 0.5) THEN
     res(n) = 273.15 + 20.0*x(n) + SIN(50.0*pi*x(n))
   ELSE IF (0.5 < x(n) .AND. x(n) <= 1) THEN 
    res(n) = 273.15 + 20.0 - 20.0*x(n) + SIN(50.0*pi*x(n))
   END IF 
  END DO

 END FUNCTION phi0

 !Forward scheme for diffusion:
 SUBROUTINE diffusion(phi_now, phi_new, bc0, bc1, k, dt, dx)
  IMPLICIT NONE
  REAL, DIMENSION(:), INTENT(INOUT) :: phi_new
  REAL, DIMENSION(:), INTENT(IN) :: phi_now
  REAL, INTENT(IN) :: bc0, bc1, k, dt, dx
  REAL :: mult
  INTEGER :: n, sz

  sz=SIZE(phi_now)

  mult = k* dt/(dx**2)
  phi_new(1) = bc0
  phi_new(sz) = bc1
  
  DO n= 2, sz-1
   phi_new(n) = phi_now(n) + mult * (phi_now(n+1) - 2.0 * phi_now(n) + phi_now(n-1))
  END DO

 END SUBROUTINE diffusion

END MODULE diffusion_mod

!Begin program:
PROGRAM ass5
 USE diffusion_mod
 IMPLICIT NONE
 REAL :: x0, x1, dx, k, dt, t, t0, t1, tp, bc0, bc1
 REAL, DIMENSION(:), ALLOCATABLE :: x, phi_now, phi00, phi_new
 REAL, DIMENSION(:,:), ALLOCATABLE :: phi_plot
 INTEGER :: nt, nx, n, ios, uni, nn
 uni=10

 nn = 0
 x0 = 0.0
 x1 = 1.0
 dx = 0.01
 k = 2.9e-5
 t0 = 0.0
 t1 = 21600.0
 dt = (dx**2)/(2.0*k)
 tp = 3600.0
 bc0 = 273.15
 bc1 = 273.15
 nx = ANINT((x1-x0)/dx) + 1
 nt = ANINT(t1/tp)
 
 ALLOCATE(x(nx), phi00(nx), phi_now(nx), phi_new(nx), phi_plot(nt,nx))

 x = xvec(x0,x1,dx)
 phi00 = phi0(x)
 phi_now = phi0(x)
 t =  t0

 DO 
  IF (t>t1) EXIT
  
  CALL diffusion(phi_now,phi_new,bc0,bc1,k,dt,dx)
  
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
 !OPEN(UNIT=uni, IOSTAT=ios, FILE='data5.dat', STATUS='new', ACTION='write')
 !DO n=1, nx
 ! WRITE(uni, *) x(n), phi00(n), phi_plot(:,n)
 !END DO 

 !CLOSE(uni)

END PROGRAM ass5

