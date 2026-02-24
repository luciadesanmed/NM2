!Assignment 6:
!Use the implicit scheme to solve the diffusion problem following the TDMA algorithm.
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

END MODULE diffusion_mod

!Begin program:
PROGRAM ass6
 USE diffusion_mod
 IMPLICIT NONE
 REAL :: x0, x1, dx, k, dt, t, t0, t1, tp, bc0, bc1, alpha
 REAL, DIMENSION(:), ALLOCATABLE :: x, phi_now, phi00, phi_new, a_vec, b_vec, c_vec, f_vec, delta
 REAL, DIMENSION(:,:), ALLOCATABLE :: phi_plot
 INTEGER :: nx, n, ios, uni, nn, j
 uni=10

 nn = 0
 x0 = 0.0
 x1 = 1.0
 dx = 0.01
 k = 2.9e-5
 t0 = 0.0
 t1 = 21600.0
 dt = 360.0
 alpha = k*dt/(dx**2)
 tp = 3600.0
 bc0 = 273.15
 bc1 = 273.15
 nx = ANINT((x1-x0)/dx) + 1
 
 ALLOCATE(x(nx), phi00(nx), phi_now(nx), phi_new(nx), phi_plot(6,nx) &
         , a_vec(nx-2), b_vec(nx-1), c_vec(nx-2), f_vec(nx-1), delta(nx-1))

 x = xvec(x0,x1,dx)
 a_vec = a_vec - alpha
 c_vec = c_vec - alpha
 b_vec = b_vec + 1.0+2.0*alpha
 f_vec = x(1:nx-1)
 delta = x(1:nx-1)
 phi00 = phi0(x)
 phi_now = phi0(x)
 t =  t0
 
 DO
  
  IF (t>t1) EXIT
  
  f_vec(1) = 0.0
  delta(1) = bc0

  DO j=1, SIZE(c_vec)  
   f_vec(j+1) = c_vec(j) / (b_vec(j+1)-a_vec(j)*f_vec(j))
   delta(j+1) = (phi_now(j+1)-a_vec(j)*delta(j)) / (b_vec(j+1)-a_vec(j)*f_vec(j))
  END DO
 
   phi_new(nx) = bc1

  DO j= SIZE(f_vec), 1, -1
   phi_new(j) = delta(j) - f_vec(j)*phi_new(j+1)
  END DO

  phi_now = phi_new
  t =  t + dt

   IF (MOD(t,tp) < dt) THEN
    nn=nn+1

    IF (nn<=6) THEN
     phi_plot(nn,:) = phi_now
    END IF

   END IF

  END DO

 !Write data for plotting:
 OPEN(UNIT=uni, IOSTAT=ios, FILE='data6.dat', STATUS='new', ACTION='write')
 DO n=1, nx
  WRITE(uni, *) x(n), phi00(n), phi_plot(1,n), phi_plot(2,n), phi_plot(3,n), phi_plot(4,n), &
          phi_plot(5,n), phi_plot(6,n)
 END DO 

 CLOSE(uni)

END PROGRAM ass6

