!Program to integrate linear ODE using the Forward scheme and the Heun's method

!Module with functions:

MODULE ode_integration
 IMPLICIT NONE
 CONTAINS

 !Function:
 FUNCTION funcion(x,y) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x, y
  REAL :: res

  res= -0.5*y + 4.0 * EXP(-0.5*x) * COS(4.0*x)

 END FUNCTION funcion

 FUNCTION analytical(x) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x
  REAL :: res

  res = EXP(-0.5*x) * SIN(4.0*x)

 END FUNCTION analytical

 FUNCTION euler(x,y,dx) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x, y, dx
  REAL :: res

  res = y + dx * funcion(x,y)
  
 END FUNCTION euler

 FUNCTION heun(x,y,dx) RESULT(res)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x, y, dx
  REAL :: res, f0, ystar, xp

  f0 = funcion(x,y)
  ystar = y + dx * f0
  xp = x + dx
  res = y + 0.5 * dx * (f0 + funcion(xp, ystar))

 END FUNCTION heun

END MODULE ode_integration

!Begin program:
PROGRAM ass1
 USE ode_integration
 IMPLICIT NONE
 REAL :: x0, x1, dx, y0, e0, y_eu, y_he, x, y_sol
 REAL, DIMENSION(:), ALLOCATABLE :: sol_e_eu, sol_e_he, x_vec
 INTEGER :: nx, n, u, ios
 u=10

 x0 = 0.0
 x1 = 10.0
 dx = 0.1
 y0 = 0.0
 e0 = 0.0
 nx = ANINT((x1-x0)/dx) + 1
 ALLOCATE(sol_e_eu(nx), sol_e_he(nx), x_vec(nx))
 y_eu = y0
 y_he = y0
 n = 1
 x = x0
  
 DO
  y_sol = analytical(x)
  sol_e_eu(n) = y_eu - y_sol
  sol_e_he(n) = y_he - y_sol
  x_vec(n) = x

  IF (x>(x1-dx/2)) EXIT
   
  n = n + 1
  y_eu = euler(x,y_eu,dx)
  y_he = heun(x,y_he,dx)
  x = x + dx
  
 END DO

 !Write data for plotting
 OPEN(UNIT=u, IOSTAT=ios, FILE='data.dat', STATUS='new', ACTION='write')
 DO n=1, nx
  WRITE(u, *) x_vec(n), sol_e_eu(n), sol_e_he(n)
 END DO

 CLOSE(u)

END PROGRAM ass1

