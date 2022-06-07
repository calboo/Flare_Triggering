MODULE initial_conditions

  USE shared_data
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

CONTAINS

  !****************************************************************************
  ! This function sets up the initial condition for the code
  ! The variables which must be set are:
  !   rho - density
  !   v{x,y,z} - Velocities in x, y, z
  !   b{x,y,z} - Magnetic fields in x, y, z
  !   energy - Specific internal energy
  !
  !****************************************************************************

  SUBROUTINE set_initial_conditions

  INTEGER :: ix, iy, iz
  REAL(num) :: beta, alpha, k, kk

! Gravity and plasma beta 
  grav=0.0_num
  beta = 1.0e-8_num 

! Define parameters for overlying field
  theta0 = pi* (theta0/180.0)
  k = 3*pi
  alpha = k*tan(theta0)/sqrt((tan(theta0)**2.0)+1.0)
  kk = sqrt(k**2.0-alpha**2.0)

! Velocities
! Static domain.

  vx = 0.0_num
  vy = 0.0_num
  vz = 0.0_num

! Define overlying Sheared Magnetic Field

  DO ix = -2, nx+2
     DO iy = -1, ny+2
         DO iz = -1, nz+2
            bx(ix,iy,iz)= alpha*exp(-kk*zc(iz))*cos(k*yc(iy))/k
         END DO
     END DO
  END DO

  DO ix = -1, nx+2
     DO iy = -2, ny+2
         DO iz = -1, nz+2
            by(ix,iy,iz)= -kk*exp(-kk*zc(iz))*cos(k*yb(iy))/k
         END DO
     END DO
  END DO

  DO ix = -1, nx+2
     DO iy = -1, ny+2
         DO iz = -2, nz+2
            bz(ix,iy,iz)= exp(-kk*zb(iz))*sin(k*yc(iy))
         END DO 
     END DO
  END DO

! Density

  DO iy= -1,ny+2 
     DO iz = -1,nz+2 
         DO ix = -1,nx+2 
             rho(ix,iy,iz) = 1.0
          END DO
     END DO
  END DO

! Energy
! The energy has been set such that the pressure will equal beta/2 everywhere.

  energy=0.5_num*(beta*1.0_num) / ((rho) * (gamma-1.0_num))

! Store Initial Values
! Initial values are stored in these arrays.

  ALLOCATE(rho0(-1:nx+2, -1:ny+2, -1:nz+2))
  ALLOCATE(bx0 (-2:nx+2, -1:ny+2, -1:nz+2))
  ALLOCATE(by0 (-1:nx+2, -2:ny+2, -1:nz+2))
  ALLOCATE(bz0 (-1:nx+2, -1:ny+2, -2:nz+2))
  ALLOCATE(energy0(-1:nx+2, -1:ny+2, -1:nz+2))

  bx0 = bx
  by0 = by
  bz0 = bz
  rho0 = rho
  energy0 = energy

! Allocate arrays for injected field coordinates

  ALLOCATE(ctv (-2:nx+2, -2:ny+2, -2: nz+2))
  ALLOCATE(ellv(-2:nx+2, -2:ny+2, -2: nz+2))
  ALLOCATE(ctx (-2:nx+2, -1:ny+2, -1: nz+2))
  ALLOCATE(ellx(-2:nx+2, -1:ny+2, -1: nz+2))
  ALLOCATE(cty (-1:nx+2, -2:ny+2, -1: nz+2))
  ALLOCATE(elly(-1:nx+2, -2:ny+2, -1: nz+2))
  ALLOCATE(ctz (-1:nx+2, -1:ny+2, -2: nz+2))
  ALLOCATE(ellz(-1:nx+2, -1:ny+2, -2: nz+2))

! Convert phif into radians phic

    phic = (phif*pi/180.0_num)

! Convert torsional amplitude into radians

    rotamp = (rotamp0*pi/180.0_num)
    
! Set up coordinate system for the injected field
! This is used to define where the injected field is
! and the magnetic field components within it.

! Note that the torsional oscillations are incorporated
! directly into the coordinate system transformations.

! The coordinates system must be defined four times
! due to the staggering of the magnetic field grid.

! Velocity aligned injected field coordinates

    DO ix = -2,nx+2
       DO iy = -2,ny+2
          DO iz = -2, 0
            xx = xb(ix)*cos(phic-rotamp*sin(2.0_num*pi*omega*time)) + &
                 yb(iy)*sin(phic-rotamp*sin(2.0_num*pi*omega*time))
            yy =-xb(ix)*sin(phic-rotamp*sin(2.0_num*pi*omega*time)) + &
                 yb(iy)*cos(phic-rotamp*sin(2.0_num*pi*omega*time))
            zz = zb(iz)-(z0 + vzf*min(t1-t0,time-t0))
            cr = sqrt(yy**2.0+zz**2.0)
            cz = xx
            ctv(ix,iy,iz) = atan2(zz,yy)
            ellv(ix,iy,iz) = (cr/r0)**2.0 + (cz/a0)**2.0
          END DO
        END DO
    END DO

! Bx aligned injected field coordinates

    DO ix = -2,nx+2
       DO iy = -1,ny+2
          DO iz = -1, nz+2
            xx = xb(ix)*cos(phic-rotamp*sin(2.0_num*pi*omega*time)) + &
                 yc(iy)*sin(phic-rotamp*sin(2.0_num*pi*omega*time))
            yy =-xb(ix)*sin(phic-rotamp*sin(2.0_num*pi*omega*time)) + &
                 yc(iy)*cos(phic-rotamp*sin(2.0_num*pi*omega*time))
            zz = zc(iz)-(z0 + vzf*min(t1-t0,time-t0))
            cr = sqrt(yy**2.0+zz**2.0)
            cz = xx
            ctx(ix,iy,iz) = atan2(zz,yy)
            ellx(ix,iy,iz) = (cr/r0)**2.0 + (cz/a0)**2.0
          END DO
        END DO
    END DO

! By aligned injected field coordinates

    DO ix = -1,nx+2
       DO iy = -2,ny+2
          DO iz = -1, nz+2
            xx = xc(ix)*cos(phic-rotamp*sin(2.0_num*pi*omega*time)) + &
                 yb(iy)*sin(phic-rotamp*sin(2.0_num*pi*omega*time))
            yy =-xc(ix)*sin(phic-rotamp*sin(2.0_num*pi*omega*time)) + &
                 yb(iy)*cos(phic-rotamp*sin(2.0_num*pi*omega*time))
            zz = zc(iz)-(z0 + vzf*min(t1-t0,time-t0))
            cr = sqrt(yy**2.0+zz**2.0)
            cz = xx
            cty(ix,iy,iz) = atan2(zz,yy)
            elly(ix,iy,iz) = (cr/r0)**2.0 + (cz/a0)**2.0
          END DO
        END DO
    END DO

! Bz aligned injected field coordinates

    DO ix = -1,nx+2
       DO iy = -1,ny+2
          DO iz = -2, nz+2
            xx = xc(ix)*cos(phic-rotamp*sin(2.0_num*pi*omega*time)) + &
                 yc(iy)*sin(phic-rotamp*sin(2.0_num*pi*omega*time))
            yy =-xc(ix)*sin(phic-rotamp*sin(2.0_num*pi*omega*time)) + &
                 yc(iy)*cos(phic-rotamp*sin(2.0_num*pi*omega*time))
            zz = zb(iz)-(z0 + vzf*min(t1-t0,time-t0))
            cr = sqrt(yy**2.0+zz**2.0)
            cz = xx
            ctz(ix,iy,iz) = atan2(zz,yy)
            ellz(ix,iy,iz) = (cr/r0)**2.0 + (cz/a0)**2.0
          END DO
        END DO
    END DO

! Define the injected field at simulation start time

  DO ix = -2, nx+2
     DO iy = -1, ny+2
        DO iz = -1, nz+2
          IF (ellx(ix,iy,iz) .LE. 1.0) THEN
            bx(ix,iy,iz) = B_e*sin(phic)*sin(ctx(ix,iy,iz))
          ELSE
            bx(ix,iy,iz) = bx0(ix,iy,iz)
          END IF
        END DO
     END DO
  END DO

  DO ix = -1, nx+2
     DO iy = -2, ny+2
         DO iz = -1, nz+2
          IF (elly(ix,iy,iz) .LE. 1.0) THEN
            by(ix,iy,iz) = - B_e*cos(phic)*sin(cty(ix,iy,iz))
          ELSE
            by(ix,iy,iz) = by0(ix,iy,iz)
          END IF
        END DO
     END DO
  END DO

  DO ix = -1, nx+2
     DO iy = -1, ny+2
         DO iz = -2, nz+2
          IF (ellz(ix,iy,iz) .LE. 1.0) THEN
            bz(ix,iy,iz) = B_e*cos(ctz(ix,iy,iz))
          ELSE
            bz(ix,iy,iz) = bz0(ix,iy,iz)
          END IF
         END DO
     END DO
  END DO

 END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
