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
  REAL(num) :: beta, alpha, k, kk,drivx,drivy

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

  ALLOCATE(ctv1 (-2:nx+2, -2:ny+2, -2: nz+2))
  ALLOCATE(ellv1(-2:nx+2, -2:ny+2, -2: nz+2))
  ALLOCATE(ctx1 (-2:nx+2, -1:ny+2, -1: nz+2))
  ALLOCATE(ellx1(-2:nx+2, -1:ny+2, -1: nz+2))
  ALLOCATE(cty1 (-1:nx+2, -2:ny+2, -1: nz+2))
  ALLOCATE(elly1(-1:nx+2, -2:ny+2, -1: nz+2))
  ALLOCATE(ctz1 (-1:nx+2, -1:ny+2, -2: nz+2))
  ALLOCATE(ellz1(-1:nx+2, -1:ny+2, -2: nz+2))

  ALLOCATE(ctv2 (-2:nx+2, -2:ny+2, -2: nz+2))
  ALLOCATE(ellv2(-2:nx+2, -2:ny+2, -2: nz+2))
  ALLOCATE(ctx2 (-2:nx+2, -1:ny+2, -1: nz+2))
  ALLOCATE(ellx2(-2:nx+2, -1:ny+2, -1: nz+2))
  ALLOCATE(cty2 (-1:nx+2, -2:ny+2, -1: nz+2))
  ALLOCATE(elly2(-1:nx+2, -2:ny+2, -1: nz+2))
  ALLOCATE(ctz2 (-1:nx+2, -1:ny+2, -2: nz+2))
  ALLOCATE(ellz2(-1:nx+2, -1:ny+2, -2: nz+2))

  ! The injected fields can be set to be either corotating,
  ! counter rotating or non-rotating.

  ! Non-rotating
  phic1 = (phif*pi/180.0_num) 
  phic2 = (phif*pi/180.0_num)

  ! Corotating
  ! phic1 = 2.0_num*pi*time
  ! phic2 = 2.0_num*pi*time

  ! Counter rotating
  ! phic1 =  2.0_num*pi*time
  ! phic2 = -2.0_num*pi*time

  ! The coordinates system must be defined four times
  ! due to the staggering of the magnetic field grid.

  ! The injected fields are offset in the direction perpendicular
  ! to the PIL, in opposite direction, by the radius of the injected
  ! field regions divided by 20.

  ! The injected fields are initially offset by 20 Mm.

  ! Velocity aligned injected field coordinates

    DO ix = -2,nx+2
       DO iy = -2,ny+2
          DO iz = -2, nz+2
            xx1 = (xb(ix)-1.0_num)*cos(phic1) + (yb(iy)-(r0/20.0_num))*sin(phic1)
            xx2 = (xb(ix)+1.0_num)*cos(phic2) + (yb(iy)+(r0/20.0_num))*sin(phic2)
            yy1 =-(xb(ix)-1.0_num)*sin(phic1) + (yb(iy)-(r0/20.0_num))*cos(phic1)
            yy2 =-(xb(ix)+1.0_num)*sin(phic2) + (yb(iy)+(r0/20.0_num))*cos(phic2)
            zz = zb(iz)-(z0 + vzf*min(t1-t0,time-t0))
            cr1 = sqrt(yy1**2.0+zz**2.0)
            cr2 = sqrt(yy2**2.0+zz**2.0)
            cz1 = xx1
            cz2 = xx2
            ctv1(ix,iy,iz) = atan2(zz,yy1)
            ctv2(ix,iy,iz) = atan2(zz,yy2)
            ellv1(ix,iy,iz) = (cr1/r0)**2.0 + (cz1/a0)**2.0
            ellv2(ix,iy,iz) = (cr2/r0)**2.0 + (cz2/a0)**2.0
          END DO
        END DO
    END DO

    ! Bx aligned injected field coordinates

    DO ix = -2,nx+2
       DO iy = -1,ny+2
          DO iz = -1, nz+2
            xx1 = (xb(ix)-1.0_num)*cos(phic1) + (yc(iy)-(r0/20.0_num))*sin(phic1)
            xx2 = (xb(ix)+1.0_num)*cos(phic2) + (yc(iy)+(r0/20.0_num))*sin(phic2)
            yy1 =-(xb(ix)-1.0_num)*sin(phic1) + (yc(iy)-(r0/20.0_num))*cos(phic1)
            yy2 =-(xb(ix)+1.0_num)*sin(phic2) + (yc(iy)+(r0/20.0_num))*cos(phic2)
            zz = zc(iz)-(z0 + vzf*min(t1-t0,time-t0))
            cr1 = sqrt(yy1**2.0+zz**2.0)
            cr2 = sqrt(yy2**2.0+zz**2.0)
            cz1 = xx1
            cz2 = xx2
            ctx1(ix,iy,iz) = atan2(zz,yy1)
            ctx2(ix,iy,iz) = atan2(zz,yy2)
            ellx1(ix,iy,iz) = (cr1/r0)**2.0 + (cz1/a0)**2.0
            ellx2(ix,iy,iz) = (cr2/r0)**2.0 + (cz2/a0)**2.0
          END DO
        END DO
    END DO

    ! By aligned injected field coordinates

    DO ix = -1,nx+2
       DO iy = -2,ny+2
          DO iz = -1, nz+2
            xx1 = (xc(ix)-1.0_num)*cos(phic1) + (yb(iy)-(r0/20.0_num))*sin(phic1)
            xx2 = (xc(ix)+1.0_num)*cos(phic2) + (yb(iy)+(r0/20.0_num))*sin(phic2)
            yy1 =-(xc(ix)-1.0_num)*sin(phic1) + (yb(iy)-(r0/20.0_num))*cos(phic1)
            yy2 =-(xc(ix)+1.0_num)*sin(phic2) + (yb(iy)+(r0/20.0_num))*cos(phic2)
            zz = zc(iz)-(z0 + vzf*min(t1-t0,time-t0))
            cr1 = sqrt(yy1**2.0+zz**2.0)
            cr2 = sqrt(yy2**2.0+zz**2.0)
            cz1 = xx1
            cz2 = xx2
            cty1(ix,iy,iz) = atan2(zz,yy1)
            cty2(ix,iy,iz) = atan2(zz,yy2)
            elly1(ix,iy,iz) = (cr1/r0)**2.0 + (cz1/a0)**2.0
            elly2(ix,iy,iz) = (cr2/r0)**2.0 + (cz2/a0)**2.0
          END DO
        END DO
    END DO

    ! Bz aligned injected field coordinates

    DO ix = -1,nx+2
       DO iy = -1,ny+2
          DO iz = -2, nz+2
            xx1 = (xc(ix)-1.0_num)*cos(phic1) + (yc(iy)-(r0/20.0_num))*sin(phic1)
            xx2 = (xc(ix)+1.0_num)*cos(phic2) + (yc(iy)+(r0/20.0_num))*sin(phic2)
            yy1 =-(xc(ix)-1.0_num)*sin(phic1) + (yc(iy)-(r0/20.0_num))*cos(phic1)
            yy2 =-(xc(ix)+1.0_num)*sin(phic2) + (yc(iy)+(r0/20.0_num))*cos(phic2)
            zz = zb(iz)-(z0 + vzf*min(t1-t0,time-t0))
            cr1 = sqrt(yy1**2.0+zz**2.0)
            cr2 = sqrt(yy2**2.0+zz**2.0)
            cz1 = xx1
            cz2 = xx2
            ctz1(ix,iy,iz) = atan2(zz,yy1)
            ctz2(ix,iy,iz) = atan2(zz,yy2)
            ellz1(ix,iy,iz) = (cr1/r0)**2.0 + (cz1/a0)**2.0
            ellz2(ix,iy,iz) = (cr2/r0)**2.0 + (cz2/a0)**2.0
          END DO
        END DO
    END DO

! Define the injected fields at simulation start time

  DO ix = -2, nx+2
     DO iy = -1, ny+2
        DO iz = -1, nz+2
          IF (ellx1(ix,iy,iz) .LE. 1.0) THEN
             bx(ix,iy,iz) = B_e*sin(phic1)*sin(ctx1(ix,iy,iz))
          ELSE IF (ellx2(ix,iy,iz) .LE. 1.0) THEN
             bx(ix,iy,iz) = B_e*sin(phic2)*sin(ctx2(ix,iy,iz))
          ELSE
             bx(ix,iy,iz) = bx0(ix,iy,iz)
          END IF
        END DO
     END DO
  END DO

  DO ix = -1, nx+2
     DO iy = -2, ny+2
         DO iz = -1, nz+2
          IF (elly1(ix,iy,iz) .LE. 1.0) THEN
             by(ix,iy,iz) = - B_e*cos(phic1)*sin(cty1(ix,iy,iz))
          ELSE IF (elly2(ix,iy,iz) .LE. 1.0) THEN
             by(ix,iy,iz) = - B_e*cos(phic2)*sin(cty2(ix,iy,iz))
          ELSE
            by(ix,iy,iz) = by0(ix,iy,iz)
          END IF
        END DO
     END DO
  END DO

  DO ix = -1, nx+2
     DO iy = -1, ny+2
         DO iz = -2, nz+2
          IF (ellz1(ix,iy,iz) .LE. 1.0) THEN
             bz(ix,iy,iz) = B_e*cos(ctz1(ix,iy,iz))
          ELSE IF (ellz2(ix,iy,iz) .LE. 1.0) THEN
             bz(ix,iy,iz) = B_e*cos(ctz2(ix,iy,iz))
          ELSE
             bz(ix,iy,iz) = bz0(ix,iy,iz)
          END IF
         END DO
     END DO
  END DO

 END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
