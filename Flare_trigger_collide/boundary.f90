!******************************************************************************
! This module contains the boundary conditions for the entire code
! Any new boundary conditions should be added here
!******************************************************************************

MODULE boundary

  USE shared_data
  USE mpiboundary

  IMPLICIT NONE

CONTAINS

  !****************************************************************************
  ! Set up any necessary variables for the chosen boundary conditions
  !****************************************************************************

  SUBROUTINE set_boundary_conditions

    any_open = .FALSE.
    IF (xbc_min == BC_OPEN .OR. xbc_max == BC_OPEN &
        .OR. ybc_min == BC_OPEN .OR. ybc_max == BC_OPEN &
        .OR. zbc_min == BC_OPEN .OR. zbc_max == BC_OPEN) any_open = .TRUE.

  END SUBROUTINE set_boundary_conditions



  !****************************************************************************
  ! Call all of the boundaries needed by the core Lagrangian solver
  !****************************************************************************

  SUBROUTINE boundary_conditions

  ! Set up coordinate system for the injected fields
  ! This is used to define where the injected fields are
  ! and the magnetic field components within them.

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

  ! The injected fields will move towards each other with velocity coll_vel.

  ! Injected field displacements
  disp1 =-1.0_num+time*coll_vel 
  disp2 = 1.0_num-time*coll_vel

  ! Velocity aligned injected field coordinates

    DO ix = -2,nx+2
       DO iy = -2,ny+2
          DO iz = -2, 0
            xx1 = (xb(ix)+disp1)*cos(phic1) + (yb(iy)-(r0/20.0_num))*sin(phic1)
            xx2 = (xb(ix)+disp2)*cos(phic2) + (yb(iy)+(r0/20.0_num))*sin(phic2)
            yy1 =-(xb(ix)+disp1)*sin(phic1) + (yb(iy)-(r0/20.0_num))*cos(phic1)
            yy2 =-(xb(ix)+disp2)*sin(phic1) + (yb(iy)+(r0/20.0_num))*cos(phic2)
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
          DO iz = -1, 0
            xx1 = (xb(ix)+disp1)*cos(phic1) + (yc(iy)-(r0/20.0_num))*sin(phic1)
            xx2 = (xb(ix)+disp2)*cos(phic2) + (yc(iy)+(r0/20.0_num))*sin(phic2)
            yy1 =-(xb(ix)+disp1)*sin(phic1) + (yc(iy)-(r0/20.0_num))*cos(phic1)
            yy2 =-(xb(ix)+disp2)*sin(phic2) + (yc(iy)+(r0/20.0_num))*cos(phic2)
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
          DO iz = -1, 0
            xx1 = (xc(ix)+disp1)*cos(phic1) + (yb(iy)-(r0/20.0_num))*sin(phic1)
            xx2 = (xc(ix)+disp2)*cos(phic2) + (yb(iy)+(r0/20.0_num))*sin(phic2)
            yy1 =-(xc(ix)+disp1)*sin(phic1) + (yb(iy)-(r0/20.0_num))*cos(phic1)
            yy2 =-(xc(ix)+disp2)*sin(phic2) + (yb(iy)+(r0/20.0_num))*cos(phic2)
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
          DO iz = -2, 0
            xx1 = (xc(ix)+disp1)*cos(phic1) + (yc(iy)-(r0/20.0_num))*sin(phic1)
            xx2 = (xc(ix)+disp2)*cos(phic2) + (yc(iy)+(r0/20.0_num))*sin(phic2)
            yy1 =-(xc(ix)+disp1)*sin(phic1) + (yc(iy)-(r0/20.0_num))*cos(phic1)
            yy2 =-(xc(ix)+disp2)*sin(phic2) + (yc(iy)+(r0/20.0_num))*cos(phic2)
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

    ! Call boundary condition subroutines

    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs
    CALL damp_boundaries

  END SUBROUTINE boundary_conditions


  !****************************************************************************
  ! Boundary conditions for magnetic field through plane
  !****************************************************************************

  SUBROUTINE bfield_bcs

    CALL bfield_mpi

    ! We have used reflective boundary conditions.
    ! At the lower boundary we define the injected magnetic fields.
    ! Where the injected fields overlap their magnetic fields
    ! are in linear superposition.

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      bx(-1,:,:) = bx(1,:,:)
      bx(-2,:,:) = bx(2,:,:)
      by( 0,:,:) = by(1,:,:)
      by(-1,:,:) = by(2,:,:)
      bz( 0,:,:) = bz(1,:,:)
      bz(-1,:,:) = bz(2,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      bx(nx+1,:,:) = bx(nx-1,:,:)
      bx(nx+2,:,:) = bx(nx-2,:,:)
      by(nx+1,:,:) = by(nx  ,:,:)
      by(nx+2,:,:) = by(nx-1,:,:)
      bz(nx+1,:,:) = bz(nx  ,:,:)
      bz(nx+2,:,:) = bz(nx-1,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      bx(:, 0,:) = bx(:,1,:)
      bx(:,-1,:) = bx(:,2,:)
      by(:,-1,:) = by(:,1,:)
      by(:,-2,:) = by(:,2,:)
      bz(:, 0,:) = bz(:,1,:)
      bz(:,-1,:) = bz(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      bx(:,ny+1,:) = bx(:,ny  ,:)
      bx(:,ny+2,:) = bx(:,ny-1,:)
      by(:,ny+1,:) = by(:,ny-1,:)
      by(:,ny+2,:) = by(:,ny-2,:)
      bz(:,ny+1,:) = bz(:,ny  ,:)
      bz(:,ny+2,:) = bz(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      bx(:,:, 0) = bx(:,:, 1)
      bx(:,:,-1) = bx(:,:, 2)
      by(:,:, 0) = by(:,:, 1)
      by(:,:,-1) = by(:,:, 2)
      bz(:,:,-1) = bz(:,:, 1)
      bz(:,:,-2) = bz(:,:, 2)
      DO ix = -2, nx+2
         DO iy = -1, ny+2
            DO iz = -1, 0
               IF (ellx1(ix,iy,iz) .LE. 1.0 .AND. ellx2(ix,iy,iz) .LE. 1.0) THEN
                  bx(ix,iy,iz) = B_e*sin(phic1)*sin(ctx1(ix,iy,iz)) &
                                 + B_e*sin(phic2)*sin(ctx2(ix,iy,iz))
               ELSE IF (ellx1(ix,iy,iz) .LE. 1.0) THEN
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
            DO iz = -1, 0
               IF (elly1(ix,iy,iz) .LE. 1.0 .AND. elly2(ix,iy,iz) .LE. 1.0) THEN
                  bx(ix,iy,iz) = - B_e*cos(phic1)*sin(cty1(ix,iy,iz)) &
                                 - B_e*cos(phic2)*sin(cty2(ix,iy,iz))
               ELSE IF (elly1(ix,iy,iz) .LE. 1.0) THEN
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
            DO iz = -2, 0
               IF (ellz1(ix,iy,iz) .LE. 1.0 .AND. ellz2(ix,iy,iz) .LE. 1.0) THEN
                  bz(ix,iy,iz) = B_e*cos(ctz1(ix,iy,iz)) &
                                 + B_e*cos(ctz2(ix,iy,iz))
               ELSE IF (ellz1(ix,iy,iz) .LE. 1.0) THEN
                  bz(ix,iy,iz) = B_e*cos(ctz1(ix,iy,iz))
               ELSE IF (ellz2(ix,iy,iz) .LE. 1.0) THEN
                  bz(ix,iy,iz) = B_e*cos(ctz2(ix,iy,iz))
               ELSE
                  bz(ix,iy,iz) = bz0(ix,iy,iz)
               END IF
            END DO
         END DO
      END DO
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      bx(:,:,nz+1) = bx(:,:,nz  )
      bx(:,:,nz+2) = bx(:,:,nz-1)
      by(:,:,nz+1) = by(:,:,nz  )
      by(:,:,nz+2) = by(:,:,nz-1)
      bz(:,:,nz+1) = bz(:,:,nz-1)
      bz(:,:,nz+2) = bz(:,:,nz-2)
    END IF

  END SUBROUTINE bfield_bcs


  !****************************************************************************
  ! Boundary conditions for specific internal energy
  !****************************************************************************

  SUBROUTINE energy_bcs

    CALL energy_mpi

    ! We reset the internal energy acroos the domain to the initial conditions
    ! to maintain a pseudo incompressible simulation.

    energy = energy0

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      energy( 0,:,:) = energy0( 0,:,:)
      energy(-1,:,:) = energy0(-1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      energy(nx+1,:,:) = energy0(nx+1,:,:)
      energy(nx+2,:,:) = energy0(nx+2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      energy(:, 0,:) = energy0(:, 0,:)
      energy(:,-1,:) = energy0(:,-1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      energy(:,ny+1,:) = energy0(:,ny+1,:)
      energy(:,ny+2,:) = energy0(:,ny+2,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      energy(:,:, 0) = energy0(:,:, 0)
      energy(:,:,-1) = energy0(:,:,-1)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      energy(:,:,nz+1) = energy0(:,:,nz+1)
      energy(:,:,nz+2) = energy0(:,:,nz+2)
    END IF

  END SUBROUTINE energy_bcs


  !****************************************************************************
  ! Boundary conditions for density
  !****************************************************************************

  SUBROUTINE density_bcs

    CALL density_mpi

    ! We reset the density acroos the domain to the initial conditions
    ! to maintain a pseudo incompressible simulation.

    rho = rho0

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      rho( 0,:,:) = rho0( 0,:,:)
      rho(-1,:,:) = rho0(-1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      rho(nx+1,:,:) = rho0(nx+1,:,:)
      rho(nx+2,:,:) = rho0(nx+2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      rho(:, 0,:) = rho0(:, 0,:)
      rho(:,-1,:) = rho0(:,-1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      rho(:,ny+1,:) = rho0(:,ny+1,:)
      rho(:,ny+2,:) = rho0(:,ny+2,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      rho(:,:, 0) = rho0(:,:, 0)
      rho(:,:,-1) = rho0(:,:,-1)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      rho(:,:,nz+1) = rho0(:,:,nz+1)
      rho(:,:,nz+2) = rho0(:,:,nz+2)
    END IF

  END SUBROUTINE density_bcs


  !****************************************************************************
  ! Boundary conditions for temperature
  !****************************************************************************

  SUBROUTINE temperature_bcs

    CALL density_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      temperature( 0,:,:) = temperature(1,:,:)
      temperature(-1,:,:) = temperature(2,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      temperature(nx+1,:,:) = temperature(nx  ,:,:)
      temperature(nx+2,:,:) = temperature(nx-1,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      temperature(:, 0,:) = temperature(:,1,:)
      temperature(:,-1,:) = temperature(:,2,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      temperature(:,ny+1,:) = temperature(:,ny  ,:)
      temperature(:,ny+2,:) = temperature(:,ny-1,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      temperature(:,:, 0) = temperature(:,:,1)
      temperature(:,:,-1) = temperature(:,:,2)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      temperature(:,:,nz+1) = temperature(:,:,nz  )
      temperature(:,:,nz+2) = temperature(:,:,nz-1)
    END IF

  END SUBROUTINE temperature_bcs



  !****************************************************************************
  ! Full timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE velocity_bcs

    ! We apply static boundary conditions for the velocity.
    ! At the lower boundary we apply vertical velocities where
    ! the injected fields emerge from the the photosphere.

    CALL velocity_mpi

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx(-2:0,:,:) = 0.0_num
      vy(-2:0,:,:) = 0.0_num
      vz(-2:0,:,:) = 0.0_num
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx(nx:nx+2,:,:) = 0.0_num
      vy(nx:nx+2,:,:) = 0.0_num
      vz(nx:nx+2,:,:) = 0.0_num
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx(:,-2:0,:) = 0.0_num
      vy(:,-2:0,:) = 0.0_num
      vz(:,-2:0,:) = 0.0_num
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx(:,ny:ny+2,:) = 0.0_num
      vy(:,ny:ny+2,:) = 0.0_num
      vz(:,ny:ny+2,:) = 0.0_num
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      vx(:,:,-2:0) = 0.0_num
      vy(:,:,-2:0) = 0.0_num 
      vz(:,:,-2:0) = 0.0_num 
      DO ix = -2,nx+2
        DO iy = -2,ny+2
          DO iz = -2, 0
             IF (ellv1(ix,iy,iz) .LE. 1.0 .AND. ellv2(ix,iy,iz) .LE. 1.0 &
                  .AND.(time .LE. t1)) THEN
                vz(ix,iy,iz) = 2.0_num*vzf
             ELSE IF ((ellv1(ix,iy,iz) .LE. 1.0) .AND. (time .LE. t1))  THEN
                vz(ix,iy,iz) = vzf
             ELSE IF ((ellv2(ix,iy,iz) .LE. 1.0) .AND. (time .LE. t1))  THEN
                vz(ix,iy,iz) = vzf
             END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      vx(:,:,nz:nz+2) = 0.0_num
      vy(:,:,nz:nz+2) = 0.0_num
      vz(:,:,nz:nz+2) = 0.0_num
    END IF

  END SUBROUTINE velocity_bcs


  !****************************************************************************
  ! Half timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE remap_v_bcs

    CALL remap_v_mpi

    ! We apply static boundary conditions for the half timestep velocity

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx1(-2:0,:,:) = 0.0_num
      vy1(-2:0,:,:) = 0.0_num
      vz1(-2:0,:,:) = 0.0_num
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx1(nx:nx+2,:,:) = 0.0_num
      vy1(nx:nx+2,:,:) = 0.0_num
      vz1(nx:nx+2,:,:) = 0.0_num
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx1(:,-2:0,:) = 0.0_num
      vy1(:,-2:0,:) = 0.0_num
      vz1(:,-2:0,:) = 0.0_num
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx1(:,ny:ny+2,:) = 0.0_num
      vy1(:,ny:ny+2,:) = 0.0_num
      vz1(:,ny:ny+2,:) = 0.0_num
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      vx1(:,:,-2:0) = 0.0_num
      vy1(:,:,-2:0) = 0.0_num 
      vz1(:,:,-2:0) = 0.0_num 
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      vx1(:,:,nz:nz+2) = 0.0_num
      vy1(:,:,nz:nz+2) = 0.0_num
      vz1(:,:,nz:nz+2) = 0.0_num
    END IF

  END SUBROUTINE remap_v_bcs


  !****************************************************************************
  ! Damped boundary conditions
  !****************************************************************************

  SUBROUTINE damp_boundaries

    ! Note that boundary damping is not applied in these simulations.

    REAL(num) :: a, d, pos, n_cells, damp_scale

    IF (.NOT.damping) RETURN
    ! number of cells near boundary to apply linearly increasing damping
    n_cells = 20.0_num 
    ! increase the damping if needed
    damp_scale = 1.0_num

    IF (proc_x_min == MPI_PROC_NULL) THEN
      d = n_cells * dxb(1)
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            pos = xb(ix) - x_min
            IF (pos < d) THEN
              a = dt * damp_scale * pos / d + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_x_max == MPI_PROC_NULL) THEN
      d = n_cells * dxb(nx)
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            pos = x_max - xb(ix)
            IF (pos < d) THEN
              a = dt * damp_scale * pos / d + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_y_min == MPI_PROC_NULL) THEN
      d = n_cells * dyb(1)
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            pos = yb(iy) - y_min
            IF (pos < d) THEN
              a = dt * damp_scale * pos / d + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_y_max == MPI_PROC_NULL) THEN
      d = n_cells * dyb(ny)
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            pos = y_max - yb(iy) 
            IF (pos < d) THEN
              a = dt * damp_scale * pos / d + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_z_min == MPI_PROC_NULL) THEN
      d = n_cells * dzb(1)
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            pos = zb(iz) - z_min
            IF (pos < d) THEN
              a = dt * damp_scale * pos / d + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (proc_z_max == MPI_PROC_NULL) THEN
      d = n_cells * dzb(nz)
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            pos = z_max - zb(iz)
            IF (pos < d) THEN
              a = dt * damp_scale * pos / d + 1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE damp_boundaries

END MODULE boundary
