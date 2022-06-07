# Flare_Triggering

To simulate the triggering of solar flares the MHD code Lare3d, https://warwick.ac.uk/fac/sci/physics/research/cfsa/people/tda/larexd/, was used. Emerging bipolar regions (BRs) of magnetic field were injected into a domain with a non-potential sheared magnetic field arcade. The overlying field changes polarity across a polarity inversion line (PIL) and the injected field emerges in the centre of the domain along the PIL. The injected field is an ellipsoid and has opposite polarity to the overlyying field and induces flaring through magnetic reconnection.

In this project which is described in this paper, https://iopscience.iop.org/article/10.3847/1538-4357/aba61a, the effects of oscillating, twisting and colliding BRs on the strength of induced flares was considered. This repository contains all the files needed to reproduce the simulations described and the figures shown in this paper. 

## Simulation files

The folders Flare_trigger_linear, Flare_trigger_twist and Flare_trigger_collide each contain four files:

- control.f90
- initial_conditions.f90
- boundary.f90
- shared_data.F90

To utilise these files in Lare3d the first three should replace files of the same name in the src directory and the fourth should replace a file of the same name in the core subdirectory which is found within the src directory.

### Flare_trigger_linear

These files are used to reproduce simulations with a single BR that oscillates both along and across the PIL. The parameters that can be changed in this simulation can all be found in the shared_data.F90 file and are:

- theta0 - the angle at which the overlying field is sheared relative to the PIL
- phif - the angle at which the injected field is oriented relative to the PIL
- vzf - the vertical velocity at which the plasma within the injected field rises
- z0 - the initial vertical coordinate of the centre of the injected field
- r0 - the radius of the injected field ellipsoid perpendicular to the PIL
- a0 - the radius of injected field ellipsoid parallel to the PIL
- t1 - the time at which the injected field stops rising
- t0 - the time at which the injected field begins rising
- B_e - the magnetic field strength within the injected field
- x_amp - the oscillation amplitude of the injected field along the PIL
- y_amp - the oscillation amplitude of the injected field across the PIL
- omega - the oscillation frequency of the injected field 


