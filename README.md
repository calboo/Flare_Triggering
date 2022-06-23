# Flare Triggering

To simulate the triggering of solar flares the MHD simulation code Lare3d was used. The Lare3d code, manual and IDL visualisation scripts can be found here:  https://warwick.ac.uk/fac/sci/physics/research/cfsa/people/tda/larexd/.

In the simulations emerging bipolar regions (BRs) of magnetic field were injected into a domain with a non-potential sheared magnetic field arcade. The overlying field changes polarity across a polarity inversion line (PIL) and the injected field emerges in the centre of the domain along the PIL. The injected field is an ellipsoid that has opposite polarity to the overlying field and induces flaring through magnetic reconnection. In this project the effects of oscillating, twisting and colliding BRs on the strength of induced flares was considered. The project is described in the paper:

*The Effects of Oscillations and Collisions of Emerging Bipolar Regions on the Triggering of Solar Flares*, https://iopscience.iop.org/article/10.3847/1538-4357/aba61a,

This repository contains all the files needed to reproduce the simulations described and the figures shown in this paper. 

## Simulation files

The folders Flare_trigger_linear, Flare_trigger_twist and Flare_trigger_collide each contain four files:

- control.f90
- initial_conditions.f90
- boundary.f90
- shared_data.F90

To utilise these files in Lare3d the first three should replace files of the same name in the src directory and the fourth should replace a file of the same name in the *core/* subdirectory which is found within the *src/* directory. Information on how to run Lare3d is given in the Lare3d manual which can be accessed from the webpage https://warwick.ac.uk/fac/sci/physics/research/cfsa/people/tda/larexd/.

### Flare_trigger_linear

These files are used to reproduce simulations with a single BR that oscillates both along and across the PIL. The parameters that can be changed in this simulation can all be found in the shared_data.F90 file, under flare triggering parameters, and are:

| Parameter | Description |
| --- | --- |
| theta0 | the angle at which the overlying field is sheared relative to the PIL|
| phif | the angle at which the injected field is oriented relative to the PIL|
| vzf | the vertical velocity at which the plasma within the injected field rises|
| z0 | the initial vertical coordinate of the centre of the injected field|
| r0 | the radius of the injected field ellipsoid perpendicular to the PIL|
| a0 | the radius of injected field ellipsoid parallel to the PIL|
| t1 | the time at which the injected field stops rising|
| t0 | the time at which the injected field begins rising|
| B_e | the magnetic field strength within the injected field|
| x_amp | the oscillation amplitude of the injected field along the PIL|
| y_amp | the oscillation amplitude of the injected field across the PIL|
| omega | the oscillation frequency of the injected field |

Note that all parameters are given in normalised units, see Lare3d manual for more info.

### Flare_trigger_twist

These files are used to reproduce simulations with a single BR that oscillates torsionally but remains in the centre of the domain. The parameters that can be changed in this simulation can all be found in the shared_data.F90 file and are the same as those for Flare_trigger_linear but with no x_amp and y_amp and instead a parameter for the amplitude of torsional oscillations, rotamp0, which is the maximum rotation of the injected field in degrees.

### Flare_trigger_collide

These files are used to reproduce simulations with a two BR that collide over the course of the simulation. The parameters that can be changed in this simulation can all be found in the shared_data.F90 file and are the same as those for Flare_trigger_linear but with no x_amp and y_amp and instead a parameter for the collision velocity coll_vel.

Additionally the injected fields can be set to be either corotating, counter rotating or non-rotating by selecting the appropriate lines in both the initial_condition.f90 and boundary.f90 files.

## Visualisation files

The visualisation script files are all located in the Visualisation_scripts folder. These scripts are named according to the figure in https://iopscience.iop.org/article/10.3847/1538-4357/aba61a that they were used to produce.

The scripts are all written in IDL and require the Start.pro script distributed with Lare3d in order to run, in addition to the SDF and IDL directories distributed with Lare3d that allow the Lare3d output files to be read within IDL. The user must initiate idl with the command idl Start.pro and can then proceed to run the chosen visualisation script. 

The user must also be sure that the simulation data is available for the simulation and the data is being read from the right location within the visualisation script.

