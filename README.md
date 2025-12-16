# SUMMERSPH - Overview

This is a Summer Project because I am very, very bored, and want to look at circumstellar discs.

# Goals:
1. make an sph code that doesnt completely suck and runs on a laptop
2. speed up nearest neighbour searching and gravitational interactions with a barnes-hut tree.
3. include artificial viscosity.
4. make some pretty pictures and do some science

Here's a couple of videos of what has been done so far

Thin ring evolution: (old version)

https://github.com/user-attachments/assets/df89fe05-3852-4479-b3c9-319c9d578eed

1D Sod Shock tube:

<img width="989" height="790" alt="Shock tube graphs" src="https://github.com/user-attachments/assets/eebe7309-01da-4f9e-b74c-2f60ed991528" />

Uniform Keplerian disc at 100yrs:

<img width="600" height="436" alt="disc at 100yr" src="https://github.com/user-attachments/assets/d4ac6eba-aede-492f-ab74-125e11c5b9a0" />


# How to run
1. go to the SUMMER_SPH.f90 and at to the bottom of the code to the "PROGRAM" section where there is a "filename equals",
   just change this to the name of the .txt initial conditions file you choose.
2. The file should have columns for the [x, y, z, vx, vy, vz, internal energy per unit mass, mass] of all particles including sinks.
3. Any line with the internal energy zeroed will be read as a sink particle.
4. I will add the need of a parameters file, as of now the code uses AU, M_sun, year units (4pi^2 = G) and the end time can be adjusted in the "simulate" subroutine.
5. compile with gfortran via the terminal line "gfortran '.\SUMMER_SPH.f90' -o run_sph"
6. run either via the .exe or again in terminal via ".\run_sph"
