# SUMMERSPH - Overview

This is a Summer Project because I am very, very bored, and it might be cool.

# Goals:
1. make an sph code that doesnt completely suck and runs on a laptop

Thats about it really, this is made with me in mind, so might be complete crap.



Here's a couple of videos of what has been done so far

Thin ring evolution:

https://github.com/user-attachments/assets/df89fe05-3852-4479-b3c9-319c9d578eed




# How to run
1. go to the SUMMER_SPH.f90 and at to the bottom of the code to the "PROGRAM" section where there is a "filename equals",
   just change this to the name of the .txt initial conditions file you choose.
2. The file should have columns for the [x, y, z, vx, vy, vz, internal energy per unit mass, mass] of all particles including sinks.
3. Any line with the internal energy zeroed will be read as a sink particle.
4. I will add the need of a parameters file, as of now the code uses AU, M_sun, year units (4pi^2 = G) and the end time can be adjusted in the 
