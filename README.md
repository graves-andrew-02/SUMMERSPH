# SUMMERSPH - Overview

This is a Summer Project because I am very, very bored, and it might be cool.

# Goals:
1. make an sph code that doesnt completely suck and runs on a laptop

Thats about it really, this is made with me in mind, so might be complete crap.



Here's a couple of videos of what has been done so far

Thin ring evolution:

https://github.com/user-attachments/assets/df89fe05-3852-4479-b3c9-319c9d578eed

NOTFINISHED!!!
1D Sod Shock tube:

<img width="578" height="455" alt="sod1" src="https://github.com/user-attachments/assets/533c0a59-b49a-448a-b363-45abea3d7a80" />
<img width="578" height="455" alt="sod2" src="https://github.com/user-attachments/assets/4a678e4f-fc02-40a0-92f0-bb1c26d81a7c" />
<img width="578" height="455" alt="sod3" src="https://github.com/user-attachments/assets/72d8fdd7-0d31-4f36-bd0f-75f2c4a6704e" />





# How to run
1. go to the SUMMER_SPH.f90 and at to the bottom of the code to the "PROGRAM" section where there is a "filename equals",
   just change this to the name of the .txt initial conditions file you choose.
2. The file should have columns for the [x, y, z, vx, vy, vz, internal energy per unit mass, mass] of all particles including sinks.
3. Any line with the internal energy zeroed will be read as a sink particle.
4. I will add the need of a parameters file, as of now the code uses AU, M_sun, year units (4pi^2 = G) and the end time can be adjusted in the 
