--------Reference--------
This code implements the algorithm proposed in:
Eyink, G.L., Gupta, A. & Zaki, T. 2019 Stochastic Lagrangian Dynamics of Vorticity. I. General Theory.

-------libary requirement------
This code needs the following libraries / functions:
1. JHU Turbulence DataBase Cluster C and Fortran Interface Library
https://github.com/idies/turblib
You need to download at least the following files:
TurbulenceService.h
TurbulenceServiceSoap.nsmap
soapC.c
soapClient.c
soapH.h
soapStub.h
stdsoap2.c
stdsoap2.h
turblib.c
turblib.h

2. pseudorandom number generator Marsenne Twister
http://coyotegulch.scottrobertladd.net/products/brahe/index.html
With the permission by the author, we provide the source codes here:
mtprng.f90
stdtypes.f90

3. The y locations of the grid points of the channel flow database
http://turbulence.pha.jhu.edu/docs/channel/y.txt
Please rename it as ygrid.txt

The main source code is strack.f90

Make sure that gcc and gfortran are installed (gcc/5.5.0 is recommended).

-------Authorization Token for JHTDB------
An authorization token for the Johns Hopkins Turbulent Databases is required to run the code.
We provide a test token in the code. If you need request a larger size of velocity data from the turbulent database, please follow the instructions on JHTDB website to ask for a token:
http://turbulence.pha.jhu.edu/authtoken.aspx

--------compile and run--------
We provide a makefile for compilation 

To compile the code, input:

make all

in your teminal. 

Then input:

mkdir checkpoint
./strack

The program will start running.

--------input and output--------
All the parameters are set in strack.f90

Input parameters:
Nsamp      ! number of particles, must be even
Nstep      ! number of backward integration time steps
timef      ! forward time of final vorticity
dt         ! backward integration time step size
x0         ! location of final vorticity 
px         ! the size of of the block over which we average, deltax=px*dx
py         ! the size of of the block over which we average, deltay=py*dy
pz         ! the size of of the block over which we average, deltaz=pz*dz

Note: px = py = pz = 0 means no coarse graining

Output files:
wallparticles.dat   each row is: wall particle index; hitting time; hitting locations; vorticity
history.dat         trajectories of first 30 particles
meanposition.dat    mean position of all particles
varposition.dat     position variance of all particles
omegamn.dat         mean vorticity of all particles
varomega.dat        vorticity variance of all particles
   
Files for restarting the program are stored in checkpoint/

--------sample results--------
With the parameters given in strack.f90, you can try running the code and visualizing omegamn.dat and varomega.dat. 
Your results should be similar to the figures in sample_results/

