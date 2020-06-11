--------Reference--------
This code implements the algorithm proposed in:
Eyink, G.L., Gupta, A. & Zaki, T. 2019 Stochastic Lagrangian Dynamics of Vorticity. 
I. General Theory. See: http://arxiv.org/abs/1912.06677

-------library requirement------
QUICK SETUP : RUN setup.sh to skip the steps below. Make sure you have git on your system.

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

2. A Fortran95 implementation of the Mersenne Twister algorithm for generating
pseudorandom sequences developed by Scott Robert Ladd. This code was formerly 
available at the official Mersenne Twister website:     
http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/fortran.html
which still provides several other Fortran implementations. With permission of Ladd, 
we provide his source codes here:
mtprng.f90
stdtypes.f90

3. The y locations of the grid points of the channel flow database: 
http://turbulence.pha.jhu.edu/docs/channel/y.txt
We also provide the y.txt file in our home directory.


The main driver code is scauchy.f90

Make sure that gcc and gfortran are installed (gcc/5.5.0 is recommended).

-------Authorization Token for JHTDB------
An authorization token for the Johns Hopkins Turbulent Databases is required to 
run the code. We provide a test token in the code, which will run with a small 
number of particles. If you need more, please follow the instructions on the 
JHTDB website to request your own token:
http://turbulence.pha.jhu.edu/authtoken.aspx

--------compile and run--------
We provide a makefile for compilation. To compile the code, input in your 
terminal window:

make all

To run the compiled program, then type in a terminal window the following commands:

mkdir checkpoint
./scauchy


--------input and output--------
All the parameters are set in scauchy.f90

Input parameters:
Nsamp      ! number of particles, which must be an even integer 
Nstep      ! number of backward integration time steps
dt         ! backward integration time step size (positive) 
timef      ! release time of particles
x0         ! release location (3-vector) of particles  
px         ! number of x-gridpoints to coarse-grain vorticity, deltax=px*dx
py         ! number of y-gridpoints to coarse-grain vorticity, deltay=py*dy
pz         ! number of x-gridpoints to coarse-grain vorticity, deltaz=pz*dz

Note: px = py = pz = 0 means no coarse graining

Output files:
wallparticles.dat   each row is: wall particle index; hitting time; hitting location; Cauchy vorticity-vector 
history.dat         full trajectories of first 30 particles
meanposition.dat    mean position of all particles
varposition.dat     position variance of all particles
omegamn.dat         mean of Cauchy vorticity-vector of all particles
varomega.dat        variance of Cauchy vorticity-vector of all particles
   
If online access to the database is interrupted or if the program 
stalls, files for restarting the program are stored in checkpoint/
You can restart the code by inputting again into the terminal window
the command
 
./scauchy

--------sample results--------
With the parameters given in scauchy.f90, you can try running the code and 
plotting output in omegamn.dat and varomega.dat. We have provided a Matlab
script 

makefigs.m 

to generate such plots. Your results should be similar to the corresponding 
figures in sample_results/

