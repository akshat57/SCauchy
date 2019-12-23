
program scauchy

!---------------------------------------------------------------------
!
!  This program is a Fortran90 implementation of a Monte Carlo Lagrangian
! algorithm to calculate realizations and statistics of stochastic Cauchy 
! vorticity vectors for spacetime events selected from an online database 
! of a turbulent channel flow. It was written by Gregory L. Eyink, Akshat 
! Gupta, Mengze Wang, and Tamer A. Zaki. This code implements the algorithm 
! proposed in:
!
! Eyink, G.L., Gupta, A. & Zaki, T.A. 2019 Stochastic Lagrangian Dynamics of 
! Vorticity. I. General Theory.   http://arxiv.org/abs/1912.06677
!
! For more information on the program and its use, consult the accompanying
! file "Readme.txt"
!
!---------------------------------------------------------------------

use stdtypes
use mtprng
implicit none


    !!  PARAMETERS TO BE SET FOR STOCHASTIC LAGRANGIAN INTEGRATION

   integer, parameter :: Nsamp = 2000 ! number of stochastic Lagrangian particles (must be an even integer)
   integer, parameter :: Nstep = 2000 ! number of backward integration time steps
   real, parameter :: dt=1d-3   ! step-size for Lagrangian time-integration
   ! time of release of particles:
   real, parameter :: timef = 25.9935
   ! spatial location of release of particles:
   ! You must choose -1<x0(2)<1, in the channel interior
   ! Two examples: (a) an ejection location:
   real, parameter, dimension(3) :: x0=(/21.094707, 0.994647, 7.563944/)
   ! (b) a sweep location:
   ! real, parameter, dimension(3) :: x0=(/0.7150, 0.9951, 0.7259/)
   ! coarse-graining intervals: to eliminate coarse-graining, set px=0,py=0,pz=0
   integer, parameter ::  px=0 ! number of x-grid spacings averaged over
   integer, parameter ::  py=0 ! number of y-grid spacings averaged over
   integer, parameter ::  pz=0 ! number of z-grid spacings averaged over
   ! parameters for repeat calls to the database and restarts from checkpoints
   integer, parameter :: checkit=20  ! number of timesteps between checkpoints for restarting
   integer, parameter :: nattempts=30  ! maximum number of attempts to call the database
   integer, parameter :: tsleep=300  ! clock time in seconds between attempted calls to the database


   !!   SPECIFY YOUR ACCESS KEY HERE  
   
   ! If you need one, please visit http://turbulence.pha.jhu.edu/
   ! (We just want to know a bit about our users!)
   character*100 :: authkey = 'edu.jhu.pha.turbulence.testing-201311' // CHAR(0)


   !! PARAMETERS TO CALL THE JHU TURBULENT CHANNEL-FLOW DATABASE (DO NOT EDIT)

   ! ---- Temporal Interpolation Options ----
   integer, parameter :: NoTInt = 0 ! No temporal interpolation
   integer, parameter :: PCHIPInt = 1 ! Piecewise cubic Hermit interpolation in time
   
   ! ---- Spatial Interpolation Flags for GetVelocity & GetVelocityAndPressure ----
   integer, parameter :: NoSInt = 0 ! No spatial interpolation
   integer, parameter :: Lag4 = 4 ! 4th order Lagrangian interpolation in space
   integer, parameter :: Lag6 = 6 ! 6th order Lagrangian interpolation in space
   integer, parameter :: Lag8 = 8 ! 8th order Lagrangian interpolation in space
   
   ! ---- Spatial Differentiation & Interpolation Flags for GetVelocityGradient & GetPressureGradient ----
   integer, parameter :: FD4NoInt = 40 ! 4th order finite differential scheme for grid values, no spatial interpolation
   integer, parameter :: FD6NoInt = 60 ! 6th order finite differential scheme for grid values, no spatial interpolation
   integer, parameter :: FD8NoInt = 80 ! 8th order finite differential scheme for grid values, no spatial interpolation
   integer, parameter :: FD4Lag4 = 44 ! 4th order finite differential scheme for grid values, 4th order Lagrangian interpolation in space
   
   
   ! ---- Spline interpolation and differentiation Flags for getVelocity,
   !      getPressure, getVelocityGradient, getPressureGradient,
   !      getVelocityHessian, getPressureHessian
   integer, parameter :: M1Q4 = 104   ! Splines with smoothness 1 (3rd order) over 4 data points. Not applicable for Hessian.
   integer, parameter :: M2Q8 = 208   ! Splines with smoothness 2 (5th order) over 8 data points.
   integer, parameter :: M2Q14 = 214  ! Splines with smoothness 2 (5th order) over 14 data points.
   
   !
   ! Choose which dataset to use in this query
   ! Currently, only the channel dataset is supported
   !
   character*100 :: dataset = 'channel' // CHAR(0)
   

   !! DATABASE SIMULATION PARAMETERS (DO NOT EDIT)

   real, parameter :: nu=5d-5  ! viscosity
   real, parameter :: dx = 0.0122656  ! x-grid spacing
   real, parameter :: dz = 0.0061328  ! z-grid spacing


   !! WORKING ARRAY DECLARATIONS

   character(len = 40) :: command
   real time
   real tdeath
   real tcross, xcross(3)
   real, dimension(3,3) :: eye
   real p0(3,1), omega0(3), ugradout9(9,1), ugrad0(3,3), omegainitial(3)
   real points(3, Nsamp), oldpoints(3,Nsamp)    ! spatial positions of particles
   integer status(Nsamp)    ! particle is off wall when status==1, on wall when status==0
   real ugrad(3,3), Xgrad(3,3)   ! temporary matrix arrays
   real history(90) ! saved trajectories of first 30 particles
   real omega(3,Nsamp), omegmn(3) ! stochastic Cauchy vorticity vectors and their mean
   real varomega(3) ! variance of the components of stochastic Cauchy vorticity vector
   real deform(3,3,Nsamp) ! deformation gradient tensors of particles
   real trace    ! trace of velocity gradient
   real, allocatable :: interior(:,:)  ! positions of interior particles
   real, allocatable :: normal(:,:)    ! normal random numbers
   real, allocatable :: dataout3(:,:)  ! results from getVelocity
   real, allocatable :: dataout9(:,:)  ! results from getVelocityGradient
   real meanposition(3), varposition(3)
   real meaninitomega, varinitialomega, temp ! initial mean and variance of vorticity
   real ygrid(512)
   integer :: it0 ! for restart
   logical :: check_exists  ! for restart
   integer :: counter, num_interior, num_interior2, Nwall
   !for stochastic interpolation of particles crossing channel wall
   real dy0, dy1, dtt, lam, mu, eta, xi, tau
   real dyy1, dyy2

   ! declarations for the Mersenne Twister prng algorithm
   real(IEEE64) :: u1, u2, n1, n2
   real(IEEE64) :: pi
   integer(INT32) :: seed
   type(mtprng_state) :: state
   
   real etime
   real timer(2)
   real UserTime, SystemTime, TotTime
       
   character(*), parameter :: fmt3 = "(3(3X,E15.8))"
   character(*), parameter :: fmt4 = "(4(3X,E15.8))"
   character(*), parameter :: fmt8 = "(8(3X,E15.8))"
   character(*), parameter :: fmt90 = "(90(3X,E15.8))"
   
   integer i,j,k,it,jt
   integer gridindex, gridflag
   
   ! For error handling
   ! Declare the return type of the turblib functions as integer.
   ! This is required for error handling
   integer :: getvelocity, getvelocitygradient
   integer, parameter :: SOAP_OK = 0  ! From stdsoap2.h
   integer rc, attempts


   !! BEGIN PROGRAM

   open(unit=20,file='input.dat',status='unknown')
   close(20)
   open(unit=20,file='input.dat',status='old')
   
   ! record initial time for diagnostic purposes
   
   TotTime    = etime(timer)
   UserTime   = timer(1)
   SystemTime = timer(2)

   ! define pi for random-number generator
   pi=datan(1d0)*4d0
   
   ! define identity matrix 
   eye=0
   do i=1,3
      eye(i,i)=1.d0
   end do
   
   !
   ! Intialize the gSOAP runtime.
   ! This is required before any WebService routines are called.
   !
   CALL soapinit()
   CALL turblibSetExitOnError(0)
   
   ! INQUIRE WHETHER CHECKPOINT EXISTS
   inquire(file='checkpoint/seed.dat', exist=check_exists)
   
   !!! INITIALIZE
   if (.not.check_exists) then
      
      !   
      !   set seed for random number generator
      !
      open(89,file='/dev/random',access='stream',form='UNFORMATTED')
      read(89) seed
      close(89)
      write(*,*) 'seed', seed
      
      
      write(20,*) 'seed=', seed
      write(20,*) 'Nsamp=', Nsamp
      write(20,*) 'Nstep=', Nstep
      write(20,*) 'timef=', timef
      write(20,*) 'dt=', dt
      write(20,*) 'x0=', x0(1),'  y0=', x0(2),'  z0=', x0(3)
      write(20,*) 'px=', px, '  py=', py,'  pz=', pz


      open (unit = 222, file = "y.txt")
      read(222,*) ygrid
      close(222)
      
      
      call mtprng_init(seed,state)
      
      
      ! Find grid points between which selected point lies
      gridflag = 0
      gridindex = 1
      do while (gridflag.ne.1)
         if ((x0(2) .gt. ygrid(gridindex)) .and. (x0(2) .le. ygrid(gridindex+1))) then
            gridflag = 1
         else
            gridindex = gridindex + 1
         end if
      end do
      write(*,*) 'gridindex = ', gridindex,ygrid(gridindex),x0(2),ygrid(gridindex+1)
      
      ! Find grid point closest to selected point
      if ((x0(2)-ygrid(gridindex)) .gt. (ygrid(gridindex+1)-x0(2))) then
         gridindex=gridindex+1
      end if
      
      ! DEFINITION OF INTERVAL FOR COARSE-GRAINING IN Y-DIRECTION
      !centering at the closest grid point, take distances to next py/2-th grid points
      if (mod(py,2).eq.0) then
         dyy1 = ygrid(gridindex) - ygrid(gridindex-py/2)
         dyy2 = ygrid(gridindex+py/2) - ygrid(gridindex)
      else
         dyy1 = (ygrid(gridindex) - ygrid(gridindex-(py-1)/2))+(ygrid(gridindex-(py-1)/2)-ygrid(gridindex-(py+1)/2))/2.0
         dyy2 = (ygrid(gridindex+(py-1)/2) - ygrid(gridindex))+(ygrid(gridindex+(py+1)/2)-ygrid(gridindex+(py-1)/2))/2.0
      end if
      
      it0=0
      num_interior=Nsamp
      allocate ( dataout3(3,Nsamp))
      allocate ( dataout9(9,Nsamp))
      
      deform=0
      do i = 1, Nsamp, 1
         status(i)=1
         oldpoints(1,i) = x0(1) - px*dx/2 + mtprng_rand_real3(state)*px*dx;
         oldpoints(2,i) = x0(2) - dyy1 + mtprng_rand_real3(state)*(dyy1 +dyy2);
         oldpoints(3,i) = x0(3) - pz*dz/2 + mtprng_rand_real3(state)*pz*dz;
         do j=1,3
            deform(j,j,i)=1.d0
            if (i .lt. 31) then
               history(j+3*(i-1))=oldpoints(j,i)
            end if
         end do
      end do
      
      
      time=timef
      
      ! CALCULATE INITIAL MEAN VORTICITY AND PARTICLE POSITION MEAN & VARIANCE
      
      ! get velocity at initial particle locations
      attempts=0
      do while (getvelocity(authkey, dataset, time, Lag6, PCHIPInt, num_interior, oldpoints, dataout3).ne.0)
         attempts = attempts + 1
         if (attempts.ge.30) then
            write(*,*) 'Fatal error: too many failures'
            CALL turblibPrintError()
            STOP
         else
            write(*,*) 'Temporary Error (#', attempts, '):'
            CALL turblibPrintError()
            CALL sleep(5)
         end if
      end do
      
      ! get velocity gradient at initial particle locations
      attempts=0
      do while (getvelocitygradient(authkey, dataset, time, FD4Lag4, PCHIPInt, num_interior, oldpoints, dataout9).ne.0)
         attempts = attempts + 1
         if (attempts.ge.30) then
            write(*,*) 'Fatal error: too many failures'
            CALL turblibPrintError()
            STOP
         else
            write(*,*) 'Temporary Error (#', attempts, '):'
            CALL turblibPrintError()
            CALL sleep(5)
         end if
      end do
      
      ! CALCULATE MEAN INITIAL VORTICITY AND POSITION
      do i=1,3
         omegmn(i)=0.0d0
         meanposition(i)=0.0d0
      end do
      
      do i = 1, Nsamp
         do j=1,3
            meanposition(j) = meanposition(j) + oldpoints(j,i)
         end do
         
         ! CALCULATE INITIAL CAUCHY INVARIANT CONTRIBUTIONS
         ugrad(1,1)=dataout9(1,i); ugrad(1,2)=dataout9(2,i);ugrad(1,3)=dataout9(3,i)
         ugrad(2,1)=dataout9(4,i); ugrad(2,2)=dataout9(5,i);ugrad(2,3)=dataout9(6,i)
         ugrad(3,1)=dataout9(7,i); ugrad(3,2)=dataout9(8,i);ugrad(3,3)=dataout9(9,i)
         
         ! Calculate and subtract the trace
         trace=0
         do jt = 1,3
            trace = trace + ugrad(jt,jt)
         end do
         ugrad=ugrad-(trace*eye)/3
         
         ! Calculate vorticity at current time
         ugrad=ugrad-transpose(ugrad)
         omega0(1)=ugrad(3,2); omega0(2)=ugrad(1,3); omega0(3)=ugrad(2,1)
         
         do j=1,3
            omega(j,i)=omega0(j)
         end do
         ! VORTICITY FOR INITIAL PARTICLES IS NOW STORED
         
         do j=1,3
            omegmn(j)=omegmn(j)+omega(j,i)
         end do
      end do
      
      do i=1,3
         omegmn(i)=omegmn(i)/Nsamp
         meanposition(i)=meanposition(i)/Nsamp
      end do
      write(20,*) omegmn(1), omegmn(2), omegmn(3)
      close(20)
      
      meaninitomega = sqrt(omegmn(1)*omegmn(1) + omegmn(2)*omegmn(2) + omegmn(3)*omegmn(3))
      !! DEFINE UNIT VECTOR IN DIRECTION OF MEAN VORTICITY, omegainitial
      do i=1,3
         omegainitial(i)=omegmn(i)/meaninitomega
      end do
      
      ! CALCULATE VARIANCE OF INITIAL VORTICITY AND POSITION
      do i=1,3
         varomega(i)=0.0d0
         varposition(i)=0.0d0
      end do
      varinitialomega=0.0d0
      
      do i = 1, Nsamp
         temp = meaninitomega - (omega(1,i)*omegainitial(1) + omega(2,i)*omegainitial(2) + omega(3,i)*omegainitial(3))
         varinitialomega = varinitialomega + temp*temp
         do j = 1,3
            varomega(j) = varomega(j) + (omega(j,i) -omegmn(j))*(omega(j,i)-omegmn(j))
            varposition(j) = varposition(j) + (oldpoints(j,i) -meanposition(j))*(oldpoints(j,i) - meanposition(j))
         end do
      end do
      varinitialomega=varinitialomega/(Nsamp-1)
      do i=1,3
         varomega(i)=varomega(i)/(Nsamp-1)
         varposition(i)=varposition(i)/(Nsamp-1)
      end do
      
      
      open(10, file='wallparticles.dat', status='new', action='write',position='append')
      open(unit = 11, file='varposition.dat', status='new', action='write',position='append')
      write(11,fmt3) varposition(1), varposition(2), varposition(3)
      open(12, file='varomega.dat', status='new', action='write',position='append')
      write(12,fmt4) varomega(1), varomega(2), varomega(3), varinitialomega
      open(15, file='meanposition.dat', status='new', action='write', position='append')
      write(15,fmt4) time, meanposition(1), meanposition(2), meanposition(3)
      open(30, file='omegamn.dat', status='new', action='write',position='append')
      write(30,fmt4) time, omegmn(1), omegmn(2), omegmn(3)
      open(40, file='history.dat', status='new', action='write',position='append')
      write(40,fmt90) (history(k), k=1,90)
      
      
      command = 'cp wallparticles.dat checkpoint'
      CALL system(command)
      command = 'cp varposition.dat checkpoint'
      CALL system(command)
      command = 'cp varomega.dat checkpoint'
      CALL system(command)
      command = 'cp meanposition.dat checkpoint'
      CALL system(command)
      command = 'cp omegamn.dat checkpoint'
      CALL system(command)
      command = 'cp history.dat checkpoint'
      CALL system(command)
      
      
      ! SAVE INITIAL VALUES
      open(21, file='checkpoint/seed.dat', form='unformatted', status='replace')
      write(21) seed
      close(21)
      open(21, file='checkpoint/omegainitial.dat', form='unformatted', status='replace')
      write(21) omegainitial
      close(21)
      
      ! SAVE CHECKPOINT
      open(22, file='checkpoint/state.dat', form='unformatted', status='replace')
      write(22) state
      close(22)
      open(23, file='checkpoint/iterate.dat', form='unformatted', status='replace')
      write(23) it0
      close(23)
      open(24, file='checkpoint/num_interior.dat', form='unformatted', status='replace')
      write(24) num_interior
      close(24)
      open(25, file='checkpoint/points.dat', form='unformatted', status='replace')
      write(25) oldpoints
      close(25)
      open(26, file='checkpoint/status.dat', form='unformatted', status='replace')
      write(26) status
      close(26)
      open(27, file='checkpoint/deform.dat', form='unformatted', status='replace')
      write(27) deform
      close(27)
      open(28, file='checkpoint/omega.dat', form='unformatted', status='replace')
      write(28) omega
      close(28)
      
      
   !!! RESTARTING
   else
      
      ! READ IN INITIAL VALUES
      open(21, file='checkpoint/seed.dat', form='unformatted', status='old')
      read(21) seed
      close(21)
      write(*,*) 'seed=', seed
      open(21, file='checkpoint/omegainitial.dat', form='unformatted', status='old')
      write(21) omegainitial
      close(21)
      
      ! READ IN CHECKPOINT
      open(22, file='checkpoint/state.dat', form='unformatted', status='old')
      read(22) state
      close(22)
      open(23, file='checkpoint/iterate.dat', form='unformatted', status='old')
      read(23) it0
      close(23)
      time=timef-dfloat(it0)*dt
      write(*,*) 'time=', time
      open(24, file='checkpoint/num_interior.dat', form='unformatted', status='old')
      read(24) num_interior
      close(24)
      write(*,*) 'num_interior=', num_interior
      open(25, file='checkpoint/points.dat', form='unformatted', status='old')
      read(25) oldpoints
      close(25)
      open(26, file='checkpoint/status.dat', form='unformatted', status='old')
      read(26) status
      close(26)
      open(27, file='checkpoint/deform.dat', form='unformatted', status='old')
      read(27) deform
      close(27)
      open(28, file='checkpoint/omega.dat', form='unformatted', status='old')
      read(28) omega
      close(28)
      
      ! OPEN DATA FILES
      command = 'cp checkpoint/wallparticles.dat .'
      CALL system(command)
      command = 'cp checkpoint/varposition.dat .'
      CALL system(command)
      command = 'cp checkpoint/varomega.dat .'
      CALL system(command)
      command = 'cp checkpoint/meanposition.dat .'
      CALL system(command)
      command = 'cp checkpoint/omegamn.dat .'
      CALL system(command)
      command = 'cp checkpoint/history.dat .'
      CALL system(command)
      
      open(10, file='wallparticles.dat', status='old', action='write', position='append')
      open(11, file='varposition.dat', status='old', action='write', position='append')
      open(12, file='varomega.dat', status='old', action='write', position='append')
      open(15, file='meanposition.dat', status='old', action='write', position='append')
      open(30, file='omegamn.dat', status='old', action='write', position='append')
      open(40, file='history.dat', status='old', action='write', position='append')
      
      
      counter=1
      allocate ( interior(3,num_interior))
      do i = 1, Nsamp
         ! fill array interior
         if (status(i) .eq. 1) then
            do j=1,3
               interior(j,counter) = oldpoints(j,i)
            end do
            counter=counter+1
         end if
      end do
      
      allocate ( dataout3(3,num_interior))
      allocate ( dataout9(9,num_interior))
      
      ! get velocity at interior particle locations
      attempts=0
      do while (getvelocity(authkey, dataset, time, Lag6, PCHIPInt, num_interior, interior, dataout3).ne.0)
        attempts = attempts + 1
        if (attempts.ge.nattempts) then
          write(*,*) 'Fatal error: too many failures'
          CALL turblibPrintError()
          STOP
        else
          write(*,*) 'Temporary Error (#', attempts, '):'
          CALL turblibPrintError()
          CALL sleep(tsleep)
        end if
      end do
      
      ! get velocity gradient at interior particle locations
      attempts=0
      do while (getvelocitygradient(authkey, dataset, time, FD4Lag4, PCHIPInt, num_interior, interior, dataout9).ne.0)
      attempts = attempts + 1
        if (attempts.ge.nattempts) then
          write(*,*) 'Fatal error: too many failures'
          CALL turblibPrintError()
          STOP
        else
          write(*,*) 'Temporary Error (#', attempts, '):'
          CALL turblibPrintError()
          CALL sleep(tsleep)
        end if
      end do
      
      deallocate (interior)
      
      
      
   ! END INITIALIZATION
   end if
   
   
   
   !!! MAIN INTEGRATION LOOP BACKWARD IN TIME
   do it=it0+1,Nstep
      ! write current time-step to screen; comment out this command to eliminate screen ouput
      write(*,*) it
      
      if (mod(num_interior,2).eq.0) then
         num_interior2=num_interior
      else
         num_interior2=num_interior+1
      end if
     
      allocate ( normal(3,num_interior2))
      
      ! Generate an array of 3*num_interior2 normal random numbers
      do i = 1, num_interior2-1, 2
         do j=1,3
            u1=mtprng_rand_real3(state)
            u2=mtprng_rand_real3(state)
            ! Box-Muller transform
            normal(j,i)=sngl(dsqrt(-2.0d0*dlog(u1))*dcos(2.0d0*pi*u2))
            normal(j,i+1)=sngl(dsqrt(-2.0d0*dlog(u1))*dsin(2.0d0*pi*u2))
         end do	
      end do 
      
      
      counter=1
      Nwall=0
      !! INTEGRATE INTERIOR PARTICLES ONE STEP BACKWARD IN TIME AND FIND WHICH NEWLY HIT THE WALL
      do i = 1, Nsamp
         
         ! Euler-Maruyama integration backward in time for interior points
         if (status(i) .eq. 1) then
            
            ! eqn(3.1)
            do j=1,3
               points(j,i) = oldpoints(j,i) - dt*dataout3(j,counter) + sqrt(2.0*nu*dt)*normal(j,counter)
            end do
            
            ! UPDATE DEFORMATION MATRIX FOR PARTICLES THAT REMAIN INTERIOR
            if((points(2,i) .gt. -1) .and. (points(2,i) .lt. 1)) then
               
               ugrad(1,1)=dataout9(1,counter); ugrad(1,2)=dataout9(2,counter); ugrad(1,3)=dataout9(3,counter)
               ugrad(2,1)=dataout9(4,counter); ugrad(2,2)=dataout9(5,counter); ugrad(2,3)=dataout9(6,counter)
               ugrad(3,1)=dataout9(7,counter); ugrad(3,2)=dataout9(8,counter); ugrad(3,3)=dataout9(9,counter)
               
               !Calculating and subtracting the trace
               trace=0d0
               do jt = 1,3
                  trace = trace + ugrad(jt,jt)
               end do
               ugrad=ugrad-(trace*eye)/3
               
               Xgrad(1,1)=deform(1,1,i); Xgrad(1,2)=deform(1,2,i); Xgrad(1,3)=deform(1,3,i)
               Xgrad(2,1)=deform(2,1,i); Xgrad(2,2)=deform(2,2,i); Xgrad(2,3)=deform(2,3,i)
               Xgrad(3,1)=deform(3,1,i); Xgrad(3,2)=deform(3,2,i); Xgrad(3,3)=deform(3,3,i)
               
               ! Euler integration backward in time for deformation matrix
               ! See eqn(3.4) in the paper
               Xgrad=matmul(Xgrad,eye+dt*ugrad)
               
               deform(1,1,i)=Xgrad(1,1); deform(1,2,i)=Xgrad(1,2); deform(1,3,i)=Xgrad(1,3)
               deform(2,1,i)=Xgrad(2,1); deform(2,2,i)=Xgrad(2,2); deform(2,3,i)=Xgrad(2,3)
               deform(3,1,i)=Xgrad(3,1); deform(3,2,i)=Xgrad(3,2); deform(3,3,i)=Xgrad(3,3)
               
            ! CALCULATE FINAL CONTRIBUTIONS FOR NEW WALL PARTICLES
            else
               
               ! stochastic estimation of hitting time
               ! See appendix A of paper
               if((points(2,i) .ne. -1) .and. (points(2,i) .ne. 1)) then
                  ! change of coordinates
                  ! dy0, dy1 are b0, b1 in appendix A
                  if( points(2,i) .le. -1) then
                     dy0=1+oldpoints(2,i)
                     dy1=1+points(2,i)
                     points(2,i) = -1
                     xcross(2) = -1
                  else
                     dy0=1-oldpoints(2,i)
                     dy1=1-points(2,i)
                     points(2,i) = 1
                     xcross(2) = 1
                  end if
                  
                  ! See eqn(A.10) in paper
                  lam=dy0**2/(2*nu*dt)
                  mu=-dy0/dy1
                  
                  ! See eqn(A.11, A.12) in paper
                  u1=mtprng_rand_real3(state)
                  u2=mtprng_rand_real3(state)
                  n1=sqrt(-2.0d0*dlog(u1))*dcos(2.0d0*pi*u2)
                  u1=mtprng_rand_real3(state)
                  eta=sngl(n1**2)
                  xi=mu*eta-sqrt(abs(4*mu*lam*eta+(mu**2)*(eta**2)))
                  xi=mu*(1+xi/(2*lam))
                  if (u1 .le. mu/(mu+xi)) then
                      tau=xi
                  else
                      tau=mu**2/xi
                  end if
                  dtt=tau*dt/(1d0+tau)
                  tcross=time-dtt
                  
                  ! stochastic estimation of xz locations
                  u1=mtprng_rand_real3(state)
                  u2=mtprng_rand_real3(state)
                  n1=dsqrt(-2.0d0*dlog(u1))*dcos(2.0d0*pi*u2)
                  n2=dsqrt(-2.0d0*dlog(u1))*dsin(2.0d0*pi*u2)
                  ! See eqn(A.18) in paper
                  xcross(1) = (points(1,i)*dtt + oldpoints(1,i)*(dt-dtt))/dt + sqrt(abs(2*sngl(nu)*tau*dt))*sngl(n1)/(1+tau)
                  ! See eqn(A.19) in paper
                  xcross(3) = (points(3,i)*dtt + oldpoints(3,i)*(dt-dtt))/dt + sqrt(abs(2*sngl(nu)*tau*dt))*sngl(n2)/(1+tau)
                  points(1,i) = xcross(1)
                  points(3,i) = xcross(3)
               end if
               
               
               tdeath=tcross
               status(i)=0
               Nwall=Nwall+1
               
               
               
               ! UPDATE DEFORMATION MATRIX FOR NEW WALL PARTICLES
               ugrad(1,1)=dataout9(1,counter); ugrad(1,2)=dataout9(2,counter); ugrad(1,3)=dataout9(3,counter)
               ugrad(2,1)=dataout9(4,counter); ugrad(2,2)=dataout9(5,counter); ugrad(2,3)=dataout9(6,counter)
               ugrad(3,1)=dataout9(7,counter); ugrad(3,2)=dataout9(8,counter); ugrad(3,3)=dataout9(9,counter)
               
               ! Calculate and subtract the trace
               trace=0d0
               do jt = 1,3
                  trace = trace + ugrad(jt,jt)
               end do
               ugrad=ugrad-(trace*eye)/3
               
               Xgrad(1,1)=deform(1,1,i); Xgrad(1,2)=deform(1,2,i); Xgrad(1,3)=deform(1,3,i)
               Xgrad(2,1)=deform(2,1,i); Xgrad(2,2)=deform(2,2,i); Xgrad(2,3)=deform(2,3,i)
               Xgrad(3,1)=deform(3,1,i); Xgrad(3,2)=deform(3,2,i); Xgrad(3,3)=deform(3,3,i)
               
               ! Euler integration backward in time for deformation matrix
               Xgrad=matmul(Xgrad,eye+dtt*ugrad)
               
               deform(1,1,i)=Xgrad(1,1); deform(1,2,i)=Xgrad(1,2); deform(1,3,i)=Xgrad(1,3)
               deform(2,1,i)=Xgrad(2,1); deform(2,2,i)=Xgrad(2,2); deform(2,3,i)=Xgrad(2,3)
               deform(3,1,i)=Xgrad(3,1); deform(3,2,i)=Xgrad(3,2); deform(3,3,i)=Xgrad(3,3)
               
               ! CALL DATABASE FOR VORTICITY AND CALCULATE CAUCHY INVARIANT CONTRIBUTION FOR NEW WALL PARTICLES
               
               ! get velocity gradient at new wall particle location
               p0(1,1)=xcross(1); p0(2,1)=xcross(2); p0(3,1)=xcross(3);
               attempts=0
               do while (getvelocitygradient(authkey, dataset, tcross, FD4Lag4, PCHIPInt,1,p0, ugradout9).ne.0)
                  attempts = attempts + 1
                  if (attempts.ge.nattempts) then
                     write(*,*) 'Fatal error: too many failures'
                     CALL turblibPrintError()
                     STOP
                  else
                     write(*,*) 'Temporary Error (#', attempts, '):'
                     CALL turblibPrintError()
                     CALL sleep(tsleep)
                  end if
               end do
               
               ugrad(1,1)=ugradout9(1,1); ugrad(1,2)=ugradout9(2,1); ugrad(1,3)=ugradout9(3,1)
               ugrad(2,1)=ugradout9(4,1); ugrad(2,2)=ugradout9(5,1); ugrad(2,3)=ugradout9(6,1)
               ugrad(3,1)=ugradout9(7,1); ugrad(3,2)=ugradout9(8,1); ugrad(3,3)=ugradout9(9,1)
               
               ! Calculate and subtract the trace
               trace=0d0
               do jt = 1,3
               trace = trace + ugrad(jt,jt)
               end do
               ugrad=ugrad-(trace*eye)/3
               
               ! calculate Cauchy vorticity vector for new wall particle
               ugrad=ugrad-transpose(ugrad)
               omega0(1)=ugrad(3,2); omega0(2)=ugrad(1,3); omega0(3)=ugrad(2,1)
               omega0=matmul(Xgrad,omega0)
               
               do j=1,3
                  omega(j,i)=omega0(j)
               end do
               write(10,fmt8) real(i), tdeath, xcross(1), xcross(2), xcross(3), omega0(1), omega0(2), omega0(3)
            end if
            
            counter = counter+1
            
         ! keep points the same for wall particles
         else
            
            do j = 1,3
               points(j,i) = oldpoints(j,i)
            end do
         
         end if
         
         ! CREATE VECTOR FOR HISTORY OF FIRST 30 PARTICLES
         do j=1,3
         if (i .lt. 31) then
            history(j+3*(i-1))=points(j,i)
         end if
         end do
         
         
      !! END ONE BACKWARD INTEGRATION STEP
      end do

      time=time-dt
      deallocate (normal)
      deallocate(dataout3)
      deallocate(dataout9)
      num_interior=num_interior-Nwall
      
      
      !! GET VELOCITIES AND VELOCITY-GRADIENTS FOR ALL INTERIOR PARTICLES
      counter=1
      allocate ( interior(3,num_interior))
      do i = 1, Nsamp
         ! fill array interior
         if (status(i) .eq. 1) then
            do j=1,3
              interior(j,counter) = points(j,i)
            end do
            counter=counter+1
         end if
      end do
      
      allocate ( dataout3(3,num_interior))
      allocate ( dataout9(9,num_interior))
      
      ! get velocity at interior particle locations
      attempts=0
      do while (getvelocity(authkey, dataset, time, Lag6, PCHIPInt, num_interior, interior, dataout3).ne.0)
         attempts = attempts + 1
         if (attempts.ge.nattempts) then
            write(*,*) 'Fatal error: too many failures'
            CALL turblibPrintError()
            STOP
         else
            write(*,*) 'Temporary Error (#', attempts, '):'
            CALL turblibPrintError()
            CALL sleep(tsleep)
         end if
      end do
      
      ! get velocity gradient at interior particle locations
      attempts=0
      do while (getvelocitygradient(authkey, dataset, time, FD4Lag4, PCHIPInt, num_interior, interior, dataout9).ne.0)
         attempts = attempts + 1
         if (attempts.ge.nattempts) then
            write(*,*) 'Fatal error: too many failures'
            CALL turblibPrintError()
            STOP
         else
            write(*,*) 'Temporary Error (#', attempts, '):'
            CALL turblibPrintError()
            CALL sleep(tsleep)
         end if
      end do
      
      deallocate (interior)

      
      ! CALCULATE ENSEMBLE AVERAGE CAUCHY VORTICITY STATISTICS AT NEW TIME
      do i=1,3
         omegmn(i)=0.0d0
         varomega(i)=0.0d0
         meanposition(i)=0.0d0
         varposition(i)=0.0d0
      end do
      meaninitomega=0.0d0
      varinitialomega=0.0d0
      counter=1
      do i = 1, Nsamp
         do j=1,3
            meanposition(j) = meanposition(j) + points(j,i)
         end do
         
         ! CALCULATE NEW CAUCHY INVARIANT CONTRIBUTION FOR INTERIOR PARTICLES
         if (status(i) .eq. 1) then
            ugrad(1,1)=dataout9(1,counter); ugrad(1,2)=dataout9(2,counter); ugrad(1,3)=dataout9(3,counter)
            ugrad(2,1)=dataout9(4,counter); ugrad(2,2)=dataout9(5,counter); ugrad(2,3)=dataout9(6,counter)
            ugrad(3,1)=dataout9(7,counter); ugrad(3,2)=dataout9(8,counter); ugrad(3,3)=dataout9(9,counter)
            
            ! Calculate and subtract the trace
            trace=0
            do jt = 1,3
             trace = trace + ugrad(jt,jt)
            end do
            ugrad=ugrad-(trace*eye)/3
            
            ! Calculate vorticity at current time
            ugrad=ugrad-transpose(ugrad)
            omega0(1)=ugrad(3,2); omega0(2)=ugrad(1,3); omega0(3)=ugrad(2,1)
            
            ! Calculate Cauchy vorticity vector at current time
            Xgrad(1,1)=deform(1,1,i); Xgrad(1,2)=deform(1,2,i); Xgrad(1,3)=deform(1,3,i)
            Xgrad(2,1)=deform(2,1,i); Xgrad(2,2)=deform(2,2,i); Xgrad(2,3)=deform(2,3,i)
            Xgrad(3,1)=deform(3,1,i); Xgrad(3,2)=deform(3,2,i); Xgrad(3,3)=deform(3,3,i)
            omega0=matmul(Xgrad,omega0)
            
            do j=1,3
               omega(j,i)=omega0(j)
            end do
            ! NEW CAUCHY VORTICITY VECTOR FOR INTERIOR PARTICLE NOW STORED
            counter=counter+1
         end if
         
         do j=1,3
            omegmn(j)=omegmn(j)+omega(j,i)
         end do
         meaninitomega = meaninitomega + omega(1,i)*omegainitial(1) + omega(2,i)*omegainitial(2) + omega(3,i)*omegainitial(3)
         
      end do
      
      meaninitomega=meaninitomega/Nsamp
      do i=1,3
         omegmn(i)=omegmn(i)/Nsamp
         meanposition(i)=meanposition(i)/Nsamp
      end do
      do i = 1, Nsamp
         temp = meaninitomega - (omega(1,i)*omegainitial(1) + omega(2,i)*omegainitial(2) + omega(3,i)*omegainitial(3))
         varinitialomega = varinitialomega + temp*temp
         do j = 1,3
            oldpoints(j,i) = points(j,i)
            varomega(j) = varomega(j) + (omega(j,i) -omegmn(j))*(omega(j,i) -omegmn(j))
            varposition(j) = varposition(j) + (points(j,i) - meanposition(j))*(points(j,i) - meanposition(j))
         end do
      end do
      varinitialomega=varinitialomega/(Nsamp-1)
      do i=1,3
         varomega(i)=varomega(i)/(Nsamp-1)
         varposition(i)=varposition(i)/(Nsamp-1)
      end do
      
      write(11,fmt3) varposition(1), varposition(2), varposition(3)
      write(12,fmt4) varomega(1), varomega(2), varomega(3), varinitialomega
      write(15,fmt4) time, meanposition(1), meanposition(2), meanposition(3)
      write(30,fmt4) time, omegmn(1), omegmn(2), omegmn(3)
      write(40,fmt90) (history(k), k=1,90)
      
      
      ! REPLACE CHECKPOINT FILES
      if (mod(it,checkit).eq.0) then
         
         open(22, file='checkpoint/state.dat', form='unformatted', status='replace')
         write(22) state
         close(22)
         open(23, file='checkpoint/iterate.dat', form='unformatted', status='replace')
         write(23) it
         close(23)
         open(24, file='checkpoint/num_interior.dat', form='unformatted', status='replace')
         write(24) num_interior
         close(24)
         open(25, file='checkpoint/points.dat', form='unformatted', status='replace')
         write(25) points
         close(25)
         open(26, file='checkpoint/status.dat', form='unformatted', status='replace')
         write(26) status
         close(26)
         open(27, file='checkpoint/deform.dat', form='unformatted', status='replace')
         write(27) deform
         close(27)
         open(28, file='checkpoint/omega.dat', form='unformatted', status='replace')
         write(28) omega
         close(28)
         
         command = 'cp wallparticles.dat checkpoint'
         CALL system(command)
         command = 'cp varposition.dat checkpoint'
         CALL system(command)
         command = 'cp varomega.dat checkpoint'
         CALL system(command)
         command = 'cp meanposition.dat checkpoint'
         CALL system(command)
         command = 'cp omegamn.dat checkpoint'
         CALL system(command)
         command = 'cp history.dat checkpoint'
         CALL system(command)
         
      end if
      
      
   ! END MAIN INTEGRATION LOOP
   end do
   
   write(71,*) 'completed'
   close(10)
   close(11)
   close(12)
   close(15)
   close(30)
   close(40)
   
   deallocate ( dataout3)
   deallocate ( dataout9)
   
   TotTime    = etime(timer) 
   UserTime   = timer(1)     
   SystemTime = timer(2)   
       
   ! print *,'Total time:', TotTime, 'secs'
   ! print *,'System time:', SystemTime, 'secs'
   ! print *,'User time:', UserTime, 'secs'


   ! Destroy the gSOAP runtime.
   ! No more WebService routines may be called.
   !
   CALL soapdestroy()
   
end program scauchy
