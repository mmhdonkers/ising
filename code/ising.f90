! Name:      Henry Jared Doster
! Course:    PHY 832
! Project:   Monte Carlo and Metropolis Test applied to Ising Model
! 
!Program Summary
!
!  Input:    Hardcode three parameters to the subroutine (temperature of electrons
!            and the length and width of the eletron array).  
!
!  Process:  Creates an array representing the arry of electrons.
!            Calls the subroutine in a nested loop from T=0 to T=arbitrary.
!            Each iteration of the outer "temperature" loop returns the magnetizaton for a
!            particular temperature.
!            The iterations of the inner loop will continue until the magnetization returned
!            from the subroutine has converged to within a tolerance range
!
!  Output:   Stores the converged magnetization as a function of temperature for each iteration
!            of the outer "temperature" loop
!
!
!Subroutine Summary
!
! Input:     Three parameters from the program
!
! Process:   For a particular temperature, runs a loop involving random selection of electrions
!            and performing a Metropolis Test
!            Many iterations of this loop cause the calculated magnetization to converge to a range
!            of values
!
! Output:    Magnetization for a particular temperature for a particular iteration.
!            Output will vary each time that the subroutine is called




program firstproj
 
  use plplot

  implicit none

!! INPUT: Final temperature (Kelvin), row and column size, step of temperature loop
  real,parameter :: T = 3d0
  integer,parameter :: rowsize = 10
  integer,parameter :: colsize = 10
  real, parameter :: step =.01

!! fortran begins indexing from 1. Start it from 0 because the rand() starts from 0
  real, dimension(0:rowsize-1, 0:colsize-1) :: spin


  real :: count = 0
  real :: mag
  integer :: time, i, j

!!!!!!!!!!! STORE MAGNETIZATION AS A FUNCTION OF ITERATIONS HARCODED SIZE OF MATRIX!!!!!!!!!!!!!!!!!!!!
  real,dimension(5) :: mag_time

!!!!!!!!!!! STORE MAGNETIZATION AS A FUNCTION OF TEMPERATURE
!  real, dimension(0:T/step) :: mag_temp

  real(8) :: rand

!! Initialize Configuration (spin up = 1, spin down = -1)

  do i=0,rowsize-1
    do j=0,colsize-1
        call random_number(rand)
        spin(i,j) = 2 * nint(rand) - 1
     end do
  end do


!! Seed random number generator

  call system_clock(time)
  call srand(time)


!! Call subroutine and run main loop

! Initialize counters  !!!!!!!!!!! THIS MAY BE DELETED
  i = 0    
  count = 0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!! MAKE A TEMPERATURE DO LOOP: RUNS THROUGH TEMPERATURES (using count = count + step)


!!!!!!!!!!++++++++++++++++++++ CONVERGING DO LOOP

  do while (count .LE. 10)
    call mainloop(spin, rowsize, colsize, T, mag)

!    magtol(i)=mag  !!!!!!!!!!!!! MAJOR SEG FAULT HERE!!!!!!!!
    print*, mag
    count = count + step
    i = i + 1
    read*,    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! END THIS DO LOOP ONCE MAGNETIZATION  HAS CONVERGED
  enddo

!!!!!!!!!!!!++++++++++++++++ END CONVERGINE DO LOOP

!! Store converged magnetization for a particular temperature
!!! END TEMPERATURE DO LOOP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end program firstproj




subroutine mainloop(spin, rowsize, colsize, T, mag)

  implicit none

!! Passed parameters, intent(in) parameters cannot be altered

  integer,intent(in) :: rowsize
  integer,intent(in) :: colsize
  real,intent(in) :: T
  real,dimension(0:rowsize-1, 0:colsize-1) :: spin
  real :: mag

 
!! Subroutine declerations
  integer :: i,j,ix,iy
  real :: totup
  real :: totdown
  real :: r1,r2,r3
  real :: expo, ediff, enew, eold
  real :: oldsp, newsp, st, sb, sl, sr, neighbors


!! Randomly choose a location in the array (the old spin)

  r1 = rand()
  r2 = rand()

  ix = (r1*rowsize)-1
  iy = (r2*colsize)-1

  oldsp = spin(ix,iy)


!! Calculate energy due to the neighbors (if statement takes into account free boundaries)

  sl = spin(modulo(ix-1,rowsize), iy)
  sr = spin(modulo(ix+1,rowsize), iy)
  st = spin(ix, modulo(iy-1,colsize))
  sb = spin(ix, modulo(iy+1,colsize))

!  if (ix .EQ. 0) sl = 0
!  if (ix .EQ. rowsize) sr = 0
!  if (iy .EQ. 0) st = 0
!  if (iy .EQ. colsize) sb = 0

  neighbors = sl + sr + st + sb


!! Calculate energy of old spin

!  if (oldsp .EQ. 1) then
!     eold = 1
!  elseif (oldsp .EQ. -1) then
!     eold = -1
!  else
!     print*, 'At least one oldsp is neither up nor down.'
!  endif
     

!! Choose new spin & calculate energy of new spin

!  newsp = -oldsp

!  if (newsp == 1) then
!     enew = 1
!  elseif (newsp == -1) then
!     enew = -1
!  endif
    


!! Calculate energy difference between old and new spins

  ediff = (-oldsp - oldsp)*neighbors



!! Calculate energy and exponential 

  expo = exp(ediff/T)


!! Metropolis test

  r3 = rand()

  if (expo > r3) then
    spin(ix,iy) = newsp
  endif


!! Calculate averages (first initialize counters)

  totup = sum(spin)
!  totup = 0
 
!  do i=0,rowsize-1
!    do j=0,colsize-1

!      if (spin(i,j) == 1) then
!        totup = totup + 1
!      elseif (spin(i,j) == -1) then
!         totdown = totdown + 1
!      else
!         print*, 'At least one spin is neither up nor down'
!      endif

!     end do
!  end do

  mag = (totup)/(rowsize*colsize)

end subroutine mainloop
