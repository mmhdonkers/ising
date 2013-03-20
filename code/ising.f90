! Name:      Jared Doster and Chiel Donkers
! Course:    International Course on Computational Physics
! Project:   Monte Carlo Ising Model: Wolff Algorithm and Metropolis Test
! 
!Program Summary (Ising)
!
!  Input:    Hardcode three parameters to the subroutine (temperature of electrons
!            and the length and width of the eletron array).  
!
!  Process:  Creates an array representing the arry of electrons.
!            Calls the subroutine in a nested loop from temp=0 to tempfinal=arbitrary.
!            Each iteration of the outer "temperature" loop returns the magnetizaton for a
!            particular temperature.
!            The iterations of the inner loop will continue until the magnetization returned
!            from the subroutine has converged.
!
!  Output:   Prints magnetization vs. temperature and magnetization vs. time into two data files
!
!
!Subroutine Summary (Mainloop)
!
! Input:     Three parameters from the program
!
! Process:   For a particular temperature, runs a loop involving random selection of electrions
!            and performing a Metropolis Test. Each time that the mainloop is called, the spins
!            of the array are randomly chosen and flipped.
!            Many iterations of this loop cause the calculated magnetization to converge to a range
!            of values.
!
! Output:    Magnetization for a particular temperature for a particular iteration.
!            Output will vary each time that the subroutine is called



program ising
 
  implicit none


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! INPUT: Final temperature (Kelvin), final time, row and column size, step of temperature loop
  integer,parameter :: tempfinal = 40     ! 10 x "Final Temperature" at end of temperature loop
  integer,parameter :: timefinal = 100000 ! Final time at end of time loop
  integer,parameter :: step = 1    
  integer,parameter :: size = 20


!! fortran begins indexing from 1. Start it from 0 because the rand() starts from 0
  integer :: spin(0:size-1,0:size-1)
  real(8) :: mag
  integer :: cputime, tempcount, timecount



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Open files for writing data !!!                                                                                                                                      
  call opentextfiles


!! Seed random number generator

  call system_clock(cputime)
  call srand(cputime)


!! Run main loop subroutine
!! In the following nested loop, the outer loop is the temperature loop (for plot of magnetization vs. temperature) 
!! and the inner loop is the converging time loop (for plot of magnetization vs. time)


  tempcount=0
  do while (tempcount .LE. tempfinal)

! Re-initialize the spin lattice for every temperature
    spin = 1       
!    weight = 10*(/-8,-4,0,4,8/)/tempcount

    timecount=0
    do while (timecount .LE. timefinal)
        call mainloop(spin, size, tempcount, mag, timecount)

! We want time to print only once (choose an arbitrary temperature)
        if (tempcount == 25) then
           WRITE(16,*) mag
        endif
     
        timecount = timecount + step
     enddo

     WRITE(15,*) abs(mag), tempcount

     tempcount = tempcount+1
  enddo

  

!!! Close text files !!!                                                                                                                                                 
  call closetextfiles
contains




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------------------------------
!! Mainloop
subroutine mainloop(spin, size, tempcount, mag, timecount)

  implicit none

!! Passed parameters, intent(in) parameters cannot be altered

  integer,intent(in) :: size, tempcount, timecount !! can remove later!!!!!!!!!!!
  integer,intent(inout) :: spin(:,:)
  real(8) :: mag

 
!! Subroutine declerations
  integer :: i,j,ix,iy
  real(8) :: temp
  real :: totup=0
  real :: totdown=0
  real :: r1,r2,r3
  real :: expo, ediff, enew=0, eold=0
  real :: oldsp, newsp, st, sb, sl, sr, neighbors


!! Tempcount is an integer, we want to rescale it to a real

  temp = tempcount/10.

!! Randomly choose a location in the array (the old spin)

  call random_number(r1)
  call random_number(r2)

  ix = (r1*size)-1
  iy = (r2*size)-1

!  oldsp = spin(ix,iy)

!! Calculate energy due to the neighbors (the if-statements takes into account the free boundaries)

  sl = spin(modulo(ix-1,size), iy)
  sr = spin(modulo(ix+1,size), iy)
  st = spin(ix, modulo(iy-1,size))
  sb = spin(ix, modulo(iy+1,size))

  expo = exp(-2 * spin(ix,iy) * (sl + sr + st + sb)/temp)

!! Calculate energy of old spin

!  if (oldsp .EQ. 1) then
!     eold = 1
!  elseif (oldsp .EQ. -1) then
!     eold = -1
!  else
!     print*, 'At least one oldsp is neither up nor down.  ', oldsp, temp, timecount
!  endif
     

!! Choose new spin & calculate energy of new spin

!  newsp = -oldsp

!  if (newsp == 1) then
!     enew = 1
!  elseif (newsp == -1) then
!     enew = -1
!  endif
    


!! Calculate energy difference between old and new spins

!  ediff = (enew - eold)*neighbors



!! Calculate energy and exponential 

!  expo = exp(ediff/temp)


!! Metropolis test

  call random_number(r3)

  if (expo > r3) then
    spin(ix,iy) = -spin(ix,iy)
  endif


!! Calculate magnetization  (quantifies how magnetic the material is)
!! (first initialize counters)

!  totup = 0
!  totdown = 0
 
!  do i=0,size-1
!    do j=0,size-1

!      if (spin(i,j) == 1) then
!        totup = totup + 1
!      elseif (spin(i,j) == -1) then
!         totdown = totdown + 1
!      else
!         print*, 'At least one spin is neither up nor down'
!      endif

!     end do
!  end do
  mag = sum(spin)/(size*size*1d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Make average of last 20 magnetizations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine mainloop
!----------------------------------------------------------------------------------
!! Open data files

subroutine opentextfiles
    integer :: OPEN_STATUS
    OPEN(UNIT=15,FILE="mag_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, mag_temp file not opened properly------------"
    endif
    OPEN(UNIT=16,FILE="mag_time.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, mag_time file not opened properly------------"
    endif
end subroutine

!---------------------------------------------------------------------------------
!! Close data files

subroutine closetextfiles
    CLOSE(UNIT=15)
    CLOSE(UNIT=16)
end subroutine


end program ising
