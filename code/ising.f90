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

!! INPUT: Final temperature (Kelvin x 10), final time, row and column size, stepsize loops
  integer,parameter :: tempfinal = 40, timefinal = 100000, size = 20
  integer,parameter :: tempstep = 1, timestep = 1


!! fortran begins indexing from 1. Start it from 0 because the rand() starts from 0
  integer :: spin(0:size-1,0:size-1)
  real(8) :: mag
  integer :: temp, time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Open files for writing data !!!                                                                                                                                      
  call opentextfiles

!! Run main loop subroutine
!! In the following nested loop, the outer loop is the temperature loop 
!!  (for plot of magnetization vs. temperature) 
!! and the inner loop is the converging time loop 
!!  (for plot of magnetization vs. time)

  do temp = 0,tempfinal,tempstep

! Re-initialize the spin lattice for every temperature
    spin = 1       
!    weight = 10*(/-8,-4,0,4,8/)/tempcount

    do time = 0,timefinal,timestep
        call metropolis(spin, size, temp/10d0, mag)

! We want time to print only once (choose an arbitrary temperature)
        if (temp == 25) then
           WRITE(16,*) mag, time
        endif
     enddo

     WRITE(15,*) abs(mag), temp/10d0
  enddo

!!! Close text files !!!                                                                                                                                                 
  call closetextfiles

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------------------------------
!! Metropolis mainloop

subroutine metropolis(spin, size, temp, mag)

  implicit none

!! Passed parameters, intent(in) parameters cannot be altered
  integer,intent(in) :: size 
  real(8),intent(in) :: temp
  integer,intent(inout) :: spin(:,:)
  real(8),intent(out) :: mag

!! Subroutine variable declerations
  integer :: ix,iy
  real(8) :: r1,r2,r3
  real(8) :: expo
  real(8) :: st, sb, sl, sr !neighbour cells

!! Randomly choose a location in the array (the old spin)
  call random_number(r1)
  call random_number(r2)

  ix = floor((r1*size))
  iy = floor((r2*size))

!! Calculate energy due to the neighbors (the modulo takes into account the boundaries)
  sl = spin(modulo(ix-1,size), iy)
  sr = spin(modulo(ix+1,size), iy)
  st = spin(ix, modulo(iy-1,size))
  sb = spin(ix, modulo(iy+1,size))

  expo = exp(-2 * spin(ix,iy) * (sl + sr + st + sb)/temp)

!! Metropolis test
  call random_number(r3)

  if (expo > r3) then
    spin(ix,iy) = -spin(ix,iy)
  endif

!! Calculate magnetization  (quantifies how magnetic the material is)
  mag = sum(spin)/(size*size*1d0)

end subroutine

!--------------------------------------------------------------------------------------------
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

!--------------------------------------------------------------------------------------------
!! Close data files

subroutine closetextfiles
    CLOSE(UNIT=15)
    CLOSE(UNIT=16)
end subroutine

end program ising
