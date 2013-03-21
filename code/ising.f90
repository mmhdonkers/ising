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

  use metropolis
  use plot
 
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! INPUT: Final temperature (Kelvin x 100), final time, row and column size, stepsize loops
  integer,parameter :: tempfinal = 400, timefinal = 100000, size = 30

!! fortran begins indexing from 1. Start it from 0 because the rand() starts from 0
  integer :: spin(0:size-1,0:size-1), temp, time
  real(8) :: weight(-2:2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call plot_init()

!!! Open files for writing data !!!                                                                                                                                      
  call opentextfiles

!! Run main loop subroutine
!! In the following nested loop, the outer loop is the temperature loop 
!!  (for plot of magnetization vs. temperature) 
!! and the inner loop is the converging time loop 
!!  (for plot of magnetization vs. time)

  do temp = 0,tempfinal

! Re-initialize the spin lattice for every temperature
    spin(:,:) = 1       
! Only 5 options for the exponent so calculate them once
    weight = [exp(-800d0/temp),exp(-400d0/temp),1d0,exp(400d0/temp),exp(800d0/temp)]

    do time = 0,timefinal
      call metropolis(spin, size, weight)
! We want time to print only once (choose an arbitrary temperature)
      if (temp == 250) then
        WRITE(16,*) sum(spin)/(size**2*1d0, time
      end if
    end do
      
    call plot_spin(spin, size)
    WRITE(15,*) abs(sum(spin)/(size**2*1d0)), temp/100d0
  end do

!!! Close text files !!!                                                                                                                                                 
  call closetextfiles
  call plot_close()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
