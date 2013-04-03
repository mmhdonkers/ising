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

  use model
  use plot
 
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Declarations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! INPUT: Row and column size
  integer,parameter :: SIZE = 32

!! INPUT: Final temperature (Kelvin x 100), final time, stepsizeloops
  integer,parameter :: TEMPFINAL = 400, TIMEFINAL = 10000

!! fortran begins indexing from 1. Start it from 0 because the rand() starts from 0
  integer :: spin(0:SIZE-1,0:SIZE-1), temp, time
  real(8) :: weight(-2:2)
  integer :: wolffspin(0:SIZE-1, 0:SIZE-1)
  integer :: swenwangspin(0:SIZE-1, 0:SIZE-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Main Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call plot_init()

!!! Open files for writing data !!!                                                                                                                                      
  call opentextfiles

!! Run main loop subroutines (both metropolis and Wolff run simultaneously)
!! In the following nested loops, the outer loop is the temperature loop 
!!  (for plot of magnetization vs. temperature) 
!! and the inner loop is the converging time loop 
!!  (for plot of magnetization vs. time)


  do temp = 100,TEMPFINAL

! Re-initialize Wolff lattice. All values in the lattice must be 1
    wolffspin(:,:) = 1
    do time = 0, 100
      call wolff(wolffspin, SIZE, temp/100d0)
    end do
    
! Re-initialize Swendswon-Wang lattice. All values in the lattice must be 1
    swenwangspin(:,:) = 1
    do time = 0, 400
      call swenwang(swenwangspin, SIZE, temp/100d0)
    end do

! Re-initialize the Metropolis lattice. All values in the lattice must be 1
    spin(:,:) = 1       
! Only 5 options for the exponent for metropolis so calculate them once
    weight = [exp(-800d0/temp),exp(-400d0/temp),1d0,exp(400d0/temp),exp(800d0/temp)]
    do time = 0,TIMEFINAL
      call metropolis(spin, SIZE, weight)

! We want time to print only once (choose an arbitrary temperature)
      if (temp == 250) then
        WRITE(16,*) sum(spin)/(SIZE**2*1d0), time
      end if
    end do

      
    call plot_spin(spin, SIZE, temp/100d0)
    WRITE(15,*) abs(sum(spin)/(SIZE**2*1d0)), temp/100d0
    WRITE(17,*) abs(sum(swenwangspin)/(SIZE**2*1d0)), temp/100d0
    WRITE(18,*) abs(sum(wolffspin)/(SIZE**2*1d0)), temp/100d0
  end do
!!! Close text files !!!                                                                                                                                                 
  call closetextfiles
  call plot_close()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine opentextfiles
    integer :: OPEN_STATUS
    OPEN(UNIT=15,FILE="metrop_mag_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, metrop_mag_temp file not opened properly------------"
    endif
    OPEN(UNIT=16,FILE="metrop_mag_time.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, metrop_mag_time file not opened properly------------"
    endif

    OPEN(UNIT=17,FILE="swenwang_mag_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, swenwang_mag_temp file not opened properly------------"
    endif
    OPEN(UNIT=18,FILE="wolff_mag_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, wolff_mag_temp file not opened properly------------"
    endif
end subroutine


subroutine closetextfiles
    CLOSE(UNIT=15)
    CLOSE(UNIT=16)

    CLOSE(UNIT=17)
    CLOSE(UNIT=18)
end subroutine

end program ising
