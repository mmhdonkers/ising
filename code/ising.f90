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
  integer,parameter :: TEMPFINAL = 350, TIMEFINAL = 100000

!! fortran begins indexing from 1. Start it from 0 because the rand() starts from 0
  integer :: spin(0:SIZE-1,0:SIZE-1), temp, time
  real(8) :: weight(-2:2), spintotal, sqspin, qspin, etot, sqetot, absspin
  integer :: wolffspin(0:SIZE-1, 0:SIZE-1), Nswc
  real(8) :: wolfftotal, sqwol, qwol, ewol, sqewol, magsupwol
  integer :: swenwangspin(0:SIZE-1, 0:SIZE-1), Nsum
  real(8) :: swenwangtotal, sqsw, qsw, esw, sqesw, magsupsw

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


  do temp = 150,TEMPFINAL

! Re-initialize Wolff lattice. All values in the lattice must be 1
    wolffspin(:,:) = 1
    wolfftotal = 0
    sqwol = 0
    qwol = 0
    ewol = 0
    sqewol = 0
    magsupwol = 0
    do time = 1, TIMEFINAL, 100
      call wolff(wolffspin, SIZE, temp/100d0, Nswc)
      wolfftotal = wolfftotal + abs(sum(wolffspin)/(SIZE**2*1d0))
      sqwol = sqwol + (sum(wolffspin)/(SIZE**2*1d0))**2
      qwol = qwol + (sum(wolffspin)/(SIZE**2*1d0))**4
      ewol = ewol + calcenergy(wolffspin, SIZE, temp/100d0)
      sqewol = sqewol + calcenergy(wolffspin, SIZE, temp/100d0)**2
      magsupwol = magsupwol + Nswc**2/(SIZE**2*1d0)
    end do

    wolfftotal = 100d0 * wolfftotal / TIMEFINAL
    sqwol = 100d0 * sqwol / TIMEFINAL
    qwol = 100d0 * qwol / TIMEFINAL
    ewol = 100d0 * ewol / TIMEFINAL
    sqewol = 100d0 * sqewol / TIMEFINAL
    magsupwol = 100d0 * magsupwol / TIMEFINAL
    
! Re-initialize Swendswon-Wang lattice. All values in the lattice must be 1
    swenwangspin(:,:) = 1
    swenwangtotal = 0
    sqsw = 0
    qsw = 0
    esw = 0
    sqesw = 0
    magsupsw = 0
    do time = 1, TIMEFINAL, 100
      call swenwang(swenwangspin, SIZE, temp/100d0, Nsum)
      swenwangtotal = swenwangtotal + abs(sum(swenwangspin)/(SIZE**2*1d0))
      sqsw = sqsw + (sum(swenwangspin)/(SIZE**2*1d0))**2
      qsw = qsw + (sum(swenwangspin)/(SIZE**2*1d0))**4
      esw = esw + calcenergy(swenwangspin, SIZE, temp/100d0)
      sqesw = sqesw + calcenergy(swenwangspin, SIZE, temp/100d0)**2
      magsupsw = magsupsw + Nsum
    end do

    swenwangtotal = 100d0 * swenwangtotal / TIMEFINAL
    sqsw = 100d0 * sqsw / TIMEFINAL
    qsw = 100d0 * qsw / TIMEFINAL
    esw = 100d0 * esw / TIMEFINAL
    sqesw = 100d0 * sqesw / TIMEFINAL
    magsupsw = 100d0 * magsupsw / (TIMEFINAL * SIZE**2 * 1d0)

! Re-initialize the Metropolis lattice. All values in the lattice must be 1
    spin(:,:) = 1
    spintotal = 0       
    sqspin = 0
    qspin = 0
    etot = 0
    sqetot = 0
    absspin = 0
! Only 5 options for the exponent for metropolis so calculate them once
    weight = [exp(-800d0/temp),exp(-400d0/temp),1d0,exp(400d0/temp),exp(800d0/temp)]
    do time = 0,TIMEFINAL
      call metropolis(spin, SIZE, weight)
      spintotal = spintotal + abs(sum(spin)/(SIZE**2*1d0))
      sqspin = sqspin + (sum(spin)/(SIZE**2*1d0))**2
      qspin = qspin + (sum(spin)/(SIZE**2*1d0))**4
      etot = etot + calcenergy(spin, SIZE, temp/100d0)
      sqetot = sqetot + calcenergy(spin, SIZE, temp/100d0)**2
      absspin = absspin + sum(abs(spin))
! We want time to print only once (choose an arbitrary temperature)
      if (temp == 250) then
        WRITE(16,*) sum(spin)/(SIZE**2*1d0), time
      end if
    end do
    spintotal = 1d0 * spintotal / TIMEFINAL
    sqspin = 1d0 * sqspin / TIMEFINAL
    qspin = 1d0 * qspin / TIMEFINAL
    etot = 1d0 * etot / TIMEFINAL
    sqetot = 1d0 * sqetot / TIMEFINAL
    absspin = 1d0 * absspin / TIMEFINAL
      
    call plot_spin(wolffspin, SIZE, temp/100d0)
    WRITE(15,*) spintotal, 1-qspin/(3*sqspin**2), (sqetot - etot**2)/(SIZE**2*1d0),sqspin - absspin**2, temp/100d0
    WRITE(17,*) swenwangtotal, 1-qsw/(3*sqsw**2), (sqesw - esw**2)/(SIZE**2*1d0), magsupsw, temp/100d0
    WRITE(18,*) wolfftotal, 1-qwol/(3*sqwol**2), (sqewol - ewol**2)/(SIZE**2*1d0), magsupwol, temp/100d0
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
    OPEN(UNIT=15,FILE="metrop_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, metrop_temp file not opened properly------------"
    endif
    OPEN(UNIT=16,FILE="metrop_mag_time.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, metrop_mag_time file not opened properly------------"
    endif

    OPEN(UNIT=17,FILE="swenwang_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, swenwang_temp file not opened properly------------"
    endif
    OPEN(UNIT=18,FILE="wolff_temp.txt",STATUS="REPLACE",IOSTAT=OPEN_STATUS)
    if (OPEN_STATUS /= 0) then
       STOP "------------Error, wolff_temp file not opened properly------------"
    endif
end subroutine


subroutine closetextfiles
    CLOSE(UNIT=15)
    CLOSE(UNIT=16)

    CLOSE(UNIT=17)
    CLOSE(UNIT=18)
end subroutine

end program ising
