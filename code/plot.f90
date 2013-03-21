module plot
  use plplot

  implicit none
  
  public plot_init, plot_close, plot_spin

contains

  subroutine plot_init() 
    call plsdev("xcairo")

    call plscol0(0, 255, 255, 255)  ! white
    call plscol0(1, 255, 0, 0)  ! red
    call plscol0(2, 0, 255, 0)  ! green
    call plscol0(3, 0, 0, 255)  ! blue
    call plscol0(4, 255, 0, 255)  ! magenta
    call plscol0(5, 0, 255, 255)  ! cyan
    call plscol0(6, 255, 255, 0)  ! yellow
    call plscol0(7, 0, 0, 0)  ! black
    call plscol0(8, 255, 77, 0)  ! orange
    call plscol0(9, 128, 128, 128)  ! gray

    call plinit()
  end subroutine

  subroutine plot_close()
    call plspause(.false.)
    call plend()
  end subroutine

  subroutine plot_spin(spin,size)
    integer,intent(in) :: size
    integer,intent(in) :: spin(0:size-1,0:size-1)

    integer :: i, j

    call plcol0(7)
    call plenv(0d0, size*1d0, 0d0, size*1d0, 0, 0)
    call pllab("x", "y", "spin")

    do i = 0, size - 1
      do j = 0, size - 1
        if (spin(i,j) .eq. -1) then
          call plcol0(1)
          call plpoin([i + 0.5d0], [j + 0.5d0], 31)
        else if (spin(i,j) .eq. 1) then
          call plcol0(2)
          call plpoin([i + 0.5d0], [j + 0.5d0], 30)
        end if
      end do
    end do

    call plspause(.false.)
  end subroutine

end module
