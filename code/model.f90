module model
  implicit none
  
  private growcluster, tryadd
  public metropolis, swenwang

contains

  subroutine metropolis(spin, SIZE, weight)
    !! Passed parameters, intent(in) parameters cannot be altered
    integer,intent(in) :: SIZE 
    real(8),intent(in) :: weight(-2:2)
    integer,intent(inout) :: spin(0:SIZE-1,0:SIZE-1)
 
    !! Subroutine variable declerations
    integer :: ix,iy
    real(8) :: r1,r2,r3
    real(8) :: expo
    real(8) :: st, sb, sl, sr !neighbour cells
 
    !! Randomly choose a location in the array (the old spin)
    call random_number(r1)
    call random_number(r2)
 
    ix = floor(r1*SIZE)
    iy = floor(r2*SIZE)
 
    !! Calculate energy due to the neighbors 
    !!  (the modulo takes into account the boundaries)
    sl = spin(modulo(ix-1,SIZE), iy)
    sr = spin(modulo(ix+1,SIZE), iy)
    st = spin(ix, modulo(iy-1,SIZE))
    sb = spin(ix, modulo(iy+1,SIZE))
 
    expo = weight(-nint(0.5 * spin(ix,iy) * (sl + sr + st + sb)))

    !! Metropolis test
    call random_number(r3)
 
    if (expo > r3) then
      spin(ix,iy) = -spin(ix,iy)
    endif
  end subroutine

  subroutine swenwang(spin, SIZE, temp)
    !! Passed parameters, intent(in) parameters cannot be altered
    integer,intent(in) :: SIZE
    real(8),intent(in) :: temp
    integer,intent(inout) :: spin(0:SIZE-1,0:SIZE-1)

    !! Subroutine variable declerations
    real(8) :: r1(0:SIZE-1,0:SIZE-1),r2(0:SIZE-1,0:SIZE-1)
    logical :: bonds(0:SIZE-1,0:SIZE-1,2)
    logical :: mark(0:SIZE-1,0:SIZE-1)

    call random_number(r1)
    call random_number(r2)

    bonds(:,:,1) = (r1 < (1 - exp(-2/temp))*spin*cshift(spin,shift=1,dim=2))
    bonds(:,:,2) = (r2 < (1 - exp(-2/temp))*spin*cshift(spin,shift=1,dim=1))

  end subroutine

  subroutine growcluster()

  end subroutine

  subroutine tryadd()

  end subroutine
end module
