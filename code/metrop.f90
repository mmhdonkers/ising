module metrop
  implicit none
  
  public metropolis

contains

  subroutine metropolis(spin, size, weight)
    !! Passed parameters, intent(in) parameters cannot be altered
    integer,intent(in) :: size 
    real(8),intent(in) :: weight(-2:2)
    integer,intent(inout) :: spin(0:size-1,0:size-1)
 
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
 
    !! Calculate energy due to the neighbors 
    !!  (the modulo takes into account the boundaries)
    sl = spin(modulo(ix-1,size), iy)
    sr = spin(modulo(ix+1,size), iy)
    st = spin(ix, modulo(iy-1,size))
    sb = spin(ix, modulo(iy+1,size))
 
    expo = weight(-nint(0.5 * spin(ix,iy) * (sl + sr + st + sb)))

    !! Metropolis test
    call random_number(r3)
 
    if (expo > r3) then
      spin(ix,iy) = -spin(ix,iy)
    endif
  end subroutine
end module
