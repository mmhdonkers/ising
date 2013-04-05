module model
  implicit none
  
  private growcluster, tryadd, backtrack
  public metropolis, swenwang, wolff, calcenergy

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

  subroutine swenwang(spin, SIZE, temp, Nsum)
    !! Passed parameters, intent(in) parameters cannot be altered
    integer,intent(in) :: SIZE
    real(8),intent(in) :: temp
    integer,intent(inout) :: spin(0:SIZE-1,0:SIZE-1), Nsum

    !! Subroutine variable declerations
    integer :: ix, iy, newspin, Nclus
    real(8) :: r1(0:SIZE-1,0:SIZE-1), r2(0:SIZE-1,0:SIZE-1), r3
    logical :: bonds(0:SIZE-1,0:SIZE-1,2), mark(0:SIZE-1,0:SIZE-1)

    call random_number(r1)
    call random_number(r2)

    bonds(:,:,1) = (r1 < (1 - exp(-2/temp))*spin*cshift(spin,shift=1,dim=2))
    bonds(:,:,2) = (r2 < (1 - exp(-2/temp))*spin*cshift(spin,shift=1,dim=1))

    mark = .FALSE.
    Nsum = 0

    do ix = 0, SIZE - 1
      do iy = 0, SIZE - 1
        Nclus = 0
        call random_number(r3)
        newspin = (2 * nint(r3) - 1) * spin(ix, iy)
        call backtrack(spin, SIZE, ix, iy, bonds, mark, newspin, Nclus)
        Nsum = Nsum + Nclus**2
      enddo
    enddo
  end subroutine

  recursive subroutine backtrack(spin, SIZE, ix, iy, bonds, mark, newspin, Nclus)
    !! Passed parameters, intent(in) parameters cannot be altered
    integer,intent(in) :: ix, iy, SIZE, newspin
    logical,intent(in) :: bonds(0:SIZE-1, 0:SIZE-1, 2)
    integer,intent(inout) :: spin(0:SIZE-1, 0:SIZE-1), Nclus
    logical,intent(inout) :: mark(0:SIZE-1, 0:SIZE-1)

    if (.NOT. mark(ix, iy)) then
      mark(ix, iy) = .TRUE.
      spin(ix, iy) = newspin
      Nclus = Nclus + 1
      if (bonds(ix, iy, 1)) call backtrack(spin, SIZE, modulo(ix + 1, SIZE), iy, bonds, mark, newspin, Nclus)
      if (bonds(ix, iy, 2)) call backtrack(spin ,SIZE, ix, modulo(iy - 1, SIZE), bonds, mark, newspin, Nclus)
      if (bonds(modulo(ix - 1, SIZE), iy, 1)) call backtrack(spin, SIZE, modulo(ix - 1, SIZE), iy, bonds, mark, newspin, Nclus)
      if (bonds(ix, modulo(iy + 1, SIZE), 2)) call backtrack(spin, SIZE, ix, modulo(iy + 1, SIZE), bonds, mark, newspin, Nclus)
    endif
  end subroutine

  subroutine wolff(spin, SIZE, temp, Nswc)
    !! Passed parameters, intent(in) parameters cannot be altered
    integer,intent(in) :: SIZE
    real(8),intent(in) :: temp
    integer,intent(inout) :: spin(0:SIZE-1,0:SIZE-1), Nswc

    !! Subroutine variable declerations
    integer :: ix, iy
    real(8) :: r1, r2
    logical :: mark(0:SIZE-1,0:SIZE-1)

    call random_number(r1)
    call random_number(r2)

    ix = floor(r1*SIZE)
    iy = floor(r2*SIZE)

    mark = .FALSE.

    Nswc = 0
    call growcluster(spin, SIZE, ix, iy, temp, mark, -spin(ix, iy), Nswc)
  end subroutine

  subroutine growcluster(spin, SIZE, ix, iy, temp, mark, newspin, Nswc)
    !! Passed parameters, intent(in) parameters cannot be altered
    integer,intent(in) :: SIZE, newspin
    real(8),intent(in) :: temp
    integer :: spin(0:SIZE-1, 0:SIZE-1), ix, iy, Nswc
    logical,intent(inout) :: mark(0:SIZE-1, 0:SIZE-1)

    Nswc = Nswc + 1
    spin(ix, iy) = newspin
    mark(ix, iy) = .TRUE.
    if (.NOT.(mark(modulo(ix + 1, SIZE), iy))) call tryadd(spin, SIZE, modulo(ix + 1, SIZE), iy, temp, mark, newspin, Nswc)
    if (.NOT.(mark(modulo(ix - 1, SIZE), iy))) call tryadd(spin, SIZE, modulo(ix - 1, SIZE), iy, temp, mark, newspin, Nswc)
    if (.NOT.(mark(ix, modulo(iy + 1, SIZE)))) call tryadd(spin, SIZE, ix, modulo(iy + 1, SIZE), temp, mark, newspin, Nswc)
    if (.NOT.(mark(ix, modulo(iy - 1, SIZE)))) call tryadd(spin, SIZE, ix, modulo(iy - 1, SIZE), temp, mark, newspin, Nswc)
  end subroutine

  subroutine tryadd(spin, SIZE, ix, iy, temp, mark, newspin, Nswc)
    !! Passed parameters, intent(in) parameters cannot be altered
    integer,intent(in) :: SIZE, newspin
    real(8),intent(in) :: temp
    integer :: spin(0:SIZE-1, 0:SIZE-1), ix, iy, Nswc
    logical,intent(inout) :: mark(0:SIZE-1, 0:SIZE-1)

    !! Subroutine variable declerations
    real(8) :: r1

    if (spin(ix, iy) /= newspin) then
      call random_number(r1)
      if (r1 < (1 - exp(-2 / temp))) then
        call growcluster(spin, SIZE, ix, iy, temp, mark, newspin, Nswc)
      end if
    end if
  end subroutine

  real(8) function calcenergy(spin, SIZE, temp) result(ener)
    !! Passed paramterers, intent(in) parameters cannot be altered
    integer,intent(in) :: SIZE, spin(0:SIZE-1, 0:SIZE-1)
    real(8),intent(in) :: temp

    !! Subroutine variable declerations
    integer :: ix, iy
    real(8) :: J

    J = 1/temp

    ener = 0
    do ix = 0, SIZE - 1
      do iy = 0, SIZE - 1
        ener = ener - J * spin(ix, iy) * spin(modulo(ix - 1, SIZE), iy)
        ener = ener - J * spin(ix, iy) * spin(modulo(ix + 1, SIZE), iy)
        ener = ener - J * spin(ix, iy) * spin(ix, modulo(iy - 1, SIZE))
        ener = ener - J * spin(ix, iy) * spin(ix, modulo(iy + 1, SIZE))
      end do
    end do
  end function
end module
