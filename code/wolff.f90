module wolff
  implicit none
  
  public wolff

contains

  subroutine wolff(wolffspin, SIZE)
 
   !! Passed parameters, intent(in) parameters cannot be altered
    integer,intent(in) :: SIZE 
    integer,intent(inout) :: wolffspin(0:SIZE-1,0:SIZE-1)
 
    !! Subroutine variable declerations
    integer :: ix,iy, xcount=0, ycount=0, xx, yy
    real(8) :: r4,r5,r6
    integer :: cluster(0:2,0:SIZE**2)  ! 1-D array for storing the xy locations of spins in a cluster
 
    !! Randomly choose a location in the array
    call random_number(r4)
    call random_number(r5)
 
    ix = floor((r4*SIZE))
    iy = floor((r5*SIZE))




!!!!! Grow cluster (Wolff cluster growth procedure)

!+++++++++++++++++++++++++++++++++++++++++++++ BEGIN LOOP ++++++++++++++++++++++++++++++++++++++++


    ! cluster must be re-initialized for each loop (initialized to any number that isn't 1 or -1
    cluster = 0 


    !! Find a neighboring spin with the same orientation (if one exists)

    select case (spin(ix,iy))
       case (spin(ix+1,iy+1))
          xx = ix+1
          yy = iy+1
          exit
       case (spin(ix+1,iy))
          xx = ix+1
          yy = iy
          exit
       case (spin(ix+1,iy-1))
          xx = ix+1
          yy = iy-1
          exit
       case (spin(ix,iy+1))
          xx = ix
          yy = iy+1
          exit
       case (spin(ix,iy-1))
          xx = ix
          yy = iy-1
          exit
       case (spin(ix-1,iy+1))
          xx = ix-1
          yy = iy+1
          exit
       case (spin(ix-1,iy))
          xx = ix-1
          yy = iy
          exit
       case (spin(ix-1,iy-1))
          xx = ix-1
          yy = iy-1
          exit
    end select



    !! possibly add spin to cluster

    call random_number(r6)

    if (1-exp(2/temp) > r6) then
       cluster(xcount,ycount) = spin(xx,yy)
       xcount = xcount+1
       ycount = ycount+1
    endif


 
    !! Flip cluster

    if (cluster == 1) then
! flip to -1
    elseif (cluster == -1) then
! flip to 1
    elseif
       print*, "Some spin values in the cluster are neither 1 nor -1"
    endif

!+++++++++++++++++++++++++++++++ END LOOP +++++++++++++++++++++++++++++++++++++++++


  end subroutine
end module

