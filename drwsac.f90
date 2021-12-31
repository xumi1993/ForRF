subroutine drsac1(datafile,data,npt1,b1,dt1,NDIM)
    ! read sac file and convert to double precision

    implicit none
    character(len=*),intent(in) :: datafile
    integer :: NDIM
    real, dimension(NDIM) :: dat_sngl
    double precision, dimension(NDIM), intent(out) :: data
    integer :: npt1, nerr
    real :: b1_sngl,dt1_sngl
    double precision :: b1,dt1

    ! read file as single precision
    call rsac1(datafile,dat_sngl,npt1,b1_sngl,dt1_sngl,NDIM,nerr)
    if (nerr > 0) then
       print *, 'Error reading sac file', trim(datafile)
       stop
    endif

    ! return double precision quantities
    b1 = dble(b1_sngl)
    dt1 = dble(dt1_sngl)
    data = dble(dat_sngl)

end subroutine drsac1

  subroutine dwsac1(datafile,data,npt1,b1,dt1)
    ! convert to single precision, then write sac file
    ! --> includes an option to add minmax values to sac file,
    !     which are used in the plotting scripts

    implicit none
    character(len=*),intent(in) :: datafile
    integer, intent(in) :: npt1
    double precision, dimension(npt1), intent(in) :: data
    double precision, intent(in) :: b1,dt1
    logical, parameter :: minmax_header = .true.

    real, dimension(npt1) :: dat_sngl,ti_sngl
    real :: b1_sngl,dt1_sngl,xmin_sngl,xmax_sngl
    integer :: nerr,i

    ! convert to single precision
    b1_sngl = real(b1)
    dt1_sngl = real(dt1)
    dat_sngl = real(data)

    if (minmax_header) then
       ! get time vector
       ti_sngl = 0.
       do i = 1,npt1
          ti_sngl(i) = b1_sngl + (i-1)*dt1_sngl
       enddo

       !call newhdr()  ! create a new header

       ! set minmax values in sac file
       xmin_sngl = minval(dat_sngl)
       xmax_sngl = maxval(dat_sngl)
       call setfhv('depmin',xmin_sngl,nerr)
       call setfhv('depmax',xmax_sngl,nerr)

       call setnhv('npts',npt1,nerr)          ! sets number of points
       !call setfhv('b',ti_sngl(1),nerr)       ! sets begin
       !call setfhv('e',ti_sngl(npt1),nerr)    ! sets end
       !call setlhv('leven',.false.,nerr)        ! sets un-even sampling
       !call setihv('iftype','itime',nerr)          ! sets file type: time file

       ! write file with headers (LQY: bug with b in wsac0())
       ! call wsac0(datafile,ti_sngl,dat_sngl,nerr)
       call wsac1(datafile,dat_sngl,npt1,b1_sngl,dt1_sngl,nerr)
    else
       call wsac1(datafile,dat_sngl,npt1,b1_sngl,dt1_sngl,nerr)
    endif
    if (nerr > 0) then
        print *, 'Error writing sac file', trim(datafile)
        stop
    endif

  end subroutine dwsac1