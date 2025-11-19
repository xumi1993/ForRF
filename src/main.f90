program main
    use sacio
    use decon_mod
    implicit none

    type(sachead) :: head_r, head_z, head_out
    double precision, allocatable :: data_r(:), data_z(:), data_rf(:)
    integer :: flag, npts
    real :: dt, tshift, f0, minderr
    integer :: maxit, ipart
    character(len=256) :: file_r, file_z, file_out

    ! Parameters
    tshift = 10.0
    f0 = 2.0
    maxit = 200
    minderr = 0.001
    ipart = 0

    ! Get arguments if needed, but for now hardcode or use local files
    file_r = "test_R.sac"
    file_z = "test_Z.sac"
    file_out = "test_RF.sac"

    print *, "Reading ", trim(file_r)
    call sacio_readsac(trim(file_r), head_r, data_r, flag)
    if (flag /= 0) then
        print *, "Error reading ", trim(file_r)
        stop
    endif

    print *, "Reading ", trim(file_z)
    call sacio_readsac(trim(file_z), head_z, data_z, flag)
    if (flag /= 0) then
        print *, "Error reading ", trim(file_z)
        stop
    endif

    npts = head_r%npts
    dt = head_r%delta

    if (head_z%npts /= npts) then
        print *, "Error: npts mismatch"
        stop
    endif

    allocate(data_rf(npts))

    print *, "Calculating Receiver Function..."
    ! deconit(utr, wtr, dt, tshift, f0, maxit, minderr, ipart, rfi)
    call deconit(data_r, data_z, dt, tshift, f0, maxit, minderr, ipart, data_rf)

    head_out = head_r
    ! Update header info if necessary
    ! For now just write it out

    print *, "Writing ", trim(file_out)
    call sacio_writesac(trim(file_out), head_out, data_rf, flag)
    if (flag /= 0) then
        print *, "Error writing ", trim(file_out)
        stop
    endif

    print *, "Done."

end program main
