program main
    use sacio
    use decon_mod
    implicit none

    type(sachead) :: head_r, head_z, head_out
    double precision, allocatable :: data_r(:), data_z(:), data_rf(:)
    integer :: flag, npts
    real :: dt, tshift, f0, minderr
    integer :: maxit, ipart
    character(len=256) :: file_r, file_z, file_out, arg
    integer :: num_args
    logical :: use_gpu
    integer(8) :: t1, t2, rate
    real :: time_sec

    ! Parameters
    tshift = 10.0
    f0 = 2.0
    maxit = 200
    minderr = 0.001
    ipart = 0
    use_gpu = .true.

    ! Parse arguments
    num_args = command_argument_count()
    if (num_args > 0) then
        call get_command_argument(1, arg)
        if (trim(arg) == 'cpu' .or. trim(arg) == 'CPU') then
            use_gpu = .false.
        endif
    endif

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

    call system_clock(t1, rate)
    if (use_gpu) then
        print *, "Calculating Receiver Function (GPU)..."
        call deconit_gpu(data_r, data_z, dt, tshift, f0, maxit, minderr, ipart, data_rf)
    else
        print *, "Calculating Receiver Function (CPU)..."
        call deconit_cpu(data_r, data_z, dt, tshift, f0, maxit, minderr, ipart, data_rf)
    endif
    call system_clock(t2)
    
    time_sec = real(t2 - t1) / real(rate)
    print *, "Calculation time: ", time_sec, " seconds"



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
