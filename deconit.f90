subroutine deconit(utr, wtr, nt, dt, tshift, f0, &
                   maxit, minderr, rfi, it, rms)
    use interpolation_mod, only: npow2, phs_shift
    use MKL_FFT

    implicit none
    Type(CLS_FFT) :: FFT 
    integer :: nt, maxit, nft, it, maxlag, i, j
    integer, dimension(1) :: i1
    real, dimension(nt), intent(in) ::utr, wtr
    real, dimension(maxit), intent(out) ::rms
    real, dimension(nt), intent(out) ::rfi
    real :: dt, f0, tshift, minderr, df, powerU, &
                        sumsq_i, sumsq, d_error, amp
    real, dimension(:), allocatable :: uflt, &
                                                   wflt, rflt, rw, p0, pflt
    complex, dimension(:), allocatable :: wf, gauss
    real, dimension(nt) :: pflt_re

    nft=nt
    call npow2(nft)
    if (.not. allocated(gauss))  allocate(gauss(nft))
    if (.not. allocated(uflt))  allocate(uflt(nft))
    if (.not. allocated(wflt))  allocate(wflt(nft))
    if (.not. allocated(rflt))  allocate(rflt(nft))
    if (.not. allocated(wf))  allocate(wf(nft))
    if (.not. allocated(rw))  allocate(rw(nft))
    if (.not. allocated(p0))  allocate(p0(nft))
    if (.not. allocated(pflt))  allocate(pflt(nft))
    ! if (.not. allocated(pflt_re))  allocate(pflt_re(nft))

    call FFT%Create(nft)
    call gaussfilter(dt, nft, f0, gauss)
    
    uflt(1:nt) = utr
    wflt(1:nt) = wtr

    ! wf = FFT%Forward(real(wflt))
    
    call filter(uflt, nft, dt, fft, gauss)
    call filter(wflt, nft, dt, fft, gauss)
    wf = FFT%Forward(real(wflt))
    rflt = uflt

    powerU = sum(uflt ** 2)

    sumsq_i = 1.0
    d_error = 100 * powerU + minderr
    maxlag = 0.5 * nft
               
    iter: do i = 1,maxit
        call correl(rflt, wflt, nft, fft, rw)
        ! if (i==1) call print_array(rw, nft)
        rw = rw / sum(wflt ** 2)
        i1 = maxloc(abs(rw(1:maxlag)))
        amp = rw(i1(1)) / dt

        ! print *, i, i1, amp
        p0(i1(1)) = p0(i1(1)) + amp
        pflt = p0
        call filter(pflt, nft, dt, fft, gauss)
        call filter(pflt, nft, dt, fft, wf)
        rflt = uflt - pflt
        sumsq = sum(rflt ** 2) / powerU
        rms(i) = sumsq
        d_error = 100 * (sumsq_i - sumsq)
        if (abs(d_error) < minderr) exit
        sumsq_i = sumsq

    end do iter
    it = i
    call filter(p0, nft, dt, fft, gauss)
    ! call print_array(p0, nft)
    call phs_shift(p0, tshift, nft, nft, dt)
    call fft%Destory()
    rfi = p0(1:nt)
    ! call print_array(p0, nft)

end subroutine deconit

program make_rf    
    use interpolation_mod, only : npow2
    use MKL_FFT
    Type(CLS_FFT) :: FFT 

    real :: start, finish
    integer, parameter :: LNPT = 15, NPT = 2**LNPT, NDIM = 200000, maxit = 400
    double precision, parameter :: f0 = 2.0, tshift = 10.0, minderr = 0.001
    double precision :: b1,dt1
    integer :: npt1, nerr, nft, i, it, npow
    real, dimension(maxit) :: rms
    double precision, dimension(NDIM) :: datar, dataz 
    real, dimension(NDIM) :: data1
    double precision, dimension(:), allocatable :: data2, corr
    real, dimension(:), allocatable :: rfi
    ! character(len=*) :: datafile
    complex, dimension(:),  allocatable :: gauss

    call drsac1('test_R.sac', datar, npt1, b1, dt1, NDIM)
    call drsac1('test_Z.sac', dataz, npt1, b1, dt1, NDIM)
    if (.not. allocated(rfi))  allocate(rfi(npt1))
    
    call deconit(real(datar), real(dataz), npt1, real(dt1), &
                 real(tshift), real(f0), &
                 maxit, real(minderr), rfi, it, rms)
    call dwsac1('rf.sac', rfi, npt1, -tshift, dt1)
    
end program make_rf