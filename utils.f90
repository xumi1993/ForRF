subroutine filter(x, nft, dt, fftcls, gauss)
    use MKL_FFT
    implicit  none
    Type(CLS_FFT) :: fftcls 

    integer:: npow, nft
    real :: dt
    ! double precision, dimension(nft), intent(in):: gauss
    real, dimension(nft), intent(inout):: x
    real, dimension(nft) :: r
    complex, dimension(nft) :: z1, gauss
    ! z1 = cmplx(x, 0.0)

    z1 = fftcls%forward(x)
    z1 = z1 * gauss *dt
    z1 = z1/nft
    x = fftcls%backward(z1)

end subroutine filter


subroutine gaussfilter(dt, nft, f0, gauss)
    implicit none
    integer nft, nft21, i
    real :: dt, f0, df
    complex, dimension(nft), intent(out) :: gauss
    real, dimension(nft) :: gauss1

    df = 1.0 / (nft * dt)
    nft21 = 0.5 * nft + 1
    do i = 1,nft21
        gauss1(i) = exp(-0.25 * (2 * 3.14159265 * df * (i-1) / f0) ** 2)/dt
        ! print *, gauss(i) 
    enddo
    do i=nft21+1,nft
        gauss1(i) = gauss1(2*nft21-i)
    enddo
    gauss = cmplx(gauss1)

end subroutine gaussfilter


subroutine correl(r, w, nft, fftcls, rw)
    use MKL_FFT
    implicit none
    Type(CLS_FFT) :: fftcls 
    integer :: nft
    complex, dimension(nft) :: r1, w1, dat
    real, dimension(nft), intent(in):: r, w
    real, dimension(nft), intent(out):: rw
    ! real, dimension()

    r1 = fftcls%forward(r)
    w1 = fftcls%forward(w)

    dat = r1* conjg(w1)
    dat = dat / nft
    rw = fftcls%backward(dat)
       
end subroutine correl


subroutine phaseshift(x, nft, dt, tshift, fftcls)
    use MKL_FFT
    implicit none
    Type(CLS_FFT) :: fftcls 
    integer :: nft, shift_i, i
    real :: dt, tshift, p
    real, dimension(nft), intent(inout) :: x
    complex, dimension(nft) :: z1

    z1 = fftcls%forward(x)
     shift_i=int((tshift/dt)+0.5)
    do i=1,nft
        p = 2*3.14159265*i*shift_i/nft
        z1(i) = z1(i) * cmplx(cos(p), -sin(p))
    enddo
    z1 = z1/nft
    x = fftcls%backward(z1) / cos(2*3.14159265*shift_i/nft)

end subroutine phaseshift

! subroutine phaseshift( x, nft, dt, shift, fftcls)
!     use MKL_FFT
!     Type(CLS_FFT) :: fftcls 
!     complex, dimension(nft) :: Xf
!     real, dimension(nft) :: p, x
!     real :: p1
!     integer shift_i
!     pi=3.141592653
!     xf = fftcls%forward(x)
!     ! do i=1, nft/2
!     !     xf(nft-i-1) = cmplx(real(xf(i)), -aimag(xf(i)))
!     ! enddo

!     shift_i=int((shift/dt)+0.5)
!     n21 = nft/2+1
!     do i=1,n21
!         p1 = 2 * pi * i * shift_i / nft
!         Xf(i) = Xf(i) * cmplx(cos(p1), -sin(p1))
!         if(i >= 1)then
!             xf(nft + 2 - i) = conjg(xf(i))
!         endif
!     enddo
!     xf(n21) = cmplx(real(xf(n21)),0.0)
!     Xf = Xf / nft
!     x = fftcls%backward(xf)
!     x = x / cos(2 * pi * shift_i / nft)

! end subroutine 


subroutine print_array(x, nt)
    real, dimension(nt) :: x
    do i=1,nt
        print *, x(i)
    enddo
    
end subroutine print_array