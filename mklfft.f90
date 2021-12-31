Module MKL_FFT !// 以下代码复制粘贴，直接使用
  Use MKL_DFTI
  private
  Type , public :: CLS_FFT
    type(DFTI_DESCRIPTOR), Pointer :: h => NULL()
    Integer :: Err, nft
  contains
    Procedure :: Create
    Procedure :: Forward
    Procedure :: Backward
    Procedure :: Destory
  End Type CLS_FFT

contains

  Subroutine Create( this , N )
    class( CLS_FFT ) :: this
    Integer , Intent( IN ) :: N
    this%nft = N
    this%Err = DftiCreateDescriptor( this%h , DFTI_SINGLE , DFTI_REAL , 1 , N )
    this%Err = DftiSetValue( this%h , DFTI_PLACEMENT , DFTI_NOT_INPLACE )
    this%Err = DftiCommitDescriptor( this%h )
  End Subroutine Create

  Function Forward( this , X ) result( F )
    class( CLS_FFT ) :: this
    Real :: X(this%nft)
    Complex :: F(this%nft/2+1 )
    this%Err = DftiComputeForward( this%h , X , F )
  End Function Forward

  Function Backward( this , X ) result( T )
    class( CLS_FFT ) :: this
    Complex :: X(this%nft)
    Real    :: T(this%nft*2-1 )
    this%Err = DftiComputeBackward( this%h , X , T )
  End Function Backward

  Subroutine Destory( this )
    class( CLS_FFT ) :: this
    this%Err = DftiFreeDescriptor( this%h )
  End Subroutine Destory

End Module MKL_FFT


! program mklfft
!     use interpolation_mod, only : npow2
!     ! implicit none
!     use MKL_FFT
!     Type(CLS_FFT) :: FFT 
!     integer, parameter :: LNPT = 15, NPT = 2**LNPT, NDIM = 200000, maxit = 200
!     double precision, parameter :: f0 = 2.0, tshift = 10.0, minderr = 0.001
!     double precision :: b1,dt1
!     integer :: npt1, nerr, nft, i, it, npow
!     double precision, dimension(maxit) :: rms
!     double precision, dimension(NDIM) :: datar, dataz 
!     real, dimension(:), allocatable :: data2
!     complex, dimension(:), allocatable :: data1
!     double precision, dimension(:), allocatable :: gauss
    
!     call drsac1('test_R.sac', datar, npt1, b1, dt1, NDIM)
!     call drsac1('test_Z.sac', dataz, npt1, b1, dt1, NDIM)
!     nft = npt1
!     call npow2(nft)
!     if (.not. allocated(data2))  allocate(data2(nft))
!     if (.not. allocated(data1))  allocate(data1(nft))
!     if (.not. allocated(gauss))  allocate(gauss(nft))

!     ! data1 = cmplx(datar)
!     call FFT%Create(nft)
!     ! data1 = FFT%Forward(real(datar))
!     ! data2  = FFT%Backward(data1/nft)
!     ! Do i = 1 , nft
!     !     Write(*,*) i , data2( i ) !// 输出频率域
!     ! end Do
!     ! call FFT%Destory()
!     ! do i=1,nft
!     !     print *, data2(i)
!     ! enddo
!     ! call savetxt('mfft.txt', data1, nft)
!     call gaussfilter(dt1, nft, f0, gauss)
!     call filter(datar, nft, dt1, FFT, cmplx(gauss, 0.))
!     call FFT%Destory()
!     ! do i=1,nft
!     !     print *, datar(i)
!     ! enddo

! end program mklfft
