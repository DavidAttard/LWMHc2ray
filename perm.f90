module perm
!! \brief Based upon http://ftp.cac.psu.edu/pub/ger/fortran/hdk/byterand.for,
!! which according to the comments in that file, is an adaptation of
!! Knuth's algorithm in his Volume 2, algorithm 3.4.2P page 145.
!! Randomly permutates integer array. Choose one among permi, permi4, permi8.
!!
!! \b Author: Kyungjin Ahn
!!
!! \b Date: 22-Nov-2009
!!

contains
!------------------------------------------------------
  subroutine permi(N, p)
    ! For random permutation of unspecified-byte integer array

    integer, intent(in) :: N
    integer, dimension(:), intent(inout) :: p

    integer :: i
    integer :: k, j, ipj, itemp, m
    real(kind=8), dimension(100) :: u


    ! Generate up to 100 U(0,1) numbers at a time.
    do i=1,N,100
       m = min(N-i+1, 100)
       call random_number(u)
       do j=1,m
          ipj = i+j-1
          k = int(u(j)*(N-ipj+1)) + ipj
          itemp = p(ipj)
          p(ipj) = p(k)
          p(k) = itemp
       end do
    end do
    return

  end subroutine permi


!------------------------------------------------------
  subroutine permi4(N, p)
    ! For random permutation of 4 byte integer array

    integer(kind=4), intent(in) :: N
    integer(kind=4), dimension(:), intent(inout) :: p

    integer(kind=4) :: i
    integer(kind=4) :: k, j, ipj, itemp, m
    real(kind=8), dimension(100) :: u


    ! Generate up to 100 U(0,1) numbers at a time.
    do i=1,N,100
       m = min(N-i+1, 100)
       call random_number(u)
       do j=1,m
          ipj = i+j-1
          k = int(u(j)*(N-ipj+1)) + ipj
          itemp = p(ipj)
          p(ipj) = p(k)
          p(k) = itemp
       end do
    end do
    return

  end subroutine permi4


!------------------------------------------------------
  subroutine permi8(N, p)
    ! For random permutation of 8 byte integer array

    integer(kind=4), intent(in) :: N
    integer(kind=8), dimension(:), intent(inout) :: p

    integer(kind=4) :: i
    integer(kind=4) :: k, j, ipj, itemp, m
    real(kind=8), dimension(100) :: u


    ! Generate up to 100 U(0,1) numbers at a time.
    do i=1,N,100
       m = min(N-i+1, 100)
       call random_number(u)
       do j=1,m
          ipj = i+j-1
          k = int(u(j)*(N-ipj+1)) + ipj
          itemp = p(ipj)
          p(ipj) = p(k)
          p(k) = itemp
       end do
    end do
    return

  end subroutine permi8

end module perm
