module rttddft_io_serial
  use m_getunit, only: getunit
  use mod_mpi_env, only: mpiinfo
  use precision, only: i32, long_int, dp

  implicit none
  
  private
  
  public :: read_array, write_array

  interface read_array
    module procedure read_array_rank4
  end interface

  interface write_array
    module procedure write_array_rank4
  end interface

contains
  subroutine read_array_rank4( file_name, first, array, mpi_env )
    character(len=*), intent(in) :: file_name
    integer(i32), intent(in) :: first
    complex(dp), intent(out) :: array(:, :, :, first:)
    type(mpiinfo), intent(in), optional :: mpi_env
    
    integer(i32) :: i, last, unit
    integer(long_int) :: size_block

    last = ubound( array, 4 )
    call getunit( unit )
    inquire( ioLength=size_block ) array(:, :, :, first)
    open( unit, file=trim(file_name), action='READ', &
      & form='UNFORMATTED', access='DIRECT', recl=size_block )
    do i = first, last
      read( unit, rec=i ) array(:, :, :, i)
    end do
    close( unit )
  end subroutine

  subroutine write_array_rank4( file_name, first, array, mpi_env )
    character(len=*), intent(in) :: file_name
    integer(i32), intent(in) :: first
    complex(dp), intent(in) :: array(:, :, :, first:)
    type(mpiinfo), intent(in), optional :: mpi_env
    
    integer(i32) :: i, last, unit
    integer(long_int) :: size_block

    last = ubound( array, 4 )
    call getunit( unit )
    inquire( ioLength=size_block ) array(:, :, :, first)
    open( unit, file=trim(file_name), action='WRITE',&
        & form='UNFORMATTED', access='DIRECT', recl=size_block)
    do i = first, last
      write( unit, rec=i ) array(:, :, :, i)
    end do
    close( unit )
  end subroutine

end module