!> Utilities for hdf5 wrapper
module hdf5_utils

#ifdef _HDF5_  
  use hdf5
#endif

  use trace, only: trace_back
  use modmpi, only: terminate
  use iso_fortran_env, only: error_unit
  use precision, only: sp, dp, str_1024

  implicit none

  !> Root group of any HDF5 file
  character(*), parameter :: hdf5_root = './'

  !> HDF5 integer kinds
#ifdef _HDF5_
  integer, parameter :: hdf5_id = HID_T
  integer, parameter :: hdf5_size = HSIZE_T
  integer, parameter :: hdf5_ssize = HSSIZE_T
#else 
  integer, parameter :: hdf5_id = sp
  integer, parameter :: hdf5_size = sp
  integer, parameter :: hdf5_ssize = sp
#endif


  !> Max. string length for write and read
  integer, parameter :: hdf5_max_string_len = str_1024
  !> Blank, to fill non needed character length
  character(1), parameter :: hdf5_blank = ' '

  contains


  !> Return HDF5 id for double
  integer(hdf5_id) function hdf5_double()
#ifdef _HDF5_
    hdf5_double = H5T_NATIVE_DOUBLE
#else
    hdf5_double = dp
#endif
  end function

  !> Return HDF5 id for real (sp)
  integer(hdf5_id) function hdf5_float()
#ifdef _HDF5_
    hdf5_float = H5T_NATIVE_REAL
#else
    hdf5_float = sp
#endif 
  end function hdf5_float


  !> Return HDF5 id for integer
  integer(hdf5_id) function hdf5_integer()
#ifdef _HDF5_    
    hdf5_integer = H5T_NATIVE_INTEGER
#else
    hdf5_integer = sp
#endif    
  end function hdf5_integer

  !> Return HDF5 id for character
  integer(hdf5_id) function hdf5_character()
#ifdef _HDF5_
    hdf5_character = H5T_NATIVE_CHARACTER
#else
    hdf5_character = sp
#endif
  end function hdf5_character


  !> Abort run if exciting was compiled exciting without HDF5.
  !> Call this routine in features that require HDF5 in the beginning.
  subroutine abort_if_not_hdf5(message)
    !> Message to
    character(*), intent(in), optional :: message

    character(*), parameter :: default_message = &
            'Error: exciting needs to be linked to HDF5 to run the current routine.'

#ifndef _HDF5_
    if (present(message)) then
      print*, message
    else
      print*, default_message
    end if 
    call trace_back()
    call terminate()
#endif 
  end subroutine abort_if_not_hdf5

  
  !> Assert if the HDF5 error flag `h5err` is zero or not. If it is not zero, the
  !> routine prints out a message that tells the routine name and the value of `h5err`.
  subroutine handle_hdf5_error(hdf5_routine_name, h5err, unit_out)
    !> Name of the hdf5 routine.
    character(*), intent(in) :: hdf5_routine_name
    !> HDF5 error flag.
    integer, intent(in) :: h5err
    !> Unit to write to. By default `error_unit`.
    integer, intent(in), optional :: unit_out

    integer :: unit_out_local

    unit_out_local = error_unit
    if (present(unit_out)) unit_out_local = unit_out

    if (h5err /= 0) then
      write(unit_out_local, '(A, A, A, I3)') 'HDF5 routine ', hdf5_routine_name, ' failed with h5err = ', h5err
      call trace_back()
      call terminate()
    end if
  end subroutine handle_hdf5_error

end module hdf5_utils  