module hdf5_group
  use hdf5_utils

  implicit none 

  private
  public :: hdf5_create_group, hdf5_link_exists, hdf5_delete_link

  
  contains


  !> Create a new `group` at `h5path`, assiciated to a link in an HDF5 file or group associated with `h5id`.
  subroutine hdf5_create_group(h5id_file, h5path, group)
    !> Identifier of the file or group, used by HDF5.
    integer(hdf5_id), intent(in) :: h5id_file
    !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
    character(*), intent(in) :: h5path
    !> Name of the group to create.
    character(*), intent(in) :: group
#ifdef _HDF5_
    integer :: h5err
    integer(hdf5_id) :: h5id_group, h5id_newgroup

    call h5gopen_f(h5id_file, h5path, h5id_group, h5err)
    call handle_hdf5_error('h5gopen_f' ,h5err)

    call h5gcreate_f(h5id_group, group, h5id_newgroup, h5err)
    call handle_hdf5_error('h5gcreate_f' ,h5err)

    call h5gclose_f(h5id_newgroup, h5err)
    call handle_hdf5_error('h5gclose_f', h5err)

    call h5gclose_f(h5id_group, h5err)
    call handle_hdf5_error('h5gclose_f', h5err)
#endif    
  end subroutine hdf5_create_group


  !> Check if `h5path` is associated with a link in an HDF5 file or group corresponding to `h5id`.
  logical function hdf5_link_exists(h5id_file, h5path)
    !> Identifier of the file or group, used by HDF5.
    integer(hdf5_id), intent(in) :: h5id_file
    !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
    character(*), intent(in) :: h5path
#ifdef _HDF5_
    integer :: h5err

    call h5lexists_f(h5id_file, h5path, hdf5_link_exists, h5err)
    call handle_hdf5_error('h5lexists_f' ,h5err)
#else 
      hdf5_link_exists = .false.
#endif     
  end function hdf5_link_exists


  !> Delete object associated with `h5path` in an HDF5 file or group corresponding to `h5id`.
  subroutine hdf5_delete_link(h5id_file, h5path)
    !> Identifier of the file or group, used by HDF5.
  integer(hdf5_id), intent(in) :: h5id_file
    !> Path to an object in an HDF5 file. Can be given as `path/to/object`.
    character(*), intent(in) :: h5path
#ifdef _HDF5_
    integer :: h5err

    call h5ldelete_f(h5id_file, h5path, h5err)
    call handle_hdf5_error('h5ldelete_f' ,h5err)
#endif    
  end subroutine hdf5_delete_link

end module hdf5_group  