module autormt_test_drivers
    use modmpi, only: mpiinfo
    ! Load test drivers here
    use init_autormt_test, only: init_autormt_test_driver
    
    private
    public :: autormt_test_driver
  
  contains
    subroutine autormt_test_driver(mpiglobal, kill_on_failure)
      !> mpi information
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program before the test driver finishes
      !> if an assertion fails 
      logical, optional :: kill_on_failure 
  
      ! Call test drivers here
      call init_autormt_test_driver(mpiglobal, kill_on_failure)
      
    end subroutine autormt_test_driver
  
  end module autormt_test_drivers