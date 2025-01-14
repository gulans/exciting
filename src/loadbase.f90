


Subroutine loadbase
    ! !USES:
    Use modinput, only: input
    Use modmain, only: nspecies,nrmt,natoms,idxas,apword,apwfr,nlorb,lofr
    Use mod_hybrids, only: ex_coef, ec_coef
    implicit none
    integer :: is,ia,ias,l,io1,ir,ilo,io2,nr
    logical :: file_exists

    inquire(file="STATE_BASE.OUT",EXIST=file_exists)
    
if (file_exists) then
    write(*,*)"Restart (HYB): loading radial basis from STATE_BASE.OUT"
    open (11, file = "STATE_BASE.OUT", status = 'old')

    !!!APW functions
    Do is = 1, nspecies
        nr = nrmt (is)
        Do ia = 1, natoms (is)
           ias = idxas (ia, is)
           Do l = 0, input%groundstate%lmaxapw
              Do io1 = 1, apword (l, is)
                do ir =1 ,nr
                    read(11,*)apwfr (ir, 1, io1, l, ias), apwfr (ir, 2, io1, l, ias)
                enddo
              enddo
            enddo
        enddo
    enddo

    !!!LO functions
    Do is = 1, nspecies
        nr = nrmt (is)
        Do ia = 1, natoms (is)
           ias = idxas (ia, is)
           Do ilo = 1, nlorb (is)
                do ir =1 ,nr
                    read(11,*)lofr (ir, 1, ilo, ias),lofr (ir, 2, ilo, ias)
                enddo
            enddo
        enddo
    enddo

    close(11)
    
else
    write(*,*)"######################################################"
    write(*,*)"##     Restart (HYB): STATE_BASE.OUT not found.     ##"
    write(*,*)"##          Generating PBE basis functions.         ##"
    write(*,*)"######################################################"
    ex_coef = 0.d0
    ec_coef = 1.d0
    
    call genapwfr()
    call genlofr(.false.)   

    ex_coef = input%groundstate%Hybrid%excoeff
    ec_coef = input%groundstate%Hybrid%eccoeff

endif

end Subroutine