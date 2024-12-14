


Subroutine storebase
    ! !USES:
    Use modinput, only: input
    Use modmain, only: nspecies,nrmt,natoms,idxas,apword,apwfr,nlorb,lofr
    implicit none
    integer :: is,ia,ias,l,io1,ir,ilo,io2,nr


    open (11, file = "STATE_BASE.OUT", status = 'replace')

    !!!APW functions
    Do is = 1, nspecies
        nr = nrmt (is)
        Do ia = 1, natoms (is)
           ias = idxas (ia, is)
           Do l = 0, input%groundstate%lmaxapw
              Do io1 = 1, apword (l, is)
                do ir =1 ,nr
                    write(11,*)apwfr (ir, 1, io1, l, ias), apwfr (ir, 2, io1, l, ias)
                enddo !ir
              enddo !io
            enddo!l
        enddo
    enddo

    !!!LO functions
    Do is = 1, nspecies
        nr = nrmt (is)
        Do ia = 1, natoms (is)
           ias = idxas (ia, is)
           Do ilo = 1, nlorb (is)
                do ir =1 ,nr
                    write(11,*)lofr (ir, 1, ilo, ias),lofr (ir, 2, ilo, ias)
                enddo
            enddo
        enddo
    enddo

    close(11)
end Subroutine