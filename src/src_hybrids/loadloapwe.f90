Subroutine loadapwloe
use modmain


integer :: is,ia,ias,io1,ilo,l,i1,i2

real(8) :: e1


write(*,*)"ielasam enerģijas no faila"
open (2, file = 'base_apw_hf.dat', status = 'old')
 Do is = 1, nspecies
         Do ia = 1, natoms (is)
read(2,*)
read(2,*)!is ia
read(2,*)
            ias = idxas (ia, is)
            Do l = 0, input%groundstate%lmaxapw
               Do io1 = 1, apword (l, is)
read(2,*)i1,i2,e1
write(*,*)i1,i2,e1
        apwe(io1, l, ias)=e1

               enddo!io1
            enddo!l
         enddo!ia
       enddo!is
close(2)

open (2, file = 'base_lo_hf.dat', status = 'old')
Do is = 1, nspecies
  do ia = 1, natoms (is)
    read(2,*)!"is ia"
    read(2,*)!is, ia
    read(2,*)!"l order energy"
    ias = idxas (ia, is)
    Do ilo = 1, nlorb (is)
      l = lorbl (ilo, is)
      Do io1 = 1, lorbord (ilo, is)
        read(2,*)i1,i2,e1
        write(*,*)i1,i2,e1
        lorbe(io1, ilo, ias)=e1
      enddo
    enddo
  enddo
enddo
end subroutine
