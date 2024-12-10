


Subroutine loadcore
    ! !USES:
    use modinput, only: input
    use mod_atoms, only: natoms, idxas, spvr, spnr, nspecies, &
    & spnst, spcore, spn, spr, spk, spl, spocc, natmmax, &
    & spnrmax
    use mod_symmetry, only: eqatoms
    use mod_muffin_tin, only: nrmt
    use mod_potential_and_density, only: veffmt
    use constants, only: y00, fourpi
    use mod_corestate, only: rhocr, rwfcr, evalcr, engy_exnl_core
    implicit none
    integer :: is,ia,ias,ir,ist
    real(8) :: t1
    ! do ist = 1, spnst(is)
    !     t1 = spocc (ist, is)
    !         If (spcore(ist, is)) Then
    !            Do ir = 1, spnr (is)
    !                rhocr (ir, ias) = rhocr (ir, ias) + t1 * &
    !              & rwfcr(ir, 1, ist, ias)**2
    !            Enddo
    !          end IF
    !     end do


    open (11, file = "STATE_CORE.OUT", status = 'old')
  
    do is = 1, nspecies
        Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist = 1, spnst (is)
                If (spcore(ist, is)) Then
                    read(11,*)evalcr (ist, ias)
                endif
            Enddo
        enddo
    enddo

    do is = 1, nspecies
        Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist = 1, spnst (is)
                If (spcore(ist, is)) Then
                    do ir = 1, spnr(is)
                        read(11,*)rwfcr(ir, 1, ist, ias)
                    enddo
                endif
            Enddo
        enddo
    enddo

    close(11)

    rhocr = 0d0
    do is = 1, nspecies
        Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            do ist=1, spnst (is)
                If (spcore(ist, is)) Then
                    t1 = spocc (ist, is)
                    Do ir = 1, spnr (is)
                        rhocr (ir, ias) = rhocr (ir, ias) + t1 * rwfcr(ir, 1, ist, ias)**2
                    End Do
                endif
            enddo
        Enddo
    enddo

    do is = 1, nspecies
        Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ir = 1, spnr (is)
                rhocr (ir, ias) = rhocr (ir, ias) / (fourpi*spr(ir,is)**2)
            End Do
        Enddo
    enddo
    
end Subroutine