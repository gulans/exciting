


Subroutine storecore
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
    use modmpi, only: mpiglobal
    integer :: is,ia,ias,ir

   if(mpiglobal%is_root) then
    open (11, file = "STATE_CORE.OUT", status = 'replace')


    do is = 1, nspecies
        Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ist = 1, spnst (is)
                If (spcore(ist, is)) Then
                    write(11,*)evalcr (ist, ias)
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
                        write(11,*)rwfcr(ir, 1, ist, ias)
                    enddo
                endif
            Enddo
        enddo
    enddo
  
    close(11)
   endif
end Subroutine
