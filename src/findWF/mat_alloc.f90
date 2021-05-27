!                     
!                     
!    The routine(s) in this file are a part of the  
!                     PMwannier                 
!    suite, developed 2020-2021, and copyrighted    
!    to the authors: Guorong Weng and Vojtech Vlcek,
!    at the University of California, Santa Barbara.
!                                                   
!                                                   
!   If you use or modify any part of this routine   
!   the header should be kept and unmodified.          
!                                                   
!                                                   
! 
subroutine mat_alloc
        use commvar
        implicit none
        integer::st
        
        if(unocc) then
          nocc_LW = nunocc_LW
          n_occ = nstates
        end if

        allocate(DOA(nx*ny*nz,noa),stat=st)
        if(st/=0) stop 'error: allocation of density of atom matrix'
        allocate(GR(nx*ny*nz,3),stat=st)
        if(st/=0) stop 'error: allocation of grid matrix'
        allocate(AWF(nx*ny*nz,noa),stat=st)
        if(st/=0) stop 'error: allocation of atomic weight function matrix'

        allocate(Q_A0(noa_LW,n_occ),stat=st)
        if(st/=0) stop 'error:allocation in Q_A0 matrix'
        allocate(orb_indx(nocc_LW),stat=st)
        if(st/=0) stop 'error:allocation in orb_indx matrix'
        allocate(forb_con(nocc_LW),stat=st)
        if(st/=0) stop 'error:allocation in found orb_con matrix'
        allocate(orb_con(n_occ),stat=st)
        if(st/=0) stop 'error:allocation in orb_con matrix'
end subroutine mat_alloc

subroutine mat_dealloc1
      use commvar
      implicit none

      If(allocated(COA)) Deallocate(COA)
      If(allocated(GR)) Deallocate(GR)
      If(allocated(Nel)) Deallocate(Nel)
      If(allocated(DOA)) Deallocate(DOA)
      If(allocated(Nvel)) Deallocate(Nvel)
end subroutine mat_dealloc1

subroutine mat_dealloc2
      use commvar

      If(allocated(Q_A0)) Deallocate(Q_A0)
      If(allocated(orb_con)) Deallocate(orb_con)
      If(allocated(forb_con)) Deallocate(forb_con)
end subroutine mat_dealloc2

subroutine mat_dealloc3
      use commvar, only: evls,LO,orb_indx

      If(allocated(evls)) Deallocate(evls)
      If(allocated(LO)) Deallocate(LO)
      If(allocated(orb_indx)) Deallocate(orb_indx)
end subroutine mat_dealloc3
