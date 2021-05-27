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
subroutine mat_alloc1
        use commvar
        implicit none
        integer::st
        
        if(unocc) then
          n_occ=n_unocc
          if(loc_wan) nocc_LW = nunocc_LW
        end if

        if(spec_state) then
          n_occ = nocc_ssi
        end if

        allocate(A(n_occ,n_occ),stat=st)
        if(st/=0) stop 'error:allocation of A matrix'
        allocate(U(n_occ,n_occ),stat=st)
        if(st/=0) stop 'error: allocation of U matrix'
        allocate(LO(nx*ny*nz,n_occ),stat=st)
        if(st/=0) stop 'error: allocation of localized orbital matrix'
        allocate(DOA(nx*ny*nz,noa),stat=st)
        if(st/=0) stop 'error: allocation of density of atom matrix'
        allocate(GR(nx*ny*nz,3),stat=st)
        if(st/=0) stop 'error: allocation of grid matrix'
        allocate(AWF(nx*ny*nz,noa),stat=st)
        if(st/=0) stop 'error: allocation of atomic weight function matrix'

        If(loc_wan) then
          allocate(Q_A0(noa_LW,n_occ,n_occ),stat=st)
          if(st/=0) stop 'error:allocation in Q_A0 matrix'
          allocate(Q_A1(noa_LW,n_occ,n_occ),stat=st)
          if(st/=0) stop 'error:allocation in Q_A1 matrix'
          allocate(orb_indx(nocc_LW),stat=st)
          if(st/=0) stop 'error:allocation in orb_indx matrix'
          allocate(forb_con(nocc_LW),stat=st)
          if(st/=0) stop 'error:allocation in found orb_con matrix'
          allocate(orb_con(n_occ),stat=st)
          if(st/=0) stop 'error:allocation in orb_con matrix'
        else
          allocate(Q_A0(noa,n_occ,n_occ),stat=st)
          if(st/=0) stop 'error:allocation in Q_A0 matrix'
          allocate(Q_A1(noa,n_occ,n_occ),stat=st)
          if(st/=0) stop 'error:allocation in Q_A1 matrix'
        end if
end subroutine mat_alloc1

subroutine mat_alloc2
        use commvar
        implicit none
        integer::st
        
        nocc_sub = nocc_LW * sub_fac
   
        allocate(A(nocc_sub,nocc_sub),stat=st)
        if(st/=0) stop 'error:allocation of A matrix'
        allocate(U(nocc_sub,nocc_sub),stat=st)
        if(st/=0) stop 'error: allocation of U matrix'
        allocate(LO(nx*ny*nz,n_occ),stat=st)
        if(st/=0) stop 'error: allocation of localized orbital matrix'
        allocate(DOA(nx*ny*nz,noa),stat=st)
        if(st/=0) stop 'error: allocation of density of atom matrix'
        allocate(GR(nx*ny*nz,3),stat=st)
        if(st/=0) stop 'error: allocation of grid matrix'
        allocate(AWF(nx*ny*nz,noa),stat=st)
        if(st/=0) stop 'error: allocation of atomic weight function matrix'
        allocate(Q_A0(noa_LW,nocc_sub,nocc_sub),stat=st)
        if(st/=0) stop 'error:allocation in Q_A0 matrix'
        allocate(Q_A1(noa_LW,nocc_sub,nocc_sub),stat=st)
        if(st/=0) stop 'error:allocation in Q_A1 matrix'
        allocate(orb_indx(nocc_LW),stat=st)
        if(st/=0) stop 'error:allocation in orb_indx matrix'
        allocate(forb_con(nocc_LW),stat=st)
        if(st/=0) stop 'error:allocation in found orb_con matrix'
        allocate(orb_con(nocc_sub),stat=st)
        if(st/=0) stop 'error:allocation in orb_con matrix'

end subroutine mat_alloc2

subroutine mat_alloc3
        use commvar
        implicit none
        integer::st

        If(allocated(Q_A0)) Deallocate(Q_A0)
        If(allocated(Q_A1)) Deallocate(Q_A1)
        If(allocated(A)) Deallocate(A)
        If(allocated(U)) Deallocate(U)
        If(allocated(orb_con)) Deallocate(orb_con)

        allocate(A(nocc_sub,nocc_sub),stat=st)
        if(st/=0) stop 'error:allocation in A matrix in subspace file'
        allocate(Q_A0(noa_LW,nocc_sub,nocc_sub),stat=st)
        if(st/=0) stop 'error:allocation in Q_A0 matrix in subspace file'
        allocate(Q_A1(noa_LW,nocc_sub,nocc_sub),stat=st)
        if(st/=0) stop 'error:allocation in Q_A1 matrix in subspace file'
        allocate(U(nocc_sub,nocc_sub),stat=st)
        if(st/=0) stop 'error:allocation in U matrix in subspace file'
        allocate(orb_con(nocc_sub),stat=st)
        if(st/=0) stop 'error:allocation in orb_con matrix in subspace file'
end subroutine mat_alloc3

subroutine mat_dealloc1
      use commvar, only: COA,GR,Nel,DOA
      implicit none

      If(allocated(COA)) Deallocate(COA)
      If(allocated(GR)) Deallocate(GR)
      If(allocated(Nel)) Deallocate(Nel)
      If(allocated(DOA)) Deallocate(DOA)
end subroutine mat_dealloc1

subroutine mat_dealloc2
      use commvar, only: Q_A0,Q_A1,A,U,orb_con,forb_con,CO

      If(allocated(Q_A0)) Deallocate(Q_A0)
      If(allocated(Q_A1)) Deallocate(Q_A1)
      If(allocated(A)) Deallocate(A)
      If(allocated(U)) Deallocate(U)
      If(allocated(CO)) Deallocate(CO)
      If(allocated(orb_con)) Deallocate(orb_con)
      If(allocated(forb_con)) Deallocate(forb_con)
end subroutine mat_dealloc2

subroutine mat_dealloc3
      use commvar, only: evls,LO,orb_indx

      If(allocated(evls)) Deallocate(evls)
      If(allocated(LO)) Deallocate(LO)
      If(allocated(orb_indx)) Deallocate(orb_indx)
end subroutine mat_dealloc3
