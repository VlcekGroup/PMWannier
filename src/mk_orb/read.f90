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
subroutine read_wanfun_sol
        use commvar
        implicit none
        integer::i,st
        
          open(111,file='wannier_trun.bin',form='unformatted',status='old',action='read')
          read(111) ch,nx
          read(111) ch,ny
          read(111) ch,nz
          read(111) ch,dx
          read(111) ch,dy
          read(111) ch,dz
          read(111) ch,nsp
          read(111) ch,nocc_LW
          write(6,*) "nocc_LW", nocc_LW
          read(111) !skip a line
          read(111) !skip a line
          read(111) !skip a line
          
          allocate(wanfun_sol(nx*ny*nz,nocc_LW),stat=st)
          if(st/=0) stop 'error:allocation of wanfun_sol matrix'
            
          do i=1,nocc_LW
            read(111) !skip a line
            read(111) wanfun_sol(:,i)
          end do 

          close(111)

        dV = dx*dy*dz
end subroutine read_wanfun_sol

subroutine read_wanfun_mol
        use commvar
        implicit none
        integer::i,st
        
          open(111,file='wannier_mol.bin',form='unformatted',status='old',action='read')
          read(111) ch,nx
          read(111) ch,ny
          read(111) ch,nz
          read(111) ch,dx
          read(111) ch,dy
          read(111) ch,dz
          read(111) ch,nsp
          read(111) ch,n_occ
          if(nocc_LW/=n_occ) stop "ERROR: number of occupied states mismatch"
          read(111) !skip a line
          read(111) !skip a line
          read(111) !skip a line

          allocate(wanfun_mol(nx*ny*nz,n_occ),stat=st)
          if(st/=0) stop 'error:allocation of wanfun_mol matrix'

          do i=1,n_occ
            read(111) !skip a line
            read(111) wanfun_mol(:,i)
          end do

          close(111)

end subroutine read_wanfun_mol

subroutine read_wf_mol
        use commvar
        implicit none
        integer::i,st

            open(111,file='wf_mol.bin',form='unformatted',status='old',action='read')
            read(111) ch,nx
            read(111) ch,ny
            read(111) ch,nz
            read(111) ch,dx
            read(111) ch,dy
            read(111) ch,dz
            read(111) ch,nsp
            read(111) ch,nstates
            allocate(evls(nstates),stat=st)
            if(st/=0) stop 'error:allocation of eigenvalue array'
            read(111) !skip 'evls'
            read(111) evls(:)
            read(111) !skip orbital
            
            allocate(rorb(nx*ny*nz),stat=st)
            if(st/=0) stop 'error:allocation of rorb matrix!'
            allocate(orb(nx*ny*nz,n_occ),stat=st)
            if(st/=0) stop 'error:allocation of canonical orbital matrix'
             
            do i=1,n_occ
              read(111)
              read(111) orb(:,i)
            end do

          close(111)

end subroutine read_wf_mol
