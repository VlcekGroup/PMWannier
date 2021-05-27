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
subroutine find_wf
      use OMP_LIB
      use commvar
      implicit none
      integer::i,j,k,l,m
      real*8::sum_itmd=0d0,OBJ_0=0d0

          !$OMP parallel
          !$OMP do schedule(dynamic)
            do j=1,n_occ
                do i=1,noa_LW
                   Q_A0(i,j) = sum(LO(:,j)*AWF(:,atom_label(i))*LO(:,j))*dV
                   !write(6,*) i,j,k
                end do
            end do
          !$OMP end do
          !$OMP end parallel

          write(6,*) "done Q_A0"; call flush(6)
   
          !$OMP parallel
          !$OMP do
            do i=1,n_occ
             orb_con(i) = sum(Q_A0(:,i))
            end do
          !$OMP end do
          !$OMP end parallel
   
          write(6,*) "done orb_con(i)"; call flush(6)

          do i=1,nocc_LW
            forb_con(i) = maxval(orb_con,n_occ)
            orb_indx(i) = maxloc(orb_con,n_occ)
            orb_con(orb_indx(i)) = 0d0
          end do

          write(6,*) "done orb_indx and forb_con"; call flush(6)
   
          !$OMP parallel shared(OBJ_0)
          !$OMP do reduction(+:OBJ_0) schedule(dynamic)
            do j=1,nocc_LW
               OBJ_0 = OBJ_0 + sum(Q_A0(:,orb_indx(j))**2d0)
            end do
          !$OMP end do
          !$OMP end parallel

          write(6,*) "done OBJ_0"; call flush(6)

          IF(allocated(AWF)) Deallocate(AWF)
          If(allocated(atom_label)) Deallocate(atom_label)

          write(6,*) 'OBJ'; call flush(6)
          write(6,*) OBJ_0; call flush(6)
   
          write(6,*) '---------------------------'
          write(6,*) 'orbital_index, contribution'; call flush(6)
            do i=1,nocc_LW
              write(6,*) orb_indx(i),forb_con(i); call flush(6)
            end do
          write(6,*) '---------------------------'

 end subroutine find_wf
