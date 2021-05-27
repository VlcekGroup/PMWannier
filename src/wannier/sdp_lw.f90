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
subroutine steepest_descent_lw
      use OMP_LIB
      use commvar
      implicit none
      integer::i,j,k,l,m,count_iter=0
      real*8::sum_itmd=0d0,OBJ_0=0d0,OBJ_1=0d0,d_OBJ
      do
        if(count_iter.eq.0) then
          call CPU_time(start)
          !$OMP parallel
          !$OMP do schedule(dynamic)
            do j=1,n_occ
              do k=1,n_occ
                do i=1,noa_LW
                   Q_A0(i,j,k) = sum(CO(:,j)*AWF(:,atom_label(i))*CO(:,k))*dV
                end do
              end do
            end do
          !$OMP end do
          !$OMP end parallel
   
          !$OMP parallel
          !$OMP do
            do i=1,n_occ
             orb_con(i) = sum(Q_A0(:,i,i))
            end do
          !$OMP end do
          !$OMP end parallel
   
          do i=1,nocc_LW
            forb_con(i) = maxval(orb_con,n_occ)
            orb_indx(i) = maxloc(orb_con,n_occ)
            orb_con(orb_indx(i)) = 0d0
          end do
   
          !$OMP parallel shared(OBJ_0)
          !$OMP do reduction(+:OBJ_0) schedule(dynamic)
            do j=1,nocc_LW
               OBJ_0 = OBJ_0 + sum(Q_A0(:,orb_indx(j),orb_indx(j))**2d0)
            end do
          !$OMP end do
          !$OMP end parallel

          IF(allocated(AWF)) Deallocate(AWF)
          If(allocated(atom_label)) Deallocate(atom_label)

          call CPU_time(finish)
          write(6,*) 'step, OBJ, dOBJ, time, and delta_t:'; call flush(6)
          write(6,*) count_iter,OBJ_0,0,finish-start,delta_t; call flush(6)
end if

          call CPU_time(start)
          
          !$OMP parallel
          !$OMP do schedule(dynamic)
           do j=1,n_occ
             do i=1,n_occ
                A(i,j) = delta_t * sum(Q_A0(:,j,i)*(Q_A0(:,j,j)-Q_A0(:,i,i))-Q_A0(:,i,j)*(Q_A0(:,i,i)-Q_A0(:,j,j)))
             end do
           end do
          !$OMP end do
          !$OMP end parallel
     
          !get the unitary matrix U and U dagger
          call r8mat_expm1 (n_occ, A, U)
     
          !get the Q_A1 matrix
          !$OMP parallel firstprivate(sum_itmd)
          !$OMP do schedule(dynamic)
           do k=1,n_occ
             do j=1,n_occ
               do i=1,noa_LW
                 do l=1,n_occ
                   do m=1,n_occ
                     sum_itmd = sum_itmd + U(m,j)*Q_A0(i,m,l)*U(l,k)
                   end do
                 end do
                 Q_A1(i,j,k) = sum_itmd
                 sum_itmd = 0d0
               end do
             end do
           end do
          !$OMP end do
          !$OMP end parallel

          !$OMP parallel
          !$OMP do
            do i=1,n_occ
               orb_con(i) = sum(Q_A1(:,i,i))
            end do
          !$OMP end do
          !$OMP end parallel
   
          do i=1,nocc_LW
            forb_con(i) = maxval(orb_con,n_occ)
            orb_indx(i) = maxloc(orb_con,n_occ)
            orb_con(orb_indx(i)) = 0d0
          end do
   
         !calculate the objective function OBJ_1
         !$OMP parallel shared(OBJ_1)
         !$OMP do reduction(+:OBJ_1) schedule(dynamic)
           do j=1,nocc_LW
              OBJ_1 = OBJ_1 + sum(Q_A1(:,orb_indx(j),orb_indx(j))**2d0)
           end do
         !$OMP end do
         !$OMP end parallel
        
         !unitary transform
         !$OMP parallel
         !$OMP do
         do j=1,nx*ny*nz
          do i=1,n_occ
              LO(j,i) = sum(U(:,i)*CO(j,:))
          end do
         end do
        !$OMP end do
        !$OMP end parallel

        d_OBJ = OBJ_1-OBJ_0
        count_iter = count_iter + 1

        call CPU_time(finish)

      write(6,*) count_iter,OBJ_1,d_OBJ,finish-start,delta_t; call flush(6)
      if(d_OBJ.gt.0) then
        if(loc_wan) then
          if(mod(count_iter,sfw)==0) then
            write(6,*) '---------------------------'
            write(6,*) 'orbital_index, contribution'; call flush(6)
            do i=1,nocc_LW
              write(6,*) orb_indx(i),forb_con(i); call flush(6)
            end do
            write(6,*) '---------------------------'
          end if
        end if
      end if

      if(count_iter.eq.max_iter) exit

      if(d_OBJ.gt.0) then
        if(d_OBJ.lt.con_thre) then
          exit
        else
          Q_A0 = Q_A1
          CO = LO
          OBJ_0 = OBJ_1
          OBJ_1 = 0d0
          if(.not.unocc) then
            open(222,file='wannier_tmp.bin',form='unformatted',status='replace',action='write')
          else
            open(222,file='wannier_unocc_tmp.bin',form='unformatted',status='replace',action='write')
          endif
          write(222) 'nx       ', nx
          write(222) 'ny       ', ny
          write(222) 'nz       ', nz
          write(222) 'dx       ', dx
          write(222) 'dy       ', dy
          write(222) 'dz       ', dz
          write(222) 'nsp      ', nsp
          write(222) 'nstates  ', n_occ
          write(222) 'evls     '
          write(222) evls(:)
          write(222) 'orbitals '
          do i=1,n_occ
            write(222) i, 1
            write(222) LO(:,i)
          end do
          close(222)
        end if
      else
         delta_t = delta_t/corr_fac
         OBJ_1 = 0d0
      end if
 end do

 end subroutine steepest_descent_lw
