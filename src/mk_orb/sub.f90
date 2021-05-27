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
subroutine coeff_prep
        use commvar
        implicit none
        integer :: st, i

        allocate(coeff(n_occ),stat=st)
        if(st/=0) stop "ERROR in allocating coeff matrix"

        do i=1,n_occ
          coeff(i) = (phase(i) * sum(wanfun_mol(:,i)*rorb(:))) / sum(abs(wanfun_mol(:,i))**2d0)
        end do
end subroutine coeff_prep

subroutine wanfun_reorder
        use commvar
        implicit none
        integer :: st, i, j
        real*8 :: maxi,mini

        allocate(ovlp(n_occ),stat=st)
        if(st/=0) stop "ERROR in allocating ovlp matrix"
        allocate(orb_indx(n_occ),stat=st)
        if(st/=0) stop "ERROR in allocating orb_indx matrix"
        allocate(mdorb(nx*ny*nz),stat=st)
        if(st/=0) stop "ERROR in allocating mdorb matrix"
        allocate(phase(n_occ),stat=st)
        if(st/=0) stop "ERROR in phase array!"

        phase = 1

        do i=1, n_occ
          do j=1, n_occ
            ovlp(j) = sum(wanfun_mol(:,i)*wanfun_sol(:,j))*dV
          enddo
          write(6,*) "Matching molecular wannier orbtial", i
          write(6,*) ovlp(:)
          if(abs(minval(ovlp,n_occ)).gt.abs(maxval(ovlp,n_occ))) then
            orb_indx(i) = minloc(ovlp,n_occ)
            phase(i) = -phase(i)
          else
            orb_indx(i) = maxloc(ovlp,n_occ)
          endif
          write(6,*) orb_indx(i)
        end do
        if(allocated(ovlp)) deallocate(ovlp)

end subroutine wanfun_reorder

subroutine mk_orb
        use commvar
        implicit none
        integer :: i, st

        mdorb = 0d0

        do i=1, n_occ
          mdorb(:) = mdorb(:) + coeff(i) * wanfun_sol(:,orb_indx(i))
        end do

        if(allocated(coeff)) deallocate(coeff)
end subroutine mk_orb

subroutine write_mdorb(iorb)
  use commvar
  implicit none
  character(30) :: fname
  integer :: iorb
   
  write(fname,'(a,i3.3,a)') 'mdorb',iorb,'.bin'

  open(222,file=trim(fname),form='unformatted',status='replace',action='write')
  write(222) 'nx       ', nx
  write(222) 'ny       ', ny
  write(222) 'nz       ', nz
  write(222) 'dx       ', dx
  write(222) 'dy       ', dy
  write(222) 'dz       ', dz
  write(222) 'nsp      ', nsp
  write(222) 'nstates  ', 1
  write(222) 'evls     '
  write(222) 0d0
  write(222) 'orbitals '

    write(222) 1, 1
    write(222) mdorb(:)

  close(222)
end subroutine write_mdorb

subroutine write_orbtxt(iorb)
  use commvar
  implicit none
  character(30) :: fname
  integer :: iorb

  write(fname,'(a,i3.3,a)') 'orb',iorb,'.txt'

  open(222,file=trim(fname),form='formatted',status='replace',action='write')
  write(222,*) 'nx       ', nx
  write(222,*) 'ny       ', ny
  write(222,*) 'nz       ', nz
  write(222,*) 'dx       ', dx
  write(222,*) 'dy       ', dy
  write(222,*) 'dz       ', dz
  write(222,*) 'nsp      ', nsp
  write(222,*) 'orb      ', evls(n_occ),1

  write(222,*) mdorb(:)

  close(222)
end subroutine write_orbtxt
