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
subroutine write_wf_trun
  use commvar
  implicit none
  integer::i
   open(222,file='wannier_trun.bin',form='unformatted',status='replace',action='write')
   write(222) 'nx       ', nx
   write(222) 'ny       ', ny
   write(222) 'nz       ', nz
   write(222) 'dx       ', dx
   write(222) 'dy       ', dy
   write(222) 'dz       ', dz
   write(222) 'nsp      ', nsp
   write(222) 'nstates  ', nocc_LW
   write(222) 'evls     '
   write(222) evls(:)
   write(222) 'orbitals '
   do i=1,nocc_LW
    write(222) i, 1
    write(222) LO(:,orb_indx(i))
   end do
  close(222)
end subroutine write_wf_trun

subroutine write_wf_unocc_trun
  use commvar
  implicit none
  integer::i
   open(222,file='wannier_unocc_trun.bin',form='unformatted',status='replace',action='write')
   write(222) 'nx       ', nx
   write(222) 'ny       ', ny
   write(222) 'nz       ', nz
   write(222) 'dx       ', dx
   write(222) 'dy       ', dy
   write(222) 'dz       ', dz
   write(222) 'nsp      ', nsp
   write(222) 'nstates  ', nocc_LW
   write(222) 'evls     '
   write(222) evls(:)
   write(222) 'orbitals '
   do i=1,nocc_LW
    write(222) i, 1
    write(222) LO(:,orb_indx(i))
   end do
  close(222)
end subroutine write_wf_unocc_trun

