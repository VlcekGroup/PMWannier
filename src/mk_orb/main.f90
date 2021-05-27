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
program wannier
      use commvar
      implicit none
      integer::iorb
      
      call read_wanfun_sol
      call read_wanfun_mol
      call read_wf_mol
      call wanfun_reorder
      do iorb=1,n_occ
        write(6,*) '------------------------'; call flush(6)
        write(6,*) 'Making Kohn-Sham Orbital',iorb; call flush(6)
        rorb(:) = orb(:,iorb)
        call coeff_prep
        call mk_orb
        call write_mdorb(iorb)
        call write_orbtxt(iorb)
        write(6,*) '------------------------'; call flush(6)
      end do

 end program wannier


