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
module commvar
      implicit none
      save
      real*8,allocatable::wanfun_sol(:,:),wanfun_mol(:,:)
      real*8,allocatable::rorb(:),orb(:,:)
      real*8,allocatable::mdorb(:)
      real*8,allocatable::ovlp(:)
      real*8,allocatable::coeff(:)
      real*8,allocatable::evls(:)
      integer,allocatable::orb_indx(:)
      integer,allocatable::phase(:)
      integer::n_occ,noa,nx,ny,nz,nsp,nstates,nthreads
      integer::noa_LW,nocc_LW
      real*8::dV,dx,dy,dz
      character(9)::ch
end module commvar
