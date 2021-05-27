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
      real*8,allocatable::CO(:,:),LO(:,:),evls(:)
      real*8,allocatable::Q_A0(:,:,:),Q_A1(:,:,:)
      real*8,allocatable::AWF(:,:),GR(:,:),Nel(:)
      real*8,allocatable::COA(:,:),DOA(:,:)
      real*8,allocatable::A(:,:),U(:,:)
      real*8,allocatable::orb_con(:),forb_con(:)
      real*8,allocatable::lumo(:)
      integer,allocatable::atom_label(:),orb_indx(:)
      integer,allocatable::Nvel(:)
      integer,allocatable::state_indx(:)
      integer::n_occ,noa,nx,ny,nz,nsp,nstates,cha,nthreads,sfw
      integer::n_unocc, nunocc_rd
      integer::max_iter,noa_LW,nocc_LW,nunocc_LW
      integer::sub_fac,nocc_sv,nocc_sub
      integer::nocc_ssi
      real*8,parameter::gamma_a=0.944863,rcut=7.18096,pi=4*atan(1.0_8)
      real*8::sqrt_2pi=dsqrt(2d0*pi),delta_t,corr_fac
      real*8::dV,dx,dy,dz
      real*8::start_all,finish_all,start,finish
      real*8::ewdw,sub_thre
      real*8::con_thre
      character(9)::ch
      logical::loc_wan
      logical::restart,restart_old
      logical::unocc
      logical::bulk
      logical::sub_wan
      logical::spec_state
end module commvar
