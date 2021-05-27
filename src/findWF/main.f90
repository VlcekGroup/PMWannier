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
program find_wanfun
      use OMP_LIB
      use commvar
      implicit none

      call getarg(1,fname)
      
      write(6,*) 'reading cnt.in';call flush(6)
      call read_cnt

      write(6,*) 'reading findWF';call flush(6)
      call read_lw

      write(6,*) 'reading wannier.bin';call flush(6)
      if(unocc) then
        call read_wannier_unocc
      else
        call read_wannier
      end if

      write(6,*) 'allocating matrices';call flush(6)
      call mat_alloc
      call OMP_set_num_threads(nthreads)

      write(6,*) 'preparing atomic weight function';call flush(6)
      call awf_prep
      call mat_dealloc1

      write(6,*) 'Searching for wannier function of interest';call flush(6)
      call find_wf
      call mat_dealloc2

      write(6,*) 'wrting wannier_trun.bin';call flush(6)
      if(unocc) then
        call write_wf_unocc_trun
      else
        call write_wf_trun
      end if

      call mat_dealloc3
     
end program find_wanfun


