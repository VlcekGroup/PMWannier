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
      use OMP_LIB
      use commvar
      implicit none
      
      call write_header
      call cpu_time(start_all)
      call read_cnt
      call read_input

      if(unocc) then
        if(bulk) then
          call read_wf_unocc_sol
        else
          call read_wf_unocc_mol
        end if
      else 
        if(spec_state) then
          call read_wf_occ_ssi
        else
          call read_wf_occ
        end if
      endif

      if(sub_wan) then
         if(restart) then
           call mat_alloc2
         else
           call mat_alloc1
         end if
      else
         call mat_alloc1
      end if

      call OMP_set_num_threads(nthreads)
      call cpu_time(start) 
      call awf_prep
      call cpu_time(finish)
      write(6,*) 'max_iter and delta t:',max_iter,delta_t; call flush(6)
      write(6,*) 'time to prepare AWF matrix:', finish-start; call flush(6)
      if(sub_wan) write(6,*) 'subspace threshold:', sub_thre; call flush(6)

      call mat_dealloc1

      if(loc_wan) then
        if(sub_wan) then
          if(restart) then
            call steepest_descent_lw_sub_re
          else
            call steepest_descent_lw_sub
          end if
        else
          call steepest_descent_lw
        end if
      else
        call steepest_descent_full
      end if
      call cpu_time(finish_all)
      write(6,*) 'CPU time:',finish_all-start_all; call flush(6)

      call mat_dealloc2

      if(unocc) then
        call write_wf_unocc
        if(loc_wan) call write_wf_unocc_lw
      else
        if(sub_wan) then
          call write_wf_full_sub
        else
          call write_wf_full
        end if
        if(loc_wan) call write_wf_trun
      endif

      call mat_dealloc3
     
 end program wannier
