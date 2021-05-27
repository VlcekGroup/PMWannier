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
subroutine write_header
  implicit none
     write(6,*)' ************************************************************** '
     write(6,*)
     write(6,*)' Pipek-Mezey Wannier Code, v1.0 (Apr/2021) '
     write(6,*)
     write(6,*)' ************************************************************** '
     write(6,*)
     call flush(6)
end subroutine write_header

