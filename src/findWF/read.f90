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
subroutine read_cnt
     use commvar
     implicit none
     character::rch,ccha
     integer::i,st,inpstat
     integer::zcha,nvele

     i=0
     open(111,file='cnt.ini',form='formatted',status='old',action='read')
       do
         read(111,*,iostat=inpstat) rch
         if(inpstat/=0) exit
         i = i+1
       enddo

       noa = i

       rewind(111)

       allocate(COA(noa,3),stat=st)
       if(st/=0) stop 'error: allcation of atomic coordinate matrix'
       allocate(Nel(noa),stat=st)
       if(st/=0) stop 'error: allocation of atomic effective electron array'
       allocate(Nvel(noa),stat=st)

       do i=1,noa
         read(111,*) ccha,COA(i,:)
         if(ccha=='Si') then
          cha = 14
          zcha = 4
         else if(ccha=='P') then
            cha = 15
            zcha = 5
         else if(ccha=='H ') then
            cha = 1
            zcha = 1
         else if(ccha=='C ') then
            cha = 6
            zcha = 4
         else if(ccha=='B ') then
            cha = 5
            zcha = 3
         elseif(ccha=='Li') then
            cha = 3
            zcha = 1
         else if(ccha=='N ') then
            cha = 7
            zcha = 5
         else if(ccha=='O ') then
            cha = 8
            zcha = 6
         else if(ccha=='F ') then
            cha = 9
            zcha = 7
         else if(ccha=='Zn') then
            cha = 30
            zcha = 12
         else if(ccha=='Ni') then
            cha = 28
            zcha = 10
         else if(ccha=='W ') then
            cha = 74
            zcha = 6
         else if(ccha=='Au') then
            cha = 79
            zcha = 11
         else if(ccha=='S ') then
            cha = 16
            zcha = 6
         else if(ccha=='He') then
            cha = 2
            zcha = 2
         else if(ccha=='Se') then
            cha = 34
            zcha = 6
         else if(ccha=='Ar') then
            cha = 18
            zcha = 8
         else if(ccha=='As') then
            cha = 33
            zcha = 5
         else if(ccha=='Cu') then
            cha = 29
            zcha = 11
         else if(ccha=='Cl') then
            cha = 17
            zcha = 7
         else if(ccha=='Br') then
            cha = 35
            zcha = 7
         else if(ccha=='I ') then
            cha = 53
            zcha = 7
         else if(ccha=='Be') then
            cha = 4
            zcha = 2
         else if(ccha=='Al') then
            cha = 13
            zcha = 3
         else if(ccha=='Rb') then
            cha = 37
            zcha = 1
         else if(ccha=='Na') then
            cha = 11
            zcha = 1
         else if(ccha=='K ') then
            cha = 19
            zcha = 1
         else if(ccha=='Ga') then
            cha = 31
            zcha = 3
         else if(ccha=='Ge') then
            cha = 32
            zcha = 4
         else if(ccha=='Kr') then
            cha = 36
            zcha = 8
         else if(ccha=='Ne') then
            cha = 10
            zcha = 8
         else if(ccha=='Ag') then
            cha = 47
            zcha = 11
         else if(ccha=='Ti') then
            cha = 22
            zcha = 4
         else if(ccha=='Xe') then
            cha = 54
            zcha = 8d0
         else if(ccha=='Mg') then
            cha = 12
            zcha = 2
         else
            write(6,*)' stopping, since we need to include also the case of atom ',ccha
            stop
         end if
         Nel(i) = cha
         Nvel(i) = zcha
     end do
     close(111)

     nvele = 0

     do i=1,noa
       nvele = nvele + Nvel(i)
     end do

     if(mod(nvele,2)==0) n_occ = nvele/2
       
end subroutine read_cnt

subroutine read_wannier
        use commvar
        implicit none
        integer::i,st
        
          open(111,file=trim(fname),form='unformatted',status='old',action='read')
          read(111) ch,nx
          read(111) ch,ny
          read(111) ch,nz
          read(111) ch,dx
          read(111) ch,dy
          read(111) ch,dz
          read(111) ch,nsp
          read(111) ch,nstates
          write(6,*) "nstates", nstates
          allocate(evls(nstates),stat=st)
          if(st/=0) stop 'error:allocation of eigenvalue array'
          read(111) !skip a line
          read(111) evls(:)
          !write(6,*) evls(:); call flush(6)
          read(111) !skip a line
          
          allocate(LO(nx*ny*nz,n_occ),stat=st)
          if(st/=0) stop 'error:allocation of canonical orbital matrix'
            
          do i=1,n_occ
            read(111) !skip a line
            read(111) LO(:,i)
          end do 

          close(111)

        dV = dx*dy*dz
end subroutine read_wannier

subroutine read_wannier_unocc
        use commvar
        implicit none
        integer::i,st

          open(111,file=trim(fname),form='unformatted',status='old',action='read')
          read(111) ch,nx
          read(111) ch,ny
          read(111) ch,nz
          read(111) ch,dx
          read(111) ch,dy
          read(111) ch,dz
          read(111) ch,nsp
          read(111) ch,nstates
          write(6,*) "nstates", nstates
          allocate(evls(nstates),stat=st)
          if(st/=0) stop 'error:allocation of eigenvalue array'
          read(111) !skip a line
          read(111) evls(:)
          !write(6,*) evls(:); call flush(6)
          read(111) !skip a line

          allocate(LO(nx*ny*nz,nstates),stat=st)
          if(st/=0) stop 'error:allocation of canonical orbital matrix'

          do i=1,nstates
            read(111) !skip a line
            read(111) LO(:,i)
          end do

          close(111)

        dV = dx*dy*dz
end subroutine read_wannier_unocc


subroutine read_lw
        use commvar
        implicit none
        integer::i,st,inpstat
        character(100) :: rch

        unocc = .false.
        nthreads = 10

        open(111,file='findWF',form='formatted',status='old',action='read')

        do
          read(111,*,iostat=inpstat) rch
          if(inpstat/=0) exit
          select case(rch)
        
          case('num_occ_lw')
            read(111,*,iostat=inpstat) nocc_LW
            if(inpstat/=0) stop "nocc_LW value is missing in findWF input"
          case('num_atom_lw')
            read(111,*,iostat=inpstat) noa_LW
            if(inpstat/=0) stop "noa_LW value is missing in findWF input"
            allocate(atom_label(noa_LW),stat=st)
            if(st/=0) stop "error: allocation of atom label array"
            read(111,*,iostat=inpstat) atom_label(:)
          case('num_of_threads')
            read(111,*,iostat=inpstat) nthreads
            if(inpstat/=0) stop "nthreads value is misssing in findWF input"
          case('unocc_state')
            read(111,*,iostat=inpstat) unocc
            if(inpstat/=0) stop "unocc value is missing in findWF input"
          case('num_unocc_lw')
            read(111,*,iostat=inpstat) nunocc_LW
            if(inpstat/=0) stop "nunocc_LW value is missing in findWF input"
          case default
            write(6,*) "unrecognized input variable", rch
            stop
          end select
        end do

        close(111)
end subroutine read_lw
