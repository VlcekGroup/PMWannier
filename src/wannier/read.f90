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
subroutine read_input
      use commvar
      implicit none
      character(100) :: rch
      integer :: inpstat,st,i

      restart = .false.
      loc_wan = .false.
      unocc   = .false.
      nthreads = 1
      max_iter = 50
      delta_t = 1
      corr_fac = 1.1
      n_unocc = 0
      nunocc_rd = 0
      ewdw    = 0d0
      bulk    = .false.
      sub_wan = .false.
      sub_fac = 3
      sub_thre = 1d0
      con_thre = 1E-8
      restart_old = .false.
      spec_state = .false.
      sfw     = 10

      open(111,file='input',form='formatted',status='old',action='read')
      do
        read(111,*,iostat=inpstat) rch
        if(inpstat/=0) exit
        select case(rch)
        case('num_of_atoms')
          read(111,*,iostat=inpstat) noa
          if(inpstat/=0) stop "noa value missing in input"
        case('num_of_occ')
          read(111,*,iostat=inpstat) n_occ
          if(inpstat/=0) stop "n_occ value missing in input"
        case('restart')
          read(111,*,iostat=inpstat) restart
          if(inpstat/=0) stop "restart value missing in input"
        case('restart_old')
          read(111,*,iostat=inpstat) restart_old
          if(inpstat/=0) stop "restart_old value missing in input"
        case('max_iter')
          read(111,*,iostat=inpstat) max_iter
          if(inpstat/=0) stop "max_iter value missing in input"
        case('delta_t')
          read(111,*,iostat=inpstat) delta_t
          if(inpstat/=0) stop "delta_t value missing in input"
        case('delta_t_correction')
          read(111,*,iostat=inpstat) corr_fac
          if(inpstat/=0) stop "corr_fac value missing in input"
        case('convergence_threshold')
          read(111,*,iostat=inpstat) con_thre
          if(inpstat/=0) stop "con_thre value missing in input"
        case('local_wannierization')
          read(111,*,iostat=inpstat) loc_wan
          if(inpstat/=0) stop "loc_wan value missing in input"
          if(loc_wan) then
            read(111,*,iostat=inpstat) noa_LW
            if(inpstat/=0) stop "noa_LW value missing in input"
            read(111,*,iostat=inpstat) nocc_LW
            if(inpstat/=0)  stop "nocc_LW value missing in input"
            allocate(atom_label(noa_LW),stat=st)
            if(st/=0) stop 'error: allocation of atom label array'
            read(111,*,iostat=inpstat) atom_label(:)
            if(inpstat/=0)  stop "atom label value missing in input"
          end if
        case('step_findwf')
          read(111,*,iostat=inpstat) sfw
          if(inpstat/=0) stop "sfw value missing in input"
        case('num_of_threads')
          read(111,*,iostat=inpstat) nthreads
          if(inpstat/=0) stop "nthreads value missing in input"
        case('unocc_state')
          read(111,*,iostat=inpstat) unocc
          if(inpstat/=0) stop "unocc value missing in input"
        case('num_unocc_lw')
          read(111,*,iostat=inpstat) nunocc_LW
          if(inpstat/=0)  stop "nunocc_LW value missing in input"
        case('subspace_wannierization')
          read(111,*,iostat=inpstat) sub_wan
          if(inpstat/=0) stop "sub_wan value missing in input"
        case('subspace_factor')
          read(111,*,iostat=inpstat) sub_fac
          if(inpstat/=0) stop "sub_fac value missing in input"
        case('subspace_threshold')
          read(111,*,iostat=inpstat) sub_thre
          if(inpstat/=0) stop "sub_thre value missing in input"
        case('num_of_unocc')
          read(111,*,iostat=inpstat) nunocc_rd
          if(inpstat/=0) stop "nunocc_rd value missing in input"
        case('bulk')
          read(111,*,iostat=inpstat) bulk
          if(inpstat/=0) stop "bulk value missing in input"
        case('ewdw')
          read(111,*,iostat=inpstat) ewdw
          if(inpstat/=0) stop "ewdw value missing in input"
        case('specific_state')
          read(111,*,iostat=inpstat) spec_state
          if(inpstat/=0) stop "specstate value missing in input"
          if(spec_state) then 
            read(111,*,iostat=inpstat) nocc_ssi
            if(inpstat/=0) stop "nocc_ssi value missing in input"
            allocate(state_indx(nocc_ssi),stat=st)
            if(st/=0) stop "error: allcation of spec_state_indx"
            read(111,*) state_indx(:)
          end if
        case default
          write(6,*) "unrecognized input variable", rch
          stop
        end select
      end do

      close(111)

 end subroutine read_input

 subroutine read_cnt
     use commvar
     implicit none
     character::rch
     character(9)::ccha
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
         else if(ccha=='H') then
            cha = 1
            zcha = 1
         else if(ccha=='C') then
            cha = 6
            zcha = 4
         else if(ccha=='B') then
            cha = 5
            zcha = 3
         elseif(ccha=='Li') then
            cha = 3
            zcha = 1
         else if(ccha=='N') then
            cha = 7
            zcha = 5
         else if(ccha=='O') then
            cha = 8
            zcha = 6
         else if(ccha=='F') then
            cha = 9
            zcha = 7
         else if(ccha=='Zn') then
            cha = 30
            zcha = 12
         else if(ccha=='Ni') then
            cha = 28
            zcha = 10
         else if(ccha=='W') then
            cha = 74
            zcha = 6
         else if(ccha=='Au') then
            cha = 79
            zcha = 11
         else if(ccha=='S') then
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
         else if(ccha=='I') then
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
         else if(ccha=='K') then
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

     n_occ = nvele/2

end subroutine read_cnt

subroutine read_wf_occ
        use commvar
        implicit none
        integer::i,st
        
        If(restart) then
          open(111,file='wannier.bin',form='unformatted',status='old',action='read')
        else
          open(111,file='wf.bin',form='unformatted',status='old',action='read')
        end if
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
          
          allocate(CO(nx*ny*nz,n_occ),stat=st)
          if(st/=0) stop 'error:allocation of canonical orbital matrix'
            
          do i=1,n_occ
            read(111) !skip a line
            read(111) CO(:,i)
          end do 

          close(111)

        dV = dx*dy*dz
end subroutine read_wf_occ

subroutine read_wf_unocc_mol
        use commvar
        implicit none
        integer::i,st
        integer::n_unoccbd=0

        If(restart) then
          open(111,file='wannier_unocc.bin',form='unformatted',status='old',action='read')
            read(111) ch,nx
            read(111) ch,ny
            read(111) ch,nz
            read(111) ch,dx
            read(111) ch,dy
            read(111) ch,dz
            read(111) ch,nsp
            read(111) ch,n_unocc
            allocate(evls(n_unocc),stat=st)
            if(st/=0) stop 'error:allocation of eigenvalue array'
            read(111) !skip 'evls'
            read(111) evls(:)
            read(111) !skip orbital
            
            allocate(CO(nx*ny*nz,n_unocc),stat=st)
            if(st/=0) stop 'error:allocation of canonical orbital matrix'
             
            do i=1,n_unocc
              read(111)
              read(111) CO(:,i)
            enddo

          close(111)

          dV=dx*dy*dz

          open(111,file='wf.bin',form='unformatted',status='old',action='read')
          read(111) 
          read(111) 
          read(111) 
          read(111) 
          read(111) 
          read(111) 
          read(111) 
          read(111) 
          read(111) 
          read(111) 
          read(111)

          allocate(lumo(nx*ny*nz),stat=st)
          if(st/=0) stop 'error:allocation of lumo vector'

          do i=1,n_occ
            read(111)
            read(111)
          end do

          read(111)
          read(111) lumo(:)

          close(111)

        else

          open(111,file='wf.bin',form='unformatted',status='old',action='read')
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

            do i=1,n_occ
              read(111)
              read(111)
            end do

            do i=n_occ+1, nstates
              if(evls(i).lt.0d0) n_unoccbd=n_unoccbd+1
              if(evls(i).gt.0d0) exit
            end do

            if(n_unoccbd==0) stop "No empty state is bound"

            write(6,*) "number of bound virtual states:", n_unoccbd; call flush(6)
            write(6,*) "number of virtual states specified:", nunocc_rd; call flush(6)
            
            if(nunocc_rd/=0) then
              n_unocc = nunocc_rd
            else
              n_unocc = n_unoccbd
            end if

            allocate(CO(nx*ny*nz,n_unocc),stat=st)
            if(st/=0) stop 'error:allocation of canonical orbital matrix'
            allocate(lumo(nx*ny*nz),stat=st)
            if(st/=0) stop 'error:allocation of lumo vector'

            do i=1,n_unocc
              read(111)
              read(111) CO(:,i)
            enddo

            lumo = CO(:,1)

          close(111)

          dV = dx*dy*dz

        end if
end subroutine read_wf_unocc_mol

subroutine read_wf_unocc_sol
        use commvar
        implicit none
        integer::i,st
        integer::n_unoccrd=0
        real*8 :: emax

        If(restart) then
          open(111,file='wannier_unocc.bin',form='unformatted',status='old',action='read')
            read(111) ch,nx
            read(111) ch,ny
            read(111) ch,nz
            read(111) ch,dx
            read(111) ch,dy
            read(111) ch,dz
            read(111) ch,nsp
            read(111) ch,n_unocc
            allocate(evls(n_unocc),stat=st)
            if(st/=0) stop 'error:allocation of eigenvalue array'
            read(111) !skip 'evls'
            read(111) evls(:)
            read(111) !skip orbital

            allocate(CO(nx*ny*nz,n_unocc),stat=st)
            if(st/=0) stop 'error:allocation of canonical orbital matrix'

            do i=1,n_unocc
              read(111)
              read(111) CO(:,i)
            enddo

          close(111)

          dV=dx*dy*dz

        else

          open(111,file='wf.bin',form='unformatted',status='old',action='read')
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

            do i=1,n_occ
              read(111)
              read(111)
            end do

            emax = evls(n_occ) + ewdw

            do i=n_occ+1, nstates
              if(evls(i).lt.emax) then
                n_unoccrd=n_unoccrd+1
              else
                exit
              endif
            end do

            if(n_unoccrd==0) stop "No virtual state is available"

            write(6,*) "number of virtual states used:", n_unoccrd; call flush(6)
            write(6,*) "number of virtual wannier functions of interest:", nunocc_LW; call flush(6)

            n_unocc = n_unoccrd

            allocate(CO(nx*ny*nz,n_unocc),stat=st)
            if(st/=0) stop 'error:allocation of canonical orbital matrix'

            do i=1,n_unocc
              read(111)
              read(111) CO(:,i)
            enddo

          close(111)

          dV = dx*dy*dz

        end if
end subroutine read_wf_unocc_sol

subroutine read_wf_occ_ssi
        use commvar
        implicit none
        integer::i,st
        real*8,allocatable :: CO_sv(:,:)

        If(restart) then
          open(111,file='wannier.bin',form='unformatted',status='old',action='read')
        else
          open(111,file='wf.bin',form='unformatted',status='old',action='read')
        end if
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

          allocate(CO_sv(nx*ny*nz,n_occ),stat=st)
          if(st/=0) stop 'error:allocation of canonical orbital matrix'
          allocate(CO(nx*ny*nz,nocc_ssi),stat=st)
          if(st/=0) stop 'error:allocation of canonical orbital matrix'

          do i=1,n_occ
            read(111) !skip a line
            read(111) CO_sv(:,i)
          end do

          close(111)

          do i=1,nocc_ssi
            CO(:,i) = CO_sv(:,state_indx(i))
          end do

          deallocate(CO_sv)
          deallocate(state_indx)

        dV = dx*dy*dz
end subroutine read_wf_occ_ssi
