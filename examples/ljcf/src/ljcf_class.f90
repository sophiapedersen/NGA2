!> Definition for a ljcf atomization class
module ljcf_class
   use precision,         only: WP,SP
   use config_class,      only: config
   use iterator_class,    only: iterator
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use hypre_str_class,   only: hypre_str
   !use ddadi_class,       only: ddadi
   use fft3d_class,       only: fft3d
   use stracker_class,    only: stracker
   use vfs_class,         only: vfs
   use tpns_class,        only: tpns
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use timer_class,       only: timer
   use pardata_class,     only: pardata
   use cclabel_class,     only: cclabel
   use irl_fortran_interface
   use interface_smoothing, only: marching_cubes2

   implicit none
   private
   
   public :: ljcf

   integer :: ierr
   
   
   ! For SILO_10.2, use f9x
   include "silo_f9x.inc"

   !> ljcf object
   type :: ljcf
      
      !> Config
      type(config) :: cfg
      
      !> Flow solver
      type(vfs)        ,public :: vf    !< Volume fraction solver
      type(tpns)        :: fs    !< Two-phase flow solver
      type(hypre_str)   :: ps    !< Structured Hypre linear solver for pressure
      !type(ddadi)       :: vs    !< DDADI solver for velocity
      type(timetracker) :: time  !< Time info
      type(cclabel)     :: ccl   !< CCLabel for local Weber number calculation

      !> Include structure tracker
      type(stracker) :: strack

      !> Ensight postprocessing
      type(surfmesh) :: smesh    !< Surface mesh for interface
      type(ensight)  :: ens_out  !< Ensight output for flow variables
      type(event)    :: ens_evt, drop_evt  !< Event trigger for Ensight output
      
      !> Simulation monitor file
      type(monitor) :: mfile    !< General simulation monitoring
      type(monitor) :: cflfile  !< CFL monitoring
      type(monitor) :: ljcf_file     !< LJCF simulation monitoring
      type(monitor) :: hitfile  !< added these lines from stracker
      type(monitor) :: cvgfile   
      
      !> Work arrays
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW      !< Residuals
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi            !< Cell-centered velocities
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU
      
      !> Iterator for VOF removal
      type(iterator) :: vof_removal_layer  !< Edge of domain where we actively remove VOF
      real(WP) :: vof_removed              !< Integral of VOF removed
      integer  :: nlayer=4                 !< Size of buffer layer for VOF removal
      
      !> Timing info
      type(monitor) :: timefile !< Timing monitoring
      type(timer)   :: tstep    !< Timer for step
      type(timer)   :: tvel     !< Timer for velocity
      type(timer)   :: tpres    !< Timer for pressure
      type(timer)   :: tvof     !< Timer for VOF
      
      !> Provide a pardata and an event tracker for saving restarts
      type(event)   :: save_evt
      type(pardata) :: df
      logical :: restarted

      !> Problem definition
      real(WP) :: djet, Vjet
      real(WP), dimension(:), allocatable :: xjet
      integer :: relax_model, nwall
      real(WP) :: gravity, liqVol, liqVolInjected, InjectionVelocity
      
   contains
      procedure :: init     !< Initialize nozzle simulation
      procedure :: step     !< Advance nozzle simulation by one time step
      procedure :: final    !< Finalize nozzle simulation
   end type ljcf
   
   !> Type for structure stats
   type :: struct_stats
      real(WP) :: vol
      real(WP) :: x_cg,y_cg,z_cg
      real(WP) :: u_avg,v_avg,w_avg
      real(WP), dimension(3,3) :: Imom
      real(WP), dimension(3) :: lengths
      real(WP), dimension(3,3) :: axes
   end type struct_stats
   
contains

   !> Perform droplet analysis
   subroutine analyse_drops(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: Pi
      use string,    only: str_medium
      use filesys,   only: makedir,isdir
      class(ljcf), intent(inout) :: this
      character(len=str_medium) :: filename,timestamp
      real(WP), dimension(:), allocatable :: dvol
      integer :: iunit,n,m,ierr
      ! Allocate droplet volume array
      allocate(dvol(1:this%strack%nstruct)); dvol=0.0_WP
      ! Loop over individual structures
      do n=1,this%strack%nstruct
         ! Loop over cells in structure and accumulate volume
         do m=1,this%strack%struct(n)%n_
            dvol(n)=dvol(n)+this%cfg%vol(this%strack%struct(n)%map(1,m),this%strack%struct(n)%map(2,m),this%strack%struct(n)%map(3,m))*&
            &                 this%vf%VF(this%strack%struct(n)%map(1,m),this%strack%struct(n)%map(2,m),this%strack%struct(n)%map(3,m))
         end do
      end do
      ! Reduce volume data
      call MPI_ALLREDUCE(MPI_IN_PLACE,dvol,this%strack%nstruct,MPI_REAL_WP,MPI_SUM,this%vf%cfg%comm,ierr)
      ! Only root process outputs to a file
      if (this%cfg%amRoot) then
         if (.not.isdir('diameter')) call makedir('diameter')
         filename='diameter_'; write(timestamp,'(es12.5)') this%time%t
         open(newunit=iunit,file='diameter/'//trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         do n=1,this%strack%nstruct
            ! Output list of diameters
            write(iunit,'(999999(es12.5,x))') (6.0_WP*dvol(n)/Pi)**(1.0_WP/3.0_WP)
         end do
         close(iunit)
      end if
   end subroutine analyse_drops

      !> Perform merge/split analysis
   subroutine analyze_merge_split(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
      use parallel, only: MPI_REAL_WP
      implicit none
      class(ljcf), intent(inout) :: this
      integer :: iunit
      logical :: file_exists
      type(struct_stats) :: stats

      ! Open the file - Created in simulation_init
      if (this%cfg%amRoot) open(iunit,file="merge_split.csv",form="formatted",status="old",position="append",action="write")
   
      analyze_merges: block
         integer :: n,nn

         ! Traverse merge events
         do n=1,this%strack%nmerge_master

            call compute_struct_stats(this%strack%merge_master(n)%newid,stats)
            if (this%cfg%amRoot) then 
               ! Write merge data to file
               this%strack%eventcount = this%strack%eventcount+1
               write(iunit,"(I0)",      advance="no")  this%strack%eventcount
               write(iunit,"(A)",       advance="no")  ', Merge,'
               do nn=1,this%strack%merge_master(n)%noldid
                  write(iunit,"(I0)",   advance="no")  this%strack%merge_master(n)%oldids(nn)
                  write(iunit,"(A)",    advance="no")  ';'
               end do
               write(iunit,"(A)",       advance="no")   ','
               write(iunit,"(I0)",      advance="no")  this%strack%merge_master(n)%newid
               write(iunit,"(A)",       advance="no")  ','
               write(iunit,"(ES12.5 )", advance="no")  this%time%t
               write(iunit,"(A)",       advance="no")   ','
               write(iunit,"(ES22.16)", advance="yes") stats%vol
               ! write SILO files for merge events
            end if 
            if (this%strack%merge_master(n)%newid.ne.1) then ! not liquid core
               !print*, "Dumping GTDA for merge event, newid=", this%strack%merge_master(n)%newid
               call dump_gtda(this,this%strack%merge_master(n)%newid,1)
            end if
            !call silo_lists_update(1)
         end do
         !call timing_stop('merge_split')
      end block analyze_merges

      analyze_splits: block 
      integer :: n,nn

         ! Traverse split events
         do n=1,this%strack%nsplit_master
            ! Write stats for each new structure after split
            do nn=1,this%strack%split_master(n)%nnewid
               call compute_struct_stats(this%strack%split_master(n)%newids(nn),stats)
               if (this%cfg%amRoot) then 

                  ! Write merge data to file
                  this%strack%eventcount = this%strack%eventcount+1
                  write(iunit,"(I0)",      advance="no")  this%strack%eventcount
                  write(iunit,"(A)",       advance="no")  ', Split,'
                  write(iunit,"(I0)",      advance="no")  this%strack%split_master(n)%oldid
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(I0)",      advance="no")  this%strack%split_master(n)%newids(nn)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES12.5 )", advance="no")  this%time%t
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES22.16)", advance="no")  stats%vol
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%x_cg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%y_cg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%z_cg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%u_avg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%v_avg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%w_avg
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%lengths(1)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="no")  stats%lengths(2)
                  write(iunit,"(A)",       advance="no")  ','
                  write(iunit,"(ES20.12)", advance="yes")  stats%lengths(3)
               end if 
               if (this%strack%split_master(n)%newids(nn).ne.1) then ! not liquid core
                  !print*, "Dumping GTDA for split event, newid=", this%strack%split_master(n)%newids(nn)
                  call dump_gtda(this,this%strack%split_master(n)%newids(nn),2)
                  !print*, "Done dumping GTDA for split event, newid=", this%strack%split_master(n)%newids(nn)
               end if
               !call silo_lists_update(2)
            end do
         end do
      
      end block analyze_splits

      if (this%cfg%amRoot) close(iunit)
   
      contains 

      subroutine compute_struct_stats(id,stats)
         implicit none 
         integer, intent(in) :: id
         type(struct_stats), intent(inout) :: stats

         integer :: n,m
         integer :: lwork,info,ierr
         integer :: ii,jj,kk
         integer  :: per_x,per_y,per_z
         real(WP) :: vol_struct
         real(WP) :: x_vol,y_vol,z_vol
         real(WP) :: u_vol,v_vol,w_vol
         real(WP), dimension(3,3) :: Imom
         real(WP) :: xtmp,ytmp,ztmp
         real(WP), dimension(3) :: lengths
         real(WP), dimension(3,3) :: axes
         
         ! Eigenvalues/eigenvectors
         real(WP), dimension(3,3) :: A
         real(WP), dimension(3) :: d
         integer , parameter :: order = 3
         real(WP), dimension(:), allocatable :: work
         real(WP), dimension(1)   :: lwork_query

         
         ! Query optimal work array size
         call dsyev('V','U',order,A,order,d,lwork_query,-1,info); lwork=int(lwork_query(1)); allocate(work(lwork))

         ! Initialize values
         vol_struct    = 0.0_WP ! Structure volume
         x_vol = 0.0_WP; y_vol = 0.0_WP; z_vol = 0.0_WP ! Center of gravity
         u_vol = 0.0_WP; v_vol = 0.0_WP; w_vol = 0.0_WP ! Average velocity inside struct

         ! Find new structure with matching newid
         do n=1,this%strack%nstruct
            ! Only deal with structure matching newid
            if (this%strack%struct(n)%id.eq.id) then
               
               ! Periodicity
               per_x = this%strack%struct(n)%per(1)
               per_y = this%strack%struct(n)%per(2)
               per_z = this%strack%struct(n)%per(3)
               
               ! Loop over cells in new structure and accumulate statistics
               do m=1,this%strack%struct(n)%n_

                  ! Indices of cells in structure
                  ii=this%strack%struct(n)%map(1,m) 
                  jj=this%strack%struct(n)%map(2,m) 
                  kk=this%strack%struct(n)%map(3,m)

                  ! Location of struct node
                  xtmp = this%strack%vf%cfg%xm(ii)-per_x*this%strack%vf%cfg%xL
                  ytmp = this%strack%vf%cfg%ym(jj)-per_y*this%strack%vf%cfg%yL
                  ztmp = this%strack%vf%cfg%zm(kk)-per_z*this%strack%vf%cfg%zL

                  ! Volume
                  vol_struct = vol_struct + this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  
                  ! Center of gravity
                  x_vol = x_vol + xtmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  y_vol = y_vol + ytmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  z_vol = z_vol + ztmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  
                  ! Average velocity inside struct
                  u_vol = u_vol + this%fs%U(ii,jj,kk)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  v_vol = v_vol + this%fs%V(ii,jj,kk)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  w_vol = w_vol + this%fs%W(ii,jj,kk)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
               end do
            end if
         end do

         ! Sum parallel stats
         call MPI_ALLREDUCE(MPI_IN_PLACE,vol_struct,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,x_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,y_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,z_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,u_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,v_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE,w_vol,1,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         
         ! Moments of inertia
         Imom=0.0_WP
         do n=1,this%strack%nstruct
            ! Only deal with structure matching newid
            if (this%strack%struct(n)%id.eq.id) then

               ! Periodicity
               per_x = this%strack%struct(n)%per(1)
               per_y = this%strack%struct(n)%per(2)
               per_z = this%strack%struct(n)%per(3)
                  
               ! Loop over cells in new structure and accumulate statistics
               do m=1,this%strack%struct(n)%n_

                  ! Indices of cells in structure
                  ii=this%strack%struct(n)%map(1,m) 
                  jj=this%strack%struct(n)%map(2,m) 
                  kk=this%strack%struct(n)%map(3,m)

                  ! Location of struct node
                  xtmp = this%strack%vf%cfg%xm(ii)-per_x*this%strack%vf%cfg%xL-x_vol/vol_struct
                  ytmp = this%strack%vf%cfg%ym(jj)-per_y*this%strack%vf%cfg%yL-y_vol/vol_struct
                  ztmp = this%strack%vf%cfg%zm(kk)-per_z*this%strack%vf%cfg%zL-z_vol/vol_struct

                  ! Moment of Inertia
                  Imom(1,1) = Imom(1,1) + (ytmp**2 + ztmp**2)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  Imom(2,2) = Imom(2,2) + (xtmp**2 + ztmp**2)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  Imom(3,3) = Imom(3,3) + (xtmp**2 + ytmp**2)*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  
                  Imom(1,2) = Imom(1,2) - xtmp*ytmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  Imom(1,3) = Imom(1,3) - xtmp*ztmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
                  Imom(2,3) = Imom(2,3) - ytmp*ztmp*this%strack%vf%cfg%vol(ii,jj,kk)*this%strack%vf%VF(ii,jj,kk)
               end do 
            end if
         end do

         ! Sum parallel stats on Imom
         do n=1,3
            call MPI_ALLREDUCE(MPI_IN_PLACE,Imom(:,n),3,MPI_REAL_WP,MPI_SUM,this%strack%vf%cfg%comm,ierr)
         end do

         ! Characteristic lengths and principle axes
         ! Eigenvalues/eigenvectors of moments of inertia tensor
         A = Imom
         n = 3
         call dsyev('V','U',n,Imom,n,d,work,lwork,info)
         ! Get rid of very small negative values (due to machine accuracy)
         d = max(0.0_WP,d)
         ! Store characteristic lengths
         lengths(1) = sqrt(5.0_WP/2.0_WP*abs(d(2)+d(3)-d(1))/vol_struct)
         lengths(2) = sqrt(5.0_WP/2.0_WP*abs(d(3)+d(1)-d(2))/vol_struct)
         lengths(3) = sqrt(5.0_WP/2.0_WP*abs(d(1)+d(2)-d(3))/vol_struct)
         ! Zero out length in 3rd dimension if 2D
         if (this%strack%vf%cfg%nx.eq.1.or.this%strack%vf%cfg%ny.eq.1.or.this%strack%vf%cfg%nz.eq.1) lengths(3)=0.0_WP
         ! Store principal axes
         axes(:,:) = A

         ! Finish computing quantities
         stats%vol     = vol_struct
         stats%x_cg    = x_vol/vol_struct
         stats%y_cg    = y_vol/vol_struct
         stats%z_cg    = z_vol/vol_struct
         stats%u_avg   = u_vol/vol_struct
         stats%v_avg   = v_vol/vol_struct
         stats%w_avg   = w_vol/vol_struct
         stats%Imom    = Imom 
         stats%lengths = lengths
         stats%axes    = axes

      end subroutine compute_struct_stats

   end subroutine analyze_merge_split

   subroutine dump_gtda(this,id,silo_case)
      !use dump_merge_split
      !use dump_struct
      !use dump_visit
      !use multiphase_fluxes
      !use interface_smoothing
      use stracker_class,    only: stracker
      use vfs_class,         only: vfs,VFlo
      use mpi_f08, only: MPI_INTEGER,MPI_LOGICAL,MPI_LOR,MPI_REAL
      use messager, only: die
      use interface_smoothing
      implicit none

      class(ljcf), intent(inout) :: this

      integer, intent(in) :: id
      integer :: i,j,k,ii,jj,kk
      integer :: c,cc,ntet,n,nn,nnn,np,nd,dim,ntri
      integer :: case,case2,v1,v2
      integer :: nAllocated,nListAllocated
      integer :: root = 0
      integer,  dimension(:), allocatable :: tmpNodeList
      integer,  dimension(:,:), allocatable :: tmpZoneIndex
      real(SP), dimension(:), allocatable :: tmpNode
      !real(WP), parameter :: VOFlo_visit=0.00001
      !real(WP), parameter :: VOFhi_visit=0.99999
      real(WP), dimension(3,8) :: verts
      real(WP), dimension(3,5) :: tverts,tverts2
      real(WP), dimension(3,4) :: mytet
      real(WP), dimension(4) :: d,plane
      real(WP), dimension(3) :: node,cen
      real(WP) :: mu
      real(WP) :: sdx2,sdy2,sdz2  

      integer, intent(in) :: silo_case
      integer, dimension(this%cfg%nproc) :: g_nStruct_array,displace,g_nMerge_array,displace_nodelist,displace_zonelist,g_nNode_array,g_nZone_array,g_nodeList_array
      real(SP),  dimension(:), allocatable :: g_xNode,g_yNode,g_zNode
      integer, dimension(:), allocatable :: g_nodeList
      integer :: g_nNodes,g_nZones
      real(SP), dimension(:), allocatable :: Vnorm,g_Vnorm
      logical  :: include_proc, g_include_proc,dir_exists
      integer :: dbfile,err,ierr,optlist,stop_
      character(len=5) :: dirname
      character(len=43) :: siloname
      character(len=26) :: folder,buffer
      real(SP), dimension(:), allocatable :: ind
      integer  :: ms_nout_time
      integer :: nZoneAllocated
   
      


      !call timing_start('gtda') 

      ! Initialize SILO LID and time lists if they do not exist yet
      !if (is_gtda .and. .not. LID_silo_exists .and. this%cfg%amRoot) then
      !   allocate(LID_silo(this%strack%nstruct)); LID_silo = this%strack%id_rmp
      !   allocate(time_silo(this%strack%nstruct)); time_silo = this%time%t

      !   nAlloc_silo = this%strack%nstruct
      !   LID_silo_exists = .true.
      !end if 

      ! Loop over split/merge/update structures
      do np=1,this%strack%nstruct
         ! Initialize stop flag
         !stop_ = 0

         ! start marching thru list, starting at first_struct
         !my_struct => first_struct
         !do m=1,this%strack%struct(np)
         ! Initialize arrays
         nAllocated    =1000
         nListAllocated=1000
         allocate(xNode   (nAllocated)); xNode = 0
         allocate(yNode   (nAllocated)); yNode = 0
         allocate(zNode   (nAllocated)); zNode = 0
         allocate(Vnorm   (nAllocated)); Vnorm = 0
         allocate(nodeList(nListAllocated)); nodeList = 0
         ! Allocate buffers for interface variables
         if (.not.allocated(zoneIndex)) then
            !print*, "Allocating zoneIndex for the first time"
            !call die("stopping")
            nZoneAllocated=1000
            allocate(zoneIndex(3,nZoneAllocated))
         else
            nZoneAllocated = size(zoneIndex,2)
         end if
         nNodes=0
         nZones=0

         ! Logical initialize = false
         include_proc = .false.
         !do while(associated(this%strack%struct(np)).and.stop_.eq.0)
            !print*, "Processor ", this%cfg%rank, " checking structure ", np, " with id ", this%strack%struct(np)%id
            !print*, "Other id" , id
            if (this%strack%struct(np)%id.eq.id) then
               ! Include these processors in MPI calls below
               include_proc = .true.
               
               ! Loop over nodes in struct
               do nd=1,this%strack%struct(np)%n_ 
               i = this%strack%struct(np)%map(1,nd)
               j = this%strack%struct(np)%map(2,nd)
               k = this%strack%struct(np)%map(3,nd)
               
               ! Output the interface using default or marching cubes
               !select case(trim(gtda_output_type))
               
               !case('Marching cubes')
                  if (this%vf%VF(i,j,k).ge.VFlo) then 
                     ! Must loop over all surrounding cells for marching cubes
                     do kk = k-1,k+1
                        do jj = j-1,j+1
                           do ii = i-1,i+1
                              ! Call marching_cubes
                              !print*, "cell", ii, jj, kk
                              !print*, "nZoneAllocated: ", nZoneAllocated, "nZones: ", nZones
                              call marching_cubes2(this%vf,ii,jj,kk)
                              !print*, "Marching cubes called for cell "
                              !if (nZoneAllocated.gt.1000) then
                              !   call die("stopping after marching cubes for testing")
                              !end if
                              ! Reallocate arrays if necessary
                              if (nAllocated-nNodes.lt.200) then
                                 allocate(tmpNode(nAllocated))
                                 tmpNode=xNode; deallocate(xNode); allocate(xNode(nAllocated+1000)); xNode(1:nAllocated)=tmpNode
                                 tmpNode=yNode; deallocate(yNode); allocate(yNode(nAllocated+1000)); yNode(1:nAllocated)=tmpNode
                                 tmpNode=zNode; deallocate(zNode); allocate(zNode(nAllocated+1000)); zNode(1:nAllocated)=tmpNode
                                 tmpNode=Vnorm; deallocate(Vnorm); allocate(Vnorm(nAllocated+1000)); Vnorm(1:nAllocated)=tmpNode
                                 deallocate(tmpNode)
                                 nAllocated=nAllocated+1000
                              end if
                              if (nListAllocated-3*nZones.lt.200) then
                                 allocate(tmpNodeList(nListAllocated))
                                 tmpNodeList=nodeList; 
                                 deallocate(nodeList); 
                                 allocate(nodeList(nListAllocated+1000)); 
                                 nodeList(1:nListAllocated)=tmpNodeList
                                 deallocate(tmpNodeList)
                                 nListAllocated=nListAllocated+1000
                              end if
                              if (nZoneAllocated-nZones.lt.200) then
                                 allocate(tmpZoneIndex(3,nZoneallocated))
                                 tmpZoneIndex=zoneIndex
                                 deallocate(zoneIndex)
                                 allocate(zoneIndex(3,nZoneAllocated+1000))
                                 Zoneindex(:,1:nZoneAllocated)=tmpZoneIndex
                                 deallocate(tmpZoneIndex)
                                 nZoneAllocated=nZoneAllocated+1000
                              end if
                           end do 
                        end do 
                     end do 

                  end if ! VOF.ge.VOFlo_Visit
               !end select 
               end do ! i,j,k
               ! Exit loop over structures
               !stop_ = 1
               !print*, "Number of nodes after looping through i,j,k", nNodes
            end if
         ! Go to next structure
         !end do
         !my_struct => my_struct%next
         !end do ! do while associated(my_struct)

         ! Create global logical - if any processors have split structures, then true
         call MPI_ALLREDUCE(include_proc,g_include_proc,1,MPI_LOGICAL,MPI_LOR,this%cfg%comm,ierr)

         ! ! Add a zero-area tri if this proc doesn't have one
         ! if (nZones.eq.0) then
         !    nZones=1
         !    nNodes=3
         !    xNode(1:3)=real(xm(imin_),SP)
         !    yNode(1:3)=real(ym(jmin_),SP)
         !    zNode(1:3)=real(zm(kmin_),SP)
         !    nodeList(1:3)=(/1,2,3/)
         !    zoneIndex(:,nZones)=(/ imin_,jmin_,kmin_ /)
         ! end if 
         !print*, "Number of zones after marching cubes", nZones
         ! Calculate normal velocity on each zone
         do n = 1,nZones
            !do c = 1,8
            i = zoneIndex(1,n)
            j = zoneIndex(2,n)
            k = zoneIndex(3,n)
            plane = getPlane(this%vf%liquid_gas_interface(i,j,k),0)
            Vnorm(i) = real(this%fs%U(i,j,k)*plane(1) &
                        + this%fs%V(i,j,k)*plane(2) &
                        + this%fs%W(i,j,k)*plane(3),SP)
            !end do
         end do 
         !print*, "Normal velocities calculated for ", nZones, " zones"
         !print*, "Number of nodes on proc ", this%cfg%rank, " is ", nNodes
         ! ---------------------- !
         ! Gather lists onto root !
         ! ---------------------- !
         if (g_include_proc) then

            ! Initialize arrays
            g_nNode_array=0
            g_nodeList_array=0
            g_nZone_array=0
            
            ! Communicate nNodes
            call MPI_AllGather(nNodes,1,MPI_INTEGER,g_nNode_array,1,MPI_INTEGER,this%cfg%comm,ierr)
            g_nNodes=sum(g_nNode_array)

            !print*, "print off global nNodes"  
            !print*, g_nNodes 
         
            ! Communicate nZones
            call MPI_AllGather(nZones,1,MPI_INTEGER,g_nZone_array,1,MPI_INTEGER,this%cfg%comm,ierr)
            g_nZones=sum(g_nZone_array)
            g_nodeList_array = 3*g_nZone_array     
            
            ! Allocate global lists
            allocate(g_nodeList(3*g_nZones)); g_nodeList=0
            allocate(g_Vnorm(g_nZones)); g_Vnorm=0
            allocate(g_xNode(g_nNodes)); g_xNode=0
            allocate(g_yNode(g_nNodes)); g_yNode=0
            allocate(g_zNode(g_nNodes)); g_zNode=0
            
            ! Compute displacements for data from each processor in global arrays
            do n=1,this%cfg%nproc
               displace(n)          = sum(g_nNode_array(1:n-1))
               displace_zonelist(n) = sum(g_nZone_array(1:n-1))
               displace_nodelist(n) = sum(g_nodeList_array(1:n-1))
            end do

            nodeList = nodeList+displace(this%cfg%rank+1)
            !print*, xNode(1:nNodes)
            !print*, yNode(1:nNodes)
            !print*, zNode(1:nNodes)

            ! Gather node lists onto root
            call MPI_GATHERV(nodeList(1:3*nZones),3*nZones,MPI_INTEGER, &
               g_nodeList,g_nodeList_array,displace_nodelist,MPI_INTEGER,root,this%cfg%comm,ierr)
            call MPI_GATHERV(Vnorm(1:nZones),nZones,MPI_REAL,    &
               g_Vnorm,g_nZone_array,displace_zonelist,MPI_REAL,root,this%cfg%comm,ierr)
            call MPI_GATHERV(xNode(1:nNodes),nNodes,MPI_REAL,    &
               g_xNode,g_nNode_array,displace,MPI_REAL,root,this%cfg%comm,ierr)
            call MPI_GATHERV(yNode(1:nNodes),nNodes,MPI_REAL,    &
               g_yNode,g_nNode_array,displace,MPI_REAL,root,this%cfg%comm,ierr)
            call MPI_GATHERV(zNode(1:nNodes),nNodes,MPI_REAL,    &
               g_zNode,g_nNode_array,displace,MPI_REAL,root,this%cfg%comm,ierr)
            !print*, "Finished gathering lists onto root"
            !print*, "Total number of nodes on root is ", g_nNodes
            !print*, "Total number of zones on root is ", g_nZones
            if (this%cfg%amRoot) then  
               if (g_nNodes.gt.0) then  
                  ! Update output counter
                  ms_nout_time=ms_nout_time+1
                  
                  ! Create directory for SILO files created this timestep
                  write(buffer,'(ES12.5)') this%time%t
                  folder = 'gtda/events_'//trim(adjustl(buffer))
                  inquire(file=folder,exist=dir_exists)
                  if (.not.dir_exists) then
                     call execute_command_line('mkdir -p '//trim(folder))
                  end if 
                  !print*, 'creating file with id ',id
                  !print*, this%strack%struct(np)%id
                  ! Filename is event type and new LID
                  write(buffer,'(I7.7)') id

                  select case(silo_case)
                  case(1) ! Write file for merge event
                     write(siloname,'(A,A,A)') trim(folder)//'/merge',trim(adjustl(buffer)),'.silo'
                  case(2) ! Write file for split event 
                     write(siloname,'(A,A,A)') trim(folder)//'/split',trim(adjustl(buffer)),'.silo'
                  case(3) ! Write file for update 
                     write(siloname,'(A,A,A)') trim(folder)//'/update',trim(adjustl(buffer)),'.silo'
                  end select    
                  !print*, 'wrote file'
                  ! Create the silo database
                  err = dbcreate(siloname, len_trim(siloname), DB_CLOBBER, DB_LOCAL,"Silo database created with NGA2", 30, DB_HDF5, dbfile)
                     if(dbfile.eq.-1) call die('Could not create Silo file!')
                  ! ierr = dbclose(dbfile)
                  !print*, 'writing interface for structure: g_nZones=', g_nZones, ' g_nNodes=', g_nNodes, ' g_nodeList size=', size(g_nodeList), ' g_Vnorm size=', size(g_Vnorm)
                  ! Write interface as unstructured mesh made of tetrahedra
                  err = dbputzl2(dbfile,"zonelist",8,g_nZones,3,g_nodeList,3*g_nZones,1,0,0 &
                        ,DB_ZONETYPE_TRIANGLE,3,g_nZones,1,DB_F77NULL,ierr) 
                  !print*, 'finished writing zonelist'
                  !print*, g_xNode(1:g_nNodes)
                  !print*, "finished printing xnodes"
                  !print*, g_yNode(1:g_nNodes)
                  !print*, "finished printing ynodes"
                  !print*, g_zNode(1:g_nNodes)
                  !print*, "finished printing znodes"
                  err = dbputum(dbfile,"Interface",9,3,g_xNode(1:g_nNodes),g_yNode(1:g_nNodes),g_zNode(1:g_nNodes) &
                        ,"xInt",4,"yInt",4,"zInt",4,DB_FLOAT,g_nNodes,g_nZones,"zonelist",8,DB_F77NULL,0,DB_F77NULL,ierr)
                  !print*, 'finished writing unstructured mesh'
                  err = dbputuv1(dbfile,"Vnorm",5,"Interface",9,g_Vnorm(1:g_nZones),g_nZones &
                        ,DB_F77NULL,0,DB_FLOAT,DB_ZONECENT,DB_F77NULL, ierr)
                  !print*, 'finished writing Vnorm variable'
                  ! Close group silo file
                  ierr = dbclose(dbfile)
                  !print*, 'closed file'
               end if 
            end if ! iroot.eq.irank 
            ! Deallocate global lists associated with this structure
            deallocate(g_nodeList,g_xNode,g_yNode,g_zNode,g_Vnorm)  
         end if 
         ! Deallocate local lists associated with this structure
         deallocate(nodeList,xNode,yNode,zNode,Vnorm)    

      end do ! np=1,nUp
      !call timing_stop('gtda')
      !print*, "Finished dumping GTDA for this event", "id=", id
   end subroutine dump_gtda

   ! Subroutine for keeping track of LIDs for SILO outputs
   !subroutine silo_lists_update(this,silo_case)
   !   !use dump_merge_split
   !   use quicksort
!
   !   integer, intent(in) :: silo_case 
   !   logical :: added = .false.
   !   
   !   if (this%cfg%amRoot) then 
   !      ! Check if we need to resize array
   !      if (nAlloc_silo-this%strack%nstruct.lt.200) then
   !         allocate(g_tmpInt(nAlloc_silo))
   !         allocate(g_tmpWP(nAlloc_silo))
   !         g_tmpInt=LID_silo; deallocate(LID_silo) ; allocate(LID_silo(nAlloc_silo+1000)) ; LID_silo = 0      ; LID_silo(1:nAlloc_silo)=g_tmpInt(1:nAlloc_silo)
   !         g_tmpWP=time_silo; deallocate(time_silo); allocate(time_silo(nAlloc_silo+1000)); time_silo = 0.0_WP; time_silo(1:nAlloc_silo)=g_tmpWP(1:nAlloc_silo)
   !         deallocate(g_tmpInt,g_tmpWP)
   !         nAlloc_silo = nAlloc_silo + 10
   !      end if
!
   !      select case(silo_case)
   !      case(1) ! Merge case
   !         !g_nStruct = g_nStruct - nUp ! Update g_nStruct from last timestep
   !         ! Remove old LIDs from list of LIDs
   !         update: do i = 1,this%strack%nstruct
   !            do j = 1,nAlloc_silo
   !               if (LID_silo(j).eq.this%strack%oldids(i)) then
   !                  LID_silo(j:this%strack%nstruct) = LID_silo(j+1:this%strack%nstruct+1)
   !                  LID_silo(this%strack%nstruct+1) = 0
   !                  time_silo(j:this%strack%nstruct) = time_silo(j+1:this%strack%nstruct+1) ! Shift lists up
   !                  time_silo(this%strack%nstruct+1) = 0.0_WP
   !                  cycle update 
   !               end if 
   !            end do 
   !         end do update 
!
   !      case(2) ! Split case 
   !         ! Remove old LIDs from list of LIDs
   !         update_remove: do i = 1,this%strack%nstruct
   !            do j = 1,nAlloc_silo
   !               if (LID_silo(j).eq.this%strack%oldids(i)) then
   !                  LID_silo(j:this%strack%nstruct) = LID_silo(j+1:this%strack%nstruct+1)
   !                  LID_silo(this%strack%nstruct+1) = 0
   !                  time_silo(j:this%strack%nstruct) = time_silo(j+1:this%strack%nstruct+1) ! Shift lists up
   !                  time_silo(this%strack%nstruct+1) = 0.0_WP
   !                  cycle update_remove 
   !               end if 
   !            end do 
   !         end do update_remove 
   !         ! Add new LIDs to list of LIDs  
   !         update_add: do i = 1,this%strack%nstruct
   !            do j = 1,nAlloc_silo
   !               if (LID_silo(j).eq.0) then
   !                  LID_silo(j)  = this%strack%id_rmp(i)
   !                  time_silo(j) = time
   !                  cycle update_add
   !               end if 
   !            end do 
   !         end do update_add
!
   !      case(3) ! Update case
   !         ! Update times when we call silo_lists_check
   !         do i = 1,this%strack%nstruct
   !            time_silo(i) = this%time%t ! update the first nUp entries
   !         end do 
   !         call quick_sort(time_silo(1:this%strack%nstruct),LID_silo(1:this%strack%nstruct))
!
   !      case(4) ! Struct removed from domain in dump_struct_remove 
   !         ! Remove old LIDs from list of LIDs
   !         update_remove2: do i = 1,this%strack%nstruct
   !            do j = 1,nAlloc_silo
   !               if (LID_silo(j).eq.up_LIDo(i)) then
   !                  LID_silo(j:this%strack%nstruct) = LID_silo(j+1:this%strack%nstruct+1)
   !                  LID_silo(this%strack%nstruct+1) = 0
   !                  time_silo(j:this%strack%nstruct) = time_silo(j+1:this%strack%nstruct+1) ! Shift lists up
   !                  time_silo(this%strack%nstruct+1) = 0.0_WP
   !                  cycle update_remove2 
   !               end if 
   !            end do 
   !         end do update_remove2 
   !      end select
!
   !   end if 
   !end subroutine silo_lists_update

   !subroutine silo_list_check
   !   use dump_merge_split
   !   use dump_struct 
   !   use dump_gtda 
!
   !   nUp = 0 ! Reset update counter
!
   !   deallocate(this%strack%id_rmp)
   !   allocate(this%strack%id_rmp(this%strack%nstruct)) 
   !   if (this%cfg%amRoot) then
   !      do i = 1,this%strack%nstruct
   !         if ((this%time%t - time_silo(i)) .ge. update_time) then 
   !            !nUp = nUp + 1 
   !            !up_LIDn(nUp) = LID_silo(i)
   !            this%strack%id_rmp(i) = LID_silo(i)
   !         else 
   !            exit ! Stop if not enough time has passed
   !         end if 
   !      end do   
   !   end if 
   !   
   !   call MPI_Bcast(nUp,1,MPI_INTEGER,iroot-1,comm,ierr)
   !   call MPI_Bcast(this%strack%id_rmp(1:this%strack%nstruct),this%strack%nstruct,MPI_INTEGER,iroot-1,comm,ierr) 
!
   !   if (this%strack%nstruct.gt.0) then
   !      ! Associate SIDs with LIDs for SILO file extraction
   !      deallocate(this%strack%id)
   !      allocate(this%strack%id(this%strack%nstruct)); this%strack%id=0 
   !      update: do n=1,this%strack%nstruct     
   !         do k=kmin_,kmax_
   !            do j=jmin_,jmax_
   !               do i=imin_,imax_ 
   !                  if (nint(this%strack%id(i,j,k)).eq.this%strack%id_rmp(n)) then 
   !                     this%strack%id(n) = this%strack%id(i,j,k) 
   !                     cycle update 
   !                  end if 
   !               end do 
   !            end do   
   !         end do           
   !      end do update 
   !   end if 
!
   !end subroutine silo_list_check

   !> Initialization of ljcf simulation
   subroutine init(this)
      implicit none
      class(ljcf), intent(inout) :: this
      !> Added parameters for GTDA dumping testing
      real(WP), dimension(3) :: center
      real(WP) :: radius
      ! Create the ljcf mesh
      create_config: block
         use sgrid_class, only: cartesian,sgrid
         use param,       only: param_read
         use parallel,    only: group
         real(WP), dimension(:), allocatable :: x,y,z
         integer, dimension(3) :: partition
         type(sgrid) :: grid
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,xlig
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1)); call param_read('X ljcf',xlig)
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-xlig
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=3,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.true.,name='ljcf')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid without walls
         this%cfg=config(grp=group,decomp=partition,grid=grid)

      end block create_config
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         call param_read('Max time',this%time%tmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays

      ! Set up walls before solvers are initialized
      create_walls: block
         use param, only: param_read,param_getsize
         integer :: i,j,k,njet
         ! Initialize liquid jet(s)
         call param_read('Jet diameter',this%djet)
         njet = param_getsize('Jet location')
         allocate(this%xjet(njet))
         call param_read('Jet location',this%xjet)
         call param_read('Froude number',this%gravity); this%gravity = 1.0_WP/this%gravity**2
         call param_read('Liquid Volume',this%liqVol)
         this%liqVolInjected = 0.0_WP
         ! Number of wall cells
         call param_read('Wall cells in domain', this%nwall, default=0)
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  if (wall(this%cfg%pgrid,i,j,k)) then
                  this%cfg%VF(i,j,k)=0.0_WP
                  end if
               end do
            end do
         end do
      end block create_walls
            
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use vfs_class, only: remap_storage,VFlo,VFhi,plicnet,r2pnet
         use mms_geom,  only: cube_refine_vol
         use param,     only: param_read
         use mpi_f08, only: MPI_INTEGER,MPI_MAX
         integer :: i,j,k,n,si,sj,sk,ierr
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call this%vf%initialize(cfg=this%cfg,reconstruction_method=plicnet,transport_method=remap_storage,name='VOF')
         this%vf%thin_thld_min=0.0_WP
         this%vf%flotsam_thld=0.0_WP
         this%vf%maxcurv_times_mesh=1.0_WP
         ! Create structure tracker
         call this%strack%initialize(vf=this%vf,phase=0,make_label=label_liquid,name='stracker_test')
         ! Initialize the interface to a ljcf
         do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
            do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
               do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[this%vf%cfg%x(i+si),this%vf%cfg%y(j+sj),this%vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  if (j.lt.this%vf%cfg%jmin) then
                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_halfdrop,0.0_WP,amr_ref_lvl)
                  else
                     ! do nothing
                  end if
                  this%vf%VF(i,j,k)=vol/this%vf%cfg%vol(i,j,k)
                  if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
                     this%vf%Lbary(:,i,j,k)=v_cent
                     this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
                  else
                     this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                     this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
                  end if
                   ! Set stracker id
                     if (vol.gt.0.0_WP) then
                        this%strack%id(i,j,k)=1
                     end if
               end do
            end do
         end do

         ! Initialize id counter to be consistent with id's
         this%strack%idcount=maxval(this%strack%id)
         call MPI_ALLREDUCE(maxval(this%strack%id),this%strack%idcount,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
         ! Update the band
         call this%vf%update_band()
         ! Perform interface reconstruction from VOF field
         call this%vf%build_interface()
         ! Set interface planes at the boundaries
         call this%vf%set_full_bcond()
         ! Now apply Neumann condition on interface at inlet to have proper round injection
         neumann_irl: block
            use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
            &                                setNumberOfPlanes,setPlane,matchVolumeFraction
            real(WP), dimension(1:4) :: plane
            type(RectCub_type) :: cell
            call new(cell)
            if (this%vf%cfg%iproc.eq.1) then
               do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
                  do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                     do i=this%vf%cfg%imino,this%vf%cfg%imin-1
                        ! Extract plane data and copy in overlap
                        plane=getPlane(this%vf%liquid_gas_interface(this%vf%cfg%imin,j,k),0)
                        call construct_2pt(cell,[this%vf%cfg%x(i  ),this%vf%cfg%y(j  ),this%vf%cfg%z(k  )],&
                        &                       [this%vf%cfg%x(i+1),this%vf%cfg%y(j+1),this%vf%cfg%z(k+1)])
                        plane(4)=dot_product(plane(1:3),[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)])
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,plane(1:3),plane(4))
                        call matchVolumeFraction(cell,this%vf%VF(i,j,k),this%vf%liquid_gas_interface(i,j,k))
                     end do
                  end do
               end do
            end if
         end block neumann_irl
         ! Create discontinuous polygon mesh from IRL interface
         call this%vf%polygonalize_interface()
         ! Calculate distance from polygons
         call this%vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call this%vf%subcell_vol()
         ! Calculate curvature
         call this%vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call this%vf%reset_volume_moments()
      end block create_and_initialize_vof
!      create_and_initialize_vof: block
!         use vfs_class, only:lvira,plicnet,remap_storage,VFhi,VFlo,r2pnet
!         use random, only: random_uniform
!         use mms_geom, only: cube_refine_vol
!         use precision, only: I4
!         use MPI, only: MPI_INTEGER,MPI_MAX
!         use param, only: param_read
!         real(WP), dimension(3,8) :: cube_vertex
!         real(WP), dimension(3) :: v_cent,a_cent
!         real(WP) :: vol,area
!         integer, parameter :: amr_ref_lvl=4
!         integer :: i,j,k,n,si,sj,sk
!         integer :: nD,nDrop
!         integer(kind=I4), allocatable, dimension(:) :: seed
!         integer(kind=I4) :: nseed
!         integer :: ierr
!         ! Create a VOF solver with r2p reconstruction
!         call this%vf%initialize(cfg=this%cfg,reconstruction_method=lvira,transport_method=remap_storage,name='VOF')
!         ! Create structure tracker
!         this%vf%thin_thld_min=0.0_WP
!         this%vf%flotsam_thld=0.0_WP
!         this%vf%maxcurv_times_mesh=1.0_WP
!         call this%strack%initialize(vf=this%vf,phase=0,make_label=label_liquid,name='stracker_test')
!         ! Initialize our bubble via r2p planes
!         call param_read('Droplet diameter',radius); radius=radius/2.0_WP
!         call param_read('Number of droplet',nDrop);
!         ! Provide seed for random number generator
!         call random_seed(size=nseed)
!         allocate(seed(nseed))
!         seed(:)=1
!         call random_seed(put=seed)
!         do nD=1,nDrop
!            center=[random_uniform(this%vf%cfg%x(this%vf%cfg%imin),this%vf%cfg%x(this%vf%cfg%imax+1)), &
!                    random_uniform(this%vf%cfg%y(this%vf%cfg%jmin),this%vf%cfg%y(this%vf%cfg%jmax+1)), &
!                    random_uniform(this%vf%cfg%z(this%vf%cfg%kmin),this%vf%cfg%z(this%vf%cfg%kmax+1))  ]
!                     
!            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
!               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
!                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
!                     ! Set cube vertices
!                     n=0
!                     do sk=0,1
!                        do sj=0,1
!                           do si=0,1
!                              n=n+1; cube_vertex(:,n)=[this%vf%cfg%x(i+si),this%vf%cfg%y(j+sj),this%vf%cfg%z(k+sk)]
!                           end do
!                        end do
!                     end do
!                     ! Call adaptive refinement code to get volume and barycenters recursively
!                     vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
!                     !!!!! levelset_sphere
!                     call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_sphere,0.0_WP,amr_ref_lvl)
!                     this%vf%VF(i,j,k)=min(1.0_WP,this%vf%VF(i,j,k)+vol/this%vf%cfg%vol(i,j,k))
!                     if (this%vf%VF(i,j,k).ge.VFlo.and.this%vf%VF(i,j,k).le.VFhi) then
!                        this%vf%Lbary(:,i,j,k)=v_cent
!                        this%vf%Gbary(:,i,j,k)=([this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]-this%vf%VF(i,j,k)*this%vf%Lbary(:,i,j,k))/(1.0_WP-this%vf%VF(i,j,k))
!                     else
!                        this%vf%Lbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
!                        this%vf%Gbary(:,i,j,k)=[this%vf%cfg%xm(i),this%vf%cfg%ym(j),this%vf%cfg%zm(k)]
!                     end if
!                     ! Set stracker id
!                     if (vol.gt.0.0_WP) then
!                        this%strack%id(i,j,k)=1
!                        !if (nD.eq.4) then    ! Test growth (merge2)
!                        !   strack%id(i,j,k)=0
!                        !end if   
!                     end if
!                  end do
!               end do
!            end do
!         end do
!         call this%vf%cfg%sync(this%vf%VF)
!         call this%vf%cfg%sync(this%strack%id)
!
!         ! Initialize id counter to be consistent with id's
!         this%strack%idcount=maxval(this%strack%id)
!         call MPI_ALLREDUCE(maxval(this%strack%id),this%strack%idcount,1,MPI_INTEGER,MPI_MAX,this%vf%cfg%comm,ierr)
!         ! Update the band
!         call this%vf%update_band()
!         ! Perform interface reconstruction from VOF field
!         call this%vf%build_interface()
!         ! Set interface planes at the boundaries
!         call this%vf%set_full_bcond()
!         ! Create discontinuous polygon mesh from IRL interface
!         call this%vf%polygonalize_interface()
!         ! Calculate distance from polygons
!         call this%vf%distance_from_polygon()
!         ! Calculate subcell phasic volumes
!         call this%vf%subcell_vol()
!         ! Calculate curvature
!         call this%vf%get_curvature()
!         ! Reset moments to guarantee compatibility with interface reconstruction
!         call this%vf%reset_volume_moments()
!      end block create_and_initialize_vof
      
      ! Create an iterator for removing VOF at edges
      create_iterator: block
         this%vof_removal_layer=iterator(this%cfg,'VOF removal',vof_removal_layer_locator)
      end block create_iterator
      
      
      ! Create a multiphase flow solver with bconds
      create_flow_solver: block
         use mathtools,       only: Pi
         use param,           only: param_read
         use tpns_class,      only: dirichlet,clipped_neumann,bcond
         use hypre_str_class, only: pcg_pfmg2
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         ! Create flow solver
         this%fs=tpns(cfg=this%cfg,name='Two-phase NS')
         ! Set fluid properties
         this%fs%rho_g=1.0_WP; call param_read('Density ratio',this%fs%rho_l)
         call param_read('Reynolds number',this%fs%visc_g); this%fs%visc_g=1.0_WP/this%fs%visc_g
         call param_read('Viscosity ratio',this%fs%visc_l); this%fs%visc_l=this%fs%visc_g*this%fs%visc_l
         call param_read('Weber number',this%fs%sigma); this%fs%sigma=1.0_WP/this%fs%sigma
         ! Define inflow boundary condition on the left
         call this%fs%add_bcond(name='inflow',type=dirichlet,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         ! Define outflow boundary condition on the right
         call this%fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         ! Define jet boundary condition on the bottom
         call this%fs%add_bcond(name='jet'    ,type=dirichlet,face='y',dir=-1,canCorrect=.false.,locator=jet_bdy)
         ! Define gravity as vector for flow solver
         this%fs%gravity(2) = this%gravity

         ! Configure pressure solver
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         this%ps%maxlevel=16
         call param_read('Pressure iteration',this%ps%maxit)
         call param_read('Pressure tolerance',this%ps%rcvg)
         ! Configure implicit velocity solver
         !this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps)!,implicit_solver=this%vs)
         ! Zero initial field
         this%fs%U=0.0_WP; this%fs%V=0.0_WP; this%fs%W=0.0_WP
         ! Apply convective velocity
         call this%fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%U(i,j,k)=1.0_WP
         end do
         ! Apply jet velocity
         call this%fs%get_bcond('jet',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            this%fs%V(i,j,k)=0 ! Start with zero velocity this%Vjet 
         end do
         ! Apply all other boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         ! Adjust MFR for global mass balance
         call this%fs%correct_mfr()
         ! Compute divergence
         call this%fs%get_div()
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      end block create_flow_solver

      ! Create CCL
      create_ccl: block
         ! Initialize CCL
         call this%ccl%initialize(pg=this%cfg%pgrid,name='ccl')
      end block create_ccl
      
      ! Handle restart/saves here
      handle_restart: block
         use param,                 only: param_read
         use string,                only: str_medium
         use filesys,               only: makedir,isdir
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane 
         character(len=str_medium) :: timestamp
         integer, dimension(3) :: iopartition
         real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
         real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
         integer :: i,j,k
         ! Create event for saving restart files
         this%save_evt=event(this%time,'Restart output')
         call param_read('Restart output period',this%save_evt%tper)
         ! Check if we are restarting
         call param_read('Restart from',timestamp,default='')
         this%restarted=.false.; if (len_trim(timestamp).gt.0) this%restarted=.true.
         ! Read in the I/O partition
         call param_read('I/O partition',iopartition)
         ! Perform pardata initialization
         if (this%restarted) then
            ! Read in the file
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,fdata='restart/data_'//trim(timestamp))
            ! Read in the planes directly and set the IRL interface
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P11',var=P11)
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P12',var=P12)
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P13',var=P13)
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P14',var=P14)
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P21',var=P21)
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P22',var=P22)
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P23',var=P23)
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_)); call this%df%pull(name='P24',var=P24)
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     ! Check if the second plane is meaningful
                     if (this%vf%two_planes.and.P21(i,j,k)**2+P22(i,j,k)**2+P23(i,j,k)**2.gt.0.0_WP) then
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),2)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),1,[P21(i,j,k),P22(i,j,k),P23(i,j,k)],P24(i,j,k))
                     else
                        call setNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k),1)
                        call setPlane(this%vf%liquid_gas_interface(i,j,k),0,[P11(i,j,k),P12(i,j,k),P13(i,j,k)],P14(i,j,k))
                     end if
                  end do
               end do
            end do
            call this%vf%sync_interface()
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
            ! Reset moments
            call this%vf%reset_volume_moments()
            ! Update the band
            call this%vf%update_band()
            ! Create discontinuous polygon mesh from IRL interface
            call this%vf%polygonalize_interface()
            ! Calculate distance from polygons
            call this%vf%distance_from_polygon()
            ! Calculate subcell phasic volumes
            call this%vf%subcell_vol()
            ! Calculate curvature
            call this%vf%get_curvature()
            ! Now read in the velocity solver data
            call this%df%pull(name='U',var=this%fs%U)
            call this%df%pull(name='V',var=this%fs%V)
            call this%df%pull(name='W',var=this%fs%W)
            call this%df%pull(name='P',var=this%fs%P)
            call this%df%pull(name='Pjx',var=this%fs%Pjx)
            call this%df%pull(name='Pjy',var=this%fs%Pjy)
            call this%df%pull(name='Pjz',var=this%fs%Pjz)
            ! Apply all other boundary conditions
            call this%fs%apply_bcond(this%time%t,this%time%dt)
            ! Compute MFR through all boundary conditions
            call this%fs%get_mfr()
            ! Adjust MFR for global mass balance
            call this%fs%correct_mfr()
            ! Compute cell-centered velocity
            call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
            ! Compute divergence
            call this%fs%get_div()
            ! Also update time
            call this%df%pull(name='t' ,val=this%time%t )
            call this%df%pull(name='dt',val=this%time%dt)
            this%time%told=this%time%t-this%time%dt
            !this%time%dt=this%time%dtmax !< Force max timestep size anyway
         else
            ! We are not restarting, prepare a new directory for storing restart files
            if (this%cfg%amRoot) then
               if (.not.isdir('restart')) call makedir('restart')
            end if
            ! Prepare pardata object for saving restart files
            call this%df%initialize(pg=this%cfg,iopartition=iopartition,filename=trim(this%cfg%name),nval=2,nvar=15)
            this%df%valname=['t ','dt']
            this%df%varname=['U  ','V  ','W  ','P  ','Pjx','Pjy','Pjz','P11','P12','P13','P14','P21','P22','P23','P24']
         end if
      end block handle_restart
      
      
      ! Create surfmesh object for interface polygon output
      create_smesh: block
         use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
         integer :: i,j,k,np,nplane
         this%smesh=surfmesh(nvar=3,name='plic')
         this%smesh%varname(1)='nplane'
         this%smesh%varname(2)='thickness'
         this%smesh%varname(3)='id'
         ! Transfer polygons to smesh
         call this%vf%update_surfmesh(this%smesh)
         ! Calculate thickness
         ! call this%vf%get_thickness()
         ! Populate nplane and thickness variables
         this%smesh%var(1,:)=1.0_WP
         np=0
         do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
            do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
               do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                  if (this%cfg%VF(i,j,k).lt.2.0_WP*epsilon(1.0_WP)) cycle ! Skip cells below VF threshold
                  do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                     if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                        np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                        this%smesh%var(2,np)=0.0_WP !this%vf%thickness(i,j,k)
                     end if
                  end do
               end do
            end do
         end do
      end block create_smesh
      
      
      ! Add Ensight output
      create_ensight: block
         use param, only: param_read
         ! Create Ensight output from cfg
         this%ens_out=ensight(cfg=this%cfg,name='ljcf')
         ! Create event for Ensight output
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period',this%ens_evt%tper)
         ! Add variables to output
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('VOF',this%vf%VF)
         call this%ens_out%add_scalar('curvature',this%vf%curv)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         call this%ens_out%add_surface('plic',this%smesh)
         call this%ens_out%add_scalar('id',this%strack%id)
         ! Output to ensight
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call this%fs%get_cfl(this%time%dt,this%time%cfl)
         call this%fs%get_max()
         call this%vf%get_max()
         ! Create simulation monitor
         this%mfile=monitor(this%fs%cfg%amRoot,'simulation_atom')
         call this%mfile%add_column(this%time%n,'Timestep number')
         call this%mfile%add_column(this%time%t,'Time')
         call this%mfile%add_column(this%time%dt,'Timestep size')
         call this%mfile%add_column(this%time%cfl,'Maximum CFL')
         call this%mfile%add_column(this%fs%Umax,'Umax')
         call this%mfile%add_column(this%fs%Vmax,'Vmax')
         call this%mfile%add_column(this%fs%Wmax,'Wmax')
         call this%mfile%add_column(this%fs%Pmax,'Pmax')
         call this%mfile%add_column(this%vf%VFint,'VOF integral')
         call this%mfile%add_column(this%vf%SDint,'SD integral')
         call this%mfile%add_column(this%vof_removed,'VOF removed')
         call this%mfile%add_column(this%vf%flotsam_error,'Flotsam error')
         call this%mfile%add_column(this%vf%thinstruct_error,'Film error')
         call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
         call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
         call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
         call this%mfile%write()
         ! Create CFL monitor
         this%cflfile=monitor(this%fs%cfg%amRoot,'cfl_atom')
         call this%cflfile%add_column(this%time%n,'Timestep number')
         call this%cflfile%add_column(this%time%t,'Time')
         call this%cflfile%add_column(this%fs%CFLst,'STension CFL')
         call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
         call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
         call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
         call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
         call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
         call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
         call this%cflfile%write()
         ! Create LJCF monitor
         this%ljcf_file=monitor(this%fs%cfg%amRoot,'ljcf')
         call this%ljcf_file%add_column(this%time%n,'Timestep number')
         call this%ljcf_file%add_column(this%time%t,'Time')
         call this%ljcf_file%add_column(this%liqVolInjected,'Liq Vol Injected')
         call this%ljcf_file%add_column(this%InjectionVelocity,'Injection Velocity')
         call this%ljcf_file%write()
      end block create_monitor
      
      
      ! Create a timing monitor
      create_timing: block
         ! Create timers
         this%tstep =timer(comm=this%cfg%comm,name='Timestep')
         this%tvof  =timer(comm=this%cfg%comm,name='VOFsolve')
         this%tvel  =timer(comm=this%cfg%comm,name='Velocity')
         this%tpres =timer(comm=this%cfg%comm,name='Pressure')
         ! Create corresponding monitor file
         this%timefile=monitor(this%fs%cfg%amRoot,'timing')
         call this%timefile%add_column(this%time%n,'Timestep number')
         call this%timefile%add_column(this%time%t,'Time')
         call this%timefile%add_column(this%tstep%time ,trim(this%tstep%name))
         call this%timefile%add_column(this%tvof%time  ,trim(this%tvof%name))
         call this%timefile%add_column(this%tvel%time  ,trim(this%tvel%name))
         call this%timefile%add_column(this%tpres%time ,trim(this%tpres%name))
      end block create_timing

      create_merge_split: block 
         integer :: i,ierr,iunit,nchar
         logical :: file_exists
         character(len=5)  :: tmpchar
         character(len=60) :: tmpchar2
         character(len=22) :: buffer
         real(WP) :: current_time
         logical :: LID_silo_exists
         integer :: nAlloc_silo, ms_nout_time 
         real(WP),dimension(:),allocatable :: ms_times,time_silo
         real(WP) :: update_time
         !character(len=str_medium) :: gtda_output_type ! Default or Marching cubes 
         integer,dimension(:),allocatable :: LID_silo  

         LID_silo_exists = .false.

         !call create_timing('gtda')

         !call param_read('Silo update frequency',update_time,max_dt*10)
         !call param_read('Gtda output type', gtda_output_type,'Marching cubes')

         if (this%cfg%amRoot) then
            ! Check if the file exists
            INQUIRE(FILE="merge_split.csv", EXIST=file_exists)
            print*, "Checking for merge_split.csv file for GTDA output. File exists: ", file_exists
            if (.not.file_exists) then
               ! Create a new file with headers if it doesn't exist
               open(newunit=iunit, file="merge_split.csv", form="formatted", status="replace", action="write")
               write(iunit, "(A)") "Event Count, Event Type, Old IDs, New ID, Time, New Vol, X, Y, Z, U, V, W, L1, L2, L3"
               close(iunit)
            end if
            ! Check for gtda folder
            inquire(file="gtda/events",exist=file_exists)
            if (file_exists) then
               ! Get number of lines in file
               open(newunit=iunit,file="gtda/events",form="formatted",iostat=ierr,status='old')
               ms_nout_time=0
               do
                  read(iunit,*,end=1)
                  ms_nout_time=ms_nout_time+1
               end do
         1       continue
               close(iunit)
               ! Read the file and keep times less than current time
               allocate(ms_times(ms_nout_time))
               open(newunit=iunit,file="gtda/events",form="formatted",iostat=ierr,status='old')
               ! Get current time (formatted correctly)
               write(buffer,'(ES12.5)') this%time%t-this%time%dt*1e-10_WP
               read(buffer,*) this%time%t
               !read(buffer,*) current_time
               do i=1,ms_nout_time
                  ! Read file
                  read(iunit,'(5A,60A)') tmpchar,tmpchar2
                  ! Extract time from string
                  nchar=len_trim(tmpchar2)
                  read(tmpchar2(1:nchar-11),*) ms_times(i)
                  ! Check if it is in the future and exit if true
                  !if (ms_times(i).ge.current_time) then
                  if (ms_times(i).ge.this%time%t) then
                     ms_nout_time=i-1
                     exit
                  end if
               end do
               close(iunit)
               ! Write new file with only past times
               open(newunit=iunit,file="gtda/events",form="formatted",iostat=ierr,status='replace')
               do i=1,ms_nout_time
                  write(tmpchar2,'(ES12.5)') ms_times(i)
                  tmpchar2='time_'//trim(adjustl(tmpchar2))//'/Visit.silo'
                  write(iunit,'(A)') trim(adjustl(tmpchar2))
               end do
               deallocate(ms_times)
               close(iunit)
            else
               print*, "GTDA output folder not found, creating folder and skipping event read/write"
               call execute_command_line('mkdir -p gtda')
            end if

         end if

         
      end block create_merge_split

      ! Initialize an event for drop size analysis
      drop_analysis: block
         use param, only: param_read
         this%drop_evt=event(time=this%time,name='Drop analysis')
         call param_read('Drop analysis period',this%drop_evt%tper)
         if (this%drop_evt%occurs()) call analyse_drops(this)
      end block drop_analysis
      
   contains
      
      
      !> Function that localizes the x- boundary
      function xm_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imin) isIn=.true.
      end function xm_locator
      
      
      !> Function that localizes the x+ boundary
      function xp_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.eq.pg%imax+1) isIn=.true.
      end function xp_locator
      
      
      !> Function that localizes region of VOF removal
      function vof_removal_layer_locator(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (i.ge.pg%imax-this%nlayer) isIn=.true.
      end function vof_removal_layer_locator
      
      
      !> Function that defines a level set function for a half droplet
      function levelset_halfdrop(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=0.5_WP*this%djet-sqrt(xyz(1)**2+(xyz(2)-this%cfg%y(this%cfg%jmin))**2+xyz(3)**2)
      end function levelset_halfdrop

      !> Function that localizes the jet(s) initial location
      function jet(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         real(WP), dimension(3) :: xyz
         logical :: isIn
         isIn=.false.
         xyz(1)=pg%xm(i); xyz(2)=pg%ym(j); xyz(3)=pg%zm(k)
         if (levelset_halfdrop(xyz,0.0_WP).gt.0.0_WP) isIn=.true.
      end function jet
      
      !> Function that localizes the walls surrounding the jets
      function wall(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         logical :: isIn
         isIn=.false.
         if (j.le.pg%jmin-1+this%nwall.and.(.not.jet(pg,i,j,k))) isIn=.true.
      end function wall
      
      !> Function that localizes the jet(s) BCs at edge of domain
      function jet_bdy(pg,i,j,k) result(isIn)
         use pgrid_class, only: pgrid
         implicit none
         class(pgrid), intent(in) :: pg
         integer, intent(in) :: i,j,k
         real(WP), dimension(3) :: xyz
         logical :: isIn
         isIn=.false.
         xyz(1)=pg%xm(i); xyz(2)=pg%y(j); xyz(3)=pg%zm(k)
         if (j.eq.pg%jmin.and.jet(pg,i,j,k)) isIn=.true.
      end function jet_bdy
      
      !> Function that identifies liquid cells
      logical function label_liquid(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (this%vf%VF(i,j,k).gt.0.00001_WP) then
            label_liquid=.true.
         else
            label_liquid=.false.
         end if
      end function label_liquid

      function levelset_sphere(xyz,t) result(G)
         implicit none
         real(WP), dimension(3),intent(in) :: xyz
         real(WP), intent(in) :: t
         real(WP) :: G
         G=radius-sqrt(sum((xyz-center)**2))
      end function levelset_sphere

   end subroutine init
   
   
   !> Take one time step
   subroutine step(this)
      use tpns_class, only: arithmetic_visc
      implicit none
      class(ljcf), intent(inout) :: this
      
      ! Reset all timers and start timestep timer
      call this%tstep%reset()
      call this%tvof%reset()
      call this%tvel%reset()
      call this%tpres%reset()
      call this%tstep%start()
      
      ! Increment time
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()

      ! Apply jet velocity
      apply_bc: block
         use tpns_class, only: bcond
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         type(bcond), pointer :: mybc
         real(WP) :: liqVolInjected_dt
         integer :: n,i,j,k
         ! Compute injection velocity
         if (this%liqVolInjected .lt. this%liqVol) then
            this%InjectionVelocity=this%gravity*this%time%t  ! Velocity increases linearly with time
         else
            this%InjectionVelocity=0.0_WP                    ! Velocity stops once volume is reached
         end if
         ! Apply injection velocity to the jet boundary condition 
         call this%fs%get_bcond('jet',mybc)
         liqVolInjected_dt = 0.0_WP
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            
            this%fs%V(i,j,k) = this%InjectionVelocity
            liqVolInjected_dt = liqVolInjected_dt + this%fs%V(i,j,k)*this%vf%VF(i,j-1,k)*this%cfg%dx(i)*this%cfg%dz(k)*this%time%dt
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,liqVolInjected_dt,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         this%liqVolInjected = this%liqVolInjected + liqVolInjected_dt
      end block apply_bc

      ! Remember old VOF
      this%vf%VFold=this%vf%VF

      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Prepare old sflaggered density (at n)
      call this%fs%get_olddensity(vf=this%vf)

      ! VOF solver step
      call this%tvof%start() ! Start VOF timer
      call this%vf%advance(dt=this%time%dt,U=this%fs%U,V=this%fs%V,W=this%fs%W)
      call this%tvof%stop() ! Stop VOF timer

      ! Advance stracker
      call this%strack%advance(make_label=label_liquid)
      call analyze_merge_split(this)
      
      ! Prepare new sflaggered viscosity (at n+1)
      call this%fs%get_viscosity(vf=this%vf,strat=arithmetic_visc)
      
      ! Perform sub-iterations
      do while (this%time%it.le.this%time%itmax)
         
         ! Start velocity timer
         call this%tvel%start()
         
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)

         ! Preliminary mass and momentum transport step at the interface
         call this%fs%prepare_advection_upwind(dt=this%time%dt)

         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)

         ! Add momentum source terms
         call this%fs%addsrc_gravity(this%resU,this%resV,this%resW)
         
         ! Assemble explicit residual
         this%resU=-2.0_WP*this%fs%rho_U*this%fs%U+(this%fs%rho_Uold+this%fs%rho_U)*this%fs%Uold+this%time%dt*this%resU
         this%resV=-2.0_WP*this%fs%rho_V*this%fs%V+(this%fs%rho_Vold+this%fs%rho_V)*this%fs%Vold+this%time%dt*this%resV
         this%resW=-2.0_WP*this%fs%rho_W*this%fs%W+(this%fs%rho_Wold+this%fs%rho_W)*this%fs%Wold+this%time%dt*this%resW   
         
         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)

         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Apply boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Stop velocity timer and start pressure timer
         call this%tvel%stop()
         call this%tpres%start()
         
         ! Solve Poisson equation
         call this%fs%update_laplacian()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         call this%fs%add_surface_tension_jump(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         ! call this%fs%add_surface_tension_jump_twoVF(dt=this%time%dt,div=this%fs%div,vf=this%vf)
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/max(epsilon(0.0_WP),this%fs%rho_U)
         this%fs%V=this%fs%V-this%time%dt*this%resV/max(epsilon(0.0_WP),this%fs%rho_V)
         this%fs%W=this%fs%W-this%time%dt*this%resW/max(epsilon(0.0_WP),this%fs%rho_W)
         
         ! Apply boundary conditions
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Stop pressure timer
         call this%tpres%stop()
         
         ! Increment sub-iteration counter
         this%time%it=this%time%it+1
         
      end do
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()
      
      ! Remove VOF at edge of domain
      remove_vof: block
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
         use parallel, only: MPI_REAL_WP
         integer :: n,i,j,k,ierr
         this%vof_removed=0.0_WP
         do n=1,this%vof_removal_layer%no_
            i=this%vof_removal_layer%map(1,n)
            j=this%vof_removal_layer%map(2,n)
            k=this%vof_removal_layer%map(3,n)
            if (n.le.this%vof_removal_layer%n_) this%vof_removed=this%vof_removed+this%cfg%vol(i,j,k)*this%vf%VF(i,j,k)
            this%vf%VF(i,j,k)=0.0_WP
         end do
         call MPI_ALLREDUCE(MPI_IN_PLACE,this%vof_removed,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr)
         call this%vf%clean_irl_and_band()
         !call silo_lists_update(4)
      end block remove_vof
      
      ! Output to ensight
      if (this%ens_evt%occurs()) then
         ! Update surface mesh
         update_smesh: block
            use irl_fortran_interface, only: getNumberOfPlanes,getNumberOfVertices
            integer :: i,j,k,np,nplane
            ! Transfer polygons to smesh
            call this%vf%update_surfmesh(this%smesh)
            ! Also populate nplane variable
            this%smesh%var(1,:)=1.0_WP
            np=0
            do k=this%vf%cfg%kmin_,this%vf%cfg%kmax_
               do j=this%vf%cfg%jmin_,this%vf%cfg%jmax_
                  do i=this%vf%cfg%imin_,this%vf%cfg%imax_
                     if (this%cfg%VF(i,j,k).lt.2.0_WP*epsilon(1.0_WP)) cycle ! Skip cells below VF threshold
                     do nplane=1,getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k))
                        if (getNumberOfVertices(this%vf%interface_polygon(nplane,i,j,k)).gt.0) then
                           np=np+1; this%smesh%var(1,np)=real(getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)),WP)
                           this%smesh%var(2,np)=0.0_WP; this%smesh%var(3,np)=real(this%strack%id(i,j,k),WP)
                        end if
                     end do
                  end do
               end do
            end do
         end block update_smesh
         call this%ens_out%write_data(this%time%t)
      end if
      
      ! Stop timestep timer
      call this%tstep%stop()

      ! Analyse droplets
         if (this%drop_evt%occurs()) call analyse_drops(this)
      
      ! Perform and output monitoring
      call this%fs%get_max()
      call this%vf%get_max()
      call this%mfile%write()
      call this%cflfile%write()
      call this%timefile%write()
      call this%ljcf_file%write()
      
      ! Finally, see if it's time to save restart files
      if (this%save_evt%occurs()) then
         save_restart: block
            use irl_fortran_interface
            use string, only: str_medium
            character(len=str_medium) :: timestamp
            real(WP), dimension(:,:,:), allocatable :: P11,P12,P13,P14
            real(WP), dimension(:,:,:), allocatable :: P21,P22,P23,P24
            integer :: i,j,k
            real(WP), dimension(4) :: plane
            ! Handle IRL data
            allocate(P11(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P12(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P13(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P14(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P21(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P22(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P23(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            allocate(P24(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
            do k=this%vf%cfg%kmino_,this%vf%cfg%kmaxo_
               do j=this%vf%cfg%jmino_,this%vf%cfg%jmaxo_
                  do i=this%vf%cfg%imino_,this%vf%cfg%imaxo_
                     ! First plane
                     plane=getPlane(this%vf%liquid_gas_interface(i,j,k),0)
                     P11(i,j,k)=plane(1); P12(i,j,k)=plane(2); P13(i,j,k)=plane(3); P14(i,j,k)=plane(4)
                     ! Second plane
                     plane=0.0_WP
                     if (getNumberOfPlanes(this%vf%liquid_gas_interface(i,j,k)).eq.2) plane=getPlane(this%vf%liquid_gas_interface(i,j,k),1)
                     P21(i,j,k)=plane(1); P22(i,j,k)=plane(2); P23(i,j,k)=plane(3); P24(i,j,k)=plane(4)
                  end do
               end do
            end do
            ! Prefix for files
            write(timestamp,'(es12.5)') this%time%t
            ! Populate df and write it
            call this%df%push(name='t'  ,val=this%time%t )
            call this%df%push(name='dt' ,val=this%time%dt)
            call this%df%push(name='U'  ,var=this%fs%U   )
            call this%df%push(name='V'  ,var=this%fs%V   )
            call this%df%push(name='W'  ,var=this%fs%W   )
            call this%df%push(name='P'  ,var=this%fs%P   )
            call this%df%push(name='Pjx',var=this%fs%Pjx )
            call this%df%push(name='Pjy',var=this%fs%Pjy )
            call this%df%push(name='Pjz',var=this%fs%Pjz )
            call this%df%push(name='P11',var=P11         )
            call this%df%push(name='P12',var=P12         )
            call this%df%push(name='P13',var=P13         )
            call this%df%push(name='P14',var=P14         )
            call this%df%push(name='P21',var=P21         )
            call this%df%push(name='P22',var=P22         )
            call this%df%push(name='P23',var=P23         )
            call this%df%push(name='P24',var=P24         )
            call this%df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
            ! Deallocate
            deallocate(P11,P12,P13,P14,P21,P22,P23,P24)
         end block save_restart
      end if

   contains

      !> Function that identifies liquid cells
      logical function label_liquid(i,j,k)
         implicit none
         integer, intent(in) :: i,j,k
         if (this%vf%VF(i,j,k).gt.0.00001_WP) then
            label_liquid=.true.
         else
            label_liquid=.false.
         end if
      end function label_liquid
      
   end subroutine step
   
   
   !> Finalize nozzle simulation
   subroutine final(this)
      implicit none
      class(ljcf), intent(inout) :: this
      
      ! Deallocate work arrays
      deallocate(this%resU,this%resV,this%resW,this%Ui,this%Vi,this%Wi)
      
   end subroutine final
   
   
end module ljcf_class