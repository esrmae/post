module read_post
!	USE shared_data
  use shared_post
  use cal_coefficient

  implicit none 

!*******************************************************************************
! We assume that this code is in a postprocessing folder, which is in the
! same folder as a number of directories containing the codes' output files
! which are labelled run1,run2,...,run(nruns) 
!*******************************************************************************

  contains

!*******************************************************************************

subroutine read_setup
	
  call default_values_main
  call default_values_post
  call input_parameters

  return
end subroutine read_setup

!*******************************************************************************

subroutine default_values_post

! Post-Processing Parameters
  !nruns=1; nstart=1; retau_av=1
  isteady=1; small=1.0e-6; height=2; itec=1
  rho=1.0e0; tstart0=0; jspec=5; read_3d=0; iflucs=1; itke=1; icorr=0
  istreak=0; ivort3d=0; nsims=1; sim_print=0; tpts=100; tmin=0; tmax=0; 
  tcplt=1; re_print=0; y_print=0; y_half=0; ibm=0; ny_solid=0; half_av=0
  sim3d=1; jvis=20; tvis=1; uns1d=0; restrt=0; tkemn=1; spec_loc=0
  ntsamps=21; reinit=0; tperiod=0.0e0; jdec=10; tperiod_plus=0.0e0
  omega_plus=0.0e0; re_val=3150.0e0; retau_val=200.0e0
  wavelen=0.0e0; wavelen_plus=0.0e0; k_plus=0.0e0; kx_plus=0.0e0; kz_plus=0.0e0
  istruct=0;itop=0; ibot=1;str_num=2000;iwxfd=0;corr_thd=0.4e0;strl_thd=150.0e0; strd_thd=10.0e0;iterstr=1
  iangle=0;iraw=0;iplot3d=0;iplot2d=0;iplot1d=0;iconstruc=0;islicexyz=10
  idecom=0;isuper=0;irlarge=0;icontr=0;ispec3d=0;kx_large=0;kz_large=0;idef_sup=1
  jmodu=10
  itrans=0; theta=0.0e0
  return
end subroutine default_values_post

!*****************************************************************************

subroutine default_values_main  
! This procedure sets all default values as in the main code
  implicit none

! Grid
  nx=64;ny=129;nz=64          
      
! Main simulation parameters
  ichann=1; ibmpfg = 0; restrt=1; roh=1.0e0; ngr=3; nbl = 1; mpi=0

! Boundary conditons parameters 
  iconbc = 0;jconbc=1;kconbc = 0;ibwest = 0;ibeast = 0;ibnoth=1;ibsoth=1         
  kbsoth=1;kbnoth=1;length=7.0e0;width=3.5e0;height=2.0e0;idomain = 0  

! Other simulation parameters
  pi=2.0e0*ASIN(1.0e0);pi2=pi*2.0e0;pi4=pi*4.0e0
  ifftw = 0; large = 1.0e6; small = 1.0e-6; very_small = 1.0e-12
                 
! Grid parameters
  cony=2.04e0

  return
end subroutine default_values_main
 
!*******************************************************************************

subroutine input_parameters

  integer :: ncheck,i,j
  character :: dummy*70,file1*40,startno*20,filedum*30,simno*20,dum6*6,dum7*7
  logical :: equals
  character(LEN=30), dimension(1000) :: charac
  real(mytype), dimension (1000) :: var

! Read the info from input-stat.par
  open(1, File='input-stat.par', Status='OLD')
  read (1,'(A70)') dummy
  read (1,'(A70)') dummy
  read (1,'(A70)') dummy
  read (1,'(A70)') dummy
  read (1,'(A70)') dummy
  read (1,'(A60)') unit1
  read (1,'(A60)') file2davg
  Write(*,*) 'reading... input-stat.par'
  ncheck=0
  do i = 1,1000
    read(1,*,end=10) charac(i),var(i)
    ncheck=ncheck+1

!   Inputed parameters
    if (charac(i) == 'nsims') then 
      nsims=int(var(i)) 
      allocate(nruns(nsims),nstart(nsims))
      nruns=1;nstart=1
    else if (charac(i) == 'nruns') then 
      allocate(nruns(1),nstart(1))
      nruns=1;nstart=1
      nruns(1)=int(var(i)) 
    else if (charac(i) == 'nstart')       then; nstart(1)=int(var(i))
    else if (charac(i) == 'isteady')      then; isteady=int(var(i))
    else if (charac(i) == 'tstart0')      then; tstart0=int(var(i))
    else if (charac(i) == 'sim_print')    then; sim_print=int(var(i))
    else if (charac(i) == 'if_les')       then; if_les=int(var(i))
    else if (charac(i) == 'tpts')         then; tpts=int(var(i))
    else if (charac(i) == 'tmin')         then; tmin=var(i)
    else if (charac(i) == 'tmax')         then; tmax=var(i)
    else if (charac(i) == 're_print')     then; re_print=int(var(i))
    else if (charac(i) == 'y_print')      then; y_print=int(var(i))
    else if (charac(i) == 'y_half')       then; y_half=int(var(i))
    else if (charac(i) == 'itec')         then; itec=int(var(i))
    else if (charac(i) == 'half_av')      then; half_av=int(var(i))
    else if (charac(i) == 'ibm')          then; ibm=int(var(i))
    else if (charac(i) == 'ny_solid')     then; ny_solid=int(var(i))
    else if (charac(i) == 'tcplt')        then; tcplt=int(var(i))
    else if (charac(i) == 'jspec')        then; jspec=int(var(i))
    else if (charac(i) == 'uns1d')        then; uns1d=int(var(i))
    else if (charac(i) == 'ntsamps')      then; ntsamps=int(var(i))
    else if (charac(i) == 'read_2d')      then; read_2d=int(var(i))
    else if (charac(i) == 'read_3d')      then; read_3d=int(var(i))
    else if (charac(i) == 'restrt')       then; restrt=int(var(i))
    else if (charac(i) == 'iflucs')       then; iflucs=int(var(i))
    else if (charac(i) == 'itke')         then; itke=int(var(i))
    else if (charac(i) == 'ivort')        then; ivort=int(var(i))
    else if (charac(i) == 'iquad')        then; iquad=int(var(i))
    else if (charac(i) == 'ishear')       then; ishear=int(var(i))
    else if (charac(i) == 'ianis')        then; ianis=int(var(i))
    else if (charac(i) == 'icorr')        then; icorr=int(var(i))
    else if (charac(i) == 'sim3d')        then; sim3d=int(var(i))
    else if (charac(i) == 'istreak')      then; istreak=int(var(i))
    else if (charac(i) == 'ivort3d')      then; ivort3d=int(var(i))
    else if (charac(i) == 'ilamq')        then; ilamq=int(var(i))
    else if (charac(i) == 'jvis')         then; jvis=int(var(i))
    else if (charac(i) == 'jdec')         then; jdec=int(var(i))
    else if (charac(i) == 'tvis')         then; tvis=int(var(i))
    else if (charac(i) == 'if_les')       then; if_les=int(var(i))
    else if (charac(i) == 'tkemn')        then; tkemn=int(var(i))
    else if (charac(i) == 'spec_loc')     then; spec_loc=int(var(i))
    else if (charac(i) == 'reinit')       then; reinit=var(i)
    else if (charac(i) == 'tperiod')      then; tperiod=var(i)
    else if (charac(i) == 'hybrid')       then; hybrid=int(var(i))
    else if (charac(i) == 'mpi')          then; mpi=int(var(i))
    else if (charac(i) == 're_val')       then; re_val=var(i)
    else if (charac(i) == 'retau_val')    then; retau_val=var(i)
    else if (charac(i) == 'tperiod_plus') then; tperiod_plus=var(i)
    else if (charac(i) == 'omega_plus')   then; omega_plus=var(i)
    else if (charac(i) == 'k_plus')   	  then; k_plus=var(i)
    else if (charac(i) == 'kx_plus')      then; kx_plus=var(i)
    else if (charac(i) == 'kz_plus')      then; kz_plus=var(i)
    else if (charac(i) == 'itrans')       then; itrans=int(var(i))
    else if (charac(i) == 'theta')        then; theta=var(i)
    else if (charac(i) == 'iangle')       then; iangle=var(i)
    else if (charac(i) == 'istruct')      then; istruct=var(i)
    else if (charac(i) == 'iensem')       then; iensem=var(i)
    else if (charac(i) == 'itop')         then; itop=var(i)
    else if (charac(i) == 'ibot')         then; ibot=var(i)
    else if (charac(i) == 'ifile1')       then; ifile1=var(i)
    else if (charac(i) == 'ifile2')       then; ifile2=var(i)
    else if (charac(i) == 'xw')           then; xw=int(var(i))
    else if (charac(i) == 'yw')           then; yw=int(var(i))
    else if (charac(i) == 'zw')           then; zw=int(var(i))
    else if (charac(i) == 'ifilt')        then; ifilt=int(var(i))
    else if (charac(i) == 'iterstr')      then; iterstr=int(var(i))
    else if (charac(i) == 'iwrite')       then; iwrite=int(var(i))
    else if (charac(i) == 'cnt_stop')     then; cnt_stop=int(var(i))
    else if (charac(i) == 'corr_thd')     then; corr_thd=var(i)
    else if (charac(i) == 'str_num')      then; str_num=int(var(i))
    else if (charac(i) == 'strl_thd')     then; strl_thd=var(i)
    else if (charac(i) == 'strd_thd')     then; strd_thd=var(i)
    else if (charac(i) == 'iosci')        then; iosci=int(var(i))
    else if (charac(i) == 'phase')        then; phase=int(var(i))
    else if (charac(i) == 'iwxfd')        then; iwxfd=int(var(i))
    else if (charac(i) == 'iraw')         then; iraw=int(var(i))
    else if (charac(i) == 'iplot3d')      then; iplot3d=int(var(i))
    else if (charac(i) == 'iplot2d')      then; iplot2d=int(var(i))
    else if (charac(i) == 'iplot1d')      then; iplot1d=int(var(i))
    else if (charac(i) == 'iconstruc')    then; iconstruc=int(var(i))
    else if (charac(i) == 'islicexyz')    then; islicexyz=int(var(i))
    else if (charac(i) == 'idecom')    	  then; idecom=int(var(i))
    else if (charac(i) == 'isuper')    	  then; isuper=int(var(i))
    else if (charac(i) == 'irlarge')      then; irlarge=int(var(i))
    else if (charac(i) == 'kx_large')     then; kx_large=int(var(i))
    else if (charac(i) == 'kz_large')     then; kz_large=int(var(i))
    else if (charac(i) == 'idef_sup')     then; idef_sup=int(var(i))
    else if (charac(i) == 'icontr')       then; icontr=int(var(i))
    else if (charac(i) == 'ispec3d')      then; ispec3d=int(var(i))
    else
      do j=1,nsims
        write(simno,*) j
        simno=adjustl(simno)
        if (charac(i) == 'nruns'//trim(simno)) then; nruns(j)=int(var(i)); goto 5
        else if (charac(i) == 'nstart'//trim(simno)) then; nstart(j)=int(var(i)); goto 5
        end if
      end do
      print*, 'Variable name is not correct in input-stat'
 5  end if
    write(*,*) charac(i),var(i)
  end do
 10  close(1)
  if (restrt == 1) restrt_file = trim(dummy)
  unit1=trim(unit1)
  Write(*,*) '################',dummy,unit1

!  if (hybrid == 1) isteady = 0

  if (nsims == 1) then; isim1 = 1
  else; isim1 = 0
  end if

  if(read_2d == 1) read_3d = 1

  if (tperiod==0.0e0) then
    if (omega_plus .ne. 0.0e0) then
      tperiod_plus=2.0e0*pi/omega_plus
      tperiod=tperiod_plus*re_val/(retau_val**2.0e0)
    end if
  end if
  if (wavelen==0.0e0) then
    if(k_plus .eq. 0.0e0 .and. kx_plus .ne. 0.0e0) k_plus=kx_plus
    if(k_plus .eq. 0.0e0 .and. kz_plus .ne. 0.0e0) k_plus=kz_plus
    if(kx_plus .ne. 0.0e0 .and. kz_plus .ne. 0.0e0) then
      write(*,*) 'oblique traveling wave! No 2D statistics! Changing isteady to 1 ...'
      isteady=1
    end if
    if (k_plus .ne. 0.0e0) then
      wavelen_plus=2.0e0*pi/k_plus
      wavelen=wavelen_plus/retau_val
    end if
  end if
!
  if(itrans==1) theta=theta*pi/180.0e0
! 
  return
end subroutine input_parameters

!*******************************************************************************

subroutine read_vals

  call stat_parameters	
  call allocate_arrays_post
  call grid_setup
  call read_0d_data
  call compute_fact
  call avg_sims_setup
  if (hybrid == 1) call calc_0d_avs

  if (isteady == 1) then
    call calc_0d_avs
    if (lpstat_1d == 1) call read_1d_avg
  else if (isteady == 0) then
    if (lp_ins_1d == 1) call read_ins_1d
    call calc_unst
  else if (isteady == 2) then
    if (lpstat_2d == 1) call read_2d_avg
    if (lpstat_2d == 1) write(*,*) 'reading_2d_avg'
    call calc_unst    
  end if

  return
end subroutine read_vals

!*******************************************************************************

subroutine stat_parameters

  integer :: ncheck,i,j,max_stats
  character :: dummy*70,file1*40,startno*20,filedum*30
  character(LEN=30), dimension(1000) :: charac
  Real(mytype), dimension (1000) :: var

  write(startno,*) nstart(1)
  startno=adjustl(startno)

! Read the info from stat-deck.par in run1 folder
  filedum='../run'//TRIM(startno)//'/stat-deck.par'
  if (isim1.ne.1) filedum='../sim1/run'//trim(startno)//'/stat-deck.par'

  open(2, file=filedum, status='OLD')
  ncheck=0
  do i = 1,1000
    read(2,*,end=11) charac(i),var(i)
    ncheck=ncheck+1
  end do
 11 close(2)

  do i = 1,ncheck
! Inputed parameters
    if (charac(i) == 'nx')                then; nx=int(var(i)) 
    else if (charac(i) == 'ny')           then; ny=int(var(i))
    else if (charac(i) == 'nz')           then; nz=int(var(i))
    else if (charac(i) == 'nstep')        then;  
    else if (charac(i) == 'itmavg')       then; itmavg=int(var(i))
    else if (charac(i) == 'ntmavg')       then; ntmavg=int(var(i))
    else if (charac(i) == 'samples')      then; samples=int(var(i))
    else if (charac(i) == 'lpstat_1d')    then; lpstat_1d=int(var(i))
    else if (charac(i) == 'lpstat_2d')    then; lpstat_2d=int(var(i))
    else if (charac(i) == 'lpstat_x')     then; lpstat_x=int(var(i))
    else if (charac(i) == 'lpstat_z')     then; lpstat_z=int(var(i))
    else if (charac(i) == 'ibin_rite')    then; ibin_rite=int(var(i))
    else if (charac(i) == 'lp_ins_1d')    then; lp_ins_1d=int(var(i))
    else if (charac(i) == 'lp_ins_2d')    then; lp_ins_2d=int(var(i))
    else if (charac(i) == 'lp_ins_3d')    then; lp_ins_3d=int(var(i))
    else if (charac(i) == 'lp_snap_1d')   then; lp_snap_1d=int(var(i))
    else if (charac(i) == 'lp_snap_2d')   then; lp_snap_2d=int(var(i))
    else if (charac(i) == 'lp_snap_3d')   then; lp_snap_3d=int(var(i))
    else if (charac(i) == 'lp_snap_x')    then; lp_snap_x=int(var(i))
    else if (charac(i) == 'lp_snap_y')    then; lp_snap_y=int(var(i))
    else if (charac(i) == 'lp_snap_z')    then; lp_snap_z=int(var(i))
    else if (charac(i) == 'lp_snap_xy')   then; lp_snap_xy=int(var(i))
    else if (charac(i) == 'lp_snap_xz')   then; lp_snap_xz=int(var(i))
    else if (charac(i) == 'lp_snap_yz')   then; lp_snap_yz=int(var(i))
    else if (charac(i) == 'nsmpl_ins_1d') then; nsmpl_ins_1d=int(var(i))
    else if (charac(i) == 'nsmpl_ins_2d') then; nsmpl_ins_2d=int(var(i))
    else if (charac(i) == 'nsmpl_ins_3d') then; nsmpl_ins_3d=int(var(i))
    else if (charac(i) == 'nsmpl_snap_1d')then; nsmpl_snap_1d=int(var(i))
    else if (charac(i) == 'nsmpl_snap_2d')then; nsmpl_snap_2d=int(var(i))
    else if (charac(i) == 'nsmpl_snap_3d')then; nsmpl_snap_3d=int(var(i))
    else if (charac(i) == 'cony')         then; cony=var(i)
    else if (charac(i) == 'length')       then; length=var(i)
    else if (charac(i) == 'width')        then; width=var(i)
    else if (charac(i) == 'height')       then; height=var(i)
    else if (charac(i) == 'nu')           then; nu=var(i)
    else if (charac(i) == 'rho')          then; rho=var(i)
    else if (charac(i) == 'total_step')   then; 
    else; print*, 'Variable name is not correct in stat-deck (run1)'
    end if
    write(*,*) charac(i),var(i)
  end do

  nxt=nx+1;nyt=ny+1;nzt=nz+1
  if(lpstat_x == 1) nxoz=nz
  if(lpstat_z == 1) nxoz=nx
  nxozt=nxoz+1

! Read nstep info from all stat-decks
  allocate(nstats(nsims))
  do i=1,nsims
    nstats(i) = nruns(i) - nstart(i) + 1
  end do
  max_stats = maxm(nstats,nsims)
  allocate(nsteps(max_stats,nsims),nsample(max_stats,nsims),ins_1d(max_stats,nsims), &
           ins_2d(max_stats,nsims),ins_3d(max_stats,nsims),snap_1d(max_stats,nsims), &
           snap_2d(max_stats,nsims),snap_3d(max_stats,nsims))
  allocate(step_tot(nsims),samp_tot(nsims),ins_1d_tot(nsims),ins_2d_tot(nsims), &
           ins_3d_tot(nsims),snap_1d_tot(nsims),snap_2d_tot(nsims),snap_3d_tot(nsims), &
	   tot_v(nsims),nmin(nsims),nmax(nsims))

! Initialise the vectors
  nsteps = 0; nsample = 0; ins_1d = 0; ins_2d = 0; ins_3d = 0; snap_1d = 0; snap_2d = 0;snap_3d = 0
  step_tot = 0; samp_tot = 0; ins_1d_tot = 0; ins_2d_tot = 0; ins_3d_tot = 0; snap_1d_tot = 0
  snap_2d_tot = 0; snap_3d_tot = 0

  do simit=1,nsims
    do i = nstart(simit), nruns(simit)
      file1 = 'stat-deck.par'
      call change_dir(file1,i)
      write(*,*) 'reading...',file1
      open(3, File=file1, Status='OLD')
      do j = 1,1000
        read(3,*,end=12) charac(j),var(j)
        if (charac(j) == 'nstep')              then; nsteps(i-nstart+1,simit)=int(var(j))
        else if (charac(j) == 'samples')       then; nsample(i-nstart+1,simit)=int(var(j))
        else if (charac(j) == 'nsmpl_ins_1d')  then; ins_1d(i-nstart+1,simit)=int(var(j))
        else if (charac(j) == 'nsmpl_ins_2d')  then; ins_2d(i-nstart+1,simit)=int(var(j))
        else if (charac(j) == 'nsmpl_ins_3d')  then; ins_3d(i-nstart+1,simit)=int(var(j))
        else if (charac(j) == 'nsmpl_snap_1d') then; snap_1d(i-nstart+1,simit)=int(var(j))
        else if (charac(j) == 'nsmpl_snap_2d') then; snap_2d(i-nstart+1,simit)=int(var(j))
        else if (charac(j) == 'nsmpl_snap_3d') then; snap_3d(i-nstart+1,simit)=int(var(j))
        else if (charac(i) == 'ntmavg')        then; ntmavg=int(var(i))
        end if
      end do
 12   close(3)
    end do

!  Finds the total number of steps/samples
    do i = 1, nstats(simit)
      step_tot(simit) = step_tot(simit) + nsteps(i,simit)
      samp_tot(simit) = samp_tot(simit) + nsample(i,simit)
      ins_1d_tot(simit) = ins_1d_tot(simit) + ins_1d(i,simit)
      ins_2d_tot(simit) = ins_2d_tot(simit) + ins_2d(i,simit)
      ins_3d_tot(simit) = ins_3d_tot(simit) + ins_3d(i,simit)
      snap_1d_tot(simit) = snap_1d_tot(simit) + snap_1d(i,simit)
      snap_2d_tot(simit) = snap_2d_tot(simit) + snap_2d(i,simit)
      snap_3d_tot(simit) = snap_3d_tot(simit) + snap_3d(i,simit)
    end do
    tot_v(simit)=step_tot(simit)/ntmavg 
  end do

  return
end subroutine stat_parameters

!*******************************************************************************

subroutine allocate_arrays_post

! Max values for array allocation
  max_tot = maxm(step_tot,nsims)
  max_v = maxm(samp_tot,nsims)
  snap_1d_max = maxm(snap_1d_tot,nsims)
  snap_2d_max = maxm(snap_2d_tot,nsims)
  snap_3d_max = maxm(snap_3d_tot,nsims)
  ins_1d_max = maxm(ins_1d_tot,nsims)
  if (nx.ge.nz) then; max_xz=nx
  else; max_xz=nz
  end if

! Grid arrays
  allocate(x(nxt),y(nyt),z(nzt),dx(-1:nx+2),dy(-1:ny+2),dz(-1:nz+2),              &
           xc(0:nxt),yc(0:nyt),zc(0:nzt),dxs(nxt),dys(nyt),dzs(nzt),xplus(nxt),   &
           yplus(nyt),zplus(nzt),xcplus(0:nxt),ycplus(0:nyt),zcplus(0:nzt))
! 0d arrays for every timestep
  allocate(time(max_tot,nsims),dt(max_tot,nsims),utau(max_tot,nsims),utau1(max_tot,nsims), &
           utau2(max_tot,nsims),dpdx_mean(max_tot,nsims),u_bulk(max_tot,nsims),            &
           retau(max_tot,nsims),re(max_tot,nsims))
  time=0.0e0

! 0d arrays for every ntmavg timesteps
  allocate(time_v(max_v,nsims),dt_v(max_v,nsims),utau_v(max_v,nsims),u_bulk_v(max_v,nsims), &
           esum_v(max_v,nsims),energy_v(max_v,nsims),retau_v(max_v,nsims),re_v(max_v,nsims))

! 1d avg arrays
  if (lpstat_1d == 1) then
    allocate(stat_avg(-1:ny+2,4,4),flucs(-1:ny+2,4,4))
    if (iquad == 1)  allocate(quad_avg(-1:ny+2,9))
    if (ivort == 1)  allocate(omega_avg(-1:ny+2,6),omegamn_avg(-1:ny+1,3),wx_avg(-1:ny+2,2,8))
    if (itke == 1)   allocate(tke_avg(-1:ny+2,4,7),kinen(-1:ny+2,4,6))
    if (ishear == 1) allocate(shear(-1:ny+2,4))
    if (ianis == 1)  allocate(anisb(-1:ny+2,4),anisII(-1:ny+2),anisIII(-1:ny+2), &
                              anisf(-1:ny+2),anisg(-1:ny+2))
    if (iangle ==1) allocate(angle(-1:ny+2,5))
  end if

  ALLOCATE(utau_av(nsims),utau1_av(nsims),utau2_av(nsims),retau_av(nsims),re_av(nsims), &
           ub_av(nsims),Cf(nsims))

! Gradient cofficients arrays
  allocate(gcdx_f(nxt,4,ngr,nbl),gcdy_f(nyt,4,ngr,nbl),gcdz_f(nzt,4,ngr,nbl),   &
           gcdx_c(0:nxt,5,ngr,nbl),gcdy_c(0:nyt,5,ngr,nbl),gcdz_c(0:nzt,5,ngr,nbl),   &
           gfdx_c(0:nxt,4,ngr,nbl),gfdy_c(nyt,4,ngr,nbl),gfdz_c(nzt,4,ngr,nbl), &
           grdsxx(nx/2+1,5,2),grdsyy(ny/2+1,5,2),grdszz(nz/2+2,5,2),            &
           gcdz_f2(nzt,4, ngr,nbl),gfdz_c2(nzt,4,ngr,nbl))

  return
end subroutine allocate_arrays_post

!*******************************************************************************

subroutine grid_setup
! This procedure defines the computational grid for the old domain
  integer :: i,j,k

  do i = 1,nxt; x(i) = 1.0e0/real(nx,kind=mytype)*real(i-1,kind=mytype); end do
  do j = 1,nyt; y(j) = 1.0e0/real(ny,kind=mytype)*real(j-1,kind=mytype); end do
  do k = 1,nzt; z(k) = 1.0e0/real(nz,kind=mytype)*real(k-1,kind=mytype); end do

! Nonuniform grid for Y
  if (cony > small) then
    do j = 1,nyt
      y(j) = tanh(cony*(2.0e0*(j-1.0e0)/ny-1.0e0))/tanh(cony)
    end do
    y=(y+1.0e0)*0.5e0
  end if

  y(1) = 0.0e0
  y(nyt) = 1.0e0

  do i = 1,nx
    xc(i) = 0.5e0*(x(i)+x(i+1))
  end do
  xc(0) = -xc(1)
  xc(nxt) = x(nxt)+xc(1)

  do j = 1,ny
    yc(j) = 0.5e0*(y(j)+y(j+1))
  end do
  yc(0) = y(1)
  yc(nyt) = y(nyt)

  do k = 1,nz
    zc(k) = 0.5e0*(z(k)+z(k+1))
  end do
  zc(0) = -zc(1)
  zc(nzt) = z(nzt)+zc(1)

! Scale the grid to the correct size
  x=x*length
  y=y*height
  z=z*width
  xc=xc*length
  yc=yc*height
  zc=zc*width
!	
! Calculates the dx grid values
  do i = 1,nx
    dx(i) = x(i+1)-x(i)
  end do
  dx(-1) = dx(1)
  dx(0) = dx(1)
  dx(nx+1) = dx(nx)
  dx(nx+2) = dx(nx)

  do j = 1,ny
    dy(j) = y(j+1)-y(j)
  end do
  dy(-1) = dy(1)
  dy(0) = dy(1)
  dy(ny+1) = dy(ny)
  dy(ny+2) = dy(ny)

  do k = 1,nz
    dz(k) = z(k+1)-z(k)
  end do
  dz(-1) = dz(1)
  dz(0) = dz(1)
  dz(nz+1) = dz(nz)
  dz(nz+2) = dz(nz)

  do i = 1,nx+1
    dxs(i) = 0.5e0*(dx(i-1)+dx(i))
  end do
  do j = 1,ny+1
    dys(j) = 0.5e0*(dy(j-1)+dy(j))
  end do
  do k = 1,nz+1
   dzs(k) = 0.5e0*(dz(k-1)+dz(k))
  end do

  return
end subroutine grid_setup

!*******************************************************************************

subroutine read_0d_data
! This reads the 0d data from STAT_0D-1.dat and STAT_0D-2.dat
  character :: file2*40,file3*40
  integer :: ii,i,j,time_step
  integer :: tpfile

! We start with STAT_0D-1.dat which has the info from every time step as:
! ii,time,dt,utau,utau1,utau2,dp/dx_mean,U_bulk
  nsteps=0
  do simit=1,nsims	
    time_step=1
    do i = nstart(simit), nruns(simit)
      file2 = 'STAT_0D-1.dat'
      call change_dir(file2,i)
      write(*,*) 'reading...',file2
      open(4, file=file2, status='OLD')
      tpfile=0
      do j = 1,1000000
        read(4,*,end=13) ii,time(time_step,simit),dt(time_step,simit),utau(time_step,simit),&
                         utau1(time_step,simit),utau2(time_step,simit),dpdx_mean(time_step,simit),&
                         u_bulk(time_step,simit)
        time_step = time_step + 1
        tpfile = tpfile + 1
      end do
 13   close(4)
      nsteps(i-nstart(simit)+1,simit) = tpfile
    end do
    if (sum(nsteps(:,simit)).ne.step_tot(simit)) then
      write(*,*) 'Some runs were ended by the 48hr time limit!'
    end if
    step_tot(simit) = sum(nsteps(:,simit))
    write(*,*) 'step: ',simit, step_tot(simit)
    max_tot_mod = maxm(step_tot,nsims)
    retau(:,simit) = utau(:,simit)/nu
    re(:,simit) = u_bulk(:,simit)/nu
!
!	Next is STAT_0D-2.dat which has the info from every ntmavg time steps as:
!	ii,time,dt,utau,u_bulk,esum,energy
!	time_step=1
!	DO i = nstart(simit), nruns(simit)
!	    file3 = 'STAT_0D-2.dat'
!	    CALL change_dir(file3,i)
!	    Write(*,*) file3
!	    Open(5, File=file3, Status='OLD')
!	    Do j = 1,1000000
 !       	READ(5,*,END=14) ii,time_v(time_step,simit),dt_v(time_step,simit),utau_v(time_step,simit),&
!				u_bulk_v(time_step,simit),esum_v(time_step,simit),energy_v(time_step,simit)
!		time_step = time_step + 1	
!            End Do
!    14	    Close(5)
!	    Writes over last time step if ins_1d is not written on last step
!	    If (nsample(i-nstart(simit)+1,simit)==ins_1d(i-nstart(simit)+1,simit)+1) time_step=time_step-1
!	END DO
!	Print*,'Total number of 1d instaneous snapshots', time_step
!	retau_v(:,simit) = utau_v(:,simit)/nu
!	re_v(:,simit) = u_bulk_v(:,simit)/nu

  end do
! get the averaged region: tmin-tmax
  nmin(:)=1;nmax(:)=step_tot(:)
  Do simit=1,nsims
    If(tmin.Ne.0.0e0) Then
    Do i=1,step_tot(simit)
      If(time(i,simit).Ge.tmin) Then
        nmin(simit)=i
        Exit
      End If
    End Do; End If
    If(tmax.Ne.0.0e0) Then
    Do i=step_tot(simit),1,-1
      If(time(i,simit).le.tmax) Then
        nmax(simit)=i
        Exit
      End If
    End Do; End If
    Write(*,*) 'simit=',simit,'nmin=',nmin,'nmax=',nmax
  End Do
! do raw average
  If(iraw == 1) Then
    Do simit=1,nsims
      write(simno,*) simit
      simno=adjustl(simno)      
      If(nsims.Eq.1) Then
	unitmp='raw_av_Cf.plt'
      Else
	unitmp='raw_av_Cf_'//simno//'.plt'
      End If
      CALL raw_avg(utau(:,simit),simit,unitmp)
    End Do
  End If

  return
end subroutine read_0d_data

!*******************************************************************************

subroutine calc_0d_avs
! This calculates some 0d statistics
  integer :: i

  utau_avg_sa = real(sum(utau_sa),kind=mytype)/real(tpts,kind=mytype)
  retau_avg_sa = utau_avg_sa/nu
  re_avg_sa = real(sum(re_sa),kind=mytype)/real(tpts,kind=mytype)
  ub_avg_sa = real(sum(u_bulk_sa),kind=mytype)/real(tpts,kind=mytype)
  Cf_sa = 2.0e0*(utau_avg_sa**2)/(ub_avg_sa**2)

  write(*,*) '****************************************'
  write(*,*) 'Re average = ', re_avg_sa
  write(*,*) 'Retau average = ', retau_avg_sa
  write(*,*) 'utau average = ', utau_avg_sa
  write(*,*) 'u bulk average = ', ub_avg_sa
  write(*,*) 'Cf = ', Cf_sa
  write(*,*) '****************************************'

  if (reinit .ge. small) then
    utau_avg_sa = reinit*nu
    retau_avg_sa =  reinit
    re_avg_sa = 1.0/nu
    !utau(:,simit) = reinit*nu
    !retau(:,simit) = reinit
    !re(:, simit) = 1.0/nu
  end if

! Scale grid to wall units
  xplus=x*retau_avg_sa
  yplus=y*retau_avg_sa
  zplus=z*retau_avg_sa
  xcplus=xc*retau_avg_sa
  ycplus=yc*retau_avg_sa
  zcplus=zc*retau_avg_sa

  return
end subroutine calc_0d_avs

!*******************************************************************************

subroutine calc_unst
  integer :: i

  if ((isteady == 0).or.(isteady == 2)) then
    allocate(flucs1d_sa(-1:ny+2,4,4,tpts))
    if (itke == 1)   allocate(kinen1d_sa(-1:ny+2,4,6,tpts))
    if (ishear == 1) allocate(shear1d_sa(-1:ny+2,4,tpts))
    if (iangle ==1) allocate(angle1d_sa(-1:ny+2,5,tpts))
    if (ianis == 1)  allocate(anisb1d_sa(-1:ny+2,4,tpts),anisII1d_sa(-1:ny+2,tpts), &
                              anisIII1d_sa(-1:ny+2,tpts),anisf1d_sa(-1:ny+2,tpts),  &
                              anisg1d_sa(-1:ny+2,tpts))
    if (ivort == 1) then; allocate(omegamn1d_sa(-1:ny+2,3,tpts)); omegamn1d_sa=0.0e0; end if
  end If
       
  allocate(xplus_sa(tpts,nx+1),yplus_sa(tpts,ny+1), &
           zplus_sa(tpts,nz+1),xcplus_sa(tpts,0:nx+1), &
           ycplus_sa(tpts,0:ny+1),zcplus_sa(tpts,0:nz+1))

  do i = 1,tpts
    xplus_sa(i,:) = retau_sa(i)*x
    xcplus_sa(i,:) = retau_sa(i)*xc
    yplus_sa(i,:) = retau_sa(i)*y
    ycplus_sa(i,:) = retau_sa(i)*yc
    zplus_sa(i,:) = retau_sa(i)*z
    zcplus_sa(i,:) = retau_sa(i)*zc
  end do

  return
end subroutine calc_unst

!*******************************************************************************

subroutine avg_sims_setup
! This subroutine sets up the avg over the simulations
  integer :: i,j,k,b,a,m
  Real(mytype) :: c,d,dum,dum2
  Real(mytype), dimension(:), allocatable :: dum0
!If tpts isn't set we put the time as in sim1
!time=amax1(time,0.e0)

  if (tmin .le. 0) then
    call minm(time,max_tot,0,tmin)
  end if
  if (tmax .le. 0) then
    call minm(time,max_tot,1,tmax)
  end if

  if (tperiod .GT. 0.0e0) then
    nperiod=ceiling(tmin/tperiod) 
    tmin=real(nperiod,kind=mytype)*tperiod    
    nperiod=floor((tmax-tmin)/tperiod)
    tmax=real(nperiod,kind=mytype)*tperiod+tmin
    write(*,*) nperiod, 'periods, of size',tperiod
    tpts=nperiod*tpts
    write(*,*) tpts, 'time points in total'
  end if
  if (isteady == 2 .and. lpstat_2d == 1) tpts=nxozt

! 1d ins arrays
  allocate(time_sa(tpts),re_sa(tpts),rev_sa(tpts),utau_sa(tpts),dpdx_sa(tpts),retau_sa(tpts),u_bulk_sa(tpts))
  if ((lp_ins_1d == 1 .and. isteady == 0).or.(lpstat_2d == 1 .and. isteady == 2)) then
    allocate(stat1d_sa(-1:ny+2,4,4,tpts),sqar1d_sa(-1:ny+2,4,4,tpts))
    if (iquad == 1) allocate(quad1d_sa(-1:ny+2,9,tpts))
    if (ivort == 1) allocate(omega1d_sa(-1:ny+2,6,tpts),wx1d_sa(-1:ny+2,2,8,tpts))
    if (itke == 1)  allocate(tke1d_sa(-1:ny+2,4,7,tpts))
  end if

  do i=1,tpts
    time_sa(i)=(real(i-1,kind=mytype)*(tmax-tmin)/real(tpts-1,kind=mytype)) + tmin
  end do
  print*,'tmin-tmax',tmin,tmax
  smin=1; smax=0

  if (sim_print.le.0) then
    smin=1; smax=nsims
  else
    do i=1,nsims
      if (sim_print == i) then
        smin=i;smax=i
      end if 
    end do
  end if
  simtot=smax-smin+1

! Avg utau,retau,etc
  re_sa=0; retau_sa=0; utau_sa=0; dpdx_sa=0; u_bulk_sa=0
  do i=1,tpts!-1
    do simit=smin,smax
!     arrays in standard time
      call above(time,max_tot,time_sa(i),simit,step_tot(simit),a)
      call below(time,max_tot,time_sa(i),simit,step_tot(simit),b)
      c=time(a,simit)-time_sa(i)
      d=time_sa(i)-time(b,simit)
      if (c .Lt. 0.0e0 .or. d .Lt. 0.0e0) then; print*,'c or d less than 0 -11111 '; end if
      if((c+d)==0.0) then; c=1; d=1; end if
      re_sa(i)=re_sa(i)+(d*re(a,simit)+c*re(b,simit))/(c+d)
      if (i == tpts) then
	print*,'check-re:'
	write(*,'(3I8,5F12.5)') a,b,max_tot,c,d,re(a,simit),re(b,simit),time_sa(i)
      endif
      retau_sa(i)=retau_sa(i)+(d*retau(a,simit)+c*retau(b,simit))/(c+d)
      utau_sa(i)=utau_sa(i)+(d*utau(a,simit)+c*utau(b,simit))/(c+d)
      u_bulk_sa(i)=u_bulk_sa(i)+(d*u_bulk(a,simit)+c*u_bulk(b,simit))/(c+d)
      dpdx_sa(i)=dpdx_sa(i)+(d*dpdx_mean(a,simit)+c*dpdx_mean(b,simit))/(c+d)	
    end do
  end do
  deallocate(re,retau,utau,dpdx_mean)

  re_sa=re_sa/real(simtot,kind=mytype)
  retau_sa=retau_sa/real(simtot,kind=mytype)
  utau_sa=utau_sa/real(simtot,kind=mytype)
  u_bulk_sa=u_bulk_sa/real(simtot,kind=mytype)
  dpdx_sa=dpdx_sa/real(simtot,kind=mytype)

! This starts the time to at t=0
  if (tstart0==1) then
    start_time=time_sa(1)
    time_sa(:)=time_sa(:)-start_time
  end if	

  return
end subroutine avg_sims_setup

!*******************************************************************************

subroutine read_1d_avg
! This reads the sum info from 'OUTPUT_1D_AVG.dat' and calculates
! the average.
  character :: file4*40,sqs*8,os*9,ts*7,ws*6
  integer :: i,j,k,n,overall_steps,rec_val,rec_sz,i2,j2,k2
  real(mytype) :: dum_val
! The following are temp arrays to read the values
  real(mytype), dimension(-1:ny+2,9) :: quad_sum
  real(mytype), dimension(-1:ny+2,6) :: omega_sum
  real(mytype), dimension(-1:ny+2,4,7) :: tke_sum
  real(mytype), dimension(-1:ny+2,4,4) :: stat_sum
  real(mytype), dimension(-1:ny+2,2,8) :: wx_sum
!
  overall_steps=0
  stat_avg=0.0e0
  if (itke == 1) tke_avg=0.0e0
  if (ivort == 1) wx_avg=0.0e0
  if (ivort == 1) omega_avg=0.0e0
  if (iquad == 1) quad_avg=0.0e0
!
! Read in the sum values from each file and sum them
  do simit=1,nsims
    do i = nstart(simit), nruns(simit)
      file4 = 'OUTPUT_1D_AVG.dat'
      call change_dir(file4,i)
      write(*,*) 'reading...',file4
      if (mpi == 0) then
        open(6, file=file4, status='OLD', form='UNFormatTED')
        call header_skip(6)
        read(6) sqs
        read(6) stat_sum(1:ny+1,1:4,1:4)
        read(6) sqs
        read(6) quad_sum(1:ny+1,1:9)
        read(6) os
        read(6) omega_sum(1:ny+1,1:4)
        read(6) ts
        read(6) tke_sum(1:ny+1,1:4,1:7)
        read(6) ws
        read(6) wx_sum(1:ny+1,1:2,1:8)
        close(6)
      else        
        open(6, file=file4, status='old', form='unformatted', &
                access='stream',carriagecontrol='none')

        do k2=1,4; do j2=1,4; do i2=-1,ny+2
          read(6) stat_sum(i2,j2,k2)       
        end do; end do; end do

        do j2=1,9; do i2=-1,ny+2
          read(6) quad_sum(i2,j2)
        end do; end do

        do j2=1,6; do i2=-1,ny+2
          read(6) omega_sum(i2,j2)
        end do; end do

        do k2=1,7; do j2=1,4; do i2=-1,ny+2
          read(6) tke_sum(i2,j2,k2)
        end do; end do; end do

        do k2=1,8; do j2=1,2; do i2=-1,ny+2
          read(6) wx_sum(i2,j2,k2)
        end do; end do; end do

        close(6)
      end if
      stat_avg(-1:ny+2,:,:) = stat_avg(-1:ny+2,:,:) + stat_sum
      if (iquad == 1) quad_avg(-1:ny+2,:) = quad_avg(-1:ny+2,:) + quad_sum
      if (ivort == 1) omega_avg(-1:ny+2,:) = omega_avg(-1:ny+2,:) + omega_sum
      if (itke == 1) tke_avg(-1:ny+2,:,:) = tke_avg(-1:ny+2,:,:) + tke_sum
      if (ivort == 1) wx_avg(-1:ny+2,:,:) = wx_avg(-1:ny+2,:,:) + wx_sum
    end do
    overall_steps=overall_steps+samp_tot(simit)
	!stat_avg_sa(:,:,:)=stat_avg_sa(:,:,:)+samp_tot(simit)*stat_avg(1:nyt,:,:)
	!If (itke == 1) tke_avg_sa(:,:,:)=tke_avg_sa(:,:,:)+samp_tot(simit)*tke_avg(1:nyt,:,:)
	!If (ivort == 1) wx_avg_sa(:,:,:)=wx_avg_sa(:,:,:)+samp_tot(simit)*wx_avg(1:nyt,:,:)
	!If (ivort == 1) omega_avg_sa(:,:)=omega_avg_sa(:,:)+samp_tot(simit)*omega_avg(1:nyt,:)
	!If (iquad == 1) quad_avg_sa(:,:)=quad_avg_sa(:,:)+samp_tot(simit)*quad_avg(1:nyt,:)
  end do

  stat_avg = stat_avg/real(overall_steps,kind=mytype)
  tke_avg = tke_avg/real(overall_steps,kind=mytype)
  wx_avg = wx_avg/real(overall_steps,kind=mytype)
  omega_avg = omega_avg/real(overall_steps,kind=mytype)
  quad_avg = quad_avg/real(overall_steps,kind=mytype)

! Update boundary values
  if (mpi == 0) then
    stat_avg(nyt,:,:) =stat_avg(ny,:,:)
    stat_avg(nyt,1:3,:) =0.0e0
    quad_avg(nyt,:) = quad_avg(ny,:)
    omega_avg(nyt,:) = omega_avg(ny,:)
    tke_avg(nyt,:,:) = tke_avg(ny,:,:)
    wx_avg(nyt,:,:) = wx_avg(ny,:,:)
!
    stat_avg(nyt+1,:,:) =stat_avg(ny,:,:)
    stat_avg(nyt+1,1:3,:) =0.0e0
    quad_avg(nyt+1,:) = quad_avg(ny,:)
    omega_avg(nyt+1,:) = omega_avg(ny,:)
    tke_avg(nyt+1,:,:) = tke_avg(ny,:,:)
    wx_avg(nyt+1,:,:) = wx_avg(ny,:,:)

    stat_avg(0,:,:) =stat_avg(1,:,:)
    stat_avg(0,1:3,:) =0.0e0
    quad_avg(0,:) = quad_avg(1,:)
    omega_avg(0,:) = omega_avg(1,:)
    tke_avg(0,:,:) = tke_avg(1,:,:)
    wx_avg(0,:,:) = wx_avg(1,:,:)
!
    stat_avg(-1,:,:) =stat_avg(1,:,:)
    stat_avg(-1,1:3,:) =0.0e0
    quad_avg(-1,:) = quad_avg(1,:)
    omega_avg(-1,:) = omega_avg(1,:)
    tke_avg(-1,:,:) = tke_avg(1,:,:)
    wx_avg(-1,:,:) = wx_avg(1,:,:)
  end if
! Move v in stat to cell centre
  if (mpi == 0) then
    do j=0,nyt; do n=1,4
      stat_avg(j,2,n)=(stat_avg(j,2,n)+stat_avg(j+1,2,n))/2.0
    end do; end do
  end if

  return
end subroutine read_1d_avg

!*******************************************************************************

subroutine read_2d_avg
! This reads the sum info from 'OUTPUT_1D_AVG.dat' and calculates
! the average.
  character :: file5*40,sqs*10,os*11,ts*9,ws*8,chtmp*40
  integer :: i,j,k,n,overall_steps,ib,i2,j2,k2,n2,rec_val,rec_sz,m,mm
  real(mytype) :: tphase
! The following are temp arrays to read the values
  real(mytype), dimension(:,:,:), allocatable :: quad2d_sum,omega2d_sum
  real(mytype), dimension(:,:,:,:), allocatable :: tke2d_sum,stat2d_sum,wx2d_sum
  allocate(quad2d_sum(-1:nxoz+2,-1:ny+2,9),omega2d_sum(-1:nxoz+2,-1:ny+2,6), &
           tke2d_sum(-1:nxoz+2,-1:ny+2,4,7),stat2d_sum(-1:nxoz+2,-1:ny+2,4,4), &
           wx2d_sum(-1:nxoz+2,-1:ny+2,2,8))

  ib=1
  overall_steps=0; mm=0
  stat1d_sa = 0.0e0
  if (iquad == 1) quad1d_sa = 0.0e0
  if (ivort == 1) omega1d_sa = 0.0e0
  if (itke == 1)tke1d_sa = 0.0e0
  if (ivort == 1) wx1d_sa = 0.0e0

! Read in the sum values fvi irom each file and sum them
  do simit=1,nsims
    do i = nstart(simit), nruns(simit)
!      file5 = 'OUTPUT_2D_AVG.dat'
      file5=file2davg
      call change_dir(file5,i)
      chtmp='OUTPUT_2D_AVG.dat'
      call change_dir(chtmp,i)
      write(*,*) file5
      if (mpi == 0) then
        open(7, File=file5, Status='OLD', Form='UNFormatTED')
        call header_skip(7)
        read(7) sqs
        read(7) stat2d_sum(1:nxoz+1,1:ny+1,4,4)
        read(7) sqs
        read(7) quad2d_sum(1:nxoz+1,1:ny+1,9)
        read(7) os
        read(7) omega2d_sum(1:nxoz+1,1:ny+1,4)
        read(7) ts
        read(7) tke2d_sum(1:nxoz+1,1:ny+1,4,7)
        read(7) ws
        read(7) wx2d_sum(1:nxoz+1,1:ny+1,2,8)
        close(7)
      else
        open(7, file=file5, status='old', form='unformatted', &
                access='stream',carriagecontrol='none')
	if(Trim(Adjustl(file5)) .ne. chtmp) then
	  read(7) tphase
	  read(7) m
	  write(*,*) 'phase',tphase,'with snapshots:',m
	end if
	if(lpstat_x == 1) then
          do n2=1,4; do k2=1,4; do i2=-1,nxoz+2; do j2=-1,ny+2
            read(7) stat2d_sum(i2,j2,k2,n2)
          end do; end do; end do; end do

          do k2=1,9; do i2=-1,nxoz+2; do j2=-1,ny+2
            read(7) quad2d_sum(i2,j2,k2)
          end do; end do; end do

          do k2=1,6; do i2=-1,nxoz+2; do j2=-1,ny+2
            read(7) omega2d_sum(i2,j2,k2)
          end do; end do; end do

          do n2=1,7; do k2=1,4; do i2=-1,nxoz+2; do j2=-1,ny+2
            read(7) tke2d_sum(i2,j2,k2,n2)
          end do; end do; end do; end do

          do n2=1,8; do k2=1,2; do i2=-1,nxoz+2; do j2=-1,ny+2
            read(7) wx2d_sum(i2,j2,k2,n2)
          end do; end do; end do; end do

        else if(lpstat_z == 1) then
          do n2=1,4; do k2=1,4; do j2=-1,ny+2; do i2=-1,nxoz+2
            read(7) stat2d_sum(i2,j2,k2,n2)
          end do; end do; end do; end do

          do k2=1,9; do j2=-1,ny+2; do i2=-1,nxoz+2
            read(7) quad2d_sum(i2,j2,k2)
          end do; end do; end do

          do k2=1,6; do j2=-1,ny+2; do i2=-1,nxoz+2
            read(7) omega2d_sum(i2,j2,k2)
          end do; end do; end do

          do n2=1,7; do k2=1,4; do j2=-1,ny+2; do i2=-1,nxoz+2
            read(7) tke2d_sum(i2,j2,k2,n2)
          end do; end do; end do; end do

          do n2=1,8; do k2=1,2; do j2=-1,ny+2; Do i2=-1,nxoz+2
            read(7) wx2d_sum(i2,j2,k2,n2)
          end do; end do; end do; end do

        end if
	Write(*,*) 'nxoz=',nxoz
!	Write(*,*) omega2d_sum(:,:,1)
        close(7)
     end if
     mm=mm+m
     do n=1,nxozt
       stat1d_sa(-1:ny+2,:,:,n) = stat1d_sa(-1:ny+2,:,:,n) + stat2d_sum(n,-1:ny+2,:,:)
       quad1d_sa(-1:ny+2,:,n) = quad1d_sa(-1:ny+2,:,n) + quad2d_sum(n,-1:ny+2,:)
       omega1d_sa(-1:ny+2,:,n) = omega1d_sa(-1:ny+2,:,n) + omega2d_sum(n,-1:ny+2,:)
       tke1d_sa(-1:ny+2,:,:,n) = tke1d_sa(-1:ny+2,:,:,n) + tke2d_sum(n,-1:ny+2,:,:)
       wx1d_sa(-1:ny+2,:,:,n) = wx1d_sa(-1:ny+2,:,:,n) + wx2d_sum(n,-1:ny+2,:,:)
     end do 
    end do
    overall_steps=overall_steps+samp_tot(simit)
  end do

! Divide by total number of time steps to get the average
  if(file5 .eq. chtmp) mm=overall_steps
  stat1d_sa(:,:,:,:) = stat1d_sa(:,:,:,:)/real(mm,kind=mytype)
  quad1d_sa(:,:,:)  = quad1d_sa(:,:,:) /real(mm,kind=mytype)
  omega1d_sa(:,:,:)  = omega1d_sa(:,:,:) /real(mm,kind=mytype)
  tke1d_sa(:,:,:,:) = tke1d_sa(:,:,:,:)/real(mm,kind=mytype)
  wx1d_sa(:,:,:,:) = wx1d_sa(:,:,:,:)/real(mm,kind=mytype)

! Move v in stat to cell centre
  if (mpi == 0) then
    stat1d_sa(0,:,:,:)=stat1d_sa(1,:,:,:)
    do j=0,nyt; do n=1,4
      stat1d_sa(j,2,n,:)=(stat1d_sa(j,2,n,:)+stat1d_sa(j+1,2,n,:))/2.0
    end do; end do
  end if

! periodic in x
!  if(mpi == 0) then
!    stat1d_sa(:,:,:,0) = 0.5e0*(stat1d_sa(:,:,:,1)+stat1d_sa(:,:,:,nx))
!    quad1d_sa(:,:,0)  = 0.5e0*(quad1d_sa(:,:,1)+quad1d_sa(:,:,nx))
!    omega1d_sa(:,:,0)  = 0.5e0*(omega1d_sa(:,:,1)+omega1d_sa(:,:,nx))
!    tke1d_sa(:,:,:,0) = 0.5e0*(tke1d_sa(:,:,:,1)+tke1d_sa(:,:,:,nx))
!    wx1d_sa(:,:,:,0) = 0.5e0*(wx1d_sa(:,:,:,1)+wx1d_sa(:,:,:,nx))

    stat1d_sa(:,:,:,nxozt) = stat1d_sa(:,:,:,1)
    quad1d_sa(:,:,nxozt)  = quad1d_sa(:,:,1)
    omega1d_sa(:,:,nxozt)  = omega1d_sa(:,:,1)
    tke1d_sa(:,:,:,nxozt) = tke1d_sa(:,:,:,1)
    wx1d_sa(:,:,:,nxozt) = wx1d_sa(:,:,:,1)
!  end if

! Move u in stat to cell centre
  if (mpi == 0) then
    do i=1,nxozt; do n=1,4
      stat1d_sa(:,1,n,i)=(stat1d_sa(:,1,n,i)+stat1d_sa(:,1,n,i+1))/2.0
    end do; end do
  end if
    
! Calculates utau/retau
  if(lpstat_x == 1) then
    time_sa(1:nxoz)=z(1:nxoz)
    time_sa(nxozt)=width+z(1)
  else if(lpstat_z == 1) then
    time_sa(1:nxoz)=x(1:nxoz)
    time_sa(nxozt)=length+x(1)
  end if
  utau_sa=0.0e0
  do i=1,nxoz
    utau_sa(i)=utau_sa(i)+sqrt(abs(nu*(gcdy_f(1,1,1,1)*stat1d_sa(2,1,1,i)+gcdy_f(1,2,1,1)*stat1d_sa(1,1,1,i))))
    utau_sa(i)=utau_sa(i)+sqrt(abs(nu*(gcdy_f(nyt,3,1,1)*stat1d_sa(ny,1,1,i)+gcdy_f(nyt,4,1,1)*stat1d_sa(ny-1,1,1,i))))
  end do
  utau_sa(tpts)=utau_sa(1)
  utau_sa=utau_sa/2.0e0
  retau_sa=utau_sa/nu

! Calculates mean vorticity
!  omegamn1d_sa=0.0e0
!  do j=1,nyt; do i=1,tpts
!    omegamn1d_sa(j,2,i) = -(gcdx_c(i,2,2,ib)*stat1d_sa(j,3,1,i+1)+ &
!                            gcdx_c(i,3,2,ib)*stat1d_sa(j,3,1,i)+ &
!                            gcdx_c(i,4,2,ib)*stat1d_sa(j,3,1,i-1))
!    omegamn1d_sa(j,3,i) = gcdx_c(i,2,2,ib)*stat1d_sa(j,2,1,i+1)+ &
!	                  gcdx_c(i,3,2,ib)*stat1d_sa(j,2,1,i)+ &
!                          gcdx_c(i,4,2,ib)*stat1d_sa(j,2,1,i-1)
!  end do; end do

  do m=1,4
    sqar1d_sa(:,m,2,:)=stat1d_sa(:,m,1,:)**2
    sqar1d_sa(:,m,3,:)=stat1d_sa(:,m,1,:)**3
    sqar1d_sa(:,m,4,:)=stat1d_sa(:,m,1,:)**4
  end do
  
  nperiod=Nint(length/wavelen)
  If(lpstat_x.Eq.1) nperiod=Nint(width/wavelen)
  tperiod=wavelen_plus/10.0e0
  Write(*,*) '>>>>>'
  Write(*,'(I5,A,F5.2)') nperiod,' waves in the domain with wave length',wavelen
! do period average
  If(tstart0 == 1) Then
    tpts=tpts/nperiod*nperiod
    Call phase_ins_1d
  End If
  isteady=0

  return
end subroutine read_2d_avg

!*******************************************************************************

subroutine read_ins_1d
! Reads the 1D Instantaneous info
  character :: file7*40
  integer :: i,j,k,time_point,ii,n,steps_used,div,b,a,m,rec_val,rec_sz,i2,j2,k2
  real(mytype) :: tdum,dum,dum2,c,d
  real(mytype), dimension(nyt) :: sums
  real(mytype), dimension(:), allocatable :: tins_1d
  real(mytype), dimension(:,:,:,:), allocatable :: stat1d,tke1d,wx1d
  real(mytype), dimension(:,:,:), allocatable :: quad1d,omega1d,vortmn
!
  allocate(stat1d(-1:ny+2,4,4,ins_1d_max),quad1d(-1:ny+2,9,ins_1d_max), &
           omega1d(-1:ny+2,6,ins_1d_max),tke1d(-1:ny+2,4,7,ins_1d_max), &
           wx1d(-1:ny+2,2,8,ins_1d_max),tins_1d(ins_1d_max),           &
           vortmn(-1:ny+2,3,ins_1d_max))

! Zero arrays
  stat1d_sa=0.0e0
  if (iquad == 1) quad1d_sa=0.0e0
  if (ivort == 1) omega1d_sa=0.0e0
  if (itke == 1) tke1d_sa=0.0e0
  if (ivort == 1) wx1d_sa=0.0e0

  do simit=1,nsims
    time_point=0
    do i = nstart(simit), nruns(simit)
      file7 = 'OUTPUT_1D_INS.dat'
      call change_dir(file7,i)
      write(*,*) file7
      if (mpi == 0) then
        open(11, file=file7, status='OLD', form='UNFormatTED')
        call header_skip(11)
!       loop over number of samples
        k = i-nstart(simit)+1
        do j=1, ins_1d(k,simit)
          time_point=time_point+1
          read(11) ii,tdum
          tins_1d(time_point) = tdum
          read(11) stat1d(1:nyt,:,:,time_point)
          read(11) quad1d(1:nyt,:,time_point)
          read(11) omega1d(1:nyt,:,time_point)
          read(11) tke1d(1:nyt,:,:,time_point)
          read(11) wx1d(1:nyt,:,:,time_point)
        end do
        close(11)
      else
     !   open(11, file=file7, status='old', form='unformatted', &
     !            access='direct',recl=mytype)
        open(11, file=file7, status='old', form='unformatted', &
                access='stream',carriagecontrol='none')
        k = i-nstart(simit)+1
        write(*,*) 'num:',ins_1d(k,simit)
        do j=1,ins_1d(k,simit)
          time_point=time_point+1
         !read(11) ii
          read(11) tdum
          tins_1d(time_point) = tdum 

          do k2=1,4; do j2=1,4; do i2=-1,ny+2
            read(11) stat1d(i2,j2,k2,time_point)
          end do; end do; end do

          do j2=1,9; do i2=-1,ny+2
            read(11) quad1d(i2,j2,time_point)
          end do; end do

          do j2=1,6; do i2=-1,ny+2
            read(11) omega1d(i2,j2,time_point)
          end do; End do

          do k2=1,7; do j2=1,4; do i2=-1,ny+2
            read(11) tke1d(i2,j2,k2,time_point)
          end do; end do; end do

          do k2=1,8; do j2=1,2; do i2=-1,ny+2
            read(11) wx1d(i2,j2,k2,time_point)
          end do; end do; end do
        end do
        close(11)
      end if
    end do
    if (tstart0==1) tins_1d(:)=tins_1d(:)-start_time
!
    do i=1,ins_1d_max
      call vorty(stat1d(:,:,:,i),vortmn(:,:,i))
    end do

!   Interpolate and sum
    do i=1,tpts
      call above(tins_1d,ins_1d_max,time_sa(i),simit,ins_1d_tot(simit),a)
      call below(tins_1d,ins_1d_max,time_sa(i),simit,ins_1d_tot(simit),b)
      c=tins_1d(a)-time_sa(i)
      d=time_sa(i)-tins_1d(b)	
      if((c+d)==0.0) then; c=1; d=1; end if
      if (c .Lt. 0.0e0 .or. d .Lt. 0.0e0) then; print*,'ins: c or d less than 0 '; stop; end if
      stat1d_sa(-1:ny+2,:,:,i)=stat1d_sa(-1:ny+2,:,:,i)+(d*stat1d(-1:ny+2,:,:,a)+c*stat1d(-1:ny+2,:,:,b))/(c+d)
      sqar1d_sa(-1:ny+2,1:3,1,i)=sqar1d_sa(-1:ny+2,1:3,1,i)+(d*vortmn(-1:ny+2,:,a)**2+c*vortmn(-1:ny+2,:,b)**2)/(c+d)
      sqar1d_sa(-1:ny+2,:,2,i)=sqar1d_sa(-1:ny+2,:,2,i)+(d*stat1d(-1:ny+2,:,1,a)**2+c*stat1d(-1:ny+2,:,1,b)**2)/(c+d)
      sqar1d_sa(-1:ny+2,:,3,i)=sqar1d_sa(-1:ny+2,:,3,i)+(d*stat1d(-1:ny+2,:,1,a)**3+c*stat1d(-1:ny+2,:,1,b)**3)/(c+d)
      sqar1d_sa(-1:ny+2,:,4,i)=sqar1d_sa(-1:ny+2,:,4,i)+(d*stat1d(-1:ny+2,:,1,a)**4+c*stat1d(-1:ny+2,:,1,b)**4)/(c+d)
      if (iquad == 1) quad1d_sa(-1:ny+2,:,i)=quad1d_sa(-1:ny+2,:,i)+(d*quad1d(-1:ny+2,:,a)+c*quad1d(-1:ny+2,:,b))/(c+d)
      if (ivort == 1) omega1d_sa(-1:ny+2,:,i)=omega1d_sa(-1:ny+2,:,i)+(d*omega1d(-1:ny+2,:,a)+c*omega1d(-1:ny+2,:,b))/(c+d)
      if (itke == 1) tke1d_sa(-1:ny+2,:,:,i)=tke1d_sa(-1:ny+2,:,:,i)+(d*tke1d(-1:ny+2,:,:,a)+c*tke1d(-1:ny+2,:,:,b))/(c+d)
      if (ivort == 1) wx1d_sa(-1:ny+2,:,:,i)=wx1d_sa(-1:ny+2,:,:,i)+(d*wx1d(-1:ny+2,:,:,a)+c*wx1d(-1:ny+2,:,:,b))/(c+d)
    end do
  end do

  stat1d_sa=stat1d_sa/real(nsims,kind=mytype)
  sqar1d_sa=sqar1d_sa/real(nsims,kind=mytype)
  if (iquad == 1) quad1d_sa=quad1d_sa/real(nsims,kind=mytype)
  if (ivort == 1) omega1d_sa=omega1d_sa/real(nsims,kind=mytype)
  if (itke == 1) tke1d_sa=tke1d_sa/real(nsims,kind=mytype)
  if (ivort == 1) wx1d_sa=wx1d_sa/real(nsims,kind=mytype)

! Update boundary values
  if (mpi == 0) then
    stat1d_sa(0,:,:,:)=stat1d_sa(1,:,:,:)
    sqar1d_sa(0,:,:,:)=sqar1d_sa(1,:,:,:)
    if (iquad == 1) quad1d_sa(0,:,:)=quad1d_sa(1,:,:)
    if (ivort == 1) omega1d_sa(0,:,:)=omega1d_sa(1,:,:)
    if (itke == 1) tke1d_sa(0,:,:,:)=tke1d_sa(1,:,:,:)
    if (ivort == 1) wx1d_sa(0,:,:,:)=wx1d_sa(1,:,:,:)

    stat1d_sa(nyt+1,:,:,:)=stat1d_sa(nyt,:,:,:)
    sqar1d_sa(nyt+1,:,:,:)=sqar1d_sa(nyt,:,:,:)
    if (iquad == 1) quad1d_sa(nyt+1,:,:)=quad1d_sa(nyt,:,:)
    if (ivort == 1) omega1d_sa(nyt+1,:,:)=omega1d_sa(nyt,:,:)
    if (itke == 1) tke1d_sa(nyt+1,:,:,:)=tke1d_sa(nyt,:,:,:)
    if (ivort == 1) wx1d_sa(nyt+1,:,:,:)=wx1d_sa(nyt,:,:,:)
  end if

! Move v in stat to cell centre
  if (mpi == 0) then
    do j=0,nyt
      stat1d_sa(j,2,1:4,:)=(stat1d_sa(j,2,1:4,:)+stat1d_sa(j+1,2,1:4,:))/2.0
      sqar1d_sa(j,2,2:4,:)=(sqar1d_sa(j,1,2:4,:)+sqar1d_sa(j+1,2,2:4,:))/2.0
    end do 
  end if
  call phase_ins_1d

  !if (hybrid == 1) then
  !  stat_avg=0.0e0; quad_avg=0.0e0; omega_avg=0.0e0; tke_avg=0.0e0; wx_avg=0.0e0
  !  do i=1,tpts
  !    stat_avg(0:nyt+1,:,:)=stat_avg(0:nyt+1,:,:)+stat1d(0:nyt+1,:,:,i)
  !    quad_avg(0:nyt+1,:)=quad_avg(0:nyt+1,:)+quad1d(0:nyt+1,:,i)
  !    omega_avg(0:nyt+1,:)=omega_avg(0:nyt+1,:)+omega1d(0:nyt+1,:,i)
  !    tke_avg(0:nyt+1,:,:)=tke_avg(0:nyt+1,:,:)+tke1d(0:nyt+1,:,:,i)
  !    wx_avg(0:nyt+1,:,:)=wx_avg(0:nyt+1,:,:)+wx1d(0:nyt+1,:,:,i)
  !  end do
  !  stat_avg(:,:,:)=stat_avg(:,:,:)/real(tpts,kind=mytype)
  !  quad_avg(:,:)=quad_avg(:,:)/real(tpts,kind=mytype)
  !  omega_avg(:,:)=omega_avg(:,:)/real(tpts,kind=mytype)
  !  tke_avg(:,:,:)=tke_avg(:,:,:)/real(tpts,kind=mytype)
  !  wx_avg(:,:,:)=wx_avg(:,:,:)/real(tpts,kind=mytype)
  !  isteady = 1	
  !end if

  return
end subroutine read_ins_1d

!*******************************************************************************

subroutine phase_ins_1d
  integer :: i,j,k,n,steps_used,div,m
  real(mytype), dimension(:,:,:,:), allocatable :: dum_stat,dum_tke,dum_wx,dum_sqar
  real(mytype), dimension(:,:,:), allocatable :: dum_quad,dum_omega
  real(mytype), dimension(:), allocatable :: dum_utau

! phase averaging
  if (tperiod .gt. 0.0e0) then
    allocate(dum_stat(-1:ny+2,4,4,tpts/nperiod),dum_quad(-1:ny+2,9,tpts/nperiod), &
             dum_omega(-1:ny+2,6,tpts/nperiod),dum_tke(-1:ny+2,4,7,tpts/nperiod),	&
             dum_wx(-1:ny+2,2,8,tpts/nperiod),dum_sqar(-1:ny+2,4,4,tpts/nperiod), &
             dum_utau(tpts/nperiod+1))

    dum_utau=0.0e0; dum_stat=0.0e0; dum_sqar=0.0e0; dum_quad=0.0e0; dum_omega=0.0e0; dum_tke=0.0e0; dum_wx=0.0e0
    do i=1,tpts
      dum_utau(mymod(i,tpts/nperiod))=dum_utau(mymod(i,tpts/nperiod))+utau_sa(i)
      dum_stat(:,:,:,mymod(i,tpts/nperiod))=dum_stat(:,:,:,mymod(i,tpts/nperiod))+stat1d_sa(:,:,:,i)
      dum_sqar(:,:,:,mymod(i,tpts/nperiod))=dum_sqar(:,:,:,mymod(i,tpts/nperiod))+sqar1d_sa(:,:,:,i)
      if (iquad == 1) dum_quad(:,:,mymod(i,tpts/nperiod))=dum_quad(:,:,mymod(i,tpts/nperiod))+quad1d_sa(:,:,i)
      if (ivort == 1) dum_omega(:,:,mymod(i,tpts/nperiod))=dum_omega(:,:,mymod(i,tpts/nperiod))+omega1d_sa(:,:,i)
      if (itke == 1) dum_tke(:,:,:,mymod(i,tpts/nperiod))=dum_tke(:,:,:,mymod(i,tpts/nperiod))+tke1d_sa(:,:,:,i)
      if (ivort == 1) dum_wx(:,:,:,mymod(i,tpts/nperiod))=dum_wx(:,:,:,mymod(i,tpts/nperiod))+wx1d_sa(:,:,:,i)
    end do
      
    deallocate(utau_sa,stat1d_sa);deallocate(sqar1d_sa)
    if (iquad == 1) deallocate(quad1d_sa)
    if (ivort == 1) deallocate(omega1d_sa,wx1d_sa)
    if (itke == 1) deallocate(tke1d_sa)	  

    tpts=tpts/nperiod
    allocate(utau_sa(tpts+1),stat1d_sa(-1:ny+2,4,4,tpts+1),sqar1d_sa(-1:ny+2,4,4,tpts+1))
    if (iquad == 1) allocate(quad1d_sa(-1:ny+2,9,tpts+1))
    if (ivort == 1) allocate(omega1d_sa(-1:ny+2,6,tpts+1),wx1d_sa(-1:ny+2,2,8,tpts+1))
    if (itke == 1) allocate(tke1d_sa(-1:ny+2,4,7,tpts+1))

    utau_sa(1:tpts)=dum_utau(1:tpts)/real(nperiod,kind=mytype)
    stat1d_sa=dum_stat/real(nperiod,kind=mytype)
    sqar1d_sa=dum_sqar/real(nperiod,kind=mytype)
    if (iquad == 1) quad1d_sa=dum_quad/real(nperiod,kind=mytype)
    if (ivort == 1) omega1d_sa=dum_omega/real(nperiod,kind=mytype)
    if (itke == 1) tke1d_sa=dum_tke/real(nperiod,kind=mytype)
    if (ivort == 1) wx1d_sa=dum_wx/real(nperiod,kind=mytype)

    !periodic
    utau_sa(tpts+1)=utau_sa(1)
    stat1d_sa(:,:,:,tpts+1)=stat1d_sa(:,:,:,1)
    sqar1d_sa(:,:,:,tpts+1)=sqar1d_sa(:,:,:,1)
    quad1d_sa(:,:,tpts+1)=quad1d_sa(:,:,1)
    omega1d_sa(:,:,tpts+1)=omega1d_sa(:,:,1)
    wx1d_sa(:,:,:,tpts+1)=wx1d_sa(:,:,:,1)
    tke1d_sa(:,:,:,tpts+1)=tke1d_sa(:,:,:,1)
    tpts=tpts+1
    deallocate(dum_utau,dum_stat,dum_sqar,dum_quad,dum_omega,dum_tke,dum_wx)
  end if
  return
end subroutine phase_ins_1d

!*******************************************************************************

subroutine header_skip(i)
! This reads in the header of the binary file with unit number i.
  character :: dummy*70
  integer :: i
	
  read(i) dummy
  read(i) dummy
  read(i) dummy
  read(i) dummy
  read(i) dummy

  return
end subroutine header_skip

!*******************************************************************************

subroutine change_dir(file1,i)
! Changes the directory to correct simulation/run
  character :: i_char*3,file1*40
  integer :: i

  write(i_char,'(I3)') i
  i_char = adjustl(i_char)
  write(simno,*) simit
  simno=adjustl(simno)
  if (isim1.eq.1) then	
    file1='../run'//trim(i_char)//'/'//trim(file1)
  else if (isim1.ne.1) then	
    file1='../sim'//trim(simno)//'/run'//trim(i_char)//'/'//trim(file1)
  end if

  return
end subroutine change_dir

!*******************************************************************************

subroutine below(iarray,nsize,val,simit,nsize_act,bel)
  implicit none

  integer :: nsize,j,simit,nsize_act,bel
  real(mytype) :: val
  real(mytype), dimension(nsize) :: iarray

  bel=1
  do j=nsize_act,1,-1
    if (iarray(j).LE.val) then; bel=j; exit; end if 
  end do

  return
end subroutine below

!*******************************************************************************

subroutine above(iarray,nsize,val,simit,nsize_act,ab)
  implicit none

  integer :: nsize,j,simit,nsize_act,ab
  real(mytype) :: val
  real(mytype), dimension(nsize) :: iarray

  ab=nsize_act!iarray(nsize,simit)
  do j=1,nsize_act
    if (iarray(j).ge.val) then ; ab=j; exit; end if 
  end do

  return
end subroutine above

!*******************************************************************************

subroutine minm(iarray,nsize,flag,minim)
  implicit none
  integer :: nsize,j,flag,nsize_act
  real(mytype) :: minim
  real(mytype), dimension(nsize,nsims) :: iarray

  if (flag==0) then
    minim=iarray(1,1)
    do j=2,nsims
      if (minim.LT.iarray(1,j)) minim=iarray(1,j)
    end do
  else 
  minim=iarray(step_tot(1),1)
    do j=2,nsims
      write(*,*) 'time-check',j,iarray(step_tot(j),j),step_tot(j),time(step_tot(j),j),nsize
      if (minim.GT.iarray(step_tot(j),j)) minim=iarray(step_tot(j),j)
    end do
  end if

  return
end subroutine minm

!*******************************************************************************

integer function maxm(iarray,nsize)
  integer :: nsize,j
  integer, dimension(nsize) :: iarray

  maxm=iarray(1)
  do j=2,nsize
    if (maxm.LT.iarray(j)) maxm=iarray(j)
  end do

  return
end function maxm

!*******************************************************************************

integer function mymod(n,m)
  integer :: n,m

  mymod=modulo(n,m)
  if (mymod==0) mymod=m     

  return
end function mymod

!*******************************************************************************

subroutine vorty(stat_in,vort_out)
!	This subroutine calculates the fluctuations of different orders in the u direction
  integer :: j,m,ib
  real(mytype) :: denoy,cy1,cy2,cy3
  real(mytype), dimension(0:ny+2,4,4), intent(in) :: stat_in
  real(mytype), dimension(1:nyt,3), intent(inout) :: vort_out

! We then calculate the mean vorticity values
  ib=1
  vort_out=0.0e0

  do j=1,ny
    vort_out(j,1) = vort_out(j,1) + (gcdy_c(j,2,2,ib)*stat_in(j+1,3,1)+ &
                                     gcdy_c(j,3,2,ib)*stat_in(j,3,1)  + &
                                     gcdy_c(j,4,2,ib)*stat_in(j-1,3,1))
    vort_out(j,3) = vort_out(j,3) - (gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)+ &
                                     gcdy_c(j,3,2,ib)*stat_in(j,1,1)  + &
                                     gcdy_c(j,4,2,ib)*stat_in(j-1,1,1))
  end do

       ! vort_out(1,1)=vort_out(1,1)+(stat_in(2,3,1)-stat_in(1,3,1))/dy(1)
        !vort_out(1,3)=vort_out(1,3)-(stat_in(2,1,1)-stat_in(1,1,1))/dy(1)
        !vort_out(ny,1)=vort_out(ny,1)+(stat_in(ny,3,1)-stat_in(ny-1,3,1))/dy(ny)
        !vort_out(ny,3)=vort_out(ny,3)-(stat_in(ny,1,1)-stat_in(ny-1,1,1))/dy(ny)

  return
end subroutine vorty

!*******************************************************************************
!
  subroutine raw_avg(raw,isim,namefile)
! do raw average between tmin and tmax
!  Implicit None
  Integer, Intent(in) :: isim
  Real(mytype), Dimension(:), Intent(in) :: raw
  Real(mytype), Allocatable, Dimension(:) :: raw_av
  Real(mytype) :: sum1
  Integer :: i,j,k
  Character*60 :: namefile

  Allocate(raw_av(1:nmax(isim)-nmin(isim)+1))
  Do i=nmax(isim),nmin(isim),-1
    j=nmax(isim)-i+1
    raw_av(j)=Sum(raw(i:nmax(isim)))/Real(j,kind=mytype)
  End Do
  Open(160,File=namefile)
  Do j=1,nmax(isim)-nmin(isim)+1
    Write(160,'(2F20.10)') time(nmax(isim)-j+1,isim),raw_av(j)
  End Do
  Close(160)
  Deallocate(raw_av)
  return
end subroutine raw_avg
!*******************************************************************************

end module read_post
