module shared_post
  implicit none

!*******************************************************************************
!  Shared Data
#ifdef DOUBLE_PREC
  integer, parameter, public :: mytype = Selected_real_kind( 12, 70 )
#else
  integer, parameter, public :: mytype = Selected_real_kind( 6, 35 )
#endif
!
! Read from input-stat
  integer :: isteady,tstart0,read_3d,iflucs,itke,icorr,ishear,ianis,istruct,tcplt,y_half,ibm,ny_solid, &
		   ivort,iquad,read_2d,nsims,simit,sim_print,tpts,re_print,y_print,half_av, &
		   if_les,max_tot_mod,ntsamps,uns1d,tkemn,spec_loc,jdec,hybrid,iangle,iraw,iplot3d,&
iplot2d,iplot1d,iconstruc,itec
  real(mytype) :: tmin,tmax,reinit,tperiod,isim1,simtot,retau_val,re_val,tperiod_plus,omega_plus,k_plus, &
	           kx_plus,kz_plus,wavelen, wavelen_plus
  integer, dimension(:), ALLOCATABLE :: nruns,nstart,nstats
  character :: simno*20

! Read from stat-deck (run1)
  integer :: nx,ny,nz,itmavg,ntmavg,samples,lpstat_1d,lpstat_2d,lpstat_x, &
             lpstat_z,ibin_rite,lp_ins_1d,lp_ins_2d,lp_ins_3d,lp_snap_1d, &
             lp_snap_2d,lp_snap_3d,lp_snap_x,lp_snap_y,lp_snap_z,lp_snap_xz, &
             lp_snap_yz,lp_snap_xy,nsmpl_ins_1d,nsmpl_ins_2d,nsmpl_ins_3d, &
             nsmpl_snap_1d,nsmpl_snap_2d,nsmpl_snap_3d,nxt,nyt,nzt,nperiod
  real(mytype) ::  cony,length,width,height,nu,small,rho

! Maximums for array allocation
  integer :: max_tot,max_v,snap_1d_max,snap_2d_max,snap_3d_max,ins_1d_max,max_xz

! Read nstep info from all stat-decks and finds total steps to average over
  integer, dimension(:,:), allocatable :: nsteps,nsample,ins_1d,ins_2d,ins_3d, &
					      snap_1d,snap_2d,snap_3d
  integer, dimension(:), allocatable :: step_tot,samp_tot,tot_v,ins_1d_tot,ins_2d_tot, & 
                                        ins_3d_tot,snap_1d_tot,snap_2d_tot,snap_3d_tot,nmin,nmax
! Read 0d data
  real(mytype), dimension(:,:), allocatable :: time,dt,utau,utau1,utau2,dpdx_mean,u_bulk, &
                                               time_v,dt_v,utau_v,u_bulk_v,esum_v,energy_v, &
                                               retau,retau_v,re,re_v,uplus,uplus_v

! Grid Setup
  real(mytype), dimension(:), allocatable :: utau_av,utau1_av,utau2_av,retau_av,re_av,ub_av,Cf
  real(mytype), dimension(:), allocatable :: x,y,z,xc,yc,zc,xplus,yplus,zplus,xcplus, &
                                             ycplus,zcplus,dx,dy,dz,dxs,dys,dzs
  real(mytype), dimension(:,:), allocatable :: xplus_sa,yplus_sa,zplus_sa,xcplus_sa, &
                                               ycplus_sa,zcplus_sa

! Read 1d avg data
  real(mytype), dimension(:,:,:), allocatable :: stat_avg,tke_avg,wx_avg,flucs,kinen
  real(mytype), dimension(:,:), allocatable :: quad_avg,omega_avg,omegamn_avg,shear,anisb,angle
  real(mytype), dimension(:), allocatable :: anisII,anisIII,anisf,anisg

! Read 2d avg data
  real(mytype), dimension(:,:,:,:), allocatable :: stat2d_avg,flucs2d_avg,wx2d_avg,tke2d_avg
  real(mytype), dimension(:,:,:), allocatable :: quad2d_avg,omega2d_avg,shear2d_avg,angle2d_avg
  integer :: nxoz,nxozt

! Read 1d instantaneous data
  real(mytype), dimension(:,:,:,:), allocatable :: stat1d_sa,sqar1d_sa,flucs1d_sa,wx1d_sa,tke1d_sa,kinen1d_sa
  real(mytype), dimension(:,:,:), allocatable :: quad1d_sa,omega1d_sa,omegamn1d_sa,shear1d_sa,anisb1d_sa,angle1d_sa
  real(mytype), dimension(:,:), allocatable :: anisII1d_sa,anisIII1d_sa,anisf1d_sa,anisg1d_sa

! Read 2d instantaneous data
  real(mytype), dimension(:,:,:,:,:), allocatable :: stat2d_sa,sqar2d_sa,flucs2d_sa,wx2d_sa,tke2d_sa,kinen2d_sa
  real(mytype), dimension(:,:,:,:), allocatable :: quad2d_sa,omega2d_sa,omegamn2d_sa,shear2d_sa,anisb2d_sa,angle2d_sa
  real(mytype), dimension(:,:,:), allocatable :: anisII2d_sa,anisIII2d_sa,anisf2d_sa,anisg2d_sa

! Spectral and corrolation calculations
  real(mytype), dimension(:,:,:), allocatable :: corr,spec
  real(mytype) :: start_time

! Write Variables
  character, dimension(40,16) :: varbls
  integer, dimension(9) :: prnt_val
  character (LEN=5), dimension(9) :: prnt_nm
  integer :: jtot,jtott,jmin,jmax,jmaxt

! Simulation averaging variables
  integer :: smin,smax,overall_steps,sims_used
  real(mytype), dimension(:), allocatable :: time_sa,re_sa,rev_sa,utau_sa,dpdx_sa,retau_sa,u_bulk_sa
  real(mytype):: re_avg_sa,retau_avg_sa,utau_avg_sa,ub_avg_sa,Cf_sa

! 3d Visualisation Variables
  integer :: sim3d,ivort3d,istreak,jvis,tvis,ilamq,restrt,jspec,ifilt,iwrite,cnt_stop,iosci,phase,islicexyz, &
	     idecom,isuper,irlarge,icontr,ispec3d,kx_large,kz_large,idef_sup,imodu,jmodu
  integer :: iensem,itop,ibot,ifile1,ifile2,struc_num,str_num,cnt_pos,cnt_neg,xw,yw,zw,iwxfd,iterstr
  character (LEN=60) :: restrt_file, unit1, unitmp, file2davg
  Real(mytype) :: corr_thd,strl_thd, strd_thd
  Real(mytype), Dimension(:,:,:,:), allocatable :: struc_avg_pos,struc_avg_neg
  Integer, dimension(:,:,:), allocatable  :: centres
  Integer, dimension(:,:), allocatable :: strucdim,strucshift
  Integer, dimension(:), allocatable :: struc_len,vort_dir,minx,xcen

! Coordinate transformation
  Integer :: itrans
  Real(mytype) :: theta

!*******************************************************************************

! Main simulation parameters
  integer :: iconbc,ibwest,ibeast,ibsoth,ibnoth,ngr,nbl,kbsoth,kbnoth,  &
             ibmpfg,igr_2pt,ibin_read,idomain,ifftw,isgsm,jconbc,kconbc,&
             ichann,mpi,if2d_noy

  real(mytype) :: roh,myu,dpmean,time0,pi,pi2,pi4,large,very_small
  real(mytype), dimension(:,:,:), allocatable ::grdsxx,grdsyy,grdszz,ap,ae,as,at
  real(mytype), dimension(:,:,:,:), allocatable :: gcdx_f,gcdy_f,gcdz_f,gcdx_c,gcdy_c, &
                                                   gcdz_f2,gfdx_c,gfdy_c,gfdz_c,gfdz_c2,gcdz_c

! Block variables
  integer::ibegbl,iendbl,jbegbl,jendbl,kbegbl,kendbl

! FFTW parameters
  real(mytype), dimension(:), allocatable :: infft1
  complex, dimension(:), allocatable :: outfft1,infft2,outfft2
  integer*8 :: planzf,planzb,planxf,planxb
  

end module shared_post
