	MODULE calc_super
!	USE shared_data
	USE shared_post
	USE cal_coefficient
	USE read_post
	IMPLICIT NONE
	Real(mytype), Dimension(:,:), Allocatable :: uu,ul,spec2dxz,vv
	Real(mytype) :: time_snap
	Character :: file2d*40,file3d*40
	CONTAINS
!*******************************************************************************
! This code performs the calculations on the values that have been read in 
! from the read_post module.
!*******************************************************************************
	SUBROUTINE struct_super
        Real(mytype), Parameter :: kx=3.0e0,kz=2.0e0
        Integer :: nsnaps,tp,n,l,m,cor_cntx,cor_cor_cntxcntz,cor_cntz,lam_cnt
	Integer :: ivisu,ishuffle,ihht,ihht2d
	Integer :: i,j,k
	Real(mytype) :: dumm
	Real(mytype), Dimension(:,:), Allocatable :: dummarr
!
	simit=1;ivisu=0;ishuffle=0;ihht=0;ihht2d=1
!	Initialisations
	Allocate(uu(-1:nx+2,-1:nz+2),vv(-1:nx+2,-1:nz+2))
!
	If(ispec3d == 1) Call read_3dspectra
!
	If(icontr == 1) Then
	  idecom=1; irlarge=0
	  If(ishear.Ne.1) Then
	    Write(*,*) 'ERROR***: ishear is not 1 for wall shearstress decomposition'
	  Stop; End If
	End If
!	read 2D xz plane
	Allocate(spec2dxz(1:nx/2+1,1:nz)); spec2dxz(1:nx/2+1,1:nz)=0.0e0
	Allocate(dummarr(-1:nx+2,-1:nz+2))
	DO n = nstart(simit), nruns(simit)
	If(ishear == 1) Then
!	file2d='OUTPUT_2D_SHEAR.dat'
	file2d='OUTPUT_2D_u_yp15.dat'
!	file2d='OUTPUT_2D_u_yh01.dat'
	Else If(lp_snap_2d == 1) Then
	file2d='OUTPUT_2D_SNAP.dat'
	Else; Goto 1112
	End IF
	Write(*,*) 'READING...', file2d
	CALL change_dir(file2d,n)
	Open(507, File=file2d, Status='OLD', Form='UNFormatTED', access='stream')
	Do l=1,1 !nsmpl_snap_2d
          Read(507) time_snap
	  Write(*,*) time_snap
	  Read(507) uu(:,:) !u
	  If(lp_snap_2d==1) Then
	    Read(507) vv(:,:) !v
	    Read(507) dummarr(:,:) !w
	    Read(507) dummarr(:,:) !p
	  End If
!         calculate averaged spectra
!	  Call spectraxzavg
	  Call spectraxz_snap(l)
! 	  visualise the dataset
	  If(ivisu==1) Call write_vtk(uu(1:nx,1:nz))
	  If(ishuffle==1) Call signal_shuffle(uu)
!
	  If(ihht == 1) Call HHT(uu)
	  If(ihht2d == 1) Call HHT2D(uu)
	End Do
!
1112	Continue
	End Do
! artificial wave
!	uu(:,:)=0.0e0
!	Do k=1,nz; Do i=1,nx
!	  uu(i,k)=sin(kx*2.0e0*pi*x(i)/length+kz*2.0e0*pi*z(k)/width)+sin(3.0e0*kx*2.0e0*pi*x(i)/length+3.0e0*kz*2.0e0*pi*z(k)/width)
!	End Do; End Do

!	Write to tecplot file
	file2d='super_2d.plt'
	Open(508, File=file2d, status='UNKNOWN')
	If(ishear == 1) Then
        WRITE(508,'(320A)') ' variables="x","z","<greek>t</greek><sub>w</sub>"'
	Else
        WRITE(508,'(320A)') ' variables="x","z","u"'
	End If
 	WRITE(508,*) 'zone t="',dumm,'" i=',nx+2,', j=',nz+2,'f=point'
!	Do k=0,nz+1; Do i=0,nx+1
!	  Write(508,'(3f20.10)') x(i),z(k),uu(i,k)
!	End Do; End Do
        Close(508)
!
        Write(*,*) 'begin to calc 2D spectra'
        If(idecom == 1) Call spectra2D
!
	If(icontr == 1) Call percent_shear
!
	Deallocate(uu,vv,dummarr)
	END SUBROUTINE struct_super

!*******************************************************************************
        Subroutine spectra2D
  	Include "fftw3.f"
!plot 2D spectra in xz plane
        Real(mytype), Allocatable, Dimension(:,:) :: spec2d
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxz
        Integer :: i,j,k,ip,jp,kp,kz
        Real(mytype) :: dumu,specavg
	Real(mytype), Allocatable, Dimension(:) :: accum,diff,lamx,lamz,kapx,kapz
	Real(mytype), Allocatable, Dimension(:,:) :: accum2,diff2
        Integer*8 :: planxz
!
	Allocate(spec2d(1:nx/2+1,1:nz),outfftxz(nx/2+1,nz),ul(1:nx,1:nz))
!
	Allocate(accum(1:Max(nx/2+1,nz/2+1)),diff(1:Max(nx/2+1,nz/2+1)),&
	lamx(1:nx/2+1),kapx(1:nx/2+1), &
	lamz(1:nz/2+1),kapz(1:nz/2+1))
	Allocate(accum2(1:nx/2+1,1:nz/2+1),diff2(1:nx/2+1,1:nz/2+1))
	accum(:)=0.0e0; lamx(:)=0.0e0; kapx(:)=0.0e0; lamz(:)=0.0e0; kapz(:)=0.0e0
!
!
	lamx(1)=0.0e0; kapx(1)=0.0e0
	Do i=2,nx/2+1; lamx(i)=length/real(i-1,mytype); kapx(i)=2.0e0*pi/lamx(i); End Do	
	lamz(1)=0.0e0; kapz(1)=0.0e0
	Do k=2,nz/2+1; lamz(k)=width/real(k-1,mytype);kapz(k)=2.0e0*pi/lamz(k); End Do
!	  
        Write(*,*) 'nx=',nx,'nz=', nz    
	ul(1:nx,1:nz)=uu(1:nx,1:nz)
	Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ul,outfftxz,FFTW_ESTIMATE)
	Call dfftw_execute(planxz)
!
	spec2d(1:nx/2+1,1:nz)=cdabs(outfftxz(1:nx/2+1,1:nz))**2/Real(nx*nz,mytype)/Real(nx*nz,mytype)
	!write averaged spectra
	spec2d(1:nx/2+1,1:nz)=spec2dxz(1:nx/2+1,1:nz)/Real(snap_2d_tot(simit),mytype)
!
	Open(Unit=314,File='spectra_2dp.plt')
	Open(Unit=315,File='spectra_2d.plt')
	Write(314,'(10A)') ' variables="k<sub>x</sub><sup>+</sup>","k<sub>z</sub><sup>+</sup>","y<sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk/u<sub><greek>t</greek></sub><sup>2</sup>", "<greek>F</greek><sub>uu</sub>/<greek>n</greek><sup>2</sup>","k<sub>x</sub>k<sub>z</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"'
	Write(315,'(10A)') ' variables="k<sub>x</sub>","k<sub>z</sub>","y",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk", "<greek>F</greek><sub>uu</sub>","k<sub>x</sub>k<sub>z</sub><greek>F</greek><sub>uu</sub>"'
!
        WRITE(314,*) 'zone i=',nx/2+1,', j=',nz/2+1,',f=point'
        WRITE(315,*) 'zone i=',nx/2+1,', j=',nz/2+1,',f=point'
	Call specxz(spec2d,kapx,kapz,accum2,diff2)
        DO k=1,nz/2+1;DO i=1,nx/2+1
	  Write(314,'(7E15.7)') kapx(i)/retau_avg_sa, kapz(k)/retau_avg_sa, yc(jp)*retau_avg_sa, &
	accum2(i,k)/utau_avg_sa**2, diff2(i,k)/(nu*nu), diff2(i,k)*kapx(i)*kapz(k)/utau_avg_sa**2
	  Write(315,'(7E15.7)') kapx(i),kapz(k),yc(jp),accum2(i,k),diff2(i,k),diff2(i,k)*kapx(i)*kapz(k)
        END DO; END DO
        Close(314); Close(315)
!
	Open(unit=348,file='spectra_snap_kxp.plt') !squashed in z direction
	Open(unit=349,file='spectra_snap_kx.plt') !squashed in z direction
	Write(348,'(10A)') ' variables="k<sub>x</sub><sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk/u<sub><greek>t</greek></sub><sup>2</sup>", "<greek>F</greek><sub>uu</sub>/(u<sub><greek>t</greek></sub><greek>n</greek>)","k<sub>x</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"'
	Write(349,'(10A)') ' variables="k<sub>x</sub>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk", "<greek>F</greek><sub>uu</sub>","k<sub>x</sub><greek>F</greek><sub>uu</sub>"'
!
	Open(unit=350,file='spectra_snap_kzp.plt') !squashed in x direction
	Open(unit=351,file='spectra_snap_kz.plt') !squashed in x direction
	Write(350,'(10A)') ' variables="k<sub>z</sub><sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk/u<sub><greek>t</greek></sub><sup>2</sup>", "<greek>F</greek><sub>uu</sub>/(u<sub><greek>t</greek></sub><greek>n</greek>)","k<sub>z</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"'
	Write(351,'(10A)') ' variables="k<sub>z</sub>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk", "<greek>F</greek><sub>uu</sub>","k<sub>z</sub><greek>F</greek><sub>uu</sub>"'
!
	Write(348,*) 'zone i=',nx/2+1,',f=point'
	Write(349,*) 'zone i=',nx/2+1,',f=point'
	Call specxy(spec2d,kapx,accum,diff)
!
	Do i=1,nx/2+1
	  Write(348,'(7E15.7)') kapx(i)/retau_avg_sa, &
	accum(i)/utau_avg_sa**2, diff(i)/(utau_avg_sa*nu), diff(i)*kapx(i)/utau_avg_sa**2
	  Write(349,'(7E15.7)') kapx(i), &
	accum(i), diff(i), diff(i)*kapx(i)
	End Do
!
	Write(350,*) 'zone i=',nz/2+1,',f=point'
	Write(351,*) 'zone i=',nz/2+1,',f=point'
	Call specyz(spec2d,kapz,accum,diff)
!
	Do k=1,nz/2+1
  	  Write(350,'(7E15.7)') kapz(k)/retau_avg_sa, &
	 accum(k)/utau_avg_sa**2, diff(k)/(utau_avg_sa*nu), diff(k)*kapz(k)/utau_avg_sa**2
  	  Write(351,'(7E15.7)') kapz(k), &
	 accum(k), diff(k), diff(k)*kapz(k)
	End Do

	Close(348); Close(349); Close(350); Close(351)
!
!	remove large scales
	Do k=1,nz; If(k<nz/2+1) Then; kz=k; Else; kz=nz-k+2; End If
	  Do i=1,nx/2+1
	  If((idef_sup.Eq.1 .And. (kz.Le.(kz_large+1)) .Or. &
	    (idef_sup.Eq.2 .And. i.Le.(kx_large+1))) .Or. &
	    (idef_sup.Eq.3 .And. kz.Le.(kz_large+1) .And. i.Le.(kx_large+1) ) .Or. &
	    (idef_sup.Eq.4 .And. (kz.Le.(kz_large+1) .Or. i.Le.(kx_large+1)) ) ) Then !kz<=kz_large
	    outfftxz(i,k)=0.0e0
	  End If
	End Do; End Do
! inverse fft
	outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
	Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz,ul,FFTW_ESTIMATE)
	Call dfftw_execute(planxz)
	ul(1:nx,1:nz)=uu(1:nx,1:nz)-Sum(uu(1:nx,1:nz))/Real(nx*nz,mytype)-ul(1:nx,1:nz)
!	Write to tecplot file
	file2d='super_2d_inverse.plt'
	Open(316, File=file2d, status='UNKNOWN')
	If(ishear == 1) Then
        WRITE(316,'(320A)') ' variables="x","z","<greek>t</greek><sub>w</sub>"'
	Else
        WRITE(316,'(320A)') ' variables="x","z","u"'
	End If
 	WRITE(316,*) 'zone t="fft-1" i=',nx,', j=',nz,'f=point'
	Do k=1,nz; Do i=1,nx
	  If(irlarge == 0) Then; Write(316,'(3f20.10)') x(i),z(k),ul(i,k)
	  Else; Write(316,'(3f20.10)') x(i),z(k),uu(1:nx,1:nz)-Sum(uu(1:nx,1:nz))/Real(nx*nz,mytype)-ul(1:nx,1:nz); End If
	End Do; End Do
        Close(316)
!
	Deallocate(spec2d,ul)
	Deallocate(accum,diff,lamx,kapx,lamz,kapz,accum2,diff2)
!
        End Subroutine spectra2D
!*******************************************************************************
        Subroutine spectraxzavg
  	Include "fftw3.f"
!plot 2D spectra in xz plane
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxz,outfftxzvv
!
        Integer :: i,j,k,ip,jp,kp,iuv
        Real(mytype) :: dumu, wl,wlx,wlz, wn,wnx,wnz
        Integer*8 :: planxz
	Allocate(outfftxz(1:nx/2+1,1:nz),outfftxzvv(1:nx/2+1,1:nz),ul(1:nx,1:nz))
!  
	iuv=1 ! 1- get spectra for uv; 0- get spectra for uu
        Write(*,*) 'nx=',nx,'nz=', nz    
!
	ul(1:nx,1:nz)=uu(1:nx,1:nz)
	Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ul,outfftxz,FFTW_ESTIMATE)
	Call dfftw_execute(planxz)
	If(iuv == 1) Then
	  ul(1:nx,1:nz)=vv(1:nx,1:nz)
	  Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ul,outfftxzvv,FFTW_ESTIMATE)
	  Call dfftw_execute(planxz)
	  Do k=1,nz; Do i=1,nx/2+1
	    spec2dxz(i,k)=spec2dxz(i,k)-Real(Conjg(outfftxz(i,k))*outfftxzvv(i,k))/Real(nx*nz,mytype)/Real(nx*nz,mytype)
	  End Do; End Do
	Else
	  spec2dxz(1:nx/2+1,1:nz)=spec2dxz(1:nx/2+1,1:nz)+Conjg(outfftxz(1:nx/2+1,1:nz))*outfftxz(1:nx/2+1,1:nz) &
	  /Real(nx*nz,mytype)/Real(nx*nz,mytype)
	End If
	Deallocate(ul,outfftxz,outfftxzvv)
!
        End Subroutine spectraxzavg
!*******************************************************************************
        Subroutine spectraxz_snap(l)
  	Include "fftw3.f"
!plot 2D spectra in xz plane
	Integer, Intent(In) :: l
        Real(mytype), Allocatable, Dimension(:,:) :: spec2d
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxz
        Integer :: i,j,k,ip,jp,kp,kz
        Real(mytype) :: dumu,specavg
	Real(mytype), Allocatable, Dimension(:) :: accum,diff,lamx,lamz,kapx,kapz
	Real(mytype), Allocatable, Dimension(:,:) :: accum2,diff2
        Integer*8 :: planxz
!
	Allocate(spec2d(1:nx/2+1,1:nz),outfftxz(nx/2+1,nz),ul(1:nx,1:nz))
!
	Allocate(accum(1:Max(nx/2+1,nz/2+1)),diff(1:Max(nx/2+1,nz/2+1)),&
	lamx(1:nx/2+1),kapx(1:nx/2+1), &
	lamz(1:nz/2+1),kapz(1:nz/2+1))
	Allocate(accum2(1:nx/2+1,1:nz/2+1),diff2(1:nx/2+1,1:nz/2+1))
	accum(:)=0.0e0; lamx(:)=0.0e0; kapx(:)=0.0e0; lamz(:)=0.0e0; kapz(:)=0.0e0
!
!
	lamx(1)=0.0e0; kapx(1)=0.0e0
	Do i=2,nx/2+1; lamx(i)=length/real(i-1,mytype); kapx(i)=2.0e0*pi/lamx(i); End Do	
	lamz(1)=0.0e0; kapz(1)=0.0e0
	Do k=2,nz/2+1; lamz(k)=width/real(k-1,mytype);kapz(k)=2.0e0*pi/lamz(k); End Do
!	  
        Write(*,*) 'nx=',nx,'nz=', nz    
	ul(1:nx,1:nz)=uu(1:nx,1:nz)
	Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ul,outfftxz,FFTW_ESTIMATE)
	Call dfftw_execute(planxz)
!
	spec2d(1:nx/2+1,1:nz)=cdabs(outfftxz(1:nx/2+1,1:nz))**2/Real(nx*nz,mytype)/Real(nx*nz,mytype)
!
	If(l==1) Then
	Open(unit=348,file='spectra_snap_kxp.plt') !squashed in z direction
	Open(unit=349,file='spectra_snap_kx.plt') !squashed in z direction
	Write(348,'(11A)') ' variables="k<sub>x</sub><sup>+</sup>","t<sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk/u<sub><greek>t</greek></sub><sup>2</sup>", "<greek>F</greek><sub>uu</sub>/(u<sub><greek>t</greek></sub><greek>n</greek>)","k<sub>x</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"'
	Write(349,'(11A)') ' variables="k<sub>x</sub>","t",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk", "<greek>F</greek><sub>uu</sub>","k<sub>x</sub><greek>F</greek><sub>uu</sub>"'
!
	Open(unit=350,file='spectra_snap_kzp.plt') !squashed in x direction
	Open(unit=351,file='spectra_snap_kz.plt') !squashed in x direction
	Write(350,'(11A)') ' variables="k<sub>z</sub><sup>+</sup>","t<sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk/u<sub><greek>t</greek></sub><sup>2</sup>", "<greek>F</greek><sub>uu</sub>/(u<sub><greek>t</greek></sub><greek>n</greek>)","k<sub>z</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"'
	Write(351,'(11A)') ' variables="k<sub>z</sub>","t",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk", "<greek>F</greek><sub>uu</sub>","k<sub>z</sub><greek>F</greek><sub>uu</sub>"'
!
	Write(348,*) 'zone i=',nx/2+1, ', j=',nsmpl_snap_2d,',f=point'
	Write(349,*) 'zone i=',nx/2+1, ', j=',nsmpl_snap_2d,',f=point'
	Write(350,*) 'zone i=',nz/2+1, ', j=',nsmpl_snap_2d,',f=point'
	Write(351,*) 'zone i=',nz/2+1, ', j=',nsmpl_snap_2d,',f=point'
	End If
!
	Call specxy(spec2d,kapx,accum,diff)
!
	Do i=1,nx/2+1
	  Write(348,'(8E15.7)') kapx(i)/retau_avg_sa, time_snap*retau_avg_sa**2*nu, &
	accum(i)/utau_avg_sa**2, diff(i)/(utau_avg_sa*nu), diff(i)*kapx(i)/utau_avg_sa**2
	  Write(349,'(8E15.7)') kapx(i), time_snap, &
	accum(i), diff(i), diff(i)*kapx(i)
	End Do
!
	Call specyz(spec2d,kapz,accum,diff)
!
	Do k=1,nz/2+1
  	  Write(350,'(8E15.7)') kapz(k)/retau_avg_sa, time_snap*retau_avg_sa**2*nu, &
	 accum(k)/utau_avg_sa**2, diff(k)/(utau_avg_sa*nu), diff(k)*kapz(k)/utau_avg_sa**2
  	  Write(351,'(8E15.7)') kapz(k), time_snap, &
	 accum(k), diff(k), diff(k)*kapz(k)
	End Do

	If(l==nsmpl_snap_2d) Then; Close(348); Close(349); Close(350); Close(351); End If
!
	Deallocate(spec2d,ul)
	Deallocate(accum,diff,lamx,kapx,lamz,kapz,accum2,diff2)
        End Subroutine spectraxz_snap
!*******************************************************************************
        Subroutine percent_shear
	Real(mytype) :: cf_total,cf_pos,cf_neg,cf_laminar
	Integer :: i,k,cnt_pos,cnt_neg,irlaminar
!
!	remove the larminar compoent from Cf
	irlaminar=1; cf_laminar=6.0e0/re_avg_sa
	If(irlaminar == 1) Then
	  Write(*,*) 'Laminar component of Cf has been removed!'
	  uu(:,:)=uu(:,:)-cf_laminar/2.0e0
	End If
!	total wall shear stress
	cf_total=Sum(uu(1:nx,1:nz))
	cf_total=cf_total/Real(nx*nz,mytype)
	cf_total=cf_total*2.0e0
!	contribution from positive super scale structure
	cf_pos=0.0e0; cf_neg=0.0e0; cnt_pos=0; cnt_neg=0
	Do k=1,nz; Do i=1,nx
	  If(ul(i,k).Ge.0.0e0) Then
	    cf_pos=cf_pos+uu(i,k)
	    cnt_pos=cnt_pos+1
	  Else
	    cf_neg=cf_neg+uu(i,k)
	    cnt_neg=cnt_neg+1
	  End If
	End Do; End Do
	cf_pos=cf_pos/Real(cnt_pos,mytype); cf_neg=cf_neg/Real(cnt_neg,mytype)
	cf_pos=cf_pos*2.0e0; cf_neg=cf_neg*2.0e0
!
	Write(*,'("*****************************************")')
	Write(*,'(" Cf_total=",E12.5," Cf_pos=",E12.5," Cf_neg=",E12.5)') &
		cf_total,cf_pos,cf_neg
	Write(*,'(" %_pos=",F10.5," %_neg=",F10.5)') &
		cf_pos*Real(cnt_pos,mytype)/cf_total/Real(nx*nz,mytype)*100.e0, &
		cf_neg*Real(cnt_neg,mytype)/cf_total/Real(nx*nz,mytype)*100.0e0
	Write(*,'("*****************************************")')
!	FIK
	If(ishear == 1 .And. step_tot(simit).Gt.1) Call fik_shear
!
        End Subroutine percent_shear
!*******************************************************************************
        Subroutine fik_shear
	Real(mytype) :: cf_total,cf_pos,cf_neg,cf_laminar
	Real(mytype),  Allocatable, Dimension(:) :: fik,accum
	Integer :: i,j,k,j1
!
!
!	FIK identity
	Allocate(fik(-1:ny+2),accum(-1:ny+2))
	fik(:)=0.e0; accum(:)=0.0e0
	Do j=0,nyt
	  fik(j)=shear(j,1)*(1.0e0-yc(j))
	End Do
!	accumulate
	Do j=0,nyt; Do j1=0,j
	  accum(j)=accum(j)+fik(j1)*dy(j1)
	End Do; End Do
!
!	write output file
	file2d='FIK.plt'
	Open(509, File=file2d, status='UNKNOWN')
        WRITE(509,'(320A)') ' variables="yc","uv","(1-y)uv","accum"'
 	WRITE(509,*) 'zone t="fik" i=',ny+2,'f=point'
	Do j=0,nyt
	  Write(509,'(3f20.10)') yc(j),shear(j,1)/utau_avg_sa**2,fik(j)/utau_avg_sa**2,&
		accum(j)/utau_avg_sa**2
	End Do
        Close(509)

	Deallocate(fik,accum)
        End Subroutine fik_shear
!
!*******************************************************************************
        Subroutine read_3dspectra
!	this subroutine reads 3d spectra file
	Real(mytype),  Allocatable, Dimension(:,:,:) :: spec_sum
	Real(mytype) :: specavg,urms,ulrms
	Integer :: i,j,jp,k,n,nyy,jmin,jmax,nread,kz
	Real(mytype), Allocatable, Dimension(:) :: accum,diff,lamx,lamz,kapx,kapz
	Real(mytype), Allocatable, Dimension(:,:) :: accum2,diff2,spec_tmp
!
	nread=1; jmin=1; jmax=ny!jspec
	Allocate(spec_sum(1:nx/2+1,1:ny,1:nz),spec_tmp(1:nx/2+1,1:ny))
	Allocate(accum2(1:nx/2+1,1:nz/2+1),diff2(1:nx/2+1,1:nz/2+1))
	Allocate(accum(1:Max(nx/2+1,nz/2+1)),diff(1:Max(nx/2+1,nz/2+1)),&
	lamx(1:nx/2+1),kapx(1:nx/2+1), &
	lamz(1:nz/2+1),kapz(1:nz/2+1))
	accum(:)=0.0e0; lamx(:)=0.0e0; kapx(:)=0.0e0; lamz(:)=0.0e0; kapz(:)=0.0e0
!
	lamx(1)=0.0e0; kapx(1)=0.0e0
	Do i=2,nx/2+1; lamx(i)=length/real(i-1,mytype); kapx(i)=2.0e0*pi/lamx(i); End Do	
	lamz(1)=0.0e0; kapz(1)=0.0e0
	Do k=2,nz/2+1; lamz(k)=width/real(k-1,mytype);kapz(k)=2.0e0*pi/lamz(k); End Do
!
	DO n = nstart(simit), nruns(simit)
	file3d='SPECTRA_3D.dat'
	CALL change_dir(file3d,n)
	Open(510, File=file3d, Status='OLD', Form='UNFormatTED',access='stream')
	Write(*,*) 'reading 3d spectra file from',file3d
	Do i=1,nread; If(i==nread) Write(*,*) 'Component ',i; 
	  Do k=1,nz; Read(510) spec_tmp(:,:)
	    If(i==nread) Then; spec_sum(:,:,k)=spec_sum(:,:,k)+spec_tmp(:,:); End If
	  End Do
        End Do; Close(510)
	End Do
	Deallocate(spec_tmp)
!
!	time average
	spec_sum(:,:,:)=(-1.0e0)**(nread==2)*spec_sum(:,:,:)/Real(samp_tot(simit),mytype)/Real(nx*nz,mytype)
	Write(*,*) 'samples=',samp_tot(simit)
!	average top and bottom
	If(half_av == 1) Then
	  nyy=ny/2; 
	  Do j=1,nyy; spec_sum(:,j,:)=(spec_sum(:,j,:)+(-1.0e0)**(nread==2)*spec_sum(:,ny+1-j,:))/2.0e0; End Do
	  Do j=1,nyy; spec_sum(:,ny+1-j,:)=(-1.0e0)**(nread==2)*spec_sum(:,j,:); End Do
	Else; nyy=ny; End If
!
	file3d='spectra_3dp.plt'; Open(unit=601,file=file3d) !3d data
	file3d='spectra_3d.plt'; Open(unit=602,file=file3d)
	Write(601,'(10A)') ' variables="k<sub>x</sub><sup>+</sup>","k<sub>z</sub><sup>+</sup>","y<sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk/u<sub><greek>t</greek></sub><sup>2</sup>", "<greek>F</greek><sub>uu</sub>/<greek>n</greek><sup>2</sup>","k<sub>x</sub>k<sub>z</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"'
	Write(602,'(10A)') ' variables="k<sub>x</sub>","k<sub>z</sub>","y",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk", "<greek>F</greek><sub>uu</sub>","k<sub>x</sub>k<sub>z</sub><greek>F</greek><sub>uu</sub>"'
	Write(601,*) 'zone i=',nx/2+1,', j=',nz/2+1,', k=',jmax-jmin+1,',f=point'
	Write(602,*) 'zone i=',nx/2+1,', j=',nz/2+1,', k=',jmax-jmin+1,',f=point'
!
	Open(unit=607,file='urmsp_decom.plt') !reconstruct the urms profile
	Open(unit=608,file='urms_decom.plt') !reconstruct the urms profile
	Write(607,'(A)') ' variables="y","u<sup>+2</sup>","u<sub>L</sub><sup>+2</sup>","u<sub>S</sub><sup>+2</sup>","y<sup>+</sup>"'
	Write(608,'(A)') ' variables="y","u<sup>2</sup>","u<sub>L</sub><sup>2</sup>","u<sub>S</sub><sup>2</sup>","y<sup>+</sup>"'
!
	Do jp=jmin,jmax
	  Call specxz(spec_sum(:,jp,:),kapx,kapz,accum2,diff2)
!
	  urms=Sum(diff2(:,:))*(2.0e0*pi/length)*(2.0e0*pi/width)
	  Write(*,*) '------Euiuj------', urms,urms/utau_avg_sa**2
	  ulrms=0.0e0
	  Do k=1,nz; If(k<nz/2+1) Then; kz=k; Else; kz=nz-k+2; End If
	    Do i=1,nx/2+1
	    If((idef_sup.Eq.1 .And. (kz.Le.(kz_large+1)) .Or. &
	      (idef_sup.Eq.2 .And. i.Le.(kx_large+1))) .Or. &
	      (idef_sup.Eq.3 .And. kz.Le.(kz_large+1) .And. i.Le.(kx_large+1) ) .Or. &
	      (idef_sup.Eq.4 .And. (kz.Le.(kz_large+1) .Or. i.Le.(kx_large+1)) ) ) Then !kz<=kz_large
	      ulrms=ulrms+diff2(i,k)*(2.0e0*pi/length)*(2.0e0*pi/width)
	    End If
	  End Do; End Do
!
	  DO k=1,nz/2+1; DO i=1,nx/2+1
!	    Write(601,'(7E15.7)') kapx(i)/retau_avg_sa, kapz(k)/retau_avg_sa, yc(jp)*retau_avg_sa, &
!	accum2(i,k)/utau_avg_sa**2, diff2(i,k)/(nu*nu), diff2(i,k)*kapx(i)*kapz(k)/utau_avg_sa**2
!	    Write(602,'(7E15.7)') kapx(i),kapz(k),yc(jp),accum2(i,k),diff2(i,k),diff2(i,k)*kapx(i)*kapz(k)
	  End Do; End Do
!
	  Write(607,'(5F20.10)') y(jp),urms/utau_avg_sa**2,ulrms/utau_avg_sa**2, &
		(urms-ulrms)/utau_avg_sa**2,y(jp)*retau_avg_sa
	  Write(608,'(5F20.10)') y(jp),urms,ulrms, &
		(urms-ulrms),y(jp)*retau_avg_sa
	End Do
!
	Open(unit=603,file='spectra_3t2_kxp.plt') !squashed in z direction
	Open(unit=604,file='spectra_3t2_kx.plt') !squashed in z direction
	Write(603,'(10A)') ' variables="k<sub>x</sub><sup>+</sup>","y<sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk/u<sub><greek>t</greek></sub><sup>2</sup>", "<greek>F</greek><sub>uu</sub>/(u<sub><greek>t</greek></sub><greek>n</greek>)","k<sub>x</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"'
	Write(604,'(10A)') ' variables="k<sub>x</sub>","y",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk", "<greek>F</greek><sub>uu</sub>","k<sub>x</sub><greek>F</greek><sub>uu</sub>"'
	Write(603,*) 'zone i=',nx/2+1,', j=',jmax-jmin+1,',f=point'
	Write(604,*) 'zone i=',nx/2+1,', j=',jmax-jmin+1,',f=point'
!
	Do jp=jmin,jmax
	  Call specxy(spec_sum(:,jp,:),kapx,accum,diff)
	  Do i=1,nx/2+1
	    Write(603,'(7E15.7)') kapx(i)/retau_avg_sa, yc(jp)*retau_avg_sa, &
	    accum(i)/utau_avg_sa**2, diff(i)/(utau_avg_sa*nu), diff(i)*kapx(i)/utau_avg_sa**2
	    Write(604,'(7E15.7)') kapx(i), yc(jp), &
	    accum(i), diff(i), diff(i)*kapx(i)
	  End Do
	End Do
!
	Open(unit=605,file='spectra_3t2_kzp.plt') !squashed in x direction
	Open(unit=606,file='spectra_3t2_kz.plt') !squashed in x direction
	Write(605,'(10A)') ' variables="k<sub>z</sub><sup>+</sup>","y<sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk/u<sub><greek>t</greek></sub><sup>2</sup>", "<greek>F</greek><sub>uu</sub>/(u<sub><greek>t</greek></sub><greek>n</greek>)","k<sub>z</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>"'
	Write(606,'(10A)') ' variables="k<sub>z</sub>","y",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk", "<greek>F</greek><sub>uu</sub>","k<sub>z</sub><greek>F</greek><sub>uu</sub>"'
	Write(605,*) 'zone i=',nz/2+1,', j=',jmax-jmin+1,',f=point'
	Write(606,*) 'zone i=',nz/2+1,', j=',jmax-jmin+1,',f=point'	
!
	Do jp=jmin,jmax
	  Call specyz(spec_sum(:,jp,:),kapz,accum,diff)
	  Do k=1,nz/2+1
  	    Write(605,'(7E15.7)') kapz(k)/retau_avg_sa, yc(jp)*retau_avg_sa, &
	   accum(k)/utau_avg_sa**2, diff(k)/(utau_avg_sa*nu), diff(k)*kapz(k)/utau_avg_sa**2
  	    Write(606,'(7E15.7)') kapz(k), yc(jp), &
	   accum(k), diff(k), diff(k)*kapz(k)
	  End Do
!
	End Do
!
	Close(601); Close(602); Close(603); Close(604); Close(605); Close(606)
	Close(607); Close(608)
!
	Deallocate(spec_sum)
!
        End Subroutine read_3dspectra
!*******************************************************************************
	Subroutine specxy(spec2d,kap,accum,diff)
	!plot 2D spectra in yz plane
	Real(mytype), Dimension(:,:), Intent(In) :: spec2d
	Real(mytype), Dimension(:), Intent(InOut) :: accum,diff,kap
	Real(mytype), Allocatable, Dimension(:) :: euu
        Real(mytype) :: specavg, dkap
	Integer :: i

	Allocate(euu(1:Max(nx/2+1,nz/2+1)))
	accum(:)=0.0e0
	dkap=(kap(2)-kap(1))
!
	accum(1)=Sum(spec2d(1,2:nz))
	euu(1)=accum(1)
	Do i=2,nx/2
	  specavg=(Sum(spec2d(i,1:nz)))*2.0e0
	  euu(i)=specavg
	  accum(i)=accum(i-1)+specavg
	End Do
	accum(nx/2+1)=accum(nx/2)+(Sum(spec2d(nx/2+1,1:nz)))
	euu(nx/2+1)=Sum(spec2d(nx/2+1,1:nz))
!
!	diff(1)=(accum(2)-accum(1))/(kap(2)-kap(1))
!	diff(nx/2+1)=(accum(nx/2+1)-accum(nx/2))/(kap(nx/2+1)-kap(nx/2))
!	Do i=2,nx/2
!	  diff(i)=(accum(i+1)-accum(i-1))/(kap(i+1)-kap(i-1))
!	End Do
!
!	another way to get spectra density funciton
	diff(1:nx/2+1)=euu(1:nx/2+1)/dkap
	Deallocate(euu)
!
	End Subroutine specxy
!*******************************************************************************
	Subroutine specyz(spec2d,kap,accum,diff)
	!plot 2D spectra in yz plane
	Real(mytype), Dimension(:,:), Intent(In) :: spec2d
	Real(mytype), Dimension(:), Intent(InOut) :: accum,diff,kap
	Real(mytype), Allocatable, Dimension(:) :: euu
        Real(mytype) :: specavg,dkap
	Integer :: k

	Allocate(euu(1:Max(nx/2+1,nz/2+1)))
	accum(:)=0.0e0; euu(:)=0.0e0
	dkap=kap(3)-kap(2)
!
        accum(1)=(Sum(spec2d(2:nx/2,1))*2.0e0+spec2d(nx/2+1,1))
	euu(1)=accum(1)
	Do k=2,nz/2
	  specavg=(Sum(spec2d(1:nx/2,k))+Sum(spec2d(2:nx/2+1,nz+2-k)))*2.0e0
	  euu(k)=specavg
	  accum(k)=accum(k-1)+specavg
	End Do
	accum(nz/2+1)=accum(nz/2)+(Sum(spec2d(1:nx/2,nz/2+1))+Sum(spec2d(2:nx/2+1,nz/2+1)))
	euu(nz/2+1)=(Sum(spec2d(1:nx/2,nz/2+1))+Sum(spec2d(2:nx/2+1,nz/2+1)))
!
!	diff(1)=(accum(2)-accum(1))/(kap(2)-kap(1))
!	diff(nz/2+1)=(accum(nz/2+1)-accum(nz/2))/(kap(nz/2+1)-kap(nz/2))
!	Do k=2,nz/2
!	  diff(k)=(accum(k+1)-accum(k-1))/(kap(k+1)-kap(k-1))
!	End Do
!
!	another way to get spectra density funciton
	diff(1:nx/2+1)=euu(1:nx/2+1)/dkap
	Deallocate(euu)
!
	End Subroutine specyz
!*******************************************************************************
	Subroutine specxz(spec2d,kapx,kapz,accum,diff)
	!plot 2D spectra in yz plane
	Real(mytype), Dimension(:,:), Intent(In) :: spec2d
	Real(mytype), Dimension(:,:), Intent(InOut) :: accum,diff
	Real(mytype), Dimension(:,:), Allocatable :: euu
	Real(mytype), Dimension(:), Intent(In) :: kapx, kapz
        Real(mytype) :: specavg,dkapx,dkapz
	Integer :: i,k,i2,k2

	Allocate(euu(1:nx/2+1,1:nz/2+1))
	accum(:,:)=0.0e0; euu(:,:)=0.0e0
!
	accum(1,1)=0.0e0; accum(1,nz/2+1)=spec2d(1,nz/2+1)
	accum(nx/2+1,1)=spec2d(nx/2+1,1); accum(nx/2+1,nz/2+1)=spec2d(nx/2+1,nz/2+1)
!
	Do i=2,nx/2 !k=1 and k=nz/2+1
          accum(i,1)=spec2d(i,1)*2.0e0
          accum(i,nz/2+1)=spec2d(i,nz/2+1)*2.0e0
	End Do
	Do k=2,nz/2 !i=1 and i=nx/2+1
          accum(1,k)=spec2d(1,k)*2.0e0
          accum(nx/2+1,k)=spec2d(nx/2+1,k)*2.0e0
	End Do
	Do k=2,nz/2; Do i=2,nx/2
	  accum(i,k)=(spec2d(i,k)+spec2d(i,nz+2-k))*2.0e0
	End Do; End Do
!
	diff(1:nx/2+1,1:nz/2+1)=accum(1:nx/2+1,1:nz/2+1)
	euu(1:nx/2+1,1:nz/2+1)=accum(1:nx/2+1,1:nz/2+1)
!
	Write(*,*) '~~~~~~~~total energy~~~~~~~~~~', Sum(diff(:,:)), Sum(diff(:,:))/utau_avg_sa**2
	Do i=2,nx/2+1; accum(i,1)=accum(i-1,1)+diff(i,1); End Do
	Do k=2,nz/2+1; accum(1,k)=accum(1,k-1)+diff(1,k); End Do
	Do k=2,nz/2+1; Do i=2,nx/2+1
	  accum(i,k)=accum(i-1,k-1)+Sum(euu(1:i-1,k))+Sum(euu(i,1:k-1))+euu(i,k)
	End Do; End Do
!
	dkapx=kapx(2)-kapx(1); dkapz=kapz(2)-kapz(1)
!	diff(:,:)=0.0e0
!	diff(1,1)=0.5e0*((accum(1,2)-accum(1,1))/dkapz+&
!		(accum(2,1)-accum(1,1))/dkapx)
!	diff(1,nz/2+1)=0.5e0*((accum(1,nz/2+1)-accum(1,nz/2))/dkapz+&
!		(accum(2,nz/2+1)-accum(1,nz/2+1))/dkapx)
!	diff(nx/2+1,1)=0.5e0*((accum(nx/2+1,2)-accum(nx/2+1,1))/dkapz+&
!		(accum(nx/2+1,1)-accum(nx/2,1))/dkapx)
!	diff(nx/2+1,nz/2+1)=0.5e0*((accum(nx/2+1,nz/2+1)-accum(nx/2+1,nz/2))/dkapz+&
!		(accum(nx/2+1,nz/2+1)-accum(nx/2,nz/2+1))/dkapx)
!
!	Do i=2,nx/2
!	  diff(i,1)=(accum(i+1,2)-accum(i-1,2)-accum(i+1,1)+accum(i-1,1))/dkapx/dkapz/2.0e0
!	  diff(i,nz/2+1)=(accum(i+1,nz/2+1)-accum(i-1,nz/2+1)-accum(i+1,nz/2)+accum(i-1,nz/2))/dkapx/dkapz/2.0e0
!	End Do
!	Do k=2,nz/2
!	  diff(1,k)=(accum(2,k+1)-accum(1,k+1)-accum(2,k-1)+accum(1,k-1))/dkapx/dkapz/2.0e0
!	  diff(nx/2+1,k)=(accum(nx/2+1,k+1)-accum(nx/2,k+1)-accum(nx/2+1,k-1)+accum(nx/2,k-1))/dkapx/dkapz/2.0e0
!	End Do
!
!	Do k=2,nz/2; Do i=2,nx/2
!	  diff(i,k)=(accum(i+1,k+1)-accum(i-1,k+1)-accum(i+1,k-1)+accum(i-1,k-1)) &
!	  	/dkapx/dkapz/4.0e0
!	End Do; End Do
!
!	another way to get spectra density funciton
	diff(1:nx/2+1,1:nz/2+1)=euu(1:nx/2+1,1:nz/2+1)/dkapx/dkapz
	Deallocate(euu)
!
	End Subroutine specxz

!*******************************************************************************
	SUBROUTINE write_vtk(uu)
! This subroutine writes out the dataset into vtk file for visualisation
	Real(mytype), Dimension(1:nx,1:nz), Intent(in) :: uu
	Integer :: i,j,k
!
	Open(637, File='super_visu.vtk', status='UNKNOWN')
	Write(637,'(A)') '# vtk DataFile Version 2.0'
	Write(637,'(A)') 'Volume example'
	Write(637,'(A)') 'ASCII'
	Write(637,'(A)') 'DATASET RECTILINEAR_GRID'
	Write(637,'(A,3I10)') 'DIMENSIONS',nx,1,nz
	Write(637,'(A,I10,A)') 'X_COORDINATES', nx, ' float'
	Do i=1,nx; Write(637,'(E15.7)') xc(i); End Do
	Write(637,'(A,I10,A)') 'Y_COORDINATES', 1, ' float'
	Do j=1,1; Write(637,'(E15.7)') yc(j); End Do
	Write(637,'(A,I10,A)') 'Z_COORDINATES', nz, ' float'
	Do k=1,nz; Write(637,'(E15.7)') zc(k); End Do
	Write(637,'(A,I10)') 'POINT_DATA',nx*nz
	Write(637,'(A)') 'SCALARS u float 1'
	Write(637,'(A)') 'LOOKUP_TABLE default'
	DO k=1,nz; DO i=1,nx
  	  Write(637,'(E15.7)') uu(i,k)
	End Do; End Do
	Close(637)
!
	END SUBROUTINE write_vtk
!
!*******************************************************************************
	SUBROUTINE signal_shuffle(uu)
! This subroutine shuffles the signal in wave space (Schlatter & Orlu 2010; Mathis et al 2011)
	Include "fftw3.f"
	Real(mytype), Dimension(-1:nx+2,-1:nz+2), Intent(in) :: uu
	Real(mytype), Dimension(:,:), Allocatable :: uxz,spec2d
	Real(mytype), Dimension(:), Allocatable :: ux,uz,pdf1d,pdf1dx, &
	  	accum,diff,lamx,kapx,arr1d
	Complex(mytype), Allocatable, Dimension(:) :: outfftx,outfftx2
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxz,outfftxz2
	Real(mytype) :: sigma,theta,tmp
	Integer :: i,j,k,ip,kp,np,ind
        Integer*8 :: planx,planxz
  	INTEGER :: i_seed
  	INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
  	INTEGER, DIMENSION(1:8) :: dt_seed

!
  	! ----- Set up random seed portably -----
  	CALL RANDOM_SEED(size=i_seed)
  	ALLOCATE(a_seed(1:i_seed))
  	CALL RANDOM_SEED(get=a_seed)
  	CALL DATE_AND_TIME(values=dt_seed)
  	a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
  	CALL RANDOM_SEED(put=a_seed)
  	DEALLOCATE(a_seed)
!
	Allocate(ux(1:nx),uz(1:nz),uxz(1:nx,1:nz),outfftx(nx/2+1),outfftx2(nx/2+1), &
	outfftxz(1:nx/2+1,1:nz),outfftxz2(1:nx/2+1,1:nz),spec2d(1:nx/2+1,1:nz),&
	arr1d(1:nx))
	Allocate(accum(1:nx/2+1),diff(1:nx/2+1),lamx(1:nx/2+1),kapx(1:nx/2+1))
	accum(:)=0.0e0; diff(:)=0.0e0
	lamx(1)=0.0e0; kapx(1)=0.0e0
	Do i=2,nx/2+1; lamx(i)=length/real(i-1,mytype); kapx(i)=2.0e0*pi/lamx(i); End Do	
!
	ip=nx/2; kp=nz/2; np=100
	Allocate(pdf1d(1:np),pdf1dx(1:np))
!
	Do k=1,nz
	  uxz(1:nx,k)=uu(1:nx,k)  
	End Do
!
!!	shuffle the signal (Schlatter & Orlu 2010)
!	Do k=1,nz
!	  arr1d(1:nx)=uxz(1:nx,k)
!	  Do i=nx,1,-1
!	    CALL RANDOM_NUMBER(theta)
!	    ind=Max(Int(theta*i),1)
!	    tmp=arr1d(i); arr1d(i)=arr1d(ind); arr1d(ind)=tmp
!	  End Do
!	  uxz(1:nx,k)=arr1d(1:nx)
!	End Do
!
!	FFT
	Call dfftw_plan_dft_r2c_1d(planx,nx,ux,outfftx,FFTW_ESTIMATE)
	Call dfftw_execute(planx)
	Call dfftw_destroy_plan(planx)	
	outfftx(:)=outfftx(:)/Real(nx,mytype)
	outfftx(1)=0.0e0
	outfftx2(:)=outfftx(:)
!	inverse FFT
	Call dfftw_plan_dft_c2r_1d(planx,nx,outfftx2,ux,FFTW_ESTIMATE)
	Call dfftw_execute(planx)
	Call dfftw_destroy_plan(planx)
!
!	FFT
	Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,uxz,outfftxz,FFTW_ESTIMATE)
	Call dfftw_execute(planxz)
	Call dfftw_destroy_plan(planxz)	
	outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)	
	outfftxz(1,1)=0.0e0
!!	scramble the signal (Mathis et al 2009)
!	Do k=1,nz; Do i=1,nx/2+1
!	  CALL RANDOM_NUMBER(theta)
!	  theta=theta*2.0e0*pi
!	  outfftxz(i,k)=outfftxz(i,k)*Cmplx(Cos(theta),Sin(theta))
!	End Do; End Do
!
!	calculate spectra
	Open(unit=639,file='signal_lx.plt') 
	Write(639,'(10A)') ' variables="<greek>l</greek><sub>x</sub><sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek>dk", "<greek>F</greek>","k<sub>x</sub><greek>F</greek>"'
!
!	calculate sigma
	spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz(1:nx/2+1,1:nz))*outfftxz(1:nx/2+1,1:nz)
	sigma=Sqrt(Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
	  	Sum(spec2d(nx/2+1,1:nz)))
	Call specxy(spec2d(:,:),kapx,accum,diff)
	Do i=2,nx/2+1
	  Write(639,'(7E15.7)') lamx(i)*retau_avg_sa, &
	    accum(i), diff(i), diff(i)*kapx(i)
	End Do
	Close(639)
!
!	inverse FFT
	Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz,uxz,FFTW_ESTIMATE)
	Call dfftw_execute(planxz)
	Call dfftw_destroy_plan(planxz)
	Open(unit=638,file='ux-ins.plt')
	Write(638,'(10A)') ' variables="x","u"'
	Do i=1,nx
	  Write(638,'(2F20.10)') x(i),uxz(i,kp)
	End Do
	Close(638)
!
!	calculate pdf
	uxz(:,:)=uxz(:,:)/sigma
	Call pdf_1d_super(uxz,uxz,0,np,pdf1dx(1:np),pdf1d(1:np))
!
	Deallocate(ux,uz,outfftx,outfftx2,pdf1d,pdf1dx,uxz,outfftxz,outfftxz2, &
	spec2d,accum,diff,lamx,kapx,arr1d)	
!
	END SUBROUTINE signal_shuffle
!
!***************************************************************************
	SUBROUTINE pdf_1d_super(u,usuper,irangefix,np,pdfx,pdf)
!	This subroutine calculates the 1D pdf (probability density function)
! This subroutine comes from pdf_1d_super in calc_post3d.f90
	Real(mytype), Dimension(1:nx,1:nz), Intent(in) :: u,usuper
	Integer, Intent(in) :: np,irangefix
	Real(mytype), Dimension(1:np), Intent(Out)  :: pdf,pdfx
!
	Integer :: i,j,k, i1,j1,k1, jp, cnt,ipos,cntp,cntn
	Real(mytype), Dimension(:,:), Allocatable :: dumu
	Real(mytype), Allocatable, Dimension(:) :: cdf
	Real(mytype) :: umin,umax,du,ulim,tmp
!
	ipos=2 !condition with positive (1) or negative (0) super scale structures or all structures (>=2)
	cntp=0;cntn=0
	Allocate(dumu(1:nx,1:nz))
	DO k=1,nz; DO i=1,nx
	  dumu(i,k)=u(i,k)
	END DO; END DO
	Allocate(cdf(1:np))
	umin=dumu(1,1); umax=dumu(1,1)
	DO k=1,nz; DO i=1,nx
	  IF(umin.GT.dumu(i,k)) umin=dumu(i,k)
	  IF(umax.LT.dumu(i,k)) umax=dumu(i,k)
	END DO; END DO
!fixed u range
	If(irangefix==1) Then; umin=-5.0e0; umax=5.0e0; End If
	Write(*,*) 'umin=',umin,'umax=',umax
	DO k=1,nz; DO i=1,nx
	  If(ipos==1) Then
	    If(usuper(i,k)<0.0e0) dumu(i,k)=umin-1.0e0
	  Else If(ipos==0) Then
	    If(usuper(i,k)>0.0e0) dumu(i,k)=umax+1.0e0
	  End If
	  If(usuper(i,k)>0.0e0) cntp=cntp+1
	END DO; END DO
	cntn=nx*nz-cntp
!
	du=(umax-umin)/REAL(np-1)
	Do i=1,np
	  pdfx(i)=du*REAL(i-1)+umin
	End Do
	cdf=0.0e0; pdf=0.0e0
	DO i1=1,np
	  cnt=0; 
	  ulim=REAL(i1-1)*du+umin
	  DO k=1,nz; DO i=1,nx
	    IF(dumu(i,k).LT.ulim.And.dumu(i,k).Ge.umin) cnt=cnt+1
	  END DO; End Do
	  If(ipos==1) Then
	    cdf(i1)=REAL(cnt)/REAL(cntp)
	  Else If(ipos==0) Then
	    cdf(i1)=REAL(cnt)/REAL(cntn)
	  Else
	    cdf(i1)=REAL(cnt)/REAL(nx*nz)
	  End If
	END DO
!Write(*,*) cntp,cntn
! get the pdf from cdf
	DO i1=2,np-1
	  pdf(i1)=(cdf(i1+1)-cdf(i1-1))/2.0e0/du
	END DO
! write data to file
	OPEN(UNIT=311,FILE='signal_pdf1d.plt')
	WRITE(311,'(320A)') ' variables="u","pdf","cdf"'
 	WRITE(311,*) 'zone i=',np-1,',f=point'
	DO i=2,np
	  WRITE(311,'(4F20.10)') pdfx(i), pdf(i), cdf(i)
	END DO
	CLOSE(311)
!! 	calculate the mean from pdf
!	tmp=Sum(u(1:nx,1:nz))/Real(nx*nz,mytype)
!	Write(*,*) 'mean=',tmp
!	tmp=0.0e0
!	Do i=1,np
!	  tmp=tmp+pdf(i)*(REAL(i-1)*du+umin)*du
!	End Do
!	Write(*,*) 'mean_pdf=',tmp
!
	DEALLOCATE(dumu,cdf)
	END SUBROUTINE pdf_1d_super
!
!*******************************************************************************
	SUBROUTINE HHT(uu)
! this subroutine calculates Hilbert-Huang transform (Huang et al 1998)
	Include "fftw3.f"
	Real(mytype), Dimension(-1:nx+2,-1:nz+2), Intent(in) :: uu
	Complex(mytype), Allocatable, Dimension(:) :: outfftx
	Real(mytype), Allocatable, Dimension(:) :: ux,uex,umx,udx,utx, &
	  ubx,uhhx,uhx,amp,pha,dpha,kapx
	Real(mytype), Allocatable, Dimension(:) :: usp,xx,aa,bb,cc,dd
	Integer, Allocatable, Dimension(:) :: indt,indb
	Real(mytype), Allocatable, Dimension(:,:) :: amp2d
	Character :: chdum*8,file3d*60
	Real(mytype) :: dum,yam,theta,alpha,beta,uol,ustar,corr1,corr2,tmp,SD
	Integer :: intdum,nsnaps
	Integer :: m,i,j,k,jp,ip,kp,it,nbs,nts,imode,nmode,nrmin
        Integer*8 :: planx
!	total number of HHT decomposed modes
	nmode=9
!	predictive model (Marusic_etal_Science2010)
	alpha=1.0e0; beta=0.5e0
!
	Allocate(outfftx(nx/2+1),ux(1:nx), &
	uex(1:nx),umx(1:nx),utx(1:nx),ubx(1:nx),uhhx(1:nx),uhx(1:nx), &
	amp(1:nx),pha(1:nx),dpha(1:nx))
	Allocate(udx(-1:nx+2),indt(1:nx),indb(1:nx))
!
	nrmin=nx/5
	Allocate(amp2d(1:nx,1:nrmin),kapx(1:nrmin))
	amp2d(:,:)=0.0e0; kapx(:)=0.0e0
	Do i=2,nrmin
	  kapx(i)=2.0e0*pi/(length/Real(i-1,mytype))
	End Do
!
	Do j=jp,jp
	  Do k=nz/2,nz/2
	  kp=k
	  ux(1:nx)=uu(1:nx,kp)
!!	  an artificial signal
	  Do i=1,nx
	    ustar=Sin(2.0e0*pi/(length/20.0e0)*x(i))
	    uol=Sin(2.0e0*pi/(length/5.e0)*x(i))
!	    ux(i)=ustar*(1.0e0+beta*uol)  +alpha*uol
!	    ux(i)=Sin(2.0e0*pi/(length/60.0e0)*x(i))*Sin(2.0e0*pi/(length/30.0e0)*x(i))
	  End Do
!Do i=1,nx/2
!  ux(i)=Sin(2.0e0*pi/(length/60.0e0)*x(i))
!End Do
!Do i=nx/2+1,nx
!  ux(i)=Sin(2.0e0*pi/(length/30.0e0)*x(i))
!End Do
!
!	  Fourier Transform
	  Call dfftw_plan_dft_r2c_1d(planx,nx,ux,outfftx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)
	  Call dfftw_destroy_plan(planx)
	  outfftx(:)=outfftx(:)/Real(nx,mytype)
	  Open(unit=642,file='HHT_spectrum_marginal.plt')
	  Write(642,'(10A)') ' variables="<greek>k</greek><sub>x</sub>",&
&"<greek>F</greek>","k<sub>x</sub><greek>F</greek>"'
	  Write(642,'(A)') 'zone t="FFT"'
	  Do i=2,nrmin
	    tmp=Conjg(outfftx(i))*outfftx(i)*2.0e0
	    Write(642,'(6E20.10)') kapx(i), &
	    tmp/(kapx(3)-kapx(2)), kapx(i)*tmp/(kapx(3)-kapx(2))
	  End Do
!	  remove the mean component
	  outfftx(1)=0.0e0
	  Call dfftw_plan_dft_c2r_1d(planx,nx,outfftx,ux,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)
	  Call dfftw_destroy_plan(planx)
!
!	  empirical mode decomposition (EMD)
	  Open(unit=640,file='HHT_modes.plt')
	  Write(640,'(A)') ' variables="x","u","u<sub>e,1</sub>","u<sub>e,2</sub>"'
	  Write(640,'(A,I3,A)') 'zone t="original"'
	  Do i=1,nx; Write(640,'(4F20.10)') x(i),ux(i),ux(i),ux(i); End Do
	  uhhx(:)=0.0e0
!
	  Do imode=1,nmode
	  ux(1:nx)=ux(1:nx)-uhhx(1:nx)
	  Open(unit=634,file='HHT_ux.plt')
	  Write(634,'(A)') ' variables="x","u","u<sub>l</sub>","u<sub>s</sub>",&
	  &"u<sub>e</sub>","u<sub>e,1</sub>","u<sub>e,2</sub>"'
	  umx(1:nx)=0.0e0; it=0; SD=1.0e0
	  uhhx(1:nx)=ux(1:nx)
	  Do While(SD>1.e-5)
	    it=it+1
	    uhhx(1:nx)=uhhx(1:nx)-umx(1:nx)
	    udx(1:nx)=uhhx(1:nx)
	    udx(-1:0)=uhhx(nx-1:nx); udx(nx+1:nx+2)=uhhx(1:2)
	    utx(:)=0.0e0; ubx(:)=0.0e0
	    ! local maximum and local minimum
	    nbs=0; nts=0
	    Do i=1,nx
	      If(udx(i-1)<udx(i).And.udx(i+1)<udx(i)) Then
		utx(i)=uhhx(i); nts=nts+1; indt(nts)=i
	      End If
	      If(udx(i-1)>udx(i).And.udx(i+1)>udx(i)) Then
		ubx(i)=uhhx(i); nbs=nbs+1; indb(nbs)=i
	      End If
	    End Do
!
	    Allocate(usp(1:nts),xx(1:nts),bb(1:nts),cc(1:nts),dd(1:nts))
!	    interpolation upper envelop
	    bb(:)=0.0e0; cc(:)=0.0e0; dd(:)=0.0e0
	    usp(1:nts)=uhhx(indt(1:nts)); xx(1:nts)=x(indt(1:nts))
	    Call spline(xx, usp, bb, cc, dd, nts) 
	    Do i=1,nx
	      utx(i) = ispline(x(i), xx, usp, bb, cc, dd, nts)
	    End Do
	    Deallocate(usp,xx,bb,cc,dd)
!	    interpolation lower envelop
	    Allocate(usp(1:nbs),xx(1:nbs),bb(1:nbs),cc(1:nbs),dd(1:nbs))
	    bb(:)=0.0e0; cc(:)=0.0e0; dd(:)=0.0e0
	    usp(1:nbs)=uhhx(indb(1:nbs)); xx(1:nbs)=x(indb(1:nbs))
	    Call spline(xx, usp, bb, cc, dd, nbs) 
	    Do i=1,nx
	      ubx(i) = ispline(x(i), xx, usp, bb, cc, dd, nbs)
	    End Do
	    Deallocate(usp,xx,bb,cc,dd)
!	    mean envelop
	    umx(:)=(utx(:)+ubx(:))/2.0e0
!	    write intermedia result
	    Write(634,'(A,I3,A)') 'zone t="',it,'"'
	    Do i=1,nx
	      Write(634,'(7F20.10)') x(i),uhhx(i),umx(i),umx(i),umx(i),utx(i),ubx(i)
	    End Do
!	    check stop criteria
	    SD=0.0e0; tmp=0.0e0
	    Do i=1,nx
	      tmp=tmp+uhhx(i)**2; SD=SD+umx(i)**2
	    End Do
	    SD=SD/tmp
	  End Do
	  Write(*,*) 'sifting time:', it
	  Close(634)
!
	  Write(640,'(A,I3,A)') 'zone t=" mode',imode,'"'
	  Do i=1,nx; Write(640,'(4F20.10)') x(i),uhhx(i),utx(i),ubx(i); End Do
!	  Hilbert spectrum
	  If(1) Then
	    uhx(1:nx)=uhhx(1:nx)
	    Call dfftw_plan_dft_r2c_1d(planx,nx,uhx,outfftx,FFTW_ESTIMATE)
	    Call dfftw_execute(planx)
	    Call dfftw_destroy_plan(planx)
	    outfftx(:)=outfftx(:)/Real(nx,mytype)
	    theta=-pi/2.0e0
	    Do i=2,nx/2+1
	      outfftx(i)=outfftx(i)*Cmplx(Cos(theta),Sin(theta))
	    End Do
	    Call dfftw_plan_dft_c2r_1d(planx,nx,outfftx,uex,FFTW_ESTIMATE)
	    Call dfftw_execute(planx)  
	    Call dfftw_destroy_plan(planx)
!	    calculate the amplitude and phase of the analytic signal
	    Do i=1,nx
	      amp(i)=Sqrt(uhx(i)**2+uex(i)**2)
	      pha(i)=Atan2(uex(i),uhx(i))
	    End Do
	    Do i=2,nx-1
	      dpha(i)=(pha(i+1)-pha(i-1))/2.0e0/dx(1)
	    End Do
	    dpha(1)=(pha(2)-pha(nx))/2.0e0/dx(1)
	    dpha(nx)=(pha(1)-pha(nx-1))/2.0e0/dx(1)
	    Do i=1,nx; Do m=2,nrmin
	      If(dpha(i).Ge.kapx(m-1).And.dpha(i).Le.kapx(m)) Then
		amp2d(i,m)=amp2d(i,m)+amp(i)**2; Cycle
	      End If
	    End Do; End Do
	  End If ! Hilbert spectum
	  End Do ! imode
!
	  Write(640,'(A,I3,A)') 'zone t="residual"'
	  Do i=1,nx; Write(640,'(4F20.10)') x(i),ux(i)-uhhx(i),0.0e0,0.0e0; End Do
	  Close(640)
!
	  Open(unit=641,file='HHT_spectrum.plt')
	  Write(641,'(10A)') ' variables="x","<greek>k</greek><sub>x</sub>",&
&"<greek>F</greek>"'
 	  Write(641,*) 'zone i=',nx,', j=',nrmin-1,',f=point'
	  Do m=2,nrmin; Do i=1,nx
	    Write(641,'(6E20.10)') x(i),kapx(m),amp2d(i,m)/Real(nx,mytype)/(kapx(3)-kapx(2))
	  End Do; End Do
	  Close(641)
!	  HHT marginal spectrum
	  Write(642,'(A)') 'zone t="HHT"'
	  Do m=2,nrmin
	    Write(642,'(6E20.10)') kapx(m), &
	    Sum(amp2d(1:nx,m))/Real(nx,mytype)/(kapx(3)-kapx(2)), &
	    kapx(m)*Sum(amp2d(1:nx,m))/Real(nx,mytype)/(kapx(3)-kapx(2))
	  End Do
	  Close(642)
!
	  End Do ! k
	End Do ! j
!
!
	Deallocate(ux,uex,udx,uhx,outfftx,utx,ubx,indt, &
	  indb,uhhx,umx,amp,pha,dpha,amp2d,kapx)	
!
	END SUBROUTINE HHT
!
!*******************************************************************************
   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
! Downloaded from http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90 
! modified for periodic boundary condition
implicit none
integer n
double precision x(n), y(n), b(n), c(n), d(n), h(n), g(n), r(n)
double precision bb(n),cc(n),dd(n),gg(n)
integer i, j, gap
double precision alpha
! solving matrix
!=================================================================================
![    b(1)	d(1)	   0	  0	   0	...	g(1)  ] [ x(1) ]  [ c(1) ]
!|    g(2)      b(2)      d(2)	  0	   0	...	0     | | x(2) |  | c(2) |
!|     0	g(3)      b(3)   d(3)      0	...	0     |x|  ... |= | c(3) |
!|    ...		...	 ...	  ...	...	...   | |  ... |  | ...  |
!|     0	 0	   0     ...     g(n-1) b(n-1)  d(n-1)| |x(n-1)|  |c(n-1)|
![    d(n)	 0         0	   0      ...    g(n)   b(n)  ] [ X(n) ]  [ c(n) ]
!=================================================================================
!                                        | |
!					 V V
!======================================================================================================================================
![2*(h(n)+h(1))	h(1)	   0	  0	   0	...		h(n)	      ] [ x(1) ]  [     (y(2)-y(1))/h(1)-(y(1)-y(n))/h(n)     ]
!|    h(1)    2*(h(1)+h(2))  h(2)	  0	   0	...		0     | | x(2) |  |     (y(3)-y(2))/h(2)-(y(2)-y(1))/h(1)     |
!|     0		h(2)  2*(h(2)+h(3))  h(3)  0	...		0     |x|  ... |= |     (y(4)-y(3))/h(3)-(y(3)-y(2))/h(2)     |
!|    ...		...	...	...	...	...		...   | |  ... |  |                ...			      |
!|     0	         0	   0    ...   h(n-2) 2*(h(n-2)+h(n-1))  h(n-1)| |x(n-1)|  |(y(n)-y(n-1))/h(n-1)-(y(n-1)-y(n-2))/h(n-2)|
![    h(n)	 0         0	   0    ...       h(n-1)   2*(h(n-1)+h(n))    ] [ X(n) ]  [   (y(1)-y(n))/h(n)-(y(n)-y(n-1))/h(n-1)   ]
!======================================================================================================================================
gap = n-1
! check input
if ( n < 2 ) then
  write(*,'(A)') 'STOP: total number of local minima (or maxima) smaller than 2! '
  Stop
end if
if ( n < 3 ) then
  write(*,'(A)') ' total number of local minima (or maxima) smaller than 3!'
  return
end if
!
! step 1: preparation
!
h(:)=0.0e0; g(:)=0.0e0; d(:)=0.0e0; r(:)=0.0e0
h(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/h(1)
do i = 2, n-1
  h(i) = x(i+1) - x(i)
  g(i) = h(i-1)
  d(i) = h(i)
  b(i) = 2.0*(h(i-1) + h(i))
  c(i+1) = (y(i+1) - y(i))/h(i)
  c(i) = c(i+1) - c(i)
end do
h(n) = length-x(n)+x(1)
!
! step 2: end conditions 
! left
b(1) = 2.0*(h(n) + h(1))
d(1) = h(1)
g(1) = h(n)
c(1) = (y(2) - y(1))/h(1) - (y(1) - y(n))/h(n)
! right
b(n) = 2.0*(h(n-1) + h(n))
d(n) = h(n)
g(n) = h(n-1)
c(n) = (y(1) - y(n))/h(n) - (y(n) - y(n-1))/h(n-1)
!
r(1) = g(1); r(n-1)=d(n-1); r(n)=b(n)
!
!
! step 3: forward elimination 
!
do i = 2, n-1
  alpha = g(i)/b(i-1)
  b(i) = b(i) - alpha*d(i-1)
  g(i) = g(i) - alpha*b(i-1)
  c(i) = c(i) - alpha*c(i-1)
  r(i) = r(i) - alpha*r(i-1)
end do
!
! step 4: back substitution
!
do i = n-2,1,-1
  alpha = d(i)/b(i+1)
  d(i) = d(i)-alpha*b(i+1)
  r(i) = r(i)-alpha*r(i+1)
  c(i) = c(i)-alpha*c(i+1)
end do
!
!
!! eliminate the (n,1) element
alpha=d(n)/b(1)
d(n)=d(n)-alpha*b(1); b(n)=r(n)-alpha*r(1)
r(n)=r(n)-alpha*r(1); c(n)=c(n)-alpha*c(1)
!! eliminate the (n,n-1) element
alpha=g(n)/b(n-1)
g(n)=g(n)-alpha*b(n-1); b(n)=r(n)-alpha*r(n-1)
r(n)=r(n)-alpha*r(n-1); c(n)=c(n)-alpha*c(n-1)
!
! step 5: calculate coefficient ci
c(n)=c(n)/r(n)
do i=n-1,1,-1
  c(i)=(c(i)-r(i)*c(n))/b(i)
end do
!
! step 5: compute spline coefficients
!
do i = 1, n-1
  b(i) = (y(i+1) - y(i))/h(i) - h(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/h(i)
  c(i) = 3.0e0*c(i)
end do
! last segment
  c(n) = 3.0e0*c(n)
  b(n) = b(n-1)+2.0e0*c(n-1)*h(n-1)+3.0e0*d(n-1)*h(n-1)**2
  d(n) = (b(1)-b(n)-2.0e0*c(n)*h(n))/3.0e0/h(n)**2
!! check gradient
!Write(*,*) b(n),c(n),d(n)
!Write(*,'(A)') 'S 1st derivative'
!Do i=1,n
!Write(*,'(2F10.5)') b(i),b(i)+2.0e0*c(i)*h(i)+3.0e0*d(i)*h(i)**2 
!End Do
!Write(*,'(A)') 'S 2nd derivative'
!Do i=1,n
!Write(*,'(2F10.5)') c(i),c(i)+3.0e0*d(i)*h(i) 
!End Do

end subroutine spline
!
!*******************************************************************************
  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
! Downloaded from http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90 
implicit none
double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
!if(u <= x(1)) then
!  ispline = y(1)
!  return
!end if
!if(u >= x(n)) then
!  ispline = y(n)
!  return
!end if
! periodic boundary
if(u < x(1)) then
  dx = u + length - x(n)
  ispline = y(n) + dx*(b(n) + dx*(c(n) + dx*d(n)))
  return
end if
if(u > x(n)) then
  dx = u - x(n)
  ispline = y(n) + dx*(b(n) + dx*(c(n) + dx*d(n)))
  return
end if

!*
!  binary search for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline
!
!*******************************************************************************
	SUBROUTINE HHT2D(uu)
! this subroutine 2D Hilbert-Huang empirical mode decomposition (Agostini & Leschziner 2014)
! envelop surface fitting using thin plate spline (TPS), based on 
! http://people.sc.fsu.edu/~jburkardt/f_src/rbf_interp_2d/rbf_interp_2d.html
! parallel using openmp
! matrix inversion uses NAG SMP library for acceleration
	Use omp_lib
        Use nag_library, Only: dsytrf, dsytrs, nag_wp
	Real(mytype), Dimension(-1:nx+2,-1:nz+2), Intent(in) :: uu
	Real(mytype), Allocatable, Dimension(:,:) :: ux,umx,uhhx,utx,  &
	  ubx,xz,xzi,a,ainv
	Integer, Allocatable, Dimension(:,:) :: indt,indb
	Real(mytype), Allocatable, Dimension(:) :: ff,coeff,ffi,pm,ap,p,    &
	  work
	Integer, Allocatable, Dimension(:) :: ipiv
	Character :: chdum*8,file3d*60
	Character*1 :: uplo
	Real(mytype) :: dum,tmp,SD
	Real(mytype) :: xmin,xmax,zmin,zmax,volume,r0,edge,xtmp,ztmp,vtmp,rtmp
	Real(mytype) :: pr,pap,rap,beta,alpha
	Integer :: intdum,ostart,oend,ifault,nullty,lwork,lda,ldb,nrhs
	Integer :: m,n,i,j,k,jp,ip,kp,it,nbs,nts,imode,nmode,ic,ns,iupper,iter
	Integer :: nd,ni,nxz
	Integer :: t_num,tid

!	total number of HHT decomposed modes
	nmode=5; edge=0.2e0
!
	Allocate(ux(-1:nx+2,-1:nz+2),umx(-1:nx+2,-1:nz+2), &
	uhhx(-1:nx+2,-1:nz+2),utx(-1:nx+2,-1:nz+2),ubx(-1:nx+2,-1:nz+2))
	Allocate(indt(1:nx*nz,1:2))
!
	Do jp=1,1
	  ux(:,:)=uu(:,:)
!	  remove the mean component
	  ux(:,:)=ux(:,:)-Sum(ux(1:nx,1:nz))/Real(nx*nz,mytype)
!!	  an artificial signal
	  Do i=1,nx; Do k=1,nz
!	    ux(i,k)=Sin(2.0e0*pi/(length/10.0e0)*x(i))*Sin(2.0e0*pi/(width/10.0e0)*zc(k)) &
!	    +2.0e0*Sin(2.0e0*pi/(length/2.0e0)*x(i))*Sin(2.0e0*pi/(width/2.0e0)*zc(k))
	  End Do; End Do
	  Do i=-1,0
!	    ux(i,:)=ux(nx+i,:); ux(nx+1-i,:)=ux(1-i,:)
!	    ux(:,i)=ux(:,nz+i); ux(:,nz+1-i)=ux(:,1-i)
	  End Do
!
!
!	  empirical mode decomposition (EMD)
	  Open(unit=641,file='HHT2D_modes.plt')
	  Write(641,'(A)') ' variables="x","z","u","u<sub>e,1</sub>","u<sub>e,2</sub>"'
 	  Write(641,'(A,I5,A,I5,A)') 'zone t="original", i=',nx,', j=',nz, &
	    ',f=point'
	  Do k=1,nz; Do i=1,nx
	    Write(641,'(5F20.10)') x(i),zc(k),ux(i,k),ux(i,k),ux(i,k)
	  End Do; End Do
	  uhhx(:,:)=0.0e0
!
	  Do imode=1,nmode
	  ux(:,:)=ux(:,:)-uhhx(:,:)
	  umx(:,:)=0.0e0; it=0; SD=1.0e0
	  uhhx(:,:)=ux(:,:)
Do ip=1,10 !While(SD>5.e-2)
	    it=it+1
	    Write(*,'(A,I5,A,I5,A)') '-------------------- imode=',imode, &
	       ',    isifting=',it, '--------------------'
	    uhhx(:,:)=uhhx(:,:)-umx(:,:)
!
	    Do iupper=0,1
!	    upper envelop
	    If(iupper==1) Then; Write(*,'(A)') '===> extract upper envelop ... '
	    Else; Write(*,'(A)') '===> extract lower envelop ... '; End If
!	    local maximum and local minimum
	    nts=0
	    If(iupper==1) Then
	      Do k=1,nz; Do i=1,nx
	        ic=1
	        Do m=-1,1
	        If(ic==1) Then; Do n=-1,1
		  If(m.Ne.0.Or.n.Ne.0) Then
		    If(uhhx(i+m,k+n).Ge.uhhx(i,k)) Then; ic=0; Cycle; End If
		  End If
	        End Do; End If
	        End Do
	        If(ic==1) Then
	          nts=nts+1; indt(nts,1)=i; indt(nts,2)=k
	        End If
	      End Do; End Do
	    Else
	      Do k=1,nz; Do i=1,nx
	        ic=1
	        Do m=-1,1
	        If(ic==1) Then; Do n=-1,1
		  If(m.Ne.0.Or.n.Ne.0) Then
		    If(uhhx(i+m,k+n).Le.uhhx(i,k)) Then; ic=0; Cycle; End If
		  End If
	        End Do; End If
	        End Do
	        If(ic==1) Then
	          nts=nts+1; indt(nts,1)=i; indt(nts,2)=k
	        End If
	      End Do; End Do
	    End If
	    Write(*,*) 'nts=',nts
!	    extract the envelop
!	    minimize the effect of non-periodic
	    xmin=edge; xmax=length-edge
	    zmin=edge; zmax=width-edge
	    ns=0
	    Do i=1,nts
	      xtmp=x(indt(i,1)); ztmp=zc(indt(i,2))
	      If(xtmp<xmin) ns=ns+1
	      If(xtmp>xmax) ns=ns+1
	      If(ztmp<zmin) ns=ns+1
	      If(ztmp>zmax) ns=ns+1
	      If(xtmp<xmin.And.ztmp<zmin) ns=ns+1
	      If(xtmp>xmax.And.ztmp<zmin) ns=ns+1
	      If(xtmp>xmax.And.ztmp>zmax) ns=ns+1
	      If(xtmp<xmin.And.ztmp>zmax) ns=ns+1	      
	    End Do
	    nd=nts+ns; ni=nx*nz
	    Allocate(xz(1:2,1:nd),ff(1:nd),coeff(1:nd))
!
	    ns=nts
	    Do i=1,nts
	      xtmp=x(indt(i,1)); ztmp=zc(indt(i,2))
	      xz(1,i)=xtmp; xz(2,i)=ztmp
	      ff(i)=uhhx(indt(i,1),indt(i,2))
	      If(xtmp<xmin) Then
		ns=ns+1; ff(ns)=ff(i); 
		xz(1,ns)=xz(1,i)+length; xz(2,ns)=xz(2,i)
	      End If
	      If(xtmp>xmax) Then
		ns=ns+1; ff(ns)=ff(i); 
		xz(1,ns)=xz(1,i)-length; xz(2,ns)=xz(2,i)
	      End If
	      If(ztmp<zmin) Then
		ns=ns+1; ff(ns)=ff(i); 
		xz(1,ns)=xz(1,i); xz(2,ns)=xz(2,i)+width
	      End If
	      If(ztmp>zmax) Then
		ns=ns+1; ff(ns)=ff(i); 
		xz(1,ns)=xz(1,i); xz(2,ns)=xz(2,i)-width
	      End If
	      If(xtmp<xmin.And.ztmp<zmin) Then
		ns=ns+1; ff(ns)=ff(i); 
		xz(1,ns)=xz(1,i)+length; xz(2,ns)=xz(2,i)+width
	      End If
	      If(xtmp>xmax.And.ztmp<zmin) Then
		ns=ns+1; ff(ns)=ff(i); 
		xz(1,ns)=xz(1,i)-length; xz(2,ns)=xz(2,i)+width
	      End If
	      If(xtmp>xmax.And.ztmp>zmax) Then
		ns=ns+1; ff(ns)=ff(i); 
		xz(1,ns)=xz(1,i)-length; xz(2,ns)=xz(2,i)-width
	      End If
	      If(xtmp<xmin.And.ztmp>zmax) Then
		ns=ns+1; ff(ns)=ff(i); 
		xz(1,ns)=xz(1,i)+length; xz(2,ns)=xz(2,i)-width
	      End If
	    End Do
!
	    xmax=maxval(xz(1,1:nd)); xmin=minval(xz(1,1:nd))
	    zmax=maxval(xz(2,1:nd)); zmin=minval(xz(2,1:nd))
!	    Write(*,'(A,4F10.5)') 'bounds:', xmin,xmax,zmin,zmax
	    volume=(xmax-xmin)*(zmax-zmin)
	    m=2
	    r0=(volume/Real(nd,mytype))**Real(m,mytype)
!-------------------------------------------------------------------------------
!	    evaluate the coefficient using thin plate spline (TPS)
!	    Call rbf_weight( m, nd, xz, r0, phi3, ff, coeff ) ! rbf_interp_2d
	    Allocate(a(1:nd,1:nd))
	    !$omp parallel default(private),shared(m,nd,xz,xzi,r0,a)
	    t_num = OMP_GET_NUM_THREADS()
	    tid=OMP_GET_THREAD_NUM()
	    If(tid==0) Write(*,'(A,I5)') "Threads NO.=",t_num
	    !$omp do
	    do i = 1, nd
	      do j = 1, i-1
	        rtmp = sqrt ( sum ( ( xz(1:m,i) - xz(1:m,j) )**2 ) )
	        a(i,j) = rtmp**2 * log ( rtmp / r0 )
	      end do
 	      a(i,i)=0.0e0
	    end do
	    !$omp enddo
	    !$omp end parallel
!	    Solve for the weights.
	    Write(*,*) 'starting inverse matrix for coefficents ...'
	    ostart=omp_get_wtime()
!	    call r8mat_solve_svd ( nd, nd, a, ff, coeff )
	    lda=nd; ldb=nd; lwork=64*nd; uplo='L'; nrhs=1
	    Allocate (work(lwork),ipiv(nd))
	    Call dsytrf(uplo,nd,a,lda,ipiv,work,lwork,ifault)
	    Call dsytrs(uplo,nd,nrhs,a,lda,ipiv,ff,ldb,ifault)
	    coeff(:)=ff(:)
	    Deallocate(work,ipiv)
!
	    oend=omp_get_wtime()
	    Write(*,*) 'running time:', oend-ostart
	    Write(*,*) 'rbf_weight done'
	    Deallocate(a)
!
	    Allocate(xzi(1:2,1:nx*nz),ffi(1:nx*nz))
	    Do k=1,nz; Do i=1,nx
	      nxz=i+(k-1)*nx; xzi(1,nxz)=x(i); xzi(2,nxz)=zc(k)
	    End Do; End Do
	    Write(*,*) 'starting interpolating ...'
!-------------------------------------------------------------------------------
!	    evaluate the values on xz plane
!	    Call rbf_interp( m, nd, xz, r0, phi3, coeff, ni, xzi, ffi ) ! rbf_interp_2d
	    ffi(:)=0.0e0
	    !$omp parallel default(private),shared(m,ni,nd,coeff,xz,xzi,r0,ffi)
	    !$omp do
	    do i = 1,ni
	      do j = 1, nd
	        rtmp = sqrt ( sum ( ( xzi(1:m,i) - xz(1:m,j) )**2 ) )
	        if ( rtmp .le. 0.0e0 ) then
	          vtmp = 0.0e0
	        else
	          vtmp = rtmp**2 * log ( rtmp / r0 )
	        end if
	        ffi(i)=ffi(i)+vtmp*coeff(j)
	      end do  
	    end do
	    !$omp enddo
	    !$omp end parallel
	    Write(*,*) 'rbf_interp done'
!-------------------------------------------------------------------------------
!
	    If(iupper==1) Then
	      Do k=1,nz; Do i=1,nx; nxz=i+(k-1)*nx; utx(i,k)=ffi(nxz);
	      End Do; End Do
	    Else
	      Do k=1,nz; Do i=1,nx; nxz=i+(k-1)*nx; ubx(i,k)=ffi(nxz);
	      End Do; End Do
	    End If
!
	    Deallocate(xz,ff,coeff,xzi,ffi)
	    End Do ! iupper
!
	    umx(:,:)=0.5e0*(utx(:,:)+ubx(:,:))
!	    check stop criteria
	    SD=0.0e0; tmp=0.0e0
	    Do k=1,nz; Do i=1,nx
	      tmp=tmp+uhhx(i,k)**2; SD=SD+umx(i,k)**2
	    End Do; End Do
	    SD=SD/tmp
	    Write(*,'(A,F10.5)') 'SD=',SD
!
	  End Do
	  Write(*,*) 'sifting time:', it
!
	  Write(641,'(A,I5,A,I5,A,I5,A)') 'zone t="mode',imode, '", i=',nx, &
	    ', j=',nz,',f=point'
	  Do k=1,nz; Do i=1,nx
	    Write(641,'(5F20.10)') x(i),zc(k),uhhx(i,k)-umx(i,k),utx(i,k),ubx(i,k)
	  End Do; End Do
!
	  End Do ! imode
!
	  Write(641,'(A,I5,A,I5,A)') 'zone t="residual", i=',nx,', j=',nz, &
	    ',f=point'
	  Do k=1,nz; Do i=1,nx
	    Write(641,'(5F20.10)') x(i),zc(k),ux(i,k)-uhhx(i,k),utx(i,k),     &
	      ubx(i,k)
	  End Do; End Do
	  Close(641)
!
	End Do ! j
!
	Deallocate(utx,ubx, &
	  uhhx,umx,indt)	
!
	END SUBROUTINE HHT2D
!
	END MODULE calc_super




























