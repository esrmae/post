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
	Integer :: i,j,k
	Real(mytype) :: dumm
	Real(mytype), Dimension(:,:), Allocatable :: dummarr
!
	simit=1
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
	file2d='OUTPUT_2D_SHEAR.dat'
	Else If(lp_snap_2d == 1) Then
	file2d='OUTPUT_2D_SNAP.dat'
	Else; Goto 1112
	End IF
	Write(*,*) 'READING...', file2d
	CALL change_dir(file2d,n)
	Open(507, File=file2d, Status='OLD', Form='UNFormatTED', access='stream')
	Do l=1,nsmpl_snap_2d
          Read(507) time_snap
	  Write(*,*) time_snap
	  Read(507) uu(:,:) !u
	  Read(507) vv(:,:) !v
	  Read(507) dummarr(:,:) !w
	  Read(507) dummarr(:,:) !p
! calculate averaged spectra
!	  Call spectraxzavg
	  Call spectraxz_snap(l)
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
!
	END MODULE calc_super




























