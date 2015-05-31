	MODULE calc_post3d
!	USE shared_data
	USE shared_post
	USE cal_coefficient
	USE read_post
	IMPLICIT NONE
	CONTAINS
!*******************************************************************************
! This code performs the calculations on the values that have been read in 
! from the read_post module.
!*******************************************************************************
	SUBROUTINE calc3d
        Integer :: nsnaps,tp,n,l,m,cor_cntx,cor_cntz,lam_cnt,iwshear,icorr2d,imodu,il2rms,iam,isk
	Character :: file3d*40
	Real(mytype) :: utau_now
        Real(mytype), dimension(8) :: dum
 	integer, dimension(3) :: intdum
	Real(mytype), Dimension(ny) :: umn_now,wmn_now,l2_mean,l2_rms,l2n_mean
	Real(mytype), Dimension(:,:,:,:), Allocatable :: dummy_3d,vx
	Real(mytype), Dimension(:,:,:), Allocatable :: l2_array, wx,AM
	Real(mytype), Dimension(:,:), Allocatable :: dummy_1d, corrm1,corrm2,corr2d
	Integer, Dimension(:,:,:), Allocatable :: pmin
	Integer :: i,j,k,j1,j2,ipldec
	Integer, Dimension(ny) :: lamn_cnt
!
!	Initialisations
	Allocate(dummy_3d(-1:nx+2,-1:ny+2,-1:nz+2,4),corr(max_xz+1,2,3),spec(max_xz,2,3))
	simit=sim3d
	cor_cntx = 0; cor_cntz = 0; lam_cnt=0; tp=0; lamn_cnt(:)=0
	spec=0.0e0; corr = 0.0e0
	l2_mean=0.0e0; l2_rms=0.0e0; l2n_mean=0.0e0
	iwshear=1
	icorr2d=0; imodu=0; jmodu=ny/2
	il2rms=0; iquad=0
	ipldec=0 ! decompose plasma force region
        iam=0 ! amplitude modulation Mathis et al 2009 (calc_am1p; calc_am1p_z; calc_am2p; calc_am2p_simplify; calc_am2p_simplify_z)
	isk=1 ! skewness decomposition (Mathis et al 2011)
	If (jvis == 0) jvis=ny
        If (ilamq == 1) ivort3d=1
	If (imodu == 1) icorr2d=1
!
	If (ilamq == 1) Allocate(l2_array(nxt,nyt,nzt),pmin(nxt,nyt,nzt),vx(nxt,nyt,nzt,3))
	If(ivort3d == 1) Allocate(wx(nx,ny,nz))
	If(icorr2d==1) Allocate(corrm1(-1:nx+2,-1:nz+2),corrm2(-1:nx+2,-1:nz+2), &
		corr2d(-1:nx+2,-1:nz+2))
	If(imodu==1) Then; Allocate(AM(1:jmodu,1:jmodu,1:3)); End If
	Call open_files3d
!
!	DO n = nstart(simit), nruns(simit)
!
	  If (restrt == 1) then
	    Open(13, File=restrt_file, Status='OLD', Form='UNFormatTED',access='stream')
	    nsnaps=1; tvis=1
	  Else 
	    file3d='OUTPUT_3D_SNAP.dat'
	    l = n-nstart(simit)+1
	    CALL change_dir(file3d,n)
	    Open(13, File=file3d, Status='OLD', Form='UNFormatTED')
	    CALL header_skip(13)
	    nsnaps=snap_3d(l,simit)
	  End If
!
	  DO m=1, nsnaps
	    tp=tp+1
!!	    Read(13) intdum(1),intdum(2),intdum(3),dum(1),dum(2),dum(3),dum(4),dum(5),dum(6),dum(7),dum(8)
!	    Read(13) intdum(:),dum(:)
!	    write(*,*) intdum, dum
!	    print *,'reading u from', restrt_file
!	    Read(13) dummy_3d(:,:,:,1) !u
!           print *,'reading v from', restrt_file
!	    Read(13) dummy_3d(:,:,:,2) !v
!           print *,'reading w from', restrt_file
!	    Read(13) dummy_3d(:,:,:,3) !w
!           print *,'reading p from', restrt_file
!	    Read(13) dummy_3d(:,:,:,4) !p	
!!
	    If(iplot3d.Eq.1) then
!	      tecplot file format
	      If(itec == 1) Then
	        Write(303,'(10f20.10)') (((xc(i),i=1,nx),j=1,jvis),k=1,nz)
	        Write(303,'(10f20.10)') (((yc(j),i=1,nx),j=1,jvis),k=1,nz)
	        Write(303,'(10f20.10)') (((zc(k),i=1,nx),j=1,jvis),k=1,nz)
	        Write(303,'(10f20.10)') (((0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1)),i=1,nx),j=1,jvis),k=1,nz)
	        Write(303,'(10f20.10)') (((0.5*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2)),i=1,nx),j=1,jvis),k=1,nz)
	        Write(303,'(10f20.10)') (((0.5*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3)),i=1,nx),j=1,jvis),k=1,nz)
	        Write(303,'(10f20.10)') (((dummy_3d(i,j,k,4),i=1,nx),j=1,jvis),k=1,nz)
	      Else
!	        paraview format
	        Write(303,'(A,3I10)') 'DIMENSIONS',nx,jvis,nz
		Write(303,'(A,I10,A)') 'X_COORDINATES', nx, ' float'
		Do i=1,nx; Write(303,'(E15.7)') xc(i); End Do
		Write(303,'(A,I10,A)') 'Y_COORDINATES', jvis, ' float'
		Do j=1,jvis; Write(303,'(E15.7)') yc(j); End Do
		Write(303,'(A,I10,A)') 'Z_COORDINATES', nz, ' float'
		Do k=1,nz; Write(303,'(E15.7)') zc(k); End Do
		Write(303,'(A,I10)') 'POINT_DATA',nx*jvis*nz
!		Write(303,'(A)') 'SCALARS u float 1'
!		Write(303,'(A)') 'LOOKUP_TABLE default'
!		DO k=1,nz; DO j=1,jvis; DO i=1,nx
!  		  Write(303,'(10F20.10)') 0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
!		End Do; End Do; End Do
                Write(303,'(A)') 'SCALARS p float 1'
                Write(303,'(A)') 'LOOKUP_TABLE default'
                DO k=1,nz; DO j=1,jvis; DO i=1,nx
                  Write(303,'(10F20.10)') dummy_3d(i,j,k,4)
                End Do; End Do; End Do
                Write(303,'(A)') 'VECTORS U float'
                DO k=1,nz; DO j=1,jvis; DO i=1,nx
                  Write(303,'(10F20.10)') 0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1)),0.5*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2)), &
                                          0.5*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
                End Do; End Do; End Do
	      End If
	      Close(303)
	    Endif
!
	    If ((tp == tvis).or.(tvis == 0)) then

	      Call dat_calc(dummy_3d,utau_now,umn_now,wmn_now)
!
! 	      remove the mean
!             DO k=0,nzt; DO j=1,ny; DO i=0,nxt
!               dummy_3d(i,j,k,1)=dummy_3d(i,j,k,1)-umn_now(j) !Sum(dummy_3d(1:nx,j,1:nz,1))/Real(nx*nz,mytype)
!               dummy_3d(i,j,k,3)=dummy_3d(i,j,k,3)-wmn_now(j) !Sum(dummy_3d(1:nx,j,1:nz,3))/Real(nx*nz,mytype)
!             End Do; End Do; End Do
!	      write wall shear stress
	      If(iwshear==1) Call wallshearplot(dummy_3d)
	      !iplot2d=1 create the 2d slice base on restart file
	      If(iplot2d.Eq.1.Or.lp_snap_2d.Eq.1) Call write_2d_slice(dummy_3d,umn_now,wmn_now)
	      !iplot1d=1 create the 2d slice base on restart file
	      If(iplot1d.Eq.1.Or.lp_snap_1d.Eq.1) Call write_1d_slice(dummy_3d)
	      If (ivort3d==1) Call vort_vis(dummy_3d,wx)
	      If (ilamq==1) then
		Call lambdaQ(dummy_3d,l2_array,pmin,vx,wx)
		Call lambdaQ_mean(l2_array,l2_mean,l2n_mean,l2_rms,lam_cnt,lamn_cnt)
	      End If
!!	      calculate the correlation of uu,vv,ww and their energy spectra
!	      If (icorr==1) Call correlation(dummy_3d,cor_cntx,cor_cntz,umn_now,wmn_now)
!	      Call pdf
!!	      2d correlation for the modulation effect from super scale structures.
!	      If(imodu==1) Then
!	        AM(:,:,:)=0.0e0
!	        Do j1=1,jmodu
!		Write(*,'(A,A,F10.5,A,F10.5)') 'modulation at: ', 'y1=',yc(j1)
!		Do j2=1,jmodu
!	          corrm1(-1:nx+2,-1:nz+2)=dummy_3d(-1:nx+2,j1,-1:nz+2,1)
!	          corrm2(-1:nx+2,-1:nz+2)=dummy_3d(-1:nx+2,j2,-1:nz+2,1)
!		  Call correlation_2D(corrm1,corrm2,corr2d)
!		  AM(j2,j1,1)=Maxval(corr2d(1:nx,1:nz))
!		  AM(j2,j1,2:3)=Maxloc(corr2d(1:nx,1:nz))
!		  Write(*,*) 'Maximum location', AM(j2,j1,2), AM(j2,j1,3), AM(j2,j1,1)
!	        End Do; End Do
!		Call write_modulation(AM)
!		Deallocate(AM)
!	      End If
!!	      2d correlation for the modulation effect from at different range of scales.
	      If(imodu==1) Call correlation_spec()
!!	      conditional average by positive and negative super scale structures
!	      Call super_condition()
!	      Call struct_log
	      If (icorr2d==1) Then
	        corrm1(-1:nx+2,-1:nz+2)=dummy_3d(-1:nx+2,26,-1:nz+2,1)
	        corrm2(-1:nx+2,-1:nz+2)=dummy_3d(-1:nx+2,26,-1:nz+2,2)
	 	Call correlation_2D(corrm1,corrm2,corr2d)
	      End If
!!	      get uv cross correlation for multiple restart files
!	      Call correlation_2D_uv
!!	      spectra calculation
!             Call spectra2D(dummy_3d,umn_now,wmn_now)
!             Call spectra_uv(dummy_3d,umn_now,wmn_now)
!---------------------------------------------------------------------
	    End If
	  END DO
	  Call write_mean(l2_array,l2_mean,l2n_mean,l2_rms,utau_now,lam_cnt,lamn_cnt,cor_cntx,cor_cntz)
	  If(il2rms==1) Call write_l2rms
	  Close(13)
!    	END DO
	If(icorr2d==1) Deallocate(corrm1,corrm2,corr2d)
!	decompose the plasma force region
	If(ipldec==1) Call plasma_decompose
!	calculate the quadrant events
	If(iquad==1) Call calc_quad
!	calculate the one (two) point amplitude modulation (AM) Mathis_etal_JFM2009 (Schlatter & Orlu 2010)
	If(iam==1) Call calc_am2p_simplify
!	calculate the decomposed skewness
	If(isk==1) Call skew_decompose
!
	Deallocate(dummy_3d)
	END SUBROUTINE calc3d

!*******************************************************************************
	SUBROUTINE dat_calc(dummy_3d,utau_now,umn,wmn)
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(in) :: dummy_3d
	Real(mytype), Intent(out) :: utau_now
	Real(mytype), Dimension(nxt) :: utaux
	Real(mytype), Dimension(ny), Intent(out) ::umn,wmn
	Real(mytype), Dimension(:,:,:,:), Allocatable :: rotfield
	Real(mytype) :: theta,urms,vrms,wrms,u2,v2,w2,uv,uw,vw,urot,wrot
!
	Integer :: i,j,k,ig,jmin,jmax,irotate,iw
!
	theta=-45.0e0/180.0e0*pi
        utau_now=0.0e0
        utaux=0.0e0
        DO k=1,nzt; DO i=1,nxt
                utau_now=utau_now+sqrt(abs(nu*(gcdy_f(1,1,1,1)*dummy_3d(i,2,k,1)+gcdy_f(1,2,1,1)*dummy_3d(i,1,k,1))))
                utau_now=utau_now+sqrt(abs(nu*(gcdy_f(nyt,3,1,1)*dummy_3d(i,ny,k,1)+gcdy_f(nyt,4,1,1)*dummy_3d(i,ny-1,k,1))))
                utaux(i)=utaux(i)+sqrt(abs(nu*(gcdy_f(1,1,1,1)*dummy_3d(i,2,k,1)+gcdy_f(1,2,1,1)*dummy_3d(i,1,k,1))))
                utaux(i)=utaux(i)+sqrt(abs(nu*(gcdy_f(nyt,3,1,1)*dummy_3d(i,ny,k,1)+gcdy_f(nyt,4,1,1)*dummy_3d(i,ny-1,k,1))))
        END DO; END DO
        utau_now=utau_now/Real(nxt)/Real(nzt)/Real(2.0e0)
        utaux=utaux/Real(nzt)/Real(2.0e0)
        Write(*,*)'utau_check', utau_avg_sa, utau_now 

	umn=0.0e0
	wmn=0.0e0
	Do k=1,nz; Do j=1,ny; Do i=1,nx
	  umn(j)=umn(j)+dummy_3d(i,j,k,1)
	  wmn(j)=wmn(j)+dummy_3d(i,j,k,3)
	End Do; End Do; End Do
	umn=umn/Real(nx)/Real(nz)
	wmn=wmn/Real(nx)/Real(nz)
!
! Write the statistics out
	Open(625, File='statistics.plt', status='UNKNOWN')
 	Write(625,'(120A,120A,120A)') ' variables="y","U","W","u<sub>rms</sub>","v<sub>rms</sub>","w<sub>rms</sub>","uv","uw","vw",', &
	  '"U1","W1","u<sub>rms</sub>1","v<sub>rms</sub>1","w<sub>rms</sub>1","uv1","uw1","vw1",', &
	 '"y<sup>+</sup>"'
	Do j=1,ny
	  urms=0.0e0;vrms=0.0e0;wrms=0.0e0
	  u2=0.0e0;v2=0.0e0;w2=0.0e0
	  uv=0.0e0;uw=0.0e0;vw=0.0e0
	  Do k=1,nz; Do i=1,nx
	    u2=u2+dummy_3d(i,j,k,1)**2
	    v2=v2+dummy_3d(i,j,k,2)**2
	    w2=w2+dummy_3d(i,j,k,3)**2
	    uv=uv+dummy_3d(i,j,k,1)*dummy_3d(i,j,k,2)
	    uw=uw+dummy_3d(i,j,k,1)*dummy_3d(i,j,k,3)
	    vw=vw+dummy_3d(i,j,k,2)*dummy_3d(i,j,k,3)
	  End Do; End Do
	  u2=u2/Real(nx*nz,mytype);v2=v2/Real(nx*nz,mytype);w2=w2/Real(nx*nz,mytype)
	  uv=uv/Real(nx*nz,mytype);uw=uw/Real(nx*nz,mytype);vw=vw/Real(nx*nz,mytype)
	  urms=Sqrt(u2-umn(j)**2);vrms=Sqrt(v2);wrms=Sqrt(w2-wmn(j)**2)
	  Write(625,'(18E20.10)') yc(j),umn(j),wmn(j),urms,vrms,wrms, &
	   uv,uw-umn(j)*wmn(j),vw, &
	   umn(j)*cos(theta)-wmn(j)*sin(theta),umn(j)*sin(theta)+wmn(j)*cos(theta), &
	   Sqrt((urms*cos(theta))**2+(wrms*sin(theta))**2-(uw-umn(j)*wmn(j))*sin(2.0e0*theta)), &
	   vrms, &
	   Sqrt((urms*sin(theta))**2+(wrms*cos(theta))**2+(uw-umn(j)*wmn(j))*sin(2.0e0*theta)), &
	   uv*cos(theta)-vw*sin(theta), &
	   (urms**2-wrms**2)*sin(theta)*cos(theta)+(uw-umn(j)*wmn(j))*(cos(theta)**2-sin(theta)**2), &
	   uv*sin(theta)+vw*cos(theta), ycplus(j)
	End Do
	Close(625)
!
! Reconstruct the statistics in the other coordinate

!
! Rotate the flow field: velocity vector will change, but not pressure.
	irotate=0
	If(irotate==1) Then
	Write(*,*) 'rotating the field ... theta=',theta
	Allocate(rotfield(-1:nx+2,-1:ny+2,-1:nz+2,1:2))
	Open(625, File='restart-rotated.rot', Status='UNKNOWN', Form='UNFormatTED',access='stream')
	Write(625) nx,ny,nz,umn(1:8)
!
	rotfield(:,:,:,1)=dummy_3d(:,:,:,1);rotfield(:,:,:,2)=0.0e0
	Do k=0,nz+1; Do j=0,ny+1; Do i=0,nx+1
	  urot=0.5e0*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	  wrot=0.5e0*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	  rotfield(i,j,k,1)=urot*cos(theta)-wrot*sin(theta)
	End Do; End Do; End Do
	! make the rotated field staged
	Do k=0,nz+1; Do j=0,ny+1; Do i=0,nx+1
	  rotfield(i,j,k,2)=0.5e0*(rotfield(i-1,j,k,1)+rotfield(i,j,k,1))
	End Do; End Do; End Do
	Write(625) rotfield(:,:,:,2)
!
	Write(625) dummy_3d(:,:,:,2)
!
	rotfield(:,:,:,1)=dummy_3d(:,:,:,3);rotfield(:,:,:,2)=0.0e0
	Do k=0,nz+1; Do j=0,ny+1; Do i=0,nx+1
	  urot=0.5e0*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	  wrot=0.5e0*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	  rotfield(i,j,k,1)=urot*sin(theta)+wrot*cos(theta)
	End Do; End Do; End Do
	! make the rotated field staged
	Do k=0,nz+1; Do j=0,ny+1; Do i=0,nx+1
	  rotfield(i,j,k,2)=0.5e0*(rotfield(i,j,k-1,1)+rotfield(i,j,k,1))
	End Do; End Do; End Do
	Write(625) rotfield(:,:,:,2)
!
	Write(625) dummy_3d(:,:,:,4)
	Close(625)
	Deallocate(rotfield)
	End If
!
! Write out a snapshot
	iw=0
	If(iw==1) Then
	jmin=jspec;jmax=jspec+1
  	Open(unit=625, file='snap_rotate.vtk')
  	Write(625,'(A)') '# vtk DataFile Version 2.0'
  	Write(625,'(A)') 'Unstructured Grid Example'
  	Write(625,'(A)') 'ASCII '
  	Write(625,'(A)') 'DATASET UNSTRUCTURED_GRID'
  	Write(625,'(A,I10,A)') 'POINTS', nx*(jmax-jmin+1)*nz,' float'
	Do k=1,nz; Do j=jmin,jmax; Do i=1,nx
  	  Write(625,'(3F10.5)') xc(i)*cos(theta)-zc(k)*sin(theta),yc(j),&
	   xc(i)*sin(theta)+zc(k)*cos(theta)
	End Do; End Do; End Do
  	Write(625,'(A,2I10)') 'CELLS ', (nx-1)*(jmax-jmin)*(nz-1),9*(nx-1)*(jmax-jmin)*(nz-1)
  	Do k=1,nz-1; Do j=1,jmax-jmin; Do i=1,nx-1
    	  ig=(k-1)*nx*(jmax-jmin+1)+(j-1)*nx+i-1
    	  Write(625,'(A,8I10)') '8',ig,ig+1,ig+nx+1,ig+nx,ig+nx*(jmax-jmin+1), &
	    ig+nx*(jmax-jmin+1)+1,ig+nx*(jmax-jmin+1)+nx+1,ig+nx*(jmax-jmin+1)+nx
  	End Do; End Do; End Do
  	Write(625,'(A,I10)') 'CELL_TYPES', (nx-1)*(jmax-jmin)*(nz-1)
  	Do i=1,(nx-1)*(jmax-jmin)*(nz-1)
    	  Write(625,*) '12'
  	End Do
  	Write(625,'(A,I10)') 'POINT_DATA', nx*(jmax-jmin+1)*nz
!  	Write(625,'(A)') 'SCALARS u float 1'
!  	Write(625,'(A)') 'LOOKUP_TABLE default'
!  	Do k=1,nz; Do j=jmin,jmax; Do i=1,nx
!    	  Write(625,'(F10.5)') dummy_3d(i,j,k,1)
!  	End Do; End Do; End Do
	Write(625,*) 'VECTORS vel float'
  	Do k=1,nz; Do j=jmin,jmax; Do i=1,nx
    	  Write(625,'(F10.5)') 0.5e0*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1)), &
		0.5e0*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2)), &
		0.5e0*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
  	End Do; End Do; End Do
	Close(625)
	End If
!
!
	END SUBROUTINE dat_calc

!*******************************************************************************
	SUBROUTINE wallshearplot(dummy_3d)
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(in) :: dummy_3d
	Real(mytype), Allocatable, Dimension(:,:) :: utau_now
!
	Integer :: i,j,k,idir
!
	Allocate(utau_now(1:nxt,1:nzt))
!
	idir=2 ! 1 for x direction; 2 for z direction
!
        DO k=1,nzt; DO i=1,nxt
          utau_now(i,k)=nu*(gcdy_f(1,1,1,1)*dummy_3d(i,2,k,1)+gcdy_f(1,2,1,1)*dummy_3d(i,1,k,1))
        END DO; END DO	
!
	If(itec==1) Then
	  Open(609, File='lowerWallShear.plt', status='UNKNOWN')
	  WRITE(609,'(120A)') ' variables="x","z","<greek>t</greek><sub>w</sub>"'
	  Write(609,'(10f20.10)') ((xc(i),i=1,nx),k=1,nz)
	  Write(609,'(10f20.10)') ((zc(i),i=1,nx),k=1,nz)
	  Write(609,'(10f20.10)') ((utau_now(i,k),i=1,nx),k=1,nz)
	Else
	  Open(609, File='lowerWallShear.vtk', status='UNKNOWN')
	  Write(609,'(A)') '# vtk DataFile Version 2.0'
	  Write(609,'(A)') 'Volume example'
	  Write(609,'(A)') 'ASCII'
	  Write(609,'(A)') 'DATASET RECTILINEAR_GRID'
	  Write(609,'(A,3I10)') 'DIMENSIONS',nx,1,nz
	  Write(609,'(A,I10,A)') 'X_COORDINATES', nx, ' float'
	  Do i=1,nx; Write(609,'(E15.7)') xc(i); End Do
	  Write(609,'(A,I10,A)') 'Y_COORDINATES', 1, ' float'
	  Do j=1,1; Write(609,'(E15.7)') yc(j); End Do
	  Write(609,'(A,I10,A)') 'Z_COORDINATES', nz, ' float'
	  Do k=1,nz; Write(609,'(E15.7)') zc(k); End Do
	  Write(609,'(A,I10)') 'POINT_DATA',nx*nz
	  Write(609,'(A)') 'SCALARS tauw float 1'
	  Write(609,'(A)') 'LOOKUP_TABLE default'
	  DO k=1,nz; DO j=1,1; DO i=1,nx
  	  Write(609,'(10F20.10)') utau_now(i,k)
	  End Do; End Do; End Do
	End If
!
	Close(609)
!
	Open(626, File='lowerWallShear_1d.plt', status='UNKNOWN')
	If(idir==1) Then
	  WRITE(626,'(120A)') ' variables="x","C<sub>f</sub>"'
	  Do i=1,nx
	    Write(626,'(2F20.10)') xc(i),Sum(utau_now(i,1:nz))/Real(nz,mytype)*2.0e0
	  End Do
	Else If(idir==2) Then
	  WRITE(626,'(120A)') ' variables="z","C<sub>f</sub>"'
	  Do k=1,nz
	    Write(626,'(2F20.10)') zc(k),Sum(utau_now(1:nx,k))/Real(nx,mytype)*2.0e0
	  End Do
	End If
	Close(626)
!
	Deallocate(utau_now)
!
	END SUBROUTINE wallshearplot
!*******************************************************************************
	SUBROUTINE vort_vis(dummy_3d,wx)
!	This subroutine calculates the correlations
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(in) :: dummy_3d
	Real(mytype), Dimension(nx,ny,nz), Intent(out) :: wx
	Real(mytype), Dimension(ny,2) :: wxj,wyj,wzj,prej
	Real(mytype) :: temp,temp2,theta,ucm,uc,ucp,vcm,vc,vcp,wcm,wc,wcp
!
	Integer :: i,j,k,iwy,iwz,irotate,iwri
	Real(mytype), Dimension(:,:,:), Allocatable :: dwdy, dvdz, wy, wz
!
	Allocate(dwdy(nx,ny,nz),dvdz(nx,ny,nz))
	Allocate(wy(nx,ny,nz),wz(nx,ny,nz))
	iwy=1; iwz=1
	iwri=0 ! write files out or not
	! if the vorticity principle vector is not aligned in x direction
	irotate=0; theta=0.0e0/180.0e0*pi ! theta should be negative
	If(irotate==1) iwz=1
!
	DO k=1,nz; DO j=1,ny; DO i=1,nx
!vortices in x direction
	  wcp=0.5e0*(dummy_3d(i,j+1,k,3)+dummy_3d(i,j+1,k+1,3))
	  wc =0.5e0*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	  wcm=0.5e0*(dummy_3d(i,j-1,k,3)+dummy_3d(i,j-1,k+1,3))
	  dwdy(i,j,k)=gcdy_c(j,2,2,1)*wcp + &
		      gcdy_c(j,3,2,1)*wc + &
		      gcdy_c(j,4,2,1)*wcm
	  vcp=0.5e0*(dummy_3d(i,j,k+1,2)+dummy_3d(i,j+1,k+1,2)) 
	  vc =0.5e0*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2)) 
	  vcm=0.5e0*(dummy_3d(i,j,k-1,2)+dummy_3d(i,j+1,k-1,2)) 
	  dvdz(i,j,k)=gcdz_c(k,2,1,1)*vcp + &
		      gcdz_c(k,3,1,1)*vc + &
		      gcdz_c(k,4,1,1)*vcm 
	  wx(i,j,k)=dwdy(i,j,k)-dvdz(i,j,k)	
	End Do; End Do; End Do
	! write binary file for wx
	If(iwri==1) Then
	Open(19,File='wx-bin-snap.dat', Status='Unknown',Form='unformatted',Access='stream')
	Write(19) wx(:,:,:)
	Close(19); End If
	If(iwy==1) Then
	DO k=1,nz; DO j=1,ny; DO i=1,nx
!!vortices in y direction
	  ! dudz
	  ucp=0.5e0*(dummy_3d(i,j,k+1,1)+dummy_3d(i+1,j,k+1,1))
	  uc =0.5e0*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	  ucm=0.5e0*(dummy_3d(i,j,k-1,1)+dummy_3d(i+1,j,k-1,1))
	  dwdy(i,j,k)=gcdz_c(k,2,1,1)*ucp + &
		      gcdz_c(k,3,1,1)*uc + &
		      gcdz_c(k,4,1,1)*ucm 
	  ! dwdx
	  wcp=0.5e0*(dummy_3d(i+1,j,k,3)+dummy_3d(i+1,j,k+1,3))
	  wc =0.5e0*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	  wcm=0.5e0*(dummy_3d(i-1,j,k,3)+dummy_3d(i-1,j,k+1,3))
	  dvdz(i,j,k)=gcdx_c(i,2,1,1)*wcp + &
		      gcdx_c(i,3,1,1)*wc + &
		      gcdx_c(i,4,1,1)*wcm 
	  wy(i,j,k)=dwdy(i,j,k)-dvdz(i,j,k)
	End Do; End Do; End Do
	! write binary file for wy
	If(iwri==1) Then
	Open(19,File='wy-bin-snap.dat', Status='Unknown',Form='unformatted',Access='stream')
	Write(19) wy(:,:,:)
	Close(19); End If
	End If
	If(iwz==1) Then
	DO k=1,nz; DO j=1,ny; DO i=1,nx
!!vortices in z direction
	  ! dvdx
	  vcp=0.5e0*(dummy_3d(i+1,j,k,2)+dummy_3d(i+1,j+1,k,2))
	  vc =0.5e0*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2))
	  vcm=0.5e0*(dummy_3d(i-1,j,k,2)+dummy_3d(i-1,j+1,k,2))
	  dwdy(i,j,k)=gcdx_c(i,2,1,1)*vcp + &
		      gcdx_c(i,3,1,1)*vc + &
		      gcdx_c(i,4,1,1)*vcm 
	  ! dudy
	  ucp=0.5e0*(dummy_3d(i,j+1,k,1)+dummy_3d(i+1,j+1,k,1))
	  uc =0.5e0*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	  ucm=0.5e0*(dummy_3d(i,j-1,k,1)+dummy_3d(i+1,j-1,k,1))
	  dvdz(i,j,k)=gcdy_c(j,2,2,1)*ucp + &
		      gcdy_c(j,3,2,1)*uc + &
		      gcdy_c(j,4,2,1)*ucm 
	  wz(i,j,k)=dwdy(i,j,k)-dvdz(i,j,k)
	End Do; End Do; End Do
	! write binary file for wz
	If(iwri==1) Then
	Open(19,File='wz-bin-snap.dat', Status='Unknown',Form='unformatted',Access='stream')
	Write(19) wz(:,:,:)
	Close(19); End If
	End If
!
	! rotate the vorticity vector field to align in the principle direction
	If(irotate==1) Then
	  Do k=1,nz; DO j=1,ny; DO i=1,nx
	    temp=wx(i,j,k)*cos(theta)-wz(i,j,k)*sin(theta)
	    temp2=wx(i,j,k)*sin(theta)+wz(i,j,k)*cos(theta)
	    wx(i,j,k)=temp; wz(i,j,k)=temp2
	  End Do; End Do; End Do
	End If
!
	If(iwri==1) Then
	If(itec==1) Then
 	  WRITE(18,*) 'zone DATAPACKING=BLOCK, i=',nx,', j=',jvis,', k=',nz
	  Write(18,'(10f20.10)') (((xc(i),i=1,nx),j=1,jvis),k=1,nz)
	  Write(18,'(10f20.10)') (((yc(j),i=1,nx),j=1,jvis),k=1,nz)
	  Write(18,'(10f20.10)') (((zc(k),i=1,nx),j=1,jvis),k=1,nz)
	  Write(18,'(10f20.10)') (((wx(i,j,k),i=1,nx),j=1,jvis),k=1,nz)
!
! 	WRITE(301,*) 'zone DATAPACKING=BLOCK, i=',nx,', j=',ny,', k=',nz
!	Write(301,'(10f20.10)') (((xc(i),i=1,nx),j=1,ny),k=1,nz)
!	Write(301,'(10f20.10)') (((yc(j),i=1,nx),j=1,ny),k=1,nz)
!	Write(301,'(10f20.10)') (((zc(k),i=1,nx),j=1,ny),k=1,nz)
!	Write(301,'(10f20.10)') (((wy(i,j,k),i=1,nx),j=1,ny),k=1,nz)
!!
! 	WRITE(302,*) 'zone DATAPACKING=BLOCK, i=',nx,', j=',ny,', k=',nz
!	Write(302,'(10f20.10)') (((xc(i),i=1,nx),j=1,ny),k=1,nz)
!	Write(302,'(10f20.10)') (((yc(j),i=1,nx),j=1,ny),k=1,nz)
!	Write(302,'(10f20.10)') (((zc(k),i=1,nx),j=1,ny),k=1,nz)
!	Write(302,'(10f20.10)') (((wz(i,j,k),i=1,nx),j=1,ny),k=1,nz)
!	Close(301)
!	Close(302)
	Else
!	paraview format
	  Write(18,'(A,3I10)') 'DIMENSIONS',nx,jvis,nz
	  Write(18,'(A,I10,A)') 'X_COORDINATES', nx, ' float'
	  Do i=1,nx; Write(18,'(E15.7)') xc(i); End Do
	  Write(18,'(A,I10,A)') 'Y_COORDINATES', jvis, ' float'
	  Do j=1,jvis; Write(18,'(E15.7)') yc(j); End Do
	  Write(18,'(A,I10,A)') 'Z_COORDINATES', nz, ' float'
	  Do k=1,nz; Write(18,'(E15.7)') zc(k); End Do
	  Write(18,'(A,I10)') 'POINT_DATA',nx*jvis*nz
	  Write(18,'(A)') 'SCALARS omegax float 1'
	  Write(18,'(A)') 'LOOKUP_TABLE default'
  	  Write(18,'(10F20.10)') (((wx(i,j,k),i=1,nx),j=1,jvis),k=1,nz)
	End If; End If
	Close(18)
!calculate vorticity and pressue fluctuation: root-mean-square
	wxj=0.0; wyj=0.0; wzj=0.0;prej=0.0
	DO k=1,nz; DO j=1,ny; DO i=1,nx
	  wxj(j,1)=wxj(j,1)+wx(i,j,k)
!	  wyj(j,1)=wyj(j,1)+wy(i,j,k)
!	  wzj(j,1)=wzj(j,1)+wz(i,j,k)
	  prej(j,1)=prej(j,1)+dummy_3d(i,j,k,4)
	  wxj(j,2)=wxj(j,2)+wx(i,j,k)**2
!	  wyj(j,2)=wyj(j,2)+wy(i,j,k)**2
!	  wzj(j,2)=wzj(j,2)+wz(i,j,k)**2
	  prej(j,2)=prej(j,2)+dummy_3d(i,j,k,4)**2
	END DO; END DO; END DO	
	  temp=Real(nx*nz,kind=mytype)
	  wxj=wxj/temp
!	  wyj=wyj/temp
!	  wzj=wzj/temp
 	  prej=prej/temp
	Open(Unit=303,File='vorticity_rms.plt')
	  WRITE(303,*) ' title="vorticity_rms"'
	  WRITE(303,'(120A)') ' variables="y<sup>+</sup>","<greek>w</greek><sub>x</sub>","<greek>w</greek><sub>y</sub>","<greek>w</greek><sub>z</sub>","p/4"'
	  Do j=1,jvis
	    Write(303,'(5f20.10)') ycplus(j),sqrt(wxj(j,2)-wxj(j,1)**2)*nu/utau_avg_sa**2,sqrt(wyj(j,2)-wyj(j,1)**2)*nu/utau_avg_sa**2,&
		sqrt(wzj(j,2)-wzj(j,1)**2)*nu/utau_avg_sa**2 ,sqrt(prej(j,2)-prej(j,1)**2)/utau_avg_sa**2/rho/4.0
	  Enddo

	Close(303)

	Deallocate(dwdy,dvdz,wy,wz)
!
	END SUBROUTINE vort_vis
!*******************************************************************************
	SUBROUTINE lambdaQ(dummy_3d,l2_array,pmin,vx,wx)
!	This subroutine calculates the correlations
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(in) :: dummy_3d
	Real(mytype), Dimension(nxt,nyt,nzt), Intent(out) :: l2_array
	Integer, Dimension(nxt,nyt,nzt), Intent(out) :: pmin 
	Real(mytype), Dimension(nxt,nyt,nzt,3), Intent(out) :: vx
	Real(mytype), Dimension(nx,ny,nz), Intent(in) :: wx
!
	Integer :: i,j,k,n,m,istep,jstep,kstep,kin,jin,iwri
	Real(mytype) :: denoy,cy1,cy2,cy3,denox,cx1,cx2,cx3,denoz,cz1,cz3,cz2, &
		SM,OM,R,Q3,PM,TEMP,A1,A2,A3,l1,l2,l3,lambda,scl
	Real(mytype), Dimension(3,3) :: S, O, SO, vmat
        Real(mytype), Dimension(3) :: v1,v2,v3
   	Real(mytype), Dimension(:,:,:,:,:), Allocatable :: vel_grad
    	Real(mytype), Dimension(:,:,:), Allocatable :: det,nk,qt
	Real(mytype), Dimension(:,:,:), Allocatable :: l1_array,l2_wxneg,l2_wxpos
!
   	ALLOCATE(vel_grad(1:nxt,1:nyt,1:nzt,3,3),l1_array(nxt,nyt,nzt),l2_wxneg(nxt,nyt,nzt),l2_wxpos(nxt,nyt,nzt)) 
        !,pmin(nxt,nyt,nzt),vx(nxt,nyt,nzt,3))
    	!ALLOCATE(nk(1:nxt,1:nyt,1:nzt),qt(1:nxt,1:nyt,1:nzt),DET(1:nxt,1:nyt,1:nzt))
!	
	istep=1;jstep=1;kstep=1
        pmin=0
	iwri=0 ! write 3d flow field or not
!
         DO j=1,nyt
          denoy=(dy(j)+dy(j-1))*(dy(j+1)+dy(j))*                               &
          (dy(j+1)+2.0e0*dy(j)+dy(j-1))
          cy1=(dy(j)+dy(j-1))*(dy(j)+dy(j-1))/denoy
          cy2=(dy(j+1)-dy(j-1))*(dy(j+1)+2.0e0*dy(j)+dy(j-1))/denoy
          cy3=-(dy(j+1)+dy(j))*(dy(j+1)+dy(j))/denoy
         
	 DO i=1,nxt
          denox=(dx(i)+dx(i-1))*(dx(i+1)+dx(i))*                             &
          (dx(i+1)+2.0e0*dx(i)+dx(i-1))
          cx1=(dx(i)+dx(i-1))*(dx(i)+dx(i-1))/denox
          cx2=(dx(i+1)-dx(i-1))*(dx(i+1)+2.0e0*dx(i)+dx(i-1))/denox
          cx3=-(dx(i+1)+dx(i))*(dx(i+1)+dx(i))/denox
        
	 DO k=1,nzt
          denoz=(dz(k)+dz(k-1))*(dz(k+1)+dz(k))*                           &
          (dz(k+1)+2.0e0*dz(k)+dz(k-1))
          cz1=(dz(k)+dz(k-1))*(dz(k)+dz(k-1))/denoz
          cz2=(dz(k+1)-dz(k-1))*(dz(k+1)+2.0e0*dz(k)+dz(k-1))/denoz
          cz3=-(dz(k+1)+dz(k))*(dz(k+1)+dz(k))/denoz
        
	  vel_grad(i,j,k,3,2)=cy1*(dummy_3d(i,j+1,k,3)+dummy_3d(i,j+1,k+1,3))                               &
          +cy2*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))+cy3*(dummy_3d(i,j-1,k,3)+dummy_3d(i,j-1,k+1,3))
        
	  vel_grad(i,j,k,2,3)=cz1*(dummy_3d(i,j,k+1,2)+dummy_3d(i,j+1,k+1,2))                               &
          +cz2*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2))+cz3*(dummy_3d(i,j,k-1,2)+dummy_3d(i,j+1,k-1,2))
          
	  vel_grad(i,j,k,1,3)=cz1*(dummy_3d(i,j,k+1,1)+dummy_3d(i+1,j,k+1,1))                               &
          +cz2*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))+cz3*(dummy_3d(i,j,k-1,1)+dummy_3d(i+1,j,k-1,1))
          
	  vel_grad(i,j,k,3,1)=cx1*(dummy_3d(i+1,j,k,3)+dummy_3d(i+1,j,k+1,3))                               &
          +cx2*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))+cx3*(dummy_3d(i-1,j,k,3)+dummy_3d(i-1,j,k+1,3))
          
	  vel_grad(i,j,k,2,1)=cx1*(dummy_3d(i+1,j,k,2)+dummy_3d(i+1,j+1,k,2))                               &
          +cx2*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2))+cx3*(dummy_3d(i-1,j,k,2)+dummy_3d(i-1,j+1,k,2))
          
	  vel_grad(i,j,k,1,2)=cy1*(dummy_3d(i,j+1,k,1)+dummy_3d(i+1,j+1,k,1))                               &
          +cy2*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))+cy3*(dummy_3d(i,j-1,k,1)+dummy_3d(i+1,j-1,k,1))
!
          vel_grad(i,j,k,1,1)=(dummy_3d(i+1,j,k,1)-dummy_3d(i,j,k,1))/dx(i)
          vel_grad(i,j,k,2,2)=(dummy_3d(i,j+1,k,2)-dummy_3d(i,j,k,2))/dy(j)
          vel_grad(i,j,k,3,3)=(dummy_3d(i,j,k+1,3)-dummy_3d(i,j,k,3))/dz(k)
	ENDDO; ENDDO; ENDDO
!
	DO K=1,nz,kstep; DO j=1,ny,jstep; DO i=1,nx,istep 
	DO M=1,3
	    DO N=1,3
	      S(M,N) = 0.5E0* (vel_grad(I,J,K,M,N) + vel_grad(I,J,K,N,M) )
	      O(M,N) = 0.5E0* (vel_grad(I,J,K,M,N) - vel_grad(I,J,K,N,M) )
	    END DO
	END DO

	!SM = SQRT ( 2.0E0*( S(1,1)**2 + S(1,2)**2 + S(1,3)**2 &
	!                   +S(2,1)**2 + S(2,2)**2 + S(2,3)**2 &
	!                   +S(3,1)**2 + S(3,2)**2 + S(3,3)**2   ) )

        !OM = SQRT ( 2.0E0*( O(1,1)**2 + O(1,2)**2 + O(1,3)**2 &
	!                   +O(2,1)**2 + O(2,2)**2 + O(2,3)**2 &
	!                   +O(3,1)**2 + O(3,2)**2 + O(3,3)**2   ) )

	!R =    vel_grad(I,J,K,1,1)*(vel_grad(I,J,K,2,2)*vel_grad(I,J,K,3,3) - vel_grad(I,J,K,3,2)*vel_grad(I,J,K,2,3) ) &
        !    -  vel_grad(I,J,K,1,2)*(vel_grad(I,J,K,2,1)*vel_grad(I,J,K,3,3) - vel_grad(I,J,K,3,1)*vel_grad(I,J,K,2,3) ) &
        !    +  vel_grad(I,J,K,1,3)*(vel_grad(I,J,K,2,1)*vel_grad(I,J,K,3,2) - vel_grad(I,J,K,3,1)*vel_grad(I,J,K,2,2) )
                             
	!QT(i,j,k) = 0.5E0*(OM**2 - SM**2)
        !NK(i,j,k) = OM/SM
	!DET(i,j,k) = (QT(i,j,k)/3.0E0)**3 + (0.5E0*R)**2

	DO M=1,3
	  DO N=1,3
! SO = S**2 + O**2 by Oliver
	    SO(M,N) =   S(M,1)*S(1,N) + S(M,2)*S(2,N) + S(M,3)*S(3,N)  &
                      + O(M,1)*O(1,N) + O(M,2)*O(2,N) + O(M,3)*O(3,N)
	  END DO
	END DO

	A1= - ( SO(1,1) + SO(2,2) + SO(3,3) )

	A2 = - (  SO(1,2)*SO(2,1) + SO(1,3)*SO(3,1) + SO(2,3)*SO(3,2) &      
                - SO(1,1)*SO(2,2) - SO(1,1)*SO(3,3) - SO(2,2)*SO(3,3) )

	A3 = - (  SO(1,1)*SO(2,2)*SO(3,3) + SO(1,3)*SO(2,1)*SO(3,2) + SO(1,2)*SO(2,3)*SO(3,1) &
               - SO(1,1)*SO(2,3)*SO(3,2) - SO(1,2)*SO(2,1)*SO(3,3) - SO(1,3)*SO(2,2)*SO(3,1) )

        call cubic ( A1, A2, A3, l1, l2, l3 )
	l1_array(i,j,k)=l1
	l2_array(i,j,k)=l2
 
        v1=0.0e0; v2=0.0e0; v3=0.0e0
        v1(1)=1.0e0; v2(1)=1.0e0; v3(1)=1.0e0
        vmat=so
        call eig_vec(l1,vmat,v1)
        vmat=so
        call eig_vec(l2,vmat,v2)
        vmat=so
        call eig_vec(l3,vmat,v3)
        !v1=v1
        !write(*,*) v1,v2,v3

          !if ((v1(1).ge.v2(1)).and.(v1(1).ge.v3(1))) then
            vx(i,j,k,:)=v1!; write(*,*) 'v1'
          !else if ((v2(1).ge.v1(1)).and.(v2(1).ge.v3(1))) then
          !  vx(i,j,k,:)=v2; write(*,*) 'v2'
          !else if ((v3(1).ge.v2(1)).and.(v3(1).ge.v1(1)))then
          !  vx(i,j,k,:)=v3; write(*,*) 'v3'
          !end if
          scl=dxs(i)/vx(i,j,k,1)
          vx(i,j,k,:)=vx(i,j,k,:)*scl

          !if (k.le.64 .and. j.le.jvis .and. i.le.64) Write(25,'(6E15.7)')xc(i),yc(j),zc(k),v1(1),v1(2),v1(3)
          !if (k.le.64 .and. j.le.jvis .and. i.le.64) Write(26,'(6E15.7)')xc(i),yc(j),zc(k),v2(1),v2(2),v2(3)
          !if (k.le.64 .and. j.le.jvis .and. i.le.64) Write(27,'(6E15.7)')xc(i),yc(j),zc(k),v3(1),v3(2),v3(3)
        l2_wxpos(i,j,k)=l2_array(i,j,k)
        l2_wxneg(i,j,k)=l2_array(i,j,k)      
        if (wx(i,j,k).ge.0.0e0) l2_wxneg(i,j,k)=0.0e0
        if (wx(i,j,k).le.0.0e0) l2_wxpos(i,j,k)=0.0e0
	ENDDO; ENDDO; ENDDO

	DO k=1,nz; DO j=1,ny; DO i=1,nx

        If (l2_array(i,j,k).lt.0.0e0) then 
	  pmin(i,j,k)=1
          DO kin=-2,2; DO jin=-2,2
            If (kin.ne.0.or.jin.ne.0) then
              If ((j+jin.ge.1).and.(j+jin.le.ny)) then
                If (dummy_3d(i,j,k,4).gt.dummy_3d(i,j+jin,k+kin,4)) pmin(i,j,k)=0  
                !If (l2_array(i,j,k).ge.l2_array(i,j+jin,rmod(k+kin,nz))) pmin(i,j,k)=0  
              End If
            End If
          ENDDO; ENDDO
        End If

	ENDDO; ENDDO; ENDDO
!
	If(iwri==1) Then
	If(itec==1) Then
 	  WRITE(23,*) 'zone DATAPACKING=BLOCK, i=',nx,', j=',jvis,', k=',nz
	  Write(23,'(10E15.7)') (((xc(i),i=1,nx),j=1,jvis),k=1,nz)
	  Write(23,'(10E15.7)') (((yc(j),i=1,nx),j=1,jvis),k=1,nz)
	  Write(23,'(10E15.7)') (((zc(k),i=1,nx),j=1,jvis),k=1,nz)
	  Write(23,'(10E15.7)') (((l2_wxneg(i,j,k)*nu**2/utau_avg_sa**4,i=1,nx),j=1,jvis),k=1,nz)
	  Write(23,'(10E15.7)') (((l2_wxpos(i,j,k)*nu**2/utau_avg_sa**4,i=1,nx),j=1,jvis),k=1,nz)
	  Write(23,'(10E15.7)') (((l2_array(i,j,k)*nu**2/utau_avg_sa**4,i=1,nx),j=1,jvis),k=1,nz)
!
! 	WRITE(23,*) ' zone i=',nx,', j=',jvis,', k=',nz,',f=point'
!	DO k=1,nz; DO j=1,jvis; DO i=1,nx
!	  write(23,'(6E15.7)')xc(i),yc(j),zc(k),l2_wxneg(i,j,k)*nu**2/utau_avg_sa**4,l2_wxpos(i,j,k)*nu**2/utau_avg_sa**4,l2_array(i,j,k)*nu**2/utau_avg_sa**4
!          If (pmin(i,j,k)==1 .and. wx(i,j,k).ge.0.0e0)  then
!            !write(25,'(6E15.7)')xc(i),yc(j),zc(k),vx(i,j,k,1),vx(i,j,k,2),vx(i,j,k,3)
!            write(31,'(4E15.7)')xc(i),yc(j),zc(k),l1_array(i,j,k)*nu**2/utau_avg_sa**4!pmin(i,j,k),I3
!          End If
!	END DO; END DO; END DO
	Else
!	paraview format
	  Write(23,'(A,3I10)') 'DIMENSIONS',nx,jvis,nz
	  Write(23,'(A,I10,A)') 'X_COORDINATES', nx, ' float'
	  Do i=1,nx; Write(23,'(E15.7)') xc(i); End Do
	  Write(23,'(A,I10,A)') 'Y_COORDINATES', jvis, ' float'
	  Do j=1,jvis; Write(23,'(E15.7)') yc(j); End Do
	  Write(23,'(A,I10,A)') 'Z_COORDINATES', nz, ' float'
	  Do k=1,nz; Write(23,'(E15.7)') zc(k); End Do
	  Write(23,'(A,I10)') 'POINT_DATA',nx*jvis*nz
	  Write(23,'(A)') 'SCALARS lambda2 float 1'
	  Write(23,'(A)') 'LOOKUP_TABLE default'
  	  Write(23,'(10F20.10)') (((l2_array(i,j,k)*nu**2/utau_avg_sa**4,i=1,nx),j=1,jvis),k=1,nz)
	End If; End If
!
	END SUBROUTINE lambdaQ
!***************************************************************************
      subroutine eig_vec(lam,mat,vec)
      Real(mytype), intent(in) :: lam
      Real(mytype), Dimension(3,3), intent(inout) :: mat
      Real(mytype), Dimension(3), intent(inout) :: vec
      Real(mytype) :: ech
      integer :: m,n

      !mat=0.0e0
      !mat(1,2)=1.0e0
      !mat(1,3)=-1.0e0
      !!mat(2,1)=1.0e0
      !mat(2,2)=1.0e0
      !mat(3,1)=-1.0e0
      !mat(3,3)=1.0e0

      do m=1,3
        mat(m,m)=mat(m,m)-lam 
      end do

!     reduce to upper echelon form
      do m=2,3
        ech = mat(m,1)/mat(1,1)
        do n=1,3
          mat(m,n) = mat(m,n)-mat(1,n)*ech
        end do
      end do
      
      ech = mat(3,2)/mat(2,2)
      do n=2,3
        mat(3,n) = mat(3,n)-mat(2,n)*ech
      end do
  
      !write(*,*) 'lam:',lam
     ! do m=1,3
      !write(*,*) (mat(m,n), n=1,3)
      !end do
      !write(*,*) '******************************'
      
!      vec(3)=vec(1)*mat(1,1)/((mat(1,2)*mat(2,3)/mat(2,2))-mat(1,3))
!      vec(2)=(-vec(3)*mat(1,3)-vec(1)*mat(1,1))/mat(1,2)

!Edit by Oliver
      vec(3)=1.0
      vec(2)=-vec(3)*mat(2,3)/mat(2,2)
      vec(1)=(-vec(3)*mat(1,3)-vec(2)*mat(1,2))/mat(1,1)

      vec=vec/sqrt(vec(1)**2+vec(2)**2+vec(3)**2)

      return
      end subroutine eig_vec
!***************************************************************************
	SUBROUTINE CUBIC ( B, C, D, l1, l2, l3 )
     implicit none

     complex  :: xx(3)
     real(mytype) ::  L(3)
     integer :: nroot, I=0,J=0,K=0
     real(mytype), intent(in) :: b,c,d
     real(mytype) :: a, phi, DD, p,qq,temp1, temp2, y1,y2,y3,u,v,y2i,y2r,  pi=3.141592654
     real(mytype), intent(out) :: l1,l2,l3
	y1=0;y2=0;y3=0;y2i=0;y2r=0
      nroot = 2
        DD = c*c-4.0E0*b*d
        if(DD .ge. 0.0E0)then
          xx(1) = cmplx((-c+sqrt(DD))/2.0E0/b, 0.0E0)
          xx(2) = cmplx((-c-sqrt(DD))/2.0E0/b, 0.0E0)
        else
          xx(1) = cmplx(-c/2.0E0/b, +sqrt(-DD)/2.0E0/b)
          xx(2) = cmplx(-c/2.0E0/b, -sqrt(-DD)/2.0E0/b)
        endif


      nroot = 3
          a = 1.0E0

      p  = c/a - b*b/a/a/3.0E0
      qq  = (2.0E0*b*b*b/a/a/a - 9.0E0*b*c/a/a + 27.0E0*d/a) / 27.0E0

      DD = p*p*p/27.0E0 + qq*qq/4.0E0

      if(DD .lt. 0.0E0)then
        phi = acos(-qq/2.0E0/sqrt(abs(p*p*p)/27.0E0))
        temp1=2.0E0*sqrt(abs(p)/3.0E0)
        y1 =  temp1*cos(phi/3.0E0)
        y2 = -temp1*cos((phi+pi)/3.0E0)
        y3 = -temp1*cos((phi-pi)/3.0E0)
      else
        temp1 = -qq/2.0E0 + sqrt(DD)
        temp2 = -qq/2.0E0 - sqrt(DD)
        u = abs(temp1)**(1.0E0/3.0E0)
        v = abs(temp2)**(1.0E0/3.0E0)
        if(temp1 .lt. 0.0E0) u=-u
        if(temp2 .lt. 0.0E0) v=-v
        y1  = u + v
        y2r = -(u+v)/2.0E0
        y2i =  (u-v)*sqrt(3.0E0)/2.0E0
      endif

      temp1 = b/a/3.0E0
      y1 = y1-temp1
      y2 = y2-temp1
      y3 = y3-temp1
      y2r=y2r-temp1

      if(DD .lt. 0.0E0)then
        xx(1) = cmplx( y1,  0.0E0)
	L(1) = y1
        xx(2) = cmplx( y2,  0.0E0)
	L(2) = y2
        xx(3) = cmplx( y3,  0.0E0)
	L(3) = y3
      elseif(DD .eq. 0.)then
        xx(1) = cmplx( y1,  0.0E0)
	L(1) = y1
        xx(2) = cmplx(y2r,  0.0E0)
	L(2) = y2r
        xx(3) = cmplx(y2r,  0.0E0)
	L(3) = y2r
      else
        xx(1) = cmplx( y1,  0.0E0)
	L(1) = y1
        xx(2) = cmplx(y2r, y2i)
	L(2) = y2r
        xx(3) = cmplx(y2r,-y2i)
	L(3) = y2r
      endif

	If ((L(1) .LE. L(2)).and.(L(2) .LE. L(3))) then
	  l1=L(3); l2=L(2); l3=L(1)
        Else If ((L(3) .LE. L(2)).and.(L(2) .LE. L(1))) then
	  l1=L(1); l2=L(2); l3=L(3)
	Else If ((L(2) .LE. L(1)).and.(L(1) .LE. L(3))) then
	  l1=L(3); l2=L(1); l3=L(2)
	Else If ((L(3) .LE. L(1)).and.(L(1) .LE. L(2))) then
	  l1=L(2); l2=L(1); l3=L(3)
	Else If ((L(1) .LE. L(3)).and.(L(3) .LE. L(2))) then
	  l1=L(2); l2=L(3); l3=L(1)
	Else If ((L(2) .LE. L(3)).and.(L(3) .LE. L(1))) then
	  l1=L(1); l2=L(3); l3=L(2)
	END IF
      RETURN
 
      END SUBROUTINE CUBIC
!*******************************************************************************
        integer function rmod(n,m)
        integer :: n,m

        rmod=modulo(n,m)
        if(rmod==0) rmod=m

        end function rmod
!*******************************************************************************
	SUBROUTINE lambdaQ_mean(l2_array,l2_mean,l2n_mean,l2_rms,lam_count,lamn_count)
!
	Real(mytype), Dimension(nxt,nyt,nzt), Intent(in) ::l2_array
	Real(mytype), Dimension(ny), Intent(inout) :: l2_mean,l2_rms,l2n_mean
	Integer, Intent(inout) :: lam_count
	Integer, Dimension(ny), Intent(inout) :: lamn_count
!
	Integer :: i,j,k
!
	DO j=1,ny; DO k=1,nz; DO i=1,nx
	  l2_mean(j) = l2_mean(j) + l2_array(i,j,k)
	  l2_rms(j) = l2_rms(j) + (l2_array(i,j,k))**2.0
	  If(l2_array(i,j,k)<0.0e0) Then
	    l2n_mean(j) = l2n_mean(j) + l2_array(i,j,k)
	    lamn_count(j) = lamn_count(j) + 1
	  End If
	END DO; END DO; END DO
	lam_count = lam_count + nx*nz 
!
	END SUBROUTINE lambdaQ_mean

!***************************************************************************
	SUBROUTINE write_mean(l2_array,l2_mean,l2n_mean,l2_rms,utau_now,lam_count,lamn_count,countx,countz)
!
	Real(mytype), Dimension(nxt,nyt,nzt), Intent(in) ::l2_array
	Real(mytype), Dimension(ny), Intent(inout) :: l2_mean,l2_rms,l2n_mean
	Real(mytype), Intent(in) :: utau_now
	Integer, Intent(in) :: lam_count,countx,countz
	Integer, Dimension(ny), Intent(in) :: lamn_count
!
	Integer :: i,j,k,m
	If (ilamq==1) then
	  l2_mean=l2_mean/REAL(lam_count)
	  Do j=1,ny
	    l2n_mean(j)=l2n_mean(j)/REAL(lamn_count(j))
	  End Do
!	  l2_rms=abs(sqrt(l2_rms/REAL(lam_count)))-l2_mean**2
!Edit by Oliver
	  l2_rms=sqrt(l2_rms/REAL(lam_count)-l2_mean**2)
	  Do j=1,ny
	    Write(28,'(6E15.7)') ycplus(j),-1.0*l2_mean(j)*nu**2/utau_avg_sa**4,l2_rms(j)*nu**2/utau_avg_sa**4, & 
		-1.0*l2n_mean(j)*nu**2/utau_avg_sa**4,-1.0*(l2_mean(j)*Real(lam_count)-l2n_mean(j)*Real(lamn_count(j))) &
		/Real(lam_count-lamn_count(j))*nu**2/utau_avg_sa**4, Real(lamn_count(j))/Real(lam_count)
	  End Do
	End If
!
	If (icorr==1) then
	  corr(:,1,:) = corr(:,1,:)/REAL(countx)
	  spec(:,1,:) = spec(:,1,:)/REAL(countx)
!
	  corr(:,2,:) = corr(:,2,:)/REAL(countz)
	  spec(:,2,:) = spec(:,2,:)/REAL(countz)
!
	  corr(nxt,1,:)=corr(1,1,:);corr(nzt,2,:)=corr(1,2,:)
!
	  Do i =1,nx+1
	    Write(57,'(4E15.7)') xc(i)-xc(1),(corr(i,1,m)/corr(1,1,m),m=1,3)
	  End Do	
	  Do k =1,nz+1
	    Write(58,'(4E15.7)') zc(k)-zc(1),(corr(k,2,m)/corr(1,2,m),m=1,3)
	  End Do		
!
	  Do i =2,nx/2+1
	    Write(56,'(4E15.7)') real(i-1)*2.0e0*pi/length, (spec(i,1,m),m=1,3)
	    Write(512,'(4E15.7)') real(i-1)*2.0e0*pi/length/retau_avg_sa, (spec(i,1,m)/utau_avg_sa**2,m=1,3)
	  End Do		
	  Do k =2,nz/2+1
	    Write(59,'(4E15.7)') real(k-1)*2.0e0*pi/width, (spec(k,2,m),m=1,3)
	    Write(513,'(4E15.7)') real(k-1)*2.0e0*pi/width/retau_avg_sa, (spec(k,2,m)/utau_avg_sa**2,m=1,3)
	  End Do
	End If	
!
	END SUBROUTINE write_mean
!***************************************************************************
	SUBROUTINE write_l2rms
!calculate the correlation of -lamd2*wx,-lamd2*wy,-lamd2*wz
!
	Real(mytype), Dimension(:,:,:), Allocatable ::l2_array,wx
	Real(mytype), Dimension(:), Allocatable :: l2_mean,l2_rms,wx_mean,wx_rms,l2wx_corr,l2wxp_corr,l2wxn_corr
!
	Integer :: i,j,k,l2count,idir
	Character :: filevort*40
!
	Allocate(l2_array(1:nxt,1:nyt,1:nzt),wx(1:nx,1:ny,1:nz), &
	  l2_mean(1:ny),l2_rms(1:ny),wx_mean(1:ny),wx_rms(1:ny),l2wx_corr(1:ny), &
	  l2wxp_corr(1:ny),l2wxn_corr(1:ny))
	l2_array(:,:,:)=0.0e0; l2_mean(:)=0.0e0; l2_rms(:)=0.0e0
	wx_mean(:)=0.0e0; wx_rms(:)=0.0e0
	l2wx_corr(:)=0.0e0; l2wxp_corr(:)=0.0e0; l2wxn_corr(:)=0.0e0
!
	idir=1 ! 1 for l2-wx; 2 for l2-wy; 3 for l2-wz
	Open(405,File='../lambda2-tmp-bin-snap1.dat',Form='unformatted',Access='stream')
	Write(*,*) 'reading ../lambda2-tmp-bin-snap1.dat ...'
	Read(405) l2_array(:,:,:)
	Read(405) wx(:,:,:)
	Close(405)
	If(idir.Ne.1) Then
	  If(idir==2) filevort='../wy-bin-snap.dat'
	  If(idir==3) filevort='../wz-bin-snap.dat'
	  Open(405,File=filevort,Form='unformatted',Access='stream')
	  Write(*,*) 'reading ',Trim(filevort), ' ...'
	  Read(405) wx(:,:,:)
	  Close(405)
	End If
!
	Do j=1,ny
	  l2count=0
	  Do k=1,nz; Do i=1,nx
	    l2_mean(j)=l2_mean(j)+l2_array(i,j,k)
	    wx_mean(j)=wx_mean(j)+wx(i,j,k)
	  End Do; End Do
	  l2_mean(j)=l2_mean(j)/Real(nx*nz,mytype)
	  wx_mean(j)=wx_mean(j)/Real(nx*nz,mytype)
	  Do k=1,nz; Do i=1,nx
	    l2_rms(j)=l2_rms(j)+(l2_array(i,j,k)-l2_mean(j))**2
	    wx_rms(j)=wx_rms(j)+(wx(i,j,k)-wx_mean(j))**2
	  End Do; End Do
	  l2_rms(j)=l2_rms(j)/Real(nx*nz,mytype)
	  l2_rms(j)=Sqrt(l2_rms(j))
	  wx_rms(j)=wx_rms(j)/Real(nx*nz,mytype)
	  wx_rms(j)=Sqrt(wx_rms(j))
	End Do	
!
	Do j=1,ny
	  l2count=0
	  Do k=1,nz; Do i=1,nx
	    l2wx_corr(j)=l2wx_corr(j)-(l2_array(i,j,k)-l2_mean(j))*Abs(wx(i,j,k)-wx_mean(j))
	    If((wx(i,j,k)-wx_mean(j))>0.0e0) Then
	    l2wxp_corr(j)=l2wxp_corr(j)-l2_array(i,j,k)*(wx(i,j,k)-wx_mean(j))
	    l2count=l2count+1
	    Else
	    l2wxn_corr(j)=l2wxn_corr(j)-l2_array(i,j,k)*(wx(i,j,k)-wx_mean(j))
	    End If
	  End Do; End Do
	  l2wx_corr(j)=l2wx_corr(j)/Real(nx*nz,mytype)
	  l2wxp_corr(j)=l2wxp_corr(j)/Real(l2count,mytype)
	  l2wxn_corr(j)=l2wxn_corr(j)/Real(nx*nz-l2count,mytype)
	End Do
!
	If(ilamq.Ne.1) Then
	  Open(28, File='lambda2-rms.plt', status='UNKNOWN')
	  WRITE(28,*) ' title="lambda2-rms.plt"'
	  WRITE(28,'(120A)') ' variables="y<sup>+</sup>","-<greek>l</greek><sub>2</sub>-mean","<greek>l</greek><sub>2</sub>-rms","<greek>w</greek><sub>x</sub>-mean","<greek>w</greek><sub>x</sub>-rms","R<sub>lw</sub>","R<sub>lw,p</sub>","R<sub>lw,n</sub>"'
	End If
	Do j=1,ny
	  Write(28,'(3E15.7)') ycplus(j),-1.0*l2_mean(j)*nu**2/utau_avg_sa**4,l2_rms(j)*nu**2/utau_avg_sa**4,wx_mean(j)*nu/utau_avg_sa**2,&
	    wx_rms(j)*nu/utau_avg_sa**2,l2wx_corr(j)/l2_rms(j)/wx_rms(j),l2wxp_corr(j)/l2_rms(j)/wx_rms(j),l2wxn_corr(j)/l2_rms(j)/wx_rms(j)
	End Do
	If(ilamq.Ne.1) Close(28)

!
	Deallocate(l2_array,wx,l2_mean,l2_rms,wx_mean,wx_rms,l2wx_corr,&
	  l2wxp_corr,l2wxn_corr)
!
	END SUBROUTINE write_l2rms
!***************************************************************************
	SUBROUTINE onepoint_correlation_lambda2
!	This subroutine calculates the correlations
	Real(mytype), Dimension(:,:,:,:), Allocatable :: dummy_3d1
	Real(mytype), Dimension(:), Allocatable :: corr1
	Real(mytype), Dimension(1:5)  :: dum,derive
!	
	Integer :: i,j,k,m
	Real(mytype), Dimension(nx,nz) :: wx1,wy1,wz1,pre1,lam2
	Real(mytype) :: dpmp,d2m2,d2mp,d2mx,d2my,d2mz
        Allocate(dummy_3d1(1:nx,1:jvis,1:nz,5))
	Allocate(corr1(4))
!	
	open(unit=306,file='vorticity-lambda2.plt')
	Read(306,*) 
	Read(306,*)
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,1) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,1) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,1) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,1) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,2) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,3)
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,4)
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,5)   	
	close(306)
	open(unit=307,file='cross-corr.plt')
	WRITE(307,'(320A)') ' variables="y+","<greek>l</greek><sub>2</sub><greek>w</greek><sub>x</sub>","<greek>l</greek><sub>2</sub><greek>w</greek><sub>y</sub>","<greek>l</greek><sub>2</sub><greek>w</greek><sub>z</sub>","<greek>l</greek><sub>2</sub>p"'
	dum=0.0
	do jspec=1,jvis !y
	write(*,*) 'calculating one-point correlation at y=',yc(jspec)
	do k=1,nz; do i=1,nx
	  dum(1)=dum(1)+dabs(dummy_3d1(i,jspec,k,1))
	  dum(2)=dum(2)+dabs(dummy_3d1(i,jspec,k,2))
	  dum(3)=dum(3)+dabs(dummy_3d1(i,jspec,k,3))
	  dum(4)=dum(4)+dummy_3d1(i,jspec,k,4)
	  dum(5)=dum(5)+dummy_3d1(i,jspec,k,5)
	enddo; enddo
	dum(:)=dum(:)/real(nx*nz,kind=mytype)
!
	derive(:)=0.0
	do k=1,nz; do i=1,nx
	  derive(1)=derive(1)+(dabs(dummy_3d1(i,jspec,k,1))-dum(1))**2
	  derive(2)=derive(2)+(dabs(dummy_3d1(i,jspec,k,2))-dum(2))**2
	  derive(3)=derive(3)+(dabs(dummy_3d1(i,jspec,k,3))-dum(3))**2
	  derive(4)=derive(4)+(dummy_3d1(i,jspec,k,4)-dum(4))**2
	  derive(5)=derive(5)+(dummy_3d1(i,jspec,k,5)-dum(5))**2
	enddo; enddo
	derive(:)=sqrt(derive(:)/real(nx*nz,kind=mytype))

	DO k=1,nz; DO i=1,nx
	    wx1(i,k) = 0.5*dabs((dummy_3d1(i,jspec,k,1)+dummy_3d1(i+1,jspec,k,1)))-dum(1)
	    wy1(i,k) = 0.5*dabs((dummy_3d1(i,jspec,k,2)+dummy_3d1(i+1,jspec,k,2)))-dum(2)
	    wz1(i,k) = 0.5*dabs((dummy_3d1(i,jspec,k,3)+dummy_3d1(i+1,jspec,k,3)))-dum(3)
	    pre1(i,k) = 0.5*(dummy_3d1(i,jspec,k,4)+dummy_3d1(i+1,jspec,k,4))-dum(4)
	    lam2(i,k) = 0.5*(dummy_3d1(i,jspec,k,5)+dummy_3d1(i+1,jspec,k,5))-dum(5)
	END DO; END DO
	DO k=1,nz;DO i=1,nx!	  
	  corr1(1) = corr1(1) + lam2(i,k)*wx1(i,k)
	  corr1(2) = corr1(2) + lam2(i,k)*wy1(i,k)
	  corr1(3) = corr1(3) + lam2(i,k)*wz1(i,k)
	  corr1(4) = corr1(4) + lam2(i,k)*pre1(i,k)
        END DO;END DO
	corr1(:)=corr1(:)/real(nz*nx,kind=mytype)
	write(307,'(5E15.7)') ycplus(jspec),(corr1(m)/derive(m)/derive(5),m=1,4)
	enddo !y
	close(307)
	Deallocate(dummy_3d1,corr1)
!
	END SUBROUTINE onepoint_correlation_lambda2
!***************************************************************************
	SUBROUTINE correlation_2D_old
!	This subroutine calculates the correlations
	Real(mytype), Dimension(:,:,:,:), Allocatable :: dummy_3d1
	Real(mytype), Dimension(:,:), Allocatable :: corr1
	Real(mytype) :: dum,derive,temp
!	
	Integer :: i,j,k,m,i1,k1,ishift,kshift,iyref
	Real(mytype), Dimension(nx,nz) :: wx1,wx1_add,wxref
	Real(mytype) :: dxdx
        Allocate(dummy_3d1(1:nx,1:jvis,1:nz,5))
	Allocate(corr1(nx,nz))
!
	open(unit=306,file='vorticity-lambda2.plt')
	Read(306,*) 
	Read(306,*)
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,1) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,1) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,1) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,1) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,2) 
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,3)
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,4)
	Read(306,'(10f20.10)') dummy_3d1(:,:,:,5)   	
	close(306)

	open(unit=308,file='cross-corr.plt')
	WRITE(308,'(320A)') ' variables="x+","z+","R<greek>w</greek><sub>x</sub><greek>w</greek><sub>x</sub>"'
! 	WRITE(308,*) 'zone i=',nx,', j=',nz,',f=point'
 	WRITE(308,*) 'zone i=',nx,', j=',60,',f=point'
	iyref=25
	DO jspec=1,60 !Y direction	
	print *, 'taking correlaton in plane YCPLUS=',ycplus(jspec)
	dum=0.0
!	x direction
	do i=1,nx; do k=1,nz
	  wx1(i,k)=dummy_3d1(i,jspec,k,1)
	  wxref(i,k)=dummy_3d1(i,iyref,k,1)
	  dum=dum+dummy_3d1(i,jspec,k,1)
	enddo; enddo
	dum=dum/real(nx*nz,kind=mytype)
	derive=0.0
	do i=1,nx; do k=1,nz
	  derive=derive+(dummy_3d1(i,jspec,k,1)-dum)**2
	enddo; enddo	
	derive=sqrt(derive/real(nx*nz,kind=mytype))
!
	wx1(:,:) = wx1(:,:)-dum
	DO kshift=0,0 !-nz/2+1,nz/2
	  DO ishift=-nx/2+1,nx/2
	    dxdx=0.0
	    do k=1,nz; do i=1,nx
	      temp=mod(i+ishift,nx)
    	      if(temp.eq.0) then
		i1=nx
	      elseif(temp.gt.0) then
		i1=temp
	      elseif(temp.lt.0) then
		i1=nx+temp
	      else
		print *, 'wrong with ishift'
	      endif
	      temp=mod(k+kshift,nz)
    	      if(temp.eq.0) then
		k1=nz
	      elseif(temp.gt.0) then
		k1=temp
	      elseif(temp.lt.0) then
		k1=temp+nz
	      else
		print *,'wrong with kshift'
	      endif
	      wx1_add(i,k)=wx1(i1,k1)
	    enddo; enddo
	    do k=1,nz; do i=1,nx
	      dxdx=dxdx+wxref(i,k)*wx1_add(i,k)
	    enddo; enddo
!	    corr1(ishift+nx/2,kshift+nz/2) = dxdx/real(nx*nz,kind=mytype)/derive**2 
	    corr1(ishift+nx/2,jspec) = dxdx/real(nx*nz,kind=mytype)/derive**2 !yz xy
	  END DO
	END DO

!	do k=1,nz; do i=1,nx
!	  write(308,'(6E15.7)') xcplus(i)-xcplus(nx/2),zcplus(k)-zcplus(nz/2),corr1(i,k)
!	enddo; enddo

	do i= 1,nx
!	  write(308,'(3E15.7)') zcplus(i)-zcplus(nz/2),ycplus(jspec),corr1(i,jspec) !yz
	  write(308,'(3E15.7)') xcplus(i)-xcplus(nx/2),ycplus(jspec),corr1(i,jspec) !xy
	enddo

	ENDDO !Y direction
	close(308)
	Deallocate(dummy_3d1,corr1)
!
	END SUBROUTINE correlation_2D_old
!***************************************************************************
	SUBROUTINE correlation_2D(corrm1,corrm2,corr2d)
!	This subroutine calculates the correlations: corrm1 and corrm2 have mean components
	Real(mytype), Dimension(-1:nx+2,-1:nz+2), Intent(In) :: corrm1,corrm2
	Real(mytype), Dimension(-1:nx+2,-1:nz+2), Intent(Out) :: corr2d
!	
	Integer :: i,j,k,m,i1,k1,ishift,kshift,iyref
	Real(mytype), Dimension(:,:), Allocatable :: corrm2s
	Real(mytype) :: dxdx,rms1,rms2,mean1,mean2
	Real(mytype) :: dum,derive,temp

	Allocate(corrm2s(-1:nx+2,-1:nz+2))
!
	corr2d=0.0e0; corrm2s=0.0e0
	mean1=0.0e0; mean2=0.0e0
	Do k=1,nz; Do i=1,nx
	  mean1=mean1+corrm1(i,k)
	  mean2=mean2+corrm2(i,k)
	End Do; End Do
	mean1=mean1/Real(nx*nz,mytype)
	mean2=mean2/Real(nx*nz,mytype)
	rms1=0.0e0; rms2=0.0e0
	Do k=1,nz; Do i=1,nx
	  rms1=rms1+(corrm1(i,k)-mean1)**2
	  rms2=rms2+(corrm2(i,k)-mean2)**2
	End Do; End Do
	rms1=Sqrt(rms1/Real(nx*nz,mytype))
	rms2=Sqrt(rms2/Real(nx*nz,mytype))

!
	DO kshift=0,0 !-50,50  !-nz/2+1,nz/2
	  DO ishift=0,0 !-100,100  !-nx/2+1,nx/2
	    dxdx=0.0
	    do k=1,nz; do i=1,nx
	      temp=mod(i+ishift,nx)
    	      if(temp.eq.0) then
		i1=nx
	      elseif(temp.gt.0) then
		i1=temp
	      elseif(temp.lt.0) then
		i1=nx+temp
	      else
		print *, 'wrong with ishift'
	      endif
	      temp=mod(k+kshift,nz)
    	      if(temp.eq.0) then
		k1=nz
	      elseif(temp.gt.0) then
		k1=temp
	      elseif(temp.lt.0) then
		k1=temp+nz
	      else
		print *,'wrong with kshift'
	      endif
	      corrm2s(i,k)=corrm2(i1,k1)-mean2
	    enddo; enddo
	    do k=1,nz; do i=1,nx
	      dxdx=dxdx+(corrm1(i,k)-mean1)*corrm2s(i,k)
	    enddo; enddo
	    corr2d(ishift+nx/2,kshift+nz/2) = dxdx/real(nx*nz,kind=mytype) !/rms1/rms2
	  END DO
	END DO
!
!! 	write file for correlation
!	Open(Unit=308,File='cross-corr.plt')
!	WRITE(308,'(320A)') ' variables="x<sup>+</sup>","z<sup>+</sup>","R"'
! 	WRITE(308,*) 'zone i=',nx,', j=',nz,',f=point'
!	Do k=1,nz; Do i=1,nx
!	  Write(308,'(3E15.7)') xcplus(i)-xcplus(nx/2),zcplus(k)-zcplus(nz/2),corr2d(i,k)
!	End Do; End Do
!	Close(308)
!
	Deallocate(corrm2s)
!
	END SUBROUTINE correlation_2D
!***************************************************************************
	SUBROUTINE correlation_2D_uv
!	This subroutine calculates the correlations for uv at different y locations using multiple restart files
	Real(mytype), Dimension(:,:,:), Allocatable :: uvw,uvw2
	Real(mytype), Dimension(:,:), Allocatable :: corrm1,corrm2,corr2d,corrm2s
	Real(mytype), Dimension(:), Allocatable :: mean1,mean2
	Integer, Dimension(:), Allocatable :: xzshift(:,:)
!	
        Real(mytype), dimension(8) :: dum
 	integer, dimension(3) :: intdum
	Integer :: i,j,k,m,i1,k1,ishift,kshift,iyref,jc1,jc2,nsp,ibt
	Real(mytype) :: dxdx,temp
	Character :: chdum*8
	Character*60, Dimension(:), Allocatable :: file3dw

	Allocate(uvw(-1:nx+2,-1:ny+2,-1:nz+2),uvw2(-1:nx+2,-1:ny+2,-1:nz+2))
	Allocate(corrm1(-1:nx+2,-1:nz+2),corrm2(-1:nx+2,-1:nz+2),corr2d(-1:nx+2,-1:nz+2),corrm2s(-1:nx+2,-1:nz+2))
	Allocate(mean1(1:nz),mean2(1:nz),xzshift(1:ny,1:2))
!
	jc1=26; jc2=26; nsp=1
	Allocate(file3dw(1:nsp))
	file3dw=(/'../restrt-start0-stand-run8-out.rst'/)
	mean1=0.0e0; mean2=0.0e0
	xzshift(:,:)=0
!
	Open(Unit=308,File='cross-corr.plt')
	WRITE(308,'(320A)') ' variables="x<sup>+</sup>","z<sup>+</sup>","R"'
	Open(Unit=633,File='loc_corruv.plt')
	WRITE(633,'(320A)') ' variables="y<sup>+</sup>","x<sup>+</sup>","z<sup>+</sup>","R<sub>max</sub>","R"'
!
	Do m=1,nsp
	  Write(chdum,'(I5)') m
	  Write(*,*) ''
          Open(Unit=28, file=file3dw(m), Status='OLD', Form='UNFormatTED',access='stream')
	  Read(28) intdum(:),dum(:)
	  Write(*,*) intdum(:),dum(:)
	  Read(28) uvw(:,:,:)
	  Read(28) uvw2(:,:,:)
	  Close(28)
	Do ibt=1,2 ! top and bottom average
!get correlation at each y location
	WRITE(633,*) 'zone t="',m,'"'
	Do j=1,ny/2
	  Write(*,*) 'j=',j
	  jc1=j; jc2=j	
	  If(ibt==1) Then
	    corrm1(:,:)=uvw(:,jc1,:)
	    mean1(1:nz)=flucs(jc1,1,1)
	  Else
	    corrm1(:,:)=uvw(:,nyt-jc1,:)
	    mean1(1:nz)=flucs(nyt-jc1,1,1)
	  End If
	  If(ibt==1) Then
	    corrm2(:,:)=uvw2(:,jc2,:)
	    mean2(1:nz)=flucs(jc2,2,1)
	  Else
	    corrm2(:,:)=uvw2(:,nyt-jc2,:)
	    mean2(1:nz)=flucs(nyt-jc2,2,1)
	  End If
!
	corr2d=0.0e0; corrm2s=0.0e0
!
	DO kshift=-10,10  !-50,50  !-nz/2+1,nz/2
	  DO ishift=-15,15  !-100,100  !-nx/2+1,nx/2
	    dxdx=0.0
	    do k=1,nz; do i=1,nx
	      temp=mod(i+ishift,nx)
    	      if(temp.eq.0) then
		i1=nx
	      elseif(temp.gt.0) then
		i1=temp
	      elseif(temp.lt.0) then
		i1=nx+temp
	      else
		print *, 'wrong with ishift'
	      endif
	      temp=mod(k+kshift,nz)
    	      if(temp.eq.0) then
		k1=nz
	      elseif(temp.gt.0) then
		k1=temp
	      elseif(temp.lt.0) then
		k1=temp+nz
	      else
		print *,'wrong with kshift'
	      endif
	      corrm2s(i,k)=corrm2(i1,k1)-mean2(k1)
	    enddo; enddo
	    do k=1,nz; do i=1,nx
	      dxdx=dxdx+(corrm1(i,k)-mean1(k))*corrm2s(i,k)
	    enddo; enddo
	    corr2d(ishift+nx/2,kshift+nz/2) = dxdx/real(nx*nz,kind=mytype)
	  END DO
	END DO
	If(ibt==2) corr2d(:,:)=-corr2d(:,:)
! find the position of the minimum of uv
	xzshift(j,1:2)=Maxloc(Abs(corr2d(1:nx,1:nz)))
	Write(633,'(5F20.10)') ycplus(j),xcplus(xzshift(j,1))-xcplus(nx/2),zcplus(xzshift(j,2))-zcplus(nz/2),Maxval(Abs(corr2d(1:nx,1:nz)))/utau_avg_sa**2,corr2d(nx/2,nz/2)/utau_avg_sa**2
!
!! 	write file for correlation
	If(j==jspec) Then
 	WRITE(308,*) 'zone t="',m,'", i=',nx,', j=',nz,',f=point'
	Do k=1,nz; Do i=1,nx
	  Write(308,'(3E15.7)') xcplus(i)-xcplus(nx/2),zcplus(k)-zcplus(nz/2),corr2d(i,k)/utau_avg_sa**2
	End Do; End Do
	End If
!
	End Do ! j
	End Do ! ibt
	End Do !m=1,nsp
!
	Close(308); Close(633)
!
	Deallocate(uvw,uvw2,corrm1,corrm2,corr2d,corrm2s,mean1,mean2,xzshift)
!
	END SUBROUTINE correlation_2D_uv
!***************************************************************************
	SUBROUTINE correlation(dummy_3d,countx,countz,umn,wmn)
!	This subroutine calculates the correlations
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(in) :: dummy_3d
	Real(mytype), Dimension(ny), Intent(in)  :: umn,wmn
	Integer, Intent(inout) :: countx,countz
!
	Integer :: i,j,k
	Real(mytype), Dimension(nx) :: ux,vx,wx,ux_add,vx_add,wx_add,dumux,dumvx,dumwx
	Real(mytype), Dimension(nz) :: uz,vz,wz,uz_add,vz_add,wz_add,dumuz,dumvz,dumwz
	Real(mytype) :: dumu,dumv,dumw
!
	WRITE(*,*) 'Taking corrolation and spectrum at y=',y(jspec),'yp=',yplus(jspec)
!
	countx=0; countz=0; corr(:,:,:)=0.0e0; spec(:,:,:)=0.0e0
!	x direction
!
	DO k=1,nz
	  DO i=1,nx
	    ux(i) = dummy_3d(i,jspec,k,1)-umn(jspec)
	    vx(i) = dummy_3d(i,jspec,k,2)
	    wx(i) = dummy_3d(i,jspec,k,3)-wmn(jspec)
	  END DO
	  DO i=0,nx-1
	    CALL vec_shift(nx,i,ux,ux_add)
	    CALL vec_shift(nx,i,vx,vx_add)
	    CALL vec_shift(nx,i,wx,wx_add)
!
	    CALL multav(nx,ux,ux_add,dumu)
	    CALL multav(nx,vx,vx_add,dumv)
	    CALL multav(nx,wx,wx_add,dumw)
!		    
	    corr(i+1,1,1) = corr(i+1,1,1) + dumu
	    corr(i+1,1,2) = corr(i+1,1,2) + dumv
	    corr(i+1,1,3) = corr(i+1,1,3) + dumw
!
	    dumux(i+1) = dumu
	    dumvx(i+1) = dumv
	    dumwx(i+1) = dumw
	  END DO
	  countx = countx + 1
! here are two ways to calculation the spectra
! 	  Using correlations
!	  CALL spec_do_corr(dumux,dumvx,dumwx,dumuz,dumvz,dumwz,1)
! 	  Using u
	  CALL spec_do_vel(ux,vx,wx,uz,vz,wz,1)
	END DO

! z direction
	DO i=1,nx
	  DO k=1,nz
	    uz(k) = dummy_3d(i,jspec,k,1)-umn(jspec)
	    vz(k) = dummy_3d(i,jspec,k,2)
	    wz(k) = dummy_3d(i,jspec,k,3)-wmn(jspec)
	  END DO
	  DO k=0,nz-1
	    CALL vec_shift(nz,k,uz,uz_add)
	    CALL vec_shift(nz,k,vz,vz_add)
	    CALL vec_shift(nz,k,wz,wz_add)
!
	    CALL multav(nz,uz,uz_add,dumu)
	    CALL multav(nz,vz,vz_add,dumv)
	    CALL multav(nz,wz,wz_add,dumw)
!		
	    corr(k+1,2,1) = corr(k+1,2,1) + dumu
	    corr(k+1,2,2) = corr(k+1,2,2) + dumv
	    corr(k+1,2,3) = corr(k+1,2,3) + dumw
!
	    dumuz(k+1) = dumu
	    dumvz(k+1) = dumv
	    dumwz(k+1) = dumw
	  END DO
	  countz = countz + 1
! here are two ways to calculation the spectra
! 	  Using correlations
!	  CALL spec_do_corr(dumux,dumvx,dumwx,dumuz,dumvz,dumwz,2)
! 	  Using u
	  CALL spec_do_vel(ux,vx,wx,uz,vz,wz,2)
	END DO

!
	END SUBROUTINE correlation
!*******************************************************************************
	SUBROUTINE spec_do_corr(infftxu,infftxv,infftxw,infftzu,infftzv,infftzw,dir)
!	This subroutine calculates the spectrum from the 3d data by using the correlation
	Real(mytype), Dimension(nx), Intent(in) :: infftxu,infftxv,infftxw
	Real(mytype), Dimension(nz), Intent(in) :: infftzu,infftzv,infftzw
!
	Integer :: i,dir
	Integer*8 :: planxu,planxv,planxw,planzu,planzv,planzw
	Complex(mytype), Dimension(nx/2+1) :: outfftxu,outfftxv,outfftxw
	Complex(mytype), Dimension(nz/2+1) :: outfftzu,outfftzv,outfftzw
!
        Integer FFTW_FORWARD
        Parameter (FFTW_FORWARD = -1)
        Integer FFTW_BACKWARD
        Parameter (FFTW_BACKWARD = +1)
        Integer FFTW_ESTIMATE
        Parameter (FFTW_ESTIMATE = 64)
!	
!	x-direction
	If (dir == 1) then
	Call dfftw_plan_dft_r2c_1d(planxu,nx,infftxu,outfftxu,FFTW_ESTIMATE)
	Call dfftw_plan_dft_r2c_1d(planxv,nx,infftxv,outfftxv,FFTW_ESTIMATE)
	Call dfftw_plan_dft_r2c_1d(planxw,nx,infftxw,outfftxw,FFTW_ESTIMATE)
	Call dfftw_execute(planxu)
	Call dfftw_execute(planxv)
	Call dfftw_execute(planxw)
	DO i=1,nx/2+1
		spec(i,1,1) = spec(i,1,1) + cdabs(outfftxu(i))
		spec(i,1,2) = spec(i,1,2) + cdabs(outfftxv(i))
		spec(i,1,3) = spec(i,1,3) + cdabs(outfftxw(i))
	END DO 
!	
	Else If (dir == 2) then
!	z-direction
	Call dfftw_plan_dft_r2c_1d(planzu,nz,infftzu,outfftzu,FFTW_ESTIMATE)
	Call dfftw_plan_dft_r2c_1d(planzv,nz,infftzv,outfftzv,FFTW_ESTIMATE)
	Call dfftw_plan_dft_r2c_1d(planzw,nz,infftzw,outfftzw,FFTW_ESTIMATE)
	Call dfftw_execute(planzu)
	Call dfftw_execute(planzv)
	Call dfftw_execute(planzw)
	DO i=1,nz/2+1
		spec(i,2,1) = spec(i,2,1) + cdabs(outfftzu(i))
		spec(i,2,2) = spec(i,2,2) + cdabs(outfftzv(i))
		spec(i,2,3) = spec(i,2,3) + cdabs(outfftzw(i))
	END DO
	End If
!
	END SUBROUTINE spec_do_corr
!*******************************************************************************
	SUBROUTINE spec_do_vel(infftxu,infftxv,infftxw,infftzu,infftzv,infftzw,dir)
!	This subroutine calculates the spectrum from the 3d data by using the original signal
	Real(mytype), Dimension(nx), Intent(in) :: infftxu,infftxv,infftxw
	Real(mytype), Dimension(nz), Intent(in) :: infftzu,infftzv,infftzw
!
	Integer :: i,dir
	Integer*8 :: planxu,planxv,planxw,planzu,planzv,planzw
	Complex(mytype), Dimension(nx/2+1) :: outfftxu,outfftxv,outfftxw
	Complex(mytype), Dimension(nz/2+1) :: outfftzu,outfftzv,outfftzw
!
        Integer FFTW_FORWARD
        Parameter (FFTW_FORWARD = -1)
        Integer FFTW_BACKWARD
        Parameter (FFTW_BACKWARD = +1)
        Integer FFTW_ESTIMATE
        Parameter (FFTW_ESTIMATE = 64)
!	
!	x-direction
	If (dir == 1) then
	Call dfftw_plan_dft_r2c_1d(planxu,nx,infftxu,outfftxu,FFTW_ESTIMATE)
	Call dfftw_plan_dft_r2c_1d(planxv,nx,infftxv,outfftxv,FFTW_ESTIMATE)
	Call dfftw_plan_dft_r2c_1d(planxw,nx,infftxw,outfftxw,FFTW_ESTIMATE)
	Call dfftw_execute(planxu)
	Call dfftw_execute(planxv)
	Call dfftw_execute(planxw)
	DO i=1,nx/2+1
		spec(i,1,1) = spec(i,1,1) + cdabs(outfftxu(i))**2/Real(nx,mytype)
		spec(i,1,2) = spec(i,1,2) + cdabs(outfftxv(i))**2/Real(nx,mytype)
		spec(i,1,3) = spec(i,1,3) + cdabs(outfftxw(i))**2/Real(nx,mytype)
	END DO 
!	
	Else If (dir == 2) then
!	z-direction
	Call dfftw_plan_dft_r2c_1d(planzu,nz,infftzu,outfftzu,FFTW_ESTIMATE)
	Call dfftw_plan_dft_r2c_1d(planzv,nz,infftzv,outfftzv,FFTW_ESTIMATE)
	Call dfftw_plan_dft_r2c_1d(planzw,nz,infftzw,outfftzw,FFTW_ESTIMATE)
	Call dfftw_execute(planzu)
	Call dfftw_execute(planzv)
	Call dfftw_execute(planzw)
	DO i=1,nz/2+1
		spec(i,2,1) = spec(i,2,1) + cdabs(outfftzu(i))**2/Real(nz,mytype)
		spec(i,2,2) = spec(i,2,2) + cdabs(outfftzv(i))**2/Real(nz,mytype)
		spec(i,2,3) = spec(i,2,3) + cdabs(outfftzw(i))**2/Real(nz,mytype)
	END DO
	End If
!
	END SUBROUTINE spec_do_vel
!*******************************************************************************
	SUBROUTINE multav(int_n,vec_a,vec_b,scal_c)
!	This subroutine calculates the average of the product of two vectors 
	INTEGER :: int_n
	REAL(mytype), DIMENSION(int_n) :: vec_a,vec_b
	REAL(mytype) :: scal_c
	INTEGER :: i
!
	scal_c = 0
	DO i=1,int_n
	    scal_c = vec_a(i)*vec_b(i) + scal_c
	END DO
	scal_c = scal_c/int_n
!
	END SUBROUTINE multav
!*******************************************************************************
	SUBROUTINE vec_shift(int_n,int_shift,vec_a,vec_b)
!	This subroutine calculates the shift in a vector
	INTEGER :: int_n, int_shift
	REAL(mytype), DIMENSION(int_n) :: vec_a,vec_b
	REAL(mytype) :: temp
	INTEGER :: i,j
!
	DO i=1,int_n
	    temp = mod((i+int_shift),int_n)
	    If (temp == 0) temp=int_n
	    vec_b(i) = vec_a(temp)
	END DO
!
	END SUBROUTINE vec_shift
!***************************************************************************
	SUBROUTINE correlation_spec()
!	This subroutine calculates the correlations in Fourier space
  	Include "fftw3.f"
        Real(mytype), Dimension(:,:,:), Allocatable :: uu
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxu,outfftxu0,outfftxus
	Real(mytype), Allocatable, Dimension(:,:) :: spec2d
	Real(mytype), Allocatable, Dimension(:,:) :: ux,uxinv,uxs,corr2d
	Real(mytype), Allocatable, Dimension(:,:,:) :: AM
	Real(mytype) :: cplxr,cplxi,cplxa,xshift,zshift,kappax,kappaz,rms0,rms, &
	  dum
	Integer :: i,j,k,jp,jmin,jmax,m,kx,kz,ishift,kshift,ismin,ismax,ksmin,ksmax,j1,nsnaps, &
	  intdum
        Integer*8 :: plan
        Character :: chdum*8,file3d*60
!
!	read data in
	Allocate(uu(-1:nx+2,-1:ny+2,-1:nz+2))
	nsnaps=1
        Do m=1,nsnaps
          Write(chdum,'(I5)') m
!         file3d='snap3d_'//Trim(Adjustl(chdum))//'.dat'
          file3d=restrt_file !'../restrt-start0-stand-run8-out.rst'
          If(nsnaps==1) file3d=restrt_file
          Open(Unit=14, file=file3d, Status='OLD', Form='UNFormatTED',access='stream')
          Read(14) intdum,intdum,intdum,dum,dum,dum,dum,dum,dum,dum,dum
          Write(*,'(A,I5,A,I5,F20.10)') 'snap:',m,'  header:',intdum,dum
          Read(14) uu(:,:,:) !u
          Close(14)
        End Do
!
	imodu=1;jmodu=ny/2
	jmin=1;jmax=jmodu;m=1
	ismin=0; ismax=0; ksmin=0; ksmax=0
	Allocate(ux(1:nx,1:nz),outfftxu(1:nx/2+1,1:nz),outfftxus(1:nx/2+1,1:nz),spec2d(1:nx/2+1,1:nz), &
	  uxinv(1:nx,1:nz),uxs(1:nx,1:nz),corr2d(1:nxt,1:nzt),outfftxu0(1:nx/2+1,1:nz))
        spec2d(:,:)=0.0e0; ux(:,:)=0.0e0; uxs(:,:)=0.0e0; outfftxus(:,:)=0.0e0
	If(imodu==1) Then
	  Allocate(AM(1:jmodu,1:jmodu,1:3))
	  AM(:,:,:)=0.0e0
	End If
!	xz plane field
	Open(Unit=611,File='u_snap_xz.vtk')
	Write(611,'(A)') '# vtk DataFile Version 2.0'
	Write(611,'(A)') 'Volume example'
	Write(611,'(A)') 'ASCII'
	Write(611,'(A)') 'DATASET RECTILINEAR_GRID'
	Write(611,'(A,3I10)') 'DIMENSIONS',nx,jmax-jmin+1,nz
	Write(611,'(A,I10,A)') 'X_COORDINATES', nx, ' float'
	Write(611,'(10E15.7)') xc(1:nx)
	Write(611,'(A,I10,A)') 'Y_COORDINATES', jmax-jmin+1, ' float'
	Write(611,'(10E15.7)') yc(jmin:jmax)
	Write(611,'(A,I10,A)') 'Z_COORDINATES', nz, ' float'
	Write(611,'(10E15.7)') zc(1:nz)
	Write(611,'(A,I10)') 'POINT_DATA',nx*(jmax-jmin+1)*nz
	Write(611,'(A)') 'SCALARS usmall float 1'
	Write(611,'(A)') 'LOOKUP_TABLE default'
!	1D correlation in x direction
	Open(612, File='xcorr_spec.plt', status='UNKNOWN')
	WRITE(612,*) ' title="xcorr_spec"'
	WRITE(612,*) ' variables="x","R<sub>11</sub>"'
!	1D correlation in z direction
	Open(613, File='zcorr_spec.plt', status='UNKNOWN')
	WRITE(613,*) ' title="zcorr_spec"'
	WRITE(613,*) ' variables="z","R<sub>11</sub>"'
!	2D correlation in xz plane
	Open(614, File='xzcorr_spec.plt', status='UNKNOWN')
	WRITE(614,*) ' title="xzcorr_spec"'
	WRITE(614,*) ' variables="x","z","R<sub>11</sub>"'
	WRITE(614,*) 'zone i=',ismax-ismin+1,', j=',ksmax-ksmin+1,',f=point'
!
	Open(615, File='modulation.plt', status='UNKNOWN')
	WRITE(615,*) ' title="modulation"'
	WRITE(615,*) ' variables="z","y","R<sub>11</sub>"'
	WRITE(615,*) 'zone i=',ksmax-ksmin+1,', j=',jmax-jmin+1,',f=point'
!
	Do j1=1,jmodu; Write(*,*) 'j1=',j1
	  ux(1:nx,1:nz)=uu(1:nx,j1,1:nz)
	  Call dfftw_plan_dft_r2c_2d(plan,nx,nz,ux,outfftxu0,FFTW_ESTIMATE)
	  Call dfftw_execute(plan)
	  outfftxu0(:,:)=outfftxu0(:,:)/Real(nx*nz,mytype)
!!	  remove large scales
	  Do k=1,nz; If(k<nz/2+1) Then; kz=k; Else; kz=nz-k+2; End If
	    Do i=1,nx/2+1
	    If((idef_sup.Eq.1 .And. (kz.Le.(kz_large+1)) .Or. &
	      (idef_sup.Eq.2 .And. i.Le.(kx_large+1))) .Or. &
	      (idef_sup.Eq.3 .And. kz.Le.(kz_large+1) .And. i.Le.(kx_large+1) ) .Or. &
	      (idef_sup.Eq.4 .And. (kz.Le.(kz_large+1) .Or. i.Le.(kx_large+1)) ) ) Then !kz<=kz_large
!	      outfftxu0(i,k)=0.0e0
	    Else
	      outfftxu0(i,k)=0.0e0
  	    End If
	  End Do; End Do
	spec2d(1:nx/2+1,1:nz)=Conjg(outfftxu0(1:nx/2+1,1:nz))*outfftxu0(1:nx/2+1,1:nz)
	rms0=Sqrt(Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
		Sum(spec2d(nx/2+1,1:nz)))
!
	Do jp=1,jmodu !1,j1 !jmin,jmax
          ux(1:nx,1:nz)=uu(1:nx,jp,1:nz)
	  Call dfftw_plan_dft_r2c_2d(plan,nx,nz,ux,outfftxu,FFTW_ESTIMATE)
	  Call dfftw_execute(plan)!
	  outfftxu(:,:)=outfftxu(:,:)/Real(nx*nz,mytype)
!
!!	remove large scales
	  Do k=1,nz; If(k<nz/2+1) Then; kz=k; Else; kz=nz-k+2; End If
	    Do i=1,nx/2+1
	    If((idef_sup.Eq.1 .And. (kz.Le.(kz_large+1)) .Or. &
	      (idef_sup.Eq.2 .And. i.Le.(kx_large+1))) .Or. &
	      (idef_sup.Eq.3 .And. kz.Le.(kz_large+1) .And. i.Le.(kx_large+1) ) .Or. &
	      (idef_sup.Eq.4 .And. (kz.Le.(kz_large+1) .Or. i.Le.(kx_large+1)) ) ) Then !kz<=kz_large
!	      outfftxu(i,k)=0.0e0
	    Else
	      outfftxu(i,k)=0.0e0
  	    End If
!	helix angle
	  outfftxu(1,1)=0.0e0
!	  If(k.Ge.nz/2+1) outfftxu(i,k)=0.0e0
	  End Do; End Do
!
	  spec2d(1:nx/2+1,1:nz)=Conjg(outfftxu(1:nx/2+1,1:nz))*outfftxu(1:nx/2+1,1:nz)
	  rms=Sqrt(Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
		Sum(spec2d(nx/2+1,1:nz)))
!
	  corr2d(:,:)=0.0e0
	  Do kshift=ksmin,ksmax; Do ishift=ismin,ismax
	  xshift=length*Real(ishift)/Real(nx); zshift=width*Real(kshift)/Real(nz)
          Do k=1,nz; 
            If(k<nz/2+1) Then; kz=k-1; Else; kz=-(nz-k+1); End If
	    kappaz=2.0*pi/width*Real(kz)
	    Do i=1,nx/2+1
            kappax=2.0e0*pi/length*Real(i-1)           
!
!!            cplxr=Real(outfftxu(i,k)); cplxi=Aimag(outfftxu(i,k))
!!            cplxa=Atan2(cplxi,cplxr)
!!            outfftxus(i,k)=Abs(outfftxu(i,k))*cmplx(cos(cplxa+kappax*xshift+kappaz*zshift), &
!!                  sin(cplxa+kappax*xshift+kappaz*zshift))
            outfftxus(i,k)=outfftxu(i,k)*cmplx(cos(kappax*xshift+kappaz*zshift), &
                  sin(kappax*xshift+kappaz*zshift))
          End Do; End Do
!
!	calculate the coviarance
	spec2d(1:nx/2+1,1:nz)=Conjg(outfftxus(1:nx/2+1,1:nz))*outfftxu0(1:nx/2+1,1:nz)
	corr2d(ishift+nx/2+1,kshift+nz/2+1)=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
		Sum(spec2d(nx/2+1,1:nz))
!
!	Write inversed flow field
	  Call dfftw_plan_dft_c2r_2d(plan,nx,nz,outfftxus,ux,FFTW_ESTIMATE)
	  Call dfftw_execute(plan)
!	  Do k=1,nz; Do i=1,nx
!	    Write(611,'(10E15.7)') ux(i,k)
!	  End Do; End Do
!
	  End Do; End Do
!
	  If(imodu==1) Then
	    AM(jp,j1,1)=Maxval(corr2d(1:nx,1:nz))
	    AM(jp,j1,2:3)=Maxloc(corr2d(1:nx,1:nz))
	  End If
!	
!	  Do ishift=ismin,ismax
!	    Write(612,*) xc(ishift+nx/2+1)-xc(nx/2+1),corr2d(ishift+nx/2+1,nz/2+1)/rms0/rms
!	  End Do
!	  Do kshift=ksmin,ksmax
!	    Write(613,*) zc(kshift+nz/2+1)-zc(nz/2+1),corr2d(nx/2+1,kshift+nz/2+1)/rms0/rms
!	  End Do
!	  Do kshift=ksmin,ksmax; Do ishift=ismin,ismax
!	    Write(614,*) xc(ishift+nx/2+1)-xc(nx/2+1),zc(kshift+nz/2+1)-zc(nz/2+1), &
!		corr2d(ishift+nx/2+1,kshift+nz/2+1)/rms0/rms  
!	  End Do; End Do
!
!!	  Do ishift=ismin,ismax
!	  Do kshift=ksmin,ksmax
!!	    Write(615,*) xc(ishift+nx/2+1)-xc(nx/2+1),yc(jp), &
!!		corr2d(ishift+nx/2+1,nz/2+1)/rms0/rms 
!	    Write(615,*) zc(kshift+nz/2+1)-zc(nz/2+1),yc(jp), &
!		corr2d(nx/2+1,kshift+nz/2+1)/rms0/rms 
!	  End Do
!
	End Do;	End Do
!
	If(imodu==1) Call write_modulation(AM)
!
	Close(611); Close(612); Close(613); Close(614); Close(615)

	Deallocate(spec2d,outfftxu,outfftxu0,outfftxus,ux,uxinv,uxs,corr2d,uu)
!
	END SUBROUTINE correlation_spec
!
!***************************************************************************
	SUBROUTINE super_condition()
!	This subroutine calculates the properties conditioned by positive and negative super structures.
  	Include "fftw3.f"
!
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxu
	Real(mytype), Allocatable, Dimension(:,:) :: ux,vx,uxinv,pdf1dx,pdf1dy,pdf2d,pdf2dtmp,uxsuper
	Real(mytype), Allocatable, Dimension(:,:,:) :: ux_3d,uxinv_3d
	Real(mytype), Allocatable, Dimension(:) :: um,ump,umn,urms,urmsp,urmsn,&
		wxm,wxmp,wxmn,wxrms,wxrmsp,wxrmsn,pdf1d,pdfx,pdfy
	Integer :: i,j,k,jp,jmin,jmax,m,kx,kz,ifile,np,irangefix,inverse,ic
	Integer, Allocatable, Dimension(:) :: pcount,ncount
        Real(mytype), dimension(8) :: dum
	Real(mytype) :: tmp
 	integer, dimension(3) :: intdum
        Integer*8 :: plan
	Character*30 :: file3d,chdum
!
	Allocate(ux(1:nx,1:nz),outfftxu(1:nx/2+1,1:nz), &
	  uxinv(1:nx,1:nz),uxinv_3d(-1:nx+2,-1:ny+2,-1:nz+2),vx(1:nx,1:nz),uxsuper(1:nx,1:nz))
	Allocate(um(1:ny),ump(1:ny),umn(1:ny),urms(1:ny),urmsp(1:ny),urmsn(1:ny),&
		wxm(1:ny),wxmp(1:ny),wxmn(1:ny),wxrms(1:ny),wxrmsp(1:ny),wxrmsn(1:ny))
	Allocate(pcount(1:ny),ncount(1:ny))
	uxinv_3d(:,:,:)=0.0e0
	m=2;inverse=0
!for pdf
	jmin=1;jmax=ny/2;np=50;irangefix=0
	Allocate(pdf1dy(1:np,1:jmax-jmin+1),pdf1dx(1:np,1:jmax-jmin+1),pdf1d(1:np))
	Allocate(pdfx(1:np),pdfy(1:np),pdf2d(1:np,1:np),pdf2dtmp(1:np,1:np))
	pdf1dy(:,:)=0.0e0;pdf1dx(:,:)=0.0e0
	pdfx(:)=0.0e0;pdfy(:)=0.0e0;pdf2d(:,:)=0.0e0
!
	um(:)=0.0e0;ump(:)=0.0e0;umn(:)=0.0e0
	urms(:)=0.0e0;urmsp(:)=0.0e0;urmsn(:)=0.0e0
	wxm(:)=0.0e0;wxmp(:)=0.0e0;wxmn(:)=0.0e0
	wxrms(:)=0.0e0;wxrmsp(:)=0.0e0;wxrmsn(:)=0.0e0
	pcount(:)=0;ncount(:)=0
!
	Do ifile=ifile1,ifile2
	Allocate(ux_3d(-1:nx+2,-1:ny+2,-1:nz+2))
        Write(chdum,'(I5)') ifile
!        restrt_file='../restrt-w0.0-k0.0-run'//Trim(Adjustl(chdum))//'-out_spectra_new.rst'
	restrt_file='../../restrt-w0.0-k0.0-run5-out-edd.rst'
	Open(616, File=restrt_file, Status='OLD', Form='UNFormatTED',access='stream')
	Read(616) intdum(:),dum(:)
	write(*,*) intdum, dum
	Do ic=1,m
	  print *,'reading component', ic, 'from', restrt_file
	  Read(616) ux_3d(:,:,:)
	  If(ic==1) Then !read data for uv pdf at jspec
	    uxinv_3d(:,:,:)=ux_3d(:,:,:)
	  Else If(ic==2) Then
	    vx(1:nx,1:nz)=0.5e0*(ux_3d(1:nx,jspec,1:nz)+ux_3d(1:nx,jspec+1,1:nz))
	  End If
	End Do
	Close(616)
	ux_3d(:,:,:)=uxinv_3d(:,:,:)
!
	If(inverse==1) Then
	Write(*,*) 'geting super scales'
	Do jp=1,ny
	  ux(1:nx,1:nz)=ux_3d(1:nx,jp,1:nz)
	  Call dfftw_plan_dft_r2c_2d(plan,nx,nz,ux,outfftxu,FFTW_ESTIMATE)
	  Call dfftw_execute(plan)!
	  outfftxu(:,:)=outfftxu(:,:)/Real(nx*nz,mytype)
!!	remove large scales
	  Do k=1,nz; If(k<nz/2+1) Then; kz=k; Else; kz=nz-k+2; End If
	    Do i=1,nx/2+1
	    If((idef_sup.Eq.1 .And. (kz.Le.(kz_large+1)) .Or. &
	      (idef_sup.Eq.2 .And. i.Le.(kx_large+1))) .Or. &
	      (idef_sup.Eq.3 .And. kz.Le.(kz_large+1) .And. i.Le.(kx_large+1) ) .Or. &
	      (idef_sup.Eq.4 .And. (kz.Le.(kz_large+1) .Or. i.Le.(kx_large+1)) ) ) Then !kz<=kz_large
!	      outfftxu(i,k)=0.0e0
	    Else
	      outfftxu(i,k)=0.0e0
  	    End If
	  End Do; End Do
	  outfftxu(1,1)=0.0e0
!
	  Call dfftw_plan_dft_c2r_2d(plan,nx,nz,outfftxu,uxinv,FFTW_ESTIMATE)
	  Call dfftw_execute(plan)
	  uxinv_3d(1:nx,jp,1:nz)=uxinv(1:nx,1:nz)
	End Do
	Cycle
	End If !inverse
!
	file3d='usuper-bin-snap'//Trim(Adjustl(chdum))//'.dat'
!	file3d='../usuper-bin-snap'//Trim(Adjustl(chdum))//'.dat'
	Open(615,File=file3d, Status='Unknown',Form='unformatted',Access='stream')
!	Write(615) uxinv_3d(:,:,:)
	Read(615) uxinv_3d(:,:,:)
	Close(615)
!
!conditioned mean and rms velocity
	Do j=1,ny
	  Do k=1,nz; Do i=1,nx
	    um(j)=um(j)+ux_3d(i,j,k)
	    urms(j)=urms(j)+ux_3d(i,j,k)**2
	    If(uxinv_3d(i,j,k).Gt.0.0e0) Then
	      ump(j)=ump(j)+ux_3d(i,j,k)
	      urmsp(j)=urmsp(j)+ux_3d(i,j,k)**2
	      pcount(j)=pcount(j)+1
	    Else
	      umn(j)=umn(j)+ux_3d(i,j,k)
	      urmsn(j)=urmsn(j)+ux_3d(i,j,k)**2
	      ncount(j)=ncount(j)+1
	    End If
	  End Do; End Do
	End Do
!
!1d pdf
	Do j=jmin,jmax
	  jp=j-jmin+1
!	  Call pdf_1d((ux_3d(1:nx,j,1:nz)-flucs(j,1,1))/utau_avg_sa,1,np,pdf1dx(1:np,jp),pdf1d(1:np))
!
	  ux(1:nx,1:nz)=(ux_3d(1:nx,j,1:nz)-flucs(j,1,1))/utau_avg_sa
	  uxsuper(1:nx,1:nz)=uxinv_3d(1:nx,j,1:nz)
	  Call pdf_1d_super(ux,uxsuper,1,np,pdf1dx(1:np,jp),pdf1d(1:np))
	  pdf1dy(1:np,jp)=pdf1dy(1:np,jp)+pdf1d(1:np)
	End Do
!2d pdf
	Do i=1,nx
	ux(i,1:nz)=((0.5e0*(ux_3d(i,jspec,1:nz)+ux_3d(i+1,jspec,1:nz)))-flucs(jspec,1,1))/utau_avg_sa
	End Do
	vx(1:nx,1:nz)=(vx(1:nx,1:nz)-flucs(jspec,2,1))/utau_avg_sa
	uxsuper(1:nx,1:nz)=uxinv_3d(1:nx,jspec,1:nz)
	Call pdf_2d_super(ux,vx,uxsuper,1,np,pdfx,pdfy,pdf2dtmp)
!	Call pdf_2d(ux,vx,0,np,pdfx,pdfy,pdf2dtmp)
	pdf2d(1:np,1:np)=pdf2d(1:np,1:np)+pdf2dtmp(1:np,1:np)
!
	Deallocate(ux_3d)
	Allocate(ux_3d(1:nx,1:ny,1:nz))
	Write(chdum,'(I5)') ifile
	file3d='wx-bin-snap'//Trim(Adjustl(chdum))//'.dat'
!	file3d='../wx-bin-snap'//Trim(Adjustl(chdum))//'.dat'
	Open(615,File=file3d, Status='Unknown',Form='unformatted',Access='stream')
	Read(615) ux_3d(:,:,:)
	Close(615)
!conditioned vorticity rms
	Do j=1,ny
	  Do k=1,nz; Do i=1,nx
	    wxm(j)=wxm(j)+ux_3d(i,j,k)
	    wxrms(j)=wxrms(j)+ux_3d(i,j,k)**2
	    If(uxinv_3d(i,j,k).Gt.0.0e0) Then
	      wxmp(j)=wxmp(j)+ux_3d(i,j,k)
	      wxrmsp(j)=wxrmsp(j)+ux_3d(i,j,k)**2
	    Else
	      wxmn(j)=wxmn(j)+ux_3d(i,j,k)
	      wxrmsn(j)=wxrmsn(j)+ux_3d(i,j,k)**2
	    End If
	  End Do; End Do
	End Do
!
	Deallocate(ux_3d)
	End Do !ifile
!
	Do j=1,ny
	  um(j)=um(j)/Real(nx*nz*(ifile2-ifile1+1))
	  ump(j)=ump(j)/Real(pcount(j));umn(j)=umn(j)/Real(ncount(j))
	  urms(j)=urms(j)/Real(nx*nz*(ifile2-ifile1+1))
	  urmsp(j)=urmsp(j)/Real(pcount(j));urmsn(j)=urmsn(j)/Real(ncount(j))
	  urms(j)=Sqrt(urms(j)-um(j)**2)
	  urmsp(j)=Sqrt(urmsp(j)-2.0e0*ump(j)*um(j)+um(j)**2)
	  urmsn(j)=Sqrt(urmsn(j)-2.0e0*umn(j)*um(j)+um(j)**2)
	  wxm(j)=wxm(j)/Real(nx*nz*(ifile2-ifile1+1))
	  wxmp(j)=wxmp(j)/Real(pcount(j));wxmn(j)=wxmn(j)/Real(ncount(j))
	  wxrms(j)=wxrms(j)/Real(nx*nz*(ifile2-ifile1+1))
	  wxrmsp(j)=wxrmsp(j)/Real(pcount(j));wxrmsn(j)=wxrmsn(j)/Real(ncount(j))
	  wxrms(j)=Sqrt(wxrms(j)-wxm(j)**2)
	  wxrmsp(j)=Sqrt(wxrmsp(j)-2.0e0*wxmp(j)*wxm(j)+wxm(j)**2)
	  wxrmsn(j)=Sqrt(wxrmsn(j)-2.0e0*wxmn(j)*wxm(j)+wxm(j)**2)
	End Do
!pdf
	Do j=jmin,jmax; jp=j-jmin+1; Do i=1,np
	  pdf1dy(i,jp)=pdf1dy(i,jp)/Real((ifile2-ifile1+1),mytype)
	End Do; End Do
	pdf2d(1:np,1:np)=pdf2d(1:np,1:np)/Real((ifile2-ifile1+1),mytype)
!
	Open(Unit=617,File='umean_super.plt')
	WRITE(617,'(320A)') ' variables="y<sup>+</sup>","U<sup>+</sup>","U<sup>+</sup><sub>p</sub>","U<sup>+</sup><sub>n</sub>","yc"'
 	WRITE(617,*) 'zone i=',ny,',f=point'
	Do j=1,ny
	  Write(617,'(5F10.5)') ycplus(j),um(j)/utau_avg_sa,ump(j)/utau_avg_sa,umn(j)/utau_avg_sa,yc(j)
	End Do
	Close(617)
	Open(Unit=618,File='urms_super.plt')
	WRITE(618,'(320A)') ' variables="y<sup>+</sup>","u","u<sub>p</sub>","u<sub>n</sub>","yc"'
 	WRITE(618,*) 'zone i=',ny,',f=point'
	Do j=1,ny
	  Write(618,'(5F10.5)') ycplus(j),urms(j)/utau_avg_sa,urmsp(j)/utau_avg_sa,urmsn(j)/utau_avg_sa,yc(j)
	End Do
	Close(618)
	Open(Unit=619,File='wxmean_super.plt')
	WRITE(619,'(320A)') ' variables="y<sup>+</sup>","<greek>W</greek><sub>x</sub><sup>+</sup>","<greek>W</greek><sub>x,p</sub><sup>+</sup>","<greek>W</greek><sub>x,n</sub><sup>+</sup>","yc"'
 	WRITE(619,*) 'zone i=',ny,',f=point'
	Do j=1,ny
	  Write(619,'(5F10.5)') ycplus(j),wxm(j)*nu/utau_avg_sa**2,wxmp(j)*nu/utau_avg_sa**2,wxmn(j)*nu/utau_avg_sa**2,yc(j)
	End Do
	Close(619)
	Open(Unit=620,File='wxrms_super.plt')
	WRITE(620,'(320A)') ' variables="y<sup>+</sup>","<greek>w</greek><sub>x</sub>","<greek>w</greek><sub>x,p</sub>","<greek>w</greek><sub>x,n</sub>","yc"'
 	WRITE(620,*) 'zone i=',ny,',f=point'
	Do j=1,ny
	  Write(620,'(5F10.5)') ycplus(j),wxrms(j)*nu/utau_avg_sa**2,wxrmsp(j)*nu/utau_avg_sa**2,wxrmsn(j)*nu/utau_avg_sa**2,yc(j)
	End Do
	Close(620)
!1d pdf
	Write(*,*) '1D pdf'
	Open(Unit=622,File='pdf1d_super.plt')
	WRITE(622,'(320A)') ' variables="u<sup>+</sup>","y<sup>+</sup>","pdf","yc"'
 	WRITE(622,*) 'zone i=',np,'j=',jmax-jmin+1,',f=point'
	Do j=jmin,jmax; jp=j-jmin+1; tmp=Maxval(pdf1dy(1:np,jp))
	  Do i=1,np
	  Write(622,'(5F20.10)') pdf1dx(i,jp),ycplus(j),pdf1dy(i,jp)/tmp,yc(j)
	End Do; End Do
	Close(622)
!2d pdf
	Write(*,*) '2D pdf'
	Open(Unit=622,File='pdf2d_super.plt')
	WRITE(622,'(320A)') ' variables="u<sup>+</sup>","v<sup>+</sup>","pdf"'
 	WRITE(622,*) 'zone i=',np,'j=',np,',f=point'
	Do j=1,np; Do i=1,np
	  Write(622,'(3F20.10)') pdfx(i),pdfy(j),pdf2d(i,j)
	End Do; End Do
	Close(622)
!
	Deallocate(ux,outfftxu,uxinv,uxinv_3d,vx,uxsuper)
	Deallocate(um,ump,umn,urms,urmsp,urmsn,wxrms,wxrmsp,wxrmsn,wxm,wxmp,wxmn)
	Deallocate(pcount,ncount)
	Deallocate(pdf1d,pdf1dx,pdf1dy)
	Deallocate(pdfx,pdfy,pdf2d,pdf2dtmp)
!
	END SUBROUTINE super_condition
!*******************************************************************************
!***************************************************************************
	SUBROUTINE struct_log
!	This subroutine finds the shapes of positive and negative super structures in yz plane
	Real(mytype), Allocatable, Dimension(:,:,:) :: ux_3d,vx_3d,wx_3d
	Real(mytype), Allocatable, Dimension(:,:) :: ux,uxsuper,uavgp,uavgn,uxsft,vavgp,vavgn,wavgp,wavgn,vx,wx
	Integer, Allocatable, Dimension(:) :: zcentp,zcentn,zcent,zcentmax
	Integer :: i,j,k,jmin,jmax,kmin,kmax,m,km,ic,ipos,ifile,cntp,cntn,cntc,cnt,zloc,zsft
	Character*30 :: file3d,chdum
        Real(mytype), dimension(8) :: dum
 	integer, dimension(3) :: intdum
!
	Allocate(ux(1:nx,1:nz),uxsuper(1:nx,1:nz),uavgp(1:ny,1:nz),uavgn(1:ny,1:nz), &
	uxsft(1:ny,1:nz),vavgp(1:ny,1:nz),vavgn(1:ny,1:nz),wavgp(1:ny,1:nz), &
	wavgn(1:ny,1:nz),vx(1:nx,1:nz),wx(1:nx,1:nz))
	Allocate(ux_3d(-1:nx+2,-1:ny+2,-1:nz+2),vx_3d(-1:nx+2,-1:ny+2,-1:nz+2), &
	wx_3d(-1:nx+2,-1:ny+2,-1:nz+2))
	zcent(:)=0;zcentp(:)=0;zcentn(:)=0
	uavgp(:,:)=0.0e0;uavgn(:,:)=0.0e0
	vavgp(:,:)=0.0e0;vavgn(:,:)=0.0e0
	wavgp(:,:)=0.0e0;wavgn(:,:)=0.0e0
	cntp=0;cntn=0;cntc=20
	Allocate(zcent(1:cntc),zcentmax(1:cntc),zcentp(1:cntc),zcentn(1:cntc))
	jmin=1;jmax=ny;kmin=1;kmax=nz;zloc=nz/2
!
	m=1;ipos=1
	Write(*,*) 'get outer structure at y+=',ycplus(jspec)
	Do ifile=ifile1,ifile2
!read super structure field
        Write(chdum,'(I5)') ifile
	file3d='usuper-bin-snap'//Trim(Adjustl(chdum))//'.dat'
!	file3d='../usuper-bin-snap'//Trim(Adjustl(chdum))//'.dat'
	Open(615,File=file3d, Status='Unknown',Form='unformatted',Access='stream')
	Read(615) ux_3d(:,:,:)
	Close(615)
	uxsuper(1:nx,1:nz)=ux_3d(1:nx,jspec,1:nz)
!
        Write(chdum,'(I5)') ifile
!        restrt_file='../restrt-w0.0-k0.0-run'//Trim(Adjustl(chdum))//'-out_spectra_new.rst'
	restrt_file='../../restrt-w0.0-k0.0-run5-out-edd.rst'
	Open(616, File=restrt_file, Status='OLD', Form='UNFormatTED',access='stream')
	Read(616) intdum(:),dum(:)
	write(*,*) intdum, dum
	Do ic=1,m
	  print *,'reading component', ic, 'from', restrt_file
	  Read(616) ux_3d(:,:,:)
	  Read(616) vx_3d(:,:,:)
	  Read(616) wx_3d(:,:,:)
	End Do
	Close(616)
!
! positive structure
Open(Unit=625,File='line-z.plt')
Write(*,*) 'selected plane',ycplus(jspec)
	Do i=1,nx
	  Do k=1,nz; Do j=1,ny
	    ux(j,k)=0.5e0*(ux_3d(i,j,k)+ux_3d(i+1,j,k))-flucs(j,1,1)
	    vx(j,k)=0.5e0*(vx_3d(i,j,k)+vx_3d(i,j+1,k))
	    wx(j,k)=0.5e0*(wx_3d(i,j,k)+wx_3d(i,j,k+1))
	  End Do; End Do
	  cnt=0
	  Do k=1,nz
	  ! find local maximum
	    If(uxsuper(i,Mod(k-1,nz)+1)*uxsuper(i,Mod(k,nz)+1)<0.0e0) Then
	      cnt=cnt+1
	      zcent(cnt)=k
	    End If
	  End Do
	  ! find the region centre for positive and negative structures
	  Do km=1,cnt-1
	    zcentmax(km)=(zcent(km)+zcent(km+1))/2
	  End Do
	  zcentmax(cnt)=Mod((zcent(cnt)+zcent(1)+nz)/2-1,nz)+1
! shift the super streaks to the centre location in z direction
	  Do km=1,cnt
	    zsft=zloc-zcentmax(km)+nz
!Write(*,*) zcentmax(km),zsft,zloc
	    Do k=1,nz
	      uxsft(1,Mod(k+zsft-1,nz)+1)=uxsuper(1,k)
	    End Do
!
	    If(uxsuper(1,zcentmax(km))>0.0e0) Then
	      Do k=1,nz
	        uavgp(1:ny,Mod(k+zsft-1,nz)+1)=uavgp(1:ny,Mod(k+zsft-1,nz)+1)+ux(1:ny,k)
	        vavgp(1:ny,Mod(k+zsft-1,nz)+1)=vavgp(1:ny,Mod(k+zsft-1,nz)+1)+vx(1:ny,k)
	        wavgp(1:ny,Mod(k+zsft-1,nz)+1)=wavgp(1:ny,Mod(k+zsft-1,nz)+1)+wx(1:ny,k)
	      End Do
	      cntp=cntp+1
	    End If
	    If(uxsuper(1,zcentmax(km))<0.0e0) Then
	      Do k=1,nz
	        uavgn(1:ny,Mod(k+zsft-1,nz)+1)=uavgn(1:ny,Mod(k+zsft-1,nz)+1)+ux(1:ny,k)
	        vavgn(1:ny,Mod(k+zsft-1,nz)+1)=vavgn(1:ny,Mod(k+zsft-1,nz)+1)+vx(1:ny,k)
	        wavgn(1:ny,Mod(k+zsft-1,nz)+1)=wavgn(1:ny,Mod(k+zsft-1,nz)+1)+wx(1:ny,k)
	      End Do
	      cntn=cntn+1
	    End If
	  End Do
Write(*,*) 'i=',i,'cnt=',cnt,'cntp=',cntp,'cntn=',cntn
!
	End Do
!
	End Do !ifile
	uavgp(:,:)=uavgp(:,:)/Real(cntp,mytype)
	uavgn(:,:)=uavgn(:,:)/Real(cntn,mytype)
	vavgp(:,:)=vavgp(:,:)/Real(cntp,mytype)
	vavgn(:,:)=vavgn(:,:)/Real(cntn,mytype)
	wavgp(:,:)=wavgp(:,:)/Real(cntp,mytype)
	wavgn(:,:)=wavgn(:,:)/Real(cntn,mytype)

Do k=1,nz
!  Write(625,'(2F20.10)') zc(k),uxsuper(1,k)
End Do
!Write(625,*) 'ZONE'
Do km=1,cnt
k=zcentmax(km) !zcent(km)
!Write(625,'(2F20.10)') zc(k),uxsuper(1,k)
Write(*,*) uxsuper(1,Mod(k-1,nz)+1),uxsuper(1,Mod(k,nz)+1)
End Do

! Write the vtk file for yz plane
	Open(Unit=624,File='usuper_yz.vtk')
	Write(624,'(A)') '# vtk DataFile Version 2.0'
	Write(624,'(A)') 'Volume example'
	Write(624,'(A)') 'ASCII'
	Write(624,'(A)') 'DATASET RECTILINEAR_GRID'
	Write(624,'(A,3I10)') 'DIMENSIONS',1,jmax-jmin+1,kmax-kmin+1
	Write(624,'(A,I10,A)') 'X_COORDINATES', 1, ' float'
	Write(624,'(10E15.7)') 0.0e0
	Write(624,'(A,I10,A)') 'Y_COORDINATES', jmax-jmin+1, ' float'
	Write(624,'(10E15.7)') yc(jmin:jmax)
	Write(624,'(A,I10,A)') 'Z_COORDINATES', kmax-kmin+1, ' float'
	Write(624,'(10E15.7)') zc(kmin:kmax)
	Write(624,'(A,I10)') 'POINT_DATA',(jmax-jmin+1)*(kmax-kmin+1)
	Write(624,'(A)') 'VECTORS up float'
	Do k=kmin,kmax; Do j=jmin,jmax
	  Write(624,'(3E15.7)') uavgp(j,k),vavgp(j,k),wavgp(j,k)
	End Do; End Do
	Write(624,'(A)') 'VECTORS un float'
	Do k=kmin,kmax; Do j=jmin,jmax
	  Write(624,'(3E15.7)') uavgn(j,k),vavgn(j,k),wavgn(j,k)
	End Do; End Do
!	Write(624,'(A)') 'SCALARS un float 1'
!	Write(624,'(A)') 'LOOKUP_TABLE default'
!	Do k=kmin,kmax; Do j=jmin,jmax
!	  Write(624,'(1E15.7)') uavgn(j,k)
!	End Do; End Do
	Close(624)
!
	Deallocate(ux_3d,vx_3d,wx_3d)
	Deallocate(ux,uxsuper,uavgp,uavgn,zcent,zcentp,zcentn,uxsft,vavgp,vavgn,wavgp,wavgn,vx,wx)
!
	END SUBROUTINE struct_log
!*******************************************************************************
!*******************************************************************************
	Subroutine spectra2D(dummy_3d1,umn,wmn)
!plot 2D spectra in xy plane
  	Include "fftw3.f"
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(in) :: dummy_3d1
	Real(mytype), Dimension(ny), Intent(in)  :: umn,wmn
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxu
	Real(mytype), Allocatable, Dimension(:,:) :: spec2d
	Real(mytype), Allocatable, Dimension(:,:) :: ux
	Real(mytype), Allocatable, Dimension(:) :: accum,diff,lamx,lamz,kapx,kapz
	Real(mytype), Allocatable, Dimension(:,:) :: accum2,diff2
	Real(mytype), Allocatable, Dimension(:,:,:) :: tmpfield
        Real(mytype) :: wl,wlx,wlz, wn,wnx,wnz, rms,rmsl,rmss, specavg
	Integer :: i,j,k,jp,jmin,jmax,ixy,ixz,iyz,ilarge3d,m,imin,imax,kz
        Integer*8 :: plan
	Character*100 :: chtmp, chtmpp

! averaged FFT in xy plane (ixy=1); xz plane (ixz=1) and yz plane (iyz=1)
! m is the variable for FFT (1 for u; 2 for v; 3 for w and 4 for p)
	ixy=1;ixz=0;iyz=0;m=1;ilarge3d=1
! FFT range in y direction
	jmin=1;jmax=ny !jspec
	imin=islicexyz;imax=islicexyz !islicexyz
!
	Allocate(accum(1:Max(nx/2+1,nz/2+1)),diff(1:Max(nx/2+1,nz/2+1)),&
	lamx(1:nx/2+1),kapx(1:nx/2+1), &
	lamz(1:nz/2+1),kapz(1:nz/2+1))
	accum(:)=0.0e0; lamx(:)=0.0e0; kapx(:)=0.0e0; lamz(:)=0.0e0; kapz(:)=0.0e0
!
	Open(Unit=309,File='spectra_2dp.plt')
	Open(Unit=310,File='spectra_2d.plt')
	If(ixz == 1) Then; chtmp='k<sub>x</sub>k<sub>z</sub><greek>F</greek><sub>uu</sub>';
	chtmpp='k<sub>x</sub>k<sub>z</sub><greek>F</greek><sub>uu</sub>/<greek>n</greek><sup>2</sup>'; 
	Else If(ixy == 1) Then; chtmp='k<sub>x</sub><greek>F</greek><sub>uu</sub>'; 
	chtmpp='k<sub>x</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>'; 
	Else If(iyz == 1) Then; chtmp='k<sub>z</sub><greek>F</greek><sub>uu</sub>'; 
	chtmpp='k<sub>z</sub><greek>F</greek><sub>uu</sub>/u<sub><greek>t</greek></sub><sup>2</sup>'; End If
	Write(309,'(10A)') ' variables="k<sub>x</sub><sup>+</sup>","k<sub>z</sub><sup>+</sup>","y<sup>+</sup>",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk/u<sub><greek>t</greek></sub><sup>2</sup>", "<greek>F</greek><sub>uu</sub>/(u<sub><greek>t</greek></sub><greek>n</greek>)",', &
'"', Trim(chtmpp), '"'
	Write(310,'(10A)') ' variables="k<sub>x</sub>","k<sub>z</sub>","y",', &
	 '"<math>r</math><greek>F</greek><sub>uu</sub>dk", "<greek>F</greek><sub>uu</sub>", "', Trim(chtmp), '"'
!
	Allocate(ux(1:nx,1:nz),outfftxu(1:nx/2+1,1:nz),spec2d(1:nx/2+1,1:nz))
        spec2d(:,:)=0.0e0; ux(:,:)=0.0e0
!
	lamx(1)=0.0e0; kapx(1)=0.0e0
	Do i=2,nx/2+1; lamx(i)=length/real(i-1,mytype); kapx(i)=2.0e0*pi/lamx(i); End Do	
	lamz(1)=0.0e0; kapz(1)=0.0e0
	Do k=2,nz/2+1; lamz(k)=width/real(k-1,mytype);kapz(k)=2.0e0*pi/lamz(k); End Do
!
	If(ixy == 1) Then ! xy view, averaged in z direction
 	  Write(309,*) 'zone i=',nx/2+1,', j=',jmax-jmin+1,',f=point'
 	  Write(310,*) 'zone i=',nx/2+1,', j=',jmax-jmin+1,',f=point'
          Do jp=jmin,jmax !y
	    Write(*,*) 'y+=', yplus(jp)
            ux(1:nx,1:nz)=dummy_3d1(1:nx,jp,1:nz,m)
	    Call dfftw_plan_dft_r2c_2d(plan,nx,nz,ux,outfftxu,FFTW_ESTIMATE)
	    Call dfftw_execute(plan)
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxu(1:nx/2+1,1:nz))*outfftxu(1:nx/2+1,1:nz)/Real(nx*nz,mytype)/Real(nx*nz,mytype)
!
	    Call specxy(spec2d,kapx,accum,diff)
!
	    Do i=1,nx/2+1
	      Write(309,'(7E15.7)') kapx(i)/retau_avg_sa, 0.0e0, yc(jp)*retau_avg_sa, &
	accum(i)/utau_avg_sa**2, diff(i)/(utau_avg_sa*nu), diff(i)*kapx(i)/utau_avg_sa**2
	      Write(310,'(7E15.7)') kapx(i), 0.0e0, yc(jp), accum(i), diff(i), diff(i)*kapx(i)
	    End Do
	  End Do
!!
	Else If(iyz == 1) Then
 	  Write(309,*) 'zone i=',nz/2+1,', j=',jmax-jmin+1,',f=point'
 	  Write(310,*) 'zone i=',nz/2+1,', j=',jmax-jmin+1,',f=point'
          Do jp=jmin,jmax !y
	    Write(*,*) 'y+=', yplus(jp)
            ux(1:nx,1:nz)=dummy_3d1(1:nx,jp,1:nz,m)
	    Call dfftw_plan_dft_r2c_2d(plan,nx,nz,ux,outfftxu,FFTW_ESTIMATE)
	    Call dfftw_execute(plan)
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxu(1:nx/2+1,1:nz))*outfftxu(1:nx/2+1,1:nz)/Real(nx*nz,mytype)/Real(nx*nz,mytype)
!
	    Call specyz(spec2d,kapz,accum,diff)
!
	    Do k=1,nz/2+1
	      Write(309,'(7E15.7)') 0.0e0, kapz(k)/retau_avg_sa, yc(jp)*retau_avg_sa, &
	 accum(k)/utau_avg_sa**2, diff(k)/(utau_avg_sa*nu), diff(k)*kapz(k)/utau_avg_sa**2
	      Write(310,'(7E15.7)') 0.0e0, kapz(k), yc(jp), accum(k), diff(k), diff(k)*kapz(k)
	    End Do
	  End Do
!!
	Else If(ixz == 1) Then
	  Allocate(tmpfield(1:nx,jmin:jmax,1:nz))
	  Open(Unit=514,File='rms_snap.plt')
	  Write(514,'(A)') ' variables="y","rms","rmsl","rmss","y<sup>+</sup>"'
	  If(ilarge3d == 1) Then !write 3d velocity file for super scale structure
	    Open(Unit=515,File='u_snap_3d_large.vtk')
	    Write(515,'(A)') '# vtk DataFile Version 2.0'
	    Write(515,'(A)') 'Volume example'
	    Write(515,'(A)') 'ASCII'
	    Write(515,'(A)') 'DATASET RECTILINEAR_GRID'
	    Write(515,'(A,3I10)') 'DIMENSIONS',imax-imin+1,jmax-jmin+1,nz
	    Write(515,'(A,I10,A)') 'X_COORDINATES', imax-imin+1, ' float'
	    Write(515,'(10E15.7)') xc(imin:imax)
	    Write(515,'(A,I10,A)') 'Y_COORDINATES', jmax-jmin+1, ' float'
	    Write(515,'(10E15.7)') yc(jmin:jmax)
	    Write(515,'(A,I10,A)') 'Z_COORDINATES', nz, ' float'
	    Write(515,'(10E15.7)') zc(1:nz)
	    Write(515,'(A,I10)') 'POINT_DATA',(imax-imin+1)*(jmax-jmin+1)*nz
	    Write(515,'(A)') 'SCALARS ularge float 1'
	    Write(515,'(A)') 'LOOKUP_TABLE default'
!
	    Open(Unit=516,File='u_snap_3d_small.vtk')
	    Write(516,'(A)') '# vtk DataFile Version 2.0'
	    Write(516,'(A)') 'Volume example'
	    Write(516,'(A)') 'ASCII'
	    Write(516,'(A)') 'DATASET RECTILINEAR_GRID'
	    Write(516,'(A,3I10)') 'DIMENSIONS',imax-imin+1,jmax-jmin+1,nz
	    Write(516,'(A,I10,A)') 'X_COORDINATES', imax-imin+1, ' float'
	    Write(516,'(10E15.7)') xc(imin:imax)
	    Write(516,'(A,I10,A)') 'Y_COORDINATES', jmax-jmin+1, ' float'
	    Write(516,'(10E15.7)') yc(jmin:jmax)
	    Write(516,'(A,I10,A)') 'Z_COORDINATES', nz, ' float'
	    Write(516,'(10E15.7)') zc(1:nz)
	    Write(516,'(A,I10)') 'POINT_DATA',(imax-imin+1)*(jmax-jmin+1)*nz
	    Write(516,'(A)') 'SCALARS usmall float 1'
	    Write(516,'(A)') 'LOOKUP_TABLE default'
	  End If
!
	  Allocate(accum2(1:nx/2+1,1:nz/2+1),diff2(1:nx/2+1,1:nz/2+1))
	  Do jp=jmin,jmax; Write(*,*) 'jp=',jp
            ux(1:nx,1:nz)=dummy_3d1(1:nx,jp,1:nz,m)
	    Call dfftw_plan_dft_r2c_2d(plan,nx,nz,ux,outfftxu,FFTW_ESTIMATE)
	    Call dfftw_execute(plan)
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxu(1:nx/2+1,1:nz))*outfftxu(1:nx/2+1,1:nz)/Real(nx*nz,mytype)/Real(nx*nz,mytype)
! calculate the toatal energy in xz plane
	    If(m == 1) Then 
	      rms=Sum((ux(1:nx,1:nz)-umn(jp))**2)/Real(nx*nz,mytype)
	    Else If(m == 2) Then
	      rms=Sum(ux(1:nx,1:nz)**2)/Real(nx*nz,mytype)
	    Else If(m == 3) Then
	      rms=Sum((ux(1:nx,1:nz)-wmn(jp))**2)/Real(nx*nz,mytype)
	    End If
	    rmss=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
		Sum(spec2d(nx/2+1,1:nz))
! super structure
	    Call specxz(spec2d,kapx,kapz,accum2,diff2)
	    rmsl=0.0e0
	    Do k=1,nz; If(k<nz/2+1) Then; kz=k; Else; kz=nz-k+2; End If
	      Do i=1,nx/2+1
	      If((idef_sup.Eq.1 .And. (kz.Le.(kz_large+1)) .Or. &
	      (idef_sup.Eq.2 .And. i.Le.(kx_large+1))) .Or. &
	      (idef_sup.Eq.3 .And. kz.Le.(kz_large+1) .And. i.Le.(kx_large+1) ) .Or. &
	      (idef_sup.Eq.4 .And. (kz.Le.(kz_large+1) .Or. i.Le.(kx_large+1)) ) ) Then !kz<=kz_large
	        rmsl=rmsl+diff2(i,k)*(2.0e0*pi/length)*(2.0e0*pi/width)
	      End If
	    End Do; End Do
	    Write(*,*) 'total energy in xz plane:', rmss,rmsl
	    Write(*,*) 'rms energy at each point:', Sqrt(rmss)
!
	    Write(*,*) 'snap total fluctuation:', rms
	    Write(*,*) 'snap rms fluctuation:', Sqrt(rms)/utau_avg_sa
	    Write(514,'(5F20.10)') yc(jp),Sqrt(rms),Sqrt(rmsl),&
	     Sqrt((rms-rmsl)),ycplus(jp)
!
	    If(ilarge3d == 1) Then ! reconstruct super scale structure field
	      outfftxu(:,:)=outfftxu(:,:)/Real(nx*nz,mytype)
!!	remove large scales
	    Do k=1,nz; If(k<nz/2+1) Then; kz=k; Else; kz=nz-k+2; End If
	      Do i=1,nx/2+1
	        If((idef_sup.Eq.1 .And. (kz.Le.(kz_large+1)) .Or. &
	        (idef_sup.Eq.2 .And. i.Le.(kx_large+1))) .Or. &
	        (idef_sup.Eq.3 .And. kz.Le.(kz_large+1) .And. i.Le.(kx_large+1) ) .Or. &
	        (idef_sup.Eq.4 .And. (kz.Le.(kz_large+1) .Or. i.Le.(kx_large+1)) ) ) Then !kz<=kz_large
	          outfftxu(i,k)=0.0e0
  	        End If
	      End Do; End Do
!	keep one fourier mode
!	      Do k=1,nz; Do i=1,nx/2+1
!	        If(i.Ne.3.Or.k.Ne.1) outfftxu(i,k)=0.0e0
!	      End Do; End Do
!
	      Call dfftw_plan_dft_c2r_2d(plan,nx,nz,outfftxu,ux,FFTW_ESTIMATE)
	      Call dfftw_execute(plan)
!
	      tmpfield(1:nx,jp,1:nz)=ux(1:nx,1:nz)
	    End If
	    If(jp == jspec) Then
 	      Write(309,*) 'zone i=',nx/2,', j=',nz/2,',f=point'
 	      Write(310,*) 'zone i=',nx/2,', j=',nz/2,',f=point'
!
	      DO k=2,nz/2+1; DO i=2,nx/2+1
	        Write(309,'(7E15.7)') kapx(i)/retau_avg_sa, kapz(k)/retau_avg_sa, yc(jp)*retau_avg_sa, &
	        accum2(i,k)/utau_avg_sa**2, diff2(i,k)/(nu*nu), diff2(i,k)*kapx(i)*kapz(k)/utau_avg_sa**2
	        Write(310,'(7E15.7)') kapx(i),kapz(k),yc(jp),accum2(i,k),diff2(i,k),diff2(i,k)*kapx(i)*kapz(k)
	      End Do; End Do; 
            End If
	  End Do; Close(514)
!      write filted structure field
	  If(ilarge3d==1) Then
	    Do k=1,nz; Do j=jmin,jmax; Do i=imin,imax
	      Write(515,'(10E15.7)') dummy_3d1(i,j,k,m)-umn(j)!-tmpfield(i,j,k)
	      Write(516,'(10E15.7)') tmpfield(i,j,k)
	    End Do; End Do; End Do
	    Deallocate(tmpfield)
	  End If
!	  
	  Deallocate(accum2,diff2)
	End If	
!
	Close(309); Close(310)
	If(ilarge3d == 1) Then; Close(515); Close(516); End If

	Deallocate(spec2d,outfftxu,ux,accum,diff,lamx,kapx,lamz,kapz)

	End Subroutine spectra2D
!*******************************************************************************
	Subroutine spectra_uv(dummy_3d1,umn,wmn)
!plot 2D spectra in xy plane
  	Include "fftw3.f"
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(in) :: dummy_3d1
	Real(mytype), Dimension(ny), Intent(in)  :: umn,wmn
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxu,outfftxv
	Real(mytype), Allocatable, Dimension(:,:) :: spec2d
	Real(mytype), Allocatable, Dimension(:,:) :: ux,vx
	Real(mytype), Allocatable, Dimension(:) :: lamx,lamz,kapx,kapz
	Real(mytype), Allocatable, Dimension(:,:) :: accum2,diff2
        Real(mytype) :: wl,wlx,wlz, wn,wnx,wnz, rms,rmsl,rmss
	Integer :: i,j,k,jp,jmin,jmax,ixy,ixz,iyz,m,kz
        Integer*8 :: plan

! averaged FFT in xy plane (ixy=1); xz plane (ixz=1) and yz plane (iyz=1)
! m is the variable for FFT (1 for u; 2 for v; 3 for w and 4 for p)
	ixy=0;ixz=1;iyz=0
! FFT range in y direction
	jmin=jspec;jmax=jspec !jspec
!
	Allocate(lamx(1:nx/2+1),kapx(1:nx/2+1), &
	lamz(1:nz/2+1),kapz(1:nz/2+1))
	lamx(:)=0.0e0; kapx(:)=0.0e0; lamz(:)=0.0e0; kapz(:)=0.0e0
!
	lamx(1)=0.0e0; kapx(1)=0.0e0
	Do i=2,nx/2+1; lamx(i)=length/real(i-1,mytype); kapx(i)=2.0e0*pi/lamx(i); End Do	
	lamz(1)=0.0e0; kapz(1)=0.0e0
	Do k=2,nz/2+1; lamz(k)=width/real(k-1,mytype);kapz(k)=2.0e0*pi/lamz(k); End Do

	Open(Unit=310,File='spectra_2d_uv.plt')
	Write(310,'(A)') ' variables="k<sub>x</sub>","k<sub>z</sub>","y",', &
	 '"E<sub>11</sub>", "k<sub>x</sub>k<sub>z</sub>E<sub>11</sub>"'
!
	Allocate(ux(1:nx,1:nz),outfftxu(1:nx/2+1,1:nz),outfftxv(1:nx/2+1,1:nz),&
	spec2d(1:nx/2+1,1:nz))
        spec2d(:,:)=0.0e0; ux(:,:)=0.0e0
!
	If(ixz == 1) Then
	  Open(Unit=514,File='rms_snap_uv.plt')
	  Write(514,'(A)') ' variables="y","rms","rmsl","rmss","y<sup>+</sup>"'
!
	  Allocate(accum2(1:nx/2+1,1:nz/2+1),diff2(1:nx/2+1,1:nz/2+1))
	  Do jp=jmin,jmax; Write(*,*) 'jp=',jp
          ux(1:nx,1:nz)=dummy_3d1(1:nx,jp,1:nz,1)
	  Call dfftw_plan_dft_r2c_2d(plan,nx,nz,ux,outfftxu,FFTW_ESTIMATE)
	  Call dfftw_execute(plan)
          ux(1:nx,1:nz)=dummy_3d1(1:nx,jp,1:nz,2)
	  Call dfftw_plan_dft_r2c_2d(plan,nx,nz,ux,outfftxv,FFTW_ESTIMATE)
	  Call dfftw_execute(plan)
!
	  spec2d(1:nx/2+1,1:nz)=-outfftxu(1:nx/2+1,1:nz)*Conjg(outfftxv(1:nx/2+1,1:nz)) &
		/Real(nx*nz,mytype)/Real(nx*nz,mytype)
! calculate the toatal energy in xz plane
	  rms=-Sum((dummy_3d1(1:nx,jp,1:nz,1)-umn(jp))*dummy_3d1(1:nx,jp,1:nz,2))
	  rmss=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
		Sum(spec2d(nx/2+1,1:nz))
! super structure
	  Call specxz(spec2d,kapx,kapz,accum2,diff2)
	  rmsl=0.0e0
	  Do k=1,nz; If(k<nz/2+1) Then; kz=k; Else; kz=nz-k+2; End If
	    Do i=1,nx/2+1
	    If((idef_sup.Eq.1 .And. (kz.Le.(kz_large+1)) .Or. &
	      (idef_sup.Eq.2 .And. i.Le.(kx_large+1))) .Or. &
	      (idef_sup.Eq.3 .And. kz.Le.(kz_large+1) .And. i.Le.(kx_large+1) ) .Or. &
	      (idef_sup.Eq.4 .And. (kz.Le.(kz_large+1) .Or. i.Le.(kx_large+1)) ) ) Then !kz<=kz_large
	      rmsl=rmsl+diff2(i,k)*(2.0e0*pi/length)*(2.0e0*pi/width)
	    End If
	  End Do; End Do
	  Write(*,*) 'uv energy at each point:', rmss, rmsl
!
	  Write(*,*) 'snap total uv:', rms
	  Write(*,*) 'snap uv:', rms/Real(nx*nz,mytype),rms/Real(nx*nz,mytype)/utau_avg_sa**2
	  Write(514,'(5F20.10)') ycplus(jp),rmss/utau_avg_sa**2,rmsl,&
	   rmss-rmsl,ycplus(jp)
!
	  If(jp == jspec) Then
 	  Write(310,*) 'zone i=',nx/2,', j=',nz-1,',f=point'
	  DO k=2,nz; DO i=2,nx/2+1
            wlx=length/Real(i-1); wlz=width/Real(k-1)
	    If(k.Ge.nz/2+1) wlz=width/Real(nz+1-k)
            wnx=2.0e0*pi/wlx; wnz=2.0e0*pi/wlz
	    Write(310,'(5E15.7)') Log10(wnx/retau_avg_sa),Log10(wnz/retau_avg_sa),0.0e0,spec2d(i,k),spec2d(i,k)*wnx*wnz/utau_avg_sa**2
	  END DO; END DO; End If
	  End Do; Close(514)
	  Deallocate(accum2,diff2)
	End If	
!
	Close(310)

	Deallocate(spec2d,outfftxu,outfftxv,ux)
	Deallocate(lamx,kapx,lamz,kapz)

	End Subroutine spectra_uv
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
!
	diff(1:nx/2+1,1:nz/2+1)=euu(1:nx/2+1,1:nz/2+1)/dkapx/dkapz
	Deallocate(euu)
!
	End Subroutine specxz
!*******************************************************************************!
	Subroutine spectra2D2(dummy_3d1,umn,wmn)
!plot 2D spectra in xy plane
  	Include "fftw3.f"
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(in) :: dummy_3d1
	Real(mytype), Dimension(ny), Intent(in)  :: umn,wmn
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxu
	Real(mytype), Allocatable, Dimension(:,:) :: spec2d
	Real(mytype), Allocatable, Dimension(:,:) :: ux
        Real(mytype) :: wl,wlx,wlz, wn,wnx,wnz
	Integer :: i,j,k,jp,jmin,jmax,ixy,ixz,iyz,m
        Integer*8 :: plan

! averaged FFT in xy plane (ixy=1); xz plane (ixz=1) and yz plane (iyz=1)
! m is the variable for FFT (1 for u; 2 for v; 3 for w and 4 for p)
	ixy=1;ixz=0;iyz=0;m=3
! FFT range in y direction
	jmin=jspec;jmax=jspec

	Open(Unit=310,File='spectra_2d.plt')
	WRITE(310,'(A)') ' variables="k<sub>x</sub>","k<sub>z</sub>","y",', &
	 '"E<sub>11</sub>", "k<sub>x</sub>k<sub>z</sub>E<sub>11</sub>"'
!
	Allocate(ux(1:nz,1:nx),outfftxu(1:nz/2+1,1:nx),spec2d(1:nz/2+1,1:nx))
        spec2d(:,:)=0.0e0; ux(:,:)=0.0e0
!
	If(ixy == 1) Then ! xy view, averaged in z direction
 	  Write(310,*) 'zone i=',nx/2,', j=',jmax-jmin+1,',f=point'
          Do jp=jmin,jmax !y
	    Write(*,*) 'y+=', yplus(jp)
	    Do k=1,nz; Do i=1,nx
              ux(k,i)=dummy_3d1(i,jp,k,m)
	    End Do; End Do
	    Call dfftw_plan_dft_r2c_2d(plan,nz,nx,ux,outfftxu,FFTW_ESTIMATE)
	    Call dfftw_execute(plan)
	    spec2d(1:nz/2+1,1:nx)=cdabs(outfftxu(1:nz/2+1,1:nx))**2/Real(nx*nz,mytype)
	    spec2d(1,1:nx)=spec2d(1,1:nx)/2.0e0
!
	    Do i=2,nx/2+1
	      wlx=length/real(i-1); wnx=2.0e0*pi/wlx
	      Write(310,'(5E15.7)') wnx, 0.0e0, yc(jp), (Sum(spec2d(1:nz/2+1,i))+Sum(spec2d(1:nz/2+1,nx+2-i)))/2.0e0/Real(nz/2+1,mytype),&
		 wnx*Sum(spec2d(1:nz/2+1,i))/Real(nz/2+1,mytype)
	    End Do
	  End Do
!!
	Else If(iyz == 1) Then
 	  Write(310,*) 'zone i=',nz/2,', j=',jmax-jmin+1,',f=point'
          Do jp=jmin,jmax !y
	    Write(*,*) 'y+=', yplus(jp)
	    Do k=1,nz; Do i=1,nx
              ux(k,i)=dummy_3d1(i,jp,k,m)
	    End Do; End Do
	    Call dfftw_plan_dft_r2c_2d(plan,nz,nx,ux,outfftxu,FFTW_ESTIMATE)
	    Call dfftw_execute(plan)
	    spec2d(1:nz/2+1,:)=cdabs(outfftxu(1:nz/2+1,:))**2/Real(nx*nz,mytype)
!
	    Do k=2,nz/2+1
	      wlz=width/real(k-1); wnz=2.0e0*pi/wlz
	      Write(310,'(5E15.7)') 0.0e0, wnz, yc(jp), &
	          Sum(spec2d(k,1:nx))/Real(nx,mytype),&
		  wnz*Sum(spec2d(k,1:nx))/Real(nx,mytype)
	    End Do
	  End Do
	Else If(ixz == 1) Then
	  Do k=1,nz; Do i=1,nx; ux(k,i)=dummy_3d1(i,jspec,k,m); End Do; End Do
	  Call dfftw_plan_dft_r2c_2d(plan,nz,nx,ux,outfftxu,FFTW_ESTIMATE)
	  Call dfftw_execute(plan)
	  Do i=1,nx; Do k=1,nz/2+1
	    spec2d(k,i)=cdabs(outfftxu(k,i))**2/Real(nx*nz,mytype)
	  End Do; End Do
! calculate the toatal energy in xz plane
	  Write(*,*) 'mean=',spec2d(1,1)
	  spec2d(1,1)=0.0e0; spec2d(1,1:nx)=spec2d(1,1:nx)/2.0e0
	  Write(*,*) 'total energy in xz plane:', Sum(spec2d(1:nz/2+1,1:nx))
	  Write(*,*) 'rms energy at each point:', Sqrt(Sum(spec2d(1:nz/2,1:nx))/Real(nx*nz,mytype)*2.0e0)
!
	  Write(*,*) 'snap total fluctuation:', Sum((ux(1:nz,1:nx)-umn(jspec))**2)*0.5e0
	  Write(*,*) 'snap rms fluctuation:', Sqrt(Sum((ux(1:nz,1:nx)-umn(jspec))**2)/Real(nx*nz,mytype))
!
 	  Write(310,*) 'zone i=',nz/2,', j=',nx-1,',f=point'
	  DO i=2,nx; DO k=2,nz/2+1
            wlx=length/Real(i-1); wlz=width/Real(k-1)
	    If(i.Ge.nx/2+1) wlx=length/Real(nx+1-i)
            wnx=2.0e0*pi/wlx; wnz=2.0e0*pi/wlz
	    Write(310,'(5E15.7)') Log10(wnx/retau_avg_sa),Log10(wnz/retau_avg_sa),0.0e0,spec2d(k,i),spec2d(k,i)*wnx*wnz/utau_avg_sa**2
	  END DO;END DO
	End If	

	Close(310)

	Deallocate(spec2d,outfftxu,ux)

	End Subroutine spectra2D2
!***************************************************************************
	SUBROUTINE pdf
!	This subroutine calculates the 2D pdf (probability density function)
! 	The first example is for u'v'
	Real(mytype), Dimension(:,:,:), Allocatable :: uvw,uvw2
!
	Integer :: i,j,k, i1,j1,k1, jp, np, cnt, m, nsp, ibt
	Real(mytype), Allocatable, Dimension(:,:) :: dumu,dumv
	Real(mytype), Allocatable, Dimension(:,:) :: cdf, pdf1
	Real(mytype) :: umin,umax,vmin,vmax, du,dv, ulim, vlim, tmp
        Real(mytype), dimension(8) :: dum
 	Integer, Dimension(3) :: intdum
	Character :: chdum*8
	Character*60, Dimension(:), Allocatable :: file3dw
!
	Allocate(uvw(-1:nx+2,-1:ny+2,-1:nz+2),uvw2(-1:nx+2,-1:ny+2,-1:nz+2))
	Allocate(dumu(1:nx,1:nz),dumv(1:nx,1:nz))
	jp=jspec; nsp=1
	Allocate(file3dw(1:nsp))
	file3dw=(/'../restrt-start0-stand-run8-out.rst'/)
	WRITE(*,*) 'Get the pdf for uv at y=',y(jp),'yp=',yplus(jp)

! write data to file
	OPEN(UNIT=311,FILE='pdf-uv.plt')
	WRITE(311,'(320A)') ' variables="u","v","pdf","cdf"'
	OPEN(UNIT=312,FILE='quadrant-scatter.plt')
	WRITE(312,'(320A)') ' variables="u","v"'

	Do m=1,nsp
	  Write(chdum,'(I5)') m
	  Write(*,*) ''
          Open(Unit=28, file=file3dw(m), Status='OLD', Form='UNFormatTED',access='stream')
	  Read(28) intdum(:),dum(:)
	  Write(*,*) intdum(:),dum(:)
	  Read(28) uvw(:,:,:)
	  Read(28) uvw2(:,:,:)
	  Close(28)

	Do ibt=1,2 ! top and bottom average
!
	DO k=1,nz; DO i=1,nx
	  dumu(i,k)=0.5e0*(uvw(i,jp,k)+uvw(i+1,jp,k))
	  dumv(i,k)=0.5e0*(uvw2(i,jp,k)+uvw2(i,jp+1,k))
 	END DO; END DO
	If(ibt==1) Then
	  DO k=1,nz; DO i=1,nx
	    dumu(i,k)=dumu(i,k)-flucs(jp,1,1)
	    dumv(i,k)=dumv(i,k)-flucs(jp,2,1)
	  END DO; END DO
	Else
	  DO k=1,nz; DO i=1,nx
	    dumu(i,k)=dumu(i,k)-flucs(nyt-jp,1,1)
	    dumv(i,k)=dumv(i,k)-flucs(nyt-jp,2,1)
	  END DO; END DO
	End If
! non-dimensionalisation
	dumu(1:nx,1:nz)=dumu(1:nx,1:nz)/utau_avg_sa
	dumv(1:nx,1:nz)=dumv(1:nx,1:nz)/utau_avg_sa
	np=50
	Allocate(cdf(1:np,1:np),pdf1(1:np,1:np))
	umin=dumu(1,1); umax=dumu(1,1)
	DO k=1,nz; DO i=1,nx
	  IF(umin.GT.dumu(i,k)) umin=dumu(i,k)
	  IF(umax.LT.dumu(i,k)) umax=dumu(i,k)
	END DO; END DO
	du=(umax-umin)/REAL(np-1)
	vmin=dumv(1,1); vmax=dumv(1,1)
	DO k=1,nz; DO i=1,nx
	  IF(vmin.GT.dumv(i,k)) vmin=dumv(i,k)
	  IF(vmax.LT.dumv(i,k)) vmax=dumv(i,k)
	END DO; END DO
	dv=(vmax-vmin)/REAL(np-1)
	cdf=0.0e0; pdf1=0.0e0
	DO k1=1,np; DO i1=1,np
	  cnt=0; 
	  ulim=REAL(i1-1)*du+umin; vlim=REAL(k1-1)*dv+vmin
	  DO k=1,nz; DO i=1,nx
	    IF(dumu(i,k).LT.ulim.AND.dumv(i,k).LT.vlim) cnt=cnt+1
	  END DO; END DO
	  cdf(i1,k1)=REAL(cnt)/REAL(nx*nz)
	END DO; END DO
! get the pdf from cdf
	DO k1=2,np-1; DO i1=2,np-1
	  pdf1(i1,k1)=(cdf(i1+1,k1+1)-cdf(i1-1,k1+1)-cdf(i1+1,k1-1)+cdf(i1-1,k1-1))/4.0e0/du/dv
!	  pdf1(i1,k1)=(cdf(i1,k1+1)-cdf(i1-1,k1+1)-cdf(i1,k1)+cdf(i1-1,k1))/4.0e0/du/dv
	END DO; END DO

 	WRITE(311,*) 'zone t="',m,'", i=',np-1,', j=',np-1,',f=point'
	DO k=2,np; DO i=2,np
	  WRITE(311,'(4F20.10)') du*REAL(i-1)+umin, dv*REAL(k-1)+vmin, pdf1(i,k), cdf(i,k)
	END DO; END DO
	
 	WRITE(312,*) 'zone t="',m,'", i=',nx,', j=',nz,',f=point'
	DO k=1,nz; DO i=1,nx
	  WRITE(312,'(2F20.10)') dumu(i,k), dumv(i,k)
	END DO; END DO
!
	Deallocate(cdf,pdf1)
!
	End Do ! ibt
	End Do !m=1,nsp

	CLOSE(311); CLOSE(312)	

	DEALLOCATE(dumu,dumv,uvw,uvw2,file3dw)
	END SUBROUTINE pdf
!*******************************************************************************
!***************************************************************************
	SUBROUTINE pdf_1d(u,irangefix,np,pdfx,pdf)
!	This subroutine calculates the 1D pdf (probability density function)
	Real(mytype), Dimension(1:nx,1:nz), Intent(in) :: u
	Integer, Intent(in) :: np,irangefix
	Real(mytype), Dimension(1:np), Intent(Out)  :: pdf,pdfx
!
	Integer :: i,j,k, i1,j1,k1, jp, cnt
	Real(mytype), Dimension(:,:), Allocatable :: dumu
	Real(mytype), Allocatable, Dimension(:) :: cdf
	Real(mytype) :: umin,umax,du,ulim, tmp
!
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
	If(irangefix==1) umin=-5.0e0; umax=5.0e0
	du=(umax-umin)/REAL(np-1)
	Do i=1,np
	  pdfx(i)=du*REAL(i-1)+umin
	End Do
	cdf=0.0e0; pdf=0.0e0
	DO i1=1,np
	  cnt=0; 
	  ulim=REAL(i1-1)*du+umin
	  DO k=1,nz; DO i=1,nx
	    IF(dumu(i,k).LT.ulim) cnt=cnt+1
	  END DO; End Do
	  cdf(i1)=REAL(cnt)/REAL(nx*nz)
	END DO
! get the pdf from cdf
	DO i1=2,np-1
	  pdf(i1)=(cdf(i1+1)-cdf(i1-1))/2.0e0/du
	END DO

! write data to file
	OPEN(UNIT=311,FILE='pdf-u.plt')
	WRITE(311,'(320A)') ' variables="u","pdf","cdf"'
 	WRITE(311,*) 'zone i=',np-1,',f=point'
	DO i=2,np
	  WRITE(311,'(4F20.10)') pdfx(i), pdf(i), cdf(i)
	END DO
	CLOSE(311)
!
	DEALLOCATE(dumu,cdf)
	END SUBROUTINE pdf_1d
!***************************************************************************
!***************************************************************************
	SUBROUTINE pdf_2d(u,v,irangefix,np,pdfx,pdfy,pdf)
!	This subroutine calculates the 2D pdf (probability density function)
! 	The first example is for u'v'
	Real(mytype), Dimension(1:nx,1:nz), Intent(in) :: u,v
	Integer, Intent(in) :: np,irangefix
	Real(mytype), Dimension(1:np,1:np), Intent(Out)  :: pdf
	Real(mytype), Dimension(1:np), Intent(Out)  :: pdfx,pdfy
!
	Integer :: i,j,k, i1,j1,k1, jp, cnt
	Real(mytype), Allocatable, Dimension(:,:) :: dumu,dumv
	Real(mytype), Allocatable, Dimension(:,:) :: cdf
	Real(mytype) :: umin,umax,vmin,vmax, du,dv, ulim, vlim, tmp
!
	Allocate(dumu(1:nx,1:nz),dumv(1:nx,1:nz))
	jp=jspec
	WRITE(*,*) 'Get the pdf for uv at y=',y(jp),'yp=',yplus(jp)
	DO k=1,nz; DO i=1,nx
	  dumu(i,k)=u(i,k)
	  dumv(i,k)=v(i,k)
	END DO; END DO
!
	Allocate(cdf(1:np,1:np))
	umin=dumu(1,1); umax=dumu(1,1)
	DO k=1,nz; DO i=1,nx
	  IF(umin.GT.dumu(i,k)) umin=dumu(i,k)
	  IF(umax.LT.dumu(i,k)) umax=dumu(i,k)
	END DO; END DO
	vmin=dumv(1,1); vmax=dumv(1,1)
	DO k=1,nz; DO i=1,nx
	  IF(vmin.GT.dumv(i,k)) vmin=dumv(i,k)
	  IF(vmax.LT.dumv(i,k)) vmax=dumv(i,k)
	END DO; END DO
!fixed u range
	If(irangefix==1) Then
	  umin=-5.0e0; umax=5.0e0
	  vmin=-1.50e0; vmax=1.50e0
	End If
	du=(umax-umin)/REAL(np-1)
	dv=(vmax-vmin)/REAL(np-1)
!
	Do i=1,np
	  pdfx(i)=du*REAL(i-1)+umin
	  pdfy(i)=dv*REAL(i-1)+vmin
	End Do
!
	cdf=0.0e0; pdf=0.0e0
	DO k1=1,np; DO i1=1,np
	  cnt=0; 
	  ulim=REAL(i1-1)*du+umin; vlim=REAL(k1-1)*dv+vmin
	  DO k=1,nz; DO i=1,nx
	    IF(dumu(i,k).LT.ulim.AND.dumv(i,k).LT.vlim) cnt=cnt+1
	  END DO; END DO
	  cdf(i1,k1)=REAL(cnt)/REAL(nx*nz)
	END DO; END DO
! get the pdf from cdf
	DO k1=2,np-1; DO i1=2,np-1
	  pdf(i1,k1)=(cdf(i1+1,k1+1)-cdf(i1-1,k1+1)-cdf(i1+1,k1-1)+cdf(i1-1,k1-1))/4.0e0/du/dv
	END DO; END DO

! write data to file
	OPEN(UNIT=311,FILE='pdf-uv.plt')
	WRITE(311,'(320A)') ' variables="u","v","pdf","cdf"'
 	WRITE(311,*) 'zone i=',np-1,', j=',np-1,',f=point'
	DO k=2,np; DO i=2,np
	  WRITE(311,'(4F20.10)') du*REAL(i-1)+umin, dv*REAL(k-1)+vmin, pdf(i,k), cdf(i,k)
	END DO; END DO
	CLOSE(311)

	DEALLOCATE(cdf,dumu,dumv)
	END SUBROUTINE pdf_2d
!*******************************************************************************
!***************************************************************************
	SUBROUTINE pdf_1d_super(u,usuper,irangefix,np,pdfx,pdf)
!	This subroutine calculates the 1D pdf (probability density function)
	Real(mytype), Dimension(1:nx,1:nz), Intent(in) :: u,usuper
	Integer, Intent(in) :: np,irangefix
	Real(mytype), Dimension(1:np), Intent(Out)  :: pdf,pdfx
!
	Integer :: i,j,k, i1,j1,k1, jp, cnt,ipos,cntp,cntn
	Real(mytype), Dimension(:,:), Allocatable :: dumu
	Real(mytype), Allocatable, Dimension(:) :: cdf
	Real(mytype) :: umin,umax,du,ulim, tmp
!
	ipos=0 !condition with positive (1) or negative (0) super scale structures or all structures (>=2)
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
	OPEN(UNIT=311,FILE='pdf-u-super.plt')
	WRITE(311,'(320A)') ' variables="u","pdf","cdf"'
 	WRITE(311,*) 'zone i=',np-1,',f=point'
	DO i=2,np
	  WRITE(311,'(4F20.10)') pdfx(i), pdf(i), cdf(i)
	END DO
	CLOSE(311)
!
	DEALLOCATE(dumu,cdf)
	END SUBROUTINE pdf_1d_super
!***************************************************************************
!***************************************************************************
	SUBROUTINE pdf_2d_super(u,v,usuper,irangefix,np,pdfx,pdfy,pdf)
!	This subroutine calculates the 2D pdf (probability density function)
! 	The first example is for u'v'
	Real(mytype), Dimension(1:nx,1:nz), Intent(in) :: u,v,usuper
	Integer, Intent(in) :: np,irangefix
	Real(mytype), Dimension(1:np,1:np), Intent(Out)  :: pdf
	Real(mytype), Dimension(1:np), Intent(Out)  :: pdfx,pdfy
!
	Integer :: i,j,k, i1,j1,k1, jp, cnt,ipos,cntp,cntn
	Real(mytype), Allocatable, Dimension(:,:) :: dumu,dumv
	Real(mytype), Allocatable, Dimension(:,:) :: cdf
	Real(mytype) :: umin,umax,vmin,vmax, du,dv, ulim, vlim, tmp
!
	ipos=1 !condition with positive (1) or negative (0) super scale structures or all structures (>=2)
	cntp=0;cntn=0
	Allocate(dumu(1:nx,1:nz),dumv(1:nx,1:nz))
	jp=jspec
	WRITE(*,*) 'Get the pdf for uv at y=',y(jp),'yp=',yplus(jp)
	DO k=1,nz; DO i=1,nx
	  dumu(i,k)=u(i,k)
	  dumv(i,k)=v(i,k)
	END DO; END DO
!
	Allocate(cdf(1:np,1:np))
	umin=dumu(1,1); umax=dumu(1,1)
	DO k=1,nz; DO i=1,nx
	  IF(umin.GT.dumu(i,k)) umin=dumu(i,k)
	  IF(umax.LT.dumu(i,k)) umax=dumu(i,k)
	END DO; END DO
	vmin=dumv(1,1); vmax=dumv(1,1)
	DO k=1,nz; DO i=1,nx
	  IF(vmin.GT.dumv(i,k)) vmin=dumv(i,k)
	  IF(vmax.LT.dumv(i,k)) vmax=dumv(i,k)
	END DO; END DO
!fixed u range
	If(irangefix==1) Then
	  umin=-5.0e0; umax=5.0e0
	  vmin=-1.50e0; vmax=1.50e0
	End If
	DO k=1,nz; DO i=1,nx
	  If(ipos==1) Then
	    If(usuper(i,k)<0.0e0) Then
	      dumu(i,k)=umin-1.0e0; dumv(i,k)=vmin-1.0e0
	    End If
	  Else If(ipos==0) Then
	    If(usuper(i,k)>0.0e0) Then
	      dumu(i,k)=umax+1.0e0; dumv(i,k)=vmax+1.0e0
	    End If
	  End If
	  If(usuper(i,k)>0.0e0) cntp=cntp+1
	END DO; END DO
	cntn=nx*nz-cntp
	du=(umax-umin)/REAL(np-1)
	dv=(vmax-vmin)/REAL(np-1)
!
	Do i=1,np
	  pdfx(i)=du*REAL(i-1)+umin
	  pdfy(i)=dv*REAL(i-1)+vmin
	End Do
!
	cdf=0.0e0; pdf=0.0e0
	DO k1=1,np; DO i1=1,np
	  cnt=0; 
	  ulim=REAL(i1-1)*du+umin; vlim=REAL(k1-1)*dv+vmin
	  DO k=1,nz; DO i=1,nx
	    IF(dumu(i,k).LT.ulim.And.dumu(i,k).GE.umin.AND.dumv(i,k).LT.vlim.AND.&
		dumv(i,k).GE.vmin) cnt=cnt+1
	  END DO; END DO
	  If(ipos==1) Then
	    cdf(i1,k1)=REAL(cnt)/REAL(cntp)
	  Else If(ipos==0) Then
	    cdf(i1,k1)=REAL(cnt)/REAL(cntn)
	  Else
	    cdf(i1,k1)=REAL(cnt)/REAL(nx*nz)
	  End If
	END DO; END DO
! get the pdf from cdf
	DO k1=2,np-1; DO i1=2,np-1
	  pdf(i1,k1)=(cdf(i1+1,k1+1)-cdf(i1-1,k1+1)-cdf(i1+1,k1-1)+cdf(i1-1,k1-1))/4.0e0/du/dv
	END DO; END DO

! write data to file
	OPEN(UNIT=311,FILE='pdf-uv.plt')
	WRITE(311,'(320A)') ' variables="u","v","pdf","cdf"'
 	WRITE(311,*) 'zone i=',np-1,', j=',np-1,',f=point'
	DO k=2,np; DO i=2,np
	  WRITE(311,'(4F20.10)') du*REAL(i-1)+umin, dv*REAL(k-1)+vmin, pdf(i,k), cdf(i,k)
	END DO; END DO
	CLOSE(311)

	DEALLOCATE(cdf,dumu,dumv)
	END SUBROUTINE pdf_2d_super
!*******************************************************************************
!***************************************************************************
	SUBROUTINE write_1d_slice(dummy_3d)
!	This subroutine write 2d slice
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(In):: dummy_3d
	Real(mytype), Dimension(:,:), Allocatable :: dummy_1d
	Real(mytype), Dimension(:), Allocatable :: xyz1
	Character :: file1d*40
	Integer :: i,j,k,l,dim1,ipsig,jps
	Real(mytype) :: dumm
!
	ipsig=0; jps=6
	If(ipsig.Eq.1) Then
	  Write(file1d,'(I)') jps
	  file1d='psigle_'//Trim(Adjustl(file1d))//'.plt'
	  Open(506, File=file1d, Status='NEW')
	  Write(506,'(A)') 'variables= "time", "u", "v", "w", "p"'
	  Write(506,'(A,I5,A)') 'zone t="pos',jps,'"'
	End If
!read OUTPUT_1D_SNAP.dat
	Write(*,*) '> create 1d snapshot'
	If(lp_snap_x.Eq.1) Then
	  dim1=nx
	  Allocate(xyz1(0:dim1+1))
	  xyz1(0:dim1+1)=xc(0:nx+1)
	Else If(lp_snap_y.Eq.1) Then
	  dim1=ny
	  Allocate(xyz1(0:dim1+1))
	  xyz1(0:dim1+1)=yc(0:ny+1)
	Else If(lp_snap_z.Eq.1) Then
	  dim1=nz
	  Allocate(xyz1(0:dim1+1))
	  xyz1(0:dim1+1)=zc(0:nz+1)
	End If
	Allocate(dummy_1d(-1:dim1+2,1:4))

	If(lp_snap_1d.Eq.1) Then
	  file1d='OUTPUT_1D_SNAP.dat'
	  CALL change_dir(file1d,nruns(simit))
	  Open(505, File=file1d, Status='OLD', Form='UNFormatTED', access='stream')
 	  If(iconstruc.Ne.0) Write(502,*) 'zone t="2d" i=',dim1+2,' j=',nsmpl_snap_1d,'f=point'
	  Do l=1,nsmpl_snap_1d
            Read(505) dumm
	    Write(*,*) dumm   
	    Read(505) dummy_1d(:,1) !u
	    Read(505) dummy_1d(:,2) !v
	    Read(505) dummy_1d(:,3) !w
	    Read(505) dummy_1d(:,4) !p
	    If(iconstruc.Eq.0) Then
 	      WRITE(502,*) 'zone t="',dumm,'" i=',dim1+2,'f=point'
	      Do i=0,dim1+1
	        Write(502,'(10f20.10)') xyz1(i),dummy_1d(i,1:4)
	      End Do	
	    Else 
	    !construct 2d slice with time dimension 
	      Do i=0,dim1+1
	        Write(502,'(10f20.10)') xyz1(i),dumm,dummy_1d(i,1:4)
	      End Do
	    End If
	    If(ipsig.Eq.1) Write(506,'(5F20.10)') dumm,dummy_1d(jps,1:4)
	  End Do
	  Close(505)
	  If(ipsig.Eq.1) Close(506)
	End If
!
	If(iplot1d.Eq.1) Then
	  If(lp_snap_x.Eq.1) Then
	    k=nz/2;j=nx/2
	    dim1=nx
	    Do i=0,dim1+1
	      dummy_1d(i,1)=0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	      dummy_1d(i,2)=0.5*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2))
	      dummy_1d(i,3)=0.5*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	      dummy_1d(i,4)=dummy_3d(i,j,k,4)
	    End Do
	  Else If(lp_snap_y.Eq.1) Then
	    i=nx/2;k=nz/2
	    dim1=ny
	    Do j=0,dim1+1
	      dummy_1d(j,1)=0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	      dummy_1d(j,2)=0.5*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2))
	      dummy_1d(j,3)=0.5*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	      dummy_1d(j,4)=dummy_3d(i,j,k,4)
	    End Do
	  Else If(lp_snap_z.Eq.1) Then
	    i=nx/2;j=ny/2
	    dim1=nz
	    Do k=0,dim1+1
	      dummy_1d(k,1)=0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	      dummy_1d(k,2)=0.5*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2))
	      dummy_1d(k,3)=0.5*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	      dummy_1d(k,4)=dummy_3d(i,j,k,4)
	    End Do
	  End If
 	  WRITE(502,*) 'zone t="slice" i=',dim1+2,'f=point'
	  Do i=0,dim1+1
	    Write(502,'(10f20.10)') xyz1(i),dummy_1d(i,1:4)
	  End Do
	End If

	Deallocate(dummy_1d,xyz1)
	Close(502)

	END SUBROUTINE write_1d_slice
!*******************************************************************************
!*******************************************************************************
	SUBROUTINE write_2d_slice(dummy_3d,umn,wmn)
!	This subroutine write 2d slice
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(In):: dummy_3d
	Real(mytype), Dimension(:), Intent(in) ::umn,wmn
	Real(mytype), Dimension(:,:,:), Allocatable :: dummy_2d
	Real(mytype), Dimension(:), Allocatable :: xyz1,xyz2
	Character :: file2d*40
	Integer :: i,j,k,l,dim1,dim2,dim1vtk,dim2vtk,dim3vtk,irmean
	Real(mytype) :: dumm

	irmean=0 !if remove mean velocity
!read OUTPUT_2D_SNAP.dat
	Write(*,*) '> create 2d snapshot'
	If(lp_snap_xy.Eq.1) Then
	  dim1=nx;dim2=ny;dim1vtk=nx;dim2vtk=ny;dim3vtk=1
	  Allocate(xyz1(0:dim1+1),xyz2(0:dim2+1))
	  xyz1(0:dim1+1)=xc(0:nx+1);xyz2(0:dim2+1)=yc(0:ny+1)
	Else If(lp_snap_yz.Eq.1) Then
	  dim1=ny;dim2=nz;dim1vtk=1;dim2vtk=ny;dim3vtk=nz
	  Allocate(xyz1(0:dim1+1),xyz2(0:dim2+1))
	  xyz1(0:dim1+1)=yc(0:ny+1);xyz2(0:dim2+1)=zc(0:nz+1)
	Else If(lp_snap_xz.Eq.1) Then
	  dim1=nx;dim2=nz;dim1vtk=nx;dim2vtk=1;dim3vtk=nz
	  Allocate(xyz1(0:dim1+1),xyz2(0:dim2+1))
	  xyz1(0:dim1+1)=xc(0:nx+1);xyz2(0:dim2+1)=zc(0:nz+1)
	End If
	Allocate(dummy_2d(-1:dim1+2,-1:dim2+2,1:4))

	If(lp_snap_2d.Eq.1) Then
	  file2d='OUTPUT_2D_SNAP.dat'
	  CALL change_dir(file2d,nruns(simit))
	  Open(504, File=file2d, Status='OLD', Form='UNFormatTED', access='stream')
	  Do l=1,nsmpl_snap_2d
            Read(504) dumm
	    Write(*,*) dumm   
	    Read(504) dummy_2d(:,:,1) !u
	    Read(504) dummy_2d(:,:,2) !v
	    Read(504) dummy_2d(:,:,3) !w
	    Read(504) dummy_2d(:,:,4) !p
!
	    Write(file2d,'(I)') l
	    If(itec==1) Then
	    file2d='plane_'//Trim(Adjustl(file2d))//'.plt'
	    Else; file2d='plane_'//Trim(Adjustl(file2d))//'.vtk'; End If
	    If(l.Ge.2) Then
	      Close(501); Open(501,File=file2d)
	      Write(501,'(A)') '# vtk DataFile Version 2.0'
	      Write(501,'(A)') 'Volume example'
	      Write(501,'(A)') 'ASCII'
	      Write(501,'(A)') 'DATASET RECTILINEAR_GRID'
	    End If
!
	    If(iconstruc.Eq.0) Then
	      If(itec==1) Then
 	        WRITE(501,*) 'zone t="',dumm,'" i=',dim1+2,', j=',dim2+2,'f=point'
	        Do j=0,dim2+1; Do i=0,dim1+1
	          Write(501,'(10f20.10)') xyz1(i),xyz2(j),dummy_2d(i,j,1:4)
	        End Do; End Do
	      Else
!	    paraview format
	        Write(501,'(A,3I10)') 'DIMENSIONS',dim1vtk,dim2vtk,dim3vtk
		Write(501,'(A,I10,A)') 'X_COORDINATES', dim1vtk, ' float'
		Do i=1,dim1vtk; Write(501,'(E15.7)') xc(i); End Do
		Write(501,'(A,I10,A)') 'Y_COORDINATES', dim2vtk, ' float'
		Do j=1,dim2vtk; Write(501,'(E15.7)') yc(j); End Do
		Write(501,'(A,I10,A)') 'Z_COORDINATES', dim3vtk, ' float'
		Do k=1,dim3vtk; Write(501,'(E15.7)') zc(k); End Do
		Write(501,'(A,I10)') 'POINT_DATA',dim1vtk*dim2vtk*dim3vtk
		Write(501,'(A)') 'SCALARS u float 1'
		Write(501,'(A)') 'LOOKUP_TABLE default'
		Write(501,'(10E20.10)') ((dummy_2d(i,j,1),i=1,dim1),j=1,dim2)
		Write(501,'(A)') 'SCALARS v float 1'
		Write(501,'(A)') 'LOOKUP_TABLE default'
		Write(501,'(10E20.10)') ((dummy_2d(i,j,2),i=1,dim1),j=1,dim2)
		Write(501,'(A)') 'SCALARS w float 1'
		Write(501,'(A)') 'LOOKUP_TABLE default'
		Write(501,'(10E20.10)') ((dummy_2d(i,j,3),i=1,dim1),j=1,dim2)
		Write(501,'(A)') 'SCALARS p float 1'
		Write(501,'(A)') 'LOOKUP_TABLE default'
		Write(501,'(10E20.10)') ((dummy_2d(i,j,4),i=1,dim1),j=1,dim2)
	      End If
	    Else
	      If(itec==1) Then
	        Do i=1,dim1
 	        WRITE(501,*) 'zone t="',dumm,x(i),'" i=',dim2+2,'f=point'
	          Do j=0,dim2+1
	            Write(501,'(10f20.10)') xyz2(j),dummy_2d(i,j,1:4)
	          End Do
	        End Do
	      Else
		Write(*,*) 'ERROR***not vtk file!'
	      End If
	    End If		  
	  End Do
	  Close(504)
	End If
!
	If(iplot2d.Eq.1 .And. lp_snap_2d.Ne.1) Then
	  If(lp_snap_xy.Eq.1) Then
	    k=islicexyz
	    dim1=nx;dim2=ny
	    Do j=0,ny+1; Do i=0,nx+1
	      dummy_2d(i,j,1)=0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	      dummy_2d(i,j,2)=0.5*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2))
	      dummy_2d(i,j,3)=0.5*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	      dummy_2d(i,j,4)=dummy_3d(i,j,k,4)
	    End Do; End Do
	    If(irmean==1) Then
	      Do j=0,ny+1; dummy_2d(:,j,1)=dummy_2d(:,j,1)-umn(j); dummy_2d(:,j,3)=dummy_2d(:,j,3)-umn(j); End Do
	    End If
	  Else If(lp_snap_yz.Eq.1) Then
	    i=islicexyz
	    dim1=ny;dim2=nz
	    Do k=0,nz+1; Do j=0,ny+1
	      dummy_2d(j,k,1)=0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	      dummy_2d(j,k,2)=0.5*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2))
	      dummy_2d(j,k,3)=0.5*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	      dummy_2d(j,k,4)=dummy_3d(i,j,k,4)
	    End Do; End Do
	    If(irmean==1) Then
	      Do j=0,ny+1; dummy_2d(j,:,1)=dummy_2d(j,:,1)-umn(j); dummy_2d(j,:,3)=dummy_2d(j,:,3)-umn(j); End Do
	    End If
	  Else If(lp_snap_xz.Eq.1) Then
	    j=islicexyz
	    dim1=nx;dim2=nz
	    Do k=0,nz+1; Do i=0,nx+1
	      dummy_2d(i,k,1)=0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1))
	      dummy_2d(i,k,2)=0.5*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2))
	      dummy_2d(i,k,3)=0.5*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3))
	      dummy_2d(i,k,4)=dummy_3d(i,j,k,4)
	    End Do; End Do
	    If(irmean==1) Then
	      Do j=islicexyz,islicexyz; dummy_2d(:,:,1)=dummy_2d(:,:,1)-umn(j); dummy_2d(:,:,3)=dummy_2d(:,:,3)-umn(j); End Do
	    End If
	  End If
	  If(itec==1) Then
 	    WRITE(501,*) 'zone t="slice" i=',dim1+2,', j=',dim2+2,'f=point'
	    Do j=0,dim2+1; Do i=0,dim1+1
	      Write(501,'(10f20.10)') xyz1(i),xyz2(j),dummy_2d(i,j,1:4)
	    End Do; End Do
	  Else
!	    paraview format
	    Write(501,'(A,3I10)') 'DIMENSIONS',dim1vtk,dim2vtk,dim3vtk
	    Write(501,'(A,I10,A)') 'X_COORDINATES', dim1vtk, ' float'
	    Do i=1,dim1vtk; Write(501,'(E15.7)') xc(i); End Do
	    Write(501,'(A,I10,A)') 'Y_COORDINATES', dim2vtk, ' float'
	    Do j=1,dim2vtk; Write(501,'(E15.7)') yc(j); End Do
	    Write(501,'(A,I10,A)') 'Z_COORDINATES', dim3vtk, ' float'
	    Do k=1,dim3vtk; Write(501,'(E15.7)') zc(k); End Do
	    Write(501,'(A,I10)') 'POINT_DATA',dim1vtk*dim2vtk*dim3vtk
	    Write(501,'(A)') 'SCALARS u float 1'
	    Write(501,'(A)') 'LOOKUP_TABLE default'
	    Write(501,'(10F20.10)') ((dummy_2d(i,j,1),i=1,dim1),j=1,dim2)
	    Write(501,'(A)') 'SCALARS v float 1'
	    Write(501,'(A)') 'LOOKUP_TABLE default'
	    Write(501,'(10F20.10)') ((dummy_2d(i,j,2),i=1,dim1),j=1,dim2)
	    Write(501,'(A)') 'SCALARS w float 1'
	    Write(501,'(A)') 'LOOKUP_TABLE default'
	    Write(501,'(10F20.10)') ((dummy_2d(i,j,3),i=1,dim1),j=1,dim2)
	    Write(501,'(A)') 'SCALARS p float 1'
	    Write(501,'(A)') 'LOOKUP_TABLE default'
	    Write(501,'(10F20.10)') ((dummy_2d(i,j,4),i=1,dim1),j=1,dim2)
	  End If
	End If

	Deallocate(dummy_2d,xyz1,xyz2)
	Close(501)

	END SUBROUTINE write_2d_slice
!*******************************************************************************
!*******************************************************************************
	SUBROUTINE open_files3d
	Use write_post
	CHARACTER (LEN=100), DIMENSION(29) :: dum
	CHARACTER (LEN=1), DIMENSION(4) :: uvwp
	CHARACTER :: filedum*40,one*1,three*3
	INTEGER :: i,dim1,dim2
	uvwp=(/'u','v','w','p'/)
!
	If (ivort3d == 1) then
	  If(itec==1) Then
	    Open(Unit=18, File='vort-ins.plt')
!	    WRITE(18,*) ' title="vort-ins.plt"'
!	    WRITE(18,*) ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","w<sub>x</sub>"'
	    WRITE(18,*) ' variables="x","y","z","wx"'

	    Open(Unit=301,File='vorty-ins.plt' )
	    WRITE(301,*) ' variables="x","y","z","wy"'
	    Open(Unit=302,File='vortz-ins.plt' )
	    WRITE(302,*) ' variables="x","y","z","wz"'
	  Else
	    Open(18, File='vort-ins.vtk', status='UNKNOWN')
	    Write(18,'(A)') '# vtk DataFile Version 2.0'
	    Write(18,'(A)') 'Volume example'
	    Write(18,'(A)') 'ASCII'
	    Write(18,'(A)') 'DATASET RECTILINEAR_GRID'
!
	    Open(301, File='vorty-ins.vtk', status='UNKNOWN')
	    Write(301,'(A)') '# vtk DataFile Version 2.0'
	    Write(301,'(A)') 'Volume example'
	    Write(301,'(A)') 'ASCII'
	    Write(301,'(A)') 'DATASET RECTILINEAR_GRID'
!
	    Open(302, File='vortz-ins.vtk', status='UNKNOWN')
	    Write(302,'(A)') '# vtk DataFile Version 2.0'
	    Write(302,'(A)') 'Volume example'
	    Write(302,'(A)') 'ASCII'
	    Write(302,'(A)') 'DATASET RECTILINEAR_GRID'
	  End If

	End If
!
	If (ilamq == 1) then
	  If(itec==1) Then
	    Open(23, File='lambdaQ-ins.plt', status='UNKNOWN')
!	    WRITE(23,*) ' title="lambdaQ-ins.plt"'
!	    WRITE(23,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","wx-neg","wx-pos","lamada2"'
	    WRITE(23,'(120A)') ' variables="x","y","z","wx-neg","wx-pos","lamada2"'
	  Else
	    Open(23, File='lambdaQ-ins.vtk', status='UNKNOWN')
	    Write(23,'(A)') '# vtk DataFile Version 2.0'
	    Write(23,'(A)') 'Volume example'
	    Write(23,'(A)') 'ASCII'
	    Write(23,'(A)') 'DATASET RECTILINEAR_GRID'
	  End If

	  Open(25, File='v1-ins.plt', status='UNKNOWN')
	  WRITE(25,*) ' title="v1-ins.plt"'
	  WRITE(25,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","v1(1)","v1(2)","v1(3)"'

	  Open(26, File='struct-ins.plt', status='UNKNOWN')
	  WRITE(26,*) ' title="v1-ins.plt"'
	  WRITE(26,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","v1(1)","v1(2)","v1(3)"'

	  Open(27, File='checkpts-ins.plt', status='UNKNOWN')
	  WRITE(27,*) ' title="v1-ins.plt"'
	  WRITE(27,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>"'

	  Open(28, File='lambda2-rms.plt', status='UNKNOWN')
	  WRITE(28,*) ' title="lambda2-rms.plt"'
	  WRITE(28,'(120A)') ' variables="y<sup>+</sup>","-<greek>l</greek><sub>2</sub>-mean","<greek>l</greek><sub>2</sub>-rms","-<greek>l</greek><sub>2-</sub>-mean","-<greek>l</greek><sub>2+</sub>-mean","<greek>f<\greek>"'

	  Open(29, File='indiv-struc.plt', status='UNKNOWN')
	  WRITE(29,*) ' title="indiv-strucplt"'
	  WRITE(29,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","<greek>l</greek><sub>2</sub>"'

	  Open(31, File='pmin.plt', status='UNKNOWN')
	  WRITE(31,*) ' title="pmin.plt"'
	  WRITE(31,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","pmin"'
	End If

	If (icorr ==1) then
	  Open(57, File='xcorr.plt', status='UNKNOWN')
	  WRITE(57,*) ' title="xcorr.plt"'
	  WRITE(57,*) ' variables="x","R<sub>11</sub>","R<sub>22</sub>","R<sub>33</sub>"'
!
	  Open(58, File='zcorr.plt', status='UNKNOWN')
	  WRITE(58,*) ' title="zcorr.plt"'
	  WRITE(58,*) ' variables="z","R<sub>11</sub>","R<sub>22</sub>","R<sub>33</sub>"'
!
	  Open(56, File='xspec.plt', status='UNKNOWN')
	  WRITE(56,*) ' title="xspec.plt"'
	  WRITE(56,*) ' variables="k<sub>x</sub>","E<sub>11</sub>","E<sub>22</sub>","E<sub>33</sub>"'
!
	  Open(59, File='zspec.plt', status='UNKNOWN')
	  WRITE(59,*) ' title="zspec.plt"'
	  WRITE(59,*) ' variables="k<sub>z</sub>","E<sub>11</sub>","E<sub>22</sub>","E<sub>33</sub>"'
!
	  Open(512, File='xspecp.plt', status='UNKNOWN')
	  WRITE(512,'(A)') ' title="xspecp.plt"'
	  WRITE(512,'(A)') ' variables="k<sub>x</sub><sup>+</sup>","E<sub>11</sub><sup>+</sup>","E<sub>22</sub><sup>+</sup>","E<sub>33</sub><sup>+</sup>"'
!
	  Open(513, File='zspecp.plt', status='UNKNOWN')
	  WRITE(513,'(A)') ' title="zspecp.plt"'
	  WRITE(513,'(A)') ' variables="k<sub>z</sub><sup>+</sup>","E<sub>11</sub><sup>+</sup>","E<sub>22</sub><sup>+</sup>","E<sub>33</sub><sup>+</sup>"'
!	  Open(305, File='spectra-contour.plt', status='UNKNOWN')
!	  WRITE(305,*) ' title="spectra-contour"'
!          WRITE(305,*) ' variables="y","k","spxR11","spxR22","spxR33"'
! 	  WRITE(305,*) 'zone i=',nx/2,', j=',ny,',f=point'
	End If

	If(iplot3d == 1) Then
	  If(itec == 1) Then; Open(Unit=303,File='3dfield.plt')
	  Else; Open(Unit=303,File='3dfield.vtk'); End If
	  If(itec == 1) Then
	    WRITE(303,'(120A)') ' variables="x","y","z","u","v","w","p"'
 	    WRITE(303,*) 'zone DATAPACKING=BLOCK, i=',nx,', j=',jvis,', k=',nz
	  Else
	    Write(303,'(A)') '# vtk DataFile Version 2.0'
	    Write(303,'(A)') 'Volume example'
	    Write(303,'(A)') 'ASCII'
	    Write(303,'(A)') 'DATASET RECTILINEAR_GRID'
	  End If
	End If
!
	If(iplot2d.Eq.1.Or.lp_snap_2d.Eq.1) Then
	  If(lp_snap_xy.Eq.1) Then
	    filedum='plane-xy.plt';dum(1)='x';dum(2)='y';dim1=nx+2;dim2=ny+2
	  Else If(lp_snap_yz.Eq.1) Then
	    filedum='plane-yz.plt';dum(1)='y';dum(2)='z';dim1=ny+2;dim2=nz+2
	  Else If(lp_snap_xz.Eq.1) Then
	    filedum='plane-xz.plt';dum(1)='x';dum(2)='z';dim1=nx+2;dim2=nz+2
	  End If
	  dum(3)='u';dum(4)='v';dum(5)='w';dum(6)='p'
  	  If(itec==1) Then
	    Open(501, File=filedum, status='UNKNOWN')
	    WRITE(501,*) ' title="',Trim(filedum),'"'
Write(*,*) 'iconstruc=',iconstruc
	    If(iconstruc.Eq.0) Then
	      WRITE(501,*) ' variables="',Trim(dum(1)),'","',Trim(dum(2)),'","',Trim(dum(3)),'","',Trim(dum(4)),'","',Trim(dum(5)),'","',Trim(dum(6)),'"'
	    Else
	      WRITE(501,*) ' variables="',Trim(dum(2)),'","',Trim(dum(3)),'","',Trim(dum(4)),'","',Trim(dum(5)),'","',Trim(dum(6)),'"'
	    End If
	  Else
	    Open(501, File='plane2D.vtk', status='UNKNOWN')
	    Write(501,'(A)') '# vtk DataFile Version 2.0'
	    Write(501,'(A)') 'Volume example'
	    Write(501,'(A)') 'ASCII'
	    Write(501,'(A)') 'DATASET RECTILINEAR_GRID'
	  End If
	End If
!
	If(iplot1d.Eq.1.Or.lp_snap_1d.Eq.1) Then
	  If(lp_snap_x.Eq.1) Then
	    filedum='line-x.plt';dum(1)='x';dim1=nx+2
	  Else If(lp_snap_y.Eq.1) Then
	    filedum='line-y.plt';dum(1)='y';dim1=ny+2
	  Else If(lp_snap_z.Eq.1) Then
	    filedum='line-z.plt';dum(1)='z';dim2=nz+2
	  End If
	  dum(2)='u';dum(3)='v';dum(4)='w';dum(5)='p'
	  Open(502, File=filedum, status='UNKNOWN')
	  WRITE(502,*) ' title="',Trim(filedum),'"'
	  If(iconstruc.Eq.0) Then
	    WRITE(502,*) ' variables="',Trim(dum(1)),'","',Trim(dum(2)),'","',Trim(dum(3)),'","',Trim(dum(4)),'","',Trim(dum(5)),'"'
	  Else
	    WRITE(502,*) ' variables="',Trim(dum(1)),'","time","',Trim(dum(2)),'","',Trim(dum(3)),'","',Trim(dum(4)),'","',Trim(dum(5)),'"'
	  End If
	End If
!
	END SUBROUTINE open_files3d
!*******************************************************************************
!*******************************************************************************
	SUBROUTINE write_modulation(AM)
!
	Real(mytype), Dimension(1:jmodu,1:jmodu,1:3), Intent(in) :: AM
	Integer :: j1,j2
!
	Open(Unit=610,File='modulation_2points.plt')
	WRITE(610,'(320A)') ' variables="y<sub>2</sub>","y<sub>1</sub>","AM","<greek>F</greek><sub>x</sub>","<greek>F</greek><sub>z</sub>"'
 	WRITE(610,*) 'zone i=',jmodu,', j=',jmodu,',f=point'
	Do j1=1,jmodu; Do j2=1,jmodu
	  Write(610,'(5E15.7)') yc(j2),yc(j1),AM(j2,j1,1:3)
	End Do; End Do
	Close(610)

!
	END SUBROUTINE write_modulation
!***************************************************************************
!*******************************************************************************
        Subroutine plasma_decompose
  	Include "fftw3.f"
! 	this subroutine decomposes the plasma force region to Fourier modes
        Real(mytype), Allocatable, Dimension(:,:) :: fbody,fmode
	Real(mytype), Allocatable, Dimension(:) :: cactu,hdisp,hdisn,dum_1d
	Complex(mytype), Allocatable, Dimension(:) :: outfftz
	Real(mytype) :: lambda,sigma,St,expo,rayl,raylstrp,raylstrn
        Integer :: i,j,k,n,nactu
        Integer*8 :: plan1,plan2
!
	nactu=24
	Allocate(outfftz(nz/2+1),fbody(1:jvis,1:nz),fmode(1:jvis,1:nz),dum_1d(1:nz))
!	plasma body force
        Allocate(cactu(1:nactu),hdisp(1:nz),hdisn(1:nz))
        cactu(:)=0.0e0; hdisp(:)=0.0e0; hdisn(:)=0.0e0
	lambda=45.0e0; sigma=0.07e0; St=1.0e0
        expo=0.0e0; rayl=0.0e0
	fmode(:,:)=0.0e0
        Do n=1,nactu
          cactu(n)=width/real(nactu,mytype)*real(n-1,mytype)
        End Do
        Do k=1,nz; Do n=1,nactu
          raylstrp=modulo(z(k)-cactu(n)+width,width)
          raylstrn=modulo(cactu(n)-z(k)+width,width)
          hdisp(k)=hdisp(k)+raylstrp/sigma**2*exp(-raylstrp**2/(2.0*sigma**2))
          hdisn(k)=hdisn(k)-raylstrn/sigma**2*exp(-raylstrn**2/(2.0*sigma**2))
        End Do; End Do
	Do j=1,ny/2
	  expo=lambda*exp(-lambda*yc(j))
	  Do k=1,nz
	  fbody(j,k)=St/lambda*sigma*sqrt(exp(1.0e0))*expo*hdisp(k)
	End Do; End Do
!
	Open(unit=628,file='spec_fbody.plt')
	Write(628,'(10A)') ' variables="<greek>l</greek><sub>z</sub><sup>+</sup>", &
	"yc","<greek>F</greek>","k<sub>z</sub><greek>F</greek>"'
	Write(628,*) 'zone i=',nz/2,', j=',jvis,',f=point'
!
	Do j=1,jvis
!	forward Fourier transform
	  dum_1d(1:nz)=fbody(j,1:nz)
	  Call dfftw_plan_dft_r2c_1d(plan1,nz,dum_1d,outfftz,FFTW_ESTIMATE)
	  Call dfftw_execute(plan1)
	  Call dfftw_destroy_plan(plan1)
	  outfftz(:)=outfftz(:)/Real(nz,mytype)
	  Do k=2,nz/2+1
	    Write(628,'(4F20.10)') width/Real(k-1,mytype)*retau_avg_sa,yc(j),Abs(outfftz(k)), &
		2.0e0*pi/(width/Real(k-1,mytype)*retau_avg_sa)*Abs(outfftz(k))
	  End Do
!	filtering (harmonics k=1,13,25,37,49,61,73 ...)
	Do k=1,nz/2+1
!	  If(Mod(k-1,12).Eq.0.And.k<=nz) Then
	  If(k==1) Then
	  Else; outfftz(k)=0.0e0; End If
	End Do
!	backward Fourier transform
	  Call dfftw_plan_dft_c2r_1d(plan2,nz,outfftz,dum_1d,FFTW_ESTIMATE)
	  Call dfftw_execute(plan2)
	  Call dfftw_destroy_plan(plan2)
	  fmode(j,1:nz)=dum_1d(1:nz)
	End Do
	Close(628)
	
!	Write the plasma body force out
	Open(627, File='fbody.vtk', status='UNKNOWN')
	Write(627,'(A)') '# vtk DataFile Version 2.0'
	Write(627,'(A)') 'Volume example'
	Write(627,'(A)') 'ASCII'
	Write(627,'(A)') 'DATASET RECTILINEAR_GRID'
	Write(627,'(A,3I10)') 'DIMENSIONS',1,jvis,nz
	Write(627,'(A,I10,A)') 'X_COORDINATES', 1, ' float'
	Do i=1,1; Write(627,'(E15.7)') xc(i); End Do
	Write(627,'(A,I10,A)') 'Y_COORDINATES', jvis, ' float'
	Do j=1,jvis; Write(627,'(E15.7)') yc(j); End Do
	Write(627,'(A,I10,A)') 'Z_COORDINATES', nz, ' float'
	Do k=1,nz; Write(627,'(E15.7)') zc(k); End Do
	Write(627,'(A,I10)') 'POINT_DATA',nz*jvis
	Write(627,'(A)') 'SCALARS fbody float 1'
	Write(627,'(A)') 'LOOKUP_TABLE default'
	Write(627,'(10E20.10)') ((fbody(j,k),j=1,jvis),k=1,nz)
	Write(627,'(A)') 'SCALARS fmode float 1'
	Write(627,'(A)') 'LOOKUP_TABLE default'
	Write(627,'(10E20.10)') ((fmode(j,k),j=1,jvis),k=1,nz)
	Close(627)
!
	Deallocate(fbody,fmode,cactu,hdisp,hdisn)

        End Subroutine plasma_decompose

!*******************************************************************************
	SUBROUTINE calc_quad
	Real(mytype), Dimension(:,:,:,:), Allocatable :: field_3d
	Real(mytype), Allocatable, Dimension(:,:) :: umean,vmean,uv_2d,q1_2d,q2_2d,q3_2d,q4_2d
	Real(mytype), Allocatable, Dimension(:) :: q1,q2,q3,q4,uv
	Real(mytype) :: utmp,vtmp,dum
	Integer :: i,j,k,idir,intdum,m,nsnaps,nactu,kp,iw
	Character :: chdum*8,file3d*60
!
	nsnaps=1; nactu=12; iw=1
	Allocate(field_3d(-1:nx+2,-1:ny+2,-1:nz+2,1:2))
	Allocate(umean(1:ny,1:nz),vmean(1:ny,1:nz),q1(1:ny),q2(1:ny),q3(1:ny),q4(1:ny),uv(1:ny))
	Allocate(q1_2d(1:ny,1:nz),q2_2d(1:ny,1:nz),q3_2d(1:ny,1:nz),q4_2d(1:ny,1:nz),uv_2d(1:ny,1:nz))
	umean(:,:)=0.0e0; vmean(:,:)=0.0e0
	q1(:)=0.0e0; q2(:)=0.0e0; q3(:)=0.0e0; q4(:)=0.0e0
	q1_2d(:,:)=0.0e0; q2_2d(:,:)=0.0e0; q3_2d(:,:)=0.0e0; q4_2d(:,:)=0.0e0
!
Do m=1,nsnaps
  Write(chdum,'(I5)') m
  file3d='snap3d_'//Trim(Adjustl(chdum))//'.dat'
  If(nsnaps==1) file3d=restrt_file
  Open(Unit=14, file=file3d, Status='OLD', Form='UNFormatTED',access='stream')
  READ(14) intdum,intdum,intdum,dum,dum,dum,dum,dum,dum,dum,dum
  Write(*,'(A,I5,A,I5,F20.10)') 'snap:',m,'  header:',intdum,dum
  READ(14) field_3d(:,:,:,1) ! u
  READ(14) field_3d(:,:,:,2) ! v
!
	Do k=1,nz; Do j=1,ny
	  umean(j,k)=stat1d_sa(j,1,1,k)  !Sum(dummy_3d(1:nx,j,1:nz))/Real(nx*nz)
	  vmean(j,k)=stat1d_sa(j,2,1,k)  !Sum(dummy_3d(1:nx,j,1:nz))/Real(nx*nz)
	End Do; End Do
!
	Do j=1,ny
	  Do k=1,nz; Do i=1,nx
	    utmp=0.5e0*(field_3d(i,j,k,1)+field_3d(i+1,j,k,1))-umean(j,k)
	    vtmp=0.5e0*(field_3d(i,j,k,2)+field_3d(i,j+1,k,2))-vmean(j,k)
	    uv(j)=uv(j)+utmp*vtmp
	    If(utmp>0.0e0) Then
	      If(vtmp>0.0e0) Then
		q1(j)=q1(j)+utmp*vtmp
	      Else If(vtmp<0.0e0) Then
		q4(j)=q4(j)+utmp*vtmp
	      End If
	    Else If(utmp<0.0e0) Then
	      If(vtmp>0.0e0) Then
		q2(j)=q2(j)+utmp*vtmp
	      Else If(vtmp<0.0e0) Then
		q3(j)=q3(j)+utmp*vtmp
	      End If
	    End If
	  End Do; End Do
	End Do
!
	Do j=1,ny; Do k=1,nz
	  Do i=1,nx
	    utmp=0.5e0*(field_3d(i,j,k,1)+field_3d(i+1,j,k,1))-umean(j,k)
	    vtmp=0.5e0*(field_3d(i,j,k,2)+field_3d(i,j+1,k,2))-vmean(j,k)
	    uv_2d(j,k)=uv_2d(j,k)+utmp*vtmp
	    If(utmp>0.0e0) Then
	      If(vtmp>0.0e0) Then
		q1_2d(j,k)=q1_2d(j,k)+utmp*vtmp
	      Else If(vtmp<0.0e0) Then
		q4_2d(j,k)=q4_2d(j,k)+utmp*vtmp
	      End If
	    Else If(utmp<0.0e0) Then
	      If(vtmp>0.0e0) Then
		q2_2d(j,k)=q2_2d(j,k)+utmp*vtmp
	      Else If(vtmp<0.0e0) Then
		q3_2d(j,k)=q3_2d(j,k)+utmp*vtmp
	      End If
	    End If
	  End Do
	End Do; End Do
!
  READ(14) field_3d(:,:,:,1) ! w (place holder)
  READ(14) field_3d(:,:,:,2) ! p (place holder)
  Close(14)
End Do
!	Write 1d files out
	Open(unit=629,file='quad-ins.plt')
	Write(629,'(10A)') ' variables="yc","uv","Q1","Q2","Q3","Q4"'
	Do j=1,ny
	  uv(j)=uv(j)/Real(nx*nz*nsnaps,mytype)
	  q1(j)=q1(j)/Real(nx*nz*nsnaps,mytype)
	  q2(j)=q2(j)/Real(nx*nz*nsnaps,mytype)
	  q3(j)=q3(j)/Real(nx*nz*nsnaps,mytype)
	  q4(j)=q4(j)/Real(nx*nz*nsnaps,mytype)
	  Write(629,'(6F20.10)') yc(j),uv(j),q1(j),q2(j),q3(j),q4(j)
	End Do
	Close(629)
!	Write 2d files out
	Open(unit=630,file='quad2d-ins.plt')
	Write(630,'(10A)') ' variables="zc","yc","uv","Q1","Q2","Q3","Q4"'
	Write(630,*) 'zone i=',nz,', j=',ny,',f=point'
!	 Write(*,'(A,5F10.5)') 'quad2d-ins:',yc(j)
	Do j=1,ny; Do k=1,nz
	  uv_2d(j,k)=uv_2d(j,k)/Real(nx*nsnaps,mytype)
	  q1_2d(j,k)=q1_2d(j,k)/Real(nx*nsnaps,mytype)
	  q2_2d(j,k)=q2_2d(j,k)/Real(nx*nsnaps,mytype)
	  q3_2d(j,k)=q3_2d(j,k)/Real(nx*nsnaps,mytype)
	  q4_2d(j,k)=q4_2d(j,k)/Real(nx*nsnaps,mytype)
	  Write(630,'(7F20.10)') zc(k),yc(j),uv_2d(j,k),q1_2d(j,k),q2_2d(j,k),q3_2d(j,k),q4_2d(j,k)
	End Do; End Do
	Close(630)
!	Write averaged 2d files out
	Open(unit=630,file='quad2d_avg-ins.plt')
	Write(630,'(10A)') ' variables="zc","yc","uv","Q1","Q2","Q3","Q4"'
	Write(630,*) 'zone i=',nz/nactu,', j=',ny,',f=point'
	Do j=1,ny; Do k=1,nz/nactu
	  Do i=1,nactu
	    kp=(i-1)*nz/nactu+k
	    uv_2d(j,k)=uv_2d(j,k)+uv_2d(j,kp)
	    q1_2d(j,k)=q1_2d(j,k)+q1_2d(j,kp)
	    q2_2d(j,k)=q2_2d(j,k)+q2_2d(j,kp)
	    q3_2d(j,k)=q3_2d(j,k)+q3_2d(j,kp)
	    q4_2d(j,k)=q4_2d(j,k)+q4_2d(j,kp)
	  End Do
	End Do; End Do
	uv_2d(:,:)=uv_2d(:,:)/Real(nsnaps,mytype)
	q1_2d(:,:)=q1_2d(:,:)/Real(nsnaps,mytype)
	q2_2d(:,:)=q2_2d(:,:)/Real(nsnaps,mytype)
	q3_2d(:,:)=q3_2d(:,:)/Real(nsnaps,mytype)
	q4_2d(:,:)=q4_2d(:,:)/Real(nsnaps,mytype)
	Do j=1,ny; Do k=1,nz/nactu
	  Write(630,'(7F20.10)') zc(k),yc(j),uv_2d(j,k),q1_2d(j,k),q2_2d(j,k),q3_2d(j,k),q4_2d(j,k)
	End Do; End Do
	Close(630)
!
	Deallocate(field_3d,umean,vmean,q1,q2,q3,q4,uv)
	Deallocate(q1_2d,q2_2d,q3_2d,q4_2d,uv_2d)
!
	END SUBROUTINE calc_quad

!*******************************************************************************
	SUBROUTINE calc_am1p
! this subroutine calculates the one point amplitude modulation based on Mathis_etal_JFM2010
	Include "fftw3.f"
	Real(mytype), Dimension(:,:,:), Allocatable :: uu
	Complex(mytype), Allocatable, Dimension(:) :: outfftx,outfftx2
	Real(mytype), Allocatable, Dimension(:) :: ux,ulx,usx,uex,uhx,uelx,am,arr1d
	Character :: chdum*8,file3d*60
	Real(mytype) :: dum,yam,theta,alpha,beta,uol,ustar,corr1,corr2,tmp
	Integer :: intdum,nsnaps
	Integer :: m,i,j,k,jp,ip,kp,ind
        Integer*8 :: planx
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
! predictive model (Marusic_etal_Science2010)
	alpha=1.0e0; beta=0.5e0
!
	nsnaps=1
	Allocate(uu(-1:nx+2,-1:ny+2,-1:nz+2))
	Allocate(outfftx(nx/2+1),outfftx2(nx/2+1),ux(1:nx),ulx(1:nx),usx(1:nx), &
	uex(1:nx),uhx(1:nx),uelx(1:nx),arr1d(1:nx))
	Allocate(am(1:ny))
!
	am(:)=0.0e0
!
	Do m=1,nsnaps
  	  Write(chdum,'(I5)') m
!  	  file3d='snap3d_'//Trim(Adjustl(chdum))//'.dat'
	  file3d='../restrt-start0-stand-run8-out.rst'
  	  If(nsnaps==1) file3d=restrt_file
  	  Open(Unit=14, file=file3d, Status='OLD', Form='UNFormatTED',access='stream')
  	  Read(14) intdum,intdum,intdum,dum,dum,dum,dum,dum,dum,dum,dum
  	  Write(*,'(A,I5,A,I5,F20.10)') 'snap:',m,'  header:',intdum,dum
  	  Read(14) uu(:,:,:) !u
  	  Close(14)
	End Do
!
!	find the y location for am calculation
	yam=15.0e0; ip=nx/2
	Do j=1,ny/2
	  If(yplus(j)<yam.And.yplus(j+1)>yam) Then
	    jp=j; Exit
	  End If
	End Do
!
	Do j=1,ny
	  Write(*,*) 'j=',j
!!	  shuffle the signal (Schlatter & Orlu 2010)
!	  Do k=1,nz
!	    arr1d(1:nx)=uu(1:nx,j,k)
!	    Do i=nx,1,-1
!	      CALL RANDOM_NUMBER(theta)
!	      ind=Max(Int(theta*i),1)
!	      tmp=arr1d(i); arr1d(i)=arr1d(ind); arr1d(ind)=tmp
!	    End Do
!	    uu(1:nx,j,k)=arr1d(1:nx)
!	  End Do
!
	  corr1=0.0e0; corr2=0.0e0
	  Do k=1,nz
	  kp=k
	  ux(1:nx)=uu(1:nx,j,kp)
!	  an artificial signal
!	  Do i=1,nx
!	    ustar=Sin(2.0e0*pi/(length/20.0e0)*x(i))
!	    uol=Sin(2.0e0*pi/(length/2.0e0)*x(i))
!	    ux(i)=ustar*(1.0e0+beta*uol)+alpha*uol
!	  End Do
!	  FFT
	  Call dfftw_plan_dft_r2c_1d(planx,nx,ux,outfftx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)
	  Call dfftw_destroy_plan(planx)	
	  outfftx(:)=outfftx(:)/Real(nx,mytype)
!!	  scramble the signal (Mathis et al 2009)
!	  Do i=1,nx/2+1
!	    CALL RANDOM_NUMBER(theta)
!	    theta=theta*2.0e0*pi
!	    outfftx(i)=outfftx(i)*Cmplx(Cos(theta),Sin(theta))
!	  End Do
!
!	  filtering
	  outfftx(1)=0.0e0
!	  inverse FFT
	  outfftx2(:)=outfftx(:)
	  Call dfftw_plan_dft_c2r_1d(planx,nx,outfftx2,ux,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)
	  Call dfftw_destroy_plan(planx)
!	  filtering
	  outfftx2(:)=outfftx(:)
	  Do i=2,nx/2+1
	    If(length/Real(i-1)>=1.0e0) outfftx2(i)=0.0e0
	  End Do
!	  inverse FFT
	  Call dfftw_plan_dft_c2r_1d(planx,nx,outfftx2,usx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)
	  Call dfftw_destroy_plan(planx)
!
	  ulx(:)=ux(:)-usx(:)
!	  Hilbert transform
	  uhx(1:nx)=usx(1:nx)
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
	  uex(:)=Sqrt(uhx(:)**2+uex(:)**2)
!	  filtering the envelop
	  Call dfftw_plan_dft_r2c_1d(planx,nx,uex,outfftx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx) 
	  Call dfftw_destroy_plan(planx)
	  outfftx(:)=outfftx(:)/Real(nx,mytype)
	  outfftx(1)=0.0e0
	  Do i=2,nx/2+1
	    If(length/Real(i-1)<1.0e0) outfftx(i)=0.0e0
	  End Do
	  Call dfftw_plan_dft_c2r_1d(planx,nx,outfftx,uelx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)  
	  Call dfftw_destroy_plan(planx)
!	  find the AM
	  Do i=1,nx
	    corr1=corr1+ulx(i)**2
	    corr2=corr2+uelx(i)**2
	    am(j)=am(j)+ulx(i)*uelx(i)
	  End Do
	  End Do ! k
	  corr1=corr1/Real(nx*nz,mytype); corr2=corr2/Real(nx*nz,mytype)
	  am(j)=am(j)/Real(nx*nz,mytype)/Sqrt(corr1)/Sqrt(corr2)
	End Do ! j
!
	Open(unit=634,file='ux-ins.plt')
	Write(634,'(10A)') ' variables="x","u","u<sub>l</sub>","u<sub>s</sub>","u<sub>e</sub>","u<sub>el</sub>"'
	Do i=1,nx
	  Write(634,'(6F20.10)') x(i),ux(i),ulx(i),usx(i),uex(i),uelx(i)
	End Do
	Close(634)
!
	Open(unit=635,file='AM_Mathis.plt')
	Write(635,'(10A)') ' variables="yc","AM","y<sup>+</sup>"'
	Do j=1,ny
	  Write(635,'(6F20.10)') yc(j),am(j),ycplus(j)
	End Do
	Close(635)
!
	Deallocate(uu,ux,ulx,usx,uex,uelx,uhx,outfftx,outfftx2,am,arr1d)	
!
	END SUBROUTINE calc_am1p

!*******************************************************************************
	SUBROUTINE calc_am1p_z
! this subroutine calculates the one point amplitude modulation based on Mathis_etal_JFM2010 for z direction
	Include "fftw3.f"
	Real(mytype), Dimension(:,:,:), Allocatable :: uu
	Complex(mytype), Allocatable, Dimension(:) :: outfftx,outfftx2
	Real(mytype), Allocatable, Dimension(:) :: ux,ulx,usx,uex,uhx,uelx,am
	Character :: chdum*8,file3d*60
	Real(mytype) :: dum,yam,theta,alpha,beta,uol,ustar,corr1,corr2
	Integer :: intdum,nsnaps
	Integer :: m,i,j,k,jp,ip,kp
        Integer*8 :: planx
! predictive model (Marusic_etal_Science2010)
	alpha=1.0e0; beta=0.5e0
!
	nsnaps=1
	Allocate(uu(-1:nx+2,-1:ny+2,-1:nz+2))
	Allocate(outfftx(nz/2+1),outfftx2(nz/2+1),ux(1:nz),ulx(1:nz),usx(1:nz),uex(1:nz),uhx(1:nz),uelx(1:nz))
	Allocate(am(1:ny))
!
	am(:)=0.0e0
!
	Do m=1,nsnaps
  	  Write(chdum,'(I5)') m
!  	  file3d='snap3d_'//Trim(Adjustl(chdum))//'.dat'
	  file3d='../restrt-start0-stand-run8-out.rst'
  	  If(nsnaps==1) file3d=restrt_file
  	  Open(Unit=14, file=file3d, Status='OLD', Form='UNFormatTED',access='stream')
  	  Read(14) intdum,intdum,intdum,dum,dum,dum,dum,dum,dum,dum,dum
  	  Write(*,'(A,I5,A,I5,F20.10)') 'snap:',m,'  header:',intdum,dum
  	  Read(14) uu(:,:,:) !u
  	  Close(14)
	End Do
!
!	find the y location for am calculation
	yam=15.0e0; ip=nx/2
	Do j=1,ny/2
	  If(yplus(j)<yam.And.yplus(j+1)>yam) Then
	    jp=j; Exit
	  End If
	End Do
!
	Do j=1,ny !jp,jp
	  Write(*,*) 'j=',j
	  corr1=0.0e0; corr2=0.0e0
	  Do i=1,nx
	  ip=i
	  ux(1:nz)=uu(ip,j,1:nz)
!	  FFT
	  Call dfftw_plan_dft_r2c_1d(planx,nz,ux,outfftx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)
	  Call dfftw_destroy_plan(planx)	
	  outfftx(:)=outfftx(:)/Real(nz,mytype)
!	  filtering
	  outfftx(1)=0.0e0
!	  inverse FFT
	  outfftx2(:)=outfftx(:)
	  Call dfftw_plan_dft_c2r_1d(planx,nz,outfftx2,ux,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)
	  Call dfftw_destroy_plan(planx)
!	  filtering
	  outfftx2(:)=outfftx(:)
	  Do k=2,nz/2+1
	    If(width/Real(k-1)>=1.0e0) outfftx2(k)=0.0e0
	  End Do
!	  inverse FFT
	  Call dfftw_plan_dft_c2r_1d(planx,nz,outfftx2,usx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)
	  Call dfftw_destroy_plan(planx)
!
	  ulx(:)=ux(:)-usx(:)
!	  Hilbert transform
	  uhx(:)=usx(:)
	  Call dfftw_plan_dft_r2c_1d(planx,nz,uhx,outfftx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)
	  Call dfftw_destroy_plan(planx)
	  outfftx(:)=outfftx(:)/Real(nz,mytype)
	  theta=-pi/2.0e0
	  Do k=2,nz/2+1
	    outfftx(k)=outfftx(k)*Cmplx(Cos(theta),Sin(theta))
	  End Do
	  Call dfftw_plan_dft_c2r_1d(planx,nz,outfftx,uex,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)  
	  Call dfftw_destroy_plan(planx)
	  uex(:)=Sqrt(uhx(:)**2+uex(:)**2)
!	  filtering the envelop
	  Call dfftw_plan_dft_r2c_1d(planx,nz,uex,outfftx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx) 
	  Call dfftw_destroy_plan(planx)
	  outfftx(:)=outfftx(:)/Real(nz,mytype)
	  outfftx(1)=0.0e0
	  Do k=2,nz/2+1
	    If(width/Real(k-1)<1.0e0) outfftx(k)=0.0e0
	  End Do
	  Call dfftw_plan_dft_c2r_1d(planx,nz,outfftx,uelx,FFTW_ESTIMATE)
	  Call dfftw_execute(planx)  
	  Call dfftw_destroy_plan(planx)
!	  find the AM
	  Do k=1,nz
	    corr1=corr1+ulx(k)**2
	    corr2=corr2+uelx(k)**2
	    am(j)=am(j)+ulx(k)*uelx(k)
	  End Do
	  End Do ! i
	  corr1=corr1/Real(nx*nz,mytype); corr2=corr2/Real(nx*nz,mytype)
	  am(j)=am(j)/Real(nx*nz,mytype)/Sqrt(corr1)/Sqrt(corr2)
	End Do ! j
!
	Open(unit=634,file='uz-ins.plt')
	Write(634,'(10A)') ' variables="zc","u","u<sub>l</sub>","u<sub>s</sub>","u<sub>e</sub>","u<sub>el</sub>"'
	Do k=1,nz
	  Write(634,'(6F20.10)') zc(k),ux(k),ulx(k),usx(k),uex(k),uelx(k)
	End Do
	Close(634)
!
	Open(unit=635,file='AM_Mathis_z.plt')
	Write(635,'(10A)') ' variables="yc","AM","y<sup>+</sup>"'
	Do j=1,ny
	  Write(635,'(6F20.10)') yc(j),am(j),ycplus(j)
	End Do
	Close(635)
!
	Deallocate(uu,ux,ulx,usx,uex,uelx,uhx,outfftx,outfftx2,am)	
!
	END SUBROUTINE calc_am1p_z

!*******************************************************************************
	SUBROUTINE calc_am2p
! this subroutine calculates the two point amplitude modulation based on Bernardini_Pirozzoli_PF2011
	Include "fftw3.f"
	Real(mytype), Dimension(:,:,:), Allocatable :: uu
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxz,outfftxz1,outfftxz2
	Real(mytype), Allocatable, Dimension(:,:) :: ux,ulx,usx,uex,uhx,uelx,spec2d
	Real(mytype), Allocatable, Dimension(:,:) :: am
	Character :: chdum*8,file3d*60
	Real(mytype) :: dum,yam,theta,alpha,beta,uol,ustar,corr0,corr1,corr2
	Integer :: intdum,nsnaps
	Integer :: m,i,j,k,jp,ip,kp,j1,j2
        Integer*8 :: planxz
!
	nsnaps=1
	Allocate(uu(-1:nx+2,-1:ny+2,-1:nz+2))
	Allocate(outfftxz(1:nx/2+1,1:nz),outfftxz1(nx/2+1,1:nz),outfftxz2(nx/2+1,1:nz),ux(1:nx,1:nz), &
ulx(1:nx,1:nz),usx(1:nx,1:nz),uex(1:nx,1:nz),uhx(1:nx,1:nz),uelx(1:nx,1:nz),spec2d(1:nx/2+1,1:nz))
	Allocate(am(1:ny,1:ny))
!
	am(:,:)=0.0e0
!
	Do m=1,nsnaps
  	  Write(chdum,'(I5)') m
!  	  file3d='snap3d_'//Trim(Adjustl(chdum))//'.dat'
	  file3d=restrt_file !'../restrt-start0-stand-run8-out.rst'
  	  If(nsnaps==1) file3d=restrt_file
  	  Open(Unit=14, file=file3d, Status='OLD', Form='UNFormatTED',access='stream')
  	  Read(14) intdum,intdum,intdum,dum,dum,dum,dum,dum,dum,dum,dum
  	  Write(*,'(A,I5,A,I5,F20.10)') 'snap:',m,'  header:',intdum,dum
  	  Read(14) uu(:,:,:) !u
  	  Close(14)
	End Do
!
!	find the y location for am calculation
	yam=15.0e0; ip=nx/2; kp=nz/2
	Do j=1,ny/2
	  If(yplus(j)<yam.And.yplus(j+1)>yam) Then
	    jp=j; Exit
	  End If
	End Do
!
	Do j1=1,jmodu
	  Write(*,*) 'j1=',j1
	  ux(1:nx,1:nz)=uu(1:nx,j1,1:nz)
!	  FFT
	  Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ux,outfftxz,FFTW_ESTIMATE)
	  Call dfftw_execute(planxz)
	  Call dfftw_destroy_plan(planxz)	
	  outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
!	  filtering
	  outfftxz(1,:)=0.0e0
	  outfftxz1(:,:)=outfftxz(:,:)
	  Do i=2,nx/2+1
	    If(length/Real(i-1)<1.0e0) outfftxz1(i,1:nz)=0.0e0
	  End Do
!	  inverse FFT
	  Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz1,ulx,FFTW_ESTIMATE)
	  Call dfftw_execute(planxz)
	  Call dfftw_destroy_plan(planxz)
!
	  Do j2=1,j1
	    ux(1:nx,1:nz)=uu(1:nx,j2,1:nz)
!	    FFT
	    Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ux,outfftxz,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)
	    Call dfftw_destroy_plan(planxz)	
	    outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
!	    filtering
	    outfftxz(1,:)=0.0e0
!	    inverse FFT
	    outfftxz2(:,:)=outfftxz(:,:)
	    Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz2,ux,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)
	    Call dfftw_destroy_plan(planxz)
!	    filtering
	    outfftxz2(:,:)=outfftxz(:,:)
	    Do i=2,nx/2+1
	      If(length/Real(i-1)<1.0e0) outfftxz2(i,1:nz)=0.0e0
	    End Do
!	    inverse FFT
	    Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz2,ulx,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)
	    Call dfftw_destroy_plan(planxz)
!
	    usx(:,:)=ux(:,:)-ulx(:,:)
!	    Hilbert transform
	    uhx(1:nx,1:nz)=usx(1:nx,1:nz)
	    Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,uhx,outfftxz,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)
	    Call dfftw_destroy_plan(planxz)
	    outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
	    theta=-pi/2.0e0
	    Do i=2,nx/2+1
	      outfftxz(i,1:nz)=outfftxz(i,1:nz)*Cmplx(Cos(theta),Sin(theta))
	    End Do
	    Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz,uex,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)  
	    Call dfftw_destroy_plan(planxz)
	    uex(:,:)=Sqrt(uhx(:,:)**2+uex(:,:)**2)
!	    filtering the envelop
	    Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,uex,outfftxz,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz) 
	    Call dfftw_destroy_plan(planxz)
	    outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
	    outfftxz(1,1:nz)=0.0e0
	    Do i=2,nx/2+1
	      If(length/Real(i-1)<1.0e0) outfftxz(i,1:nz)=0.0e0
	    End Do
	    Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz,uelx,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)  
	    Call dfftw_destroy_plan(planxz)
!	    find the AM  
!	    calculate the coviarance
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz(1:nx/2+1,1:nz))*outfftxz1(1:nx/2+1,1:nz)
	    corr0=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
		Sum(spec2d(nx/2+1,1:nz))
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz(1:nx/2+1,1:nz))*outfftxz(1:nx/2+1,1:nz)
	    corr1=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
	  	Sum(spec2d(nx/2+1,1:nz))
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz1(1:nx/2+1,1:nz))*outfftxz1(1:nx/2+1,1:nz)
	    corr2=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
	  	Sum(spec2d(nx/2+1,1:nz))
	    am(j1,j2)=corr0/Sqrt(corr1*corr2)
	  End Do !j2
	End Do !j1
!
!	Open(unit=634,file='uxz-ins.plt')
!	Write(634,'(10A)') ' variables="x","z","u","u<sub>l</sub>","u<sub>s</sub>","u<sub>e</sub>","u<sub>el</sub>"'
!	Write(634,'(A,I5,A,I5,A)') 'zone i=',nx,', j=',nz,',f=point'
!	Do k=1,nz; Do i=1,nx
!	  Write(634,'(7F20.10)') x(i),zc(k),ux(i,k),ulx(i,k),usx(i,k),uex(i,k),uelx(i,k)
!	End Do; End Do
!	Close(634)
!
	Open(unit=635,file='AM_Mathis.plt')
	Write(635,'(10A)') ' variables="yc","AM","y<sup>+</sup>"'
	Do j=1,ny
	  Write(635,'(6F20.10)') yc(j),am(j,j),ycplus(j)
	End Do
	Close(635)
!
	Open(Unit=636,File='AM_2p.plt')
	WRITE(636,'(320A)') ' variables="y<sub>1</sub><sup>+</sup>","y<sub>2</sub><sup>+</sup>","AM"'
 	WRITE(636,*) 'zone i=',jmodu,', j=',jmodu,',f=point'
	Do j1=1,jmodu; Do j2=1,jmodu
	  Write(636,'(5E15.7)') ycplus(j1),ycplus(j2),am(j1,j2)
	End Do; End Do
	Close(636)
!
	Deallocate(uu,ux,ulx,usx,uex,uelx,uhx,outfftxz,outfftxz1,outfftxz2,am,spec2d)	
!
	END SUBROUTINE calc_am2p

!*******************************************************************************
	SUBROUTINE calc_am2p_simplify
! this subroutine calculates the two point amplitude modulation based on Bernardini_Pirozzoli_PF2011
! simplified version of calc_am2p to reduced computing time
! filtered in x direction
	Include "fftw3.f"
	Real(mytype), Dimension(:,:,:), Allocatable :: uu
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxz,outfftxz1,outfftxz2
	Real(mytype), Allocatable, Dimension(:,:) :: ux,ulx,usx,uex,uhx,uelx,spec2d
	Real(mytype), Allocatable, Dimension(:,:) :: am
	Character :: chdum*8,file3d*60
	Real(mytype) :: dum,yam,theta,alpha,beta,uol,ustar,corr0,corr1,corr2, &
	xshift,zshift,kappax,kappaz
	Integer :: intdum,nsnaps
	Integer :: m,i,j,k,jp,ip,kp,j1,j2
        Integer*8 :: planxz
!
	nsnaps=1; jmodu=ny/2
	Allocate(uu(-1:nx+2,-1:ny+2,-1:nz+2))
	Allocate(outfftxz(1:nx/2+1,1:nz),outfftxz1(nx/2+1,1:nz),outfftxz2(nx/2+1,1:nz),ux(1:nx,1:nz), &
ulx(1:nx,1:nz),usx(1:nx,1:nz),uex(1:nx,1:nz),uhx(1:nx,1:nz),uelx(1:nx,1:nz),spec2d(1:nx/2+1,1:nz))
	Allocate(am(1:ny,1:ny))
!
	am(:,:)=0.0e0
!
	Do m=1,nsnaps
  	  Write(chdum,'(I5)') m
!  	  file3d='snap3d_'//Trim(Adjustl(chdum))//'.dat'
	  file3d=restrt_file !'../restrt-start0-stand-run8-out.rst'
  	  If(nsnaps==1) file3d=restrt_file
  	  Open(Unit=14, file=file3d, Status='OLD', Form='UNFormatTED',access='stream')
  	  Read(14) intdum,intdum,intdum,dum,dum,dum,dum,dum,dum,dum,dum
  	  Write(*,'(A,I5,A,I5,F20.10)') 'snap:',m,'  header:',intdum,dum
  	  Read(14) uu(:,:,:) !u
  	  Close(14)
	End Do
!
!	find the y location for am calculation
	yam=15.0e0; ip=nx/2; kp=nz/2
	Do j=1,ny/2
	  If(yplus(j)<yam.And.yplus(j+1)>yam) Then
	    jp=j; Exit
	  End If
	End Do
!
	Do j1=1,jmodu
	  Write(*,*) 'j1=',j1
	  ux(1:nx,1:nz)=uu(1:nx,j1,1:nz)
!	  FFT
	  Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ux,outfftxz,FFTW_ESTIMATE)
	  Call dfftw_execute(planxz)
	  Call dfftw_destroy_plan(planxz)	
	  outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
!	  filtering
	  outfftxz(1,:)=0.0e0
	  outfftxz1(:,:)=outfftxz(:,:)
	  Do i=2,nx/2+1
	    If(length/Real(i-1,mytype)<1.0e0) outfftxz1(i,1:nz)=0.0e0
	  End Do
!	  inverse FFT
	  outfftxz(:,:)=outfftxz1(:,:)
!*	  Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz,ulx,FFTW_ESTIMATE)
!*	  Call dfftw_execute(planxz)
!*	  Call dfftw_destroy_plan(planxz)
!
	  spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz1(1:nx/2+1,1:nz))*outfftxz1(1:nx/2+1,1:nz)
	  corr2=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
	  	Sum(spec2d(nx/2+1,1:nz))
!
	  Do j2=1,jmodu
	    ux(1:nx,1:nz)=uu(1:nx,j2,1:nz)
!	    FFT
	    Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ux,outfftxz,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)
	    Call dfftw_destroy_plan(planxz)	
	    outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
!	    filtering
	    outfftxz(1,:)=0.0e0
	    Do i=2,nx/2+1
	      If(length/Real(i-1,mytype)>=1.0e0) outfftxz(i,1:nz)=0.0e0
	    End Do
	    outfftxz2(:,:)=outfftxz(:,:)
!	    inverse FFT
	    Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz,uhx,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)
	    Call dfftw_destroy_plan(planxz)
!
!	    Hilbert transform
	    theta=-pi/2.0e0
	    Do i=2,nx/2+1
	      outfftxz2(i,1:nz)=outfftxz2(i,1:nz)*Cmplx(Cos(theta),Sin(theta))
	    End Do
	    Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz2,uex,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)  
	    Call dfftw_destroy_plan(planxz)
	    uex(:,:)=Sqrt(uhx(:,:)**2+uex(:,:)**2)
!	    filtering the envelop
	    Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,uex,outfftxz,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz) 
	    Call dfftw_destroy_plan(planxz)
	    outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
	    outfftxz(1,1:nz)=0.0e0
	    Do i=2,nx/2+1
	      If(length/Real(i-1,mytype)<1.0e0) outfftxz(i,1:nz)=0.0e0
	    End Do
!	    phase shift of uelx in streamwise direction considering the streamwise leaning of the super scale structures
	    theta=14.0e0/180.0*pi
	    xshift=-(yc(j1)-yc(j2))/Tan(theta)
            Do i=1,nx/2+1; Do k=1,nz
              kappax=2.0e0*pi/length*Real(i-1,mytype)
              outfftxz(i,k)=outfftxz(i,k)*Cmplx(Cos(kappax*xshift), &
                  Sin(kappax*xshift))
            End Do; End Do
!
	    outfftxz2(:,:)=outfftxz(:,:)
!*	    Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz2,uelx,FFTW_ESTIMATE)
!*	    Call dfftw_execute(planxz)  
!*	    Call dfftw_destroy_plan(planxz)
!	    find the AM  
!	    calculate the coviarance
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz(1:nx/2+1,1:nz))*outfftxz1(1:nx/2+1,1:nz)
	    corr0=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
		Sum(spec2d(nx/2+1,1:nz))
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz(1:nx/2+1,1:nz))*outfftxz(1:nx/2+1,1:nz)
	    corr1=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
	  	Sum(spec2d(nx/2+1,1:nz))
	    am(j1,j2)=corr0 !/Sqrt(corr1*corr2)
!
!Write(*,*) 'corr0=',corr0,corr1,corr2
!corr0=0.0e0
!Do k=1,nz; Do i=1,nx
!  corr0=corr0+uelx(i,k)*ulx(i,k)
!End Do; End Do
!Write(*,*) 'corr1=',corr0/Real(nx*nz,mytype)
!
	  End Do !j2
	End Do !j1
!
	Open(Unit=636,File='AM_2p.plt')
	WRITE(636,'(320A)') ' variables="y<sub>1</sub><sup>+</sup>","y<sub>2</sub><sup>+</sup>","AM"'
 	WRITE(636,*) 'zone i=',jmodu,', j=',jmodu,',f=point'
	Do j1=1,jmodu; Do j2=1,jmodu
	  Write(636,'(5E15.7)') ycplus(j1),ycplus(j2),am(j1,j2)
	End Do; End Do
	Close(636)
!
	Deallocate(uu,ux,ulx,usx,uex,uelx,uhx,outfftxz,outfftxz1,outfftxz2,am,spec2d)	
!
	END SUBROUTINE calc_am2p_simplify

!*******************************************************************************
	SUBROUTINE calc_am2p_simplify_z
! this subroutine calculates the two point amplitude modulation based on Bernardini_Pirozzoli_PF2011
! simplified version of calc_am2p to reduced computing time
! filtered in z direction
	Include "fftw3.f"
	Real(mytype), Dimension(:,:,:), Allocatable :: uu
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxz,outfftxz1,outfftxz2
	Real(mytype), Allocatable, Dimension(:,:) :: ux,ulx,usx,uex,uhx,uelx,spec2d
	Real(mytype), Allocatable, Dimension(:,:) :: am
	Character :: chdum*8,file3d*60
	Real(mytype) :: dum,yam,theta,alpha,beta,uol,ustar,corr0,corr1,corr2
	Integer :: intdum,nsnaps
	Integer :: m,i,j,k,jp,ip,kp,j1,j2,kz
        Integer*8 :: planxz
!
	nsnaps=1; jmodu=ny/2
	Allocate(uu(-1:nx+2,-1:ny+2,-1:nz+2))
	Allocate(outfftxz(1:nx/2+1,1:nz),outfftxz1(nx/2+1,1:nz),outfftxz2(nx/2+1,1:nz),ux(1:nx,1:nz), &
ulx(1:nx,1:nz),usx(1:nx,1:nz),uex(1:nx,1:nz),uhx(1:nx,1:nz),uelx(1:nx,1:nz),spec2d(1:nx/2+1,1:nz))
	Allocate(am(1:ny,1:ny))
!
	am(:,:)=0.0e0
!
	Do m=1,nsnaps
  	  Write(chdum,'(I5)') m
!  	  file3d='snap3d_'//Trim(Adjustl(chdum))//'.dat'
	  file3d=restrt_file !'../restrt-start0-stand-run8-out.rst'
  	  If(nsnaps==1) file3d=restrt_file
  	  Open(Unit=14, file=file3d, Status='OLD', Form='UNFormatTED',access='stream')
  	  Read(14) intdum,intdum,intdum,dum,dum,dum,dum,dum,dum,dum,dum
  	  Write(*,'(A,I5,A,I5,F20.10)') 'snap:',m,'  header:',intdum,dum
  	  Read(14) uu(:,:,:) !u
  	  Close(14)
	End Do
!
!	find the y location for am calculation
	yam=15.0e0; ip=nx/2; kp=nz/2
	Do j=1,ny/2
	  If(yplus(j)<yam.And.yplus(j+1)>yam) Then
	    jp=j; Exit
	  End If
	End Do
!
	Do j1=1,jmodu
	  Write(*,*) 'j1=',j1
	  ux(1:nx,1:nz)=uu(1:nx,j1,1:nz)
!	  FFT
	  Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ux,outfftxz,FFTW_ESTIMATE)
	  Call dfftw_execute(planxz)
	  Call dfftw_destroy_plan(planxz)	
	  outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
!	  filtering
	  outfftxz(1,:)=0.0e0
	  outfftxz1(:,:)=outfftxz(:,:)
	  Do k=2,nz
	    If(k>nz/2+1) kz=nz+2-k
	    If(width/Real(kz-1)<0.5e0) outfftxz1(1:nx/2,k)=0.0e0
	  End Do
!
	  spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz1(1:nx/2+1,1:nz))*outfftxz1(1:nx/2+1,1:nz)
	  corr2=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
	  	Sum(spec2d(nx/2+1,1:nz))
!
	  Do j2=1,jmodu
	    ux(1:nx,1:nz)=uu(1:nx,j2,1:nz)
!	    FFT
	    Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ux,outfftxz,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)
	    Call dfftw_destroy_plan(planxz)	
	    outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
!	    filtering
	    outfftxz(1,:)=0.0e0
	    Do k=2,nz
	      If(k>nz/2+1) kz=nz+2-k
	      If(width/Real(kz-1)>=0.5e0) outfftxz(1:nx/2,k)=0.0e0
	    End Do
	    outfftxz2(:,:)=outfftxz(:,:)
!	    inverse FFT
	    Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz,uhx,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)
	    Call dfftw_destroy_plan(planxz)
!
!	    Hilbert transform
	    theta=-pi/2.0e0
	    Do i=2,nx/2+1
	      outfftxz2(i,1:nz)=outfftxz2(i,1:nz)*Cmplx(Cos(theta),Sin(theta))
	    End Do
	    Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz2,uex,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz)  
	    Call dfftw_destroy_plan(planxz)
	    uex(:,:)=Sqrt(uhx(:,:)**2+uex(:,:)**2)
!	    filtering the envelop
	    Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,uex,outfftxz,FFTW_ESTIMATE)
	    Call dfftw_execute(planxz) 
	    Call dfftw_destroy_plan(planxz)
	    outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
	    outfftxz(1,1:nz)=0.0e0
	    Do k=2,nz
	      If(k>nz/2+1) kz=nz+2-k
	      If(width/Real(kz-1)<0.5e0) outfftxz(1:nx/2,k)=0.0e0
	    End Do
	    outfftxz2(:,:)=outfftxz(:,:)
!
!	    find the AM  
!	    calculate the coviarance
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz(1:nx/2+1,1:nz))*outfftxz1(1:nx/2+1,1:nz)
	    corr0=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
		Sum(spec2d(nx/2+1,1:nz))
	    spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz(1:nx/2+1,1:nz))*outfftxz(1:nx/2+1,1:nz)
	    corr1=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
	  	Sum(spec2d(nx/2+1,1:nz))
	    am(j1,j2)=corr0 !/Sqrt(corr1*corr2)
!
	  End Do !j2
	End Do !j1
!
	Open(Unit=636,File='AM_2p.plt')
	WRITE(636,'(320A)') ' variables="y<sub>1</sub><sup>+</sup>","y<sub>2</sub><sup>+</sup>","AM"'
 	WRITE(636,*) 'zone i=',jmodu,', j=',jmodu,',f=point'
	Do j1=1,jmodu; Do j2=1,jmodu
	  Write(636,'(5E15.7)') ycplus(j1),ycplus(j2),am(j1,j2)
	End Do; End Do
	Close(636)
!
	Deallocate(uu,ux,ulx,usx,uex,uelx,uhx,outfftxz,outfftxz1,outfftxz2,am,spec2d)	
!
	END SUBROUTINE calc_am2p_simplify_z

!*******************************************************************************
	SUBROUTINE skew_decompose
! this subroutine calculates the terms in decomposed skewness of streamwise velocity (Mathis et al 2011)
! u^3=uL^3 + 3*uL^2*uS + 3*uL*uS^2 + uS^3
	Include "fftw3.f"
	Real(mytype), Dimension(:,:,:), Allocatable :: uu
	Complex(mytype), Allocatable, Dimension(:,:) :: outfftxz,outfftxz2
	Real(mytype), Allocatable, Dimension(:,:) :: ux,ulx,usx,spec2d
	Real(mytype), Allocatable, Dimension(:,:) :: sk
	Character :: chdum*8,file3d*60
	Real(mytype) :: dum,yam,urms
	Integer :: intdum,nsnaps
	Integer :: m,i,j,k,jp,ip,kp
        Integer*8 :: planxz
!
	nsnaps=1
	Allocate(uu(-1:nx+2,-1:ny+2,-1:nz+2))
	Allocate(outfftxz(1:nx/2+1,1:nz),outfftxz2(nx/2+1,1:nz),ux(1:nx,1:nz), &
	ulx(1:nx,1:nz),usx(1:nx,1:nz),sk(1:ny,1:5),spec2d(1:nx/2+1,1:nz))
	sk(:,:)=0.0e0
!
	Do m=1,nsnaps
  	  Write(chdum,'(I5)') m
!  	  file3d='snap3d_'//Trim(Adjustl(chdum))//'.dat'
	  file3d=restrt_file !'../restrt-start0-stand-run8-out.rst'
  	  If(nsnaps==1) file3d=restrt_file
  	  Open(Unit=14, file=file3d, Status='OLD', Form='UNFormatTED',access='stream')
  	  Read(14) intdum,intdum,intdum,dum,dum,dum,dum,dum,dum,dum,dum
  	  Write(*,'(A,I5,A,I5,F20.10)') 'snap:',m,'  header:',intdum,dum
  	  Read(14) uu(:,:,:) !u
  	  Close(14)
	End Do
!
!	find the y location for am calculation
	yam=15.0e0; ip=nx/2; kp=nz/2
	Do j=1,ny/2
	  If(yplus(j)<yam.And.yplus(j+1)>yam) Then
	    jp=j; Exit
	  End If
	End Do
!
	Do jp=1,ny
	  Write(*,*) 'j=',jp
	  ux(1:nx,1:nz)=uu(1:nx,jp,1:nz)
!	  FFT
	  Call dfftw_plan_dft_r2c_2d(planxz,nx,nz,ux,outfftxz,FFTW_ESTIMATE)
	  Call dfftw_execute(planxz)
	  Call dfftw_destroy_plan(planxz)	
	  outfftxz(:,:)=outfftxz(:,:)/Real(nx*nz,mytype)
!	  filtering
	  outfftxz(1,1)=0.0e0
!	  variance
	  spec2d(1:nx/2+1,1:nz)=Conjg(outfftxz(1:nx/2+1,1:nz))*outfftxz(1:nx/2+1,1:nz)
	  urms=Sum(spec2d(2:nx/2,1:nz))*2.0e0+Sum(spec2d(1,2:nz))+&
	  	Sum(spec2d(nx/2+1,1:nz))
!
	  outfftxz2(:,:)=outfftxz(:,:)
!	  inverse FFT
	  outfftxz(:,:)=outfftxz2(:,:)
	  Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz,ux,FFTW_ESTIMATE)
	  Call dfftw_execute(planxz)
	  Call dfftw_destroy_plan(planxz)
!
	  Do i=2,nx/2+1
	    If(length/Real(i-1,mytype)<1.0e0) outfftxz2(i,1:nz)=0.0e0
	  End Do
!	  inverse FFT for large scales
	  Call dfftw_plan_dft_c2r_2d(planxz,nx,nz,outfftxz2,ulx,FFTW_ESTIMATE)
	  Call dfftw_execute(planxz)
	  Call dfftw_destroy_plan(planxz)
!	  small scales
	  usx(:,:)=ux(:,:)-ulx(:,:)
!
	  Do i=1,nx; Do k=1,nz
	    sk(jp,1)=sk(jp,1)+ux(i,k)**3 ! u^3
	    sk(jp,2)=sk(jp,2)+ulx(i,k)**3 ! uL^3
	    sk(jp,3)=sk(jp,3)+3.0e0*ulx(i,k)**2*usx(i,k) ! 3*uL^2*uS
	    sk(jp,4)=sk(jp,4)+3.0e0*ulx(i,k)*usx(i,k)**2 ! 3*uL*uS^2 
	    sk(jp,5)=sk(jp,5)+usx(i,k)**3 ! uS^3
	  End Do; End Do
	  sk(jp,1:5)=sk(jp,1:5)/(Sqrt(urms))**3/Real(nx*nz,mytype)
!
	End Do !jp
!
	Open(unit=635,file='skewness_decompose.plt')
	Write(635,'(10A)') ' variables="yc","S","u<sub>L</sub><sup>3</sup>",&
&"3u<sub>L</sub><sup>2</sup>u<sub>S</sub>","3u<sub>L</sub>u<sub>S</sub><sup>2</sup>",&
&"u<sub>S</sub><sup>3</sup>","y<sup>+</sup>"'
	Do j=1,ny
	  Write(635,'(7F20.10)') yc(j),sk(j,1:5),ycplus(j)
	End Do
	Close(635)
!
	Deallocate(uu,ux,ulx,usx,outfftxz,outfftxz2,spec2d)	
!
	END SUBROUTINE skew_decompose

END MODULE calc_post3d

















