	MODULE calc_struct
	USE calc_post3d
	IMPLICIT NONE
	CONTAINS
!*******************************************************************************
! This code performs the calculations on the values that have been read in 
! from the read_post module.
!*******************************************************************************
	SUBROUTINE struct_av
! this subroutine calculate the averaged lam2 structures (considered phase average of controlled cases)
        Integer :: nsnaps,tp,n,l,m,cor_cntx,cor_cor_cntxcntz,cor_cntz,lam_cnt,i3dfield,loop,str
	Character :: file3d*40, chdum*8
	Real(mytype) :: utau_now,rngsch
        Real(mytype), dimension(8) :: dum
 	integer, dimension(3) :: intdum
	Real(mytype), Dimension(ny) :: umn_now,wmn_now,l2_mean,l2_rms
	Real(mytype), Dimension(:,:,:,:), Allocatable :: dummy_3d,vx
	Real(mytype), Dimension(:,:,:), Allocatable :: l2_array, wx,wy,wz
	Real(mytype), Allocatable, Dimension(:) :: cocomax
	Integer, Dimension(:,:,:), Allocatable :: pmin 
	Integer :: i,j,k,i2,j2,k2,ifile,kmin,kmax,jmin,jmax
 	integer, dimension(0:20) :: struc_st
!
!	Initialisations
	Allocate(dummy_3d(-1:nx+2,-1:ny+2,-1:nz+2,4),corr(max_xz,2,3),spec(max_xz,2,3), &
			l2_array(nxt,nyt,nzt),pmin(nxt,nyt,nzt),vx(nxt,nyt,nzt,3),wx(nx,ny,nz))
	simit=sim3d
	cor_cntx = 0; cor_cor_cntxcntz = 0;lam_cnt=0; tp=0
	spec=0.0e0; corr = 0.0e0
	l2_mean=0.0e0; l2_rms=0.0e0
	i3dfield=0
	If (jvis == 0) jvis=ny
        If (ilamq == 1) ivort3d=1
!
	Call open_files3d
!
!	DO n = nstart(simit), nruns(simit)
!
	  If (restrt == 1) then
	    Open(13, File=restrt_file, Status='OLD', Form='UNFormatTED',access='stream')
	    nsnaps=1; tvis=1
	  Else 
	    file3d='OUTPUT_3D_SNAP.dat'
!	    l = n-nstart(simit)+1
	    n=nstart(simit)
	    CALL change_dir(file3d,n)
	    Open(13, File=file3d, Status='OLD', Form='UNFormatTED',access='stream')
!	    CALL header_skip(13)
!	    nsnaps=snap_3d(l,simit)
	    nsnaps=nsmpl_snap_3d
	  End If
!
	  DO m=1, nsnaps
	    write(*,*) 'Reading 3d snapshot',m
	    tp=tp+1
!	    Read(13) intdum(1),intdum(2),intdum(3),dum(1),dum(2),dum(3),dum(4),dum(5),dum(6),dum(7),dum(8)
	    Read(13) intdum(:),dum(:)
	    write(*,*) intdum, dum
	    print *,'reading u from', restrt_file
	    Read(13) dummy_3d(:,:,:,1) !u
            print *,'reading v from', restrt_file
	    Read(13) dummy_3d(:,:,:,2) !v
            print *,'reading w from', restrt_file
	    Read(13) dummy_3d(:,:,:,3) !w
            print *,'reading p from', restrt_file
	    Read(13) dummy_3d(:,:,:,4) !p

	    Write(chdum,'(I5)') m
	    file3d='3dfield-snap'//Trim(Adjustl(chdum))//'.snp'
            If(restrt.Ne.1) Then
              Open(Unit=13+m, file=file3d, Status='NEW', Form='UNFormatTED',access='stream')
	      Write(13+m) intdum(:),dum(:)
	      Write(13+m) dummy_3d(:,:,:,1)
	      Write(13+m) dummy_3d(:,:,:,2)
	      Write(13+m) dummy_3d(:,:,:,3)
	      Write(13+m) dummy_3d(:,:,:,4)
	      Close(13+m)
            End If	
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!only left veloctiy fluctuation
!	    Do j=-1,ny+2
!	      dummy_3d(:,j,:,1)=dummy_3d(:,j,:,1)-flucs(j,1,1)
!	      dummy_3d(:,j,:,3)=dummy_3d(:,j,:,3)-flucs(j,3,1)
!	      dummy_3d(:,j,:,4)=dummy_3d(:,j,:,4)-flucs(j,4,1)
!	    End Do
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	    If(i3dfield.Eq.1) then
	      Open(Unit=303,File='3dfield.plt')
	      WRITE(303,'(120A)') ' variables="x","y","z","u","v","w","p"'
 	      WRITE(303,*) 'zone t="3dfiled", f=point , i=',nx,', j=',jvis,', k=',nz
	      Do k=1,nz; Do j=1,jvis; Do i=1,nx
	        Write(303,'(7f20.10)') xc(i),yc(j),zc(k),0.5*(dummy_3d(i,j,k,1)+dummy_3d(i+1,j,k,1)),0.5*(dummy_3d(i,j,k,2)+dummy_3d(i,j+1,k,2)),0.5*(dummy_3d(i,j,k,3)+dummy_3d(i,j,k+1,3)),dummy_3d(i,j,k,4)
	      End Do; End Do; End Do
	      Close(303)
	    Endif

	    If ((tp == tvis).or.(tvis == 0)) then

	      Call dat_calc(dummy_3d,utau_now,umn_now,wmn_now)
!
	      If (ivort3d==1) Call vort_vis(dummy_3d,wx)
	      If (ilamq==1) then
		Call lambdaQ(dummy_3d,l2_array,pmin,vx,wx)
		Call lambdaQ_mean(l2_array,l2_mean,l2_rms,lam_cnt)
	        Write(chdum,'(I5)') m
	        file3d='lambda2-tmp-bin-snap'//Trim(Adjustl(chdum))//'.dat'
		Open(405,File=file3d, Status='NEW',Form='unformatted',Access='stream')
		Write(405) l2_array(:,:,:)
		!Write(405) pmin(:,:,:)
		Write(405) wx(:,:,:)
		Close(405)
	      End If

	      If(iensem==1) Then
	        Allocate(centres(str_num,50,3),strucdim(str_num,9),strucshift(str_num,2)) !subroutine centre_calc, strc_len<50
	        Allocate(struc_len(str_num),vort_dir(str_num),minx(str_num),xcen(str_num))
	        Allocate(struc_avg_pos(1:xw,1:yw,1:zw,11),struc_avg_neg(1:xw,1:yw,1:zw,11))
	        struc_avg_pos(:,:,:,:)=0.0e0; struc_avg_neg(:,:,:,:)=0.0e0
	        struc_num=0; centres(:,:,:)=0; strucdim(:,:)=0
	        struc_st(:)=0;vort_dir(:)=0
	        cnt_pos=0; cnt_neg=0
	        Open(unit=170,file='structure-info.dat')
	        WRITE(170,'(120A)') ' variables="str","xc","yc","zc","len","<greek>l</greek><sub>2</sub>","<greek>w</greek><sub>x</sub>"'
	        ! get the field averaged wx
	        If(iwxfd.Eq.1) Then
	          Write(*,*) 'get phase averaged wx from instanteous flow field'
	          omega1d_sa(:,4,phase)=0.0e0
	          Do ifile=ifile1,ifile2
	            Write(chdum,'(I5)') (ifile-1)*16+(phase-1)/((tpts-1)/16)+1
	            file3d='lambda2-tmp-bin-snap'//Trim(Adjustl(chdum))//'.dat'
		    Open(405,File=Trim(unit1)//file3d,Form='unformatted',Access='stream')
		    Read(405) l2_array(:,:,:)
		    !Read(405) pmin(:,:,:)
		    Read(405) wx(:,:,:)
	            Close(405)
	            Do i=1,nx; Do k=1,nz	
	              omega1d_sa(:,4,phase)=omega1d_sa(:,4,phase)+wx(i,:,k)
	            End Do; End Do
	          End Do
	          omega1d_sa(:,4,phase)=omega1d_sa(:,4,phase)/Real(nx*nz*(ifile2-ifile1+1))
	        End If	  
		!get wx'      
	        struc_st(ifile1-1)=0
	        Do ifile=ifile1,ifile2
	          Write(chdum,'(I5)') (ifile-1)*16+(phase-1)/((tpts-1)/16)+1
	          file3d='lambda2-tmp-bin-snap'//Trim(Adjustl(chdum))//'.dat'
	          Write(*,*) 'Reading lambda2 field from ->',Trim(unit1)//file3d
	          Write(*,*) '****************************************'
		  Open(405,File=Trim(unit1)//file3d,Form='unformatted',Access='stream')
		  Read(405) l2_array(:,:,:)
	          print *,'reading l2_array from  ', file3d
!		  Read(405) pmin(:,:,:)
!	          print *,'reading pmin from  ', file3d
		  Read(405) wx(:,:,:)
	          Do j=1,ny
	            If(iosci.Eq.1) Then
	              wx(:,j,:)=wx(:,j,:)-omega1d_sa(j,4,phase)
	            Else
	              wx(:,j,:)=wx(:,j,:)-omega_avg(j,4)
	            End If
	          End Do
	          print *,'reading wx from  ', file3d
		  Close(405)
!
		  !find local minimum point for lambda2 value, pmin
		  !goto 1123
	     	  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		  pmin(:,:,:)=0
	          rngsch=strd_thd/retau_avg_sa
	          kmin=-Floor(rngsch/dzs(1)); kmax=Ceiling(rngsch/dzs(1))
	  	  DO k=1,nz; DO j=1,ny; DO i=1,nx
		    If (l2_array(i,j,k).lt.0.0e0) then 
		      pmin(i,j,k)=1
		      jmin=j-1;jmax=j+1
		      Do While((yc(j)-yc(Max(jmin,1))).Lt.rngsch.And.jmin.Ne.1)
		        jmin=Max(jmin-1,1)
		      End Do
		      Do While((yc(Min(jmax,ny))-yc(j)).Lt.rngsch.And.jmax.Ne.ny)
		        jmax=Min(jmax+1,ny)
		      End Do
		      DO k2=kmin,kmax; DO j2=jmin-j,jmax-j
		        If (k2.ne.0.or.j2.ne.0) then
		          If ((j+j2.ge.1).and.(j+j2.le.ny)) then
		  !         If (l2_array(i,j,k).ge.l2_array(i,j+j2,rmod(k+k2,nz)).And. &
		  !            wx(i,j,k)*wx(i,j+j2,rmod(k+k2,nz)).Ge.0.0e0) pmin(i,j,k)=0  
	  	            If (l2_array(i,j,k).ge.l2_array(i,j+j2,rmod(k+k2,nz))) pmin(i,j,k)=0  
		          End If
		        End If
		      ENDDO; ENDDO
		    End If
		  End Do; End Do; End Do
		  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		  !analysis the local minimum pressure and lambda2
!		  Open(Unit=520, File='pminima.plt')
!		  Write(520,*) ' variables="x","y","z"'
!		  Write(520,*) ' zone t="pminima"'
!		  Open(Unit=521, File='l2minima.plt')
!		  Write(521,*) ' variables="x","y","z"'
!		  Write(521,*) ' zone t="l2minima"'
!		  Open(Unit=522, File='l2_wx.plt')
!		  Write(522,*) ' variables="x","y","z","l2","W","wx"'
		  Open(Unit=523, File='wx_phase.plt')
		  Write(523,*) ' variables="y","W","wx"'
		  Write(523,*) ' zone t="wx_phase"'
!		  DO k=1,nz; DO j=1,jvis; DO i=1,nx
!		    If(pmin(i,j,k).Eq.1.And.wx(i,j,k).Gt.0.0e0) Write(521,'(3F20.10)') xc(i),yc(j),zc(k)
!		    If(pmin(i,j,k).Eq.1.And.wx(i,j,k).Le.0.0e0) Write(520,'(3F20.10)') xc(i),yc(j),zc(k)
!		  End Do; End Do; End Do
!		  Write(522,*) ' zone t="bottom", f=point, i=',nx,' ,j=',jvis,' ,k=',nz
!		  DO k=1,nz; DO j=1,jvis; DO i=1,nx
!		    Write(522,'(6F20.10)') xc(i),yc(j),zc(k),l2_array(i,j,k),dummy_3d(i,j,k,3),wx(i,j,k)
!		  End Do; End Do; End Do
!		  Write(522,*) ' zone t="top", f=point, i=',nx,' ,j=',jvis,' ,k=',nz
!		  DO k=1,nz; DO j=ny-jvis+1,ny; DO i=1,nx
!		    Write(522,'(6F20.10)') xc(i),yc(j),zc(k),l2_array(i,j,k),dummy_3d(i,j,k,3),wx(i,j,k)
!		  End Do; End Do; End Do
		  DO j=1,ny
		    If(iosci.Eq.1) Write(523,'(3F20.10)') yc(j),flucs1d_sa(j,3,1,phase),omega1d_sa(j,4,phase)
		  End Do
!		  Close(520);Close(521);Close(522);Close(523)
		  Write(*,*) 'done minima calc'
		  1123 Continue
		  !Stop
!===============================================================================
		  !identify all the satisifed structures
		  vx(:,:,:,:)=0.0e0
                  Call struct_id(l2_array,pmin,wx)
	          struc_st(ifile)=struc_num
	        End Do
		!do raw average of the structures
	        struc_avg_pos(:,:,:,1)=struc_avg_pos(:,:,:,1)/Real(cnt_pos)
	        struc_avg_neg(:,:,:,1)=struc_avg_neg(:,:,:,1)/Real(cnt_neg)
!
		! take correlation for several times to delete low correlation structure
		Do loop=1,iterstr
		  Write(*,*) 'selecting structures at iteration',loop
	    	  struc_avg_pos(:,:,:,3)=0.0e0; struc_avg_neg(:,:,:,3)=0.0e0
	    	  struc_avg_pos(:,:,:,11)=0.0e0; struc_avg_neg(:,:,:,11)=0.0e0
		  strucshift(:,1:2)=0
		  cnt_neg=0; cnt_pos=0
		  Do ifile=ifile1,ifile2
		    Write(*,*) 'struc(ifile)=',struc_st(:)
		    Allocate(cocomax(1:struc_st(ifile)-struc_st(ifile-1)))
	            Write(chdum,'(I5)') (ifile-1)*16+(phase-1)/((tpts-1)/16)+1
	            file3d='lambda2-tmp-bin-snap'//Trim(Adjustl(chdum))//'.dat'
	            Write(*,*) '****************************************'
		    Open(405,File=Trim(unit1)//file3d,Form='unformatted',Access='stream')
		    Read(405) l2_array(:,:,:)
		    !Read(405) pmin(:,:,:)
		    Read(405) wx(:,:,:)
	            Do j=1,ny
	              If(iosci.Eq.1) Then
	                wx(:,j,:)=wx(:,j,:)-omega1d_sa(j,4,phase)
	              Else
	                wx(:,j,:)=wx(:,j,:)-omega_avg(j,4)
	              End If
	            End Do
		    Close(405)
!
	    	    Call struct_corr(l2_array,wx,struc_avg_pos(:,:,:,1),struc_avg_neg(:,:,:,1),struc_st(ifile-1)+1,struc_st(ifile),cocomax(1:struc_st(ifile)-struc_st(ifile-1)))
	      	    j=struc_st(ifile-1); j2=struc_st(ifile-1)
	      	    Do str=struc_st(ifile-1)+1,struc_num
	              If(str.Le.struc_st(ifile)) Then
			If(cocomax(str-struc_st(ifile-1)).Ge.corr_thd*Real(loop)) Then
			   j2=j2+1; Else; Cycle; End If
		      End If
		      j=j+1
	              Do i=1,struc_len(str)
	                struc_len(j)=struc_len(str)
	  	        minx(j)=minx(str)
	     	        xcen(j)=xcen(str)
	                strucdim(j,1:9)=strucdim(str,1:9)
	  	        centres(j,i,1:3)=centres(str,i,1:3)
	                vort_dir(j)=vort_dir(str)
			strucshift(j,1:2)=strucshift(str,1:2)
	              End Do
	            End Do	
	            struc_num=struc_num-struc_st(ifile)+j2
		    Do i=ifile+1,ifile2
		      struc_st(i)=struc_st(i)-struc_st(ifile)+j2
		    End Do
		    struc_st(ifile)=j2
		    Write(*,*) 'struc_num=',struc_num
		    Deallocate(cocomax)	    		    
	            Call struct_average_lam(l2_array,wx,struc_st(ifile-1)+1,struc_st(ifile))
		  End Do
	    	  struc_avg_pos(:,:,:,3)=struc_avg_pos(:,:,:,3)/Real(cnt_pos)
	    	  struc_avg_neg(:,:,:,3)=struc_avg_neg(:,:,:,3)/Real(cnt_neg)
	    	  struc_avg_pos(:,:,:,1)=struc_avg_pos(:,:,:,3)
	    	  struc_avg_neg(:,:,:,1)=struc_avg_neg(:,:,:,3)
	    	  struc_avg_pos(:,:,:,11)=struc_avg_pos(:,:,:,11)/Real(cnt_pos)
	    	  struc_avg_neg(:,:,:,11)=struc_avg_neg(:,:,:,11)/Real(cnt_neg)
		End Do ! loop
!===============================================================================
	        !get the final average
	        struc_avg_pos(:,:,:,2)=0.0e0; struc_avg_pos(:,:,:,4:10)=0.0e0
	        struc_avg_neg(:,:,:,2)=0.0e0; struc_avg_neg(:,:,:,4:10)=0.0e0
	        cnt_pos=0; cnt_neg=0
	        Do ifile=ifile1,ifile2
	          Write(chdum,'(I5)') (ifile-1)*16+(phase-1)/((tpts-1)/16)+1
	          file3d='./3dfield-snap'//Trim(Adjustl(chdum))//'.snp'
	          Open(406, File=Trim(unit1)//file3d, Status='OLD', Form='UNFormatTED',access='stream')
	          Write(*,*) 'Reading 3D field',Trim(unit1)//file3d
	          Write(*,*) '****************************************'
	          Read(406) intdum(:),dum(:)
	          Write(*,*) intdum, dum
	          Read(406) dummy_3d(:,:,:,1) !u
	          print *,'reading u from  ', file3d
	          Read(406) dummy_3d(:,:,:,2) !v
	          print *,'reading v from  ', file3d
	          Read(406) dummy_3d(:,:,:,3) !w
	          print *,'reading w from  ', file3d
	          Read(406) dummy_3d(:,:,:,4) !p
	          print *,'reading p from  ', file3d
	          Close(406)
!
                  Call struct_average_vel(dummy_3d,struc_st(ifile-1)+1,struc_st(ifile))
	        End Do	
	        struc_avg_pos(:,:,:,2)=struc_avg_pos(:,:,:,2)/real(cnt_pos)
	        struc_avg_neg(:,:,:,2)=struc_avg_neg(:,:,:,2)/real(cnt_neg)        
	        struc_avg_pos(:,:,:,4:10)=struc_avg_pos(:,:,:,4:10)/real(cnt_pos)
	        struc_avg_neg(:,:,:,4:10)=struc_avg_neg(:,:,:,4:10)/real(cnt_neg)
	      ! write averaged structure
                Open(407, File='averaged-struc-pos.plt', status='UNKNOWN')
                WRITE(407,*) ' title="averaged-struc.plt"'
                WRITE(407,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","<greek>l</greek><sub>2</sub><sub>crd</sub>","C","<greek>l</greek><sub>2</sub><sub>ref</sub>","Ua","Wa","Pa","u","v","w","p","<greek>w</greek><sub>x</sub>"'
                write(407,*) ' zone t="postive",i=',xw,', j=',yw,', k=',zw,',f=point'
                do k=-zw/2,zw-zw/2-1; do j=1,yw; do i=-xw/2,xw-xw/2-1
	          i2=i+xw/2+1; j2=j; k2=k+zw/2+1
                  write(407,'(14E15.7)') real(i)*dx(1)*retau_avg_sa,yc(j)*retau_avg_sa,real(k)*dz(1)*retau_avg_sa,struc_avg_pos(i2,j2,k2,1)*nu**2/utau_avg_sa**4,struc_avg_pos(i2,j2,k2,2),struc_avg_pos(i2,j2,k2,3)*nu**2/utau_avg_sa**4,(struc_avg_pos(i2,j2,k2,l),l=4,11)
                end do ;end do ;end do
                Open(408, File='averaged-struc-neg.plt', status='UNKNOWN')
                WRITE(408,*) ' title="averaged-struc.plt"'
                WRITE(408,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","<greek>l</greek><sub>2</sub><sub>crd</sub>","C","<greek>l</greek><sub>2</sub><sub>ref</sub>","Ua","Wa","Pa","u","v","w","p","<greek>w</greek><sub>x</sub>"'
                write(408,*) ' zone t="negative",i=',xw,', j=',yw,', k=',zw,',f=point'
                do k=-zw/2,zw-zw/2-1; do j=1,yw; do i=-xw/2,xw-xw/2-1
	          i2=i+xw/2+1; j2=j; k2=k+zw/2+1
                  write(408,'(14E15.7)') real(i)*dx(1)*retau_avg_sa,yc(j)*retau_avg_sa,real(k)*dz(1)*retau_avg_sa,struc_avg_neg(i2,j2,k2,1)*nu**2/utau_avg_sa**4,struc_avg_neg(i2,j2,k2,2),struc_avg_neg(i2,j2,k2,3)*nu**2/utau_avg_sa**4,(struc_avg_neg(i2,j2,k2,l),l=4,11)
                end do ;end do ;end do
!	        Deallocate (centres,strucdim,strucshift,struc_len,vort_dir,minx,xcen,struc_avg_pos,struc_avg_neg)
	        Close(407);Close(408);Close(170)
	      End If
!-----------------------------------------------------------------------------------
	      if(ivort3d==1.and.0) then
!write middle data file for wx, wy, wz and lamda2, so that can calculate the correlation between -lamda2 and others
		open(unit=306,file='vorticity-lambda2.plt')
	        WRITE(306,'(120A)') ' variables="x","y","z","wx","wy","wz","p","lam2"'
 	        WRITE(306,*) 'zone DATAPACKING=BLOCK, i=',nx,', j=',jvis,', k=',nz
	        Write(306,'(10f20.10)') (((xc(i),i=1,nx),j=1,jvis),k=1,nz)
	        Write(306,'(10f20.10)') (((yc(j),i=1,nx),j=1,jvis),k=1,nz)
	        Write(306,'(10f20.10)') (((zc(k),i=1,nx),j=1,jvis),k=1,nz)
	        Write(306,'(10f20.10)') (((wx(i,j,k)*nu/utau_avg_sa**2,i=1,nx),j=1,jvis),k=1,nz)
	        Write(306,'(10f20.10)') (((wy(i,j,k)*nu/utau_avg_sa**2,i=1,nx),j=1,jvis),k=1,nz)
	        Write(306,'(10f20.10)') (((wz(i,j,k)*nu/utau_avg_sa**2,i=1,nx),j=1,jvis),k=1,nz)
	        Write(306,'(10f20.10)') (((dummy_3d(i,j,k,4)/utau_avg_sa**2/rho,i=1,nx),j=1,jvis),k=1,nz)
	        Write(306,'(10f20.10)') (((-1.0*l2_array(i,j,k)*nu**2/utau_avg_sa**4,i=1,nx),j=1,jvis),k=1,nz)		
		close(306)
	      endif
	    End If
	  END DO
	  Close(13)
!    	END DO
!
	END SUBROUTINE struct_av

!*******************************************************************************
	SUBROUTINE struct_id(l2_array,pmin,wx)
!	This subroutine calculates the correlations
	Real(mytype), Dimension(nxt,nyt,nzt), Intent(in) :: l2_array
	Integer, Dimension(nxt,nyt,nzt), Intent(inout) :: pmin 
!	Real(mytype), Dimension(nxt,nyt,nzt,3), Intent(in) :: vx
	Real(mytype), Dimension(nx,ny,nz), Intent(in) :: wx

        Integer :: i,j,k,m,jin,kin,i1,j1,k1,struc_end,i2,j2,k2,struc_start,neg, &
                   jstr,jend,str,iext,is,ks,ishift,kshift,iloop,&
                   cnt_n,cnt_p, itmp,imax,idiscd
        Real(mytype), dimension(3) :: lenl,maxlen
        Integer, dimension(3) :: maxall
        Integer, Allocatable, dimension(:) :: sort
	Real(mytype), Allocatable, Dimension(:) :: cocomax
        Real(mytype) :: dysum,dzsum,xs,zs,wx_avg,l2_avg, tmp
	Real(mytype), Dimension(:,:,:,:), allocatable :: vx_dum
	Real(mytype), Dimension(1:xw,1:yw,1:zw) :: struc_pos,struc_neg
	Character*10 :: chord,chdir

	! centres stores i,j,k for 100 structures max lenth=50
        allocate(vx_dum(nxt,nyt,nzt,3))
	struc_start=struc_num+1
	! find structures in the lower channel
	if(ibot.eq.1) then
          do i=1,nx; do j=1,jvis; do k=1,nz          
            if(pmin(i,j,k)==1) then
              struc_num=struc_num+1
              struc_len(struc_num)=0
!              write(*,*) 'Structure ',struc_num
              !forwards
              if (wx(i,j,k) .lt. 0.0e0) then; vort_dir(struc_num)=-1
              else; vort_dir(struc_num)=1; end if
              i2=i;j2=j;k2=k
              call centre_calc(l2_array,pmin,vx_dum,wx,i2,j2,k2,1,struc_len(struc_num))
              !backwards
              i2=i;j2=j;k2=k
              call centre_calc(l2_array,pmin,vx_dum,wx,i2,j2,k2,-1,struc_len(struc_num))
!              write(*,*) 'structure',struc_num,'length=', struc_len(struc_num)
!              write(*,*) '***********************************'
	      !discard low lambda2 value structu!res
	      wx_avg=0.0e0; l2_avg=0.0e0
	      do j1=1,struc_len(struc_num)
!	        Write(*,*) centres(struc_num,j1,1),centres(struc_num,j1,2),centres(struc_num,j1,3)
	        wx_avg=wx_avg+wx(centres(struc_num,j1,1),centres(struc_num,j1,2),centres(struc_num,j1,3))
	        l2_avg=l2_avg+l2_array(centres(struc_num,j1,1),centres(struc_num,j1,2),centres(struc_num,j1,3))
	      end do
	      wx_avg=wx_avg/Real(struc_len(struc_num))
	      l2_avg=l2_avg/Real(struc_len(struc_num))
!              if (struc_len(struc_num).lt.15.or.l2_avg.ge.(-1.0e0)) then
              if (struc_len(struc_num).lt.Ceiling(strl_thd/retau_avg_sa/dx(1)).Or.struc_len(struc_num).Gt.50) then
                struc_num=struc_num-1 
	      else
                write(*,*) 'Structure ',struc_num,' with voricity ',vort_dir(struc_num),' wx_avg=',wx_avg,' l2_avg=',l2_avg
                do str=1,struc_len(struc_num)
!Write(*,*) str,wx(centres(struc_num,str,1),centres(struc_num,str,2),centres(struc_num,str,3)),l2_array(centres(struc_num,str,1),centres(struc_num,str,2),centres(struc_num,str,3))
                  pmin(centres(struc_num,str,1),centres(struc_num,str,2),centres(struc_num,str,3))=0
                end do
              end if
            end if 
              if (struc_num==str_num) goto 15 
!              if (struc_num==10) goto 15
          end do; end do; end do
        end if
	! =======================================
	! find structures in the upper channel
	if(itop.eq.1) then   
          do i=1,nx; do j=ny,ny-jvis+1,-1; do k=1,nz       
            if(pmin(i,j,k)==1) then
              struc_num=struc_num+1
              struc_len(struc_num)=0
!              write(*,*) 'Structure ',struc_num
              !forwards
              if (wx(i,j,k) .lt. 0.0e0) then; vort_dir(struc_num)=-1
              else; vort_dir(struc_num)=1; end if
              i2=i;j2=j;k2=k
              call centre_calc(l2_array,pmin,vx_dum,wx,i2,j2,k2,1,struc_len(struc_num))
              !backwards
              i2=i;j2=j;k2=k
              call centre_calc(l2_array,pmin,vx_dum,wx,i2,j2,k2,-1,struc_len(struc_num))
!              write(*,*) 'structure',struc_num,'length=', struc_len(struc_num)
!              write(*,*) '***********************************'
	      !discard low lambda2 value structures
	      wx_avg=0.0e0; l2_avg=0.0e0
	      do j1=1,struc_len(struc_num)
!	        Write(*,*) centres(struc_num,j1,1),centres(struc_num,j1,2),centres(struc_num,j1,3)
	        wx_avg=wx_avg+wx(centres(struc_num,j1,1),centres(struc_num,j1,2),centres(struc_num,j1,3))
	        l2_avg=l2_avg+l2_array(centres(struc_num,j1,1),centres(struc_num,j1,2),centres(struc_num,j1,3))
	      end do
	      wx_avg=wx_avg/Real(struc_len(struc_num))
	      l2_avg=l2_avg/Real(struc_len(struc_num))
!              if (struc_len(struc_num).lt.15.or.l2_avg.ge.(-1.0e0)) then
              if (struc_len(struc_num).lt.Ceiling(strl_thd/retau_avg_sa/dx(1)).Or.struc_len(struc_num).Gt.50) then
                struc_num=struc_num-1 
	      else
                write(*,*) 'Structure ',struc_num,' with voricity ',vort_dir(struc_num),' wx_avg=',wx_avg,' l2_avg=',l2_avg
                do str=1,struc_len(struc_num)
!Write(*,*) str,wx(centres(struc_num,str,1),centres(struc_num,str,2),centres(struc_num,str,3)),l2_array(centres(struc_num,str,1),centres(struc_num,str,2),centres(struc_num,str,3))
                  pmin(centres(struc_num,str,1),centres(struc_num,str,2),centres(struc_num,str,3))=0
                end do
              end if
            end if 
              if (struc_num==str_num) goto 15 
!              if (struc_num==10) goto 15
          end do; end do; end do
        end if

15        write(*,*) 'done ', struc_num, ' structures  >>>><<<<<<  itop=',itop, '  ibot=',ibot
        maxlen=0.0e0; maxall=0
        do i=struc_start,struc_num
          call middle_calc(centres(i,1:struc_len(i),1:3),lenl,struc_len(i),maxall,minx(i),strucdim(i,1:9))
          if(lenl(1).gt.maxlen(1)) maxlen(1)=lenl(1)
          if(lenl(2).gt.maxlen(2)) maxlen(2)=lenl(2)
          if(lenl(3).gt.maxlen(3)) maxlen(3)=lenl(3)
!          write(*,*) 'length=', struc_len(i),lenl(1),lenl(2),lenl(3)          
!          write(*,*) '***********************************'
	  do j=1,struc_len(i)
	    if(strucdim(i,2).eq.centres(i,j,1)) then
	      xcen(i)=j
	      exit
	    end if
	  end do
	end do
!        write(*,*) 'max length=', maxlen(1),maxlen(2),maxlen(3)
!        write(*,*) 'max pts=', maxall(1)/struc_num,maxall(2)/struc_num,maxall(3)/struc_num
        maxall(1)=xw  !maxall(1)/struc_num+20
        maxall(2)=yw  !maxall(2)/struc_num+5
        maxall(3)=zw  !2*maxall(3)/struc_num+4
	iext=0
!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
	cnt_n=0; cnt_p=0
	Do str=struc_start,struc_num
	  If(vort_dir(str).Eq.1) cnt_p=cnt_p+1
	  If(vort_dir(str).Eq.-1) cnt_n=cnt_n+1
	End Do
	! filt some structures
	If(ifilt.Eq.1.And.(cnt_n.Gt.(cnt_stop*(ibot+itop)).Or.cnt_p.Gt.(cnt_stop*(ibot+itop)))) Then
	  Allocate(cocomax(1:(struc_num-struc_start+1)),sort(1:struc_num-struc_start+1))
	  Do iloop=1,1
	    ! average the structure
	    struc_pos(:,:,:)=0.0e0; struc_neg(:,:,:)=0.0e0
	    Do str=struc_start,struc_num
	      do k=-iext-maxall(3)/2,iext+maxall(3)-maxall(3)/2-1; do j=1,maxall(2)+2*iext; do i=-iext-maxall(1)/2,iext+maxall(1)-maxall(1)/2-1
	        i1=i+iext+maxall(1)/2+1; j1=j; k1=k+iext+maxall(3)/2+1
	        i2=i+centres(str,xcen(str),1); j2=j; k2=k+centres(str,xcen(str),3)
	        if(strucdim(str,5).gt.ny/2) then
	          j2=ny-j+1; k2=-k+centres(str,xcen(str),3)
	        end if
	        i2=rmod(i2+nx,nx); k2=rmod(k2+nz,nz)
	        if(vort_dir(str).eq.1) then
	          struc_pos(i1,j1,k1)=struc_pos(i1,j1,k1)+l2_array(i2,j2,k2)  !lambda2
	        end if
	        if(vort_dir(str).eq.-1) then
	          struc_neg(i1,j1,k1)=struc_neg(i1,j1,k1)+l2_array(i2,j2,k2)  !lambda2
	        end if
	      end do; end do; end do
	    End Do
	    struc_pos(:,:,:)=struc_pos(:,:,:)/Real(cnt_p)
	    struc_neg(:,:,:)=struc_neg(:,:,:)/Real(cnt_n)
	    ! take the correlation
	    Call struct_corr(l2_array,wx,struc_pos,struc_neg,struc_start,struc_num,cocomax)
	    ! find the structure with higher correlation
!
	    !check negative structure
	    If(cnt_n.Gt.(cnt_stop*(ibot+itop))) Then
	      imax=cnt_n
	      i=0
	      Do str=struc_start,struc_num
	        If(vort_dir(str).Eq.-1) Then
	          i=i+1
	          sort(i)=str
	        End If	      
	      End Do
	      ! sort structure as decrecent order
	      Do j=imax-1,1,-1
	        Do i=1,j
	          If(cocomax(sort(i)).Lt.cocomax(sort(i+1))) Then
	            itmp=sort(i); sort(i)=sort(i+1); sort(i+1)=itmp
	          End If	      
	        End Do
	        !Write(*,*) j,'sort(:)', sort(:)
	        !Write(*,*) cocomax(sort(:))
	        !Write(*,*) vort_dir(sort(:))
	      End Do

	      j=struc_start-1
	      Do str=struc_start,struc_num
	        idiscd=0
	        Do k=cnt_stop*(ibot+itop)+1,imax
                  If(str.Eq.sort(k)) Then
	            idiscd=1
	            Exit
	          End If
	        End Do
	        If(idiscd.Eq.0) Then
	          j=j+1
	          Do i=1,struc_len(str)
	            struc_len(j)=struc_len(str)
	  	    minx(j)=minx(str)
	     	    xcen(j)=xcen(str)
	            strucdim(j,1:9)=strucdim(str,1:9)
	  	    centres(j,i,1:3)=centres(str,i,1:3)
	            vort_dir(j)=vort_dir(str)
	          End Do
	        Else
!	          Write(*,*) 'throwing low correlation structure ->', str, vort_dir(str)
	        End If
	      End Do	
	      struc_num=j      
	    End If
	    !check positive structure
	    If(cnt_p.Gt.(cnt_stop*(ibot+itop))) Then
	      imax=cnt_p
	      i=0
	      Do str=struc_start,struc_num
	        If(vort_dir(str).Eq.1) Then
	          i=i+1
	          sort(i)=str
	        End If	      
	      End Do
	      ! sort structure as decrecent order
	      Do j=imax-1,1,-1
	        Do i=1,j
	          If(cocomax(sort(i)).Lt.cocomax(sort(i+1))) Then
	            itmp=sort(i); sort(i)=sort(i+1); sort(i+1)=itmp
	          End If	      
	        End Do
	        !Write(*,*) j,'sort(:)', sort(:)
	        !Write(*,*) cocomax(sort(:))
	        !Write(*,*) vort_dir(sort(:))
	      End Do

	      j=struc_start-1
	      Do str=struc_start,struc_num
	        idiscd=0
	        Do k=cnt_stop*(ibot+itop)+1,imax
                  If(str.Eq.sort(k)) Then
	            idiscd=1
	            Exit
	          End If
	        End Do
	        If(idiscd.Eq.0) Then
	          j=j+1
	          Do i=1,struc_len(str)
	            struc_len(j)=struc_len(str)
	  	    minx(j)=minx(str)
	     	    xcen(j)=xcen(str)
	            strucdim(j,1:9)=strucdim(str,1:9)
	  	    centres(j,i,1:3)=centres(str,i,1:3)
	            vort_dir(j)=vort_dir(str)
	          End Do
	        Else
!	          Write(*,*) 'throwing low correlation structure ->', str, vort_dir(str)
	        End If
	      End Do	
	      struc_num=j 
	    End If
	  End Do
	  Deallocate(cocomax, sort)
	End If
!->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->->
        do i=struc_start,struc_num
	  Write(chord,'(I5)') i
	  If(vort_dir(i).Eq.1) Then; chdir='pos'; Else; chdir='neg'; End If
!	  if(iwrite.eq.1.And.vort_dir(i).Eq.1) then
	  if(iwrite.eq.1) then
            Open(26, File=trim(adjustl(chdir))//'_'//trim(adjustl(chord))//'-struct-ins.plt', status='UNKNOWN')
            WRITE(26,*) ' title="v1-ins.plt"'
            WRITE(26,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","v1(1)","v1(2)","v1(3)"'
	    WRITE(26,*) ' zone t=" ',i ,' " ' 
	  end if
          do j=1,struc_len(i)
            i2=centres(i,j,1); j2=centres(i,j,2); k2=centres(i,j,3)
!            write(*,*) i2,j2,k2
	    xs=xc(i2)-xc(centres(i,xcen(i),1))
	    if(strucdim(i,1).gt.strucdim(i,3)) then
	      if(xs.gt.xc(nx/2)) xs=xs-xc(nx)
	      if(xs.lt.-xc(nx/2)) xs=xs+xc(nx)
	    end if 
	    zs=zc(k2)-zc(centres(i,xcen(i),3))
	    if(strucdim(i,7).gt.strucdim(i,9)) then
	      if(zs.gt.zc(nz/2)) zs=zs-zc(nz)
	      if(zs.lt.-zc(nz/2)) zs=zs+zc(nz)
	    end if 
            If (iwrite.eq.1) Then
!	     write(26,'(6E15.7)') xs*retau_avg_sa,yc(j2)*retau_avg_sa,zs*retau_avg_sa,vx_dum(i2,j2,k2,1),vx_dum(i2,j2,k2,2),vx_dum(i2,j2,k2,3)
	     write(26,'(6E15.7)') xs,yc(j2),zs,vx_dum(i2,j2,k2,1),vx_dum(i2,j2,k2,2),vx_dum(i2,j2,k2,3)
!	     write(26,'(6E15.7)') xc(i2),yc(j2),zc(k2),vx_dum(i2,j2,k2,1),vx_dum(i2,j2,k2,2),vx_dum(i2,j2,k2,3)
	    End If
          end do
	  If (iwrite.eq.1) Close(26)      
        end do    
!
        do str=struc_start,struc_num
          if (iwrite.eq.1) then
            Write(chord,'(I5)') str
	    If(vort_dir(str).Eq.1) Then; chdir='pos'; Else; chdir='neg'; End If  
            Open(29, File=trim(adjustl(chdir))//'_'//trim(adjustl(chord))//'-indiv-struc.plt', status='UNKNOWN')
            WRITE(29,*) ' title="indiv-struc.plt"'
            WRITE(29,'(120A)') ' variables="x<sup>+</sup>","y<sup>+</sup>","z<sup>+</sup>","<greek>l</greek><sub>2</sub>"'
!            write(*,*) 'writing ',str
 	    write(29,*) ' zone t="',str,'",i=',maxall(1)+2*iext,', j=',maxall(2)+2*iext,', k=',maxall(3)+2*iext,',f=point'
!	    write(*,*) 'dimension=',-iext-maxall(1)/2,iext+maxall(1)-maxall(1)/2,1,strucdim(str,5)+iext,-iext-maxall(3)/2,iext+maxall(3)-maxall(3)/2
            do k=-iext-maxall(3)/2,iext+maxall(3)-maxall(3)/2-1; do j=1,maxall(2)+2*iext; do i=-iext-maxall(1)/2,iext+maxall(1)-maxall(1)/2-1
              i2=i+centres(str,xcen(str),1); j2=j; k2=k+centres(str,xcen(str),3)
	      if(strucdim(str,5).gt.ny/2) j2=ny-j+1
	      i2=rmod(i2+nx,nx); k2=rmod(k2+nz,nz)
!              write(29,'(5E15.7)') real(i)*dx(1)*retau_avg_sa,yc(j2)*retau_avg_sa,real(k)*dz(1)*retau_avg_sa,l2_array(i2,j2,k2)*nu**2/utau_avg_sa**4
              write(29,'(5E15.7)') real(i)*dx(1),yc(j2),real(k)*dz(1),l2_array(i2,j2,k2) !*nu**2/utau_avg_sa**4
            end do ;end do ;end do 
	    Close(29)
          end if
	  ! write structure information
	  Write(170,'(5I8,2E15.7)') str,centres(str,xcen(str),1),centres(str,xcen(str),2),centres(str,xcen(str),3),struc_len(str),l2_array(centres(str,xcen(str),1),centres(str,xcen(str),2),centres(str,xcen(str),3)),wx(centres(str,xcen(str),1),centres(str,xcen(str),2),centres(str,xcen(str),3))
        end do
! get the averaged structure
	write(*,*) 'rude averaged structure'
	do str=struc_start,struc_num
	  if(vort_dir(str).eq.1) cnt_pos=cnt_pos+1
	  if(vort_dir(str).eq.-1) cnt_neg=cnt_neg+1
	  do j=1,struc_len(str)
!	write(*,*) j
!	write(*,*) "**************************"
	    i2=centres(str,j,1); j2=centres(str,j,2); k2=centres(str,j,3)
!        write(*,*) i2,j2,k2
            if(strucdim(str,1).gt.strucdim(str,3)) then
	      if(centres(str,xcen(str),1).gt.nx/2) then
                if(i2.lt.nx/2) i2=i2+nx
	      else
	        if(i2.gt.nx/2) i2=i2-nx
	      end if
            end if
            if(strucdim(str,7).gt.strucdim(str,9)) then
	      if(centres(str,xcen(str),3).gt.nz/2) then
                if(k2.lt.nz/2) k2=k2+nz
	      else
	        if(k2.gt.nz/2) k2=k2-nz
	      end if
            end if
	    i2=i2-centres(str,xcen(str),1); k2=k2-centres(str,xcen(str),3)
!	    if(itop.eq.1.and.strucdim(str,5).gt.ny/2) then
	    if(strucdim(str,5).gt.ny/2) then
	      j2=ny-j2+1; k2=-k2
	    end if
	    i2=i2+iext+maxall(1)/2+1; k2=k2+iext+maxall(3)/2+1
	    if(i2.ge.1.and.i2.le.maxall(1)+2*iext.and.k2.ge.1.and.k2.le.maxall(3)+2*iext) then
	      if(vort_dir(str).eq.1) then
	        struc_avg_pos(i2,j2,k2,2)=struc_avg_pos(i2,j2,k2,2)+1.0e0
	      else
	        struc_avg_neg(i2,j2,k2,2)=struc_avg_neg(i2,j2,k2,2)+1.0e0
	      end if
	    end if
	  end do
	  do k=-iext-maxall(3)/2,iext+maxall(3)-maxall(3)/2-1; do j=1,maxall(2)+2*iext; do i=-iext-maxall(1)/2,iext+maxall(1)-maxall(1)/2-1
	    i1=i+iext+maxall(1)/2+1; j1=j; k1=k+iext+maxall(3)/2+1
	    i2=i+centres(str,xcen(str),1); j2=j; k2=k+centres(str,xcen(str),3)
!	    if(itop.eq.1.and.strucdim(str,5).gt.ny/2) then
	    if(strucdim(str,5).gt.ny/2) then
	      j2=ny-j+1; k2=-k+centres(str,xcen(str),3)
	    end if
	    i2=rmod(i2+nx,nx); k2=rmod(k2+nz,nz)
	    if(vort_dir(str).eq.1) then
	      struc_avg_pos(i1,j1,k1,1)=struc_avg_pos(i1,j1,k1,1)+l2_array(i2,j2,k2)  !lambda2
	    end if
	    if(vort_dir(str).eq.-1) then
	      struc_avg_neg(i1,j1,k1,1)=struc_avg_neg(i1,j1,k1,1)+l2_array(i2,j2,k2)  !lambda2
	    end if
	  end do; end do; end do
	end do
	Write(*,*) 'cnt_pos=',cnt_pos,'cnt_neg=',cnt_neg
!Write(*,*) 'stopped at the end of struct_id'
!STOP
	END SUBROUTINE struct_id
!*******************************************************************************
	SUBROUTINE struct_corr(l2_array,wx,struc_pos,struc_neg,struc_start,struc_end,cocomax)
!	This subroutine calculates the correlations
	Real(mytype), Dimension(nxt,nyt,nzt), Intent(in) :: l2_array
	Real(mytype), Dimension(nx,ny,nz), Intent(in) :: wx
	Real(mytype), Dimension(:), Intent(out) :: cocomax
	Integer, Intent(in) :: struc_start,struc_end
	Real(mytype), Dimension(1:xw,1:yw,1:zw), Intent(in) :: struc_pos,struc_neg

        Integer :: i,j,k,m,jin,kin,i1,j1,k1,i2,j2,k2,neg, &
                   str,iext,is,ks,ishift,kshift
        Real(mytype), dimension(3) :: lenl,maxlen
        Integer, dimension(3) :: maxall
	Integer :: coragx,coragy,coragz
        Real(mytype) :: dysum,dzsum,xs,zs,coco,coco1_pos,coco1_neg,coco2
	Character*10 :: chord,chdir

	! take the correlation
	write(*,*) 'take correlation'
	! display window for structures
        maxall(1)=xw 
        maxall(2)=yw 
        maxall(3)=zw 
	! range to take cross correlation, 150 wall units in streamwise and 50 wall units in spanwise direction
        coragx=Int(150.0e0/retau_avg_sa/dx(1))
        coragy=yw 
        coragz=Int(50.0e0/retau_avg_sa/dz(1))
	iext=0
	!shifted 50 wall units in streamwise and 30 wall units in spanwise
	ishift=Int(50.0e0/retau_avg_sa/dx(1)); 
	kshift=Int(30.0e0/retau_avg_sa/dz(1))  
	Write(*,*) 'correlation at range:', 'x+=',(coragx)*dx(1)*retau_avg_sa,'/  y+=',yc(coragy)*retau_avg_sa,'/  z+=', (coragz)*dz(1)*retau_avg_sa
	cocomax(:)=0.0e0
	coco1_pos=0.0e0; coco1_neg=0.0e0
	do k=-iext-coragz/2,iext+coragz-coragz/2-1; do j=1,coragy+2*iext; do i=-iext-coragx/2,iext+coragx-coragx/2-1
	  i1=i+iext+maxall(1)/2+1; j1=j; k1=k+iext+maxall(3)/2+1
	  coco1_pos=coco1_pos+struc_pos(i1,j1,k1)**2
	end do; end do; end do
	coco1_pos=coco1_pos/real((coragx+2*iext)*(coragy+2*iext)*(coragz+2*iext))
	do k=-iext-coragz/2,iext+coragz-coragz/2-1; do j=1,coragy+2*iext; do i=-iext-coragx/2,iext+coragx-coragx/2-1
	  i1=i+iext+maxall(1)/2+1; j1=j; k1=k+iext+maxall(3)/2+1
	  coco1_neg=coco1_neg+struc_neg(i1,j1,k1)**2
	end do; end do; end do
	coco1_neg=coco1_neg/real((coragx+2*iext)*(coragy+2*iext)*(coragz+2*iext))
	do str=struc_start,struc_end
	  if(vort_dir(str).eq.1) then ! positive structure
	    do ks=-kshift,kshift; do is=-ishift,ishift
	      coco2=0.0e0; coco=0.0e0
	      do k=-iext-coragz/2,iext+coragz-coragz/2-1; do j=1,coragy+2*iext; do i=-iext-coragx/2,iext+coragx-coragx/2-1
	        i1=i+iext+maxall(1)/2+1; j1=j; k1=k+iext+maxall(3)/2+1
	        i2=i+centres(str,xcen(str),1)+is; j2=j; k2=k+centres(str,xcen(str),3)+ks
!	        if(itop.eq.1.and.strucdim(str,5).gt.ny/2) then
	        if(strucdim(str,5).gt.ny/2) then
	          j2=ny-j+1; k2=-k+centres(str,xcen(str),3)-ks
	        end if
	        i2=rmod(i2+nx,nx); k2=rmod(k2+nz,nz)
	        coco2=coco2+l2_array(i2,j2,k2)**2
	        coco=coco+struc_pos(i1,j1,k1)*l2_array(i2,j2,k2)
	      end do; end do; end do
	      coco2=coco2/real((coragx+2*iext)*(coragy+2*iext)*(coragz+2*iext))
	      coco=coco/real((coragx+2*iext)*(coragy+2*iext)*(coragz+2*iext))
	      coco=coco/sqrt(coco1_pos)/sqrt(coco2)
	      if(coco.ge.cocomax(str-struc_start+1)) then
	        cocomax(str-struc_start+1)=coco
	        strucshift(str,1)=is
	        strucshift(str,2)=ks
	      end if
	    end do; end do
!	    Write(*,*) 'cross-correlation_max=',cocomax(str-struc_start+1),'ishift=',strucshift(str,1),'kshift=',strucshift(str,2)
	    Write(*,*) 'cross-correlation_max=',cocomax(str-struc_start+1),strucshift(str,1),strucshift(str,2)
	  end if
	  if(vort_dir(str).eq.-1) then  ! negative structure
	    do ks=-kshift,kshift; do is=-ishift,ishift
	      coco2=0.0e0; coco=0.0e0
	      do k=-iext-coragz/2,iext+coragz-coragz/2-1; do j=1,coragy+2*iext; do i=-iext-coragx/2,iext+coragx-coragx/2-1
	        i1=i+iext+maxall(1)/2+1; j1=j; k1=k+iext+maxall(3)/2+1
	        i2=i+centres(str,xcen(str),1)+is; j2=j; k2=k+centres(str,xcen(str),3)+ks
!	        if(itop.eq.1.and.strucdim(str,5).gt.ny/2) then
	        if(strucdim(str,5).gt.ny/2) then
	          j2=ny-j+1; k2=-k+centres(str,xcen(str),3)-ks
	        end if
	        i2=rmod(i2+nx,nx); k2=rmod(k2+nz,nz)
	        coco2=coco2+l2_array(i2,j2,k2)**2
	        coco=coco+struc_neg(i1,j1,k1)*l2_array(i2,j2,k2)
	      end do; end do; end do
	      coco2=coco2/real((coragx+2*iext)*(coragy+2*iext)*(coragz+2*iext))
	      coco=coco/real((coragx+2*iext)*(coragy+2*iext)*(coragz+2*iext))
	      coco=coco/sqrt(coco1_neg)/sqrt(coco2)
	      if(coco.ge.cocomax(str-struc_start+1)) then
	        cocomax(str-struc_start+1)=coco
	        strucshift(str,1)=is
	        strucshift(str,2)=ks
	      end if
	    end do; end do
!	    Write(*,*) 'cross-correlation_max=',cocomax(str-struc_start+1),'ishift=',strucshift(str,1),'kshift=',strucshift(str,2)
	    Write(*,*) 'cross-correlation_max=',cocomax(str-struc_start+1),strucshift(str,1),strucshift(str,2)
	  end if
	end do

	END SUBROUTINE struct_corr
!*******************************************************************************
!*******************************************************************************
	SUBROUTINE struct_average_lam(l2_array,wx,struc_start,struc_end)
!	This subroutine do the averaged lambda2 and wx' field
	Real(mytype), Dimension(nxt,nyt,nzt), Intent(in) :: l2_array
	Real(mytype), Dimension(nx,ny,nz), Intent(in) :: wx
	Integer, Intent(in) ::struc_start,struc_end

        Integer :: i,j,k,m,jin,kin,i1,j1,k1,i2,j2,k2,neg, &
                   str,iext,is,ks,ishift,kshift
        Integer, dimension(3) :: maxall
        Real(mytype) :: dysum,dzsum,xs,zs
	Character*10 :: chord,chdir

        maxall(1)=xw 
        maxall(2)=yw 
        maxall(3)=zw 
	iext=0
	do str=struc_start,struc_end
	    if(vort_dir(str).eq.1) cnt_pos=cnt_pos+1
	    if(vort_dir(str).eq.-1) cnt_neg=cnt_neg+1
	    do k=-iext-maxall(3)/2,iext+maxall(3)-maxall(3)/2-1; do j=1,maxall(2)+2*iext; do i=-iext-maxall(1)/2,iext+maxall(1)-maxall(1)/2-1
	      i1=i+iext+maxall(1)/2+1; j1=j; k1=k+iext+maxall(3)/2+1
	      i2=i+(centres(str,xcen(str),1)+strucshift(str,1)); j2=j; k2=k+(centres(str,xcen(str),3)+strucshift(str,2))
!	      if(itop.eq.1.and.strucdim(str,5).gt.ny/2) then
	      if(strucdim(str,5).gt.ny/2) then
	        j2=ny-j+1; k2=-k+(centres(str,xcen(str),3)-strucshift(str,2))
	      end if
	      i2=rmod(i2+nx,nx); k2=rmod(k2+nz,nz)
	      if(vort_dir(str).eq.1) then
	        struc_avg_pos(i1,j1,k1,3)=struc_avg_pos(i1,j1,k1,3)+l2_array(i2,j2,k2)  !lambda2
	        struc_avg_pos(i1,j1,k1,11)=struc_avg_pos(i1,j1,k1,11)+wx(i2,j2,k2) !wx
	      end if
	      if(vort_dir(str).eq.-1) then
	        struc_avg_neg(i1,j1,k1,3)=struc_avg_neg(i1,j1,k1,3)+l2_array(i2,j2,k2)  !lambda2
	        struc_avg_neg(i1,j1,k1,11)=struc_avg_neg(i1,j1,k1,11)+wx(i2,j2,k2) !wx
	      end if
	    end do; end do; end do
	end do
	Write(*,*) 'cnt_pos=',cnt_pos,'cnt_neg=',cnt_neg
	END SUBROUTINE struct_average_lam
!*******************************************************************************
	SUBROUTINE struct_average_vel(dummy_3d,struc_start,struc_end)
!	This subroutine do the averaged velocity field
	Real(mytype), Dimension(-1:nx+2,-1:ny+2,-1:nz+2,4), Intent(in) :: dummy_3d
!	Real(mytype), Dimension(nxt,nyt,nzt,3), Intent(in) :: vx
	Integer, Intent(in) ::struc_start,struc_end

        Integer :: i,j,k,m,jin,kin,i1,j1,k1,i2,j2,k2,neg, &
                   str,iext,is,ks,ishift,kshift
        Integer, dimension(3) :: maxall
        Real(mytype) :: dysum,dzsum,xs,zs
	Character*10 :: chord,chdir

        maxall(1)=xw 
        maxall(2)=yw 
        maxall(3)=zw 
	iext=0
	do str=struc_start,struc_end
	    if(vort_dir(str).eq.1) cnt_pos=cnt_pos+1
	    if(vort_dir(str).eq.-1) cnt_neg=cnt_neg+1
	    do j=1,struc_len(str)
	      i2=centres(str,j,1); j2=centres(str,j,2); k2=centres(str,j,3)
!        write(*,*) i2,k2
              if(strucdim(str,1).gt.strucdim(str,3)) then
	        if(centres(str,xcen(str),1).gt.nx/2) then
                  if(i2.lt.nx/2) i2=i2+nx
	        else
	          if(i2.gt.nx/2) i2=i2-nx
	        end if
              end if
              if(strucdim(str,7).gt.strucdim(str,9)) then
	        if(centres(str,xcen(str),3).gt.nz/2) then
                  if(k2.lt.nz/2) k2=k2+nz
	        else
	          if(k2.gt.nz/2) k2=k2-nz
	        end if
              end if
	      i2=i2-(centres(str,xcen(str),1)+strucshift(str,1)); k2=k2-(centres(str,xcen(str),3)+strucshift(str,2))
!	      if(itop.eq.1.and.strucdim(str,5).gt.ny/2) then
	      if(strucdim(str,5).gt.ny/2) then
	        j2=ny-j2+1; k2=-k2
	      end if
	      i2=i2+iext+maxall(1)/2+1; k2=k2+iext+maxall(3)/2+1
!        write(*,*) i2,k2
	      if(i2.ge.1.and.i2.le.maxall(1)+2*iext.and.k2.ge.1.and.k2.le.maxall(3)+2*iext) then
	        if(vort_dir(str).eq.1) then
	          struc_avg_pos(i2,j2,k2,2)=struc_avg_pos(i2,j2,k2,2)+1.0e0
	        else
	          struc_avg_neg(i2,j2,k2,2)=struc_avg_neg(i2,j2,k2,2)+1.0e0
	        end if
	      end if
	    end do
!
	    do k=-iext-maxall(3)/2,iext+maxall(3)-maxall(3)/2-1; do j=1,maxall(2)+2*iext; do i=-iext-maxall(1)/2,iext+maxall(1)-maxall(1)/2-1
	      i1=i+iext+maxall(1)/2+1; j1=j; k1=k+iext+maxall(3)/2+1
	      i2=i+(centres(str,xcen(str),1)+strucshift(str,1)); j2=j; k2=k+(centres(str,xcen(str),3)+strucshift(str,2))
!	      if(itop.eq.1.and.strucdim(str,5).gt.ny/2) then
	      if(strucdim(str,5).gt.ny/2) then
	        j2=ny-j+1; k2=-k+(centres(str,xcen(str),3)-strucshift(str,2))
	      end if
	      i2=rmod(i2+nx,nx); k2=rmod(k2+nz,nz)
	      if(vort_dir(str).eq.1) then
!	        struc_avg_pos(i1,j1,k1,3)=struc_avg_pos(i1,j1,k1,3)+l2_array(i2,j2,k2)  !lambda2
	        struc_avg_pos(i1,j1,k1,4)=struc_avg_pos(i1,j1,k1,4)+0.5*(dummy_3d(i2,j2,k2,1)+dummy_3d(i2+1,j2,k2,1)) !umean
	        If(iosci.Eq.1) Then
	          struc_avg_pos(i1,j1,k1,7)=struc_avg_pos(i1,j1,k1,7)+(0.5*(dummy_3d(i2,j2,k2,1)+dummy_3d(i2+1,j2,k2,1))-flucs1d_sa(j2,1,1,phase)) !u
	        Else
	          struc_avg_pos(i1,j1,k1,7)=struc_avg_pos(i1,j1,k1,7)+(0.5*(dummy_3d(i2,j2,k2,1)+dummy_3d(i2+1,j2,k2,1))-flucs(j2,1,1)) !u
	        End If
!	        if(itop.eq.1.and.strucdim(str,5).gt.ny/2) then     
	        if(strucdim(str,5).gt.ny/2) then
	          struc_avg_pos(i1,j1,k1,8)=struc_avg_pos(i1,j1,k1,8)-0.5*(dummy_3d(i2,j2,k2,2)+dummy_3d(i2,j2+1,k2,2)) !v
	          struc_avg_pos(i1,j1,k1,5)=struc_avg_pos(i1,j1,k1,5)-0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3)) !wmean
	          If(iosci.Eq.1) Then
	            struc_avg_pos(i1,j1,k1,9)=struc_avg_pos(i1,j1,k1,9)-(0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3))-flucs1d_sa(j2,3,1,phase)) !w
	          Else
	            struc_avg_pos(i1,j1,k1,9)=struc_avg_pos(i1,j1,k1,9)-(0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3))-flucs(j2,3,1)) !w
	          End If
	        else
	          struc_avg_pos(i1,j1,k1,8)=struc_avg_pos(i1,j1,k1,8)+0.5*(dummy_3d(i2,j2,k2,2)+dummy_3d(i2,j2+1,k2,2)) !v
	          struc_avg_pos(i1,j1,k1,5)=struc_avg_pos(i1,j1,k1,5)+0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3)) !wmean         
	          If(iosci.Eq.1) Then
	            struc_avg_pos(i1,j1,k1,9)=struc_avg_pos(i1,j1,k1,9)+(0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3))-flucs1d_sa(j2,3,1,phase)) !w
	          Else
	            struc_avg_pos(i1,j1,k1,9)=struc_avg_pos(i1,j1,k1,9)+(0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3))-flucs(j2,3,1)) !w
	          End If
	        end if
	        struc_avg_pos(i1,j1,k1,6)=struc_avg_pos(i1,j1,k1,6)+dummy_3d(i2,j2,k2,4) !pmean
	        struc_avg_pos(i1,j1,k1,10)=struc_avg_pos(i1,j1,k1,10)+(dummy_3d(i2,j2,k2,4)-flucs(j2,4,1)) !p
!	        struc_avg_pos(i1,j1,k1,11)=struc_avg_pos(i1,j1,k1,11)+wx(i2,j2,k2) !wx
	      end if
	      if(vort_dir(str).eq.-1) then
!	        struc_avg_neg(i1,j1,k1,3)=struc_avg_neg(i1,j1,k1,3)+l2_array(i2,j2,k2)  !lambda2
	        struc_avg_neg(i1,j1,k1,4)=struc_avg_neg(i1,j1,k1,4)+0.5*(dummy_3d(i2,j2,k2,1)+dummy_3d(i2+1,j2,k2,1)) !umean
	        If(iosci.Eq.1) Then
	          struc_avg_neg(i1,j1,k1,7)=struc_avg_neg(i1,j1,k1,7)+(0.5*(dummy_3d(i2,j2,k2,1)+dummy_3d(i2+1,j2,k2,1))-flucs1d_sa(j2,1,1,phase)) !u
	        Else
	          struc_avg_neg(i1,j1,k1,7)=struc_avg_neg(i1,j1,k1,7)+(0.5*(dummy_3d(i2,j2,k2,1)+dummy_3d(i2+1,j2,k2,1))-flucs(j2,1,1)) !u
	        End If
!	        if(itop.eq.1.and.strucdim(str,5).gt.ny/2) then     
	        if(strucdim(str,5).gt.ny/2) then
	          struc_avg_neg(i1,j1,k1,8)=struc_avg_neg(i1,j1,k1,8)-0.5*(dummy_3d(i2,j2,k2,2)+dummy_3d(i2,j2+1,k2,2)) !v
	          struc_avg_neg(i1,j1,k1,5)=struc_avg_neg(i1,j1,k1,5)-0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3)) !wmean
	          If(iosci.Eq.1) Then
	            struc_avg_neg(i1,j1,k1,9)=struc_avg_neg(i1,j1,k1,9)-(0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3))-flucs1d_sa(j2,3,1,phase)) !w
	          Else
	            struc_avg_neg(i1,j1,k1,9)=struc_avg_neg(i1,j1,k1,9)-(0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3))-flucs(j2,3,1)) !w
	          End If
	        else
	          struc_avg_neg(i1,j1,k1,8)=struc_avg_neg(i1,j1,k1,8)+0.5*(dummy_3d(i2,j2,k2,2)+dummy_3d(i2,j2+1,k2,2)) !v
	          struc_avg_neg(i1,j1,k1,5)=struc_avg_neg(i1,j1,k1,5)+0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3)) !wmean	          
	          If(iosci.Eq.1) Then
	            struc_avg_neg(i1,j1,k1,9)=struc_avg_neg(i1,j1,k1,9)+(0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3))-flucs1d_sa(j2,3,1,phase)) !w
	          Else
	            struc_avg_neg(i1,j1,k1,9)=struc_avg_neg(i1,j1,k1,9)+(0.5*(dummy_3d(i2,j2,k2,3)+dummy_3d(i2,j2,k2+1,3))-flucs(j2,3,1)) !w
	          End If
	        end if
	        struc_avg_neg(i1,j1,k1,6)=struc_avg_neg(i1,j1,k1,6)+dummy_3d(i2,j2,k2,4) !pmean
	        struc_avg_neg(i1,j1,k1,10)=struc_avg_neg(i1,j1,k1,10)+(dummy_3d(i2,j2,k2,4)-flucs(j2,4,1)) !p
!	        struc_avg_neg(i1,j1,k1,11)=struc_avg_neg(i1,j1,k1,11)+wx(i2,j2,k2) !wx
	      end if
	    end do; end do; end do
	end do
	Write(*,*) 'cnt_pos=',cnt_pos,'cnt_neg=',cnt_neg
	END SUBROUTINE struct_average_vel
!*******************************************************************************
	SUBROUTINE centre_calc(l2_array,pmin,vx,wx,i2,j2,k2,fb,struc_len1)
	Integer, Dimension(nxt,nyt,nzt), Intent(inout) :: pmin
	Real(mytype), Dimension(nxt,nyt,nzt), Intent(in) :: l2_array
	Real(mytype), Dimension(nxt,nyt,nzt,3), Intent(inout) :: vx
        Integer, intent(inout) :: i2,j2,k2,struc_len1
        Integer, intent(in) :: fb !fb: 1=forwards -1=backwards
	Real(mytype), Dimension(nx,ny,nz), Intent(in) :: wx

        Integer :: jin,kin,k1,j1,i1,struc_end1,neg, &
                   kmin,kmax,jmin,jmax,wx_neg,itmp,jtmp,ktmp
        real(mytype) :: dzsum, grd_l2, grd
!if(wx(i2,j2,k2).le.0.0e0) goto 1122
        struc_end1=0; wx_neg=1
        if (wx(i2,j2,k2).lt.0.0e0) wx_neg=-1

	  If(struc_len1.Lt.2) Then
            vx(i2,j2,k2,1)=real(fb)*dxs(i2)
            vx(i2,j2,k2,2)=0.0e0
            vx(i2,j2,k2,3)=0.0e0
	  Else
	    vx(i2,j2,k2,1)=-vx(centres(struc_num,2,1),centres(struc_num,2,2),centres(struc_num,2,3),1)
	    vx(i2,j2,k2,2)=-vx(centres(struc_num,2,1),centres(struc_num,2,2),centres(struc_num,2,3),2)
	    vx(i2,j2,k2,3)=-vx(centres(struc_num,2,1),centres(struc_num,2,2),centres(struc_num,2,3),3)
	  End If
            if (struc_len1 == 0) then
              struc_len1=struc_len1+1
              centres(struc_num,struc_len1,1) = i2
              centres(struc_num,struc_len1,2) = j2
              centres(struc_num,struc_len1,3) = k2
            end if    
            do while(struc_end1==0.And.struc_len1.Lt.50)   
              !pmin(i2,j2,k2)=0
              !vx(i2,j2,k2,:)=fb*vx(i2,j2,k2,:)
              !Write(26,'(6E15.7)')xc(i2),yc(j2),zc(k2),vx(i2,j2,k2,1),vx(i2,j2,k2,2),vx(i2,j2,k2,3)
              call dist_calc(vx(i2,j2,k2,:),struc_len1,i2,j2,k2,jmin,jmax,kmin,kmax)
              struc_end1=1
	      grd_l2=1.0e8; itmp=0; jtmp=0; ktmp=0
              do k1=kmin,kmax; do j1=jmin,jmax
		i1=fb
                if ((j2+j1).ge.1 .and. (j2+j1).le.ny) then
                  !write(*,*) 'checking: ',rmod(i2+i1,nx),j2+j1,rmod(k2+k1,nz)
                  !if (j2+j1.le.jvis ) Write(27,'(6E15.7)') xc(rmod(i2+i1,nx)),yc(j2+j1),zc(rmod(k2+k1,nz))
	          if (pmin(rmod(i2+i1,nx),j2+j1,rmod(k2+k1,nz)) == 1 .and. &
                      real(wx_neg)*wx(rmod(i2+i1,nx),j2+j1,rmod(k2+k1,nz)) .ge. 0.0e0) then
	            grd=Abs(l2_array(i2,j2,k2)-l2_array(rmod(i2+i1,nx),j2+j1,rmod(k2+k1,nz)))/ &
			Sqrt((xc(i2)-xc(rmod(i2+i1,nx)))**2+(yc(j2)-yc(j2+j1))**2+(zc(k2)-zc(rmod(k2+k1,nz)))**2)
	            If(grd.Lt.grd_l2) Then
	              grd_l2=grd
	              itmp=i1; jtmp=j1; ktmp=k1
	            End If
                    struc_end1=0
                  end if
                else
               ! write(*,*) 'did not check: ',rmod(i2+i1,nx),j2+j1,rmod(k2+k1,nz)
                end if
              end do; end do
	      If(struc_end1.Eq.0) Then
!               set new vector from current vortex direction
	        i1=itmp; j1=jtmp; k1=ktmp
                vx(rmod(i2+i1,nx),j2+j1,rmod(k2+k1,nz),1)=xc(rmod(i2+i1,nx))-xc(i2)
                vx(rmod(i2+i1,nx),j2+j1,rmod(k2+k1,nz),2)=yc(j2+j1)-yc(j2)
                vx(rmod(i2+i1,nx),j2+j1,rmod(k2+k1,nz),3)=zc(rmod(k2+k1,nz))-zc(k2)
                i2=rmod(i2+i1,nx);j2=j2+j1;k2=rmod(k2+k1,nz)
!               catch overlap
                if (vx(i2,j2,k2,1).gt.length/2.0e0) vx(i2,j2,k2,1)=length-vx(i2,j2,k2,1)
                if (vx(i2,j2,k2,3).gt.width/2.0e0) vx(i2,j2,k2,3)=width-vx(i2,j2,k2,3)
                struc_len1=struc_len1+1
                centres(struc_num,struc_len1,1) = i2
                centres(struc_num,struc_len1,2) = j2
                centres(struc_num,struc_len1,3) = k2
	      End If
         16 end do  
1122 continue
	END SUBROUTINE centre_calc
!*******************************************************************************

!*******************************************************************************
	SUBROUTINE dist_calc(vec,struc_len1,i2,j2,k2,jmin,jmax,kmin,kmax)
	Real(mytype), dimension(3), intent(in) :: vec
        Integer, intent(in) :: i2,j2,k2,struc_len1
        Integer, intent(out) :: jmin,jmax,kmin,kmax

        integer :: neg
        real(mytype) :: dispy,dispz,dysum,dzsum,lim

!	If(struc_len1.Eq.1) Then
!          dispy=4.0e0*dzs(1)  ! search angle
!          dispz=2.0e0*dzs(1)
!	Else
          dispy=2.0e0*dzs(1)  ! search angle
          dispz=2.0e0*dzs(1)
!	End If

!       y direction
        jmin=0; neg=1; dysum=0.0e0
        if (vec(2).lt.0.0e0) neg=-1
        lim=abs(vec(2))-real(neg)*dispy
        if (lim.lt.0.0e0) neg=-1*neg
        do while(dysum.lt.abs(lim)) 
          jmin=jmin+neg
          if ((j2+jmin).ge.1 .and. (j2+jmin).le.ny) then
            dysum=dysum+dys(j2+jmin)
          else
            dysum=abs(lim)
          end if
        end do  

        jmax=0; neg=1; dysum=0
        if (vec(2).lt.0.0e0) neg=-1
        lim=abs(vec(2))+real(neg)*dispy
        if (lim.lt.0.0e0) neg=-1*neg
        do while(dysum.lt.abs(lim)) 
          jmax=jmax+neg
          if ((j2+jmax).ge.1 .and. (j2+jmax).le.ny) then
            dysum=dysum+dys(j2+jmax)
          else
            dysum=abs(lim)
          end if
        end do    
!        if (vec(2).lt.0.0e0) then; jmax=jmax+1
!        else; jmin=jmin-1; end if
        
!       k direction
        kmin=0; neg=1; dzsum=0.0e0
        if (vec(3).lt.0.0e0) neg=-1
        lim=abs(vec(3))-real(neg)*dispz
        if (lim.lt.0.0e0) neg=-1*neg
        do while(dzsum.lt.abs(lim)) 
          kmin=kmin+neg
          dzsum=dzsum+dzs(rmod(k2+kmin,nz))
        end do  

        kmax=0; neg=1; dzsum=0
        if (vec(3).lt.0.0e0) neg=-1
        lim=abs(vec(3))+real(neg)*dispz
        if (lim.lt.0.0e0) neg=-1*neg
        do while(dzsum.lt.abs(lim)) 
          kmax=kmax+neg
          dzsum=dzsum+dzs(rmod(k2+kmax,nz))
        end do  
!        if (vec(3).lt.0.0e0) then; kmax=kmax+1
!        else; kmin=kmin-1; end if
!Write(*,*) 'i2,j2,k2,xc,yc,zc',i2,j2,k2,xc(i2),yc(j2),zc(k2)
!Write(*,*) 'dispy,dispz',dispy,dispz
!Write(*,*) 'vec',vec(1),vec(2),vec(3)
!Write(*,*) 'kmin,kmax,zc(k2+kmin),zc(k2+kmax)',kmin,kmax,zc(k2+kmin),zc(k2+kmax)


       END SUBROUTINE dist_calc
!*******************************************************************************
	SUBROUTINE middle_calc(centres,lenl,struc_len1,maxall,minx,strucdim1)
        Integer :: struc_len1,minx
        Integer, dimension(struc_len1,4), Intent(in) :: centres
        Integer, dimension(3), Intent(inout) :: maxall
	Integer, dimension(1:9), Intent(out) :: strucdim1
        Real(mytype), dimension(3), Intent(out) :: lenl

        Integer, dimension(3) :: maxn,maxl,minl
        Integer :: k

        maxn(1)=nx; maxn(2)=ny; maxn(3)=nz

!         calulates upper and lower extremeties of structure
          do k=1,3
            if (minval(centres(1:struc_len1,k))==1) then
              minl(k)=minval(centres(1:struc_len1,k),MASK=centres(1:struc_len1,k).gt.maxn(k)/2)
              if (minl(k).gt.maxn(k)) minl(k)=1
            else 
              minl(k)=minval(centres(1:struc_len1,k))
            end if      
	    strucdim1(3*k-2)=minl(k)

            if (maxval(centres(1:struc_len1,k))==maxn(k)) then
              maxl(k)=maxval(centres(1:struc_len1,k),MASK=centres(1:struc_len1,k).lt.maxn(k)/2)
              if (maxl(k).lt.1) maxl(k)=maxn(k)
            else 
              maxl(k)=maxval(centres(1:struc_len1,k))
            end if    
	    strucdim1(3*k)=maxl(k)
	    if(minl(k).gt.maxl(k)) then
	      strucdim1(3*k-1)=rmod((minl(k)+maxl(k)+maxn(k))/2,maxn(k))          
	    else
	      strucdim1(3*k-1)=rmod((minl(k)+maxl(k))/2,maxn(k))          
	    end if
	  end do

!          write(*,*) minl(1),maxl(1),minl(2),maxl(2),minl(3),maxl(3)
!	  write(*,*) strucdim1(1:9)
!         calculates length,height and width of structure
          if (minl(1).lt.maxl(1)) then
            lenl(1)=xc(maxl(1))-xc(minl(1))
            !if ((maxl(1)-minl(1)+1).gt.maxall(1)) maxall(1)=maxl(1)-minl(1)+1
            maxall(1)=maxall(1)+maxl(1)-minl(1)+1
          else
            lenl(1)=xc(maxl(1))-xc(1)+xc(nx+1)-xc(minl(1))
            !if ((maxl(1)+nx-minl(1)+1).gt.maxall(1)) maxall(1)=maxl(1)+nx-minl(1)+1
            maxall(1)=maxall(1)+maxl(1)+nx-minl(1)+1
          end if

          lenl(2)=yc(maxl(2))-yc(minl(2))
          !if ( maxl(2).gt.maxall(2)) maxall(2)=maxl(2)
          maxall(2)=maxall(2)+maxl(2)
 
          if (minl(3).lt.maxl(3)) then
            lenl(3)=zc(maxl(3))-zc(minl(3))
            !if ((maxl(3)-minl(3)+1).gt.maxall(3)) maxall(3)=maxl(3)-minl(3)+1
            maxall(3)=maxall(3)+maxl(3)-minl(3)+1
          else
            lenl(3)=zc(maxl(3))-zc(3)+zc(nz+1)-zc(minl(3))
            !if ((maxl(3)+nz-minl(3)+1).gt.maxall(3)) maxall(3)=maxl(3)+nz-minl(3)+1
            maxall(3)=maxall(3)+maxl(3)+nz-minl(3)+1
          end if

        minx=minl(1)

        return
	END SUBROUTINE middle_calc

	END MODULE calc_struct

















