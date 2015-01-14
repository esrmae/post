
	MODULE write_post
!	USE shared_data
	USE shared_post
	IMPLICIT NONE
!*******************************************************************************
! 	This code writes all the information to the correct files.
!*******************************************************************************
	CONTAINS
!*******************************************************************************
	SUBROUTINE write_vals
!
        If(tstart0 == 1) time_sa=time_sa/time_sa(tpts)
	CALL y_print_calc
	CALL open_files_0d
	CALL write_0d_data
	If (isteady == 1) then
	    CALL open_files_steady1d
	    If (lpstat_1d == 1) CALL write_1d_avg
	Else If (isteady == 0) then
	    CALL open_files_unst
	    CALL open_1dpts_unst
	    CALL write_1d_unst
	    If (uns1d.ne.0) CALL write_1dpts_unst
	    CALL write_rates_unst
	Else If (isteady == 2) then
	    CALL open_files_steady2d
	    If (lpstat_2d == 1) CALL write_2d_avg
	End If
	CALL close_files
!
	END SUBROUTINE write_vals
!*******************************************************************************

	SUBROUTINE open_file1d_point(file_num,filename,char_arr,num_vals,ivals)
	INTEGER :: file_num,num_vals,i,ivals
	CHARACTER :: filename*40
	CHARACTER (LEN=100), DIMENSION(29) :: char_arr
!
	If (tcplt == 1) then
		filename=TRIM(filename)//'.plt'
	Else
		filename=TRIM(filename)//'.dat'	
	End If
	Open(Unit=file_num, File=TRIM(filename), Status='UNKNOWN')
    If (tcplt==1) then
	Write(file_num,*) 'title="',TRIM(filename),'"'
	Write(file_num,'(A)',advance='no') ' variables='
	DO i=1,num_vals-1
	Write(file_num,'(3A)',advance='no') '"',TRIM(char_arr(i)),'",'
	END DO 
	Write(file_num,'(3A)') '"',TRIM(char_arr(num_vals)),'"'
	Write(file_num,'(A,I6)') ' zone f=point, i =', ivals
    End If
!
	END SUBROUTINE open_file1d_point
!*******************************************************************************

	SUBROUTINE open_file1d(file_num,filename,char_arr,num_vals,ivals)
	INTEGER :: file_num,num_vals,i,ivals
	CHARACTER :: filename*40
	CHARACTER (LEN=100), DIMENSION(29) :: char_arr
!
	If (tcplt == 1) then
		filename=TRIM(filename)//'.plt'
	Else
		filename=TRIM(filename)//'.dat'	
	End If
	Open(Unit=file_num, File=TRIM(filename), Status='UNKNOWN')
    If (tcplt==1) then
	Write(file_num,*) 'title="',TRIM(filename),'"'
	Write(file_num,'(A)',advance='no') ' variables='
	DO i=1,num_vals-1
	Write(file_num,'(3A)',advance='no') '"',TRIM(char_arr(i)),'",'
	END DO 
	Write(file_num,'(3A)') '"',TRIM(char_arr(num_vals)),'"'
	Write(file_num,'(A,I6)') ' zone i =', ivals, ',f=point'
    End If
!
	END SUBROUTINE open_file1d
!*******************************************************************************

	SUBROUTINE open_file2d(file_num,filename,char_arr,num_vals,ivals,jvals)
	INTEGER :: file_num,num_vals,i,ivals,jvals
	CHARACTER :: filename*40
	CHARACTER (LEN=100), DIMENSION(29) :: char_arr
!
	filename=TRIM(filename)//'.plt'
	Open(Unit=file_num, File=TRIM(filename), Status='UNKNOWN')
    !If (tcplt==1) then
	Write(file_num,*) 'title="',TRIM(filename),'"'
	Write(file_num,'(A)',advance='no') ' variables='
	DO i=1,num_vals-1
	Write(file_num,'(3A)',advance='no') '"',TRIM(char_arr(i)),'",'
	END DO 
	Write(file_num,'(3A)') '"',TRIM(char_arr(num_vals)),'"'
	Write(file_num,*) 'zone t="z1",i=',ivals,',j=',jvals,',f=point'
    !End If
!
	END SUBROUTINE open_file2d
!*******************************************************************************

	SUBROUTINE y_print_calc
!
	   jmin = 0
	   jmax = nyt
	   jmaxt = nyt
	If (ibm == 1) then
	   jmin = 1+ny_solid
	   jmax = ny-ny_solid
	   jmaxt = nyt-ny_solid
	End If
!
	If (y_half == 1) then
	   jmax = nyt/2+1
	   jmaxt = nyt/2+1
	End If
	   jtot = jmax-jmin+1
	   jtott = jmaxt-jmin+1
!
	END SUBROUTINE y_print_calc
!*******************************************************************************
	SUBROUTINE open_files_0d
	CHARACTER (LEN=100), DIMENSION(29) :: dum
	CHARACTER (LEN=1), DIMENSION(4) :: uvwp
	CHARACTER :: filedum*40
	INTEGER :: i
	uvwp=(/'u','v','w','p'/)
	dum=''
!       Open 0 DIMENSIONAL files
	dum(1)='time'
	dum(2)='u<sub><greek>t</greek></sub>'
	filedum='utau-xza'
	CALL open_file1d(35,filedum,dum,2,tpts)
!
	dum(1)='time'
	dum(2)='Re<sub><greek>t</greek></sub>'
	dum(3)='Re'
	filedum='retau-xza'
	CALL open_file1d(36,filedum,dum,3,tpts)
!
	dum(1)='time'
	dum(2)='dp/dx<sub>mean</sub>'
	filedum='dpdx_mean-xza'
	CALL open_file1d(37,filedum,dum,2,tpts)
!
!	dum(1)='y'
!	dum(2)='u'
!	dum(3)='v'
!	dum(4)='w'
!	dum(5)='<greek>w</greek><sub>x</sub>'
!	dum(6)='<greek>w</greek><sub>y</sub>'
!	dum(7)='<greek>w</greek><sub>z</sub>'
!	filedum='rmsofrms'
!	CALL open_file1d(38,filedum,dum,7,ny)

!	dum(1)='y'
!	dum(2)='u'
!	dum(3)='v'
!	dum(4)='w'
!	filedum='avgflucs'
!	CALL open_file1d(39,filedum,dum,4,ny)
!
!
	If (lp_snap_1d == 1) then
	    If (lp_snap_y == 1) then
		DO i=1,3
		    dum(1)='time'
		    dum(2)=uvwp(i)//'<sub>(ny/32)</sub>'
		    dum(3)=uvwp(i)//'<sub>(ny/16)</sub>'
		    dum(4)=uvwp(i)//'<sub>(ny/8)</sub>'
		    dum(5)=uvwp(i)//'<sub>(3ny/16)</sub>'
		    dum(6)=uvwp(i)//'<sub>(ny/4)</sub>'
		    dum(7)=uvwp(i)//'<sub>(3ny/8)</sub>'
		    dum(8)=uvwp(i)//'<sub>(ny/2)</sub>'
		    filedum='freqy'//uvwp(i)
		    CALL open_file1d(29+i,filedum,dum,8,tpts)
		END DO
!
		dum(1)='time'
		dum(2)='p<sub>(1)</sub>'
		dum(3)='p<sub>(ny/6)</sub>'
		dum(4)='p<sub>(ny/3)</sub>'
		dum(5)='p<sub>(ny/2)</sub>'
		dum(6)='p<sub>(2ny/3)</sub>'
		dum(7)='p<sub>(5ny/6)</sub>'
		dum(8)='p<sub>(ny)</sub>'
		filedum='freqyp'
		CALL open_file1d(33,filedum,dum,8,tpts)
	    End If
	End if
!
	END SUBROUTINE open_files_0d
!*******************************************************************************
	SUBROUTINE open_files_steady1d
!       1 DIMENSIONAL TIME AVERAGE
	CHARACTER (LEN=100), DIMENSION(29) :: dum
	CHARACTER (LEN=1), DIMENSION(4) :: uvwp
	CHARACTER :: filedum*40,one*1,three*3
	INTEGER :: i
	uvwp=(/'u','v','w','p'/)
!
!	Opens files	
!       1 DIMENSIONAL TIME AVERAGE
	If (iflucs == 1) then
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='u'
	    dum(3)='v'
	    dum(4)='w'
	    dum(5)='p'
	    dum(6)='yc'
	    filedum='mean-txza'
	    CALL open_file1d(41,filedum,dum,6,jtot)
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='u<sup>+</sup>'
	    dum(3)='v<sup>+</sup>'
	    dum(4)='w<sup>+</sup>'
	    dum(5)='p<sup>+</sup>'
	    dum(6)='yc'
	    filedum='meanp-txza'
	    CALL open_file1d(48,filedum,dum,6,jtot)
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='u<sub>rms</sub>'
	    dum(3)='v<sub>rms</sub>'
	    dum(4)='w<sub>rms</sub>'
	    dum(5)='p<sub>rms</sub>'
	    dum(6)='yc'
	    filedum='flucs-txza'
	    CALL open_file1d(42,filedum,dum,6,jtot)
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='u<sub>rms</sub><sup>+</sup>'
	    dum(3)='v<sub>rms</sub><sup>+</sup>'
	    dum(4)='w<sub>rms</sub><sup>+</sup>'
	    dum(5)='p<sub>rms</sub><sup>+</sup>'
	    dum(6)='yc'
	    filedum='flucsp-txza'
	    CALL open_file1d(53,filedum,dum,6,jtot)
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='u<sub>skew</sub>'
	    dum(3)='v<sub>skew</sub>'
	    dum(4)='w<sub>skew</sub>'
	    dum(5)='p<sub>skew</sub>'
	    dum(6)='yc'
	    filedum='skew-txza'
	    CALL open_file1d(45,filedum,dum,6,jtot)
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='u<sub>kurt</sub>'
	    dum(3)='v<sub>kurt</sub>'
	    dum(4)='w<sub>kurt</sub>'
	    dum(5)='p<sub>kurt</sub>'
	    dum(6)='yc'
	    filedum='kurt-txza'
	    CALL open_file1d(46,filedum,dum,6,jtot)
	End If
!
	If (ishear == 1) then
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='<greek>t</greek><sub>turbx</sub>'
	    dum(3)='<greek>t</greek><sub>lamx</sub>'
	    dum(4)='<greek>t</greek><sub>turbz</sub>'
	    dum(5)='<greek>t</greek><sub>lamz</sub>'
	    dum(6)='yc'
	    filedum='shear-txza'
	    CALL open_file1d(44,filedum,dum,6,jtot)
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='<greek>t</greek><sub>turbx</sub><sup>+</sup>'
	    dum(3)='<greek>t</greek><sub>lamx</sub><sup>+</sup>'
	    dum(4)='<greek>t</greek><sub>SGS</sub><sup>+</sup>'
	    dum(5)='<greek>t</greek><sub>total</sub><sup>+</sup>'
	    dum(6)='yc'
	    filedum='shearp-txza'
	    CALL open_file1d(54,filedum,dum,6,jtot)
	End If
!
	If (iangle == 1) then
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='<greek>b</greek>'
	    dum(3)='<greek>a</greek>'
	    dum(4)='<greek>e</greek>'
	    dum(5)='<greek>g</greek>'
	    dum(6)='a<sub>1</sub>'
	    dum(7)='yc'
	    filedum='angle-txza'
	    CALL open_file1d(401,filedum,dum,7,jtot)
	End If
!
	If (ivort == 1) then
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='<greek>w</greek><sub>x(rms)</sub>'
	    dum(3)='<greek>w</greek><sub>y(rms)</sub>'
	    dum(4)='<greek>w</greek><sub>z(rms)</sub>'
	    dum(5)='<greek>w</greek><sub>x(mean)</sub>'
	    dum(6)='<greek>w</greek><sub>y(mean)</sub>'
	    dum(7)='<greek>w</greek><sub>z(mean)</sub>'
	    dum(8)='yc'
	    filedum='rmsvort-txza'
	    CALL open_file1d(43,filedum,dum,8,jtot)
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='<greek>w</greek><sub>x(rms)</sub><sup>+</sup>'
	    dum(3)='<greek>w</greek><sub>y(rms)</sub><sup>+</sup>'
	    dum(4)='<greek>w</greek><sub>z(rms)</sub><sup>+</sup>'
	    dum(5)='<greek>w</greek><sub>x(mean)</sub><sup>+</sup>'
	    dum(6)='<greek>w</greek><sub>y(mean)</sub><sup>+</sup>'
	    dum(7)='<greek>w</greek><sub>z(mean)</sub><sup>+</sup>'
	    dum(8)='yc'
	    filedum='rmsvortp-txza'
	    CALL open_file1d(55,filedum,dum,8,jtot)
	End If
!
	If (iquad == 1) then
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='Q1'
	    dum(3)='Q2'
	    dum(4)='Q3'
	    dum(5)='Q4'
	    dum(6)='uv'
	    dum(7)='vw'
	    dum(8)='uw'	    
	    dum(9)='yc'
	    filedum='quad-txza'
	    CALL open_file1d(47,filedum,dum,9,jtot)
	End If
!
	If (itke ==1) then
	    DO i=1,3
	        Write(three,'(I3)') i
		one=TRIM(adjustl(three))
	        dum(1)='yc<sup>+</sup>'
		dum(2)='P<sub>'//one//one//'</sub>'
		dum(3)='T<sub>'//one//one//'</sub>'
		dum(4)='<greek>P</greek><sub>'//one//one//'</sub>'
		dum(5)='<greek>F</greek><sub>'//one//one//'</sub>'
		dum(6)='D<sub>'//one//one//'</sub>'
		dum(7)='<greek>e</greek><sub>'//one//one//'</sub>'
		dum(8)='yc'
		filedum='tke'//uvwp(i)//uvwp(i)//'-txza'
		CALL open_file1d(48+i,filedum,dum,8,jtot)
	    END DO
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='P<sub>12</sub>'
	    dum(3)='T<sub>12</sub>'
	    dum(4)='<greek>P</greek><sub>12</sub>'
	    dum(5)='<greek>F</greek><sub>12</sub>'
	    dum(6)='D<sub>12</sub>'
	    dum(7)='<greek>e</greek><sub>12</sub>'
	    dum(8)='yc'
	    filedum='tkeuv-txza'
	    CALL open_file1d(52,filedum,dum,8,jtot)
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='P<sub>k</sub>'
	    dum(3)='T<sub>k</sub>'
	    dum(4)='<greek>P</greek><sub>k</sub>'
	    dum(5)='<greek>F</greek><sub>k</sub>'
	    dum(6)='D<sub>k</sub>'
	    dum(7)='<greek>e</greek><sub>k</sub>'
	    dum(8)='yc'
	    filedum='tkek-txza'
	    CALL open_file1d(66,filedum,dum,8,jtot)
!
	    DO i=1,3
	        Write(three,'(I3)') i
		one=TRIM(adjustl(three))
	        dum(1)='yc<sup>+</sup>'
		dum(2)='P<sub>'//one//one//'</sub><sup>+</sup>'
		dum(3)='T<sub>'//one//one//'</sub><sup>+</sup>'
		dum(4)='<greek>P</greek><sub>'//one//one//'</sub><sup>+</sup>'
		dum(5)='<greek>F</greek><sub>'//one//one//'</sub><sup>+</sup>'
		dum(6)='D<sub>'//one//one//'</sub><sup>+</sup>'
		dum(7)='<greek>e</greek><sub>'//one//one//'</sub><sup>+</sup>'
		dum(8)='yc'
		filedum='tke'//uvwp(i)//uvwp(i)//'p-txza'
		CALL open_file1d(61+i,filedum,dum,8,jtot)
	    END DO
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='P<sub>12</sub><sup>+</sup>'
	    dum(3)='T<sub>12</sub><sup>+</sup>'
	    dum(4)='<greek>P</greek><sub>12</sub><sup>+</sup>'
	    dum(5)='<greek>F</greek><sub>12</sub><sup>+</sup>'
	    dum(6)='D<sub>12</sub><sup>+</sup>'
	    dum(7)='<greek>e</greek><sub>12</sub><sup>+</sup>'
	    dum(8)='yc'
	    filedum='tkeuvp-txza'
	    CALL open_file1d(65,filedum,dum,8,jtot)
!
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='P<sub>k</sub><sup>+</sup>'
	    dum(3)='T<sub>k</sub><sup>+</sup>'
	    dum(4)='<greek>P</greek><sub>k</sub><sup>+</sup>'
	    dum(5)='<greek>F</greek><sub>k</sub><sup>+</sup>'
	    dum(6)='D<sub>k</sub><sup>+</sup>'
	    dum(7)='<greek>e</greek><sub>k</sub><sup>+</sup>'
	    dum(8)='yc'
	    filedum='tkekp-txza'
	    CALL open_file1d(67,filedum,dum,8,jtot)
	End If
!
	If (ianis==1) then
	    dum(1)='yc<sup>+</sup>'
	    dum(2)='b<sub>11</sub>'
	    dum(3)='b<sub>22</sub>'
	    dum(4)='b<sub>33</sub>'
	    dum(5)='b<sub>12</sub>'
	    dum(6)='II'
	    dum(7)='III'
	    dum(8)='F'
	    dum(9)='G'
	    dum(10)='yc'
	    filedum='anis-txza'
	    CALL open_file1d(61,filedum,dum,10,jtot)
	End If
!
	END SUBROUTINE open_files_steady1d
!*******************************************************************************
	SUBROUTINE open_files_steady2d
!       1 DIMENSIONAL TIME AVERAGE
	CHARACTER (LEN=100), DIMENSION(29) :: dum
	CHARACTER (LEN=1), DIMENSION(4) :: uvwp
	CHARACTER :: filedum*40,one*1,three*3
	INTEGER :: i
	uvwp=(/'u','v','w','p'/)
!
!	Opens files	
!       1 DIMENSIONAL TIME AVERAGE
	dum(1)='xc<sup>+</sup>'
	dum(2)='yc<sup>+</sup>'
	dum(3)='u'
	dum(4)='v'
	dum(5)='w'
	filedum='umean2d-txza'
	CALL open_file2d(110,filedum,dum,5,nx,ny)
!
	dum(1)='xc<sup>+</sup>'
	dum(2)='yc<sup>+</sup>'
	dum(3)='u<sup>+</sup>'
	dum(4)='v<sup>+</sup>'
	dum(5)='w<sup>+</sup>'
	filedum='umeanp2d-txza'
	CALL open_file2d(111,filedum,dum,5,nx,ny)
!
	If (iflucs == 1) then
	    dum(1)='xc<sup>+</sup>'
	    dum(2)='yc<sup>+</sup>'
	    dum(3)='u<sub>rms</sub>'
	    dum(4)='v<sub>rms</sub>'
	    dum(5)='w<sub>rms</sub>'
	    dum(6)='p<sub>rms</sub>'
	    filedum='flucs2d-txza'
	    CALL open_file2d(112,filedum,dum,6,nx,ny)
!
	    dum(1)='xc<sup>+</sup>'
	    dum(2)='yc<sup>+</sup>'
	    dum(3)='u<sub>rms</sub><sup>+</sup>'
	    dum(4)='v<sub>rms</sub><sup>+</sup>'
	    dum(5)='w<sub>rms</sub><sup>+</sup>'
	    dum(6)='p<sub>rms</sub><sup>+</sup>'
	    filedum='flucsp2d-txza'
	    CALL open_file2d(113,filedum,dum,6,nx,ny)
!
	    dum(1)='xc<sup>+</sup>'
	    dum(2)='yc<sup>+</sup>'
	    dum(3)='u<sub>skew</sub>'
	    dum(4)='v<sub>skew</sub>'
	    dum(5)='w<sub>skew</sub>'
	    dum(6)='p<sub>skew</sub>'
	    filedum='skew2d-txza'
	    CALL open_file2d(114,filedum,dum,6,nx,ny)
!
	    dum(1)='xc<sup>+</sup>'
	    dum(2)='yc<sup>+</sup>'
	    dum(3)='u<sub>kurt</sub>'
	    dum(4)='v<sub>kurt</sub>'
	    dum(5)='w<sub>kurt</sub>'
	    dum(6)='p<sub>kurt</sub>'
	    filedum='kurt2d-txza'
	    CALL open_file2d(115,filedum,dum,6,nx,ny)
	End If
!
	If (ishear == 1) then
	    dum(1)='xc<sup>+</sup>'
	    dum(2)='yc<sup>+</sup>'
	    dum(3)='<greek>t</greek><sub>turbx</sub>'
	    dum(4)='<greek>t</greek><sub>lamx</sub>'
	    dum(5)='<greek>t</greek><sub>turbz</sub>'
	    dum(6)='<greek>t</greek><sub>lamz</sub>'
	    filedum='shear2d-txza'
	    CALL open_file2d(116,filedum,dum,6,nx,ny)
	End If
	If (iangle == 1) then
	    dum(1)='xc<sup>+</sup>'
	    dum(2)='yc<sup>+</sup>'
	    dum(3)='<greek>b</greek>'
	    dum(4)='<greek>a</greek>'
	    dum(5)='<greek>e</greek>'
	    dum(6)='<greek>g</greek>'
	    dum(7)='a<sub>1</sub>'
	    filedum='angle2d-txza'
	    CALL open_file2d(402,filedum,dum,7,nx,ny)
	End If
!
	END SUBROUTINE open_files_steady2d
!*******************************************************************************
	SUBROUTINE open_files_unst
	INTEGER :: IC,ifin
	CHARACTER (LEN=100), DIMENSION(29) :: dum
	CHARACTER (LEN=1), DIMENSION(4) :: uvwp
	CHARACTER :: filedum*40,one*1,three*3,re_val*30,y_val*15
	INTEGER :: i,j,p
	INTEGER, DIMENSION(:), ALLOCATABLE :: twrite
	uvwp=(/'u','v','w','p'/)
	If (re_print==0) then; re_val='Re'
	Else If (re_print==1) then; re_val='Re<sub><greek>t</greek></sub>'
	Else; re_val='time'
	End If 
	If (y_print==0) then; y_val='yc'
	Else; y_val='yc<sup>+</sup>'
	End If
!
	dum(1)='time'
	dum(2)=y_val
	dum(3)='u'
	dum(4)='v'
	dum(5)='w'
	dum(6)='p'
	filedum='mean-unst-ty'
	CALL open_file2d(70,filedum,dum,6,tpts,jtot)
!
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='u'
	dum(4)='v'
	dum(5)='w'
	dum(6)='p'
	filedum='mean-unst-rey'
	CALL open_file2d(71,filedum,dum,6,tpts,jtot)
!
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='u<sub>mean</sub><sup>+</sup>'
	dum(4)='v<sub>mean</sub><sup>+</sup>'
	dum(5)='w<sub>mean</sub><sup>+</sup>'
	dum(6)='p<sub>mean</sub><sup>+</sup>'
	filedum='meanp-unst-rey'
	CALL open_file2d(90,filedum,dum,6,tpts,jtot)
!
	!dum(1)='time'
	!dum(2)=re_val
	!CALL open_file1d(76,filedum,dum,2,tpts)
!
!
	If (iflucs == 1) then
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='u<sub>rms</sub>'
	dum(4)='v<sub>rms</sub>'
	dum(5)='w<sub>rms</sub>'
	dum(6)='p<sub>rms</sub>'
	filedum='rms-unst-rety'
	CALL open_file2d(78,filedum,dum,6,tpts,jtot)
!
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='u<sup>+</sup><sub>rms</sub>'
	dum(4)='v<sup>+</sup><sub>rms</sub>'
	dum(5)='w<sup>+</sup><sub>rms</sub>'
	dum(6)='p<sup>+</sup><sub>rms</sub>'
	filedum='rmsp-unst-rey'
	CALL open_file2d(82,filedum,dum,6,tpts,jtot)
!
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='u<sup>+</sup><sub>skew</sub>'
	dum(4)='v<sup>+</sup><sub>skew</sub>'
	dum(5)='w<sup>+</sup><sub>skew</sub>'
	dum(6)='p<sup>+</sup><sub>skew</sub>'
	filedum='skew-unst-rey'
	CALL open_file2d(72,filedum,dum,6,tpts,jtot)
!
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='u<sup>+</sup><sub>kurt</sub>'
	dum(4)='v<sup>+</sup><sub>kurt</sub>'
	dum(5)='w<sup>+</sup><sub>kurt</sub>'
	dum(6)='p<sup>+</sup><sub>kurt</sub>'
	filedum='kurt-unst-rey'
	CALL open_file2d(73,filedum,dum,6,tpts,jtot)
	End If
!
	If (ishear == 1) then
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='<greek>t</greek><sub>turbx</sub>'
	dum(4)='<greek>t</greek><sub>lamx</sub>'
	dum(5)='<greek>t</greek><sub>turby</sub>'
	dum(6)='<greek>t</greek><sub>turbz</sub>'
	filedum='shear-unst-rey'
	CALL open_file2d(83,filedum,dum,6,tpts,jtot)
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='<greek>t</greek><sub>turbx</sub>'
	dum(4)='<greek>t</greek><sub>lamx</sub>'
	dum(5)='-vw'
	dum(6)='dWdy'
	filedum='shearp-unst-rey'
	CALL open_file2d(84,filedum,dum,6,tpts,jtot)
	End If
!
	If (iangle == 1) then
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='<greek>b</greek>'
	dum(4)='<greek>a</greek>'
	dum(5)='<greek>e</greek>'
	dum(6)='<greek>g</greek>'
	dum(7)='a<sub>1</sub>'
	filedum='angle-unst-rey'
	CALL open_file2d(403,filedum,dum,7,tpts,jtot)
	End If
!
!
	If (itke == 1) then
	    DO i=1,3
	        Write(three,'(I3)') i
		one=TRIM(adjustl(three))
	        dum(1)=re_val
	        dum(2)=y_val
		dum(3)='P<sub>'//one//one//'</sub>'
		dum(4)='T<sub>'//one//one//'</sub>'
		dum(5)='<greek>P</greek><sub>'//one//one//'</sub>'
		dum(6)='<greek>F</greek><sub>'//one//one//'</sub>'
		dum(7)='D<sub>'//one//one//'</sub>'
		dum(8)='<greek>e</greek><sub>'//one//one//'</sub>'
		filedum='tke'//uvwp(i)//uvwp(i)//'-unst-rey'
		CALL open_file2d(85+i,filedum,dum,8,tpts,jtot)
	        Write(three,'(I3)') i
		one=TRIM(adjustl(three))
	        dum(1)=re_val
	        dum(2)=y_val
		dum(3)='P<sub>'//one//one//'</sub>'
		dum(4)='T<sub>'//one//one//'</sub>'
		dum(5)='<greek>P</greek><sub>'//one//one//'</sub>'
		dum(6)='<greek>F</greek><sub>'//one//one//'</sub>'
		dum(7)='D<sub>'//one//one//'</sub>'
		dum(8)='<greek>e</greek><sub>'//one//one//'</sub>'
		filedum='tke'//uvwp(i)//uvwp(i)//'p-unst-rey'
		CALL open_file2d(185+i,filedum,dum,8,tpts,jtot)
	    END DO
!
	    dum(1)=re_val
	    dum(2)=y_val
	    dum(3)='P<sub>12</sub>'
	    dum(4)='T<sub>12</sub>'
	    dum(5)='<greek>P</greek><sub>12</sub>'
	    dum(6)='<greek>F</greek><sub>12</sub>'
	    dum(7)='D<sub>12</sub>'
	    dum(8)='<greek>e</greek><sub>12</sub>'
	    filedum='tkeuv-unst-rey'
	    CALL open_file2d(89,filedum,dum,8,tpts,jtot)
!
	    dum(1)=re_val
	    dum(2)=y_val
	    dum(3)='P<sub>12</sub>'
	    dum(4)='T<sub>12</sub>'
	    dum(5)='<greek>P</greek><sub>12</sub>'
	    dum(6)='<greek>F</greek><sub>12</sub>'
	    dum(7)='D<sub>12</sub>'
	    dum(8)='<greek>e</greek><sub>12</sub>'
	    filedum='tkeuvp-unst-rey'
	    CALL open_file2d(189,filedum,dum,8,tpts,jtot)
!
	    dum(1)=re_val
	    dum(2)=y_val
	    dum(3)='P<sub>k</sub>'
	    dum(4)='T<sub>k</sub>'
	    dum(5)='<greek>P</greek><sub>k</sub>'
	    dum(6)='<greek>F</greek><sub>k</sub>'
	    dum(7)='D<sub>k</sub>'
	    dum(8)='<greek>e</greek><sub>k</sub>'
	    filedum='tkek-unst-rey'
	    CALL open_file2d(92,filedum,dum,8,tpts,jtot)
!
	    dum(1)=re_val
	    dum(2)=y_val
	    dum(3)='P<sub>k</sub>'
	    dum(4)='T<sub>k</sub>'
	    dum(5)='<greek>P</greek><sub>k</sub>'
	    dum(6)='<greek>F</greek><sub>k</sub>'
	    dum(7)='D<sub>k</sub>'
	    dum(8)='<greek>e</greek><sub>k</sub>'
	    filedum='tkekp-unst-rey'
	    CALL open_file2d(192,filedum,dum,8,tpts,jtot)
!
	End If
!
	If (ianis==1) then
	    dum(1)=re_val
	    dum(2)=y_val
	    dum(3)='b<sub>11</sub>'
	    dum(4)='b<sub>22</sub>'
	    dum(5)='b<sub>33</sub>'
	    dum(6)='b<sub>12</sub>'
	    dum(7)='II'
	    dum(8)='III'
	    dum(9)='F'
	    dum(10)='G'
	    filedum='anis-unst-rey'
	    CALL open_file2d(91,filedum,dum,10,tpts,jtot)
	End If
!
	If (iquad == 1) then
	dum(1)=re_val
	dum(2)='yc'
	dum(3)='Q1'
	dum(4)='Q2'
	dum(5)='Q3'
	dum(6)='Q4'
	dum(7)='uv'
	dum(8)='vw'
	dum(9)='uw'
	filedum='quad-unst-rey'
	CALL open_file2d(98,filedum,dum,9,tpts,ny)
!
	dum(1)=re_val
	dum(2)='yc'
	dum(3)='Q1'
	dum(4)='Q2'
	dum(5)='Q3'
	dum(6)='Q4'
	dum(7)='uv'
	dum(8)='vw'
	dum(9)='uw'
	filedum='quadp-unst-rey'
	CALL open_file2d(198,filedum,dum,9,tpts,ny)
	End If
!
	If (ivort == 1) then
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='w<sub>x</sub>-mean'
	dum(4)='w<sub>y</sub>-mean'
	dum(5)='w<sub>z</sub>-mean'
	filedum='vortmean-unst-rey'
	CALL open_file2d(148,filedum,dum,5,tpts,jtot)

	dum(1)=re_val
	dum(2)=y_val
	dum(3)='w<sub>x</sub><sup>+</sup>-mean'
	dum(4)='w<sub>y</sub><sup>+</sup>-mean'
	dum(5)='w<sub>z</sub><sup>+</sup>-mean'
	filedum='vortmeanp-unst-rey'
	CALL open_file2d(149,filedum,dum,5,tpts,jtot)

	dum(1)=re_val
	dum(2)=y_val
	dum(3)='w<sub>x</sub>-rms'
	dum(4)='w<sub>y</sub>-rms'
	dum(5)='w<sub>z</sub>-rms'
	filedum='vortrms-unst-rey'
	CALL open_file2d(150,filedum,dum,5,tpts,jtot)

	dum(1)=re_val
	dum(2)=y_val
	dum(3)='w<sub>x</sub><sup>+</sup>-rms'
	dum(4)='w<sub>y</sub><sup>+</sup>-rms'
	dum(5)='w<sub>z</sub><sup>+</sup>-rms'
	filedum='vortrmsp-unst-rey'
	CALL open_file2d(151,filedum,dum,5,tpts,jtot)
	End If
!
	If (if_les == 1) Then
	 Open (101,file='nut-rey-unst-2d.dat')
!
        Write(101,*) 'title="nut-rey-unst"'
        Write(101,'(A140)') 'variables="Re","yc","<greek>n</greek><sub>SGS", "<greek>n</greek><sub>SGS</sub><sup>+","<greek>n</greek><sub>SGSi</sub><sup>+","c_s"'
        Write(101,*) 'zone t="z1",i=',tpts,',j=',jtot,'f=point'
	End If
!
	dum(1)=re_val
	dum(2)=y_val
	dum(3)='u<sub>mean</sub>'
	dum(4)='u<sub>rms</sub>'
	dum(5)='v<sub>rms</sub>'
	dum(6)='w<sub>rms</sub>'
	dum(7)='p<sub>rms</sub>'
	filedum='rates-unst-rey'
	CALL open_file2d(102,filedum,dum,7,tpts-1,jtot)

!	If (re_print==0) then; re_val='Re'
!	Else If (re_print==1) then; re_val='Re<sub><greek>t</greek></sub>'
!	Else; re_val='time'
!	End If 
	END SUBROUTINE open_files_unst
!*******************************************************************************
	SUBROUTINE open_1dpts_unst 
!		
	CHARACTER (LEN=100), DIMENSION(29) :: dum
	CHARACTER :: filedum*40,one*1,three*3,re_val*30,y_val*15
	INTEGER :: i,j,p
	INTEGER, DIMENSION(:), ALLOCATABLE :: twrite
	If (re_print==0) then; re_val='Re'
	Else If (re_print==1) then; re_val='Re<sub><greek>t</greek></sub>'
	Else; re_val='time'
	End If 
	If (y_print==0) then; y_val='yc'
	Else; y_val='yc<sup>+</sup>'
	End If
!
    If ((uns1d == 1).or.(uns1d == 3)) then
!
	dum(1)=re_val
	dum(2)='ny/32'
	dum(3)='ny/16'
	dum(4)='ny/8'
	dum(5)='3*ny/16'
	dum(6)='ny/4'
	dum(7)='3*ny/8'
	dum(8)='ny/2'
	dum(9)='time2'
!
!	Non wall units
	filedum='umean-unst-1d-t'
	CALL open_file1d(201,filedum,dum,9,tpts)
!
	If (iflucs == 1) then
	   filedum='urms-unst-1d-t'
	   CALL open_file1d(202,filedum,dum,9,tpts)
	   filedum='vrms-unst-1d-t'
	   CALL open_file1d(203,filedum,dum,9,tpts)
	   filedum='wrms-unst-1d-t'
	   CALL open_file1d(204,filedum,dum,9,tpts)
!
	   filedum='uskew-unst-1d-t'
	   CALL open_file1d(218,filedum,dum,9,tpts)
	   filedum='vskew-unst-1d-t'
	   CALL open_file1d(219,filedum,dum,9,tpts)
	   filedum='wskew-unst-1d-t'
	   CALL open_file1d(220,filedum,dum,9,tpts)
!
	   filedum='ukurt-unst-1d-t'
	   CALL open_file1d(221,filedum,dum,9,tpts)
	   filedum='vkurt-unst-1d-t'
	   CALL open_file1d(222,filedum,dum,9,tpts)
	   filedum='wkurt-unst-1d-t'
	   CALL open_file1d(223,filedum,dum,9,tpts)
	End If
!
	If (ishear == 1) then
	   filedum='shear-unst-1d-t'
	   CALL open_file1d(205,filedum,dum,9,tpts)
	End If
!
	If (iangle == 1) then
	   filedum='angle-unst-1d-t'
	   CALL open_file1d(404,filedum,dum,9,tpts)
	End If
!
	If (ivort == 1) then
	   filedum='vortx-unst-1d-t'
	   CALL open_file1d(206,filedum,dum,9,tpts)
	   filedum='vorty-unst-1d-t'
	   CALL open_file1d(207,filedum,dum,9,tpts)
	   filedum='vortz-unst-1d-t'
	   CALL open_file1d(208,filedum,dum,9,tpts)
	End If
!
	If (ianis == 1) then
	   filedum='anisb11-unst-1d-t'
	   CALL open_file1d(209,filedum,dum,9,tpts)
	   filedum='anisb22-unst-1d-t'
	   CALL open_file1d(210,filedum,dum,9,tpts)
	   filedum='anisb33-unst-1d-t'
	   CALL open_file1d(211,filedum,dum,9,tpts)
	   filedum='anisb12-unst-1d-t'
	   CALL open_file1d(212,filedum,dum,9,tpts)
	   filedum='anisII-unst-1d-t'
	   CALL open_file1d(213,filedum,dum,9,tpts)
	   filedum='anisIII-unst-1d-t'
	   CALL open_file1d(214,filedum,dum,9,tpts)
	   filedum='anisF-unst-1d-t'
	   CALL open_file1d(215,filedum,dum,9,tpts)
	   filedum='anisG-unst-1d-t'
	   CALL open_file1d(216,filedum,dum,9,tpts)
	End If
!
	If (if_les == 1) then
	   filedum='cles-unst-1d-t'
	   CALL open_file1d(217,filedum,dum,9,tpts)
	End If
!
	If (itke == 1) then
	    filedum='tkeuu-P-unst-1d-t'
	    CALL open_file1d(218,filedum,dum,9,tpts)
	    filedum='tkeuu-T-unst-1d-t'
	    CALL open_file1d(219,filedum,dum,9,tpts)
	    filedum='tkeuu-Pi-unst-1d-t'
	    CALL open_file1d(220,filedum,dum,9,tpts)
	    filedum='tkeuu-Phi-unst-1d-t'
	    CALL open_file1d(221,filedum,dum,9,tpts)
	    filedum='tkeuu-D-unst-1d-t'
	    CALL open_file1d(222,filedum,dum,9,tpts)
	    filedum='tkeuu-Epsilon-unst-1d-t'
	    CALL open_file1d(223,filedum,dum,9,tpts)
!
	    filedum='tkevv-P-unst-1d-t'
	    CALL open_file1d(224,filedum,dum,9,tpts)
	    filedum='tkevv-T-unst-1d-t'
	    CALL open_file1d(225,filedum,dum,9,tpts)
	    filedum='tkevv-Pi-unst-1d-t'
	    CALL open_file1d(226,filedum,dum,9,tpts)
	    filedum='tkevv-Phi-unst-1d-t'
	    CALL open_file1d(227,filedum,dum,9,tpts)
	    filedum='tkevv-D-unst-1d-t'
	    CALL open_file1d(228,filedum,dum,9,tpts)
	    filedum='tkevv-Epsilon-unst-1d-t'
	    CALL open_file1d(229,filedum,dum,9,tpts)
!
	    filedum='tkeww-P-unst-1d-t'
	    CALL open_file1d(230,filedum,dum,9,tpts)
	    filedum='tkeww-T-unst-1d-t'
	    CALL open_file1d(231,filedum,dum,9,tpts)
	    filedum='tkeww-Pi-unst-1d-t'
	    CALL open_file1d(232,filedum,dum,9,tpts)
	    filedum='tkeww-Phi-unst-1d-t'
	    CALL open_file1d(233,filedum,dum,9,tpts)
	    filedum='tkeww-D-unst-1d-t'
	    CALL open_file1d(234,filedum,dum,9,tpts)
	    filedum='tkeww-Epsilon-unst-1d-t'
	    CALL open_file1d(235,filedum,dum,9,tpts)
!
	    filedum='tkeuv-P-unst-1d-t'
	    CALL open_file1d(236,filedum,dum,9,tpts)
	    filedum='tkeuv-T-unst-1d-t'
	    CALL open_file1d(237,filedum,dum,9,tpts)
	    filedum='tkeuv-Pi-unst-1d-t'
	    CALL open_file1d(238,filedum,dum,9,tpts)
	    filedum='tkeuv-Phi-unst-1d-t'
	    CALL open_file1d(239,filedum,dum,9,tpts)
	    filedum='tkeuv-D-unst-1d-t'
	    CALL open_file1d(240,filedum,dum,9,tpts)
	    filedum='tkeuv-Epsilon-unst-1d-t'
	    CALL open_file1d(241,filedum,dum,9,tpts)
	End If
!
!	Wall units
	filedum='umeanp-unst-1d-t'
	CALL open_file1d(251,filedum,dum,9,tpts)
!
	If (iflucs == 1) then
	   filedum='urmsp-unst-1d-t'
	   CALL open_file1d(252,filedum,dum,9,tpts)
	   filedum='vrmsp-unst-1d-t'
	   CALL open_file1d(253,filedum,dum,9,tpts)
	   filedum='wrmsp-unst-1d-t'
	   CALL open_file1d(254,filedum,dum,9,tpts)
	End If
!
	If (ishear == 1) then
	   filedum='shearp-unst-1d-t'
	   CALL open_file1d(255,filedum,dum,9,tpts)
	End If
!
	If (ivort == 1) then
	   filedum='vortxp-unst-1d-t'
	   CALL open_file1d(256,filedum,dum,9,tpts)
	   filedum='vortyp-unst-1d-t'
	   CALL open_file1d(257,filedum,dum,9,tpts)
	   filedum='vortzp-unst-1d-t'
	   CALL open_file1d(258,filedum,dum,9,tpts)
	End If
!
	If (if_les == 1) then
	   filedum='clesp-unst-1d-t'
	   CALL open_file1d(267,filedum,dum,9,tpts)
	End If
!
	If (itke == 1) then
	    filedum='tkeuup-P-unst-1d-t'
	    CALL open_file1d(268,filedum,dum,9,tpts)
	    filedum='tkeuup-T-unst-1d-t'
	    CALL open_file1d(269,filedum,dum,9,tpts)
	    filedum='tkeuup-Pi-unst-1d-t'
	    CALL open_file1d(270,filedum,dum,9,tpts)
	    filedum='tkeuup-Phi-unst-1d-t'
	    CALL open_file1d(271,filedum,dum,9,tpts)
	    filedum='tkeuup-D-unst-1d-t'
	    CALL open_file1d(272,filedum,dum,9,tpts)
	    filedum='tkeuup-Epsilon-unst-1d-t'
	    CALL open_file1d(273,filedum,dum,9,tpts)
!
	    filedum='tkevvp-P-unst-1d-t'
	    CALL open_file1d(274,filedum,dum,9,tpts)
	    filedum='tkevvp-T-unst-1d-t'
	    CALL open_file1d(275,filedum,dum,9,tpts)
	    filedum='tkevvp-Pi-unst-1d-t'
	    CALL open_file1d(276,filedum,dum,9,tpts)
	    filedum='tkevvp-Phi-unst-1d-t'
	    CALL open_file1d(277,filedum,dum,9,tpts)
	    filedum='tkevvp-D-unst-1d-t'
	    CALL open_file1d(278,filedum,dum,9,tpts)
	    filedum='tkevvp-Epsilon-unst-1d-t'
	    CALL open_file1d(279,filedum,dum,9,tpts)
!
	    filedum='tkewwp-P-unst-1d-t'
	    CALL open_file1d(280,filedum,dum,9,tpts)
	    filedum='tkewwp-T-unst-1d-t'
	    CALL open_file1d(281,filedum,dum,9,tpts)
	    filedum='tkewwp-Pi-unst-1d-t'
	    CALL open_file1d(282,filedum,dum,9,tpts)
	    filedum='tkewwp-Phi-unst-1d-t'
	    CALL open_file1d(283,filedum,dum,9,tpts)
	    filedum='tkewwp-D-unst-1d-t'
	    CALL open_file1d(284,filedum,dum,9,tpts)
	    filedum='tkewwp-Epsilon-unst-1d-t'
	    CALL open_file1d(285,filedum,dum,9,tpts)
!
	    filedum='tkeuvp-P-unst-1d-t'
	    CALL open_file1d(286,filedum,dum,9,tpts)
	    filedum='tkeuvp-T-unst-1d-t'
	    CALL open_file1d(287,filedum,dum,9,tpts)
	    filedum='tkeuvp-Pi-unst-1d-t'
	    CALL open_file1d(288,filedum,dum,9,tpts)
	    filedum='tkeuvp-Phi-unst-1d-t'
	    CALL open_file1d(289,filedum,dum,9,tpts)
	    filedum='tkeuvp-D-unst-1d-t'
	    CALL open_file1d(290,filedum,dum,9,tpts)
	    filedum='tkeuvp-Epsilon-unst-1d-t'
	    CALL open_file1d(291,filedum,dum,9,tpts)
	End If
!
    End If
!
    If ((uns1d == 2).or.(uns1d == 3)) then
	ALLOCATE(twrite(ntsamps))
!	twrite(1)=Floor(0.015e0*(tpts-1))+1
!	Do p = 2,ntsamps
!	      twrite(p) = twrite(1)+Floor((0.0625e0*Real(p-1))*(tpts-1))+1
!	End Do
	Do p = 1,ntsamps
	  twrite(p) = ((p-1)*(tpts-1)/(ntsamps-1)) + 1
	End Do
!
!	For non wall units
	dum(1)='yc'
	DO j=1,ntsamps
	    WRITE(dum(j+1),'(E10.5)') time_sa(twrite(j))
 	    dum(j+1) = 't='//TRIM(adjustl(dum(j+1)(1:10)))
	END DO
!
	filedum='umean-unst-1d-y'
	CALL open_file1d(301,filedum,dum,ntsamps+1,jtot)
	filedum='wmean-unst-1d-y'
	CALL open_file1d(450,filedum,dum,ntsamps+1,jtot)
	filedum='vmean-unst-1d-y'
	CALL open_file1d(451,filedum,dum,ntsamps+1,jtot)
!
	If (iflucs == 1) then
	   filedum='urms-unst-1d-y'
	   CALL open_file1d(302,filedum,dum,ntsamps+1,jtot)
	   filedum='vrms-unst-1d-y'
	   CALL open_file1d(303,filedum,dum,ntsamps+1,jtot)
	   filedum='wrms-unst-1d-y'
	   CALL open_file1d(304,filedum,dum,ntsamps+1,jtot)
	End If
!
	If (ishear == 1) then
	   filedum='shear-unst-1d-y'
	   CALL open_file1d(305,filedum,dum,ntsamps+1,jtot)
	End If
!
	If (iangle == 1) then
	   filedum='angle-unst-1d-y'
	   CALL open_file1d(405,filedum,dum,ntsamps+1,jtot)
	End If
!
	If (ivort == 1) then
	   filedum='vortxmean-unst-1d-y'
	   CALL open_file1d(398,filedum,dum,ntsamps+1,jtot)
	   filedum='vortx-unst-1d-y'
	   CALL open_file1d(306,filedum,dum,ntsamps+1,jtot)
	   filedum='vorty-unst-1d-y'
	   CALL open_file1d(307,filedum,dum,ntsamps+1,jtot)
	   filedum='vortz-unst-1d-y'
	   CALL open_file1d(308,filedum,dum,ntsamps+1,jtot)
	End If
!
	If (if_les == 1) then
	   filedum='cles-unst-1d-y'
	   CALL open_file1d(317,filedum,dum,ntsamps+1,jtot)
	End If
!
	If (itke == 1) then
	    filedum='tkeuu-P-unst-1d-y'
	    CALL open_file1d(318,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuu-T-unst-1d-y'
	    CALL open_file1d(319,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuu-Pi-unst-1d-y'
	    CALL open_file1d(320,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuu-Phi-unst-1d-y'
	    CALL open_file1d(321,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuu-D-unst-1d-y'
	    CALL open_file1d(322,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuu-Epsilon-unst-1d-y'
	    CALL open_file1d(323,filedum,dum,ntsamps+1,jtot)
!
	    filedum='tkevv-P-unst-1d-y'
	    CALL open_file1d(324,filedum,dum,ntsamps+1,jtot)
	    filedum='tkevv-T-unst-1d-y'
	    CALL open_file1d(325,filedum,dum,ntsamps+1,jtot)
	    filedum='tkevv-Pi-unst-1d-y'
	    CALL open_file1d(326,filedum,dum,ntsamps+1,jtot)
	    filedum='tkevv-Phi-unst-1d-y'
	    CALL open_file1d(327,filedum,dum,ntsamps+1,jtot)
	    filedum='tkevv-D-unst-1d-y'
	    CALL open_file1d(328,filedum,dum,ntsamps+1,jtot)
	    filedum='tkevv-Epsilon-unst-1d-y'
	    CALL open_file1d(329,filedum,dum,ntsamps+1,jtot)
!
	    filedum='tkeww-P-unst-1d-y'
	    CALL open_file1d(330,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeww-T-unst-1d-y'
	    CALL open_file1d(331,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeww-Pi-unst-1d-y'
	    CALL open_file1d(332,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeww-Phi-unst-1d-y'
	    CALL open_file1d(333,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeww-D-unst-1d-y'
	    CALL open_file1d(334,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeww-Epsilon-unst-1d-y'
	    CALL open_file1d(335,filedum,dum,ntsamps+1,jtot)
!
	    filedum='tkeuv-P-unst-1d-y'
	    CALL open_file1d(336,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuv-T-unst-1d-y'
	    CALL open_file1d(337,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuv-Pi-unst-1d-y'
	    CALL open_file1d(338,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuv-Phi-unst-1d-y'
	    CALL open_file1d(339,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuv-D-unst-1d-y'
	    CALL open_file1d(340,filedum,dum,ntsamps+1,jtot)
	    filedum='tkeuv-Epsilon-unst-1d-y'
	    CALL open_file1d(341,filedum,dum,ntsamps+1,jtot)
	End If
!
!	For wall units
	dum(1)='yc'
	DO j=1,ntsamps
	    WRITE(dum(2*j),'(F10.3)') time_sa(twrite(j))
	    WRITE(dum(2*j),*) 'yc<sup>+</sup> at '//TRIM(adjustl(dum(2*j)(1:10)))
	    WRITE(dum(2*j+1),*) time_sa(twrite(j))
 	    dum(2*j+1) = 't='//TRIM(adjustl(dum(2*j+1)(1:10)))
	END DO
!
	filedum='umeanp-unst-1d-y'
	CALL open_file1d(351,filedum,dum,2*ntsamps+1,jtot)
!
	If (iflucs == 1) then
	   filedum='urmsp-unst-1d-y'
	   CALL open_file1d(352,filedum,dum,2*ntsamps+1,jtot)
	   filedum='vrmsp-unst-1d-y'
	   CALL open_file1d(353,filedum,dum,2*ntsamps+1,jtot)
	   filedum='wrmsp-unst-1d-y'
	   CALL open_file1d(354,filedum,dum,2*ntsamps+1,jtot)
!
	   filedum='uskewp-unst-1d-y'
	   CALL open_file1d(392,filedum,dum,2*ntsamps+1,jtot)
	   filedum='vskewp-unst-1d-y'
	   CALL open_file1d(393,filedum,dum,2*ntsamps+1,jtot)
	   filedum='wskewp-unst-1d-y'
	   CALL open_file1d(394,filedum,dum,2*ntsamps+1,jtot)
!
	   filedum='ukurtp-unst-1d-y'
	   CALL open_file1d(395,filedum,dum,2*ntsamps+1,jtot)
	   filedum='vkurtp-unst-1d-y'
	   CALL open_file1d(396,filedum,dum,2*ntsamps+1,jtot)
	   filedum='wkurtp-unst-1d-y'
	   CALL open_file1d(397,filedum,dum,2*ntsamps+1,jtot)
	End If
!
	If (ishear == 1) then
	   filedum='shearp-unst-1d-y'
	   CALL open_file1d(355,filedum,dum,2*ntsamps+1,jtot)
	End If
!
	If (ivort == 1) then
	   filedum='vortxmeanp-unst-1d-y'
	   CALL open_file1d(399,filedum,dum,2*ntsamps+1,jtot)
	   filedum='vortxp-unst-1d-y'
	   CALL open_file1d(356,filedum,dum,2*ntsamps+1,jtot)
	   filedum='vortyp-unst-1d-y'
	   CALL open_file1d(357,filedum,dum,2*ntsamps+1,jtot)
	   filedum='vortzp-unst-1d-y'
	   CALL open_file1d(358,filedum,dum,2*ntsamps+1,jtot)
	End If
!
	If (ianis == 1) then
	   filedum='anisb11p-unst-1d-y'
	   CALL open_file1d(359,filedum,dum,2*ntsamps+1,jtot)
	   filedum='anisb22p-unst-1d-y'
	   CALL open_file1d(360,filedum,dum,2*ntsamps+1,jtot)
	   filedum='anisb33p-unst-1d-y'
	   CALL open_file1d(361,filedum,dum,2*ntsamps+1,jtot)
	   filedum='anisb12p-unst-1d-y'
	   CALL open_file1d(362,filedum,dum,2*ntsamps+1,jtot)
	   filedum='anisIIp-unst-1d-y'
	   CALL open_file1d(363,filedum,dum,2*ntsamps+1,jtot)
	   filedum='anisIIIp-unst-1d-y'
	   CALL open_file1d(364,filedum,dum,2*ntsamps+1,jtot)
	   filedum='anisFp-unst-1d-y'
	   CALL open_file1d(365,filedum,dum,2*ntsamps+1,jtot)
	   filedum='anisGp-unst-1d-y'
	   CALL open_file1d(366,filedum,dum,2*ntsamps+1,jtot)
	End If
!
	If (if_les == 1) then
	   filedum='clesp-unst-1d-y'
	   CALL open_file1d(367,filedum,dum,2*ntsamps+1,jtot)
	End If
!
	If (itke == 1) then
	    filedum='tkeuup-P-unst-1d-y'
	    CALL open_file1d(368,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuup-T-unst-1d-y'
	    CALL open_file1d(369,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuup-Pi-unst-1d-y'
	    CALL open_file1d(370,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuup-Phi-unst-1d-y'
	    CALL open_file1d(371,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuup-D-unst-1d-y'
	    CALL open_file1d(372,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuup-Epsilon-unst-1d-y'
	    CALL open_file1d(373,filedum,dum,2*ntsamps+1,jtot)
!
	    filedum='tkevvp-P-unst-1d-y'
	    CALL open_file1d(374,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkevvp-T-unst-1d-y'
	    CALL open_file1d(375,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkevvp-Pi-unst-1d-y'
	    CALL open_file1d(376,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkevvp-Phi-unst-1d-y'
	    CALL open_file1d(377,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkevvp-D-unst-1d-y'
	    CALL open_file1d(378,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkevvp-Epsilon-unst-1d-y'
	    CALL open_file1d(379,filedum,dum,2*ntsamps+1,jtot)
!
	    filedum='tkewwp-P-unst-1d-y'
	    CALL open_file1d(380,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkewwp-T-unst-1d-y'
	    CALL open_file1d(381,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkewwp-Pi-unst-1d-y'
	    CALL open_file1d(382,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkewwp-Phi-unst-1d-y'
	    CALL open_file1d(383,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkewwp-D-unst-1d-y'
	    CALL open_file1d(384,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkewwp-Epsilon-unst-1d-y'
	    CALL open_file1d(385,filedum,dum,2*ntsamps+1,jtot)
!
	    filedum='tkeuvp-P-unst-1d-y'
	    CALL open_file1d(386,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuvp-T-unst-1d-y'
	    CALL open_file1d(387,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuvp-Pi-unst-1d-y'
	    CALL open_file1d(388,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuvp-Phi-unst-1d-y'
	    CALL open_file1d(389,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuvp-D-unst-1d-y'
	    CALL open_file1d(390,filedum,dum,2*ntsamps+1,jtot)
	    filedum='tkeuvp-Epsilon-unst-1d-y'
	    CALL open_file1d(391,filedum,dum,2*ntsamps+1,jtot)
	End If
!
    End If
!
	END SUBROUTINE open_1dpts_unst
!*******************************************************************************
	SUBROUTINE close_files !Needs to be sorted out
!
	Close(35);Close(36);Close(37);Close(38)
!
	If (isteady == 1) then
	Close(41)
	If (iflucs == 1) Close(42)
	If (ivort == 1) Close(43)
	If (iquad == 1) Close(47)
 	Close(48);
	If (itke == 1) then
	Close(49);Close(50);Close(51);Close(52)
	End If
	If (icorr == 1) then
	Close(56);Close(57);Close(58);Close(59)
	Close(512);Close(513)
	End If
	If (ishear == 1) then
	Close(44);Close(54)
	End If
	If (iangle == 1) then
	Close(401)
	End If
!
	Else If (isteady == 0) then
!
	Close(70);Close(71);Close(72);Close(73);Close(74);Close(75);Close(76)
	Close(77);Close(101);Close(102);Close(103);Close(104);Close(105);
	Close(106)
	If (iflucs == 1) then
	Close(78);Close(79);Close(80);Close(81);Close(82);Close(83);Close(84)
	Close(85);
	End If
	If (itke == 1) then
	Close(86);Close(87);Close(88);Close(89); Close(94);Close(95);Close(96)
	Close(97)
	End If

	Close(90);Close(91);Close(92);Close(93);

	If (iquad == 1) then
	Close(98);Close(99)
	End If
	If (ishear == 1) then
	Close(83);Close(84);Close(205);Close(255);Close(305);Close(355)
	End If
	If (iangle == 1) then
	Close(403);Close(404);Close(405)
	End If
!
	End If
!
	END SUBROUTINE close_files
!*******************************************************************************
	SUBROUTINE write_0d_data
	INTEGER :: i,j2,j,m
	INTEGER, DIMENSION (7) :: jins
        REAL(mytype), DIMENSION(6) :: rmsofrms
        REAL(mytype), DIMENSION(3) :: avgflucs
        REAL(mytype), DIMENSION(tpts) :: val,sqr    
!	This writes the 1d time averaged files
	DO i=1,tpts
		Write(35,'(2E15.7)')time_sa(i),utau_sa(i)
		Write(36,'(3E15.7)')time_sa(i),utau_sa(i)/nu, re_sa(i)
		Write(37,'(2E15.7)')time_sa(i),dpdx_sa(i)
	END DO
	Open (381,file='re-check.dat')
	DO i=1,tpts
	Write (381,*)1.0e0/re_sa(i)/utau_sa(i)**2,re_sa(i),utau_sa(i)**2
	Enddo

	!DO j=1,ny
         !       val=sqrt(abs(flucs1d_sa(j,1,2,:)))
         !       sqr=val**2
         !       rmsofrms(1)=((SUM(sqr)/real(tpts))-(SUM(val)/real(tpts))**2)*3150/200
         !       val=sqrt(abs(flucs1d_sa(j,2,2,:)))
         !!       sqr=val**2
         !       rmsofrms(2)=((SUM(sqr)/real(tpts))-(SUM(val)/real(tpts))**2)*3150/200
         !       val=sqrt(abs(flucs1d_sa(j,3,2,:)))
         !       sqr=val**2
         !       rmsofrms(3)=((SUM(sqr)/real(tpts))-(SUM(val)/real(tpts))**2)*3150/200

        !        val=omega1d_sa(j,1,:)
        !        sqr=val**2
        !        rmsofrms(4)=((SUM(sqr)/real(tpts))-(SUM(val)/real(tpts))**2)*3150/200/200
        !        val=omega1d_sa(j,2,:)
        !        sqr=val**2
        !        rmsofrms(5)=((SUM(sqr)/real(tpts))-(SUM(val)/real(tpts))**2)*3150/200/200
        !        val=omega1d_sa(j,3,:)
        !        sqr=val**2
        !        rmsofrms(6)=((SUM(sqr)/real(tpts))-(SUM(val)/real(tpts))**2)*3150/200/200
	!	Write(38,'(7E15.7)')yc(j),(rmsofrms(m),m=1,6)

        !        avgflucs(1)=sqrt(abs(SUM(flucs1d_sa(j,1,2,:))/real(tpts)))
        !        avgflucs(2)=sqrt(abs(SUM(flucs1d_sa(j,2,2,:))/real(tpts)))
        !        avgflucs(3)=sqrt(abs(SUM(flucs1d_sa(j,3,2,:))/real(tpts)))
	!	Write(39,'(4E15.7)')yc(j),(avgflucs(m),m=1,3)                
	!END DO

	Close (381)
	END SUBROUTINE write_0d_data
!*******************************************************************************
	SUBROUTINE write_1d_avg
	INTEGER :: i,j,k,m,n

!	This writes the 1d time averaged files
	DO j =jmin,jmax
	    If (iflucs == 1) then   
	       Write(41,'(6E15.7)')ycplus(j),(flucs(j,m,1),m=1,4),yc(j)
	       Write(48,'(6E15.7)')ycplus(j),(flucs(j,m,1)/utau_avg_sa,m=1,3),flucs(j,m,1)/utau_avg_sa**2/rho,yc(j) 
	       Write(42,'(6E15.7)')ycplus(j),(sqrt(abs(flucs(j,m,2))),m=1,4),yc(j)
	       Write(53,'(6E15.7)')ycplus(j),(sqrt(abs(flucs(j,m,2)))/utau_avg_sa,m=1,3),sqrt(abs(flucs(j,4,2)))/utau_avg_sa**2/rho,yc(j)
               Write(45,'(6E15.7)')ycplus(j),(flucs(j,m,3)/abs(flucs(j,m,2))**1.5, m=1,4),yc(j)
               Write(46,'(6E15.7)')ycplus(j),((flucs(j,m,4)/flucs(j,m,2)**2.0)-3.0, m=1,4),yc(j)
	    End If
	    If (ishear == 1) then
		Write(44,'(6E15.7)')ycplus(j),(shear(j,m),m=1,4),yc(j)
		Write(54,'(6E15.7)')ycplus(j),(shear(j,m)/utau_avg_sa**2,m=1,4),yc(j)
	    End If
	    If (iangle == 1) then
		Write(401,'(7E15.7)')ycplus(j),(angle(j,m),m=1,5),yc(j)
	    End If
!
            If (ivort == 1) then
                Write(43,'(8E15.7)') ycplus(j),omega_avg(j,1),omega_avg(j,2),omega_avg(j,3), omega_avg(j,4),omega_avg(j,5),omega_avg(j,6),yc(j)
                Write(55,'(8E15.7)') ycplus(j),omega_avg(j,1)*nu/utau_avg_sa**2,omega_avg(j,2)*nu/utau_avg_sa**2,   &
                                     omega_avg(j,3)*nu/utau_avg_sa**2, omega_avg(j,4)*nu/utau_avg_sa**2,omega_avg(j,5)*nu/utau_avg_sa**2,omega_avg(j,6)*nu/utau_avg_sa**2,yc(j)
	    End If
	END DO

  	DO j = jmin,jmax
            If (iquad == 1) then
		Write(47,'(9E15.7)') ycplus(j),(quad_avg(j,m),m=1,6),quad_avg(j,7)-stat_avg(j,1,1)*stat_avg(j,3,1),yc(j)
!		Write(47,'(9E15.7)') ycplus(j),(quad_avg(j,m)-stat_avg(j,1,1)*Sqrt(flucs(j,2,2)),m=1,2),(quad_avg(j,m)+stat_avg(j,1,1)*Sqrt(flucs(j,2,2)),m=3,4),(quad_avg(j,m),m=5,6),quad_avg(j,7)-stat_avg(j,1,1)*stat_avg(j,3,1),yc(j)
	    End If
	    If (itke == 1) then
	        Write(49,'(8E15.7)') ycplus(j),(kinen(j,1,m),m=1,6),yc(j)
	        Write(50,'(8E15.7)') ycplus(j),(kinen(j,2,m),m=1,6),yc(j)
	        Write(51,'(8E15.7)') ycplus(j),(kinen(j,3,m),m=1,6),yc(j)
	        Write(52,'(8E15.7)') ycplus(j),(kinen(j,4,m),m=1,6),yc(j)
	        Write(66,'(8E15.7)') ycplus(j),((kinen(j,1,m)+kinen(j,2,m)+kinen(j,3,m))/2.0,m=1,6),yc(j)
!
	        Write(62,'(8E15.7)') ycplus(j),(kinen(j,1,m)*nu/utau_avg_sa**4,m=1,6),yc(j)
	        Write(63,'(8E15.7)') ycplus(j),(kinen(j,2,m)*nu/utau_avg_sa**4,m=1,6),yc(j)
	        Write(64,'(8E15.7)') ycplus(j),(kinen(j,3,m)*nu/utau_avg_sa**4,m=1,6),yc(j)
	        Write(65,'(8E15.7)') ycplus(j),(kinen(j,4,m)*nu/utau_avg_sa**4,m=1,6),yc(j)
	        Write(67,'(8E15.7)') ycplus(j),((kinen(j,1,m)+kinen(j,2,m)+kinen(j,3,m))*nu/2.0/utau_avg_sa**4,m=1,6),yc(j)
	    End If
	    If (ianis == 1) Write(61,'(10E15.7)')ycplus(j),(anisb(j,m),m=1,4),anisII(j),anisIII(j), &
						anisf(j),anisg(j),yc(j)
   	END DO
!
	If(itrans==1) Call write_1d_avg_transform
!
	END SUBROUTINE write_1d_avg
!*******************************************************************************
        SUBROUTINE write_1d_avg_transform
        Real(mytype) :: urms,vrms,wrms,u2,v2,w2,uv,uw,vw,urot,wrot,umn,wmn
        Real(mytype), Allocatable, Dimension(:) :: urmsavg, vrmsavg, wrmsavg, &
	 uvavg, uwavg, vwavg, &
	 umnavg, wmnavg, urmsavg1, vrmsavg1, wrmsavg1, &
	 uvavg1, uwavg1, vwavg1, umnavg1, wmnavg1
        INTEGER :: i,j,k,m,n,iccum,isty
!       This writes the 1d time averaged files with coordinate transformation
! Write the statistics out
        Open(625, File='statistics.plt', status='UNKNOWN')
        Write(625,'(120A,120A,120A)') ' variables="y","U","W","u<sub>rms</sub>","v<sub>rms</sub>","w<sub>rms</sub>","uv","uw","vw",', &
          '"U1","W1","u<sub>rms</sub>1","v<sub>rms</sub>1","w<sub>rms</sub>1","uv1","uw1","vw1",', &
         '"y<sup>+</sup>"'
!
	isty=1
	Allocate(urmsavg(-1:ny+2),vrmsavg(-1:ny+2),wrmsavg(-1:ny+2), &
		uvavg(-1:ny+2),uwavg(-1:ny+2),vwavg(-1:ny+2),&
		urmsavg1(-1:ny+2),vrmsavg1(-1:ny+2),wrmsavg1(-1:ny+2), &
		uvavg1(-1:ny+2),uwavg1(-1:ny+2),vwavg1(-1:ny+2), &
		umnavg(-1:ny+2),wmnavg(-1:ny+2),umnavg1(-1:ny+2),wmnavg1(-1:ny+2))
	umnavg(:)=0.0e0;wmnavg(:)=0.0e0
	urmsavg(:)=0.0e0;vrmsavg(:)=0.0e0;wrmsavg(:)=0.0e0
	uvavg(:)=0.0e0;uwavg(:)=0.0e0;vwavg(:)=0.0e0
	umnavg1(:)=0.0e0;wmnavg1(:)=0.0e0
	urmsavg1(:)=0.0e0;vrmsavg1(:)=0.0e0;wrmsavg1(:)=0.0e0
	uvavg1(:)=0.0e0;uwavg1(:)=0.0e0;vwavg1(:)=0.0e0
!
	If(isty==2) Then; iccum=tpts; Else; iccum=1; End If
        Do i=1,iccum; Do j=jmin,jmax
	  If(isty==2) Then 
            umn=stat1d_sa(j,1,1,i);wmn=stat1d_sa(j,3,1,i)
            u2=stat1d_sa(j,1,2,i);v2=stat1d_sa(j,2,2,i);w2=stat1d_sa(j,3,2,i)
            uv=quad1d_sa(j,5,i);vw=quad1d_sa(j,6,i);uw=quad1d_sa(j,7,i)
	  Else
            umn=stat_avg(j,1,1);wmn=stat_avg(j,3,1)
            u2=stat_avg(j,1,2);v2=stat_avg(j,2,2);w2=stat_avg(j,3,2)
            uv=quad_avg(j,5);vw=quad_avg(j,6);uw=quad_avg(j,7)
	  End If
          urms=Sqrt(u2-umn**2);vrms=Sqrt(v2);wrms=Sqrt(w2-wmn**2)
!
          umnavg(j)=umnavg(j)+umn
	  wmnavg(j)=wmnavg(j)+wmn
          umnavg1(j)=umnavg1(j)+umn*cos(theta)-wmn*sin(theta)
	  wmnavg1(j)=wmnavg1(j)+umn*sin(theta)+wmn*cos(theta)
	  urmsavg(j)=urmsavg(j)+urms
	  vrmsavg(j)=vrmsavg(j)+vrms
	  wrmsavg(j)=wrmsavg(j)+wrms
	  uvavg(j)=uvavg(j)+uv
	  uwavg(j)=uwavg(j)+uw-umn*wmn
	  vwavg(j)=vwavg(j)+vw
	  urmsavg1(j)=urmsavg1(j)+Sqrt((urms*cos(theta))**2+(wrms*sin(theta))**2-(uw-umn*wmn)*sin(2.0e0*theta))
          vrmsavg1(j)=vrmsavg1(j)+vrms
          wrmsavg1(j)=wrmsavg1(j)+Sqrt((urms*sin(theta))**2+(wrms*cos(theta))**2+(uw-umn*wmn)*sin(2.0e0*theta))
	  uvavg1(j)=uvavg1(j)+uv*cos(theta)-vw*sin(theta)
	  uwavg1(j)=uwavg1(j)+(urms**2-wrms**2)*sin(theta)*cos(theta)+(uw-umn*wmn)*(cos(theta)**2-sin(theta)**2)
	  vwavg1(j)=vwavg1(j)+uv*sin(theta)+vw*cos(theta)
        End Do; End Do
! average within a period
	umnavg(:)=umnavg(:)/Real(iccum,kind=mytype)
	wmnavg(:)=wmnavg(:)/Real(iccum,kind=mytype)
	umnavg1(:)=umnavg1(:)/Real(iccum,kind=mytype)
	wmnavg1(:)=wmnavg1(:)/Real(iccum,kind=mytype)
	urmsavg(:)=urmsavg(:)/Real(iccum,kind=mytype)
	vrmsavg(:)=vrmsavg(:)/Real(iccum,kind=mytype)
	wrmsavg(:)=wrmsavg(:)/Real(iccum,kind=mytype)
	urmsavg1(:)=urmsavg1(:)/Real(iccum,kind=mytype)
	vrmsavg1(:)=vrmsavg1(:)/Real(iccum,kind=mytype)
	wrmsavg1(:)=wrmsavg1(:)/Real(iccum,kind=mytype)
	uvavg(:)=uvavg(:)/Real(iccum,kind=mytype)
	uwavg(:)=uwavg(:)/Real(iccum,kind=mytype)
	vwavg(:)=vwavg(:)/Real(iccum,kind=mytype)
	uvavg1(:)=uvavg1(:)/Real(iccum,kind=mytype)
	uwavg1(:)=uwavg1(:)/Real(iccum,kind=mytype)
	vwavg1(:)=vwavg1(:)/Real(iccum,kind=mytype)
!
	Do j=jmin,jmax
	  Write(625,'(18E20.10)') yc(j),umnavg(j),wmnavg(j),urmsavg(j),vrmsavg(j), &
	wrmsavg(j),uvavg(j),uwavg(j),vwavg(j), &
	umnavg1(j),wmnavg1(j),urmsavg1(j),vrmsavg1(j), &
	wrmsavg1(j),uvavg1(j),uwavg1(j),vwavg1(j),ycplus(j)
	End Do
        Close(625)
!
	Deallocate(urmsavg, vrmsavg, wrmsavg, &
	 uvavg, uwavg, vwavg, umnavg, wmnavg, &
	 urmsavg1, vrmsavg1, wrmsavg1, &
	 uvavg1, uwavg1, vwavg1, umnavg1, wmnavg1)
!
        END SUBROUTINE write_1d_avg_transform
!*******************************************************************************
	SUBROUTINE write_2d_avg
	INTEGER :: i,j,k,m,n
!	This writes the 1d time averaged files
!
	DO j=1,ny
	    DO i=1,nx
		Write(110,'(5E15.7)')xcplus(i),ycplus(j),(flucs2d_avg(i,j,m,1),m=1,3)
		Write(111,'(5E15.7)')xcplus(i),ycplus(j),(flucs2d_avg(i,j,m,1)*re_avg_sa/retau_avg_sa,m=1,3)
	    	If (iflucs == 1) then    
	       	    Write(112,'(6E15.7)')xcplus(i),ycplus(j),(sqrt(abs(flucs2d_avg(i,j,m,2))),m=1,4)
	            Write(113,'(6E15.7)')xcplus(i),ycplus(j),(sqrt(abs(flucs2d_avg(i,j,m,2)))*re_avg_sa/retau_avg_sa,m=1,4)
                    Write(114,'(6E15.7)')xcplus(i),ycplus(j),(flucs2d_avg(i,j,m,3)/(abs(flucs2d_avg(i,j,m,2))**1.5),m=1,4)
                    Write(115,'(6E15.7)')xcplus(i),ycplus(j),((flucs2d_avg(i,j,m,4)/(flucs2d_avg(i,j,m,2)**2.0))-3.0,m=1,4)
	        End If
	        If (ishear == 1) Write(116,'(6E15.7)')xcplus(i),ycplus(j),(shear2d_avg(i,j,m),m=1,4)
	        If (iangle == 1) Write(402,'(7E15.7)')xcplus(i),ycplus(j),(angle2d_avg(i,j,m),m=1,5)
	    END DO
	END DO
!
	END SUBROUTINE write_2d_avg
!*******************************************************************************
	SUBROUTINE write_1dpts_unst
	INTEGER :: i,j,m,fileval,j2,p,i2,nsnap_target,n
	REAL(mytype) :: re_val,y_val
	INTEGER, DIMENSION (7) :: jins
	INTEGER, DIMENSION (:), ALLOCATABLE :: twrite
	Real(mytype), DIMENSION (:), ALLOCATABLE :: snap_re_target
!
    If ((uns1d == 1).or.(uns1d == 3)) then
        jins(1) = ny/32
        jins(2) = ny/16
        jins(3) = ny/8
        jins(4) = 3*ny/16
        jins(5) = ny/4
        jins(6) = 3*ny/8
        jins(7) = ny/2
	Do i = 1,tpts
	    If (re_print==0) Then; re_val=re_sa(i)
	    Else If (re_print==1) Then; re_val=retau_sa(i)
	    Else; re_val=time_sa(i)
	    End If
	    Write(201,'(9E15.7)') re_val,(stat1d_sa(jins(j2),1,1,i), j2=1,7),time_sa(i)
	    Write(251,'(9E15.7)') re_val,(stat1d_sa(jins(j2),1,1,i)/utau_sa(i), j2=1,7),time_sa(i)
	    If (iflucs == 1) then
	    	Write(202,'(9E15.7)') re_val,(sqrt(abs(flucs1d_sa(jins(j2),1,2,i))), j2=1,7),time_sa(i)
	    	Write(252,'(9E15.7)') re_val,(sqrt(abs(flucs1d_sa(jins(j2),1,2,i)))/utau_sa(i), j2=1,7),time_sa(i)
	    	Write(203,'(9E15.7)') re_val,(sqrt(abs(flucs1d_sa(jins(j2),2,2,i))), j2=1,7),time_sa(i)
	    	Write(253,'(9E15.7)') re_val,(sqrt(abs(flucs1d_sa(jins(j2),2,2,i)))/utau_sa(i), j2=1,7),time_sa(i)
	    	Write(204,'(9E15.7)') re_val,(sqrt(abs(flucs1d_sa(jins(j2),3,2,i))), j2=1,7),time_sa(i)
	    	Write(254,'(9E15.7)') re_val,(sqrt(abs(flucs1d_sa(jins(j2),3,2,i)))/utau_sa(i), j2=1,7),time_sa(i)

	        Write(218,'(9E15.7)') re_val,(flucs1d_sa(jins(j2),1,3,i)/abs(flucs1d_sa(jins(j2),1,2,i))**1.5, j2=1,7),time_sa(i)
	        Write(219,'(9E15.7)') re_val,(flucs1d_sa(jins(j2),2,3,i)/abs(flucs1d_sa(jins(j2),2,2,i))**1.5, j2=1,7),time_sa(i)
	        Write(220,'(9E15.7)') re_val,(flucs1d_sa(jins(j2),3,3,i)/abs(flucs1d_sa(jins(j2),3,2,i))**1.5, j2=1,7),time_sa(i)

	        Write(221,'(9E15.7)') re_val,((flucs1d_sa(jins(j2),1,4,i)/(flucs1d_sa(jins(j2),1,2,i)**2.0)-3.0), j2=1,7),time_sa(i)
	        Write(222,'(9E15.7)') re_val,((flucs1d_sa(jins(j2),2,4,i)/(flucs1d_sa(jins(j2),2,2,i)**2.0)-3.0), j2=1,7),time_sa(i)	
	        Write(223,'(9E15.7)') re_val,((flucs1d_sa(jins(j2),3,4,i)/(flucs1d_sa(jins(j2),3,2,i)**2.0)-3.0), j2=1,7),time_sa(i)	
	    End If
	    If (ishear == 1) then
		Write(205,'(9E15.7)') re_val,(shear1d_sa(jins(j2),4,i), j2=1,7),time_sa(i)
		Write(255,'(9E15.7)') re_val,(shear1d_sa(jins(j2),4,i)/utau_sa(i)**2, j2=1,7),time_sa(i)
	    End If
	    If (iangle == 1) then
		Write(404,'(9E15.7)') re_val,(angle1d_sa(jins(j2),4,i), j2=1,7),time_sa(i)
	    End If
	    If (ivort == 1) then
	    	Write(206,'(9E15.7)') re_val,(omega1d_sa(jins(j2),1,i), j2=1,7),time_sa(i)
	    	Write(256,'(9E15.7)') re_val,(omega1d_sa(jins(j2),1,i)*nu/utau_sa(i)**2, j2=1,7),time_sa(i)
	    	Write(207,'(9E15.7)') re_val,(omega1d_sa(jins(j2),2,i), j2=1,7),time_sa(i)
	    	Write(257,'(9E15.7)') re_val,(omega1d_sa(jins(j2),2,i)*nu/utau_sa(i)**2, j2=1,7),time_sa(i)
	    	Write(208,'(9E15.7)') re_val,(omega1d_sa(jins(j2),3,i), j2=1,7),time_sa(i)
	    	Write(258,'(9E15.7)') re_val,(omega1d_sa(jins(j2),3,i)*nu/utau_sa(i)**2, j2=1,7),time_sa(i)
	    End If
	    If (ianis == 1) then
	    	Write(209,'(9E15.7)') re_val,(anisb1d_sa(jins(j2),1,i), j2=1,7),time_sa(i)
	    	Write(210,'(9E15.7)') re_val,(anisb1d_sa(jins(j2),2,i), j2=1,7),time_sa(i)
	    	Write(211,'(9E15.7)') re_val,(anisb1d_sa(jins(j2),3,i), j2=1,7),time_sa(i)
	    	Write(212,'(9E15.7)') re_val,(anisb1d_sa(jins(j2),4,i), j2=1,7),time_sa(i)
	    	Write(213,'(9E15.7)') re_val,(anisII1d_sa(jins(j2),i), j2=1,7),time_sa(i)
	    	Write(214,'(9E15.7)') re_val,(anisIII1d_sa(jins(j2),i), j2=1,7),time_sa(i)
	    	Write(215,'(9E15.7)') re_val,(anisf1d_sa(jins(j2),i), j2=1,7),time_sa(i)
	    	Write(216,'(9E15.7)') re_val,(anisg1d_sa(jins(j2),i), j2=1,7),time_sa(i)
	    End If
	    If (iquad == 1) then
		If (if_les == 1) Write(217,'(9E15.7)') re_val,(quad1d_sa(jins(j2),9,i), j2=1,7),time_sa(i)
		If (if_les == 1) Write(267,'(9E15.7)') re_val,(quad1d_sa(jins(j2),9,i)/nu, j2=1,7),time_sa(i)
	    End If
	    If (itke == 1) then
		DO n = 1,4
		    DO m=1,6
			Write(217+(n-1)*6+m,'(9E15.7)') re_val,(kinen1d_sa(jins(j2),n,m,i), j2=1,7),time_sa(i)
			Write(267+(n-1)*6+m,'(9E15.7)') re_val,(kinen1d_sa(jins(j2),n,m,i)*nu/utau_sa(i)**4, j2=1,7),time_sa(i)
		    END DO
		END DO
	    End If
	End Do
    End If
!
    If ((uns1d == 2).or.(uns1d == 3)) then
	ALLOCATE(twrite(ntsamps))
	Do p = 1,ntsamps
	      twrite(p) = ((p-1)*(tpts-1)/(ntsamps-1)) + 1
	End Do
!
!	if_spec_loc_time = 1
      If (spec_loc == 1) Then
	Print*, 'Snapshots in time at specific locations'
	!nsnap_target = 21
	!If (nsnap_target.gt. ntsamps) Then
	!Print*,'Reduce the number of target snapshots in subroutine write_1dpts_unst',ntsamps,nsnap_target
	!End if
	Allocate (snap_re_target(ntsamps))
	Do i=1,ntsamps-2
	snap_re_target(i+1)=4.0e3+(i-1)*1.0e3
	End Do
!
	Do i=1,tpts
	Do j= 2,ntsamps-1
	If (re_sa(i).le. snap_re_target(j).and.re_sa(i+1).gt.snap_re_target(j)) Then
	twrite(j)=i
	Exit
	End if
	End Do
	End Do
	twrite(ntsamps)=tpts
	Deallocate (snap_re_target)
      End If
	WRITE(*,*) 'Location of time snap shots in y -direction', twrite
!
	Do j = jmin,jmax
	    Write(301,'(30E15.7)') yc(j),(stat1d_sa(j,1,1,twrite(i2)), i2=1,ntsamps)
	    Write(451,'(30E15.7)') yc(j),(stat1d_sa(j,2,1,twrite(i2)), i2=1,ntsamps)
	    Write(450,'(30E15.7)') yc(j),(stat1d_sa(j,3,1,twrite(i2)), i2=1,ntsamps)
	    Write(351,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),stat1d_sa(j,1,1,twrite(i2))/utau_sa(twrite(i2)), i2=1,ntsamps)
	    If (iflucs == 1) then
	    	Write(302,'(30E15.7)') yc(j),(sqrt(abs(flucs1d_sa(j,1,2,twrite(i2)))), i2=1,ntsamps)
	    	Write(352,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),sqrt(abs(flucs1d_sa(j,1,2,twrite(i2))))/utau_sa(twrite(i2)), i2=1,ntsamps)
	    	Write(303,'(30E15.7)') yc(j),(sqrt(abs(flucs1d_sa(j,2,2,twrite(i2)))), i2=1,ntsamps)
	    	Write(353,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),sqrt(abs(flucs1d_sa(j,2,2,twrite(i2))))/utau_sa(twrite(i2)), i2=1,ntsamps)
	    	Write(304,'(30E15.7)') yc(j),(sqrt(abs(flucs1d_sa(j,3,2,twrite(i2)))), i2=1,ntsamps)
	    	Write(354,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),sqrt(abs(flucs1d_sa(j,3,2,twrite(i2))))/utau_sa(twrite(i2)), i2=1,ntsamps)

	        Write(392,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),flucs1d_sa(j,1,3,twrite(i2))/abs(flucs1d_sa(j,1,2,twrite(i2)))**1.5, i2=1,ntsamps)
	        Write(393,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),flucs1d_sa(j,2,3,twrite(i2))/abs(flucs1d_sa(j,2,2,twrite(i2)))**1.5, i2=1,ntsamps)
	        Write(394,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),flucs1d_sa(j,3,3,twrite(i2))/abs(flucs1d_sa(j,3,2,twrite(i2)))**1.5, i2=1,ntsamps)

	        Write(395,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),(flucs1d_sa(j,1,4,twrite(i2))/(flucs1d_sa(j,1,2,twrite(i2))**2.0)-3.0), i2=1,ntsamps)
	        Write(396,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),(flucs1d_sa(j,2,4,twrite(i2))/(flucs1d_sa(j,2,2,twrite(i2))**2.0)-3.0), i2=1,ntsamps)	
	        Write(397,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),(flucs1d_sa(j,3,4,twrite(i2))/(flucs1d_sa(j,3,2,twrite(i2))**2.0)-3.0), i2=1,ntsamps)	
	    End If
	    If (ishear == 1) then
		Write(305,'(30E15.7)') yc(j),(shear1d_sa(j,1,twrite(i2)), i2=1,ntsamps)
		Write(355,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),shear1d_sa(j,1,twrite(i2))/utau_sa(twrite(i2))**2, i2=1,ntsamps)
	    End If
	    If (iangle == 1) then
		Write(405,'(30E15.7)') yc(j),(angle1d_sa(j,1,twrite(i2)), i2=1,ntsamps)
	    End If
	    If (ivort == 1) then
	    	Write(398,'(30E15.7)') yc(j),(omega1d_sa(j,4,twrite(i2)), i2=1,ntsamps)
	    	Write(399,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j), &
					omega1d_sa(j,4,twrite(i2))*nu/utau_sa(twrite(i2))**2, i2=1,ntsamps)
	    	Write(306,'(30E15.7)') yc(j),(omega1d_sa(j,1,twrite(i2)), i2=1,ntsamps)
	    	Write(356,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j), &
					omega1d_sa(j,1,twrite(i2))*nu/utau_sa(twrite(i2))**2, i2=1,ntsamps)
	    	Write(307,'(30E15.7)') yc(j),(omega1d_sa(j,2,twrite(i2)), i2=1,ntsamps)
		Write(357,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j), &
					omega1d_sa(j,2,twrite(i2))*nu/utau_sa(twrite(i2))**2, i2=1,ntsamps)
	    	Write(308,'(30E15.7)') yc(j), (omega1d_sa(j,3,twrite(i2)), i2=1,ntsamps)
	    	Write(358,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),omega1d_sa(j,3,twrite(i2))*nu/utau_sa(twrite(i2))**2, i2=1,ntsamps)
	    End If
	    If (ianis == 1) then
	    	Write(359,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),anisb1d_sa(j,1,twrite(i2)), i2=1,ntsamps)
	    	Write(360,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),anisb1d_sa(j,2,twrite(i2)), i2=1,ntsamps)
	    	Write(361,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),anisb1d_sa(j,3,twrite(i2)), i2=1,ntsamps)
	    	Write(362,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),anisb1d_sa(j,4,twrite(i2)), i2=1,ntsamps)
	    	Write(363,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),anisII1d_sa(j,twrite(i2)), i2=1,ntsamps)
	    	Write(364,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),anisIII1d_sa(j,twrite(i2)), i2=1,ntsamps)
	    	Write(365,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),anisf1d_sa(j,twrite(i2)), i2=1,ntsamps)
	    	Write(366,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),anisg1d_sa(j,twrite(i2)), i2=1,ntsamps)
	    End If
	    If (iquad == 1) then
		If (if_les == 1) Write(317,'(30E15.7)') yc(j),(quad1d_sa(j,9,twrite(i2)), i2=1,ntsamps)
		If (if_les == 1) Write(367,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),quad1d_sa(j,9,twrite(i2))/nu, i2=1,ntsamps)
	    End If
	    If (itke == 1) then
		DO n = 1,4
		    DO m=1,6
			Write(317+(n-1)*6+m,'(30E15.7)') yc(j),(kinen1d_sa(j,n,m,twrite(i2)), i2=1,ntsamps)
			Write(367+(n-1)*6+m,'(60E15.7)') yc(j),(ycplus_sa(twrite(i2),j),kinen1d_sa(j,n,m,twrite(i2))*nu/utau_sa(twrite(i2))**4, i2=1,ntsamps)
		    END DO
		END DO
	    End If
	End Do
    End If
!
	END SUBROUTINE write_1dpts_unst
!*******************************************************************************
	SUBROUTINE write_rates_unst
	INTEGER :: i,j,m
	REAL(mytype) :: re_val,y_val
!
	Do j = jmin,jmax
	Do i = 1,tpts-1
	    If (re_print==0) Then; re_val=re_sa(i)
	    Else If (re_print==1) Then; re_val=retau_sa(i)
	    Else; re_val=time_sa(i)
	    End If
	    If (y_print==0) Then; y_val=yc(j)
	    Else; y_val=ycplus_sa(i,j)
	    End If
  	    Write(102,'(7E15.7)') re_val,y_val,(stat1d_sa(j,1,1,i+1)-stat1d_sa(j,1,1,i))/stat1d_sa(j,1,1,i), &
				(sqrt(abs(flucs1d_sa(j,1,2,i+1)))-sqrt(abs(flucs1d_sa(j,1,2,i))))/sqrt(abs(flucs1d_sa(j,1,2,i))), &
				(sqrt(abs(flucs1d_sa(j,2,2,i+1)))-sqrt(abs(flucs1d_sa(j,2,2,i))))/sqrt(abs(flucs1d_sa(j,2,2,i))), &
				(sqrt(abs(flucs1d_sa(j,3,2,i+1)))-sqrt(abs(flucs1d_sa(j,3,2,i))))/sqrt(abs(flucs1d_sa(j,3,2,i))), &
				(sqrt(abs(flucs1d_sa(j,4,2,i+1)))-sqrt(abs(flucs1d_sa(j,4,2,i))))/sqrt(abs(flucs1d_sa(j,4,2,i)))
	End Do;End Do
!
	END SUBROUTINE write_rates_unst
!*******************************************************************************
SUBROUTINE write_1d_unst
  INTEGER :: i,j,m,fileval,j2
  REAL(mytype) :: re_val,y_val
  INTEGER, DIMENSION (7) :: jins
!
  Do i = 1,tpts
    If (re_print==0) Then; re_val=re_sa(i)
    Else If (re_print==1) Then; re_val=retau_sa(i)
    Else; re_val=time_sa(i)
    End If
    !Write(76,'(2E15.7)') time_sa(i),re_val  
  End Do
!        
  Do j = jmin,jmax
    Do i = 1,tpts
      If (re_print==0) Then; re_val=re_sa(i)
      Else If (re_print==1) Then; re_val=retau_sa(i)
      Else; re_val=time_sa(i)
      End If
      If (y_print==0) Then; y_val=yc(j)
      Else; y_val=ycplus_sa(i,j)
      End If

      Write(70,'(6E15.7)') time_sa(i),y_val,stat1d_sa(j,1,1,i), &
                           stat1d_sa(j,2,1,i),stat1d_sa(j,3,1,i),stat1d_sa(j,4,1,i)
      Write(71,'(6E15.7)') re_val,y_val,stat1d_sa(j,1,1,i), &
                           stat1d_sa(j,2,1,i),stat1d_sa(j,3,1,i),stat1d_sa(j,4,1,i)
!
      Write(90,'(6E15.7)') re_val,y_val,stat1d_sa(j,1,1,i)/utau_sa(i),stat1d_sa(j,2,1,i)/utau_sa(i), &
                           stat1d_sa(j,3,1,i)/utau_sa(i),stat1d_sa(j,4,1,i)/utau_sa(i)**2/rho
!
      If (iflucs == 1) then
        Write(78,'(6E15.7)') re_val,y_val,sqrt(abs(flucs1d_sa(j,1,2,i))),sqrt(abs(flucs1d_sa(j,2,2,i))), &
                             sqrt(abs(flucs1d_sa(j,3,2,i))),sqrt(abs(flucs1d_sa(j,4,2,i)))
        Write(82,'(6E15.7)') re_val,y_val,sqrt(abs(flucs1d_sa(j,1,2,i)))/utau_sa(i),sqrt(abs(flucs1d_sa(j,2,2,i)))/utau_sa(i), &
                             sqrt(abs(flucs1d_sa(j,3,2,i)))/utau_sa(i),sqrt(abs(flucs1d_sa(j,4,2,i)))/utau_sa(i)
        Write(72,'(6E15.7)') re_val,y_val,(flucs1d_sa(j,m,3,i)/abs(flucs1d_sa(j,m,2,i))**1.5, m=1,4)
        Write(73,'(6E15.7)') re_val,y_val,((flucs1d_sa(j,m,4,i)/(flucs1d_sa(j,m,2,i)**2.0)-3.0), m=1,4)
     End If
!
     If (ishear == 1) then
       Write(83,'(6E15.7)') re_val,y_val,shear1d_sa(j,1,i),shear1d_sa(j,2,i), &
                            shear1d_sa(j,3,i),shear1d_sa(j,4,i)
       Write(84,'(6E15.7)') re_val,y_val,shear1d_sa(j,1,i)/utau_sa(i)**2,shear1d_sa(j,2,i)/utau_sa(i)**2, &
                            shear1d_sa(j,3,i)/utau_sa(i)**2,shear1d_sa(j,4,i)/utau_sa(i)**2
     End If
     If (iangle == 1) then
       Write(403,'(7E15.7)') re_val,y_val,(angle1d_sa(j,m,i),m=1,5)
     End If
!
     If (itke == 1) then
       Write(86,'(8E15.7)') re_val,y_val,(kinen1d_sa(j,1,m,i),m=1,6)
       Write(186,'(8E15.7)') re_val,y_val,(kinen1d_sa(j,1,m,i)*nu/utau_sa(i)**4,m=1,6)
       Write(87,'(8E15.7)') re_val,y_val,(kinen1d_sa(j,2,m,i),m=1,6)
       Write(187,'(8E15.7)') re_val,y_val,(kinen1d_sa(j,2,m,i)*nu/utau_sa(i)**4,m=1,6)
       Write(88,'(8E15.7)') re_val,y_val,(kinen1d_sa(j,3,m,i),m=1,6)
       Write(188,'(8E15.7)') re_val,y_val,(kinen1d_sa(j,3,m,i)*nu/utau_sa(i)**4,m=1,6)
       Write(89,'(8E15.7)') re_val,y_val,(kinen1d_sa(j,4,m,i),m=1,6)
       Write(189,'(8E15.7)') re_val,y_val,(kinen1d_sa(j,4,m,i)*nu/utau_sa(i)**4,m=1,6)
       Write(92,'(8E15.7)') re_val,y_val,((kinen1d_sa(j,1,m,i)+kinen1d_sa(j,2,m,i)+kinen1d_sa(j,3,m,i))/2.0,m=1,6)
       Write(192,'(8E15.7)') re_val,y_val,((kinen1d_sa(j,1,m,i)+kinen1d_sa(j,2,m,i)+kinen1d_sa(j,3,m,i))*nu/2.0/utau_sa(i)**4,m=1,6)
    End If
!
    If (ianis == 1) Write(91,'(10E15.7)') re_val,y_val,(anisb1d_sa(j,m,i),m=1,4),anisII1d_sa(j,i), &
                                          anisIII1d_sa(j,i),anisf1d_sa(j,i),anisg1d_sa(j,i)
!
    If (if_les == 1) Then
      Write(101,'(6E15.7)') re_val,y_val,quad1d_sa(j,7,i),quad1d_sa(j,7,i)/nu, quad1d_sa(j,7,i)/nu,quad1d_sa(j,9,i)
    End If
!
    If (ivort == 1) then
!      Write(148,'(5E15.7)') re_val,y_val,omegamn1d_sa(j,1,i),omegamn1d_sa(j,2,i),omegamn1d_sa(j,3,i)
!      Write(149,'(5E15.7)') re_val,y_val,omegamn1d_sa(j,1,i)*nu/utau_sa(i)**2,omegamn1d_sa(j,2,i)*nu/utau_sa(i)**2,&
!                            omegamn1d_sa(j,3,i)*nu/utau_sa(i)**2
      Write(148,'(5E15.7)') re_val,y_val,omega1d_sa(j,4,i),omega1d_sa(j,5,i),omega1d_sa(j,6,i)
      Write(149,'(5E15.7)') re_val,y_val,omega1d_sa(j,4,i)*nu/utau_sa(i)**2,omega1d_sa(j,5,i)*nu/utau_sa(i)**2,&
                            omega1d_sa(j,6,i)*nu/utau_sa(i)**2
      Write(150,'(5E15.7)') re_val,y_val,omega1d_sa(j,1,i),omega1d_sa(j,2,i),omega1d_sa(j,3,i)
      Write(151,'(5E15.7)') re_val,y_val,omega1d_sa(j,1,i)*nu/utau_sa(i)**2,omega1d_sa(j,2,i)*nu/utau_sa(i)**2,&
                            omega1d_sa(j,3,i)*nu/utau_sa(i)**2
    End If
!
    If (iquad == 1) then
      Write(98,'(11E15.7)') re_val,y_val,(quad1d_sa(j,m,i),m=1,7)
      Write(198,'(11E15.7)') re_val,y_val,(quad1d_sa(j,m,i)/utau_sa(i)**2,m=1,7)
    End If
  End Do; End Do
!
END SUBROUTINE write_1d_unst
!*******************************************************************************
!
END MODULE write_post
