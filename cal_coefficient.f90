!	BRIEF DESCRIPTION
!*******************************************************************************
!
!	Coeffieient names are as following
!
!	'gcgx_f' means gradient in 'x' direction when variable is at 'cell center'
!	and gradient is reuired at 'cell face'
!
!
!	'gfgy_f' means gradient in 'y' direction when variable is at 'cell face'
!	and gradient is reuired at 'cell face'
!
!
!	'gcgz_c' means gradient in 'z' direction when variable is at 'cell center'
!	and gradient is reuired at 'cell center'
!
!
!	'gfgx_f' means gradient in 'x' direction when variable is at 'cell face'
!	and gradient is reuired at 'cell face'
!	
!       Added gcdz_f to consider nonuniform grid in z-Direction
!
! 	Defined gfdx_c,gfdy_c,gfdz_c in order to have a symmetry between 
!       coeffiecints otherwise for second order accuracy there is no need 
!	to define this
!
!	***SCHEME1*** Determines gradient at cell center when the variables are
!	at cell face
!	***SCHEME2*** Determines gradient at cell face when the variables are 
!	at cell cell center
!       ***SCHEME3*** Determines gradient at cell center when the variables are
!	at cell center
!	***SCHEME4*** Determines gradient at cell face when the variable are at
!	cell face (which is not required in the present code)
!
	Module cal_coefficient
!	USE shared_data
	use shared_post
!
	Implicit None
	Integer , Private :: isec
	Contains
!
!******************************************************************************
	Subroutine COMPUTE_FACT
!******************************************************************************
!------ This procedure computes geometrical factors to calculate gradients
!------ using 3rd order interpolation.
	Implicit None
!
	Call GFDS_C_COMP
	Call GCDS_F_COMP
	Call GCDS_C_COMP
	!Call COFF_CALC_OUT
!
	Return
        End Subroutine COMPUTE_FACT
!
!******************************************************************************
	Subroutine GFDS_C_COMP
!******************************************************************************
!	This procedure calculates gradient at cel center with variable laocated
!	at cell face using two points
	Implicit None
	Integer :: ib,iblock
!
        iblock = 0
        ib=iblock+1
!
!	GRADIENT AT CELL CENTER WHEN VARIABLE IS AT CELL FACE (TWO-POINT)
!
	Call SCHEME1(gfdx_c,dx,0,nx,1,ib)
	Call SCHEME1(gfdy_c,dy,1,ny,1,ib)
	Call SCHEME1(gfdz_c,dz,1,nz,1,ib)
	Call SCHEME1(gfdz_c2,dz,1,nz,0,ib)
If (if2d_noy == 1) gfdy_c = 0.0e0
!
	Return
        End Subroutine GFDS_C_COMP 
!
!******************************************************************************
	Subroutine GCDS_F_COMP
!*****************************************************************************
!	THis procedure calculates the gradient at cell face with variable located
!	at cell center using four points
	Implicit None
	Integer :: i,ib,iblock,j,m
!
        iblock = 0
        ib=iblock+1
!
!	GRADIENTS AT CELL FACE WHEN VARIABLE IS AT CELL CENTER(FOUR_POINT)
!
	If (igr_2pt.eq.0) Then
	Call SCHEME2(dx,1,nx,1,1,ib)
	Call SCHEME2(dy,1,ny,2,1,ib)
	!Call SCHEME2(dz,1,nz,3,1,ib)
	Else If (igr_2pt.eq.1) Then
	Do i = 1,nxt
	gcdx_f(i,1,1,ib) = 0.0e0
	gcdx_f(i,2,1,ib) = 1.0/dxs(i)
	gcdx_f(i,3,1,ib) = -1.0/dxs(i)
	gcdx_f(i,4,1,ib) = 0.0e0
	End Do
	Do i = 1,nyt
        gcdy_f(i,1,1,ib) = 0.0e0
        gcdy_f(i,2,1,ib) = 1.0/dys(i)
        gcdy_f(i,3,1,ib) = -1.0/dys(i)
        gcdy_f(i,4,1,ib) = 0.0e0
        End Do
	End If
!
	Do i = 1,nz
	gcdz_f(i,1,1,ib) = 0.0e0
	gcdz_f(i,2,1,ib) = 1.0/dzs(i)
	gcdz_f(i,3,1,ib) = -1.0/dzs(i)
	gcdz_f(i,4,1,ib) = 0.0e0
	gcdz_f2(i,1,1,ib) = 0.0e0
	gcdz_f2(i,2,1,ib) = 1.0/dzs(i+1)
	gcdz_f2(i,3,1,ib) = -1.0/dzs(i+1)
	gcdz_f2(i,4,1,ib) = 0.0e0
	End Do
!
!	Pressure coeffiients
	Do i = 1,nxt
        do m = 1,4
          gcdx_f(i,m,2,ib) = gcdx_f(i,m,1,ib)
        End Do
	End Do
	Do i = 1,nyt
        do m = 1,4
         gcdy_f(i,m,2,ib) = gcdy_f(i,m,1,ib)
        End Do
	End Do
	Do i = 1,nzt
        do m = 1,4
         gcdz_f(i,m,2,ib) = gcdz_f(i,m,1,ib)
        End Do
        End Do
!
!	Boundary gradeints
!
	If (igr_2pt.eq.0) Then
	If (iconbc == 1)Call SCHEME2_BC(dx,1,nx,ib,1)
	If (jconbc /= 0)Call SCHEME2_BC(dy,1,ny,ib,2)
	!If (kconbc.eq.1)Call SCHEME2_BC(dz,1,nz,ib,3)
	Else If (igr_2pt.eq.1) Then
	If (iconbc == 1) then
	i=1
        gcdx_f(i,1,1,ib) = 0.0e0
        gcdx_f(i,2,1,ib) = 2.0/dx(i)
        gcdx_f(i,3,1,ib) = -2.0/dx(i)
        gcdx_f(i,4,1,ib) = 0.0e0
	gcdx_f(i,1:4,2,ib) = 0.0e0
	i=nxt
	gcdx_f(i,1,1,ib) = 0.0e0
        gcdx_f(i,2,1,ib) = 2.0/dx(nx)
        gcdx_f(i,3,1,ib) = -2.0/dx(nx)
        gcdx_f(i,4,1,ib) = 0.0e0
	gcdx_f(i,1:4,2,ib) = 0.0e0
	End If
	If (jconbc /= 0) Then
	j=1
        gcdy_f(j,1,1,ib) = 0.0e0
        gcdy_f(j,2,1,ib) = 2.0/dy(j)
        gcdy_f(j,3,1,ib) = -2.0/dy(j)
        gcdy_f(j,4,1,ib) = 0.0e0
	gcdy_f(j,1:4,2,ib) = 0.0e0
        j=nyt
        gcdy_f(j,1,1,ib) = 0.0e0
        gcdy_f(j,2,1,ib) = 2.0/dy(ny)
        gcdy_f(j,3,1,ib) = -2.0/dy(ny)
        gcdy_f(j,4,1,ib) = 0.0e0
	gcdy_f(j,1:4,2,ib) = 0.0e0
	End If
        End If
!
!	Gradients at the block sur=face
!
	If (ibmpfg /= 0) Then
	gcdx_f(1:nxt,1:4,1:2,2) = gcdx_f(1:nxt,1:4,1:2,1)
	gcdy_f(1:nyt,1:4,1:2,2) = gcdy_f(1:nyt,1:4,1:2,1)
	gcdz_f(1:nzt,1:4,1:2,2) = gcdz_f(1:nzt,1:4,1:2,1)
        Call SCHEME2_BC_BLOCK
	End If
If (if2d_noy == 1) gcdy_f = 0.0e0
	Return
        End Subroutine GCDS_F_COMP 
!
!******************************************************************************
	Subroutine GCDS_C_COMP
!*****************************************************************************
!	This subroutine calculates the gradient at cell center with variable located
!	cell center using five point formulation.
	Implicit None
	Integer :: i,ib,iblock
!
        iblock = 0
        ib=iblock+1
!
!	GRADIENTS AT CELL CENTER WHEN VARIABLE IS AT CELL CENTER(FIVE_POINT)
!
	Call SCHEME3(gcdx_c,dx,1,nx,1,ib)
	Call SCHEME3(gcdy_c,dy,1,ny,2,ib)
	Call SCHEME3(gcdz_c,dz,1,nz,3,ib)
	Call SCHEME3_LES(grdsxx,dx,1,nx/2,1,ib)
	Call SCHEME3_LES(grdsyy,dy,1,ny/2,1,ib)
	Call SCHEME3_LES(grdszz,dz,1,nz/2,2,ib)
	If (iconbc == 1)Call SCHEME3_BC(gcdx_c,grdsxx,dx,1,nx,1,ib) 
	If (jconbc /= 0)Call SCHEME3_BC(gcdy_c,grdsyy,dy,1,ny,1,ib) 
	!If (kconbc.ne.0) Call SCHEME3_BC(gcdz_c,grdszz,dz,1,nz,2,ib) 
!
!	Gradients near block surface
!
	If (ibmpfg /= 0) Call SCHEME3_BC_BLOCK
!
If (if2d_noy == 1) gcdy_c = 0.0e0

	Return
        End Subroutine GCDS_C_COMP 
!
!******************************************************************************
	Subroutine COFF_CALC_OUT
!*****************************************************************************
!	This procedure writes the gradient calculated in this module to cmpfct.dat file
	Implicit None
	Integer :: i,j,k,m
!
        Write(11,*)'i,gcdx_f(i,1,1,1)'
        Do i = 1,nxt
          Write(11,*) i,(gcdx_f(i,m,1,1),m=1,4)
	End Do
        Write(11,*)'i,gcdx_f(i,1,2,1)'
        Do i = 1,nxt
          Write(11,*) i,(gcdx_f(i,m,2,1),m=1,4)
	End Do
	If (ibmpfg.eq.1) Then
        Write(11,*)'i,gcdx_f(i,1,1,2)'
        Do i = 1,nxt
          Write(11,*) i,(gcdx_f(i,m,1,2),m=1,4)
	End Do
        Write(11,*)'i,gcdx_f(i,1,2,2)'
        Do i = 1,nxt
          Write(11,*) i,(gcdx_f(i,m,2,2),m=1,4)
	End Do
	End If
        Write(11,*)'i,gcdx_c(i,1,1,1)'
        Do i = 1,nx
          Write(11,*) i,(gcdx_c(i,m,1,1),m=1,5)
	End Do
        Write(11,*)'i,gcdx_c(i,1,2,1)'
        Do i = 1,nx
          Write(11,*) i,(gcdx_c(i,m,2,1),m=1,5)
	End Do
        Write(11,*)'j,Y(j),dy(j)'
        DO  j=1,nyt
          Write(11,*) j,y(j),dy(j)
	End Do
        Write (11,*)'j,gcdy_f(j,1,1,1)'
        Do j = 1,nyt
          Write(11,*) j,(gcdy_f(j,m,1,1),m=1,4)
	End Do
        WRITe(11,*)'j,gcdy_f(j,1,2,1)'
        Do j = 1,nyt
          Write(11,*) j,(gcdy_f(j,m,2,1),m=1,4)
	End Do
	If (ibmpfg.eq.1) Then
        Write (11,*)'j,gcdy_f(j,1,1,2)'
        Do j = 1,nyt
          Write(11,*) j,(gcdy_f(j,m,1,2),m=1,4)
	End Do
        WRITe(11,*)'j,gcdy_f(j,1,2,2)'
        Do j = 1,nyt
          Write(11,*) j,(gcdy_f(j,m,2,2),m=1,4)
	End Do
	End If
        Write(11,*)'j,gcdy_c(j,1,1,1)'
        Do j = 1,ny
          Write(11,*) j,(gcdy_c(j,m,1,1),m=1,5)
	End Do
        Write(11,*)'j,gcdy_c(j,1,2,1)'
        Do j = 1,ny
          Write(11,*) j,(gcdy_c(j,m,2,1),m=1,5)
	End Do
!
        !Write(11,*)'i,gcdz_f(i,1,1,1)'
        !Do i = 1,nxt
        !  Write(11,*) i,(gcdz_f(i,m,1,1),m=1,4)
	!End Do
        !Write(11,*)'i,gcdz_f(i,1,2,1)'
        !Do i = 1,nxt
        !  Write(11,*) i,(gcdz_f(i,m,2,1),m=1,4)
	!End Do
	!Close(11)
!
        Return
        End Subroutine COFF_CALC_OUT
!
!******************************************************************************
	Subroutine THIRD (ij,a,b,c,d,iuv,ik,ib)
!*****************************************************************************
!------ This subroutine computes gradient coefficients using 3rd order
	Implicit None
	Integer,Intent( In ) :: ij,iuv,ik,ib
  real(mytype),Intent( In ) :: a,b,c,d
  real(mytype) :: temp0,temp1,temp2,temp3,temp4,temp5,deno
!
        temp0=c*c*d*d
        temp3=temp0*(c-d)
        temp4=b*b*c*c*(b-c)
        temp5=a*a*b*b*(a-b)
        temp1=(d+c)*temp4+(c*c-b*b)*temp0
        temp2=(a*a-b*b)*temp0-(d+c)*temp5
        deno=temp1*(b*a*(b-a)*temp3-c*d*(d-c)*temp5)-temp2*(c*d*(d-c)*temp4-   &
	c*b*(c-b)*temp3)
!
        If (iuv == 1) Then
          gcdx_f(ij,1,ik,ib) = temp1*temp3*b*b/deno
          gcdx_f(ij,2,ik,ib) = (temp2*temp3*c*c-temp1*temp3*a*a)/deno
          gcdx_f(ij,3,ik,ib) = (-temp1*temp5*d*d-temp2*temp4*d*d-temp2*temp3*    &
	  b*b)/deno
          gcdx_f(ij,4,ik,ib) = (temp1*temp5*c*c+temp2*temp4*c*c)/deno
        Else If (iuv == 2) Then
          gcdy_f(ij,1,ik,ib) = temp1*temp3*b*b/deno
          gcdy_f(ij,2,ik,ib) = (temp2*temp3*c*c-temp1*temp3*a*a)/deno
          gcdy_f(ij,3,ik,ib) = (-temp1*temp5*d*d-temp2*temp4*d*d-temp2*temp3*     &
	  b*b)/deno
          gcdy_f(ij,4,ik,ib) = (temp1*temp5*c*c+temp2*temp4*c*c)/deno
	Else If (iuv == 3) Then
          gcdz_f(ij,1,ik,ib) = temp1*temp3*b*b/deno
          gcdz_f(ij,2,ik,ib) = (temp2*temp3*c*c-temp1*temp3*a*a)/deno
          gcdz_f(ij,3,ik,ib) = (-temp1*temp5*d*d-temp2*temp4*d*d-temp2*temp3*     &
          b*b)/deno
          gcdz_f(ij,4,ik,ib) = (temp1*temp5*c*c+temp2*temp4*c*c)/deno
        End If
!
        Return
        End Subroutine THIRD 
!
!******************************************************************************
	Subroutine THIRD3 (ij,a,b,c,d,iuv,nabv,ik,ib)
!*****************************************************************************
!------ This subroutine computes gradient coefficients using 3rd order
!       NABV is the number of points above the face under consideration.
!
	Implicit None
!
	Integer, Intent( In ) :: ij,nabv,iuv,ik,ib
  real(mytype), Intent( In ) :: a,b,c,d
  real(mytype) :: deno
!
        If (NABV == 1) Then
          deno=a*(b+c)*(b-c)+b*(c+a)*(c-a)+c*(a+b)*(a-b)
          If (iuv == 1) Then
            gcdx_f(ij,1,ik,ib) = (b+c)*(b-c)/deno
            gcdx_f(ij,2,ik,ib) = (c+a)*(c-a)/deno
            gcdx_f(ij,3,ik,ib) = (a+b)*(a-b)/deno
            gcdx_f(ij,4,ik,ib) = 0.0e0
          Else If (iuv == 2) Then
            gcdy_f(ij,1,ik,ib) = (b+c)*(b-c)/deno
            gcdy_f(ij,2,ik,ib) = (c+a)*(c-a)/deno
            gcdy_f(ij,3,ik,ib) = (a+b)*(a-b)/deno
            gcdy_f(ij,4,ik,ib) = 0.0e0
          Else If (iUV == 3) Then
            gcdz_f(ij,1,ik,ib) = (b+c)*(b-c)/deno
            gcdz_f(ij,2,ik,ib) = (c+a)*(c-a)/deno
            gcdz_f(ij,3,ik,ib) = (a+b)*(a-b)/deno
            gcdz_f(ij,4,ik,ib) = 0.0e0
          End If
        Else If (nabv == 2) Then
          deno=b*(c+d)*(c-d)+c*(d+b)*(d-b)+d*(b+c)*(b-c)
          If (iuv == 1) Then
            gcdx_f(ij,1,ik,ib) = 0.0e0
            gcdx_f(ij,2,ik,ib) = (c+d)*(c-d)/deno
            gcdx_f(ij,3,ik,ib) = (d+b)*(d-b)/deno
            gcdx_f(ij,4,ik,ib) = (b+c)*(b-c)/deno
          Else If (iuv == 2) Then
            gcdy_f(ij,1,ik,ib) = 0.0e0
            gcdy_f(ij,2,ik,ib) = (c+d)*(c-d)/deno
            gcdy_f(ij,3,ik,ib) = (d+b)*(d-b)/deno
            gcdy_f(ij,4,ik,ib) = (b+c)*(b-c)/deno
          Else If (iuv == 3) Then
            gcdz_f(ij,1,ik,ib) = 0.0e0
            gcdz_f(ij,2,ik,ib) = (c+d)*(c-d)/deno
            gcdz_f(ij,3,ik,ib) = (d+b)*(d-b)/deno
            gcdz_f(ij,4,ik,ib) = (b+c)*(b-c)/deno
          End If
        End If
!
        Return
        End Subroutine THIRD3
!
!******************************************************************************
Subroutine COMP_PRES_COFF
!*****************************************************************************
!------ This procedure computes coefficients a's before solving Poisson eq.'
!       for pressure.
Implicit None
Integer :: i,j,k
!
   Do k = 1,nzt
   Do j = 1,nyt
   Do i = 1,nxt
      ae(i,j,k) = dy(j)*dz(k)/dxs(i)
      as(i,j,k) = dx(i)*dz(k)/dys(j)
      at(i,j,k) = dx(i)*dy(j)/dzs(k)
      ap(i,j,k) = 1.0e0
   End Do; End Do; End Do
If (if2d_noy == 1) as = 0.0e0
!
as(:,1,:) = 0.0e0
Do k = 1,nzt
Do i = 1,nxt
   ae(i,nyt,k) = 0.0e0
   as(i,nyt,k) = 0.0e0
   at(i,nyt,k) = 0.0e0
End Do; End Do
!
!       Periodic B.C.
!       AE, AS, AT are O.K at the inlet plane.
If (iconbc == 0) Then
   ae(nxt,:,:) = ae(1,:,:)
!
!       Convective B.C.
Else If (iconbc == 1) Then
   ae(1,:,:) = 0.0e0
   ae(nxt,:,:) = 0.0e0
End If
!
!------ For cutting the links
!
Do k = 1,nz
Do j = 1,ny
Do i = 1,nx
   ap(i,j,k) = ae(i,j,k)+ae(i+1,j,k)+as(i,j,k)+as(i,j+1,k)+at(i,j,k)+at(i,j,k+1)
   ap(i,j,k) = 1.0e0/ap(i,j,k)
End Do; End Do; End Do
!
!       Periodic B.C.
If (iconbc == 0) ap(nxt,:,:) = ap(1,:,:)
!
        Return
        End Subroutine COMP_PRES_COFF
!
!******************************************************************************
	Subroutine SCHEME1(temp,ds,mm,nn,jj,ib)
!******************************************************************************
!	THIS Subroutine CALCULATES THE GRADIENT AT CELL CENTER USING TWO POINTS 
!	WHEN THE VARABLES ARE LOCATED AT CELL FACE.
!       jj IS AN ARGUMENT WHICH IS NEEDED TO CALCULATE GRADIENT AT NEXT CELL OR
!	THE PREVIOUS ONE
	Implicit None
	Integer, Intent( In ) :: nn,mm,jj,ib
  real(mytype), Dimension(mm:nn+1,4,ngr,nbl), Intent( InOut ) :: temp
  real(mytype), Dimension(-1:nn+2), Intent( In ) :: ds
	Integer :: i,j
!
	If (jj == 1) Then
          Do i = mm,nn
            temp(i,1,1,ib) = 0.0e0
            temp(i,2,1,ib) = 1.0e0/ds(i)
            temp(i,3,1,ib) = -1.0e0/ds(i)
            temp(i,4,1,ib) = 0.0e0
          End Do
	Else If (jj == 0) Then
	  Do i = mm,nn
            temp(i,1,1,ib) = 0.0e0
            temp(i,2,1,ib) = 1.0e0/ds(i-1)
            temp(i,3,1,ib) = -1.0e0/ds(i-1)
            temp(i,4,1,ib) = 0.0e0
          End Do
	Else If (jj == 2) Then
	  Do i = mm,nn
            temp(i,1,1,ib) = 0.0e0
            temp(i,2,1,ib) = 1.0e0/ds(i+1)
            temp(i,3,1,ib) = -1.0e0/ds(i+1)
            temp(i,4,1,ib) = 0.0e0
          End Do
	Else
	Print*,'jj MUST BE BETWEEN 0 AND 2 IN SCHEME1 Subroutine'
	STOP
	End If
!
	If (ibmpfg.eq.1) Then
        Do j=1,4
        Do i = mm,nn
        temp(i,j,1,2) = temp(i,j,1,1)
        End Do
        End Do
        End If
!
	Return
	End Subroutine SCHEME1
!
!******************************************************************************
	Subroutine SCHEME2(ds,mm,nn,iuvw,ik,ib)
!******************************************************************************
!	THIS Subroutine CALCULATES THE GRADIENT AT CELL FACE WITH VARIABLES AT THE
!	CELL CENTER(4-POINT)
	Implicit None
	Integer, Intent( In ) :: nn,mm,ib,iuvw,ik
  real(mytype), Dimension(-1:nn+2), Intent( In ) :: ds
	Integer :: i,m
  real(mytype) :: a,b,c,d
!
	Do i = mm,nn+1
          a=ds(i)+0.5e0*ds(i+1)
          b = 0.5e0*ds(i)
          c = -0.5e0*ds(i-1)
          d = -ds(i-1)-0.5e0*ds(i-2)
!
          Call THIRD (i,a,b,c,d,iuvw,ik,ib)
        End Do
!
	Return
	End Subroutine SCHEME2
!
!******************************************************************************
	Subroutine SCHEME2_BC(ds,mm,nn,ib,ij)
!******************************************************************************
! 	THIS Subroutine DETERMINES THE GRADIENT AT CELL FACE AT BOUNDARIES WITH
!	VARibLES AT CELL CENTER(ibwest=1=VARIABLE AT BOUNDARY, ibwest = 0=AT CELL CENTER)
!	Ij=1 X_DIRECTION, Ij=2 Y-dIRECTION,Ij=3 Z DIRECTION
	Implicit None
	Integer, Intent( In ) :: nn,mm,ib,ij
  real(mytype), Dimension(-1:nn+2), Intent( In ) :: ds
	Integer :: i,m
  real(mytype) :: a,b,c,d
!
        Do i = mm,mm+1
          a=ds(i)+0.5e0*ds(i+1)
          b = 0.5e0*ds(i)
          c = -0.5e0*ds(i-1)
          d = -ds(i-1)-0.5e0*ds(i-2)
	  If (ij == 1.AND.i == 1) Then
            If (ibwest == 1) Then
              c = 0.0e0
              d = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,1,1,ib)
            Else If (ibwest == 0) Then
              d = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,1,1,ib)
            End If
            do m = 1,4
              gcdx_f(i,m,2,ib) = 0.0e0
            End Do
	  Else If (ij == 1.AND.i == 2) Then
            If (ibwest == 1) Then
              d = -ds(i-1)
              Call THIRD (i,a,b,c,d,ij,1,ib)
            End If
            Call THIRD3 (i,a,b,c,d,ij,1,2,ib) 
	  Else If (ij == 2.AND.i == 1) Then
            If (ibsoth == 1) Then
              c = 0.0e0
              d = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,1,1,ib)
            Else If (ibsoth == 0) Then
              d = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,1,1,ib)
            End If
            do m = 1,4
              gcdy_f(i,m,2,ib) = 0.0e0
            End Do
          Else If (ij == 2.AND.i == 2) Then
            If (ibsoth == 1) Then
              d = -ds(i-1)
             Call THIRD (i,a,b,c,d,ij,1,ib)
            End If
	    Call THIRD3 (i,a,b,c,d,ij,1,2,ib)
	  Else If (ij == 3.AND.i == 1) Then
            If (kbsoth == 1) Then
              c = 0.0e0
              d = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,1,1,ib)
            Else If (kbsoth == 0) Then
              d = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,1,1,ib)
            End If
            do m = 1,4
              gcdz_f(i,m,2,ib) = 0.0e0
            End Do
          Else If (ij == 3.AND.i == 2) Then
            If (kbsoth == 1) Then
              d = -ds(i-1)
             Call THIRD (i,a,b,c,d,ij,1,ib)
            End If
            Call THIRD3 (i,a,b,c,d,ij,1,2,ib)
	  Else
	    Print*,'IJ VALUE MUST BE BETWEEN 1 AND 3 IN SCHEM2_BC Subroutine'
	    STOP
	  End If
	End Do
!
	Do i = nn,nn+1
          a=ds(i)+0.5e0*ds(i+1)
          b = 0.5e0*ds(i)
          c = -0.5e0*ds(i-1)
          d = -ds(i-1)-0.5e0*ds(i-2)
!
	  If (ij == 1.AND.i == nn) Then 
            If (ibeast == 1) Then
              a=ds(i)
              Call THIRD (i,a,b,c,d,ij,1,ib)
            End If
            Call THIRD3 (i,a,b,c,d,ij,2,2,ib)
	  Else If (ij == 1.AND.i == (nn+1)) Then
	    If (ibeast == 1) Then
              a = 0.0e0
              b = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,2,1,ib)
            Else If (ibeast == 0) Then
              a = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,2,1,ib)
            End If
            do m = 1,4
              gcdx_f(i,m,2,ib) = 0.0e0
            End Do
	  Else If (ij == 2.AND.i == nn) Then 
	   If (ibnoth == 1) Then
             a=ds(i)
             Call THIRD (i,a,b,c,d,ij,1,ib)
           End If
           Call THIRD3 (i,a,b,c,d,ij,2,2,ib)
	  Else If (ij == 2 .And. i == nn+1) Then
            If (ibnoth == 1) Then
              a = 0.0e0
              b = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,2,1,ib)
            Else If (ibnoth == 0) Then
              a = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,2,1,ib)
            End If
            do m = 1,4
              gcdy_f(i,m,2,ib) = 0.0e0
            End Do
	  Else If (ij == 3.AND.i == nn) Then
           If (kbnoth == 1) Then
             a=ds(i)
             Call THIRD (i,a,b,c,d,ij,1,ib)
           End If
           Call THIRD3 (i,a,b,c,d,ij,2,2,ib)
          Else If (ij == 3.AND.i == (nn+1)) Then
            If (kbnoth == 1) Then
              a = 0.0e0
              b = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,2,1,ib)
            Else If (kbnoth == 0) Then
              a = 0.0e0
              Call THIRD3 (i,a,b,c,d,ij,2,1,ib)
            End If
            do m = 1,4
              gcdz_f(i,m,2,ib) = 0.0e0
            End Do
	  Else
	    Print*,'IJ VALUE MUST BE BETWEEN 1 AND 3 IN SCHEME2_BC'
	  End If
	End Do
!
	Return
	End Subroutine SCHEME2_BC
!
!******************************************************************************
        Subroutine SCHEME3(temp,ds,mm,nn,jj,ib)
!******************************************************************************
!       THIS Subroutine CALCULATES THE GRADIENT AT CELL CENTER USING THREE POINTS
!       WHEN THE VARABLES ARE LOCATED AT CELL CENTER.
!       Jj=1 X DIRECTION ,Jj=2 Y DIRECTION, Jj=3 Z DIRECTION 
        Implicit None
        Integer, Intent( In ) :: nn,mm,ib,jj
  real(mytype), Dimension(mm-1:nn+1,5,ngr,nbl), Intent( InOut ) :: temp
  real(mytype), Dimension(-1:nn+2), Intent( In ) :: ds
        Integer :: i,m
  real(mytype) :: c,d,deno
!
	Do i = mm,nn
          c = -(ds(i-1)+ds(i))*0.5e0
          d=(ds(i)+ds(i+1))*0.5e0
          deno=c*d*(c-d)
          temp(i,1,1,ib) = 0.0e0
          temp(i,2,1,ib) = c*c/deno
          temp(i,3,1,ib) = (d+c)*(d-c)/deno
          temp(i,4,1,ib) = -d*d/deno
          temp(i,5,1,ib) = 0.0e0
          do m = 1,5
            temp(i,m,2,ib) = temp(i,m,1,ib)
          End Do
        End Do
!
        Return
        End Subroutine SCHEME3
!
!******************************************************************************
        Subroutine SCHEME3_LES(temp,ds,mm,nn,jj,ib)
!******************************************************************************
!       THIS Subroutine CALCULATES THE GRADIENT AT CELL CENTER USING THREE POINTS
!       WHEN THE VARABLES ARE LOCATED AT CELL CENTER *****FOR LES PART***a**
!       jj ARGUMENT IS ONLY USED BECAUSE THE Dimension OF grdszz IS DIFFERENT FROM 
!	OTHRS TWO ARRAY Jj=2 FOR ZDIRECTION AND 1 FOR OTHERS
        Implicit None
        Integer, Intent( In ) :: nn,mm,ib,jj
  real(mytype), Dimension(mm:nn+jj,5,2), Intent( InOut ) :: temp
  real(mytype), Dimension(-1:(nn*2)+2), Intent( In ) :: ds
        Integer :: i,m
  real(mytype) :: c,d,deno
!
        Do i = mm,nn
          c = -(ds(2*I-3)+ds(2*I-2)+ds(2*I-1)+ds(2*I))*0.5e0
          d=(ds(2*I-1)+ds(2*I)+ds(2*I+1)+ds(2*I+2))*0.5e0
          deno=c*d*(c-d)
          temp(i,1,1) = 0.0e0
          temp(i,2,1) = c*c/deno
          temp(i,3,1) = (d+c)*(d-c)/deno
          temp(i,4,1) = -d*d/deno
          temp(i,5,1) = 0.0e0
!
!	********
!	grdsxx(i,1-5,**2**) IS FOR BUMP AND AT THE MOMENT WE ARE TAKING IT EQUAL TO
!	grdsxx(i,1-5,1)
!
          do m = 1,5
            temp(i,m,2) = temp(i,m,1)
          End Do
        End Do
!
        Return
        End Subroutine SCHEME3_LES
!
!******************************************************************************
        Subroutine SCHEME3_BC(temp1,temp2,ds,mm,nn,jj,ib)
!******************************************************************************
!       THIS Subroutine CALCULATES THE GRADIENT AT CELL CENTER USING THREE POINTS
!       WHEN THE VARABLES ARE LOCATED AT CELL CENTER FOR CONVECTIVE BOUNDARY CONDITIONS.
!       jj ARGUMENT IS ONLY USED BECAUSE THE Dimension OF grdszz IS DIFFERENT FROM 
!	OTHRS TWO ARRAY jj=2 FOR Z DIRECTION AND 1 FOR OTHERS
        Implicit None
        Integer, Intent( In ) :: nn,mm,ib,jj
  real(mytype), Dimension(mm-1:nn+1,5,ngr,nbl), Intent( InOut ) :: temp1
  real(mytype), Dimension(mm:nn/2+jj,5,2), Intent( InOut ) :: temp2
  real(mytype), Dimension(-1:nn+2), Intent( In ) :: ds
        Integer :: i,m
  real(mytype) :: c,d,deno
!
        i=mm
        c=(ds(i)+ds(i+1))*0.5e0
        d=ds(i+1)+0.5e0*(ds(i)+ds(i+2))
        deno=c*d*(c-d)
        temp1(i,1,1,ib) = c*c/deno
        temp1(i,2,1,ib) = -d*d/deno
        temp1(i,3,1,ib) = (d+c)*(d-c)/deno
        temp1(i,4,1,ib) = 0.0e0
        temp1(i,5,1,ib) = 0.0e0
        c = -ds(i)*0.5e0
        d=(ds(i)+ds(i+1))*0.5e0
        deno=c*d*(c-d)
        temp1(i,1,2,ib) = 0.0e0
        temp1(i,2,2,ib) = c*c/deno
        temp1(i,3,2,ib) = (d+c)*(d-c)/deno
        temp1(i,4,2,ib) = -d*d/deno !0.0e0        !Neumann B.c.
        temp1(i,5,2,ib) = 0.0e0
!
        deno=ds(i)*0.5e0
        temp1(i-1,1,2,ib) = 0.0e0
        temp1(i-1,2,2,ib) = 1.0e0/deno
        temp1(i-1,3,2,ib) = -1.0e0/deno
        temp1(i-1,4,2,ib) = 0.0e0
        temp1(i-1,5,2,ib) = 0.0e0
!
        i=mm
        c=(ds(1)+ds(2)+ds(3)+ds(4))*0.5e0
        d=ds(3)+ds(4)+0.5e0*(ds(1)+ds(2)+ds(5)+ds(6))
        deno=c*d*(c-d)
        temp2(i,1,ib) = c*c/deno
        temp2(i,2,ib) = -d*d/deno
        temp2(i,3,ib) = (d+c)*(d-c)/deno
        temp2(i,4,ib) = 0.0e0
        temp2(i,5,ib) = 0.0e0
!
        i=nn
        c = -(ds(i-1)+ds(i))*0.5e0
        d = -ds(i-1)-(ds(i)+ds(i-2))*0.5e0
        deno=c*d*(c-d)
        temp1(i,1,1,ib) = 0.0e0
        temp1(i,2,1,ib) = 0.0e0
        temp1(i,3,1,ib) = (d+c)*(d-c)/deno
        temp1(i,4,1,ib) = -d*d/deno
        temp1(i,5,1,ib) = c*c/deno
        c = -(ds(i-1)+ds(i))*0.5e0
        d=ds(i)*0.5e0
        deno=c*d*(c-d)
        temp1(i,1,2,ib) = 0.0e0
        temp1(i,2,2,ib) = c*c/deno !0.0e0     !Neumann B.c.
        temp1(i,3,2,ib) = (d+c)*(d-c)/deno
        temp1(i,4,2,ib) = -d*d/deno
        temp1(i,5,2,ib) = 0.0e0
!
        deno=ds(i)*0.5
        temp1(i+1,1,2,ib) = 0.0e0
        temp1(i+1,2,2,ib) = 0.0e0
        temp1(i+1,3,2,ib) = 1.0e0/deno
        temp1(i+1,4,2,ib) = -1.0e0/deno
        temp1(i+1,5,2,ib) = 0.0e0

!
        i=nn/2
        c = -(ds(nn-3)+ds(nn-2)+ds(nn-1)+ds(nn))*0.5e0
        d = -(ds(nn-2)+ds(nn-3))-(ds(nn)+ds(nn-1)+ds(nn-4)+ds(nn-5))*0.5e0
        deno=c*d*(c-d)
        temp2(i,1,ib) = 0.0e0
        temp2(i,2,ib) = 0.0e0
        temp2(i,3,ib) = (d+c)*(d-c)/deno
        temp2(i,4,ib) = -d*d/deno
        temp2(i,5,ib) = c*c/deno
!
        Return
        End Subroutine SCHEME3_BC
!
!******************************************************************************
	Subroutine SCHEME2_BC_BLOCK
!******************************************************************************
! 	THIS Subroutine DETERMINES THE GRADIENT AT CELL FACE AT BOUNDARIES WITH
!	VARibLES AT CELL CENTER(ibwest=1=VARIABLE AT BOUNDARY, ibwest = 0=AT CELL CENTER)
!	Ij=1 X_DIRECTION, Ij=2 Y-dIRECTION,Ij=3 Z DIRECTION
	Implicit None
	Integer :: i1,i2,ib,i,m
  real(mytype) :: a,b,c,d
        ib=2
!
!	Left surface Gradient
	If (ibegbl.gt.1) Then
        i2=ibegbl-1
	Do i = i2,i2+1
          a=dx(i)+0.5e0*dx(i+1)
          b = 0.5e0*dx(i)
          c = -0.5e0*dx(i-1)
          d = -dx(i-1)-0.5e0*dx(i-2)
	  If (i == i2) Then 
              a=dx(i)
              Call THIRD (i,a,b,c,d,1,1,ib)
	  Write(11,140)
	!Write(11,152)
	  Write(11,153) i,(gcdx_f(i,m,1,2),m=1,4)
	  Else If (i == i2+1) Then
              a = 0.0e0
              b = 0.0e0
              Call THIRD3 (i,a,b,c,d,1,2,1,ib)
            do m = 1,4
!	Neumann boundary condition for pressure gradient at boundary
              gcdx_f(i,m,2,ib) = 0.0e0    
            End Do
	!Write(11,152)
	  Write(11,151) i,(gcdx_f(i,m,1,2),m=1,4),(gcdx_f(i,m,2,2),m=1,4)
	  End If
        End Do
	End If
!
!	Lower surface Gradient
!
	If (jbegbl.gt.1) Then 
        i2=jbegbl-1
	Do i = i2,i2+1
          a=dy(i)+0.5e0*dy(i+1)
          b = 0.5e0*dy(i)
          c = -0.5e0*dy(i-1)
          d = -dy(i-1)-0.5e0*dy(i-2)
	  If (i == i2) Then 
             a=dy(i)
             Call THIRD (i,a,b,c,d,2,1,ib)
	!Write(11,152)
	  Write(11,141)
	  Write(11,153) i,(gcdy_f(i,m,1,2),m=1,4)
	  Else If (i == i2+1) Then
              a = 0.0e0
              b = 0.0e0
              Call THIRD3 (i,a,b,c,d,2,2,1,ib)
!	Neumann boundary condition for pressure gradient at boundary
            do m = 1,4
              gcdy_f(i,m,2,ib) = 0.0e0
            End Do
	!Write(11,152)
	  Write(11,151) i,(gcdy_f(i,m,1,2),m=1,4),(gcdy_f(i,m,2,2),m=1,4)
	  End If
        End Do
	End If
!
!	Lower index surface in z-direction
!
	!If (kbegbl.gt.1) Then 
        !i2=kbegbl-1
	!Do i = i2,i2+1
        !  a=dz(i)+0.5e0*dz(i+1)
        !  b = 0.5e0*dz(i)
        !  c = -0.5e0*dz(i-1)
        !  d = -dz(i-1)-0.5e0*dz(i-2)
	!  If (i == i2) Then
        !     a=dz(i)
        !     Call THIRD (i,a,b,c,d,3,1,ib)
        !  Else If (i == i2+1) Then
        !      a = 0.0e0
        !      b = 0.0e0
        !      Call THIRD3 (i,a,b,c,d,3,2,1,ib)
!	Neumann boundary condition for pressure gradient at boundary
        !    do m = 1,4
        !      gcdz_f(i,m,2,ib) = 0.0e0
        !    End Do
	!  End If
        !End Do
	!End If
!
!	Gradients for Right surface of the Block
!
	If (iendbl.lt.nxt) Then
        i1=iendbl
        Do i = i1,i1+1
          a=dx(i)+0.5e0*dx(i+1)
          b = 0.5e0*dx(i)
          c = -0.5e0*dx(i-1)
          d = -dx(i-1)-0.5e0*dx(i-2)
	!If (ij.eq.2.and.i.eq.ind) print *,'abcd',a,b,c,d
	  If (i == i1) Then
            c = 0.0e0
            d = 0.0e0
            Call THIRD3 (i,a,b,c,d,1,1,1,ib)
!	Neumann boundary condition for pressure gradient at boundary
            do m = 1,4
              gcdx_f(i,m,2,ib) = 0.0e0
            End Do
	!Write(11,150)
	  Write(11,142)
	  Write(11,151) i,(gcdx_f(i,m,1,2),m=1,4),(gcdx_f(i,m,2,2),m=1,4)
	  Else If (i == i1+1) Then
            d = -dx(i-1)
            Call THIRD (i,a,b,c,d,1,1,ib)
	!Write(11,150)
	  Write(11,153) i,(gcdx_f(i,m,1,2),m=1,4)
          End If
        End Do
	End If
!
!	Gradients for Upper surface of the Block
!
	If (jendbl.lt.nyt) Then
        i1=jendbl
        Do i = i1,i1+1
          a=dy(i)+0.5e0*dy(i+1)
          b = 0.5e0*dy(i)
          c = -0.5e0*dy(i-1)
          d = -dy(i-1)-0.5e0*dy(i-2)
	  If (i == i1) Then
            c = 0.0e0
            d = 0.0e0
            Call THIRD3 (i,a,b,c,d,2,1,1,ib)
!	Neumann boundary condition for pressure gradient at boundary
            do m = 1,4
              gcdy_f(i,m,2,ib) = 0.0e0
            End Do
	!Write(11,150)
	  Write(11,143)
	  Write(11,151) i,(gcdy_f(i,m,1,2),m=1,4),(gcdy_f(i,m,2,2),m=1,4)
          Else If (i == i1+1) Then
            d = -dy(i-1)
            Call THIRD (i,a,b,c,d,2,1,ib)
	!Write(11,150)
	  Write(11,153) i,(gcdy_f(i,m,1,2),m=1,4)
          End If
        End Do
	End If
!
!	Gradients for Upper index surface of the Block in spanwise direction
!
	!If (kendbl.lt.nzt) Then
        !i1=kendbl
        !Do i = i1,i1+1
        !  a=dz(i)+0.5e0*dz(i+1)
        !  b = 0.5e0*dz(i)
        !  c = -0.5e0*dz(i-1)
        !  d = -dz(i-1)-0.5e0*dz(i-2)
	!  If (i == i1) Then
        !    c = 0.0e0
        !    d = 0.0e0
        !    Call THIRD3 (i,a,b,c,d,3,1,1,ib)
!	Neumann boundary condition for pressure gradient at boundary
        !    do m = 1,4
        !      gcdz_f(i,m,2,ib) = 0.0e0
        !    End Do
        !  Else If (i == i1+1) Then
        !     d = -dz(i-1)
        !     Call THIRD (i,a,b,c,d,3,1,ib)
        !    End If
        !End Do
	!End If
!
  140   Format(' Gradients at cell face on  Left surface of block'//)
  141   Format(' Gradients at cell face on  Lower surface of block'//)
  142   Format(' Gradients at cell face on  Right surface of block'//)
  143   Format(' Gradients at cell face on  Upper surface of block'//)
!  150   Format(' Gradients at cell face on  Right/Upper surface of block'//)
!  152   Format(' Gradient at cell face  on Left/Lower surface of block'//)
  151   Format('Location='i,4x,'Graients for velocity are'/,4E14.6,/ &
               'Gradients for pressure are'/,4E14.6//)
  153   Format('Location='i,4x,'Graients for velocity are'/,4E14.6//)
!
	Return
	End Subroutine SCHEME2_BC_BLOCK
!
!******************************************************************************
        Subroutine SCHEME3_BC_BLOCK
!******************************************************************************
!       This procedure calculates the cell centered gradient near the block surface
!	when the variables asre located at cell center.
!	ONE SIDED DIFFERENCING
	Implicit None
	Integer :: i,ib,m
  real(mytype) :: a,b,c,d,deno
        ib=2
!
!	Gradient near the left surface of the block
!
	If (ibegbl.gt.1) Then
        i=ibegbl-1
!
        c=(dx(i)+dx(i+1))*0.5e0
        d=dx(i+1)+0.5e0*(dx(i)+dx(i+2))
        deno=c*d*(c-d)
        gcdx_c(i,1,1,ib) = c*c/deno
        gcdx_c(i,2,1,ib) = -d*d/deno
        gcdx_c(i,3,1,ib) = (d+c)*(d-c)/deno
        gcdx_c(i,4,1,ib) = 0.0e0
        gcdx_c(i,5,1,ib) = 0.0e0
        c = -dx(i)*0.5e0
        d=(dx(i)+dx(i+1))*0.5e0
        deno=c*d*(c-d)
        gcdx_c(i,1,2,ib) = 0.0e0
        gcdx_c(i,2,2,ib) = c*c/deno
        gcdx_c(i,3,2,ib) = (d+c)*(d-c)/deno
        gcdx_c(i,4,2,ib) = 0.0e0        !Neumann B.c.
        gcdx_c(i,5,2,ib) = 0.0e0
!
	Write(11,140)
	Write(11,151) i,(gcdx_c(i,m,1,1),m=1,5), (gcdx_c(i,m,2,1),m=1,5)
	End If
!
!	Gradient near the lower surface of the block
!
	If (jbegbl.gt.1) Then
        i=jbegbl-1
!
        c=(dy(i)+dy(i+1))*0.5e0
        d=dy(i+1)+0.5e0*(dy(i)+dy(i+2))
        deno=c*d*(c-d)
        gcdy_c(i,1,1,ib) = c*c/deno
        gcdy_c(i,2,1,ib) = -d*d/deno
        gcdy_c(i,3,1,ib) = (d+c)*(d-c)/deno
        gcdy_c(i,4,1,ib) = 0.0e0
        gcdy_c(i,5,1,ib) = 0.0e0
        c = -dy(i)*0.5e0
        d=(dy(i)+dy(i+1))*0.5e0
        deno=c*d*(c-d)
        gcdy_c(i,1,2,ib) = 0.0e0
        gcdx_c(i,2,2,ib) = c*c/deno
        gcdy_c(i,3,2,ib) = (d+c)*(d-c)/deno
        gcdy_c(i,4,2,ib) = 0.0e0        !Neumann B.c.
        gcdy_c(i,5,2,ib) = 0.0e0
!
	Write(11,141)
	Write(11,151) i,(gcdy_c(i,m,1,1),m=1,5), (gcdy_c(i,m,2,1),m=1,5)
	End If
!
!	Gradient near the right surface of the block
!
	If (iendbl.lt.nxt) Then
        i=iendbl
!
        c = -(dx(i-1)+dx(i))*0.5e0
        d = -dx(i-1)-(dx(i)+dx(i-2))*0.5e0
        deno=c*d*(c-d)
        gcdx_c(i,1,1,ib) = 0.0e0
        gcdx_c(i,2,1,ib) = 0.0e0
        gcdx_c(i,3,1,ib) = (d+c)*(d-c)/deno
        gcdx_c(i,4,1,ib) = -d*d/deno
        gcdx_c(i,5,1,ib) = c*c/deno
        c = -(dx(i-1)+dx(i))*0.5e0
        d=dx(i)*0.5e0
        deno=c*d*(c-d)
        gcdx_c(i,1,2,ib) = 0.0e0
        gcdx_c(i,2,2,ib) = 0.0e0     !Neumann B.c.
        gcdx_c(i,3,2,ib) = (d+c)*(d-c)/deno
        gcdx_c(i,4,2,ib) = -d*d/deno
        gcdx_c(i,5,2,ib) = 0.0e0
!
	Write(11,142)
	Write(11,151) i,(gcdx_c(i,m,1,1),m=1,5), (gcdx_c(i,m,2,1),m=1,5)
	End If
!
!	Gradient near the upper surface of the block
!
	If (jendbl.lt.nyt) Then
        i=jendbl
!
        c = -(dy(i-1)+dy(i))*0.5e0
        d = -dy(i-1)-(dy(i)+dx(i-2))*0.5e0
        deno=c*d*(c-d)
        gcdy_c(i,1,1,ib) = 0.0e0
        gcdy_c(i,2,1,ib) = 0.0e0
        gcdy_c(i,3,1,ib) = (d+c)*(d-c)/deno
        gcdy_c(i,4,1,ib) = -d*d/deno
        gcdy_c(i,5,1,ib) = c*c/deno
        c = -(dy(i-1)+dy(i))*0.5e0
        d=dy(i)*0.5e0
        deno=c*d*(c-d)
        gcdy_c(i,1,2,ib) = 0.0e0
        gcdy_c(i,2,2,ib) = 0.0e0     !Neumann B.c.
        gcdy_c(i,3,2,ib) = (d+c)*(d-c)/deno
        gcdy_c(i,4,2,ib) = -d*d/deno
        gcdy_c(i,5,2,ib) = 0.0e0
!
	Write(11,143)
	Write(11,151) i,(gcdy_c(i,m,1,1),m=1,5), (gcdy_c(i,m,2,1),m=1,5)
	End If
!
  140   Format(' Gradients at cell center near  Left surface of block'//)
  141   Format(' Gradients at cell center near  Lower surface of block'//)
  142   Format(' Gradients at cell center near  Right surface of block'//)
  143   Format(' Gradients at cell center near  Upper surface of block'//)
  151   Format('Location='i,4x,/'Graients for velocity are'/,5E14.6,/ &
              'Gradients for pressure are'/,5E14.6//)
!
	Return
	End Subroutine SCHEME3_BC_BLOCK
!
	End Module cal_coefficient
