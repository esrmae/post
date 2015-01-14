	MODULE calc_post
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
	SUBROUTINE calc_perform
!
	If ((isteady == 1).and.(lpstat_1d == 1)) then
	    Call calc_steady
	Else If (isteady == 0) then 
	    Call calc_unsteady   
            If (hybrid == 1) Call hybrid_calc
	End If
!
	END SUBROUTINE calc_perform
!*******************************************************************************
	SUBROUTINE calc_steady
        REAL(mytype), DIMENSION(:,:,:), ALLOCATABLE :: sqar_avg
	Allocate(sqar_avg(-1:ny+2,4,4))

        sqar_avg=0.0e0
        !sqar_avg(:,3,1)=omega_avg(:,4)**2
        sqar_avg(:,:,2)=stat_avg(:,:,1)**2
        sqar_avg(:,:,3)=stat_avg(:,:,1)**3
        sqar_avg(:,:,4)=stat_avg(:,:,1)**4

	If (iflucs == 1) CALL flucs_calc(stat_avg(-1:ny+2,:,:),sqar_avg(-1:ny+2,:,:), &
                                         flucs(-1:ny+2,:,:))
	If (ivort == 1) then
         ! CALL vortmn1d(stat_avg(-1:ny+2,:,:),sqar_avg(-1:ny+2,:,:),omega_avg(-1:ny+2,:))
          CALL vort_calc(sqar_avg(-1:ny+2,:,:),omega_avg(-1:ny+2,:))
        End If
	If (itke == 1) CALL tke_calc(stat_avg(-1:ny+2,:,:),quad_avg(-1:ny+2,:), &
                                     tke_avg(-1:ny+2,:,:),kinen(-1:ny+2,:,:))
	If (ishear == 1) CALL shear_calc(stat_avg(-1:ny+2,:,:),quad_avg(-1:ny+2,:), &
                                         shear(-1:ny+2,:))
	If (ianis == 1) CALL anis_calc(stat_avg(-1:ny+2,:,:),quad_avg(-1:ny+2,:), &
                                       anisb(-1:ny+2,:),anisII(-1:ny+2), &
                                       anisIII(-1:ny+2),anisf(-1:ny+2),anisg(-1:ny+2))
	If (iflucs == 1) CALL flucs_div(flucs(-1:ny+2,:,:))
	If (iangle == 1) CALL angle_calc(stat_avg(-1:ny+2,:,:),quad_avg(-1:ny+2,:),sqar_avg(-1:ny+2,:,:),angle(-1:ny+2,:))

	If ((half_av == 1)) CALL avg_walls(stat_avg(-1:ny+2,:,:),flucs(-1:ny+2,:,:), &
                                           tke_avg(-1:ny+2,:,:),wx_avg(-1:ny+2,:,:), &
			                   omega_avg(-1:ny+2,:),kinen(-1:ny+2,:,:),  &
                                           anisb(-1:ny+2,:),anisII(-1:ny+2),anisIII(-1:ny+2), &
                                           anisf(-1:ny+2),anisg(-1:ny+2),shear(-1:ny+2,:))
!
	END SUBROUTINE calc_steady
!*******************************************************************************
	SUBROUTINE calc_unsteady
	Integer :: i
!
	If (iflucs == 1) then
	  Do i=1,tpts
	    Call flucs_calc(stat1d_sa(-1:ny+2,:,:,i),sqar1d_sa(-1:ny+2,:,:,i), &
                            flucs1d_sa(-1:ny+2,:,:,i))
          End Do
        End If
	If (ivort == 1) then
	  Do i=1,tpts
	    Call vort_calc(sqar1d_sa(-1:ny+2,:,:,i),omega1d_sa(-1:ny+2,:,i))
          End Do
        End If
	If (itke == 1) then
	  Do i=1,tpts
	    Call tke_calc(stat1d_sa(-1:ny+2,:,:,i),quad1d_sa(-1:ny+2,:,i), &
			  tke1d_sa(-1:ny+2,:,:,i),kinen1d_sa(-1:ny+2,:,:,i))
          End Do
        End If
	If (ishear == 1) then
	  Do i=1,tpts
	    Call shear_calc(stat1d_sa(-1:ny+2,:,:,i),quad1d_sa(-1:ny+2,:,i), &
			    shear1d_sa(-1:ny+2,:,i))
          End Do
        End If
	If (ianis == 1) then
	  Do i=1,tpts
	    Call anis_calc(stat1d_sa(-1:ny+2,:,:,i),quad1d_sa(-1:ny+2,:,i), &
			   anisb1d_sa(-1:ny+2,:,i),anisII1d_sa(-1:ny+2,i), &
			   anisIII1d_sa(-1:ny+2,i),anisf1d_sa(-1:ny+2,i),anisg1d_sa(-1:ny+2,i))
          End Do
        End If
	If (iangle == 1) then
	  Do i=1,tpts
	    Call angle_calc(stat1d_sa(-1:ny+2,:,:,i),quad1d_sa(-1:ny+2,:,i), &
			    sqar1d_sa(-1:ny+2,:,:,i),angle1d_sa(-1:ny+2,:,i))
          End Do
        End If
	If (iflucs == 1) then
	  Do i=1,tpts
	    Call flucs_div(flucs1d_sa(-1:ny+2,:,:,i))
          End Do
        End If

	If (half_av == 1) then
	  Do i=1,tpts
	    Call avg_walls(stat1d_sa(-1:ny+2,:,:,i),flucs1d_sa(-1:ny+2,:,:,i), &
			   tke1d_sa(-1:ny+2,:,:,i),wx1d_sa(-1:ny+2,:,:,i), &
			   omega1d_sa(-1:ny+2,:,i),kinen1d_sa(-1:ny+2,:,:,i), &
			   anisb1d_sa(-1:ny+2,:,i),anisII1d_sa(-1:ny+2,i), &
			   anisIII1d_sa(-1:ny+2,i),anisf1d_sa(-1:ny+2,i), &
			   anisg1d_sa(-1:ny+2,i),shear1d_sa(-1:ny+2,:,i))
          End Do
        End If
!
	END SUBROUTINE calc_unsteady
!*******************************************************************************
	SUBROUTINE flucs_calc(stat_in,sqar_in,flucs_out)
!	This subroutine calculates the fluctuations of different orders in the u direction
	Integer :: j,m
	Real(mytype), dimension(-1:ny+2,4,4), intent(in) :: stat_in
	Real(mytype), dimension(-1:ny+2,4,4), intent(in) :: sqar_in
	Real(mytype), dimension(-1:ny+2,4,4), intent(out) :: flucs_out
!
!	We then calculate the tke values
        DO m=1,4 !for u,v,w,p
          DO j=0,nyt
            flucs_out(j,m,1) = stat_in(j,m,1)
            flucs_out(j,m,2) = stat_in(j,m,2)-sqar_in(j,m,2)
            flucs_out(j,m,3) = stat_in(j,m,3)-sqar_in(j,m,3) &
                               -3.0*stat_in(j,m,1)*flucs_out(j,m,2)
            flucs_out(j,m,4) = stat_in(j,m,4)-sqar_in(j,m,4) & 
                               -6.0*sqar_in(j,m,2)*flucs_out(j,m,2) &
                               -4.0*stat_in(j,m,1)*flucs_out(j,m,3)
          END DO
        END DO
!
        END SUBROUTINE flucs_calc
!*******************************************************************************
        SUBROUTINE flucs_div(flucs_out)
!       This subroutine sorts out the fluctuations divide by zero
        Integer :: j,m
        Real(mytype), dimension(-1:ny+2,4,4), intent(out) :: flucs_out

        DO m=1,4
            If (flucs_out(0,m,2)==0.0e0) flucs_out(0,m,3)=flucs_out(1,m,3)*small**1.5e0/flucs_out(1,m,2)**1.5e0
            If (flucs_out(nyt,m,2)==0.0e0) flucs_out(nyt,m,3)=flucs_out(ny,m,3)*small**1.5e0/flucs_out(ny,m,2)**1.5e0
            If (flucs_out(0,m,2)==0.0e0) flucs_out(0,m,4)=flucs_out(1,m,4)*small**2.0e0/flucs_out(1,m,2)**2.0e0
            If (flucs_out(nyt,m,2)==0.0e0) flucs_out(nyt,m,4)=flucs_out(ny,m,4)*small**2.0e0/flucs_out(ny,m,2)**2.0e0
        END DO

        DO m=1,4
          DO j=0,nyt
            If (flucs_out(j,m,2)==0.0e0) flucs_out(j,m,2)=small
          END DO
        END DO
!
        END SUBROUTINE flucs_div
!*******************************************************************************
	SUBROUTINE vortmn1d(stat_in,sqar_out,omega_in)
!	This subroutine calculates the mean vorticity for 1d only
	Integer :: j,m,ib
	Real(mytype), dimension(-1:ny+2,4,4), intent(in) :: stat_in
	Real(mytype), dimension(-1:ny+2,6), intent(in) :: omega_in
	Real(mytype), dimension(-1:ny+2,4,4), intent(inout) :: sqar_out
        ib=1
!
!	We then calculate the mean vorticity values
!
        DO j=0,nyt
          sqar_out(j,1,1)=gcdy_c(j,2,2,ib)*stat_in(j+1,3,1)+ &
                          gcdy_c(j,3,2,ib)*stat_in(j,3,1)+ &
                          gcdy_c(j,4,2,ib)*stat_in(j-1,3,1)
          sqar_out(j,3,1)=omega_in(j,6)!-(gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)+ &
                            !gcdy_c(j,3,2,ib)*stat_in(j,1,1)+ &
                            !gcdy_c(j,4,2,ib)*stat_in(j-1,1,1))
          sqar_out(j,1,1)=sqar_out(j,1,1)**2
          sqar_out(j,3,1)=sqar_out(j,3,1)**2
        END DO

        END SUBROUTINE    
!*******************************************************************************
	SUBROUTINE vort_calc(sqar_in,omega_out)
!	This subroutine calculates the fluctuations of different orders in the u direction
	Integer :: j,m,ib
	Real(mytype), dimension(-1:ny+2,4,4), intent(in) :: sqar_in
	Real(mytype), dimension(-1:ny+2,6), intent(inout) :: omega_out
!
!	We then calculate the mean vorticity values
!
	If(mpi.Eq.0) Then
          DO j=0,nyt
            omega_out(j,1)=sqrt(abs(omega_out(j,1)))
            omega_out(j,2)=sqrt(abs(omega_out(j,2)))
            omega_out(j,3)=sqrt(abs(omega_out(j,3)-omega_out(j,4)**2))
          END DO
	Else
          DO m=1,3
            DO j=0,nyt
              omega_out(j,m)=sqrt(abs(omega_out(j,m)-omega_out(j,m+3)**2))
            END DO
          END DO
	End If
!
	END SUBROUTINE vort_calc
!*******************************************************************************
	SUBROUTINE tke_calc(stat_in,quad_in,tke_in,kinen_out)
!	This subroutine calculates the tke values
	Integer :: j,m,n,i,ib
	Real(mytype), dimension(-1:ny+2,4,4), intent(in) :: stat_in   
	Real(mytype), dimension(-1:ny+2,9), intent(in) :: quad_in
	Real(mytype), dimension(-1:ny+2,4,7), intent(in) :: tke_in     
	Real(mytype), dimension(-1:ny+2,4,6), intent(out) :: kinen_out 
	Real(mytype), dimension(-1:ny+2) :: ddy,d2dy2
	ib = 1
!
	! Calculates P11=-2*uv*dU/dy & P12=-vv*dU/dy
	DO j=0,nyt
	    	 ddy(j)=gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)+ &
			gcdy_c(j,3,2,ib)*stat_in(j,1,1)+ &
			gcdy_c(j,4,2,ib)*stat_in(j-1,1,1)
	    kinen_out(j,1,1)=-2*quad_in(j,5)*ddy(j)
	    kinen_out(j,4,1)=-stat_in(j,2,2)*ddy(j)
	END DO
!
! 	Calculates Tii=-ddy(ui^2u2)
	DO m=1,3 !for u,v,w
	    DO j=0,nyt
		kinen_out(j,m,2)= -(gcdy_c(j,2,2,ib)*tke_in(j+1,m,2)+ &
				gcdy_c(j,3,2,ib)*tke_in(j,m,2)+ &
				gcdy_c(j,4,2,ib)*tke_in(j-1,m,2))
	   END DO
	END DO
! 	Calculates T12=-ddy(uv^2)
	DO j=0,nyt
	    kinen_out(j,4,2)= -(gcdy_c(j,2,2,ib)*tke_in(j+1,2,1)+ &
			  gcdy_c(j,3,2,ib)*tke_in(j,2,1)+ &
			  gcdy_c(j,4,2,ib)*tke_in(j-1,2,1))
	END DO
!
! 	Calculates PI22=-(2/rho)ddy(pv)
	DO j=0,nyt
	    kinen_out(j,2,3)=-(2/rho)*(gcdy_c(j,2,2,ib)*tke_in(j+1,2,3)+ &
			  	gcdy_c(j,3,2,ib)*tke_in(j,2,3)+ &
			  	gcdy_c(j,4,2,ib)*tke_in(j-1,2,3))
	END DO
! 	Calculates PI12=-(1/rho)ddy(pu)
	DO j=0,nyt
	    kinen_out(j,4,3)=-(1/rho)*(gcdy_c(j,2,2,ib)*tke_in(j+1,1,3)+ &
			  	gcdy_c(j,3,2,ib)*tke_in(j,1,3)+ &
			  	gcdy_c(j,4,2,ib)*tke_in(j-1,1,3))
	END DO
!
! 	Calculates PHIii=(2/rho)p(du_i/dx_i)
	DO m=1,3 !for u,v,w
	    DO j=0,nyt
	    	kinen_out(j,m,4)=(2/rho)*tke_in(j,m,4)
	    END DO
	END DO
	!because a bug in MPI main code causes the value of dudx with a opposite sign
	If(mpi.Eq.1) kinen_out(:,1,4)=-kinen_out(:,1,4)
! 	Calculates PHI12=(1/rho)(p*du/dy+p*dv/dx)
	DO j=0,nyt
	    kinen_out(j,4,4)=(1/rho)*(tke_in(j,4,3)+tke_in(j,4,4))
	END DO
!
! 	Calculates Dii=nu*d^2/dy^2(u_i^2)
	DO m=1,3 !for u,v,w
	    DO j=0,nyt
		ddy(j)=(gcdy_c(j,2,2,ib)*stat_in(j+1,m,2)+ &
			  gcdy_c(j,3,2,ib)*stat_in(j,m,2)+ &
			  gcdy_c(j,4,2,ib)*stat_in(j-1,m,2))
	    END DO
!
	    d2dy2(0)=ddy(1)/(0.5e0*dys(1))
	    kinen_out(0,m,5)=nu*d2dy2(0)
            d2dy2(1)=(ddy(2)- ddy(1))/dys(2)
            kinen_out(1,m,5)=nu*d2dy2(1)
	    DO j=2,ny-1
		d2dy2(j)=(gcdy_c(j,2,2,ib)*ddy(j+1)+ &
			  gcdy_c(j,3,2,ib)*ddy(j)+ &
			  gcdy_c(j,4,2,ib)*ddy(j-1))
	    	kinen_out(j,m,5)=nu*d2dy2(j)
	    END DO
            d2dy2(ny)=(ddy(ny)-ddy(ny-1))/dys(ny)
            kinen_out(ny,m,5)=nu*d2dy2(ny)
	    d2dy2(nyt)=(-ddy(ny))/(0.5e0*dys(nyt))
	    kinen_out(nyt,m,5)=nu*d2dy2(nyt)
	END DO
! 	Calculates D12=nu*d^2/dy^2(uv)
	DO j=0,nyt
	    ddy(j)=(gcdy_c(j,2,2,ib)*quad_in(j+1,5)+ &
			  gcdy_c(j,3,2,ib)*quad_in(j,5)+ &
			  gcdy_c(j,4,2,ib)*quad_in(j-1,5))
	END DO
	d2dy2(0)=(ddy(1)- ddy(0))/dy(1)
	kinen_out(0,4,5)=nu*d2dy2(0)
	DO j=1,ny
	    d2dy2(j)=(gcdy_c(j,2,2,ib)*ddy(j+1)+ &
			  gcdy_c(j,3,2,ib)*ddy(j)+ &
			  gcdy_c(j,4,2,ib)*ddy(j-1))
	    kinen_out(j,4,5)=nu*d2dy2(j)
	END DO
	d2dy2(nyt)=(ddy(nyt)- ddy(ny))/dy(ny)
	kinen_out(nyt,4,5)=nu*d2dy2(nyt)
!
! 	Calculates epsilonii=-(2*nu)(du_i/dx_k)^2 (sum over k)
	DO m=1,3 !for u,v,w
	    DO j=0,nyt
	    	kinen_out(j,m,6)=-2*nu*(tke_in(j,m,5)+tke_in(j,m,6)+tke_in(j,m,7))
	    END DO
	END DO
! 	Calculates epsilon12=-(2*rho)(du_1/dx_k*du_2/dx_k)(sum over k)
	DO j=0,nyt
	    kinen_out(j,4,6)=-2*nu*(tke_in(j,4,5)+tke_in(j,4,6)+tke_in(j,4,7))
	END DO
!
!	If means are not subtracted in code, this sorts it out:
	If (tkemn == 1) then
	   DO j=0,nyt
!		T11
		kinen_out(j,1,2)= kinen_out(j,1,2)+2.0*(gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)*quad_in(j+1,5)+ &
				gcdy_c(j,3,2,ib)*stat_in(j,1,1)*quad_in(j,5)+ &
				gcdy_c(j,4,2,ib)*stat_in(j-1,1,1)*quad_in(j-1,5))
!		T12
		kinen_out(j,4,2)= kinen_out(j,4,2)+(gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)*stat_in(j+1,2,2)+ &
				gcdy_c(j,3,2,ib)*stat_in(j,1,1)*stat_in(j,2,2)+ &
				gcdy_c(j,4,2,ib)*stat_in(j-1,1,1)*stat_in(j-1,2,2))
!		PI12
		kinen_out(j,4,3)= kinen_out(j,4,3)+(1/rho)*(gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)*stat_in(j+1,4,1)+ &
				gcdy_c(j,3,2,ib)*stat_in(j,1,1)*stat_in(j,4,1)+ &
				gcdy_c(j,4,2,ib)*stat_in(j-1,1,1)*stat_in(j-1,4,1))
!		PHI12
		kinen_out(j,4,4)= kinen_out(j,4,4)-(1/rho)*stat_in(j,4,1)*(gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)+ &
				gcdy_c(j,3,2,ib)*stat_in(j,1,1) + gcdy_c(j,4,2,ib)*stat_in(j-1,1,1))	
!		D11 (finds d/dy)
	    	ddy(j)=(gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)**2 + gcdy_c(j,3,2,ib)*stat_in(j,1,1)**2+ &
			  gcdy_c(j,4,2,ib)*stat_in(j-1,1,1)**2)
!		E11
		kinen_out(j,1,6)= kinen_out(j,1,6)+2.0*nu*(gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)+ &
				gcdy_c(j,3,2,ib)*stat_in(j,1,1) + gcdy_c(j,4,2,ib)*stat_in(j-1,1,1))**2 		
!	     For unsteady case
	     IF(isteady == 0) Then
!	       P33
	       kinen_out(j,3,1) = kinen_out(j,3,1)+2.0*(gcdy_c(j,2,2,ib)*stat_in(j+1,3,1)*quad_in(j+1,6)+ &
				gcdy_c(j,3,2,ib)*stat_in(j,3,1)*quad_in(j,6)+ &
				gcdy_c(j,4,2,ib)*stat_in(j-1,3,1)*quad_in(j-1,6))
	     END IF
	   END DO
!
!		D11 (finds d^2/dy^2)
	    d2dy2(0)=ddy(1)/(0.5e0*dys(1))
	    kinen_out(0,1,5)=kinen_out(0,1,5)-nu*d2dy2(0)
            d2dy2(1)=(ddy(2)- ddy(1))/dys(2)
            kinen_out(1,1,5)=kinen_out(1,1,5)-nu*d2dy2(1)

            d2dy2(ny)=(ddy(ny)-ddy(ny-1))/dys(ny)
            kinen_out(ny,1,5)=kinen_out(ny,1,5)-nu*d2dy2(ny)
	    d2dy2(nyt)=(-ddy(ny))/(0.5e0*dys(nyt))
	    kinen_out(nyt,1,5)=kinen_out(nyt,1,5)-nu*d2dy2(nyt)

	   DO j=2,ny-1
                d2dy2(j) = (gcdy_c(j,2,2,ib)*ddy(j+1)+ &
                            gcdy_c(j,3,2,ib)*ddy(j) + gcdy_c(j,4,2,ib)*ddy(j-1))
                 kinen_out(j,1,5)= kinen_out(j,1,5)-nu*d2dy2(j)
	   END DO
!
!	   PIij=-1/rho*(ui*dp/dxj+uj*dp/dxi)=Phiij+phiij : velocity pressure gradient term	
!	   DO j=0,nyt; DO m=1,4
!	     kinen_out(j,m,3)=kinen_out(j,m,3)+kinen_out(j,m,4)
!	   END DO; END DO
	End If
!
	END SUBROUTINE tke_calc
!*******************************************************************************
	SUBROUTINE anis_calc(stat_in,quad_in,anisb_out,anisII_out, &
			     anisIII_out,anisf_out,anisg_out)
!	This subroutine calculates the anisotrophy values bij,I,II,F,G
	Integer :: j,ib
	Real(mytype) :: u,v,w,uvw,den
	Real(mytype), dimension(-1:ny+2,4,4), intent(in) :: stat_in   
	Real(mytype), dimension(-1:ny+2,9), intent(in) :: quad_in
	Real(mytype), dimension(-1:ny+2,4), intent(out) :: anisb_out
	Real(mytype), dimension(-1:ny+2), intent(out) :: anisII_out,anisIII_out,anisf_out,anisg_out
	ib = 1
!
	DO j=0,nyt
	    u = stat_in(j,1,2)-stat_in(j,1,1)**2
	    v = (stat_in(j,2,2)+stat_in(j+1,2,2))/2.0-stat_in(j,2,1)**2
	    w = stat_in(j,3,2) - stat_in(j,3,1)**2
            uvw=u+v+w
            if (uvw == 0.0e0) uvw=small
!	    Calculates b11	
	    anisb_out(j,1)=u/uvw-1.0/3.0
!	    Calculates b22
	    anisb_out(j,2)=v/uvw-1.0/3.0
!	    Calculates b33
	    anisb_out(j,3)=w/uvw-1.0/3.0
!	    Calculates b12
	    anisb_out(j,4)=(quad_in(j,5)-stat_in(j,1,1)*stat_in(j,2,1))/uvw
!
!	    Calculates II and III
	    anisII_out(j)=-0.5*(anisb_out(j,1)**2+anisb_out(j,2)**2+anisb_out(j,3)**2+ &
				2.0*anisb_out(j,4)**2)
	    anisIII_out(j)=(1.0/3.0)*(anisb_out(j,1)**3+anisb_out(j,2)**3+anisb_out(j,3)**3+ &
				3.0*anisb_out(j,1)*anisb_out(j,4)**2+3.0*anisb_out(j,2)*anisb_out(j,4)**2)
!
!	    Calculates F and G
	    anisf_out(j) = 1 + 9.0*anisII_out(j) + 27.0*anisIII_out(j)
            den=(anisII_out(j)/3.0)**3
            if (den==0.0e0) den=small
	    anisg_out(j) = -(anisIII_out(j)/2.0)**2/den
        END DO

        anisb_out(0,:)=anisb_out(1,:)
        anisII_out(0)=anisII_out(1)
        anisIII_out(0)=anisIII_out(1)
        anisf_out(0)=anisf_out(1)
        anisg_out(0)=anisg_out(1)

        anisb_out(nyt,:)=anisb_out(ny,:)
        anisII_out(nyt)=anisII_out(ny)
        anisIII_out(nyt)=anisIII_out(ny)
        anisf_out(nyt)=anisf_out(ny)
        anisg_out(nyt)=anisg_out(ny)
!
	END SUBROUTINE anis_calc
!*******************************************************************************
	SUBROUTINE shear_calc(stat_in,quad_in,shear_out)
	Integer :: i,j,ib
	Real(mytype), dimension(-1:ny+2,4,4), intent(in) :: stat_in   
	Real(mytype), dimension(-1:ny+2,9), intent(in) :: quad_in
	Real(mytype), dimension(-1:ny+2,4), intent(out) :: shear_out
	Real(mytype), dimension(-1:ny+2) :: dudy,dwdy
	ib=1
!
!	Calculates dudy
	DO j=0,nyt
	    dudy(j)=gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)+ &
			gcdy_c(j,3,2,ib)*stat_in(j,1,1)+ &
			gcdy_c(j,4,2,ib)*stat_in(j-1,1,1)
!	Calculates tau_turbulent
	    shear_out(j,1)=-quad_in(j,5)
!	Calculates tau_laminar
	    shear_out(j,2)=nu*dudy(j)
!!	Calculates tau_sgs
!	    shear_out(j,3)=quad_in(j,7)*dudy(j)
!!	Calculates tau_total
!	    shear_out(j,4)=-quad_in(j,5)+(nu+quad_in(j,7))*dudy(j)
!	Calculates spanwise turbulent shear stress -vw
	    shear_out(j,3)=-quad_in(j,6)  
!	Calculates dWdy
	    dwdy(j)=gcdy_c(j,2,2,ib)*stat_in(j+1,3,1)+ &
			gcdy_c(j,3,2,ib)*stat_in(j,3,1)+ &
			gcdy_c(j,4,2,ib)*stat_in(j-1,3,1)
	    shear_out(j,4)=nu*dwdy(j)	    	    
	END DO
!        DO j=0,nyt
!          write(*,*) yc(j),shear_out(j,1),shear_out(j,2),shear_out(j,2)
!        END DO
!
	END SUBROUTINE shear_calc
!*******************************************************************************
	SUBROUTINE angle_calc(stat_in,quad_in,sqar_in,angle_out)
	Integer :: i,j,ib
	Real(mytype), dimension(-1:ny+2,4,4), intent(in) :: stat_in,sqar_in
	Real(mytype), dimension(-1:ny+2,9), intent(in) :: quad_in
	Real(mytype), dimension(-1:ny+2,5), intent(out) :: angle_out
	Real(mytype) :: tmp,tmp1,tmp2,tmp3
	ib=1
!
	DO j=0,nyt
	!mean velocity angle -> beta=atan(W/U)
	  If(j.Le.ny/2) Then
	    tmp=(stat_in(j,3,1)-stat_in(0,3,1))/stat_in(j,1,1)
	  Else
	    tmp=(stat_in(j,3,1)-stat_in(nyt,3,1))/stat_in(j,1,1)
	  End If
	  angle_out(j,1)=Atan(tmp)
	!mean velocity gradient angle -> alpha=atan(dWdy/dUdy)
	  tmp1=gcdy_c(j,2,2,ib)*stat_in(j+1,3,1)+ &
			gcdy_c(j,3,2,ib)*stat_in(j,3,1)+ &
			gcdy_c(j,4,2,ib)*stat_in(j-1,3,1)
	  tmp2=gcdy_c(j,2,2,ib)*stat_in(j+1,1,1)+ &
			gcdy_c(j,3,2,ib)*stat_in(j,1,1)+ &
			gcdy_c(j,4,2,ib)*stat_in(j-1,1,1)
	  angle_out(j,2)=Atan(tmp1/tmp2)
	! turbulent shear stress angle -> eta=atan(vw/uv)
	  tmp=quad_in(j,6)/quad_in(j,5)
	  angle_out(j,3)=Atan(tmp)
	! turbulent intensity angle -> gamma=atan(2uw/(u^2-w^2))
	  tmp1=stat_in(j,1,2)-sqar_in(j,1,2)
	  tmp2=stat_in(j,3,2)-sqar_in(j,3,2)
	  tmp3=quad_in(j,7)-stat_in(j,1,1)*stat_in(j,3,1)
	  tmp=2.0e0*tmp3/(tmp1-tmp2)
	  angle_out(j,4)=Atan(tmp)
	! structure parameter -> a1=sqrt(uv^2+vw^2)/(u^2+v^2+w^2)
	  tmp1=quad_in(j,5)**2+quad_in(j,6)**2
	  tmp2=0.0e0
	  DO i=1,3
	    tmp2=tmp2+(stat_in(j,i,2)-sqar_in(j,i,2))
	  END DO
	  angle_out(j,5)=Atan(Sqrt(tmp1)/tmp2)
	END DO
	! update edge
	DO i=1,5
	  angle_out(0,i)=angle_out(1,i)
	  angle_out(nyt,i)=angle_out(nyt-1,i)
	END DO
	angle_out(:,1:4)=angle_out(:,1:4)/pi*180.0e0
!
	END SUBROUTINE angle_calc
!*******************************************************************************
	SUBROUTINE avg_walls(stat_out,flucs_out,tke_out,wx_out,omega_out, &
			     kinen_out,anisb_out,anisII_out,anisIII_out, &
			     anisf_out,anisg_out,shear_out)
!
!	This subroutine calculates the avg between the walls
	Integer :: j
	Real(mytype), dimension(-1:ny+2,4,4), intent(inout) :: stat_out  
	Real(mytype), dimension(-1:ny+2,4,4), intent(inout) :: flucs_out 
	Real(mytype), dimension(-1:ny+2,4,7), intent(inout) :: tke_out  
	Real(mytype), dimension(-1:ny+2,2,8), intent(inout) :: wx_out
	Real(mytype), dimension(-1:ny+2,6), intent(inout) :: omega_out
	Real(mytype), dimension(-1:ny+2,4,6), intent(inout) :: kinen_out 
	Real(mytype), dimension(-1:ny+2,4), intent(inout) :: anisb_out
	Real(mytype), dimension(-1:ny+2), intent(inout) :: anisII_out,anisIII_out, &
						anisf_out,anisg_out
	Real(mytype), dimension(-1:ny+2,4), intent(inout) :: shear_out
!
!	Sets bottom of plane as sum
   DO j=0,ceiling(real(ny)/2.0)
     stat_out(j,:,:)=stat_out(j,:,:)+stat_out(ny-j+1,:,:)
     If (iflucs == 1) then
       flucs_out(j,:,1:2)=flucs_out(j,:,1:2)+flucs_out(ny-j+1,:,1:2)
       flucs_out(j,1,3)=flucs_out(j,1,3)+flucs_out(ny-j+1,1,3)
       flucs_out(j,2,3)=flucs_out(j,2,3)-flucs_out(ny-j+1,2,3)
       flucs_out(j,3:4,3)=flucs_out(j,3:4,3)+flucs_out(ny-j+1,3:4,3)
       flucs_out(j,:,4)=flucs_out(j,:,4)+flucs_out(ny-j+1,:,4)
     End If
     If (itke == 1) tke_out(j,:,:)=tke_out(j,:,:)+tke_out(ny-j+1,:,:)
     If (ivort == 1) then
       wx_out(j,:,:)=wx_out(j,:,:)+wx_out(ny-j+1,:,:)
       omega_out(j,1:5)=omega_out(j,1:5)+omega_out(ny-j+1,1:5)
       omega_out(j,6)=omega_out(j,6)-omega_out(ny-j+1,6)
     End If
     If (itke == 1) kinen_out(j,:,:)=kinen_out(j,:,:)+kinen_out(ny-j+1,:,:)
     If (ianis == 1) then
       anisb_out(j,:)=anisb_out(j,:)+anisb_out(ny-j+1,:)
       anisII_out(j)=anisII_out(j)+anisII_out(ny-j+1)
       anisIII_out(j)=anisIII_out(j)+anisIII_out(ny-j+1)
       anisf_out(j)=anisf_out(j)+anisf_out(ny-j+1)
       anisg_out(j)=anisg_out(j)+anisg_out(ny-j+1)
     End If
     If (ishear == 1) shear_out(j,:)=shear_out(j,:)-shear_out(ny-j+1,:)
   END DO
!
!	Sets top of plane equal to bottom
   DO j=0,ceiling(real(ny)/2.0)
     stat_out(ny-j+1,:,:)=stat_out(j,:,:)
     If (iflucs == 1) then
       flucs_out(ny-j+1,:,1:2)=flucs_out(j,:,1:2)
       flucs_out(ny-j+1,1,3)=flucs_out(j,1,3)
       flucs_out(ny-j+1,2,3)=-flucs_out(j,2,3)
       flucs_out(ny-j+1,3:4,3)=flucs_out(j,3:4,3)
       flucs_out(ny-j+1,:,4)=flucs_out(j,:,4)
     End If
     If (itke == 1) tke_out(ny-j+1,:,:)=tke_out(j,:,:)
     If (ivort == 1) then
       wx_out(ny-j+1,:,:)=wx_out(j,:,:)
       omega_out(ny-j+1,1:5)=omega_out(j,1:5)
       omega_out(ny-j+1,6)=-omega_out(j,6)
     End If
     If (itke == 1) kinen_out(ny-j+1,:,:)=kinen_out(j,:,:)
     If (ianis == 1) then
       anisb_out(ny-j+1,:)=anisb_out(j,1:)
       anisII_out(ny-j+1)=anisII_out(j)
       anisIII_out(ny-j+1)=anisIII_out(j)
       anisf_out(ny-j+1)=anisf_out(j)
       anisg_out(ny-j+1)=anisg_out(j)
    End If
    If (ishear == 1) shear_out(ny-j+1,:)=-shear_out(j,:)
   END DO
!
!	Average by dividing by 2
   stat_out=stat_out/2.0e0
   If (iflucs == 1) flucs_out=flucs_out/2.0e0
   If (itke == 1) then
     tke_out=tke_out/2.0e0
     kinen_out=kinen_out/2.0e0
   End If
   If (ivort == 1) then
     wx_out=wx_out/2.0e0
     omega_out=omega_out/2.0e0
   End If
   If (ianis == 1) then
	    anisb_out=anisb_out/2.0e0
	    anisII_out=anisII_out/2.0e0
	    anisIII_out=anisIII_out/2.0e0
	    anisf_out=anisf_out/2.0e0
	    anisg_out=anisg_out/2.0e0
   End If
   If (ishear == 1) shear_out=shear_out/2.0e0
!	
	END SUBROUTINE avg_walls
!*******************************************************************************
  subroutine hybrid_calc
!
! This averages over tpts
  integer :: i,j,k
  do k=1,4; do j=1,4; do i=-1,ny+2
    stat_avg(i,j,k)=SUM(stat1d_sa(i,j,k,1:tpts))/real(tpts,kind=mytype)
  end do; end do; end do
  do k=1,4; do j=1,4; do i=-1,ny+2
    flucs(i,j,k)=SUM(flucs1d_sa(i,j,k,1:tpts))/real(tpts,kind=mytype)
  end do; end do; end do

  if (iquad == 1)  then
    do j=1,9; do i=-1,ny+2
      quad_avg(i,j)=SUM(quad1d_sa(i,j,1:tpts))/real(tpts,kind=mytype)
    end do; end do
  end if

  if (ivort == 1)  then
    do j=1,6; do i=-1,ny+2
      omega_avg(i,j)=SUM(omega1d_sa(i,j,1:tpts))/real(tpts,kind=mytype)
    end do; end do
    do k=1,8; do j=1,2; do i=-1,ny+2
      wx_avg(i,j,k)=SUM(wx1d_sa(i,j,k,1:tpts))/real(tpts,kind=mytype)
    end do; end do; end do
  end if

  if (itke == 1) then
    do k=1,7; do j=1,4; do i=-1,ny+2
      tke_avg(i,j,k)=SUM(tke1d_sa(i,j,k,1:tpts))/real(tpts,kind=mytype)
    end do; end do; end do
    do k=1,6; do j=1,4; do i=-1,ny+2
      kinen(i,j,k)=SUM(kinen1d_sa(i,j,k,1:tpts))/real(tpts,kind=mytype)
    end do; end do; end do
  end if

  if (ishear == 1) then
    do j=1,4; do i=-1,ny+2
      shear(i,j)=SUM(shear1d_sa(i,j,1:tpts))/real(tpts,kind=mytype)
    end do; end do
  end if

  if (iangle == 1) then
    do j=1,5; do i=-1,ny+2
      angle(i,j)=SUM(angle1d_sa(i,j,1:tpts))/real(tpts,kind=mytype)
    end do; end do
  end if

  if (ianis == 1) then
    do j=1,4; do i=-1,ny+2
      anisb(i,j)=SUM(anisb1d_sa(i,j,1:tpts))/real(tpts,kind=mytype)
    end do; end do
    do i=-1,ny+2
      anisII(i)=SUM(anisII1d_sa(i,1:tpts))/real(tpts,kind=mytype)
    end do
    do i=-1,ny+2
      anisIII(i)=SUM(anisIII1d_sa(i,1:tpts))/real(tpts,kind=mytype)
    end do
    do i=-1,ny+2
      anisf(i)=SUM(anisf1d_sa(i,1:tpts))/real(tpts,kind=mytype)
    end do
    do i=-1,ny+2
      anisg(i)=SUM(anisg1d_sa(i,1:tpts))/real(tpts,kind=mytype)
    end do
  end if

  isteady=1

  return
end subroutine hybrid_calc
!*******************************************************************************

	END MODULE calc_post
