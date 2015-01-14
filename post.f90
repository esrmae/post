program post
  use shared_post
  use read_post
  use cal_coefficient
  use calc_post3d
  use calc_struct
  use calc_post
  use calc_super
  use write_post

  implicit none

!*******************************************************************************
! 	The main program for the post processing code.
!*******************************************************************************

  call read_setup
  write (*,*) 'DONE READ SETUP'
  call read_vals
  write (*,*) 'DONE READING'
  call calc_perform
  write (*,*) 'DONE CALCS'
  if (read_3d == 1) then 
   CALL calc3d
   write (*,*) 'DONE CALCS_3D'
  end if
  if (istruct ==1) then
    call struct_av
    write(*,*) 'DONE STRUCT_AVERAGE'
  end if
  if (isuper ==1) then
    call struct_super
    write(*,*) 'DONE CALC_SUPER'
  end if
  simit=sim_print
  call write_vals
  write (*,*) 'DONE WRITING'

  write (*,*) 'DONE FLOW VIS'	

!*******************************************************************************

  stop
end program post
