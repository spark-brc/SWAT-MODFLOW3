      program main

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine links sets up the SWAT-MODFLOW-RT3D linking subroutines
!!    and "events" and then calls SWAT-MODFLOW through "events"
!!
        use mf_rt_link, only: rt_active !MODFLOW-RT3D Linkage
        use smrt_parm, only: mf_active !SWAT-MODFLOW linkage
        implicit none
        
        ! Read the SWAT-MODFLOW input files
        call smrt_read_link
 
        ! Run SWAT's main subroutine
        call swat_main

        ! Close swat-modflow (close files, deallocate variables, etc.)s
        call smrt_close
        call mf_close
        call rt_close
      end program main
