      subroutine smrt_rt3d_stress
!!    Ryan Bailey
!!    Colorado State University 2012

!!
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine calls RT3D's stress period reader
!!      
        use mf_rt_link, only: rt_kstp !MODFLOW-RT3D linkage
        implicit none
          
        rt_kstp = 0
        call rt_stress

      end subroutine smrt_rt3d_stress
