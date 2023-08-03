      subroutine smrt_rt3d_run

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine calls conversions from SWAT to RT3D, then runs RT3D
!!    for one day, then the groundwater results are converted back from RT3D
!!    to SWAT in smrt_conversion2swat
    

      !Convert SWAT variables units into RT3D variables and units
      call smrt_conversion2rt3d
          
      !Prep. RT3D
      call rtmf_prepare !rtb
          
      !Call RT3D
      call rt_run !rtb
          
      !There is no call to conversion back to SWAT here because this is handled
      !in smrt_converstion2swat which is called after this from smrt_mf_run       !tcw

      end subroutine smrt_rt3d_run
