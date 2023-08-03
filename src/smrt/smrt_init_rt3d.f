      subroutine smrt_init_rt3d

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine initializes RT3D linking subroutins of
!!    SWAT-MODFLOW-RT3D

!     Set up RT3D data and allocate arrays
      print *, 'RT3D is being used'
      call rt_read
      
      !write out information to log file
      write(6008,*) 'swatmf_init: RT3D has been initialized'
      
      end subroutine smrt_init_rt3d
