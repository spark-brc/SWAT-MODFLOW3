*     RT3D - Reactive Transport in 3 Dimensions
*     This is a modified version from Clement (1997) that includes dynamic allocation of arrays

*     Author: Ryan Bailey (Colorado State University)
*     Date: 10-24-2012

*     read UZF-RT3D data for the groundwater reactive transport simulation
      subroutine rt_read

      use rt_global

      implicit none
      integer   N,L
      common   /GCGIDX/L(19) !for implicit shceme     

*     initialize the UZF-RT3D simulation: open files, allocate arrays
      call rt_init

      return 
      end
