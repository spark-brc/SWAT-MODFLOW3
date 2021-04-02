      subroutine smrt_dhru2rtgrid3D (smrtVar, rt3dVar)

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts SWAT-MODFLOW-based disaggregated hrus (dhru) to 3D RT3D-based grids
!!    where the index of the 3rd dimension being modified (number of chemical species) is defined by "dim"
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dhru2grid   |none          |Array of percent areas of SWAT-MODFLOW-dhrus
!!                |              |contributing to each MODFLOW-grid per watershed
!!                |              |format: (MODFLOW grid ID, hru)
!!    smrtVar     |unknown       |SWAT-MODFLOW variable (dhru)
!!    dim         |none          |the location in the 4th dimension of the RT3D grid
!!                |              |in which to store the converted data from smrtVar
!!    rt3dVar     |unkwon        |3D RT3D variable (NCOL,NROW,ncomp) to be
!!                |              |populated with the contents of the SWAT-MODFLOW variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    rt3dVar     |unkwon        |3D RT3D variable (NCOL,NROW,ncomp) to be
!!                |              |populated with the contents of the SWAT-MODFLOW variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |MODFLOW row loop index
!!    j           |none          |MODFLOW column loop index
!!    l           |none          |SWAT-MODFLOW dhru loop index
!!    ctr         |none          |index for what total grid ID should be used
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use GLOBAL, only: NCOL,NROW,DELR,DELC !MODFLOW
      use rt_global
      use smrt_parm
      implicit none
      
!     Define local variables
      real smrtVar(dhru)
      real rt3dVar(NCOL,NROW)
      integer mult_TF, i, j, k, l, ctr, dim, dhruID, dum
      logical cellUsed
      real perc_area
      
!     Initialize variables
      cellUsed = .false.
      rt3dVar = 0.

!     Convert SWAT-MODFLOW-dhru variables into RT3D grid variables by multiplying each
!     contributing dhru variable by its percent area contributing to each grid
!     Note: data is only placed in the 3rd dimension of the array specified by the variable "dim"
!     contributing dhru variable by its percent area contributing to each grid
      ctr = 1
      do j=1,NROW
        do i=1,NCOL
          
          !loop through the dhrus that contribute to the grid cell
          do k=1,size(d2g_map(ctr)%dhru_id)
            dhruID = d2g_map(ctr)%dhru_id(k)
            perc_area = d2g_map(ctr)%dhru_perc(k)
            
            if(cellUsed)then
              !if this is the second or more time referencing this cell, add to its contents
              rt3dVar(i,j) = rt3dVar(i,j)+smrtVar(dhruID)*
     &                             perc_area
            else
              !if this is the first time referencing this cell, overwrite its contents
              rt3dVar(i,j) = smrtVar(dhruID)*perc_area
              cellUsed = .true.
            endif

          enddo
        
          !reset the counters for the loop
          cellUsed = .false.
          ctr = ctr + 1

        enddo !go to next cell
      enddo

      return
      end