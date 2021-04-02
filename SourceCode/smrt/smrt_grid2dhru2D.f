      subroutine smrt_grid2dhru2D (mfVar2,smrtVar)

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts 2D MODFLOW-based grids to SWAT-MODFLOW-based disaggregated HRUs (dhrus)
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mfVar2      |unknown       |2D MODFLOW variable (NCOL, NROW)
!!    smrtVar     |unknown       |SWAT-MODFLOW variable (list of dhrus) to be populated
!!                |              |with the contents of the MODFLOW variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    smrtVar     |unknown       |SWAT-MODFLOW variable (list of dhrus) now populated
!!                |              |with the contents of the MODFLOW variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |SWAT-MODFLOW dhru loop index
!!    j           |none          |g2d_size loop index
!!    row         |none          |row index for MODFLOW
!!    col         |none          |column index for MODFLOW
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use GLOBAL, only: NCOL,NROW !MODFLOW
      use smrt_parm
      implicit none
      
!     Define local variables
      real mfVar2(NCOL,NROW)
      real smrtVar(dhru)
      integer i, j, row, col
      real perc_area
      
!     Initialize variables
      smrtVar = 0.
      
!     Convert MODFLOW-grid variables into SWAT-MODFLOW-dhru variables by multiplying each
!     contributing grid variable by its percent area contributing to each dhru
      
      !loop through all of the DHRUs
      do i=1,dhru
        !loop through all of the cells that contribute to the DHRU
        do j=1,size(g2d_map(i)%cell_row)
          row = g2d_map(i)%cell_row(j)
          col = g2d_map(i)%cell_col(j)
          perc_area = g2d_map(i)%cell_perc(j)
          if(col.ne.0 .and. row.ne.0)then
            smrtVar(i) = smrtVar(i) + mfVar2(col,row)*perc_area
          endif
        enddo
      enddo
      
      return
      end