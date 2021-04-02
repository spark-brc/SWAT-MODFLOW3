      subroutine smrt_grid2dhru4D (rt4dVar,smrtVar,location,dim)

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts 4D RT3D-based grids where the index of the 4th dimension being 
!!    modified (number of chemical species) is defined by "dim" to SWAT-MODFLOW-based disaggregated hrus (dhru)
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dim         |none          |the location in the 4th dimension of the RT3D grid
!!                |              |in which to retrieve the data, see below
!!    rt4dVar     |unknown       |4D RT3D variable (NCOL,NROW,location,dim)
!!    smrtVar     |unknown       |SWAT-MODFLOW variable (list of dhru) to be populated
!!                |              |with the contents of the MODFLOW variable
!!    location    |none          |an array(NCOL, NROW) with the value equal to which
!!                |              |layer(NLAY) the MODFLOW value is to be taken from
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
!!    lay         |none          |layer index for MODFLOW
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use GLOBAL, only: NCOL,NROW,NLAY !MODFLOW
      use rt_global
      use smrt_parm
      implicit none
      
!     Define local variables
!     real sum
      double precision rt4dVar(NCOL,NROW,NLAY,ncomp)
      real smrtVar(dhru)
      real location(NCOL,NROW)
      integer i, j, row, col, lay, dim
      real perc_area
      
!     Initialize variables
      smrtVar = 0.
      
!     Convert MODFLOW-grid variables into SWAT-MODFLOW-dhru variables by multiplying each
!     contributing grid variable by its percent area contributing to each dhru
!     For 3D variable arrays taking specified layer value only
      do i=1,dhru
        do j=1,size(g2d_map(i)%cell_row)
          row = g2d_map(i)%cell_row(j)
          col = g2d_map(i)%cell_col(j)
          perc_area = g2d_map(i)%cell_perc(j)
          if(col.ne.0 .and. row.ne.0)then
            lay = location(col,row)
            smrtVar(i) = smrtVar(i) + 
     &          rt4dVar(col,row,lay,dim)*perc_area
          endif
        enddo
      enddo
      
      return
      end