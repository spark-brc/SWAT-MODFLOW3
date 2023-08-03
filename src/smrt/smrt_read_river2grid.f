      subroutine smrt_read_river2grid()

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine reads in the file containing the information to calculate
!!    river conductance for MODFLOW's river cells
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    subtot        |none          |number of subbasins in watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    grid2riv_id   |none          |a list per MODFLOW river cells containing 
!!                  |              |the SWAT river IDs within that river cell
!!    grid2riv_len  |m             |a list per MODFLOW river cells containing 
!!                  |              |the SWAT river lengths within that river cell
!!    riv_nsubs     |none          |a list per MODFLOW river cells containing 
!!                  |              |the number of subbasins that have spatial areas within
!!                  |              |that river cell
!!    river_cells   |LENUNI        |an array of the properties of the river grid cells,
!!                  |              |are the thicknesses of the river bed of the grid cells
!!    nriver_cells  |n/a           |# of river cells that are linked to SWAT sub-basins
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i             |none          |counter index for the number of SWAT subbasins
!!    j             |none          |counter index for reading contributing
!!                  |              |MODFLOW rive grids
!!    temp          |none          |index of current SWAT subbasin
!!    ngrid_current |none          |index of the number of contributing MODFLOW
!!                  |              |river grids
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use GLOBAL, only:LENUNI,ITMUNI !MODFLOW
      use smrt_parm !smrt linkage
      implicit none
      
!     Initialize local variables
      integer i, j, temp, nsub_current      
      
!     Read in the id and percent area of each SWAT HRU contributing to each
!     MODFLOW grid
      open (6005,file="swatmf_river2grid.txt")
      print *, 'Reading Subbasin to Grid mapping...'
      print *

      !# of river cells that are linked to SWAT sub-basins
      read(6005,*) nrivcells_subs
      
!     Allocate the variable sizes, which needs to be called here before reading the file continues
      call smrt_allocate()
      
!     Read which SWAT river reaches are within each MODFLOW river cell, and their associated river lengths
      do i=1,nrivcells_subs
        read(6005,*) temp,grid2riv_gridID(i),nsub_current ! cell #, then the number of subbasins associated with the cell
        if (nsub_current.gt.0) then
          read(6005,*) (grid2riv_id(i,j), j=1,nsub_current) ! list of subbasin ID numbers which are within the current cell
          read(6005,*) (grid2riv_len(i,j),j=1,nsub_current) ! list of length of river within the current cell
          riv_nsubs(i) = nsub_current
        else
          read(6005,*)
          read(6005,*)
        endif
        
      enddo
      close(6005)

      !write out information to log file
      write(6008,*) 'swatmf_init: swatmf_river2grid.txt has been read'

      return
      end