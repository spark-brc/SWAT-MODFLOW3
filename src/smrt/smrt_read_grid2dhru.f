      subroutine smrt_read_grid2dhru

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine reads in the file containing the information to convert 
!!    MODFLOW-based grid variables to SWAT-MODFLOW-based disaggregated HRU (dhru) variables
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    dhru         |none          |the total number of disaggregated HRUs in the 
!!                 |              |entire watershed (not just the current subbasin)
!!    NCOL         |none          |the current number of columns of MODFLOW grids
!!    NROW         |none          |the current number of rowss of MODFLOW grids
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    g2d_map      |none          |Array that stores the attributes (row, col, percent area)
!!                                |of the MODFLOW grid cells that intersect each DHRU.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name         |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i            |none         |counter index for the number of SWAT-MODFLOW dhru
!!    j            |none         |counter index for reading in contributing
!!                 |             |areas and looping through all grids
!!    dhruID       |none         |index of current SWAT-MODFLOW dhru
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use GLOBAL, only: NCOL, NROW, NLAY !MODFLOW
      use smrt_parm !SWAT-MODFLOW linkage
      implicit none
      
      !Initialize local variables
      integer i,j,dhruID,num_dhru_cells


      !Read in the ID and percent area of each MODFLOW grid contributing to each SWAT-MODFLOW dhru
      open (6004,file="swatmf_grid2dhru.txt")
      print *, 'Reading Grid to DHRU mapping...'

      !The first line of this file is the total number of DHRUs in the watershed (ie, how many will be read in)
      read(6004,*) dhru

      !allocate the main mapping array (grid --> DHRU)
      allocate(g2d_map(dhru))


      !loop through each DHRU, reading in the information for each intersected grid cell
      do i=1,dhru
        
        !DHRU ID, # of cells contributing to the DHRU
        read(6004,*) dhruID,num_dhru_cells
        
        if(num_dhru_cells.gt.0) then
        
        !allocate the second dimension of the array, according to the number of contributing cells
        allocate(g2d_map(i)%cell_row(num_dhru_cells))
        allocate(g2d_map(i)%cell_col(num_dhru_cells))
        allocate(g2d_map(i)%cell_perc(num_dhru_cells))

        !now, fill in the second dimension of each array
        read(6004,*) (g2d_map(i)%cell_row(j),j=1,num_dhru_cells)
        read(6004,*) (g2d_map(i)%cell_col(j),j=1,num_dhru_cells)
        read(6004,*) (g2d_map(i)%cell_perc(j),j=1,num_dhru_cells)

        else
          read(6004,*)
          read(6004,*)
          read(6004,*)
        endif

      enddo
      close(6004)

      !write out information to log file
      write(6008,*) 'swatmf_init: swatmf_grid2dhru.txt has been read'

      return
      end