      subroutine smrt_read_dhru2grid

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine reads in the file containing the information to convert 
!!    SWAT-MODFLOW-based disaggregated hrus (dhru) variables to MODFLOW-based grid variables
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    NCOL         |none          |the current number of columns of MODFLOW grids
!!    NROW         |none          |the current number of rowss of MODFLOW grids
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    d2g_map      |none          |Array that stores the attributes (ID, percent area)
!!                                |of the DHRUs that intersect each MODFLOW grid cell.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i            |none          |counter index for the number of MOFLOW grids
!!    j            |none          |counter index for reading in contributing
!!                 |              |areas and looping through all hrus
!!    temp         |none          |index of current MODFLOW grid
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use GLOBAL, only: NCOL,NROW,NLAY !MODFLOW
      use smrt_parm
      implicit none
      
      !Initialize local variables
      integer i,j,cellID,num_cell_dhrus,ncells
      

      !The first line of this file is the total number of MODFLOW grids (active and inactive) = NROW * NCOL
      open (6002,file="swatmf_dhru2grid.txt")
      print *, 'Reading DHRU to Grid mapping...'
      
      !The first line of this file is the total number of MODFLOW grid cells (ie, how many will be read in)
      read(6002,*) ncells
      
      !allocate the main mapping array (grid --> DHRU)
      allocate(d2g_map(ncells))
      
      !loop through the MODFLOW grid cells, storing attributes of the DHRUs that intersect each cell
      do i=1,ncells

        !MODFLOW cell ID, # of DHRUs contributing to the cell
        read(6002,*) cellID,num_cell_dhrus

        if(num_cell_dhrus.gt.0) then
          !allocate the second dimension of the array, according to the number of contributing cells
          allocate(d2g_map(i)%dhru_id(num_cell_dhrus))
          allocate(d2g_map(i)%dhru_perc(num_cell_dhrus))
          read(6002,*) (d2g_map(i)%dhru_id(j),j=1,num_cell_dhrus) ! list of dhru ID numbers which contribute to this grid cell
          read(6002,*) (d2g_map(i)%dhru_perc(j),j=1,num_cell_dhrus) ! list of % areas of that dhru contributing to this grid cell
        endif
        
      enddo

      close(6002)
      
      !write out information to log file
      write(6008,*) 'swatmf_init: swatmf_dhru2grid.txt has been read'

      return
      end
