      subroutine smrt_read_dhru2hru

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine reads in the file containing the information to convert 
!!    SWAT-MODFLOW-based disaggregated HRU (dhru) variables to SWAT-based normal/aggregated HRU variables
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    nhru         |none          |the current number of SWAT hrus
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    h2d_map      |none          |Array that stores the IDs of the DHRUs that make 
!!                                |up each HRU.
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name         |units        |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i            |none         |counter index for reading in contributing
!!                 |             |areas and looping through all grids
!!    j            |none         |counter index for the number of SWAT-MODFLOW dhru
!!    hruID        |none         |index of current SWAT-MODFLOW dhru
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use smrt_parm
      implicit none
      
!     Initialize local variables
      integer num_hru,i,j,hruID,num_hru_dhrus,subID


      !open the input file containing the HRU --> DHRU mapping data
      open (6003,file="swatmf_dhru2hru.txt")
      print *, 'Reading DHRU to HRU mapping...'

      !The first line of this file is the total number of HRUs in the watershed (ie, how many will be read in)
      read(6003,*) num_hru !should equal SWAT's nhru
      
      !allocate the main mapping array (grid --> DHRU)
      allocate(h2d_map(dhru))
      

      !loop throug the HRUs
      do i=1,num_hru
        
        ! the HRU's global ID within the watershed, the number of dhrus contributing to this HRU, 
        ! and the subbasin number for this HRU.
        read(6003,*) hruID, num_hru_dhrus, subID
        
        if(num_hru_dhrus.gt.0)then
          allocate(h2d_map(i)%dhru_id(num_hru_dhrus))
          allocate(h2d_map(i)%dhru_perc(num_hru_dhrus))
          read(6003,*) (h2d_map(i)%dhru_id(j),j=1,num_hru_dhrus) ! list of dhru ID numbers which contribute to this hru
          read(6003,*) (h2d_map(i)%dhru_perc(j),j=1,num_hru_dhrus) ! list of % areas of that dhru contributing to this hru
        else
          read(6003,*)
          read(6003,*)
        endif
        
      enddo

      close(6003)
      
      !write out information to log file
      write(6008,*) 'swatmf_init: swatmf_dhru2hru.txt has been read'

      return
      end