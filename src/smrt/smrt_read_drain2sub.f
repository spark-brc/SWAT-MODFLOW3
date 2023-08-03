      subroutine smrt_read_drain2sub() !rtb drain

      !Ryan Bailey (CSU) 2017 08-09

      !This subroutine reads in the file containing the information to link
      !MODFLOW Drain cells to SWAT sub-basins. This is only used if the MODFLOW
      !Drain package is turned on.

!     Import variables
      use GLOBAL, only:LENUNI,ITMUNI,NROW,NCOL !MODFLOW
      use smrt_parm !smrt linkage
      implicit none
      
!     Initialize local variables
      integer i,drn_row,drn_col,sub_basin     
      
!     Read in row, column, and associated sub-basin for each DRAIN cell
      open (6007,file="swatmf_drain2sub.txt")
      print *, 'Reading Drain to Subbasin Mapping...'
      print *
      read(6007,*) ndrn_subs ! # of MODLFOW drain cells in the SWAT domain
      read(6007,*)

      !allocate and populate the DRAIN_SUBS array
      allocate(drn_subs(NROW,NCOL))
      drn_subs = 0
      do i=1,ndrn_subs
        read(6007,*) drn_row,drn_col,sub_basin
        drn_subs(drn_row,drn_col) = sub_basin
      enddo
      close(6007)

      !write out information to log file
      write(6008,*) 'swatmf_init: swatmf_drain2sub.txt has been read'

      return
      end