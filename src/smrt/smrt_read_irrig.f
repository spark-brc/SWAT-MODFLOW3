      subroutine smrt_read_irrig()

            !This subroutine reads the file containing the information that links MODFLOW Well package cells (those
      !that provide irrigation water) to SWAT sub-basins

!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name            |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mf_sub_irrigate |none          |array that holds the list of MODFLOW cells that
                                     !provide irrigation to SWAT.
                                     ! 1 = SWAT sub-basin in which cell is located
                                     ! 2 = Row of MODFLOW cell
                                     ! 3 = Col of MODFLOW cell
                                     ! 4 = SWAT HRU in which cell is located
!!    ncell_irrigate  |              |number of MODFLOW cells that provide irrigation 
                                     !to SWAT
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


      !import variables
      use parm
      use smrt_parm !smrt linkage
      use GLOBAL,only:NROW,NCOL,NLAY !MODFLOW

      implicit none
      
      !initialize local variables
      integer m,n,sub_id,cell_row,cell_col,cell_lay,hru_id,nhru_irrigate
      
      !write task to screen
      print *, 'Reading MODFLOW Cells that provide Irrigation Water...'

      !open the file and read in number of MODFLOW cells that provide irrigation to SWAT
      open (6006,file="swatmf_irrigate.txt")
      read(6006,*)
      read(6006,*) ncell_irrigate
      
      !allocate the main arrays
      allocate(mf_sub_irrigate(ncell_irrigate,5))
      allocate(mf_sub_irrig_depth(msub))
      allocate(mf_sub_irrig_flag(msub))
      allocate(mf_convey_ratio(msub))
      allocate(mf_runoff_ratio(msub))
      mf_sub_irrigate = 0.
      mf_sub_irrig_depth = 0.
      mf_sub_irrig_flag = 0
      mf_convey_ratio = 0.
      mf_runoff_ratio = 0.
      
      !allocate the size of the applied irrigation array, if SWAT HRUs are providing irrigation water
      !to MODFLOW pumping cells
      if(mf_irrigation_swat.eq.1) then
        allocate(irrig_volume(nhru))
        allocate(irrig_depth(nhru))
        irrig_volume = 0.
        irrig_depth = 0.
      endif

      !read conveyance efficiency and surface runoff ratio for each sub-basin
      !only read if irrigation depends on specified MODFLOW pumping (in WELL package)
      read(6006,*)
      if(mf_irrigation_swat.eq.0) then
        do m=1,msub
          read(6006,*) mf_convey_ratio(m),mf_runoff_ratio(m)
        enddo
        read(6006,*)
      endif

      !loop through cells and read in information
      do n=1,ncell_irrigate
        
        !read sub-basin #, cell Row, cell Column, and HRU ID (not used)
        if(mf_irrigation_swat.eq.0) then
          read(6006,*) sub_id,cell_row,cell_col,hru_id
        else
          read(6006,*) sub_id,cell_row,cell_col,cell_lay,hru_id
        endif

        !store these values in the main array
        mf_sub_irrigate(n,1) = sub_id
        mf_sub_irrigate(n,2) = cell_row
        mf_sub_irrigate(n,3) = cell_col
        mf_sub_irrigate(n,4) = hru_id
        if(mf_irrigation_swat) mf_sub_irrigate(n,5) = cell_lay

      enddo

      !read the HRUs that are to receive irrigation water (no other HRUs will receive
      !pumped irrigation water from MODFLOW)
      if(mf_irrigation_swat.eq.0) then
      allocate(mf_hru_irrig_flag(nhru))
      allocate(mf_nsub_hru_irrig(msub))
      allocate(mf_sub_hru_irrig(msub,nhru))
      mf_hru_irrig_flag = 0
      mf_nsub_hru_irrig = 0
      mf_sub_hru_irrig = 0
      read(6006,*)
      read(6006,*) nhru_irrigate
      do n=1,nhru_irrigate
        read(6006,*) sub_id,hru_id
        mf_hru_irrig_flag(hru_id) = 1
        mf_nsub_hru_irrig(sub_id) = mf_nsub_hru_irrig(sub_id) + 1
        mf_sub_hru_irrig(sub_id,mf_nsub_hru_irrig(sub_id)) = hru_id
      enddo

      endif

      !open up files for pumped irrigation output (filled up in 'smrt_conversion2swat.f')
      if(mf_irrigation_swat.eq.0) then
        open(300009,file='swatmf_out_pumped_sub')
        write(300009,*) 'MODFLOW Pumped, SWAT Received'
        open(300010,file='swatmf_out_pumped_hru')
        write(300010,*) 'MODFLOW pumped groundwater applied to HRUs'
      else
        !open up files for pumping rates specified by SWAT HRUs (auto-irrigation)
        open(300011,file='swatmf_out_pumping')
        write(300011,*) 'MODFLOW Pump Rates specified by AutoIrrigation'
        write(300011,*) 'Units: according to ITMUNI and LENUNI in .dis'
      endif

      !close the file
      close(6006)

      !write out information to log file
      write(6008,*) 'swatmf_init: swatmf_irrigate.txt has been read'

      return
      end