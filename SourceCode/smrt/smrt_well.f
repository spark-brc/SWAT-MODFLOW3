      subroutine smrt_well

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the SWAT applied irrigation volumes to MODLFOW
!!    pumping rates in the WELL package.

        !Import variables
        use parm, only: day_total,nhru
        use GLOBAL, only: ITMUNI,LENUNI
        use GWFWELMODULE, ONLY:NWELLS,WELL
        use smrt_parm
        implicit none
        
        !Define local variables
        integer i,j,n,hru_id,cell_row,cell_col,cell_lay,num_mf_wells
        real    volume,pump_rate
        integer, dimension (:), allocatable :: hru_wells
        
        write(300011,*)
        write(300011,*) 'Day:',day_total
        write(300011,*) 'daily pumping rates provided to MODFLOW'



        !determine the number of WELL cells in each HRU ---------------------------------------------------------------
        allocate(hru_wells(nhru))
        hru_wells = 0
        do i=1,ncell_irrigate
          hru_id = mf_sub_irrigate(i,4)
          hru_wells(hru_id) = hru_wells(hru_id) + 1
        enddo
        

        !loop through the MODFLOW WELL cells that are to receive irrigation water from SWAT; --------------------------
        !calculate pumping rate 
        num_mf_wells = NWELLS !copy in the constant value from MODFLOW's files
        do i=1,ncell_irrigate
          
          !increment the number of cells in the MODFLOW WELL array
          num_mf_wells = num_mf_wells + 1

          !row and column of the WELL cell
          cell_row = mf_sub_irrigate(i,2)
          cell_col = mf_sub_irrigate(i,3)
          cell_lay = mf_sub_irrigate(i,5)

          !total volume --> pumping rate
          hru_id = mf_sub_irrigate(i,4)
          volume = irrig_volume(hru_id) !m3/day
          pump_rate = volume / hru_wells(hru_id) !depends on the number of wells in the HRU

          !convert to MODFLOW units - time
          if(ITMUNI.eq.1) then !seconds
            pump_rate = pump_rate / 86400.
          elseif(ITMUNI.eq.2) then !minutes
            pump_rate = pump_rate / 1440.
          elseif(ITMUNI.eq.3) then !hours
            pump_rate = pump_rate / 24.
          elseif(ITMUNI.eq.5) then !years
            pump_rate = pump_rate * 365.
          endif

          !convert to MODFLOW units - length
          if(LENUNI.eq.1) then !feet
            pump_rate = pump_rate * 35.287552
          elseif(LENUNI.eq.3) then !cm
            pump_rate = pump_rate * 1000000.
          endif

          !add to WELL array
          WELL(1,num_mf_wells) = cell_lay
          WELL(2,num_mf_wells) = cell_row
          WELL(3,num_mf_wells) = cell_col
          WELL(4,num_mf_wells) = pump_rate * (-1) !make it extraction
          
          !print out to file
          write(300011,*) (WELL(j,num_mf_wells),j=1,4)

        enddo

      end subroutine smrt_well
