      subroutine smrt_dhru2grid2D (smrtVar, mfVar2, mult_TF)

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts SWAT-MODFLOW-based disaggregated HRUs (dhrus) to 2D MODFLOW-based grids.
!!    Additionally, this subroutine multiplies a MODFLOW 2D variable by the grid area of each
!!    MODFLOW grid if mult_TF is 1, divides a MODFLOW 2D variable by the grid area of each
!!    MODFLOW grid if mult_TF is 2, or does nothing additional if mult_TF is 0.
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    smrtVar     |unknown       |smrt variable (list per dhru)
!!    mfVar2      |unknown       |2D MODFLOW variable (NCOL, NROW) to be
!!                |              |populated with the contents of the SWAT-MODFLOW variable
!!    mult_TF     |n/a           |integer, if "1" will multiply array
!!                |              |by the grid area, if "2" will divide
!!                |              |array by the grid area, if "0" will not
!!                |              |convert the variable based on MODFLOW area
!!    DELR        |LENUNI        |MODFLOW variable, the thickness of the row per column
!!    DELC        |LENUNI        |MODFLOW variable, the thickness of the column per row
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mfVar2      |unknown       |2D MODFLOW variable (NCOL, NROW) now populated
!!                |              |with the contents of the SWAT-MODFLOW variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    newVar2     |unknown       |2D MODFLOW variable (NCOL, NROW) that is
!!                |              |temporarily populated with the contents of 
!!                |              |the SM variable before checking for remaining
!!                |              |grid area not contributed to by a dhru
!!    areaUsed    |%             |2D MODFLOW variable (NCOL, NROW) that is
!!                |              |used to track how much of the grid area has
!!                |              |been contributed to by the dhrus so that any
!!                |              |remaining grid area not contributed to by a dhru
!!                |              |uses the original MODFLOW value   
!!    i           |n/a           |MODFLOW row loop index
!!    j           |n/a           |MODFLOW column loop index
!!    ctr         |n/a           |MODFLOW grid id index
!!    k           |n/a           |d2g_size loop index
!!    cellUsed    |n/a           |a true/false for whether the MODFLOW grid
!!                |              |has interacted with a SWAT-MODFLOW dhru and should
!!                |              |therefore be converted for area/unit reasons
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use GLOBAL, only: NCOL,NROW,DELR,DELC !MODFLOW
      use smrt_parm
      implicit none
      
!     Define local variables
      real smrtVar(dhru)
      real mfVar2(NCOL,NROW)
      real newVar2(NCOL,NROW)
      real areaRemain(NCOL,NROW)
      integer mult_TF, i, j, k, ctr, dhruID
      logical cellUsed
      real perc_area
      
!     Initialize local variables
      newVar2 = 0.
      areaRemain = 1.
      cellUsed = .false.
      
!     Convert SWAT-MODFLOW dhru variables into MODFLOW-grid variables by multiplying each
!     contributing dhru variable by its percent area contributing to each grid
      
      !loop through the MODFLOW grid cells
      ctr = 1
      do j=1,NROW
        do i=1,NCOL

          !loop through the DHRUs that contribute area to the MODFLOW grid cell
          do k=1,size(d2g_map(ctr)%dhru_id)
            
            !ID and percent area of the current DHRU
            dhruID = d2g_map(ctr)%dhru_id(k)
            perc_area = d2g_map(ctr)%dhru_perc(k)

            !Convert the DHRU information into grid information
            newVar2(i,j) = newVar2(i,j)+ (smrtVar(dhruID)*perc_area)
            
            !Track how much of the area of the grid has been contributed to (%)
            areaRemain(i,j) = areaRemain(i,j) - perc_area
            cellUsed = .true.

          enddo
          
          !If the grid cell was contributed to by dhrus, and
          !the units need to convert the area, do so here
          if(cellUsed .and. mult_TF.eq.1) then
            !Multiply variable by cell area
            newVar2(i,j) = newVar2(i,j)*DELR(i)*DELC(j)
          else if(cellUsed .and. mult_TF.eq.2) then
            !Divide variable by cell area
            newVar2(i,j) = newVar2(i,j)/(DELR(i)*DELC(j))
          endif
          
          !Store the converted dhru results into the MODFLOW variable
          if(areaRemain(i,j) > 0) then
            !If the dhrus do not completely cover the grid cell, retain a portion of the original grid data
            mfVar2(i,j) = newVar2(i,j) + (mfVar2(i,j)*areaRemain(i,j))
          else
            mfVar2(i,j) = newVar2(i,j)
          endif
          
          !reset the counters for the loop
          cellUsed = .false.
          ctr = ctr + 1
        enddo
      enddo

      return
      end