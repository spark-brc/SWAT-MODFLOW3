      subroutine smrt_dhru2hru (smrtVar, swatVar, mult_TF)

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts SWAT-MODFLOW-based disaggregated HRUs (dhru) to SWAT HRUs.
!!    Additionally, this subroutine multiplies a SWAT variable by the area, in
!!    km**2, of each SWAT hru if mult_TF is true 1, divides a SWAT variable by
!!    the area, in km**2, of each SWAT hru if mult_TF is false 2, or does nothing
!!    additional if mult_TF is 0.
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    smrtVar     |unknown       |SWAT-MODFLOW variable (list of dhrus)
!!    swatVar     |unknown       |SWAT variable (list of nhru) to be populated
!!                |              |with the contents of the SWAT-MODFLOW variable
!!    hru_km(:)   |km**2         |area of HRU in square kilometers
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    swatVar     |unknown       |SWAT variable (list of nhru) to be populated
!!                |              |with the contents of the SWAT-MODFLOW variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |SWAT HRU loop index
!!    j           |none          |SWAT-MODFLOW dhru loop index
!!    cellUsed    |n/a           |a true/false for whether the SWAT hru
!!                |              |has interacted with a SWAT-MODFLOW dhru and should
!!                |              |therefore be converted for area/unit reasons
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use parm, only: nhru, hru_km !SWAT
      use smrt_parm
      implicit none
      
!     Define local variables
      real smrtVar(dhru)
      real swatVar(nhru)
      integer mult_TF, i, j, dhruID
      logical cellUsed
      real perc_area
      
!     Initialize local variables
      cellUsed = .false.
      
!     Convert SWAT-MODFLOW-dhru variables into SWAT-HRU variables by multiplying each
!     contributing dhru variable by its percent area contributing to each HRU
      do i=1,nhru
        
        !loop through the DHRUs that are within the current HRU
        do j=1,size(h2d_map(i)%dhru_id)
          
          !ID and percent contributing area for the current DHRU
          dhruID = h2d_map(i)%dhru_id(j)
          perc_area = h2d_map(i)%dhru_perc(j)

          if(cellUsed)then
            !if this is the second or more time referencing this cell, add to its contents
            swatVar(i) = swatVar(i) + (smrtVar(dhruID)*perc_area)
          else
            !if this is the first time referencing this cell, overwrite its contents
            swatVar(i) = smrtVar(dhruID)*perc_area
            cellUsed = .true.
          endif

        enddo
        
        !If the hru was contributed to by dhrus, and !the units need to convert the area, do so here
        if(cellUsed .and. mult_TF.eq.1)then
          !Multiply variable by cell area, in km
          swatVar(i)=swatVar(i)*hru_km(i)
        else if(cellUsed .and. mult_TF.eq.2)then
          !Divide variable by cell area, in km
          swatVar(i)=swatVar(i)/(hru_km(i))
        endif
        
        !reset the counters for the loop
        cellUsed = .false.
      enddo
      
      return
      end