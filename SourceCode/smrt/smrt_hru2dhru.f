      subroutine smrt_hru2dhru (swatVar, smrtVar)

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts a SWAT HRU-variable to a SWAT-MODFLOW disaggregated HRUs (dhrus) variable
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    swatVar     |unknown       |SWAT variable (list of nhru)
!!    smrtVar     |unknown       |SWAT-MODFLOW variable (list of dhrus) to be populated
!!                |              |with the contents of the SWAT variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    smrtVar     |unknown       |SWAT-MODFLOW variable (list of dhrus) to be populated
!!                |              |with the contents of the SWAT variable
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    i           |none          |SWAT HRU loop index
!!    j           |none          |SWAT-MODFLOW dhru loop index
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use parm, only: nhru !SWAT
      use smrt_parm
      implicit none
      
!     Define local variables
      real swatVar(nhru)
      real smrtVar(dhru)
      integer i, j
      
!     Initialize variables
      smrtVar = 0.
      
!     Convert SWAT-HRU variables to SWAT-MODFLOW-dhru variables
!     The conversion does not involve a weighted average because the SWAT variables
!     are calculated such that the variable's value is the same within each DHRU
!     The Water Table in SWAT is a depth. The depth is the same everywhere in the HRU.
!     Thus, any dhru will also have the same value. Hence, do not apply a weighted average.
      do i=1,nhru
        do j=1,size(h2d_map(i)%dhru_id)
          smrtVar(h2d_map(i)%dhru_id(j)) = swatVar(i) !each dhru gets the hru value
        enddo
      enddo
      
      return
      end