      subroutine units(array,unit_in,unit_out,magnitude,asize,leapyr)
!!    ~ ~ ~ Author ~ ~ ~
!!    Tyler Wible, Masters student
!!    Colorado State University 2012-2014
!!    Comment initials "tcw"
!!
!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts the provided array to different units based on 
!!    provided flags for the incoming and outgoing units (based on a modified 
!!    set of MODFLOW's unit flags listed below)
!!
!!    0 = undefined
!!    1 = seconds
!!    2 = minutes
!!    3 = hours
!!    4 = days
!!    5 = years    Note: a year is assumed to be equal to 365 days
!!
!!    11 = feet
!!    12 = meters
!!    13 = centimeters
!!    14 = millimeters
!!    15 = kilometers
!!
!!    21 = kilograms
!!    22 = grams
!!    23 = milligrams
!!    24 = micrograms
!!
!!    31 = square feet
!!    32 = square meters
!!    33 = hectacre (ha)
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    array(:)    |unit_in units |the variable to have its units converted
!!                |              |currently in the units of unit_in that has
!!                |              |"asize" elements in the array
!!    unit_in     |none          |integer representing the incoming unit of the 
!!                |              |parameter following the codes listed above
!!    unit_out    |none          |integer representing the outgoing unit of the 
!!                |              |parameter following the codes listed above
!!    magnitude   |none          |the exponent of the units to be converted
!!                |              |example: if converting from square feet (ft^2)
!!                |              |to square meeters (m^2), magnitude should have
!!                |              |a value of 2
!!    asize       |              |the number of elements in array(:)
!!    leapyr      |none          |leap year flag, this is only used if the in/out
!!                |              |units are time units involving years
!!                |              |0  leap year
!!                |              |1  regular year
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    array       |unit_out units|the originally provided variable array now 
!!                |              |in the units of unit_out
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
      
!     Initialize local variables
      implicit none
      integer asize, leapyr
      real array(asize)
      real conversion
      integer unit_in, unit_out, magnitude, i
      conversion = 1.
      
!     Determine the time unit conversion
      if (unit_in.eq.1 .and. unit_out.eq.2) then
        conversion = 1./60.        !!Convert seconds to minutes
      else if (unit_in.eq.1 .and. unit_out.eq.3) then
        conversion = 1./3600.      !!Convert seconds to hours
      else if (unit_in.eq.1 .and. unit_out.eq.4) then
        conversion = 1./86400.     !!Convert seconds to days
      else if (unit_in.eq.1 .and. unit_out.eq.5 .and. leapyr.ne.0) then
        conversion = 1./31536000.  !!Convert seconds to non-leap years
      else if (unit_in.eq.1 .and. unit_out.eq.5 .and. leapyr.eq.0) then
        conversion = 1./31622400.  !!Convert seconds to leap years
      
      else if (unit_in.eq.2 .and. unit_out.eq.1) then
        conversion = 60.           !!Convert minutes to seconds
      else if (unit_in.eq.2 .and. unit_out.eq.3) then
        conversion = 1./60.        !!Convert minutes to hours
      else if (unit_in.eq.2 .and. unit_out.eq.4) then
        conversion = 1./1440.      !!Convert minutes to days
      else if (unit_in.eq.2 .and. unit_out.eq.5 .and. leapyr.ne.0) then
        conversion = 1./525600.    !!Convert minutes to non-leap years
      else if (unit_in.eq.2 .and. unit_out.eq.5 .and. leapyr.eq.0) then
        conversion = 1./527040.    !!Convert minutes to leap years
      
      else if (unit_in.eq.3 .and. unit_out.eq.1) then
        conversion = 3600.         !!Convert hours to seconds
      else if (unit_in.eq.3 .and. unit_out.eq.2) then
        conversion = 60.           !!Convert hours to minutes
      else if (unit_in.eq.3 .and. unit_out.eq.4) then
        conversion = 1./24.        !!Convert hours to days
      else if (unit_in.eq.3 .and. unit_out.eq.5 .and. leapyr.ne.0) then
        conversion = 1./8760.      !!Convert hours to non-leap years
      else if (unit_in.eq.3 .and. unit_out.eq.5 .and. leapyr.eq.0) then
        conversion = 1./8784.      !!Convert hours to leap years
      
      else if (unit_in.eq.4 .and. unit_out.eq.1) then
        conversion = 86400.        !!Convert days to seconds
      else if (unit_in.eq.4 .and. unit_out.eq.2) then
        conversion = 1440.         !!Convert days to minutes
      else if (unit_in.eq.4 .and. unit_out.eq.3) then
        conversion = 24.           !!Convert days to hours
      else if (unit_in.eq.4 .and. unit_out.eq.5 .and. leapyr.ne.0) then
        conversion = 1./365.       !!Convert days to non-leap years
      else if (unit_in.eq.4 .and. unit_out.eq.5 .and. leapyr.eq.0) then
        conversion = 1./366.       !!Convert days to leap years
      
      else if (unit_in.eq.5 .and. unit_out.eq.1 .and. leapyr.ne.0) then
        conversion = 31536000.     !!Convert non-leap years to seconds
      else if (unit_in.eq.5 .and. unit_out.eq.2 .and. leapyr.ne.0) then
        conversion = 525600.       !!Convert non-leap years to minutes
      else if (unit_in.eq.5 .and. unit_out.eq.3 .and. leapyr.ne.0) then
        conversion = 8760.         !!Convert non-leap years to hours
      else if (unit_in.eq.5 .and. unit_out.eq.4 .and. leapyr.ne.0) then
        conversion = 365.          !!Convert non-leap years days
      
      else if (unit_in.eq.5 .and. unit_out.eq.1 .and. leapyr.eq.0) then
        conversion = 31622400.     !!Convert leap years to seconds
      else if (unit_in.eq.5 .and. unit_out.eq.2 .and. leapyr.eq.0) then
        conversion = 527040.       !!Convert leap years to minutes
      else if (unit_in.eq.5 .and. unit_out.eq.3 .and. leapyr.eq.0) then
        conversion = 8784.         !!Convert leap years to hours
      else if (unit_in.eq.5 .and. unit_out.eq.4 .and. leapyr.eq.0) then
        conversion = 366.          !!Convert leap years days
      
      
      
!     Determine the distance unit conversion
      else if (unit_in.eq.11 .and. unit_out.eq.12) then
        conversion = 1./3.28084    !!Convert feet to meters
      else if (unit_in.eq.11 .and. unit_out.eq.13) then
        conversion = 1./0.0328084  !!Convert feet to centimeters
      else if (unit_in.eq.11 .and. unit_out.eq.14) then
        conversion = 1./0.00328084 !!Convert feet to millimeters
      else if (unit_in.eq.11 .and. unit_out.eq.15) then
        conversion = 1./3280.84    !!Convert feet to kilometers
      
      else if (unit_in.eq.12 .and. unit_out.eq.11) then
        conversion = 3.28084       !!Convert meters to feet
      else if (unit_in.eq.12 .and. unit_out.eq.13) then
        conversion = 100.          !!Convert meters to centimeters
      else if (unit_in.eq.12 .and. unit_out.eq.14) then
        conversion = 1000.         !!Convert meters to millimeters
      else if (unit_in.eq.12 .and. unit_out.eq.15) then
        conversion = 1./1000.      !!Convert meters to kilometers
      
      else if (unit_in.eq.13 .and. unit_out.eq.11) then
        conversion = 0.0328084     !!Convert centimeters to feet
      else if (unit_in.eq.13 .and. unit_out.eq.12) then
        conversion = 1./100.       !!Convert centimeters to meters
      else if (unit_in.eq.13 .and. unit_out.eq.14) then
        conversion = 10.           !!Convert centimeters to millimeters
      else if (unit_in.eq.13 .and. unit_out.eq.15) then
        conversion = 1./100000.    !!Convert centimeters to kilometers
      
      else if (unit_in.eq.14 .and. unit_out.eq.11) then
        conversion = 0.00328084    !!Convert millimeters to feet
      else if (unit_in.eq.14 .and. unit_out.eq.12) then
        conversion = 1./1000.      !!Convert millimeters to meters
      else if (unit_in.eq.14 .and. unit_out.eq.13) then
        conversion = 1./10.        !!Convert millimeters to centimeters
      else if (unit_in.eq.14 .and. unit_out.eq.15) then
        conversion = 1./1000000.   !!Convert millimeters to kilometers
      
      else if (unit_in.eq.15 .and. unit_out.eq.11) then
        conversion = 3280.84       !!Convert kilometers to feet
      else if (unit_in.eq.15 .and. unit_out.eq.12) then
        conversion = 1000.         !!Convert kilometers to meters
      else if (unit_in.eq.15 .and. unit_out.eq.13) then
        conversion = 100000.       !!Convert kilometers to centimeters
      else if (unit_in.eq.15 .and. unit_out.eq.14) then
        conversion = 1000000.      !!Convert kilometers to milimeters
      
      
      
!     Determine the weight unit conversion
      else if (unit_in.eq.21 .and. unit_out.eq.22) then
        conversion = 1000.         !!Convert kilograms to grams
      else if (unit_in.eq.21 .and. unit_out.eq.23) then
        conversion = 1000000.      !!Convert kilograms to milligrams
      else if (unit_in.eq.21 .and. unit_out.eq.24) then
        conversion = 1000000000.   !!Convert kilograms to micrograms
      
      else if (unit_in.eq.22 .and. unit_out.eq.21) then
        conversion = 1./1000.      !!Convert grams to kilograms
      else if (unit_in.eq.22 .and. unit_out.eq.23) then
        conversion = 1000.         !!Convert grams to milligrams
      else if (unit_in.eq.22 .and. unit_out.eq.24) then
        conversion = 1000000.      !!Convert grams to micrograms
      
      else if (unit_in.eq.23 .and. unit_out.eq.21) then
        conversion = 1./1000000.   !!Convert milligrams to kilograms
      else if (unit_in.eq.23 .and. unit_out.eq.22) then
        conversion = 1./1000.      !!Convert milligrams to grams
      else if (unit_in.eq.23 .and. unit_out.eq.24) then
        conversion = 1000.         !!Convert milligrams to micrograms
      
      else if (unit_in.eq.24 .and. unit_out.eq.21) then
        conversion = 1./1000000000.!!Convert micrograms to kilograms
      else if (unit_in.eq.24 .and. unit_out.eq.22) then
        conversion = 1./1000000.   !!Convert micrograms to grams
      else if (unit_in.eq.24 .and. unit_out.eq.23) then
        conversion = 1./1000.      !!Convert micrograms to milligrams
      
      
      
!     Determine the area unit conversion
      else if (unit_in.eq.31 .and. unit_out.eq.32) then
        conversion = 1./10.7639    !!Convert square feet to square meters
      else if (unit_in.eq.31 .and. unit_out.eq.33) then
        conversion = 1./107639.    !!Convert square feet to hectacres
      
      else if (unit_in.eq.32 .and. unit_out.eq.31) then
        conversion = 10.7639       !!Convert square meters to square feet
      else if (unit_in.eq.32 .and. unit_out.eq.33) then
        conversion = 0.0001        !!Convert square meters to hectacres
      
      else if (unit_in.eq.33 .and. unit_out.eq.31) then
        conversion = 107639.       !!Convert hectacres to square feet
      else if (unit_in.eq.33 .and. unit_out.eq.32) then
        conversion = 10000.        !!Convert hectacres to square meters
      
        
      else
        conversion = 1.            !!Conversion Error, conversion not accounted for
      endif
      
      
!     Convert line by line the array into its new units
      do i=1,asize
        array(i) = array(i) * (conversion**magnitude)
      enddo
      
      return
      end