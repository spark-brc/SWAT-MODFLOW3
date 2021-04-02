      subroutine smrt_conversion2mf

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the necessary SWAT variables (per hru) into SWAT-MODFLOW variables
!!    (per disaggregated hru, dhru) then converts these into the proper units for MODFLOW
!!
!!    ~ ~ ~ Variables Used~ ~ ~
!!    name            |units                 |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    LENUNI          |unit_in               |the MODFLOW variable for which length units are being used
!!    ITMUNI          |unit_out              |the MODFLOW variable for which time units are being used
!!    sepbtm(:)       |mm H2O                |percolation from bottom of soil profile for
!!                    |                      |the day in HRU
!!    etremain(:)     |mm H2O                |remaining et to be passed from SWAT to
!!                    |                      |MODFLOW-UZF, etremain(j) = pet_day-etday
!!    etremain_dis(:) |mm H2O           (in) |remaining et to be passed from SWAT to
!!                    |LENUNI**3/ITMUNI (out)|MODFLOW-UZF per disaggregated hru (dhru)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name            |units                 |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    rchrg_dhru(:)   |mm H2O           (in) |recharge to the water table
!!                    |                      |populated from rchrg and smrt_hru2dhru conversion
!!                    |LENUNI**3/ITMUNI (out)|for the day in HRU, now in MODFLOW units
!!    etremain_dis(:) |mm H2O           (in) |remaining et to be passed from SWAT to
!!                    |LENUNI**3/ITMUNI (out)|MODFLOW-UZF per disaggregated hru (dhru)
!!                    |                      |populated from etremain and smrt_hru2dhru conversion
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name          |units                 |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mf_lengthUnit |MODFLOW length unit   |the modified inteteger to represent
!!                  |                      |the MODFLOW unit of length for units.f
!!    mf_timeUnit   |MODFLOW time unit     |the modified integer to represent
!!                  |                      |the MODFLOW unit of time for units.f
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
      
!     Import variables
      use parm, only:nhru,hru_km,rchrg,gw_delaye,sepbtm,day_total,
     &               etremain,stsol_rd,leapyr,hru_dafr,wshddayo  !SWAT
      use GLOBAL, only:LENUNI,ITMUNI !MODFLOW
      use smrt_parm !smrt linkage
      implicit none
      
!     Define local variables
      integer mf_lengthUnit,mf_timeUnit,n,j,k,dhruID,cell_row,cell_col
      real    hru_area,rchrg1
      mf_lengthUnit = LENUNI + 10
      mf_timeUnit = ITMUNI
      
      !calculate recharge to the water table, based on the percolation from the soil profile
      do n=1,nhru
        rchrg1 = rchrg(n)
        rchrg(n) = 0.
        rchrg(n) = ((1.-gw_delaye(n))*sepbtm(n)) + (gw_delaye(n)*rchrg1)
        if (rchrg(n) < 1.e-6) rchrg(n) = 0.
        wshddayo(107) = wshddayo(107) + rchrg(n) * hru_dafr(n) !watershed amount
      enddo
      
      !calculate average deep percolation for the current month (rtb avg)
      smrt_DP_tot_mo = smrt_DP_tot_mo + rchrg
      smrt_DP_avg_mo = smrt_DP_tot_mo / day_count_mo

      !calculate average recharge for the current year (rtb avg)
      smrt_DP_tot_yr = smrt_DP_tot_yr + rchrg
      smrt_DP_avg_yr = smrt_DP_tot_yr / day_count_yr

      !write out deep percolation values to file
      if(out_SWAT_recharge.eq.1 .and.
     &   day_total.eq.outswatmf(swatmf_out_ctr)) then
        write(30001,*)
        write(30001,*) 'Day:',day_total
        write(30001,*) 'Recharge (mm), Soil Perc. (mm) for each HRU'
        do n=1,nhru
          write(30001,*) rchrg(n),sepbtm(n)
        enddo
        write(30001,*)
      endif

!     Convert SWAT hru based variables into SWAT-MODFLOW dis-aggregated hrus for more accurate surface/groundwater interaction
      call smrt_hru2dhru(rchrg, rchrg_dhru) !mm
      call smrt_hru2dhru(etremain, etremain_dhru) !mm
      call smrt_hru2dhru(stsol_rd, stsol_rd_dhru) !mm
      
!     Convert SWAT-MODFLOW variables into MODFLOW units
      call units(rchrg_dhru, 14, mf_lengthUnit, 1, dhru, leapyr)! to convert units (mm/day to LENUNI/day)
      call units(rchrg_dhru, mf_timeUnit, 4, 1, dhru, leapyr)! to convert time units (LENUNI/days to LENUNI/ITMUNI)
      call units(etremain_dhru, 14, mf_lengthUnit, 1, dhru, leapyr)! to convert length units (mm/day to LENUNI/day)
      call units(etremain_dhru, mf_timeUnit, 4, 1, dhru, leapyr)! to convert time units (LENUNI/days to LENUNI/ITMUNI)
      call units(stsol_rd_dhru, 14, mf_lengthUnit, 1, dhru, leapyr)! to convert length units (mm/day to LENUNI/day)
      call units(stsol_rd_dhru, mf_timeUnit, 4, 1, dhru, leapyr)! to convert time units (LENUNI/days to LENUNI/ITMUNI)


      !Calculate applied irrigation volume for each HRU ---------------------------------------------------------------
      !This will be converted to pumping rates for MODFLOW well cells, when MODFLOW is run
      if(mf_irrigation_swat.eq.1) then
        !loop through the HRUs
        do n=1,nhru

          !find area of HRU (m2)
          hru_area = hru_km(n) * 1000000.

          !find volume of applied irrigation water
          irrig_volume(n) = (irrig_depth(n)/1000.) * hru_area !m3

        enddo
      endif
      irrig_depth = 0. !zero out, to prepare for SWAT's next day calculations

      return
      end