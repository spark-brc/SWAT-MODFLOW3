      subroutine smrt_conversion2rt3d

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the necessary SWAT variables (per hru) into SWAT-MODFLOW variables
!!    (per disaggregated hru, dhru) then converts these into the proper units for RT3D
!!
!!    ~ ~ ~ Variables Used~ ~ ~
!!    name            |units                 |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    LENUNI          |unit_in               |the MODFLOW variable for which length units are being used
!!    ITMUNI          |unit_out              |the MODFLOW variable for which time units are being used
!!    nitraten(:)     |mg N/L                |nitrate concentration in reach
!!    disolvp(:)      |mg P/L                |dissolved phosphorus concentration in reach !rtb P
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name            |units                 |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    rchrg_no3_conc_dhru(:)  |mg/L                |concentration of nitrate in recharge water
!!    rchrg_p_conc_dhru(:)   |mg/L                |concentration of phosphorus in recharge water
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
      
      !Import variables
      use parm, only:leapyr,msub,nitraten,disolvp,gw_delaye,rchrg,
     &               percn,percp,rchrg_n,rchrg_p,nhru,iida,iyr,hru_ha,
     &               rchrg_no3_conc,rchrg_p_conc,day_total !SWAT
      use GLOBAL, only:LENUNI,ITMUNI,NCOL,NROW !MODFLOW
      use GWFRIVMODULE, only:RIVR,MXRIVR !MODFLOW
      use mf_rt_link, only: rt_criv,crch !MODFLOW-RT3D linkage

      use smrt_parm !smrt linkage
      implicit none
      
      !Define local variables
      integer mf_lengthUnit,mf_timeUnit,dum
      integer i,j,n,col,row,lay,subbasin_id,mf_gridID
      real conc_recharge(ncol,nrow)
      real rivlen,rivNitraten,rivPdsolv
      real no3_mass,p_mass,rchrgn1,rchrgp1,hru_area_m2,rch_volume

      mf_lengthUnit = LENUNI + 10
      mf_timeUnit = ITMUNI
      
      !compute nitrate and phosphate in recharge water for current day
      do n=1,nhru
        
        !compute mass of nitrate reaching water table
        rchrgn1 = 0.
        rchrgn1 = rchrg_n(n)
        if (rchrgn1 < 1.e-6) rchrgn1 = 0.0
        rchrg_n(n) = 0.
        rchrg_n(n) = ((1.- gw_delaye(n)) * percn(n)) + 
     &               (gw_delaye(n) * rchrgn1)
        no3_mass = rchrg_n(n) * hru_ha(n) * 1000. !g of nitrate-n

        !compute mass of phosphorus reaching water table
        rchrgp1 = 0.
        rchrgp1 = rchrg_p(n)
        if (rchrgp1 < 1.e-6) rchrgp1 = 0.0
        rchrg_p(n) = 0.
        rchrg_p(n) = ((1.- gw_delaye(n)) * percp(n)) + 
     &               (gw_delaye(n) * rchrgp1)
        p_mass = rchrg_p(n) * hru_ha(n) * 1000. !g of phosphorus

        !compute volume of recharge water for the HRU
        hru_area_m2 = hru_ha(n) * 10000. !m2
        rch_volume = (rchrg(n)/1000.) * hru_area_m2 !m3 of groundwater

        !compute concentration in recharge water
        if(rch_volume.gt.0) then
          rchrg_no3_conc(n) = no3_mass / rch_volume !g/m3 = mg/L
          rchrg_p_conc(n) = p_mass / rch_volume !g/m3 = mg/L
        else
          rchrg_no3_conc(n) = 0.
          rchrg_p_conc(n) = 0.
        endif

      enddo

      !write out deep percolation values to file
      if(day_total.eq.outswatmf(swatmf_out_ctr)) then
      write(30010,*) 'NO3 recharge (mg/L) for each HRU for day',
     &                iida,'year',iyr
      do i=1,nhru
        write(30010,*) rchrg_no3_conc(i)
      enddo
      write(30010,*)
      write(30014,*) 'Dissolved P recharge (mg/L) for each HRU 
     &                for day',iida,'year',iyr
      do i=1,nhru
        write(30014,*) rchrg_p_conc(i)
      enddo
      write(30014,*)
      endif


      !Convert SWAT hru based variables into SWAT-MODFLOW dis-aggregated hrus for more accurate surface/groundwater interaction
      call smrt_hru2dhru(rchrg_no3_conc, rchrg_no3_conc_dhru)
      call smrt_hru2dhru(rchrg_p_conc, rchrg_p_conc_dhru)
      

      !Convert SWAT-MODFLOW variables into RT3D variables
      call smrt_dhru2rtgrid3D(rchrg_no3_conc_dhru, conc_recharge)
      do i=1,nrow
        do j=1,ncol
          CRCH(j,i,1) = conc_recharge(j,i) !no3 concentration mapped to the 1st index of CRCH
        enddo
      enddo
      call smrt_dhru2rtgrid3D(rchrg_p_conc_dhru, conc_recharge)
      do i=1,nrow
        do j=1,ncol
          CRCH(j,i,2) = conc_recharge(j,i) !p concentration mapped to the 2nd index of CRCH
        enddo
      enddo


      !write out cell-by-cell values to a file (for testing/plotting)
      if(day_total.eq.outswatmf(swatmf_out_ctr)) then
      write(30009,*) 'NO3 recharge (mg/L) for each cell for day',
     &               iida,'year',iyr
      do i=1,nrow
        write(30009,100) (CRCH(j,i,1),j=1,ncol)
      enddo
      write(30009,*)
      write(30013,*) 'P recharge (mg/L) for each cell for day',
     &               iida,'year',iyr
      do i=1,nrow
        write(30013,100) (CRCH(j,i,2),j=1,ncol)
      enddo
      write(30013,*)
      endif


      !loop through the MODFLOW river cells that are intersected with SWAT's subbasins (as listed in the swatmf_river2grid.txt file)
      !for each cell, map the NO3 and P river water concentration to the RT3D array
      do i=1,nrivcells_subs
        
        !reset averaged river properties
        rivlen = 0.
        rivNitraten = 0.
        rivPdsolv = 0.

        !loop through the SWAT rivers (one for each subbasin)
        do j=1,riv_nsubs(i)
          
          !get the subbasin ID
          subbasin_id = grid2riv_id(i,j)
          
          !Get the river's properties to be based on a weighted average with: 
          !weights = subbasin's river segment length / total river length in grid cell (rivlen)
          rivlen = rivlen + grid2riv_len(i,j)
          rivNitraten = rivNitraten + 
     &                     (nitraten(subbasin_id) * grid2riv_len(i,j))
          rivPdsolv = rivPdsolv + 
     &                     (disolvp(subbasin_id) * grid2riv_len(i,j))
        enddo
        
        !prevent divide by zero problems
        if(rivlen.eq.0) rivlen = 1. 
        
        !take weighted average of river nitrate concentration
        rivNitraten = rivNitraten / rivlen !mg/L
        rivPdsolv = rivPdsolv / rivlen !mg/L

        !add the concentration values to the correct MODFLOW river cell
        do j=1,MXRIVR
          row = RIVR(2,j)
          col = RIVR(3,j)
          mf_gridID = ((row-1)*NCOL) + col !grid ID  
          if(mf_gridID.eq.grid2riv_gridID(i)) then
            rt_criv(j,1) = rivNitraten !note that nitrate is hard coded to the 1st entry for RT3D
            rt_criv(j,2) = rivPdsolv !note that phosphorus is hard coded to the 2nd entry for RT3D
          endif
        enddo

      enddo !go to the next cell


100   format (1000(f12.4))

      return
      end