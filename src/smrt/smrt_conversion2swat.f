      subroutine smrt_conversion2swat

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine converts all the MODFLOW variables, to pass into SWAT, into the proper units and SWAT variables
!!
!!    ~ ~ ~ Variables Used ~ ~ ~
!!    name            |units            |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    gw_q(:)         |mm H2O           |groundwater contribution to streamflow from
!!                    |                 |HRU on current day
!!    gw_qdeep(:)     |mm H2O           |groundwater contribution to streamflow from deep aquifer from
!!                    |                 |HRU on current day
!!    gwht(:)         |m                |groundwater height
!!    hru_fr(:)       |none             |fraction of subbasin area in HRU
!!    hru_ha(:)       |ha               |area of HRU in hectares
!!    hrutot(:)       |none             |number of HRUs in subbasin
!!    leapyr          |none             |leap year flag:
!!                    |                 |0  leap year
!!                    |                 |1  regular year
!!    msubo           |none             |max number of variables in output.sub
!!    nhru            |none             |number of HRUs in watershed
!!    no3gw(:)        |kg N/ha          |nitrate loading to reach in groundwater
!!                    |                 |in HRU for the day
!!    rchstor(:)      |m^3 H2O          |water stored in reach
!!    rttlc(:)        |m^3 H2O          |transmission losses from reach on day
!!    sub_fr(:)       |km2/km2          |fraction of total watershed area contained
!!                    |                 |in subbasin
!!    sub_km(:)       |km^2             |area of subbasin in square kilometers
!!    LENUNI          |MODFLOW length   |the MODFLOW variable for which length units are being used
!!    ITMUNI          |MODFLOW time     |the MODFLOW variable for which time units are being used
!!    NROW            |n/a              |the MODFLOW variable for number of rows of grids
!!    NCOL            |n/a              |the MODFLOW variable for number of columns of grids
!!    NLAY            |n/a              |the MODFLOW variable for number of layers of grids
!!    HNEW(:,:,:)     |LENUNI           |MODFLOW variable for the new head in the aquifer
!!    BOTM(:,:,:)     |LENUNI           |MODFLOW variable for the elevation of the 
!!                    |                 |bottom of each grid cell layer
!!    RIV(6,:)        |LENUNI           |MODFLOW variable for properties of the river cells
!!                    |n/a              |(1,:) cell layer index
!!                    |n/a              |(2,:) cell row index
!!                    |n/a              |(3,:) cell column index
!!                    |LENUNI           |(4,:) cell river stage (bottom elevation + depth)
!!                    |LENUNI/ITMUNI    |(5,:) cell conductance
!!                    |LENUNI           |(6,:) cell river bottom elevation
!!    grid2seg_id     |none             |a list per MODFLOW river grids (cols) containing 
!!                    |                 |the SWAT river segment IDs within that grid cell
!!    grid2seg_len    |m                |a list per MODFLOW river grids (cols) containing 
!!                    |                 |the SWAT river lengths within that grid cell
!!    nrivcells_subs  |n/a              |number of river grid cells in MODFLOW that are linked with SWAT sub-basins
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name            |units                 |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    gwht(:)         |m                     |groundwater height SWAT variable, populated by MODFLOW's HNEW
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!
!!    ~ ~ ~ Local Definitions ~ ~ ~
!!    name            |units                 |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mf_lengthUnit   |MODFLOW length unit   |the modified integer to represent
!!                    |                      |the MODFLOW unit of length for units.f
!!    mf_timeUnit     |MODFLOW time unit     |the modified integer to represent
!!                    |                      |the MODFLOW unit of time for units.f
!!    IL              |none                  |the layer index for getting info. from RIVR
!!    IC              |none                  |the column index for getting info. from RIVR
!!    IR              |none                  |the row index for getting info. from RIVR
!!    sum_rivrate(1)  |LENUNI**3/ITMUNI      |holder for the sum of river rates contributing to a
!!                    |                      |given subbasin to be mapped to SWAT
!!    wtlocation(:,:) |n/a                   |smrt variable of the rows and columns of 
!!                    |                      |MODFLOW grid cells with a value indicating
!!                    |                      |which layer the water table is located in
!!    gw_tot          |m**3/day              |The total amount of groundwater discharge from the current
!!                    |                      |grid cell to the subbasin outlet
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ End Specifications ~ ~ ~ ~ ~ ~
      
!     Import variables
      use parm, only: gw_q,gw_qdeep,gwht,hru_fr,hru_ha,hru1, !SWAT
     &                hrutot,day_total,
     &                leapyr,msub,nhru,rchstor,rttlc,sub_fr,
     &                sub_km,varoute,sw_gw_q,drn_q,mhru,hru_dafr, !rtb drain
     &                qdr,hrumono,vadose_thick_hru,hru_km,
     &                wshddayo,
     &                no3gw,sw_gw_no3,minpgw,sw_gw_p
      use GLOBAL,only:LENUNI,ITMUNI,NROW,NCOL,NLAY,HNEW,BOTM, !MODFLOW
     &                DELR,DELC,IBOUND,IUNIT
      use GWFUPWMODULE, only: SC2UPW  !MODFLOW
      use GWFRIVMODULE, only:RIVR,MXRIVR !MODFLOW
      use GWFWELMODULE, only:WELL,NWELLS
      use GWFDRNMODULE !MODFLOW Drain Package !rtb drain
      use GWFSTRMODULE, only:STRM !MODFLOW Stream Package
      use GWFEVTMODULE, only:EVTvol
      use GWFBASMODULE, only:VBVL
      use GWFUPWMODULE, only:LAYTYPUPW
      use mf_rt_link, only: rt_active,rt_rivmass,rt_drnmass !MODFLOW-RT3D Linkage
      use rt_global, only: CNEW !RT3D
      use smrt_parm !smrt linkage
      implicit none
      
!     define local variables
      integer found,mf_gridID
      integer mf_lengthUnit,mf_timeUnit,nsubs,dum
      integer IL,IC,IR,i,j,k,h,m,n,w,subID,dhruID,hru_id,
     &        well_row,well_col,cell_row,cell_col,sub_id,
     &        num_layer,ii,jj,kk,icomp
      DOUBLE PRECISION HHNEW,CHRIV,RRBOT,CCRIV
      real grid_rivlen, HRIV, CRIV, RBOT, RATE, gw_tot, discharge,
     &     loading,gw_sw_mm,sw_gw_mm,
     &     basin_area,sub_area,
     &     no3_gw_loading,p_gw_loading,
     &     sub_gw_sw,sub_sw_gw,
     &     sub_gw_sw_no3,sub_sw_gw_no3,sub_drn_sw_no3,
     &     sub_gw_sw_p,sub_sw_gw_p,sub_drn_sw_p,
     &     thickness,cell_area,
     &     cell_gw_volume,aquifer_gw_volume,sy,aquifer_gw_m,hru_gw,
     &     storativity,spec_storage
      real grams_N, Ntot
      real wtlocation(NCOL,NROW)
      real sum_rivrate(1)
      real rt_cnew_location(NCOL,NROW)
      real sub_gw(msub)
      real sub_gwno3(msub)
      real sub_gwsolp(msub)
      real gwht_dhru(dhru),vadose_thick_dhru(dhru)
      real sub_irrig_vol(msub)
      real Q_rate,Q_array(1),sum_mf_pump,sum_sub_pump,irrig_area
      real gw_change,gw_conhead,gw_et


      !rtb drain (DRAIN package --> flow to SWAT streams)
      real sub_drn(msub),sub_drn_no3(msub),sub_drn_p(msub)
      real drnrate(1),no3_drn_loading,p_drn_loading
      integer drn_lay,drn_row,drn_col,dhru_row,dhru_col,
     &        drn_hru,sub_basin
      real    drn_elev,cell_head,drn_cond,drn_flow,sub_drn_sw

      !rtb stream (STREAM package --> flow to SWAT streams)
      real sub_str(msub)
      real strrate(1)
      integer str_seg,str_rch,str_sub,nstrm_swat
      real sub_str_sw,str_flow

      mf_lengthUnit = LENUNI + 10
      mf_timeUnit = ITMUNI
      wtlocation = 1 !Set default water table location
      rt_cnew_location = 1 !Set UZF-R3D concentration location to the 1st layer, rtb
      subID = 0
      gwht_dhru = 0.
      vadose_thick_dhru = 0.
      gw_tot = 0
      h = 1 !set starting hru index
      
      !Zero out gwmod values. gw_q values will be populated using output from MODFLOW
      gw_q = 0.
      gw_qdeep = 0.
      sw_gw_q = 0.
      drn_q = 0. !rtb drain
      rttlc = 0.

      !zero out array for subbasin groundwater discharge and nitrate mass loading !rtb
      sub_gw = 0. 
      sub_gwno3 = 0.
      sub_gwsolp = 0.

      !zero out array for subbasin drainage water (from DRAIN package) !rtb drain
      sub_drn = 0.
      sub_str = 0.
      sub_drn_no3 = 0.
      sub_drn_p = 0.

      !get the total area of the basin
      basin_area = 0.
      do i=1,msub
        basin_area = basin_area + sub_km(i) !km2
      enddo

      !print out the area (km2) of each subbasin, for the first day of the simulation
      if(area_print.eq.1) then
        open(99999,file='sub_km')
        do i=1,msub
          write(99999,*) sub_km(i)
        enddo
        area_print = 0
      endif

      !output files (for MODFLOW, SWAT values)
      if(out_MODFLOW_gwsw.eq.1 .and.
     &   day_total.eq.outswatmf(swatmf_out_ctr)) then
        write(30005,*) !skip line, to start values for current day
        write(30005,*) 'Day:',day_total
        write(30005,*) 'Layer, Row, Column, Flow Rate'
      endif

      if(out_SWAT_gwsw.eq.1 .and.
     &   day_total.eq.outswatmf(swatmf_out_ctr)) then
        write(30006,*) !skip line, to start values for current day
        write(30006,*) 'Day:',day_total
        write(30006,*) 'Daily GW/SW Exchange for each Subbasin'
        write(30006,*) 'Subbasin, RIVER pkg, DRAIN pkg, STREAM pkg'
      endif

      if(rt_active.eq.1 .and.
     &   day_total.eq.outswatmf(swatmf_out_ctr)) then
        write(30007,*) !skip line, to start values for current day
        write(30007,*) 'Day:',day_total
        write(30007,*) 'Daily GW/SW NO3 Exchange for each River Cell'
        write(30008,*) !skip line, to start values for current day
        write(30008,*) 'Day:',day_total
        write(30008,*) 'Daily GW/SW NO3 Exchange for each Subbasin'
        write(30008,*) 'Subbasin, RIVER pkg, DRAIN pkg'
        write(30011,*) !skip line, to start values for current day
        write(30011,*) 'Day:',day_total
        write(30011,*) 'Daily GW/SW P Exchange for each River Cell'
        write(30012,*) !skip line, to start values for current day
        write(30012,*) 'Day:',day_total
        write(30012,*) 'Daily GW/SW P Exchange for each Subbasin'
        write(30012,*) 'Subbasin, RIVER pkg, DRAIN pkg'
      endif

      !loop through the MODFLOW river cells -------------------------------------------------------------------------------------
      if(IUNIT(4).gt.0) then !only proceed if River package is active
      do i=1,MXRIVR
        
        !determine if the River cell is linked with a SWAT sub-basin
        IR = RIVR(2,i) !row of MODFLOW river cell
        IC = RIVR(3,i) !column of MODFLOW river cell
        mf_gridID = ((IR-1)*NCOL) + IC !grid ID of MODFLOW river cell
        found = 0
        do m=1,nrivcells_subs
          
        !only proceed if there is a match
        if(mf_gridID.eq.grid2riv_gridID(m)) then
      
        !reset averaged river properties
        grid_rivlen = 0.
        sum_rivrate(1) = 0.
        RATE = 0.
        
        !Use original MODFLOW code to calculate the rate of loss/gain from/to rivers (LENUNI^3/ITMUNI)
        !the below code is borrowed from the MODFLOW subroutine LMT7RIV7 to calculate RATE
        IL=RIVR(1,i)
        IR=RIVR(2,i)
        IC=RIVR(3,i)
        !--GET RIVER PARAMETERS FROM RIVER LIST.
        HRIV=RIVR(4,i)
        CRIV=RIVR(5,i)
        RBOT=RIVR(6,i)
        HHNEW=HNEW(IC,IR,IL)
        CHRIV=CRIV*HRIV
        CCRIV=CRIV
        RRBOT=RBOT
        !--COMPARE HEAD IN AQUIFER TO BOTTOM OF RIVERBED.
        IF(HHNEW.GT.RRBOT) RATE = CHRIV-CCRIV*HHNEW !--AQUIFER HEAD > BOTTOM THEN RATE=CRIV*(HRIV-HNEW).
        IF(HHNEW.LE.RRBOT) RATE = CRIV*(HRIV-RBOT) !--AQUIFER HEAD < BOTTOM THEN RATE=CRIV*(HRIV-RBOT)
        if(out_MODFLOW_gwsw.eq.1 .and.
     &             day_total.eq.outswatmf(swatmf_out_ctr)) then
          write(30005,*) IL,IR,IC,RATE !write out for analysis !rtb
        endif
        smrt_GWSW_MF_tot_mo(i) = smrt_GWSW_MF_tot_mo(i) + RATE
        smrt_GWSW_MF_avg_mo(i) = smrt_GWSW_MF_tot_mo(i) / day_count_mo
        smrt_GWSW_MF_tot_yr(i) = smrt_GWSW_MF_tot_yr(i) + RATE
        smrt_GWSW_MF_avg_yr(i) = smrt_GWSW_MF_tot_yr(i) / day_count_yr
        !end MODFLOW code

        !retrieve groundwater nitrate mass loading
        if(rt_active.eq.1) then
          no3_gw_loading = rt_rivmass(i,1) !grams of no3-n, from MODFLOW river cell
          no3_gw_loading = no3_gw_loading / 1000. !kg of no3-n
          p_gw_loading = rt_rivmass(i,2) !grams of dissolved P, from MODFLOW river cell
          p_gw_loading = p_gw_loading / 1000. !kg of P
          if(day_total.eq.outswatmf(swatmf_out_ctr)) then
            write(30007,*) IL,IR,IC,no3_gw_loading
            write(30011,*) IL,IR,IC,p_gw_loading
          endif
        endif
        
        !The river cell might be in multiple subbasins. Split the groundwater discharge and mass loadings
        !to these subbasins based on stream length (usually, the cell is only in one subbasin).  !rtb

        !Get total stream segment length within current river cell !rtb
        do j=1,riv_nsubs(m) !rtb
          grid_rivlen = grid_rivlen + grid2riv_len(m,j)
        enddo
        
        !loop through the subbasins that have spatial areas within the current river cell  !rtb
        nsubs = riv_nsubs(m) !rtb
        do j=1,nsubs !rtb
          
          !Calculate the volume of groundwater discharge to the current subbasin, and add to
          !the total groundwater discharge for the subbasin  !rtb
          subID = grid2riv_id(m,j) !rtb
          
          !change to positive value (to indicate discharge from the aquifer to the stream) !rtb
          discharge = RATE * (-1) !rtb

          !convert to units used by SWAT !rtb
          call units(discharge, mf_lengthUnit, 12, 3, 1, leapyr)! convert length unit (LENUNI**3/ITMUNI to m**3/ITMUNI)
          call units(discharge, 4, mf_timeUnit, 1, 1, leapyr)! convert time unit (m**3/ITMUNI to m**3/day)
          !discharge = discharge * 86400. !ft3 for the day !rtb
          !discharge = discharge * (0.3048**3) !m3 for the day !rtb

          !portion of river cell discharge for the subbasin, based on stream length !rtb
          discharge = discharge * (grid2riv_len(m,j)/grid_rivlen) !rtb
            
          !add river cell discharge to the total groundwater discharge for the subbasin !rtb
          sub_gw(subID) = sub_gw(subID) + discharge !rtb
          smrt_GWSW_SWAT_tot_mo(subID) = smrt_GWSW_SWAT_tot_mo(subID) + 
     &                                   discharge
          smrt_GWSW_SWAT_avg_mo(subID) = smrt_GWSW_SWAT_tot_mo(subID) / 
     &                                   day_count_mo
          smrt_GWSW_SWAT_tot_yr(subID) = smrt_GWSW_SWAT_tot_yr(subID) + 
     &                                   discharge
          smrt_GWSW_SWAT_avg_yr(subID) = smrt_GWSW_SWAT_tot_yr(subID) / 
     &                                   day_count_yr

          !add river cell nitrate mass loading to the total groundwater mass loading for the subbasin !rtb
          if(rt_active.eq.1) then
            loading = no3_gw_loading * (-1) !positive now signifies loading to the stream
            loading = loading * (grid2riv_len(m,j)/grid_rivlen)
            sub_gwno3(subID) = sub_gwno3(subID) + loading
            loading = p_gw_loading * (-1) !positive now signifies loading to the stream
            loading = loading * (grid2riv_len(m,j)/grid_rivlen)
            sub_gwsolp(subID) = sub_gwsolp(subID) + loading
          endif

        enddo !rtb

          endif
        enddo !end loop to determine if MODFLOW River cell has a match with a cell that is linked with SWAT sub-basins     
      
      enddo !go to next MODFLOW River Cell

      !rtb drain
      !Determine drain water added to each sub-basin (from MODFLOW's drain package) -----------------------------------
      if(IUNIT(3).GT.0 .and. mf_drain_subs.eq.1) then
        
        !loop through the DRAIN cells
        do i=1,NDRAIN

          !get the cell index for the current DRAIN cell
          drn_lay = DRAI(1,i)
          drn_row = DRAI(2,i)
          drn_col = DRAI(3,i)

          !only proceed if the MODFLOW cell is active
          if(IBOUND(drn_col,drn_row,drn_lay).gt.0) then

            drn_elev = DRAI(4,i) !elevation of drain

            !only proceed if groundwater head is above the drain
            cell_head = HNEW(drn_col,drn_row,drn_lay)
            if(cell_head.gt.drn_elev) then
              drn_cond = DRAI(5,i) !conductance between drain and aquifer
              
              !calculate flow rate to the drain (and change units to m3/day)
              drn_flow = drn_cond * (cell_head - drn_elev) !L3/day
              drnrate(1) = drn_flow
              call units(drnrate, mf_lengthUnit, 12, 3, 1, leapyr)! convert length unit (LENUNI**3/ITMUNI to m**3/ITMUNI)
              call units(drnrate, 4, mf_timeUnit, 1, 1, leapyr)! convert time unit (m**3/ITMUNI to m**3/day)

              !calculate loading of N and P to the drain (kg)
              if(rt_active.eq.1) then
                no3_drn_loading = rt_drnmass(i,1)
                no3_drn_loading = no3_drn_loading / 1000. !g to kg
                p_drn_loading = rt_drnmass(i,2)
                p_drn_loading = p_drn_loading / 1000. !g to kg
              endif

              !find the sub-basin associated with this DRAIN cell
              sub_basin = drn_subs(drn_row,drn_col)
              if(sub_basin.gt.0) then
                
                !water
                sub_drn(sub_basin) = sub_drn(sub_basin) + drnrate(1)

                if(rt_active.eq.1) then
                !nitrate
                loading = no3_drn_loading * (-1) !positive now signifies loading to the stream
                sub_drn_no3(sub_basin) = sub_drn_no3(sub_basin)+loading

                !phosphorus
                loading = p_drn_loading * (-1) !positive now signifies loading to the stream
                sub_drn_p(sub_basin) = sub_drn_p(sub_basin)+loading

                if(isnan(loading)) then
                  dum = 10
                endif

                endif

              endif

            endif
          endif
        enddo !go to the next DRAIN cell

      endif !if DRAIN package is used      
      


      !rtb Stream
      !Determine stream water added to each sub-basin (from MODFLOW's Stream package) ---------------------------------
      !(this happens if a stream that is not in connection with SWAT, e.g. canal, discharges water to the river network
      if(IUNIT(18).GT.0) then

        !loop through stream segments that discharge water to SWAT sub-basins
        open(88889,file='swatmf_str2sub.txt')
        read(88889,*) nstrm_swat
        read(88889,*)
        do i=1,nstrm_swat
          read(88889,*) str_rch,str_sub
          
          !get flow rate from the Stream reach, and convert to SWAT units
          str_flow = STRM(9,str_rch) !L3/day
          strrate(1) = str_flow
          call units(strrate, mf_lengthUnit, 12, 3, 1, leapyr)! convert length unit (LENUNI**3/ITMUNI to m**3/ITMUNI)
          call units(strrate, 4, mf_timeUnit, 1, 1, leapyr)! convert time unit (m**3/ITMUNI to m**3/day)

          !store in the correct subbasin
          sub_str(str_sub) = sub_str(str_sub) + strrate(1)

        enddo
        close(88889)
      endif



      !Add discharge/seepage to SWAT routing arrays -------------------------------------------------------------------

      !Add subbasin groundwater discharge (River package) and drain water discharge (Drain package) 
      !to the total volume of water in the stream !rtb drain
      do j=1,msub !rtb
        varoute(2,j) = varoute(2,j) + sub_gw(j) !rtb WATER (m3)
        varoute(2,j) = varoute(2,j) + sub_drn(j) !rtb drain WATER (m3)
        varoute(2,j) = varoute(2,j) + sub_str(j) !rtb stream WATER (m3)
        varoute(6,j) = varoute(6,j) + sub_gwno3(j) !rtb NO3-N from aquifer-river interaction (kg)
        varoute(6,j) = varoute(6,j) + sub_drn_no3(j) !rtb NO3-N from drains (kg)
        varoute(7,j) = varoute(7,j) + sub_gwsolp(j) !rtb P from aquifer-river interaction (kg)
        varoute(7,j) = varoute(7,j) + sub_drn_p(j) !rtb P from drains (kg)
        if(out_SWAT_gwsw.eq.1) then
          
          write(30006,*) j,sub_gw(j),sub_drn(j),sub_str(j) !rtb drain

          if(rt_active.eq.1) then
            if(day_total.eq.outswatmf(swatmf_out_ctr)) then
              write(30008,*) j,sub_gwno3(j),sub_drn_no3(j)
              write(30012,*) j,sub_gwsolp(j),sub_drn_p(j)
            endif
          endif

        endif
      enddo  
      
      
      
      !Calculate water balance terms (gw_q, sw_gw, drn_q gw) ----------------------------------------------------------
      !rtb drain

      !calculate gw_q, sw_gw, and drn_q
      sub_gw_output = 0.
      do j=1,msub
        
        sub_gw_sw = sub_gw(j)
        sub_drn_sw = sub_drn(j) !rtb drain
        sub_gw_sw_no3 = sub_gwno3(j)
        sub_gw_sw_p = sub_gwsolp(j)
        sub_drn_sw_no3 = sub_drn_no3(j)
        sub_drn_sw_p = sub_drn_p(j)
        sub_area = sub_km(j)

        !change to mm/day
        sub_gw_sw = sub_gw_sw / (sub_km(j) * 1000000.) !m of water
        sub_gw_sw = sub_gw_sw * 1000. !mm of water
        sub_drn_sw = sub_drn_sw / (sub_km(j) * 1000000.) !m of water
        sub_drn_sw = sub_drn_sw * 1000. !mm of water
        sub_str_sw = sub_str_sw / (sub_km(j) * 1000000.) !m of water
        sub_str_sw = sub_str_sw * 1000. !mm of water

        !store for output.sub output
        sub_gw_output(j) = sub_gw_sw !mm of water over the subbasin
        
        !change to kg/ha
        sub_gw_sw_no3 = sub_gw_sw_no3 / sub_km(j) !kg/km2
        sub_gw_sw_no3 = sub_gw_sw_no3 / 100. !kg/ha
        sub_drn_sw_no3 = sub_drn_sw_no3 / sub_km(j) !kg/km2
        sub_drn_sw_no3 = sub_drn_sw_no3 / 100. !kg/ha
        sub_gw_sw_p= sub_gw_sw_p / sub_km(j) !kg/km2
        sub_gw_sw_p = sub_gw_sw_p / 100. !kg/ha
        sub_drn_sw_p = sub_drn_sw_p / sub_km(j) !kg/km2
        sub_drn_sw_p = sub_drn_sw_p / 100. !kg/ha

        if(sub_gw_sw.gt.0) then !groundwater entering stream
          
          !calculate % of subbasin groundwater contributing to entire basin
          wshddayo(104) = wshddayo(104) + 
     &                      (sub_gw_sw * (sub_area/basin_area))

          !calculate % of subbasin groundwater NO3 contributing to entire basin
          wshddayo(110) = wshddayo(110) + 
     &                      (sub_gw_sw_no3 * (sub_area/basin_area))
          
          !calculate % of subbasin groundwater P contributing to entire basin
          wshddayo(118) = wshddayo(118) + 
     &                      (sub_gw_sw_p * (sub_area/basin_area))
          
          !now, calculate gw_q for each HRU
          h = hru1(j)
          do k=1,hrutot(j)     
            !scale the value by the % hru area in the subbasin area
            gw_q(h) = sub_gw_sw * hru_fr(h) !mm of water
            no3gw(h) = sub_gw_sw_no3 * hru_fr(h) !kg/ha of NO3
            minpgw(h) = sub_gw_sw_p * hru_fr(h) !kg/ha of P
            h = h + 1
          enddo
        
        else !stream water entering groundwater
          sub_sw_gw = sub_gw_sw * (-1)
          sub_sw_gw_no3 = sub_gw_sw_no3 * (-1)
          sub_sw_gw_p = sub_gw_sw_p * (-1)

          !calculate % of subbasin sw_gw contributing to entire basin
          wshddayo(114) = wshddayo(114) + 
     &                      (sub_sw_gw * (sub_area/basin_area))
          
          !calculate % of subbasin sw_gw NO3 contributing to entire basin
          wshddayo(116) = wshddayo(116) + 
     &                      (sub_sw_gw_no3 * (sub_area/basin_area))
          
          !calculate % of subbasin sw_gw P contributing to entire basin
          wshddayo(119) = wshddayo(119) + 
     &                      (sub_sw_gw_p * (sub_area/basin_area))
          
          !now, calculate gw_q for each HRU
          h = hru1(j)
          do k=1,hrutot(j)     
            !scale the value by the % hru area in the subbasin area
            sw_gw_q(h) = sub_sw_gw * hru_fr(h) !mm of water
            sw_gw_no3(h) = sub_sw_gw_no3 * hru_fr(h) !kg/ha of NO3
            sw_gw_p(h) = sub_sw_gw_p * hru_fr(h) !kg/ha of P
            h = h + 1
          enddo
        endif

        !drainage discharge to streams
        !calculate % of subbasin drainage contributing to entire basin
        wshddayo(117) = wshddayo(117) + 
     &                    (sub_drn_sw * (sub_area/basin_area))
        
        !calculate % of subbasin drain NO3 contributing to entire basin
        wshddayo(120) = wshddayo(120) + 
     &                    (sub_drn_sw_no3 * (sub_area/basin_area))
         
        !calculate % of subbasin drain P contributing to entire basin
        wshddayo(121) = wshddayo(121) + 
     &                    (sub_drn_sw_p * (sub_area/basin_area))
    
        !now, calculate drn_q for each HRU
        h = hru1(j)
        do k=1,hrutot(j)     
          !scale the value by the % hru area in the subbasin area
          drn_q(h) = sub_drn_sw * hru_fr(h) !mm of water
          h = h + 1
        enddo
      
      enddo
      endif


      !calculate gw (total groundwater volume in the watershed) -----------------------------------------------------------------
      aquifer_gw_volume = 0.
      if(IUNIT(1)) then
        num_layer = mf_NTOP
      else
        num_layer = NLAY
      endif
      do i=1,NROW
        do j=1,NCOL
          do k=1,num_layer
            thickness = 0
            if(ibound(j,i,k).ne.0) then
            !Determine saturated thickness of the current layer
            if(HNEW(j,i,k).gt.BOTM(j,i,0)) then !head above the top of layer 1
              thickness = HNEW(j,i,k) - BOTM(j,i,0)
            else
              if(HNEW(j,i,k).lt.BOTM(j,i,k-1) .and. !BOTM(J,I,0) contains the top of layer 1
     &             HNEW(j,i,k).gt.BOTM(j,i,k)) then 
                thickness = HNEW(j,i,k) - BOTM(j,i,k)
              else
                thickness = BOTM(j,i,k-1) - BOTM(j,i,k)
              endif
            endif
            cell_area = DELR(j) * DELC(i)
            
            ! forrtl: severe (157): Program Exception - access violation !spark commented out
            !if(LAYTYPUPW(k).eq.0) then !confined aquifer
            !  spec_storage = mf_SC1(j,i,k) / cell_area
            !  storativity = spec_storage * thickness
            !  cell_gw_volume = thickness * cell_area * storativity
            !else !unconfined aquifer
            !  if(IUNIT(23).gt.0 .or. IUNIT(1).gt.0) then
            !    sy = mf_SC2(j,i,k) / cell_area !specific yield from cell storage capacity
            !  else
            !    sy = SC2UPW(j,i,k) / cell_area !specific yield from cell storage capacity
            !  endif
            !  cell_gw_volume = thickness * cell_area * sy
            !endif
            
            if(IUNIT(23).gt.0 .or. IUNIT(1).gt.0) then
              sy = mf_SC2(j,i,k) / cell_area !specific yield from cell storage capacity
            else
              sy = SC2UPW(j,i,k) / cell_area !specific yield from cell storage capacity
            endif
            cell_gw_volume = thickness * cell_area * sy            
            
            
            gw_available(j,i,k) = cell_gw_volume !save for SWAT auto-irrigation (see irrigate.f)
            aquifer_gw_volume = aquifer_gw_volume + cell_gw_volume
            endif
          enddo
        enddo
      enddo
      aquifer_gw_m = aquifer_gw_volume/(basin_area*1000000.) !convert m3 to m
      wshddayo(115) = aquifer_gw_m * 1000. !convert m to mm of water

      !include gw, sw_gw, and drn_q values into the overall water yield for the watershed !rtb drain
      wshddayo(6) = wshddayo(6) + wshddayo(104) - wshddayo(114) + 
     &              wshddayo(117)

      !change in groundwater volume (m3)
      gw_change = (VBVL(3,1) - VBVL(4,1)) * (-1)
      call units(gw_change, mf_lengthUnit, 12, 3, 1, leapyr) !length
      call units(gw_change, 4, mf_timeUnit, 1, 1, leapyr) !time
      gw_change = gw_change / (basin_area*1000000.) !convert to m
      gw_change = gw_change * 1000. !convert m to mm of water
      wshddayo(122) = gw_change

      !total groundwater added/removed via constant head boundaries
      gw_conhead = (VBVL(3,2) - VBVL(4,2))
      call units(gw_conhead, mf_lengthUnit, 12, 3, 1, leapyr) !length
      call units(gw_conhead, 4, mf_timeUnit, 1, 1, leapyr) !time
      gw_conhead = gw_conhead / (basin_area*1000000.) !convert to m
      gw_conhead = gw_conhead * 1000. !convert m to mm of water
      wshddayo(123) = gw_conhead

      !total groundwater removed via ET (EVT package)
      gw_et = 0.
      if(iunit(5).gt.0) then
        do i=1,NROW
          do j=1,NCOL
            do k=1,NLAY
              gw_et = gw_et + (EVTvol(j,i,k)*(-1)) !make a positive value
            enddo
          enddo
        enddo
      endif
      gw_et = gw_et / (basin_area*1000000.) !convert to m
      gw_et = gw_et * 1000. !convert m to mm of water
      wshddayo(124) = gw_et

      !Add gw_q and sw_gw to total water yield (to be output in output.std)
      hru_gw = 0.
      do h=1,mhru
        hru_gw = hru_gw + (gw_q(h) * hru_dafr(h))
      enddo

      !HRU monthly summations
      do h=1,mhru
        hrumono(6,h) = hrumono(6,h) + gw_q(h)
        hrumono(70,h) = 0. !no gw_qdeep
        hrumono(74,h) = hrumono(74,h) + sw_gw_q(h)
      enddo



      !Provide irrigation water to SWAT HRUs, from MODFLOW Well package cells ---------------------------------------------------
      if(IUNIT(2).gt.0 .and. mf_irrigation.eq.1 .and.
     &                       mf_irrigation_swat.eq.0) then
        if(ncell_irrigate.gt.0) then !only if there are designated irrigation wells

          !allocate array to hold irrigation volumes for each sub-basin
          sub_irrig_vol = 0.
          mf_sub_irrig_flag = 0

          !global sums, for output check
          sum_mf_pump = 0.
          sum_sub_pump = 0.

          !loop through the active MODFLOW Well package cells. Determine if they are irrigation wells,
          !and sum of water applied to each sub-basin for the current day.
          do w=1,NWELLS

            well_row = WELL(2,w)
            well_col = WELL(3,w)
            sum_mf_pump = sum_mf_pump + (WELL(4,w)*-1)
            
            !loop through the cells designated as irrigation cells, and determine if there is a match
            do n=1,ncell_irrigate
              sub_id = mf_sub_irrigate(n,1)
              cell_row = mf_sub_irrigate(n,2)
              cell_col = mf_sub_irrigate(n,3)
              if(cell_row.eq.well_row .and. cell_col.eq.well_col) then
                mf_sub_irrig_flag(sub_id) = 1 !set flag indicating irrigation for the sub-basin
                Q_rate = WELL(4,w) * (-1) !get pumping rate for the cell
                goto 10
              endif
            enddo
            
            !convert pumping rate to m3/day
 10         continue
            Q_array(1) = Q_rate
            call units(Q_array,ITMUNI,4,1,1,1) !convert to days
            call units(Q_array,mf_lengthUnit,12,3,1,1) !convert to m3

            !add volume (in m3) to the sub-basin total
            sub_irrig_vol(sub_id) = sub_irrig_vol(sub_id) + Q_array(1)

          enddo !go to the next WELL cell

          !sum up volume of irrigation water provided to SWAT sub-basins
          do m=1,msub
            sum_sub_pump = sum_sub_pump + sub_irrig_vol(m)
          enddo
          write(300009,*) sum_mf_pump,sum_sub_pump

          !multiply by conveyance efficiency (if there are losses between the pumping well and the field of application
          do m=1,msub
            sub_irrig_vol(m) = sub_irrig_vol(m) * mf_convey_ratio(m)
          enddo

          !determine depth of irrigation water for the HRUs in each sub-basin
          do m=1,msub
              
            if(mf_nsub_hru_irrig(m).gt.0) then
            
            !calculate total area in sub-basin that receives irrigation water
            irrig_area = 0.
            do n=1,mf_nsub_hru_irrig(m)
              hru_id = mf_sub_hru_irrig(m,n)
              irrig_area = irrig_area + hru_km(hru_id) !km2 of land area
            enddo
              
            !calculate depth of irrigation water
            irrig_area = irrig_area * 1000000. !m2 of land area
            mf_sub_irrig_depth(m) = sub_irrig_vol(m) / irrig_area !m of water
            mf_sub_irrig_depth(m) = mf_sub_irrig_depth(m) * 1000. !mm of water

            endif
          enddo

        endif
      endif
      mf_hru_irrig = 0. !keeps track of total applied irrigation water to HRUs (filled up in "percmain.f")



      !Determine water table location (layer) and the thickness of the vadose zone ----------------------------------------------
      !(this will be used to determine if the water table is in SWAT's soil profile, and
      !(thus whether the soil processes are affected by the water table)

      !Loop through MODFLOW variables and pull out information
      do j=1, NROW
        do i=1, NCOL
          do k=1, NLAY
            !Loop through and find which layer the water table is in
            if(HNEW(i,j,k).LT.BOTM(i,j,k-1) .and. !BOTM(J,I,0) contains the top of layer1
     &           HNEW(i,j,k).GT.BOTM(i,j,k)) then 
              wtlocation(i,j) = k

              !distance between the ground surface and the water table
              vadose_thick(i,j) = BOTM(i,j,0) - HNEW(i,j,k)
            endif
          enddo
        enddo
      enddo
      
      !Fill in SWAT variable "gwht" with water table elevation from MODFLOW
      call smrt_grid2dhru3D(HNEW, gwht_dhru, wtlocation)
      call units(gwht_dhru, mf_lengthUnit, 12, 1, dhru, leapyr)!convert length unit (LENUNI to m)
      call smrt_dhru2hru(gwht_dhru, gwht, 0) !map dhru to hru
      
      !Convert MODFLOW variable (vadose_thick) into a SWAT variable (first to DHRUs, then to HRUs)
      !(now, each HRU has a thickness of the vadose zone - this is to know if the water table
      ! is within the soil profile)
      call smrt_grid2dhru2D(vadose_thick, vadose_thick_dhru)
      call units(vadose_thick_dhru, mf_lengthUnit, 14, 1, dhru, leapyr)! to convert length units (LENUNI to mm)
      call units(vadose_thick_dhru, mf_timeUnit, 4, 1, dhru, leapyr)! to convert time units (ITMUNI to days)
      call smrt_dhru2hru(vadose_thick_dhru, vadose_thick_hru, 0)
       
      !increment counter (if necessary)
      if(day_total.eq.outswatmf(swatmf_out_ctr)) then
        swatmf_out_ctr = swatmf_out_ctr + 1
      endif 

  100 format(1000(e12.4))
      
      return
      end