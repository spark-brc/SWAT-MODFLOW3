      subroutine smrt_allocate()

      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine allocates array sizes for variables used in the linkage between SWAT, MODFLOW, and RT3D

      use parm, only: msub,nhru,idaf !SWAT
      use smrt_parm !smrt linkage
      use mf_rt_link
      use rt_global, only: ncomp
      use GLOBAL, only:NCOL,NROW,NLAY,IBOUND !MODFLOW
      use GWFBASMODULE, only:HNOFLO
      use GWFRIVMODULE, only:MXRIVR 
      implicit none
      
      integer i,j,k

!     Allocate the size of the river-to-grid variables
      allocate(grid2riv_id(nrivcells_subs,msub))
      allocate(grid2riv_len(nrivcells_subs,msub))
      allocate(riv_nsubs(nrivcells_subs))
      allocate(grid2riv_gridID(nrivcells_subs))
      grid2riv_id = 0.
      grid2riv_len = 0.
      riv_nsubs = 0.
      grid2riv_gridID = 0
      
!     Allocate the size of HRU variables before smrt_upfluxToSoil runs
      allocate(upflux_no3(nhru))
      allocate(vadose_thick(NCOL,NROW))
      upflux_no3 = 0.
      vadose_thick = 0.

!     Allocate the size of disaggregated hru variables before MODFLOW and RT3D runs
      allocate(etremain_dhru(dhru))
      allocate(stsol_rd_dhru(dhru))
      allocate(rchrg_dhru(dhru))
      allocate(rchrg_no3_conc_dhru(dhru))
      allocate(rchrg_p_conc_dhru(dhru))
      allocate(upflux_no3_dhru(dhru))
      allocate(qtile_hru(nhru)) !rtb tile
      allocate(qtile_cells(NCOL,NROW))
      etremain_dhru = 0.
      stsol_rd_dhru = 0.
      rchrg_dhru = 0.
      rchrg_no3_conc_dhru = 0.
      rchrg_p_conc_dhru = 0.
      upflux_no3_dhru = 0.
      qtile_hru = 0. !rtb tile
      qtile_cells = 0.
      mf_tile_active = 0 !rtb tile
      
      !allocate and populate the month_days array (for stored variable averages) (rtb avg)
      allocate(month_days(12))
      month_days(1) = 31
      month_days(2) = 59
      month_days(3) = 90
      month_days(4) = 120
      month_days(5) = 151
      month_days(6) = 181
      month_days(7) = 212
      month_days(8) = 243
      month_days(9) = 273
      month_days(10) = 304
      month_days(11) = 334
      month_days(12) = 365

      allocate(month_days_leap(12))
      month_days_leap(1) = 31
      month_days_leap(2) = 60
      month_days_leap(3) = 91
      month_days_leap(4) = 121
      month_days_leap(5) = 152
      month_days_leap(6) = 183
      month_days_leap(7) = 213
      month_days_leap(8) = 244
      month_days_leap(9) = 274
      month_days_leap(10) = 305
      month_days_leap(11) = 335
      month_days_leap(12) = 366

      !initialize counters
      day_count_mo = 1
      day_count_yr = 1

      !initialize the month (based on the "idaf" variable: starting Julian day of the SWAT simulation)
      do i=1,12
        if(idaf.lt.month_days(i)) then
          smrt_month_counter = i
          goto 10
        endif
      enddo

      !arrays for variable averages
10    allocate(smrt_RECH_tot_mo(ncol,nrow))
      allocate(smrt_RECH_avg_mo(ncol,nrow))
      allocate(smrt_RECH_tot_yr(ncol,nrow))
      allocate(smrt_RECH_avg_yr(ncol,nrow))
      smrt_RECH_tot_mo = 0.
      smrt_RECH_avg_mo = 0.
      smrt_RECH_tot_yr = 0.
      smrt_RECH_avg_yr = 0.

      allocate(smrt_HNEW_tot_mo(ncol,nrow,nlay))
      allocate(smrt_HNEW_avg_mo(ncol,nrow,nlay))
      allocate(smrt_HNEW_tot_yr(ncol,nrow,nlay))
      allocate(smrt_HNEW_avg_yr(ncol,nrow,nlay))
      do k=1,nlay
        do j=1,ncol
          do i=1,nrow
            if(ibound(j,i,k).eq.1) then
              smrt_HNEW_tot_mo(j,i,k) = 0.
              smrt_HNEW_avg_mo(j,i,k) = 0.
              smrt_HNEW_tot_yr(j,i,k) = 0.
              smrt_HNEW_avg_yr(j,i,k) = 0.
            else
              smrt_HNEW_tot_mo(j,i,k) = HNOFLO
              smrt_HNEW_avg_mo(j,i,k) = HNOFLO
              smrt_HNEW_tot_yr(j,i,k) = HNOFLO
              smrt_HNEW_avg_yr(j,i,k) = HNOFLO
            endif
          enddo
        enddo
      enddo
      
      allocate(smrt_DP_tot_mo(nhru))
      allocate(smrt_DP_avg_mo(nhru))
      allocate(smrt_DP_tot_yr(nhru))
      allocate(smrt_DP_avg_yr(nhru))
      smrt_DP_tot_mo = 0.
      smrt_DP_avg_mo = 0.
      smrt_DP_tot_yr = 0.
      smrt_DP_avg_yr = 0.

      allocate(smrt_GWSW_MF_tot_mo(MXRIVR))
      allocate(smrt_GWSW_MF_avg_mo(MXRIVR))
      allocate(smrt_GWSW_MF_tot_yr(MXRIVR))
      allocate(smrt_GWSW_MF_avg_yr(MXRIVR))
      smrt_GWSW_MF_tot_mo = 0.
      smrt_GWSW_MF_avg_mo = 0.
      smrt_GWSW_MF_tot_yr = 0.
      smrt_GWSW_MF_avg_yr = 0.

      allocate(smrt_GWSW_SWAT_tot_mo(msub))
      allocate(smrt_GWSW_SWAT_avg_mo(msub))
      allocate(smrt_GWSW_SWAT_tot_yr(msub))
      allocate(smrt_GWSW_SWAT_avg_yr(msub))
      smrt_GWSW_SWAT_tot_mo = 0.
      smrt_GWSW_SWAT_avg_mo = 0.
      smrt_GWSW_SWAT_tot_yr = 0.
      smrt_GWSW_SWAT_avg_yr = 0.

      if(rt_active.eq.1) then
        allocate(smrt_csolute_tot_mo(ncol,nrow,nlay,2))
        allocate(smrt_csolute_avg_mo(ncol,nrow,nlay,2))
        allocate(smrt_csolute_tot_yr(ncol,nrow,nlay,2))
        allocate(smrt_csolute_avg_yr(ncol,nrow,nlay,2))
      endif
      smrt_csolute_tot_mo = 0.
      smrt_csolute_avg_mo = 0.
      smrt_csolute_tot_yr = 0.
      smrt_csolute_avg_yr = 0.

      area_print = 1

      allocate(gw_available(ncol,nrow,nlay))
      gw_available = 0.

      allocate(sub_gw_output(msub))
      sub_gw_output = 0.
      
      return
      end
