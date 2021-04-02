      subroutine smrt_read_link

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine reads in the file containing the information to link SWAT and MODFLOW
!!
!!    ~ ~ ~ Variables Modified ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    mf_active    |none          |index whether or not to use MODFLOW to
!!                 |              |calculate groundwater flow processes
!!                 |              |instead of SWAT's gwmod
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

!     Import variables
      use parm, only: nhru,gw_delaye
      use GLOBAL, only: mf_ran,mf_interval !MODFLOW
      use mf_rt_link, only: rt_active !MODFLOW-RT3D Linkage
      use smrt_parm !SWAT-MODFLOW linkage
      implicit none
      
      integer n


      !Read in flags for linking SWAT and MODFLOW -----------------------------------------------------------
      open(6001,file='swatmf_link.txt')
      read(6001,*) mf_active
      read(6001,*) mf_irrigation
      read(6001,*) mf_irrigation_swat
      read(6001,*) mf_drain_subs
      read(6001,*) rt_active
      mf_ran = .false.
      mf_interval = 1 !this is set to 1 day (constant for all SWAT-MODFLOW simulations)
      read(6001,*) mf_obs_flag

      !open up log file
      if(mf_active.eq.1) then
        open(6008,file='swatmf_log')
        write(6008,*) 'Progress of SWAT-MODFLOW simulation'
        write(6008,*)
      endif

      !write information to log file
      write(6008,*) 'swatmf_link.txt:    file flags have been read'
      if(mf_active.eq.1) then
        write(6008,*) 'swatmf_link.txt:    MODFLOW is active'
      endif
      if(mf_irrigation.eq.1) then
        write(6008,*) 'swatmf_link.txt:    MODFLOW irrigation is active'
      endif
      if(mf_irrigation_swat.eq.1) then
        write(6008,*) 'swatmf_link.txt:    SWAT irrigation is active'
      endif
      if(mf_drain_subs.eq.1) then
        write(6008,*) 'swatmf_link.txt:    DRAIN cells are active'
      endif
      if(rt_active.eq.1) then
        write(6008,*) 'swatmf_link.txt:    RT3D (N,P) is active'
      endif
      if(mf_obs_flag.eq.1) then
        write(6008,*) 'swatmf_link.txt:    Observation cell file needed'
      endif


      !Optional output files (for SWAT-MODFLOW information) -------------------------------------------------
        
      !SWAT Recharge to Water Table (by SWAT HRU)
      read(6001,*)
      read(6001,*) out_SWAT_recharge 
      if(out_SWAT_recharge.eq.1) then
        open(30001,file='swatmf_out_SWAT_recharge') !rtb
        write(30001,*)'SWAT recharge to water table (mm) (for each HRU)'
      endif

      !MODFLOW Recharge (by MODFLOW cell)
      read(6001,*) out_MF_recharge 
      if(out_MF_recharge.eq.1) then
        open(30002,file='swatmf_out_MF_recharge') !rtb
        write(30002,*) 'MODFLOW Recharge (L3/T) (for each cell)'
        write(30002,*) '--Calculated from SWAT HRU recharge--'
      endif

      !SWAT Channel Depth (by SWAT subbasin)
      read(6001,*) out_SWAT_channel 
      if(out_SWAT_channel.eq.1) then
        open(30003,file='swatmf_out_SWAT_channel') !rtb
        write(30003,*) 'SWAT channel depth (m) (for each subbasin)'
      endif

      !MODFLOW River Stage (by MODFLOW River Cell)
      read(6001,*) out_MF_channel
      if(out_MF_channel.eq.1) then
        open(30004,file='swatmf_out_MF_riverstage') !rtb
        write(30004,*) 'MODFLOW River Stage (L) (for each River Cell)'
        write(30004,*) '--Calculated from SWAT Channel Depth--'
      endif

      !GW/SW exchange (by MODFLOW River cell)
      read(6001,*) out_MODFLOW_gwsw 
      if(out_MODFLOW_gwsw.eq.1) then
        open(30005,file='swatmf_out_MF_gwsw') !rtb
        write(30005,*) 'Groundwater/Surface Water exchange (L3/T)'
        write(30005,*) 'for each MODFLOW River Cell'
      write(30005,*) 'Positive: River water seeps to the aquifer'
      write(30005,*) 'Negative: Groundwater flows from aquifer to river'
      endif

      !GW/SW exchange (by SWAT subbasin)
      read(6001,*) out_SWAT_gwsw 
      if(out_SWAT_gwsw.eq.1) then
        open(30006,file='swatmf_out_SWAT_gwsw') !rtb
        write(30006,*) 'Groundwater/Surface Water exchange (m3/day)'
        write(30006,*) 'for each SWAT Subbasin'
      write(30006,*) 'Positive: Volume entering stream from the aquifer'
      write(30006,*) 'Negative: Volume seeps from stream to the aquifer'
      endif

      !read flag for printing out monthly and annual average output
      read(6001,*) swatmf_out_avg

      write(6008,*) 'swatmf_link.txt:    output flags have been read'

      !output if RT3D is active
      if(rt_active.eq.1) then
      
      !GW/SW NO3 exchange (by MODFLOW River Cell) --------------------------------------- NITRATE
      open(30007,file='swatmf_out_RT_rivno3') !rtb
      write(30007,*) 'Groundwater/Surface Water NO3 exchange (kg/day)'
      write(30007,*) 'for each MODFLOW River Cell'
      write(30007,*) 'Positive: Mass seeps to the aquifer'
      write(30007,*) 'Negative: Mass from aquifer to river'

      !GW/SW NO3 exchange (by SWAT subbasin)
      open(30008,file='swatmf_out_SWAT_rivno3') !rtb
      write(30008,*) 'GW/SW and Drain NO3 exchange (kg/day)'
      write(30008,*) 'for each SWAT Subbasin'
      write(30008,*) 'Positive: Mass entering stream from the aquifer'
      write(30008,*) 'Negative: Mass seeps from stream to the aquifer'
      
      !NO3 recharge concentration to MODFLOW grid
      open(30009,file='swatmf_out_RT_rechno3') !rtb
      write(30009,*) 'RT3D NO3-N Recharge Conc. (mg/L) for each cell'
      write(30009,*) '--Calculated from SWAT HRU recharge--'

      !NO3 recharge concentration from each SWAT HRU
      open(30010,file='swatmf_out_SWAT_rechno3') !rtb
      write(30010,*) 'SWAT NO3-N Recharge Conc. (mg/L)'

      !GW/SW P exchange (by MODFLOW River Cell) ----------------------------------------- PHOSPHORUS
      open(30011,file='swatmf_out_RT_rivP') !rtb
      write(30011,*) 'GW/SW and Drain P exchange (kg/day)'
      write(30011,*) 'for each MODFLOW River Cell'
      write(30011,*) 'Positive: Mass seeps to the aquifer'
      write(30011,*) 'Negative: Mass from aquifer to river'

      !GW/SW P exchange (by SWAT subbasin)
      open(30012,file='swatmf_out_SWAT_rivP') !rtb
      write(30012,*) 'Groundwater/Surface Water P exchange (kg/day)'
      write(30012,*) 'for each SWAT Subbasin'
      write(30012,*) 'Positive: Mass entering stream from the aquifer'
      write(30012,*) 'Negative: Mass seeps from stream to the aquifer'
      
      !P recharge concentration to MODFLOW grid
      open(30013,file='swatmf_out_RT_rechP') !rtb
      write(30013,*) 'RT3D P Recharge Conc. (mg/L) for each cell'
      write(30013,*) '--Calculated from SWAT HRU recharge--'

      !P recharge concentration from each SWAT HRU
      open(30014,file='swatmf_out_SWAT_rechP') !rtb
      write(30014,*) 'SWAT Dissolved P Recharge Conc. (mg/L)'

      endif


      !Read output control for SWAT-MODFLOW variables -------------------------------------------------------
      read(6001,*)
      read(6001,*) n_outswatmf
      allocate(outswatmf(n_outswatmf))
      do n=1,n_outswatmf
        read(6001,*) outswatmf(n)
      enddo
      swatmf_out_ctr = 1
      write(6008,*) 'swatmf_link.txt:    output control has been read'

      !open files for variable average output (rtb avg) -----------------------------------------------------
      if(swatmf_out_avg) then
        
        !recharge (modflow)
        open(30020,file='swatmf_out_MF_recharge_monthly')
        write(30020,*) 'Monthly Averaged Recharge Values for MODFLOW'
        write(30020,*)

        open(30021,file='swatmf_out_MF_recharge_yearly')
        write(30021,*) 'Yearly Averaged Recharge Values for MODFLOW'
        write(30021,*)

        !groundwater head
        open(30022,file='swatmf_out_MF_head_monthly')
        write(30022,*) 'Monthly Averaged Head Values from MODFLOW'
        write(30022,*)

        open(30023,file='swatmf_out_MF_head_yearly')
        write(30023,*) 'Yearly Averaged Head Values from MODFLOW'
        write(30023,*)

        !recharge to water table (swat)
        open(30024,file='swatmf_out_SWAT_recharge_monthly')
        write(30024,*) 'Monthly Averaged Recharge Values from SWAT'
        write(30024,*)

        open(30025,file='swatmf_out_SWAT_recharge_yearly')
        write(30025,*) 'Yearly Averaged Recharge Values from SWAT'
        write(30025,*)

        !gw/sw exchange (modflow)
        open(30026,file='swatmf_out_MF_gwsw_monthly')
        write(30026,*) 'Monthly Averaged GW/SW Rates for MODFLOW'
        write(30026,*)

        open(30027,file='swatmf_out_MF_gwsw_yearly')
        write(30027,*) 'Yearly Averaged GW/SW Rates for MODFLOW'
        write(30027,*)

        !gw/sw exchange (swat)
        open(30028,file='swatmf_out_SWAT_gwsw_monthly')
        write(30028,*) 'Monthly Averaged GW/SW Rates for SWAT'
        write(30028,*)

        open(30029,file='swatmf_out_SWAT_gwsw_yearly')
        write(30029,*) 'Yearly Averaged GW/SW Rates for SWAT'
        write(30029,*)

        !solute concentration in groundwater (RT3D)
        if(rt_active.eq.1) then
          open(30030,file='swatmf_out_RT_cno3_monthly')
          write(30030,*) 'Monthly Averaged GW Nitrate Concentration'
          write(30030,*)

          open(30031,file='swatmf_out_RT_cno3_yearly')
          write(30031,*) 'Yearly Averaged GW Nitrate Concentration'
          write(30031,*)
          
          open(30032,file='swatmf_out_RT_cp_monthly')
          write(30032,*) 'Monthly Averaged GW Phosphorus Concentration'
          write(30032,*)

          open(30033,file='swatmf_out_RT_cp_yearly')
          write(30033,*) 'Yearly Averaged GW Phosphorus Concentration'
          write(30033,*)    
        endif

      endif

      write(6008,*)


      return
      end
