      module smrt_parm
!!    ~ ~ ~ Author ~ ~ ~
!!    Tyler Wible, Masters student
!!    Colorado State University 2012-2014
!!    Comment initials "tcw"
!!
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This module contains variables used for linking SWAT variables to MODFLOW and 
!!    RT3D as well as MODFLOW variables to SWAT and RT3D and UZF-RT3D to SWAT
      
      !flag for whether modflow is active -----------------------------------------------------------------------------
      integer :: mf_active

      !Declare global variables for SWAT-MODFLOW-RT3D linkage
      real, dimension (:,:), allocatable :: d2h_id, d2h_area
      integer :: d2h_size
      real, dimension (:,:), allocatable :: grid2riv_id, grid2riv_len
      integer, dimension (:), allocatable :: riv_nsubs,cell_dhrus

      !HRU --> DHRU mapping
      type hru_dhru_mapping
        integer, allocatable :: dhru_id(:)
        real, allocatable :: dhru_perc(:)
      endtype hru_dhru_mapping
      type (hru_dhru_mapping), dimension(:), allocatable :: h2d_map

      !DHRU --> Grid mapping
      type dhru_grid_mapping
        integer, allocatable :: dhru_id(:)
        real, allocatable :: dhru_perc(:)
      endtype dhru_grid_mapping
      type (dhru_grid_mapping), dimension(:), allocatable :: d2g_map

      !grid --> DHRU mapping
      type grid_dhru_mapping
        integer, allocatable :: cell_row(:)
        integer, allocatable :: cell_col(:)
        real, allocatable :: cell_perc(:)
      endtype grid_dhru_mapping
      type (grid_dhru_mapping), dimension(:), allocatable  :: g2d_map

      !Specific Yield variables
      real, dimension (:,:,:), allocatable :: mf_SC2,mf_SC1
      integer :: mf_NTOP
      
      !Variables for irrigation: MODFLOW Well package supplies irrigation water to SWAT HRUs, or vice versa
      integer :: mf_irrigation,mf_irrigation_swat
      real, dimension (:,:), allocatable :: mf_sub_irrigate
      real, dimension (:), allocatable :: mf_sub_irrig_depth,
     &                                    mf_convey_ratio,
     &                                    mf_runoff_ratio
      real, dimension (:,:,:), allocatable :: gw_available
      integer, dimension (:), allocatable :: mf_sub_irrig_flag
      integer, dimension (:), allocatable :: mf_hru_irrig_flag,
     &                                       mf_nsub_hru_irrig
      integer, dimension (:,:), allocatable :: mf_sub_hru_irrig
      integer :: ncell_irrigate
      real :: mf_hru_irrig
      
      !SWAT-MODFLOW Disaggregated HRU (dhru) variables
      real, dimension (:), allocatable :: etremain_dhru,rchrg_dhru,
     &                                    stsol_rd_dhru
      real, dimension (:), allocatable :: rchrg_no3_conc_dhru, 
     &                                    upflux_no3_dhru,
     &                                    rchrg_p_conc_dhru
      integer :: dhru
      
      !Variables for DRAIN package linkage to SWAT sub-basins
      integer, dimension (:,:), allocatable :: drn_subs !rtb drain
      integer :: ndrn_subs,mf_drain_subs

      !Variable for Tile Drainage linkage to MODFLOW !rtb tile
      real, dimension (:), allocatable :: qtile_hru
      real, dimension (:,:), allocatable :: qtile_cells
      integer, dimension (:), allocatable :: hru_ncells
      integer :: mf_tile_active

!!    MODFLOW/SWAT River segment variables
      integer :: nrivcells_subs
      integer, dimension (:), allocatable :: grid2riv_gridID
      
!!    SWAT variables for smrt_upfluxToSoil as populated by RT3D
      real, dimension (:), allocatable :: upflux_no3
      real, dimension (:,:), allocatable :: vadose_thick

!!    Flags for writing output data (SWAT, MODFLOW)
      integer :: out_SWAT_recharge,out_MF_recharge,
     &           out_SWAT_channel,out_MF_channel,
     &           out_MODFLOW_gwsw,out_SWAT_gwsw

!!    Variables and Arrays for storing SWAT-MODFLOW variable averages (monthly, annual) (rtb avg)
      integer :: day_count_mo,day_count_yr,smrt_month_counter,
     &           swatmf_out_avg
      
      integer, dimension (:), allocatable :: month_days,
     &                                       month_days_leap
      real, dimension (:,:), allocatable :: smrt_RECH_tot_mo,
     &                                      smrt_RECH_avg_mo,
     &                                      smrt_RECH_tot_yr,
     &                                      smrt_RECH_avg_yr
      real, dimension (:,:,:), allocatable :: smrt_HNEW_tot_mo,
     &                                        smrt_HNEW_avg_mo,
     &                                        smrt_HNEW_tot_yr,
     &                                        smrt_HNEW_avg_yr
      real, dimension (:), allocatable :: smrt_DP_tot_mo,
     &                                    smrt_DP_avg_mo,
     &                                    smrt_DP_tot_yr,
     &                                    smrt_DP_avg_yr
      real, dimension (:), allocatable :: smrt_GWSW_MF_tot_mo,
     &                                    smrt_GWSW_MF_avg_mo,
     &                                    smrt_GWSW_MF_tot_yr,
     &                                    smrt_GWSW_MF_avg_yr
      real, dimension (:), allocatable :: smrt_GWSW_SWAT_tot_mo,
     &                                    smrt_GWSW_SWAT_avg_mo,
     &                                    smrt_GWSW_SWAT_tot_yr,
     &                                    smrt_GWSW_SWAT_avg_yr
      real, dimension (:,:,:,:), allocatable :: smrt_csolute_tot_mo,
     &                                          smrt_csolute_avg_mo,
     &                                          smrt_csolute_tot_yr,
     &                                          smrt_csolute_avg_yr

!!    Array for storing MODFLOW observation cells
      integer :: num_MF_obs,mf_obs_flag
      integer, dimension (:,:), allocatable :: MF_obs

!!    Array for storing SWAT-MODFLOW output times
      integer :: n_outswatmf,swatmf_out_ctr
      integer, dimension (:), allocatable :: outswatmf

!!    Array for applied irrigation water from SWAT HRUs, that will be used to populate
!!    pumping rates in the WELL package
      real, dimension (:), allocatable :: irrig_depth,irrig_volume

      !flag for printing the subbasin area in km2
      integer :: area_print
      
      !arrays for subbasin groundwater output (for output.sub file)
      real, dimension (:), allocatable :: sub_gw_output

      end module smrt_parm
