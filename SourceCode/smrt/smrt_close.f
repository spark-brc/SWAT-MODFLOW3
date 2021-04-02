      subroutine smrt_close

!!    ~ ~ ~ Purpose ~ ~ ~
!!    This subroutine deallocates the variables from smrt_parm and some
!!    additional MODFLOW and RT3D variables which had their deallocate calls moved
      
!     Import variables
      use GWFRIVMODULE, only: RIVAUX, RIVR !MODFLOW
      use smrt_parm !smrt linkage
      implicit none
      
!     Deallocate MODFLOW variables
      deallocate(RIVAUX)
      deallocate(RIVR)
      
!     Deallocate smrt variables
      deallocate(g2d_map)
      deallocate(d2g_map)
      deallocate(h2d_map)
      deallocate(grid2riv_id)
      deallocate(grid2riv_len)
      deallocate(riv_nsubs)
      deallocate(upflux_no3)
      
      deallocate(etremain_dhru)
      deallocate(stsol_rd_dhru)
      deallocate(rchrg_dhru)
      deallocate(rchrg_no3_conc_dhru)
      deallocate(rchrg_p_conc_dhru)
      deallocate(upflux_no3_dhru)
      
      return
      end