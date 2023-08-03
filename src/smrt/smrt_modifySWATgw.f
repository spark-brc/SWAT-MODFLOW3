      subroutine smrt_modifySWATgw

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine calls conversions from SWAT to MODFLOW, then runs MODFLOW
!!    for one day, then converts back the groundwater results from MODFLOW to
!!    SWAT
!!      
        use parm, only: deepst,gw_q,gw_qdeep,hru_sub,ihru,rttlc,shallst,
     &                  gwseep !SWAT
        use smrt_parm !sm linkage
        implicit none
        
        !define local variables
        integer i,j
        integer dhruID,nsubs,subID

        !zero out groundwater variables (groundwater is handled by MODFLOW)
        do i=1,size(h2d_map(ihru)%dhru_id)
          dhruID = h2d_map(ihru)%dhru_id(i)
          if(dhruID.ne.0)then
            gw_q(ihru) = 0
            shallst(ihru) = 0 !also zero the shallow groundwater storage because this is now held in MODFLOW
            gw_qdeep(ihru) = 0 !also zero the deep groundwater flow to streams because MODFLOW only passes back 1 baseflow value
            deepst(ihru) = 0 !also zero the deep groundwater storage because this is now held in MODFLOW
            gwseep = 0.
            exit
          endif
        enddo
        
        !zero out reach loss (stream-groundwater is handled by MODFLOW)
        do i=1,nrivcells_subs
          nsubs = riv_nsubs(i)
          do j=1,nsubs
            subID = grid2riv_id(i,j)
            if(subID.eq.hru_sub(ihru))then
              rttlc(subID) = 0. !also zero the reach losses because this is calculated in MODFLOW and passed directly into varoute
              exit
            endif
          enddo
        enddo

      end subroutine smrt_modifySWATgw
