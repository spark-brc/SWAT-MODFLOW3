      subroutine smrt_recharge

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the necessary SWAT variables for the MODFLOW recharge
!!    (RCH) package's variables
!!      
        !Import variables
        use parm, only:day_total
        use GLOBAL, only: NCOL,NROW
        use GWFRCHMODULE, only: RECH
        use smrt_parm
        implicit none
        
        !Define local variables
        integer i,j


        !convert to MODFLOW Recharge array
        call smrt_dhru2grid2D(rchrg_dhru,RECH,1)

        !calculate average recharge for the current month (rtb avg)
        smrt_RECH_tot_mo = smrt_RECH_tot_mo + RECH
        smrt_RECH_avg_mo = smrt_RECH_tot_mo / day_count_mo

        !calculate average recharge for the current year (rtb avg)
        smrt_RECH_tot_yr = smrt_RECH_tot_yr + RECH
        smrt_RECH_avg_yr = smrt_RECH_tot_yr / day_count_yr

        !output Recharge values to a file
        if(out_MF_recharge.eq.1 .and.
     &    day_total.eq.outswatmf(swatmf_out_ctr)) then
          write(30002,*)
          write(30002,*) 'Day:',day_total
          write(30002,*) 'daily recharge values provided to MODFLOW'
          do i=1,nrow
            write(30002,100) (RECH(j,i),j=1,ncol)
          enddo
        endif

  100 format(1000(e17.10))

      end subroutine smrt_recharge 
