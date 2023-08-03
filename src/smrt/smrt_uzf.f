      subroutine smrt_uzf

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the necessary SWAT variables for the MODFLOW
!!    unsaturated zone flow (UZF1) package's variables
!!      
        use parm, only: day_total
        use smrt_parm
        use GLOBAL, only: NCOL,NROW
        use GWFUZFMODULE, only: FINF,PETRATE
        implicit none 

        !Define local variables
        integer i,j

        call smrt_dhru2grid2D(rchrg_dhru,FINF,1)
        call smrt_dhru2grid2D(etremain_dhru,PETRATE,1)

        !output Recharge values to a file
        if(out_MF_recharge.eq.1 .and.
     &    day_total.eq.outswatmf(swatmf_out_ctr)) then
          write(30002,*)
          write(30002,*) 'Day:',day_total
          write(30002,*) 'daily infiltration values provided to UZF'
          do i=1,nrow
            write(30002,100) (FINF(j,i),j=1,ncol)
          enddo
        endif

  100 format(1000(e17.10))

      end subroutine 
