      subroutine smrt_evt

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine converts the SWAT remaining ET (PET - Actual ET) to rates for MODFLOW's EVT package

        !Import variables
        use parm, only: day_total,nhru
        use GLOBAL, only: NCOL,NROW,BOTM
        use GWFEVTMODULE, ONLY:EVTR,EXDP,SURF
        use smrt_parm
        implicit none
        
        !Define local variables
        integer i,j

        !MODFLOW variables
        !EVTR: maximum ET rate
        !EXDP: extinction depth (depth at which ET ceases)
        !SURF: surface from which extinction depth is measured

        !SWAT variables
        !stsol_rd: current rooting depth (mm)
        !etremain: remaining ET (potential ET - actual ET) (mm)

        !map DHRU etremain values to MODFLOW grid cells (and multiply by cell area to get flow rate)        
        call smrt_dhru2grid2D(etremain_dhru,EVTR,1)

        !Extinction depth is set in the *.evt MODFLOW input file
        !EXDP

        !store ground surface elevation in SURF
        do i=1,NROW
          do j=1,NCOL
            SURF(j,i) = BOTM(j,i,0)  
          enddo
        enddo

      end subroutine smrt_evt
