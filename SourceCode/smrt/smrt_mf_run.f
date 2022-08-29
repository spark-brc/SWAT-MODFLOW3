      subroutine smrt_mf_run
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine calls conversions from SWAT to MODFLOW, then runs MODFLOW
!!    for one day, then converts back the groundwater results from MODFLOW to
!!    SWAT
!!      
        use parm, only: etremain,percn,percp,sepbtm,nhru, !SWAT
     &                  leapyr,iida,iyr,msub
        use GLOBAL, only: mf_interval,mf_ran,NCOL,NROW,NLAY,IBOUND !MODFLOW
        use GWFRIVMODULE, only:MXRIVR,RIVR
        use GWFBASMODULE, only:HNOFLO 
        use mf_rt_link, only: rt_active !MODFLOW-RT3D Linkage
        use smrt_parm

        implicit none
        
        integer i,j,k,IL,IC,IR,m,mf_gridID
        integer days_in_month(12) !rtb avg
        
        
        !Convert SWAT units into MODFLOW units
        call smrt_conversion2mf
          
        !run modflow (water table, river discharge/seepage, pumping)
        call mf_run
        mf_ran = .true.
          
        !Convert MODFLOW units back to SWAT units
        if(mf_ran) then
          call smrt_conversion2swat
        endif


        !Print out SWAT-MODFLOW variable averages (rtb avg) -----------------------------------------------------------
        if(swatmf_out_avg) then
        if(leapyr.eq.1) then
          days_in_month = month_days
        else
          days_in_month = month_days_leap
        endif
        
        !check to see if end of month has been reached
        if(iida.eq.days_in_month(smrt_month_counter)) then
          
          !print out monthly averages for MODFLOW recharge
          write(30020,*), 'month:',smrt_month_counter,'year:',iyr
          do i=1,NROW
            write(30020,100) (smrt_RECH_avg_mo(j,i),j=1,NCOL)
          enddo
          write(30020,*)
          
          !print out monthly averages for SWAT deep percolation
          write(30024,*), 'month:',smrt_month_counter,'year:',iyr
          do i=1,nhru
            write(30024,*) smrt_DP_avg_mo(i)
          enddo
          write(30024,*)

          !MODFLOW GW/SW
          write(30026,*), 'month:',smrt_month_counter,'year:',iyr
          write(30026,*) 'Layer, Row, Column, Rate'
          do i=1,MXRIVR !spark - added Layer, row, column info
            !determine if the River cell is linked with a SWAT sub-basin
            IR = RIVR(2,i) !row of MODFLOW river cell
            IC = RIVR(3,i) !column of MODFLOW river cell
            mf_gridID = ((IR-1)*NCOL) + IC !grid ID of MODFLOW river cell
            do m=1,nrivcells_subs
              !only proceed if there is a match
              if(mf_gridID.eq.grid2riv_gridID(m)) then
                IL=RIVR(1,i)  
                IR=RIVR(2,i)
                IC=RIVR(3,i)
                write(30026,*) IL, IR, IC, smrt_GWSW_MF_avg_mo(i)
              end if
            end do
          enddo !spark
          write(30026,*)

          !SWAT GW/SW
          write(30028,*), 'month:',smrt_month_counter,'year:',iyr
          do i=1,msub
            write(30028,*) i,smrt_GWSW_SWAT_avg_mo(i)
          enddo
          write(30028,*)

          !solute concentration
          if(rt_active.eq.1) then
          write(30030,*), 'month:',smrt_month_counter,'year:',iyr
          do k=1,NLAY
            write(30030,*),'layer:',k
            do i=1,NROW
              write(30030,100) (smrt_csolute_avg_mo(j,i,k,1),j=1,NCOL)
            enddo
          enddo
          write(30030,*)
          write(30032,*), 'month:',smrt_month_counter,'year:',iyr
          do k=1,NLAY
            write(30032,*),'layer:',k
            do i=1,NROW
              write(30032,100) (smrt_csolute_avg_mo(j,i,k,2),j=1,NCOL)
            enddo
          enddo
          write(30032,*)
          endif

          !print out monthly averages for groundwater head
          write(30022,*), 'month:',smrt_month_counter,'year:',iyr
          do k=1,NLAY
            write(30022,*),'layer:',k
            do i=1,NROW
              write(30022,101) (smrt_HNEW_avg_mo(j,i,k),j=1,NCOL)
            enddo
          enddo
          write(30022,*)

          !reset total monthly values
          smrt_RECH_tot_mo = 0.
          smrt_DP_tot_mo = 0.
          smrt_GWSW_MF_tot_mo = 0.
          smrt_GWSW_SWAT_tot_mo = 0.
          smrt_csolute_tot_mo = 0.
          do k=1,nlay
            do j=1,ncol
              do i=1,nrow
                if(ibound(j,i,k).eq.1) then
                  smrt_HNEW_tot_mo = 0.
                endif
              enddo
            enddo
          enddo

          !print out yearly averages if last month of year
          if(smrt_month_counter.eq.12) then
            
            !MODFLOW recharge
            write(30021,*), 'year:',iyr
            do i=1,NROW
              write(30021,100) (smrt_RECH_avg_yr(j,i),j=1,NCOL)
            enddo
            write(30021,*)

            !SWAT deep percolation
            write(30025,*), 'year:',iyr
            do i=1,nhru
              write(30025,*) smrt_DP_avg_yr(i)
            enddo
            write(30025,*)

            !MODFLOW GW/SW
            write(30027,*), 'year:',iyr
            write(30027,*) 'Layer, Row, Column, Rate'
            do i=1,MXRIVR !spark - added Layer, row, column info
              !determine if the River cell is linked with a SWAT sub-basin
              IR = RIVR(2,i) !row of MODFLOW river cell
              IC = RIVR(3,i) !column of MODFLOW river cell
              mf_gridID = ((IR-1)*NCOL) + IC !grid ID of MODFLOW river cell
              do m=1,nrivcells_subs
              !only proceed if there is a match
                if(mf_gridID.eq.grid2riv_gridID(m)) then
                  IL=RIVR(1,i)  
                  IR=RIVR(2,i)
                  IC=RIVR(3,i)
                  write(30027,*) IL, IR, IC, smrt_GWSW_MF_avg_yr(i)
                endif
              enddo
            enddo
            write(30027,*)

            !SWAT GW/SW
            write(30029,*), 'year:',iyr
            do i=1,msub
              write(30029,*) smrt_GWSW_SWAT_avg_yr(i)
            enddo
            write(30029,*)

            !head
            write(30023,*), 'year:',iyr
            do k=1,NLAY
              write(30023,*),'layer:',k
              do i=1,NROW
                write(30023,101) (smrt_HNEW_avg_yr(j,i,k),j=1,NCOL)
              enddo
            enddo
            write(30023,*)

            !solute concentration
            if(rt_active.eq.1) then
            write(30031,*), 'year:',iyr
            do k=1,NLAY
              write(30031,*),'layer:',k
              do i=1,NROW
                write(30031,100) (smrt_csolute_avg_yr(j,i,k,1),j=1,NCOL)
              enddo
            enddo
            write(30031,*)
            write(30033,*), 'year:',iyr
            do k=1,NLAY
              write(30033,*),'layer:',k
              do i=1,NROW
                write(30033,100) (smrt_csolute_avg_yr(j,i,k,2),j=1,NCOL)
              enddo
            enddo
            write(30033,*)
            endif

            !reset total annual values
            smrt_RECH_tot_yr = 0.
            smrt_DP_tot_yr = 0.
            smrt_GWSW_MF_tot_yr = 0.
            smrt_GWSW_SWAT_tot_yr = 0.
            smrt_csolute_tot_yr = 0.
            do k=1,nlay
              do j=1,ncol
                do i=1,nrow
                  if(ibound(j,i,k).eq.1) then
                    smrt_HNEW_tot_yr = 0.
                  endif
                enddo
              enddo
            enddo

            day_count_yr = 0
            smrt_month_counter = 0 !start in January
          endif
          smrt_month_counter = smrt_month_counter + 1
          day_count_mo = 0 !reset
        endif

        !increment each counter by one day
        day_count_mo = day_count_mo + 1
        day_count_yr = day_count_yr + 1
      endif

  100 format(1000(f12.4))
  101 format(1000(1pe13.5)) !spark modified format for high altitude

      end subroutine smrt_mf_run
