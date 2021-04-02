*     This subroutine runs the RT3D simulator for the time period indicated. Typically,
*     the stress period and flow time step are the same as for SWAT (1 day).

*     This subroutine is called for each flow time step of MODFLOW

      subroutine rt_run

      use rt_global
      use mf_rt_link
      use smrt_parm
      use GLOBAL,only:ITMUNI,NROW,NCOL,NLAY,MF_leapyr,HNEW !MODFLOW
      use GWFBASMODULE,only:HDRY

      implicit none
      integer   i,j,k,N,L,dum
      real      time_update
      common   /GCGIDX/L(19) !for implicit shceme
           

      !print a message to the screen, indicating that RT3D is running
      print *, '            RT3D is running'

*     This subroutine is called within the MODFLOW flow time step loop. Thus, it
*     calculates the reactive transport of chemical species for 1 day (to match
*     the SWAT step size)
      
      !Calculate the dispersion coefficients for the current flow time step
      if(TRNOP(2)) call dsp_coeff

      dt = rt_delt !(calculated in mf_run)
      t1 = t2
      t2 = t2 + dt
      
      !zero out mass transfer arrays
      rt_rivmass = 0.
      rt_drnmass = 0.
      
      !loop through the transport time steps for the flow time step
      do ntrans=1,5000

        !advance once transport step
        call btn_ts

        if(time2.ge.4233600) then
          dum = 10
        endif

        !for each mobile species, compute change in concentration due to advection, 
        !dispersion, and source-sink mixing. this is done implicitly
        do icomp=1,ncomp
          if(phase(icomp).eq.1) then !species is mobile
            call rt_implicit
          endif
        enddo
          
        !compute change in concentration due to chemical reactions
        if(trnop(4)) call rct_solve

        !check for NAN values
        do k=1,nlay
          do i=1,nrow
            do j=1,ncol
              do icomp=1,ncomp
                if(isnan(cnew(j,i,k,icomp))) then
                  dum = 10
                endif
              enddo
            enddo
          enddo
        enddo

                 
        !compute mass budgets and output results for the transport step
        do icomp=1,ncomp
          call btn_budget !mass budgets
          call btn_output !save outputs   
        enddo 

        !calculate portion of flow time step accomplished
        time_update = time2 / t2

        print *, '                    Transport step',ntrans,time_update

        !if end of flow time step, exit
        if(time2.ge.t2) goto 100

      enddo !go to next transport time step


      !store concentration values for average outputs
100   call units(dt, ITMUNI, 4, 1, 1, MF_leapyr) !tcw
      sum_tracker = sum_tracker + dt !in days
      
      !calculate average concentration for the current month (rtb avg)
      smrt_csolute_tot_mo = smrt_csolute_tot_mo + cnew
      smrt_csolute_avg_mo = smrt_csolute_tot_mo / day_count_mo

      !calculate average recharge for the current year (rtb avg)
      smrt_csolute_tot_yr = smrt_csolute_tot_yr + cnew
      smrt_csolute_avg_yr = smrt_csolute_tot_yr / day_count_yr
      
      do k=1,nlay
        do j=1,ncol
          do i=1,nrow
            if(ICBUND(j,i,k,1).eq.0 .or. HNEW(j,i,k).eq.HDRY) then
              do icomp=1,ncomp
                smrt_csolute_avg_mo(j,i,k,icomp) = 0.
                smrt_csolute_avg_yr(j,i,k,icomp) = 0.
              enddo
            endif
          enddo
        enddo
      enddo

 200  format(1000(e12.4))

      return
	end
