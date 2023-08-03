*     Subroutine:   rt3d_finish
*     Purpose:      finalize output and deallocate memory. stop rt3d.
      
      subroutine rt_close

      use rt_global
      use mf_rt_link
      use GLOBAL,only:NROW,NCOL,NLAY !MODFLOW

*	Calculate the total CPU time for the simulation
	!call cpu_time(time_end)
      !time_total = time_end - time_begin
	!WRITE(5000,*) 'CPU Time = ', time_total

*     Print out restart file (final concentrations, water content)  
      call btn_restart

*     Write out termination message      
      write(IOUT,1200)
      write(IOUT,1225)
      write(IOUT,1200)
	write(*,1200)
      write(*,1225)
	write(*,*) 'PROGRAM TERMINATED SUCCESSFULLY'
 
      stop !rt3d is now stopped
 

 1199 format(1000(f12.4)) 
 1200 format(1X,' ----- ')
 1225 format(1X,'| R T |'
     &      /1X,'| 3 D | END OF MODEL OUTPUT')
	
      return
      end
