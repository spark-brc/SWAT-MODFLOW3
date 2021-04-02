*     Subroutine:   rt3d_stress
*     Purpose:      read and prepare information for the current stress period
      
      subroutine rt_stress

      use rt_global
      use mf_rt_link


*     Write out timing information to the screen and output file      
      write(IOUT,*)
      write(IOUT,*)
      write(IOUT,*)
      write(IOUT,51) 
      write(IOUT,52) rt_kper

*     Read and prepare Source-Sink Mixing input
      if(TRNOP(3)) call ssm_read


 50   format(/1X,'STRESS PERIOD NO.',I5)
 51   format(120('*'))
 52   format(15X,'STRESS PERIOD NO.',I8,3X,'OF',I8)
      
      return
      end