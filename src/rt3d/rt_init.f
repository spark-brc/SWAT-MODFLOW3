*     Subroutine:   rt3d_init
*     Purpose:      Read and prepare data (allocate) required for the rt3d simulation
      
      subroutine rt_init

      use rt_global
      use mf_rt_link
      use GLOBAL,only:NROW,NCOL,NLAY !MODFLOW

*     Calculate time at beginning of simulation
	call cpu_time(time_begin)

*     Open up the input and output files
	call rt_openio

*     Write program title to output file
      write(IOUT,*) '**************************************************'
      write(IOUT,*) '                    UZF-RT3D'
      write(IOUT,*) '    Variably-Saturated Reactive Transport of' 
      write(IOUT,*) '     Multiple Chemical Species in Soil and'
      write(IOUT,*) '              Groundwater Systems'
      write(IOUT,*) '**************************************************'
      write(IOUT,*)

*     Define problem dimension and simulation options
      call btn_define

*     Allocate storage space for data arrays and read file input
      call btn_alloc !grid and concentration arrays
      call fmi_alloc !head and flow data
      if(TRNOP(1)) call adv_alloc !advection
      if(TRNOP(2)) call dsp_alloc !dispersion
      if(TRNOP(3)) call ssm_alloc !source-sink mixing
      if(TRNOP(4)) call rct_alloc !chemical reactions
      if(TRNOP(5)) call gcg_alloc !implicit solver (GCG)
      if(TRNOP(7)) call vst_alloc !variably-saturated transport
      call btn_read !basic transport file

*     Prepare time-stepping information      
      t1 = 0. !beginning time of the flow time step
      t2 = 0. !ending time of the flow time step
      nps = 1
      kflow = 1 !cumulative flow time step

      return
      end
