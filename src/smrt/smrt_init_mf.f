      subroutine smrt_init_mf

!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This subroutine initializes MODFLOW linking subroutines of SWAT-MODFLOW-RT3D

      use GLOBAL,only:IUNIT !rtb drain
      use parm, only:nhru,gw_delaye
      use smrt_parm
      use mf_rt_link
      implicit none
      
      integer i,j,k,n,dhruID,hru_read
      real    delay

      print *, 'MODFLOW is being used' !rtb
        
      !set up MODFLOW data and allocate arrays
      call mf_read !rtb
      write(6008,*) 'swatmf_init: MODFLOW files have been read'

      !read SWAT-MODFLOW linkage files
      call smrt_read_grid2dhru
      call smrt_read_dhru2grid
      call smrt_read_dhru2hru
      call smrt_read_river2grid
        
      !call for pumping-irrigation linkage
      if(mf_irrigation.eq.1 .or.
     &   mf_irrigation_swat.eq.1) then
        call smrt_read_irrig
      endif

      !call DRAIN cells linkage if the MODFLOW DRAIN package is active
      if(IUNIT(3).GT.0 .and. mf_drain_subs.eq.1) then !rtb drain
        call smrt_read_drain2sub              
      endif

      !read groundwater delay values
      read(6001,*)
      read(6001,*) hru_read
      if(hru_read.eq.0) then
        read(6001,*) delay
        gw_delaye = Exp(-1./(delay + 1.e-6))
      else
        do n=1,nhru
          read(6001,*) delay
          gw_delaye(n) = Exp(-1./(delay + 1.e-6))
        enddo
      endif
      close(6001)

      !Additional calculations required for later subroutines.
      !determine the number of MODFLOW grid cells in each HRU
      allocate(hru_ncells(nhru))
      hru_ncells = 0
      do n=1,nhru !loop through the HRUs
        do j=1,size(h2d_map(n)%dhru_id) 
            
          !ID of the current DHRU
          dhruID = h2d_map(n)%dhru_id(j)

          !loop through the MODFLOW cells that intersect the DHRU
          do k=1,size(g2d_map(dhruID)%cell_row)
            hru_ncells(n) = hru_ncells(n) + 1
          enddo
        enddo
      enddo !go to the next HRU

      !read RT3D information if active
      if(rt_active) call smrt_init_rt3d

      !write out information to log file
      write(6008,*) 'swatmf_init: initialization finished'


      return
      end subroutine smrt_init_mf
