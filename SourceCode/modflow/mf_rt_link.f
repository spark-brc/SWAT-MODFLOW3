      module mf_rt_link
!!    ~ ~ ~ Author ~ ~ ~
!!    Tyler Wible, Masters student
!!    Colorado State University 2012-2014
!!    Comment initials "tcw"
!!
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    This module contains variables used for linking MODFLOW and 
!!    UZF-RT3D
!!      
!C     Data for linking between MODFLOW and UZF-RT3D ---------------------------------------------------

        !Flags for source-sink packages (turned on/off by MODFLOW)
        logical::FWEL,FDRN,FRCH,FEVT,FRIV,FGHB,FCNH

        !Size of flow and head arrays
        integer::maxcnh,maxwel,maxdrn,maxriv,maxghb,
     &           rt_ncnh,rt_nwel,rt_ndrn,rt_nriv,rt_nghb
        
        !Arrays populated by MODFLOW simulation results (flows and heads)
        real, allocatable, save :: QX(:,:,:),QY(:,:,:),QZ(:,:,:)
        real, allocatable, save :: QSTO(:,:,:),DH(:,:,:)
        real, allocatable, save :: MC(:,:,:),UZQ(:,:,:),UZQSTO(:,:,:),
     &                             GWQOUT(:,:,:),RT_GWET(:,:,:)
        integer, allocatable, save :: RT_IRCH(:,:),RT_IEVT(:,:)
        real, allocatable, save :: CNH(:,:),WEL(:,:),DRN(:,:),RIV(:,:),
     &                             GHB(:,:)
        real, allocatable, save :: RT_RECH(:,:),RT_EVTR(:,:)
        
        !Arrays populated by SWAT simulation results
        !solute concentration in river cells, solute concentration in deep
        !percolation water.
        real, allocatable, save :: rt_criv(:,:),crch(:,:,:)

        !Array populated by RT3D for SWAT simulation
        real, allocatable, save :: rt_rivmass(:,:),rt_drnmass(:,:)

        !flag for MODFLOW model
        integer :: rt_active !rtb
        
        !stress period, flow time step, time step
        integer :: rt_kper,rt_kstp,rt_delt

      endmodule mf_rt_link
