C ************************************************************************************************************
C RXNS.f
C This file holds subroutines that calculate the change in concentration of simulation species due to
C kinetic chemical reactions (mass-balance ordinary differential equations). 
C ************************************************************************************************************ 

C ************************************************************************************************************
C RXNS
C Mass balance ordinary differential equations (ODE) for user-defined kinetic chemical reactions 
C ************************************************************************************************************ 
      subroutine rxns(J,I,K,CONC_RG,K_RG,K_SLOPE)

      use rt_global
      use VST_ARRAYS
      use mf_rt_link

      implicit  none
      logical   nitrogen_complex,nitrogen_simple,nitrate_only,O2DIFF
      integer   J,I,K,KK,K_SLOPE
      real      SOURCE,SINK
      real      O2,NO3,NH4,NO2,DOC,N2,
     &          kO2H,kDEN,kNIT1,kNIT2,KO2,KNO3,KNH4,KNO2,KDOC,IO2,
     &          O2H,DENH,NIT1,NIT2,
     &          kdenh,
     &          SPEC_CONC,D1,D2,D3,k1,k2,k3,KSP1,KSP2,KSP3,ISP1,ISP2,
     &          LAMBDA,RED2,RED3,
     &          UNSATZ,AC,DiffMedia,DiffAir,VP,BP,O2_SAT,REAERATION,
     &          SUM,MCMAX,POROS_AVG
      double precision CONC_RG,K_RG
      dimension CONC_RG(NCOMP),K_RG(4,NCOMP)

      nitrogen_complex = .false.
      nitrogen_simple = .false.
      nitrate_only = .true.

      !only nitrate transport is simulated --------------------------------------------------------
      if(nitrate_only) then 
        
        !Retrieve the species concentration
        no3 = conc_rg(1)

        !Retrieve the parameter values
        kden = rc(1)
        kno3 = rc(2)
        
        !Write the chemical reaction rate laws     
        denh = kden*no3*(no3/(kno3+no3))

        !Change in species mass
        k_rg(k_slope,1) = -denh

      !O2-NH4-NO3 system is simulated -------------------------------------------------------------
      elseif(nitrogen_complex) then 

C       Retrieve the species concentration
        o2 = conc_rg(1)
        no3 = conc_rg(2)
        nh4 = conc_rg(3)
        no2 = conc_rg(4)
        doc = conc_rg(5)

C       Retrieve the parameters (either spatially-constant or spatially-variable)
        ko2 = rc(1)
        kno3 = rc(2)
        knh4 = rc(3)
        kno2 = rc(4)
        kdoc = rc(5)
        io2 = rc(6)
        ko2h = vrc(J,I,K,1)
        kden = vrc(J,I,K,2)
        knit1 = vrc(J,I,K,3)
        knit2 = vrc(J,I,K,4)

C       Write the chemical reaction rate laws
        o2h = kO2H*O2*(O2/(KO2+O2))*(DOC/(KDOC+DOC))
        denh = kden*NO3*(NO3/(KNO3+NO3))*(DOC/(KDOC+DOC))*(IO2/(IO2+O2))
        nit1 = kNIT1*NH4*(NH4/(KNH4+NH4))*(O2/(KO2+O2))
        nit2 = kNIT2*NO2*(NO2/(KNO2+NO2))*(O2/(KO2+O2))

C       Change in species mass
        k_rg(k_slope,1) = -o2h
        k_rg(k_slope,2) = (nit2*1.348) - denh !Nitrate
        k_rg(k_slope,3) = -nit1 !Ammonium
        k_rg(k_slope,4) = (nit1*2.556) - nit2 !Nitrite
        k_rg(k_slope,5) = -(o2h*0.9375) - (denh*0.605) !Dissolved Organic Carbon
      

      !Simple three-species system is simulated ---------------------------------------------------
      elseif(nitrogen_simple) then 
        
        NH4 = CONC_RG(1)
        NO2 = CONC_RG(2)
        NO3 = CONC_RG(3)

        k1 = RC(1)
        k2 = RC(2)
        k3 = RC(3)

        K_RG(K_SLOPE,1) = -(K1*NH4)
        K_RG(K_SLOPE,2) = (K1*NH4) - (K2*NO2)
        K_RG(K_SLOPE,3) = (K2*NO2) - (K3*NO3)

      endif


      return
      end
