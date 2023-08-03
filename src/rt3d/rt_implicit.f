*     Subroutine:   rt_implicit
*     Purpose:      solve the advection-dispersion reaction using implicit method (solver)
      
      subroutine rt_implicit

      use rt_global
      use GLOBAL,     only:NCOL,NROW,NLAY

      
      !update matrix if nonlinear sorption or multicomponent
      if (TRNOP(4).AND.ISOTHM.GT.1) UPDLHS=.TRUE.
      if (NCOMP.GT.1) UPDLHS=.TRUE.
      
      call btn_imp_update(X(LCCNEWN),IX(LCIBND)) !update concentrations
      
      !go through outer iterations until convergence is achieved
      do ito=1,rt_mxiter
        
        !prepare A and RHS matrices
        call btn_imp(X(LCA),X(LCRHS)) !initialize matrices
        if(TRNOP(1).AND.MIXELM.EQ.0) call adv_imp(X(LCA)) !advection
        if(TRNOP(2)) call dsp_imp(X(LCA),X(LCRHS)) !dispersion
        if(TRNOP(3)) call ssm_imp(X(LCA),X(LCRHS)) !source-sink mixing
     
        !solve the matrix
        call gcg_solve(X(LCCNCG),IX(LCLRCH),X(LCA),X(LCCNEWN),X(LCRHS),
     &                 IX(LCIBND),X(LCQ),X(LCWK))
        
        !if converged, go to next outer iteration
        if(icnvg.eq.1) GOTO 110

      enddo !go to the next iteration

      !if convergence not achieved after rt_mxiter, then stop RT3D
      if(impsol.EQ.1 .AND. icnvg.EQ.0) then
        write(*,*) 'STOP. GCG SOLVER FAILED TO CONVERGE.'
        stop
      endif

110   continue

      !calculate mass budgets for implicit schemes
      if(TRNOP(1).AND.MIXELM.EQ.0) call adv_imp_budget !advection
      if(TRNOP(2)) call dsp_imp_budget !dispersion
      if(TRNOP(3)) call ssm_imp_budget !source-sink mixing


      return
      end
