      subroutine alexi(ia,ja,ierr,ibad,iter)
c     **************************************************************
c     *
c     *
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      include 'ALEXI.inc'
      save

      xlai0=xlai
      fc0=fc
c  Initialize run
      call runinit(badinput,ibad,writeme)
      if (badinput) then
         ierr=2
         go to 3000
      endif

      delhmx=0.1 	! Convergence parameter for H [W/m2]
      fcnew=fc
      ftheta=fthetanew

1000  continue		! Return to here to try new value of FC

c  Initialize state variables
      psima=psi0
      hn=hn0
      tc1=taobs1
      tc2=taobs2
      ta1=taobs1
      ta2=taobs2
      ts1=taobs1
      ts2=taobs2
      tac1=taobs1
      tac2=taobs2

c  Initialize counters & flags
      iter=0
      ierr=0
      converged=.FALSE.
      stopiter=.FALSE.

c  BEGIN FLUX CONVERGENCE LOOP
c-----------------------------------------------------------------
c  Guess an initial value for Hn and iterate until the system of
c  equations is solved

      do while (.not.converged)

         fc=fcnew
         call updatefc

c  Use ALEXI algorithm to find HN, the scaling flux.
c  Relax nudge factor if not converging.
         call findhn (hn,hnnew,converged,stopiter,ierr)
!         write(6,100) iter,hn,hnnew,tc2,ts2
 100     format("ITER HN HNNEW TC TS",i5,4f9.2)
         if (stopiter) go to 2000
         delthn=hnnew-hn
         if (abs(delthn).gt.delhmx) then
            if(iter.lt.300)then
               f=0.25/(1+int(iter/20))
               dh=f*delthn
            else
               ierr=1
               go to 2000
            endif
            hn=hn+dh
            iter=iter+1
         else
            converged=.TRUE.
         end if

      enddo	! while not converged

c  END FLUX CONVERGENCE LOOP
c-----------------------------------------------------------------
2000  continue	! Escape from flux convergence loop


c  SOLUTION DID NOT CONVERGE:
c ------------------------------------------------------------------
      if (.not.converged) then
c 		Move fcnew closer to 0.5 and try again.
c  		Bail if fcnew reaches 0.5 and it still
c               doesn't converge.

c          fcnew=fcnew+sign(0.01,0.5-fcnew)
c          if (abs(fcnew-0.5).le.1e-2) then
c             nfail=nfail+1
c          else
c             go to 1000
c          endif

c  SOLUTION CONVERGED:
c ------------------------------------------------------------------
      else
c        call checksoln(g2,rnsoil2,ts2,trad2,converged,ierr)
         if (converged) then
c           write(6,*) "CONVERGED"
           nconv=nconv+1
           ierr=0
c  		Save these values as best first guess for next pixel

c           hn0=hn
c           psi0=psima
         endif
      endif

3000  continue
      call ALEXI_errorcode(writeme,ierr)

      xlai=xlai0

      return
      end

      subroutine findhn (hn,hnnew,converged,stopiter,ierr)
c     **************************************************************
c     *
c     *
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      include 'ALEXI_parm.inc'
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &       jzmax,jz1,nlev
      common/stability/zdlamx,zdlamn
      common/view/theta,ftheta
      logical converged,stopiter
      save

      h1=hn*t1/thrise
      h2=hn*t2/thrise
c      write(6,*) "H2 = ", h2, hn, t2, thrise
1000  continue      	! Return to here to try new value of FC

c  Compute net radiation and soil heat flux at times t1 and t2

      call getnetrad(ts1,tc1,xlwdn1,sdn1,albedo1,taubtv1,taubtn1,
     &               rnet1,rnsoil1,rndiv1,swup1,xlwup1)
      call getnetrad(ts2,tc2,xlwdn2,sdn2,albedo2,taubtv2,taubtn2,
     &               rnet2,rnsoil2,rndiv2,swup2,xlwup2)
!      write(6,*) "RNET = ", rnet1, rnet2
!      write(6,*) "ts = ", ts1, ts2
!      write(6,*) "tc = ", tc1, tc2
!      write(6,*) "xlwdn = ", xlwdn1, xlwdn2
!      write(6,*) "albedo = ", albedo1, albedo2
!      write(6,*) "swup = ", swup1, swup2
!      write(6,*) "xlwup2 = ", xlwup1, xlwup2
!      write(6,*) "taubtv1 = ", taubtv1, taubtv2
!      write(6,*) "taubtn = ", taubtn1, taubtn1
      if(rnet2.lt.0.0)then
         ierr=9
         stopiter=.TRUE.
         go to 2000
      endif

      call getsoilheat(rnsoil1,g1,tloc1)
      call getsoilheat(rnsoil2,g2,tloc2)

c  Find air temperature and component fluxes at times t1 and t2

1500  continue
      zdlamx=.4
      zdlamn=-.5
c      write(6,*),trad1,trad2,ta1,ta2,w1,w2,g1,g2
      call findta (trad1,ta1,ts1,tc1,rnet1,rnsoil1,rndiv1,h1,hs1,hc1,
     &             xle1,xles1,xlec1,g1,ra1,rs1,rx1,w1,rhocp1,tac1,
     &             stopiter,ierr,1)
      if (isnan(ta1).eqv..TRUE.) then
       stopiter=.TRUE.
      endif


      if(stopiter)go to 2000

c      zdlamx=-.2		! System should not be stable at time t2
      zdlamx=0.0
      zdlamn=-10.0
      call findta (trad2,ta2,ts2,tc2,rnet2,rnsoil2,rndiv2,h2,hs2,hc2,
     &             xle2,xles2,xlec2,g2,ra2,rs2,rx2,w2,rhocp2,tac2,
     &             stopiter,ierr,2)
c      write(6,*) trad2, ta2, ts2, tc2
c      write(6,*)xlai,fc,height
c      write(6,*)ra2,rs2,rx2,w2
c      write(6,*)class
      if (isnan(ta2)) then
       stopiter=.TRUE.
      endif


      if(stopiter)go to 2000

c  Trap bad solutions for ts and tc and alter fc accordingly.
c  Appears to be a high fc-triggered condition

      if (ts2.lt.tc2.and.fc.gt.0.1) then
c         stopiter=.TRUE.
c         ierr=11
c         go to 2000
      endif

c  Ta2 must be > Ta1 for PBL model.
c  Ta2 < Ta1 appears to be a low fc, high hn-triggered condition.
c  Setting ta2=ta1+2 yields low hn guess for next iter.

      if(ta2.lt.ta1+0.5)then
         ta2=ta1+0.5
      endif

c  Now, are Ta1, Ta2 consistent with H1 -> H2 heat input into the PBL?
c  Find new value of Hn and return for convergence check.

      call growPBL(hnnew,stopiter,ierr)
c      write(6,*) "hnnew = ", hnnew, thrise
c  Limit the next guess for H2 to be <= RNET2-G2.
c  Otherwise we may wander into Never-never Land.

      h2new=hnnew*t2/thrise
!      write(6,*) "h2new = ", hnnew, t2, thrise
      h2max=rnet2-g2
      if(h2new.gt.h2max)hnnew=h2max*thrise/t2
      if (hnnew<=0.0)then
        stopiter=.TRUE.
        ierr=30
      endif

2000  continue		! Escape from flux computation

      return
      end

      subroutine growPBL(hnnew,stopiter,ierr)
c     **************************************************************
c     *
c     *
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      include 'ALEXI_parm.inc'
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &       jzmax,jz1,nlev
      common/pbl2/tabtheta(mt),tabz2(mt),dtheta,itmax,thpblz1
      logical stopiter,notdone
      save

c  Find lookup table index IT corresponding with current
c  change in air temperature

      th1=ta1+273.15
      th2=ta2+273.15
      dta=th2-th1
      it=int(dta/dtheta)
c      write(6,*) th1, th2, dta, it, itmax
      if (it.gt.8000) then
c       write(6,*) "inside stop", it, itmax
c         stopiter=.TRUE.
c         ierr=3
         go to 2000
      endif

      thetasum=tabtheta(it)
      z2=tabz2(it)
c      write(6,*) "it = ", it, thetasum, z2
c  Shift profile in temp by constant amount THSHIFT such that at
c  height Z1 it has value TH1.
c  Must subtract this shift (integrated over Z) off precomputed
c  integral in TABTHETA lookup table.

      thshift=thpblz1-th1
      thetasum=thetasum-thshift*(z2-z1)
      hnnew=thrise*(rhocp1+rhocp2)*(z2*th2-z1*th1-thetasum)/
     &      (t2*t2-t1*t1)
 2000 continue		! Escape PBL growth routine

      return
      end


      subroutine growPBL2(hnnew,stopiter,ierr)
c     **************************************************************
c     *
c     *
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      include 'ALEXI_parm.inc'
      common/data1/t1,tloc1,trad1,taobs1,ea1,w1,pres1,ta1,ts1,tc1,
     &       h1,hs1,hc1,xle1,xles1,xlec1,g1,rnet1,rnsoil1,rndiv1,
     &       sdn1,par1,xlwdn1,ra1,rs1,rx1,th1,z1,rhocp1,tac1,
     &       albedo1,taubtv1,taubtn1,zen1,tb1,xlwup1,swup1
      common/data2/t2,tloc2,trad2,taobs2,ea2,w2,pres2,ta2,ts2,tc2,
     &       h2,hs2,hc2,xle2,xles2,xlec2,g2,rnet2,rnsoil2,rndiv2,
     &       sdn2,par2,xlwdn2,ra2,rs2,rx2,th2,z2,rhocp2,tac2,
     &       albedo2,taubtv2,taubtn2,zen2,tb2,xlwup2,swup2
      common/pbl/thrise,zpbli(mli),thpbli(mli),zpbl(ml),thpbl(ml),
     &       jzmax,jz1,nlev
      logical stopiter,notdone
      save

c             Convert ta1 and ta2 to potential temperatures
      th1=((ta1+273.15)*(1000./pres1)**0.286)
      th2=((ta2+273.15)*(1000./pres2)**0.286)
      th1=ta1+273.15
      th2=ta2+273.15
      if (th2.lt.th1) then
         hnnew=5
         go to 2000
      endif
      tlast=th1
      zlast=z1
      thetasum=0.
      notdone=.TRUE.
      j=jz1

      do while (notdone)
         if (j.gt.jzmax) then
            stopiter=.TRUE.
            ierr=3
            go to 2000
         endif
         dt=thpbl(j)-thpbl(j-1)
         dz=zpbl(j)-zpbl(j-1)
         t=tlast+dt
         z=zlast+dz
         if (tlast.lt.th2.and.th2.le.t) then
            z2=z
            notdone=.FALSE.
         endif
         thetasum=thetasum+t*dz
         j=j+1
         tlast=t
         zlast=z
      enddo 		! while not done

      hnnew=thrise*(rhocp1+rhocp2)*(z2*th2-z1*th1-thetasum)/
     &      (t2*t2-t1*t1)
 2000 continue		! Escape PBL growth routine

      return
      end

      subroutine findta (trad,ta,ts,tc,rn,rnsoil,rndiv,h,hs,hc,
     &                   xle,xles,xlec,g,ra,rs,rx,wind,rhocp,tac,
     &                   stopiter,ierr,it)
c     **************************************************************
c     *
c     *
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/view/theta,ftheta
      logical stopiter
      save

      niter=20		! niter>10 keeps LE & H imbalance < 10Wm-2 in FIFE
c      niter=1
      ts=(trad-ftheta*tc)/(1.-ftheta)

c  Check system energy balance.  Does assumed value of H require negative LE?
c  If so, set all LE components to 0 and absorb residual into soil heat.

      xle=rn-h-g
      if(xle.lt.0.)then
         xle=0.
         g=rn-h-xle
         xles=0.
         xlec=0.
         hs=rnsoil-g
         hc=h-hs
         do i=1,niter
            call getresistance (h,wind,ta,ts,tc,ra,rs,rx,rhocp)
            ta=trad-h*ra/rhocp-ftheta*hc*(rx*0.5)/rhocp-
     &                  (1.-ftheta)*(h-hc)*rs/rhocp
            tac=h*ra/rhocp+ta
            tc=hc*(rx*0.5)/rhocp+tac
            ts=(trad-ftheta*tc)/(1.-ftheta)
         enddo
         go to 2000	! That's that - get outa here.
      endif

c  LE is positive: proceed to find canopy and soil temperatures and fluxes

c		     C A N O P Y   P R O P E R T I E S
c                    ---------------------------------
      call getresistance (h,wind,ta,ts,tc,ra,rs,rx,rhocp)
      tac=h*ra/rhocp+ta

c 		TC Iteration loop
      do i=1,niter
         esat=6.108*10**(7.5*tac/(237.3+tac))		! [mb]
         xlam=(2.501-0.00237*tac)*1e6			! latent heat o' vaporization [J/kg]
         s=2.17e-3*xlam*esat/((273.15+tac)*(273.15+tac))! [mb/K]
         xlec=rndiv*1.3*fg*(s/(s+0.66))
c        if (xlec.gt.rndiv) xlec=rndiv
         hc=rndiv-xlec
         ta=trad-h*ra/rhocp-ftheta*hc*(rx*0.5)/rhocp-
     &                  (1.-ftheta)*(h-hc)*rs/rhocp
         call getresistance (h,wind,ta,ts,tc,ra,rs,rx,rhocp)
         tac=h*ra/rhocp+ta
         tc=hc*(rx*0.5)/rhocp+tac
         ts=(trad-ftheta*tc)/(1.-ftheta)
      enddo
c 		End TC Iteration loop

c		       S O I L   P R O P E R T I E S
c                      -----------------------------

c  Now check the soil and canopy energy balances.
c  If LEs < 0  it means that the soil is probably dry and the canopy is
c  not transpiring at full capacity.
c  Throttle back LEc and get new Hc, then get consistent Ta and Tc.

      if (rs.gt.0) then
        hs=rhocp*(ts-tac)/rs
      else
        hs=h-hc
      endif
!      write(6,*) rnsoil, hs, g, ts, tac
      xles=rnsoil-hs-g
      if (it.eq.1) go to 2000		! We will allow condensation at time T1

c		| Check soil and canopy ET for positivity:
      if (xles.gt.-1.e-3)then
c              	| Good job.  A solution for the soil and canopy
c              	| energy fluxes has been reached.
         go to 2000
      endif

c  		| Hs is too large because Hc is too small because
c		| LEc is too large.
c      write(6,*)'Im screwing with the soil budget'
      xles=0.
      hs=rnsoil-g
      hc=h-hs
      xlec=rndiv-hc
      if (xlec.lt.-1.e-3)then
c  		| We're here because Hc is too big.   It can't
c		| be larger than rndiv.
c        write(6,*)'Im screwing with the canopy budget'
         xlec=0
         hc=rndiv
      endif
c		| Update Ta, Tc and Ts
      ta=trad-h*ra/rhocp-ftheta*hc*(rx*0.5)/rhocp-
     &                  (1.-ftheta)*(h-hc)*rs/rhocp
      tac=h*ra/rhocp+ta
      tc=hc*(rx*0.5)/rhocp+tac
      ts=(trad-ftheta*tc)/(1.-ftheta)
      if (rs.gt.0) then
        hs=rhocp*(ts-tac)/rs
      else
        hs=h-hc
      endif

2000  continue		! Escape TA computation

c		       S Y S T E M   B U D G E T
c                      -------------------------

c  Does system energy budgets balance?

c      write(6,*),ra,rs,rx,fc,xlai
      xlesum=xles+xlec
      xleresid=xle-xlesum
      hsum=hs+hc
      hresid=hsum-h
      resid=rn-h-xle-g
c      if (abs(xleresid).gt.0.5) write (6,*) 'Imbalance in LE!',xleresid
c      if (abs(hresid).gt.0.5)   write (6,*) 'Imbalance in H!',hresid
c      if (resid.ge.1e-3) write (6,*) 'Imbalance in energy budget!',resid

      return
      end

      subroutine findta_parallel (trad,ta,ts,tc,rn,rnsoil,rndiv,
     &                   h,hs,hc,xle,xles,xlec,g,ra,rs,rx,wind,
     &                   rhocp,tac,stopiter,ierr)
c     **************************************************************
c     *
c     *  If TS computed through series model is unreasonable in
c     *  comparison with TRAD (TS-TRAD>15C or TS-TRAD<0), move
c     *  into a parallel resistance form.  BK has found this to be
c     *  more robust at high cover, where this condition typically
c     *  occurs.
c     *
c     *  Martha Anderson
c     *  Created:  01/14/02
c     *
c     ***************************************************************
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump
      common/view/theta,ftheta
      logical stopiter
      save

      tdiff=ts-trad
      if(tdiff.gt.15)ts=trad+15.0
      if(tdiff.lt.0) ts=trad
      tc=(trad-ts*(1-ftheta))/ftheta
      hc=rhocp*(tc-ta)/(ra+0.5*rx)
      hs=(ts-ta)*rhocp/(ra+rs)
      xlec=rndiv-hc
      xles=rnsoil-hs-g

c  If LEs < 0  it means that the soil is probably dry and the canopy is
c  not transpiring at full capacity.
c  Throttle back LEc and get new Hc, then get consistent Ta and Tc.

      if (xles.gt.-1.e-3)then
c              	| Good job.  A solution for the soil and canopy
c              	| energy fluxes has been reached.
         go to 2000
      endif
c  		| Hs is too large because Hc is too small because
c		| LEc is too large.
c      write(6,*)'Im screwing with the soil budget'
      xles=0.
      hs=rnsoil-g

      ts=ta+hs*(ra+rs)/rhocp
      tc=(trad-ts*(1-ftheta))/ftheta
      hc=rhocp*(tc-ta)/(ra+0.5*rx)
      xlec=rndiv-hc

      if (xlec.lt.-1.e-3)then
c  		| We're here because Hc is too big.   It can't
c		| be larger than rndiv. Use parallel form to
c               | avoid iteration
c         write(6,*)'Im screwing with the canopy budget'
         xlec=0
         hc=rndiv
         tc=(hc*(ra+0.5*rx)/rhocp)+ta
         ts=(trad-ftheta*tc)/(1.-ftheta)
         hs=rhocp*(ts-ta)/(rs+ra)
         g=rnsoil-hs
      endif

2000  continue

      h=hc+hs
      xle=xlec+xles

      return
      end

      subroutine getsoilheat (rnsoil,g,tloc)
c     **************************************************************
c     *
c     *  This ties G to DISP/HEIGHT.  At very low covers, wind speed
c     *  near soil is very high and most of RN at soil is convected
c     *  away from soil as H. DISP/HEIGHT exhibits this behavior too
c     *  (approaches 0 rapidly as fc->0, for fc<0.05)
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      common/sun/strt,end
      common/cover/iclass,beta,gvs,xndvi,fpar,dstom

c      a=0.31
c      b=74000
      a=0.35
!      if (fpar.gt.0.10) then
!       a=0.55
!      endif
      b=100000
      tnoon=0.5*(strt+end)
c      a=xndvi
c      if (t.lt.strt) t=strt
c      if (t.gt.end) t=end
      t=(tloc-tnoon)*3600.
      f=a*cos(2*3.14159*(t+10800.)/b)
c      f=0.35
      g=f*rnsoil

      return
      end

      subroutine getsoilheat_old (rnsoil,g)
c     **************************************************************
c     *
c     *  This ties G to DISP/HEIGHT.  At very low covers, wind speed
c     *  near soil is very high and most of RN at soil is convected
c     *  away from soil as H. DISP/HEIGHT exhibits this behavior too
c     *  (approaches 0 rapidly as fc->0, for fc<0.05)
c     *
c     *  Martha Anderson
c     *  Created:  02/14/01
c     *
c     ***************************************************************
      common/canopy/refhtw,disp,z0,height,xl,xlai,fc,fg,clump

      fmin=0.15
      fmax=0.31
      if (height.eq.0.0) then
        f=fmin
      else
        f=disp/height
        f=max(f,fmin)
        f=min(f,fmax)
      endif
      g=f*rnsoil
c  TESTING:
      g=0.31*rnsoil

      return
      end

      subroutine checksoln (rnsoil,g,ts,trad,converged,ierr)
c     **************************************************************
c     *
c     *  Check F=G/RNSOIL.  If it has been modified such that F<0.1
c     *  (including F negative), then it is likely that we may have
c     *  a mismatch (spatial or temporal) between TRAD and FC.  Best
c     *  to flag this point BAD.
c     *
c     *  Martha Anderson
c     *  Created:  01/28/02
c     *
c     ***************************************************************
      logical converged

      f=g/rnsoil
      if (f.lt.0.1) then
         converged=.FALSE.
         ierr=12
      endif
      tdiff=ts-trad
      if (tdiff.gt.15) then
         converged=.FALSE.
         ierr=13
      endif

      return
      end
