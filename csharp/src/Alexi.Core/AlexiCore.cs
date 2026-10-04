namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_0.1/src/ALEXI.f — the two-source flux model core
/// (Massman &amp; Weil 1996 / ALEXI surface-energy-balance solver).
/// Field names verbatim; all state flows through PixelState/DailyState.
///
/// Notes on source quirks (ported, not fixed):
/// - `ftheta = fthetanew` in alexi(): fthetanew is an uninitialized implicit
///   local; on the original Linux build it reads as 0.0 (BSS). The port sets
///   s.ftheta = 0f.
/// - `nconv = nconv + 1` in alexi() increments a dead static local in ALEXI.f
///   (the driver's nconv is a different variable); not reproduced.
/// - The fc-retry loop (`go to 1000`) and checksoln call are commented out in
///   the source; the port keeps them as comments.
/// - growPBL: `it = int(dta/dtheta)` may be 0 or negative; Fortran would read
///   tabtheta(0) out of bounds (garbage). The port clamps it to [1, Mt] and
///   documents it.
/// </summary>
public static class AlexiCore
{
    /// <summary>
    /// alexi(ia,ja,ierr,ibad,iter) — solve the two-source flux model for one
    /// clear-sky pixel. ia/ja are unused in the source body. On return,
    /// s.converged carries the solution flag; ierr/ibad carry diagnostics.
    /// </summary>
    public static void Alexi(PixelState s, DailyState d, TextWriter stdout,
        ref int ierr, ref int ibad, out int iter)
    {
        float xlai0 = s.xlai;
        s.fc0 = s.fc;
        iter = 0; // Fortran: uninitialized local; C# requires definite assignment

        AlexiUtl.RunInit(s, d, stdout);
        if (s.badinput)
        {
            ierr = 2;
            AlexiUtl.AlexiErrorcode(s, ierr, stdout);
            s.xlai = xlai0;
            return; // go to 3000
        }

        const float delhmx = 0.1f; // convergence parameter for H [W/m2]
        float fcnew = s.fc;
        s.ftheta = 0f; // fthetanew: uninitialized local in source (BSS zero)

        // 1000 continue — return here to try a new FC (retry loop commented out in source)
        float psima = s.psi0;
        s.psima = psima;
        float hn = s.hn0;
        s.tc1 = s.taobs1;
        s.tc2 = s.taobs2;
        s.ta1 = s.taobs1;
        s.ta2 = s.taobs2;
        s.ts1 = s.taobs1;
        s.ts2 = s.taobs2;
        s.tac1 = s.taobs1;
        s.tac2 = s.taobs2;

        iter = 0;
        ierr = 0;
        s.converged = false;
        s.stopiter = false;

        float hnnew = 0f;
        while (!s.converged)
        {
            s.fc = fcnew;
            AlexiUtl.UpdateFc(s);

            Findhn(s, d, ref hn, ref hnnew, ref s.converged, ref s.stopiter, ref ierr);
            if (s.stopiter)
                break; // go to 2000

            float delthn = hnnew - hn;
            if (MathF.Abs(delthn) > delhmx)
            {
                if (iter < 300)
                {
                    float f = 0.25f / (1 + iter / 20); // Fortran integer division
                    float dh = f * delthn;
                    hn = hn + dh;
                    iter = iter + 1;
                }
                else
                {
                    ierr = 1;
                    break; // go to 2000
                }
            }
            else
            {
                s.converged = true;
            }
        }
        // 2000: escape from flux convergence loop

        if (!s.converged)
        {
            // fc-retry block is commented out in the source — nothing to do.
        }
        else
        {
            // checksoln call is commented out in the source.
            // nconv++ in the source targets a dead static local — not reproduced.
            ierr = 0;
        }

        // 3000:
        AlexiUtl.AlexiErrorcode(s, ierr, stdout);
        s.xlai = xlai0;
    }

    /// <summary>
    /// findhn(hn,hnnew,converged,stopiter,ierr) — iterate the scaling flux HN
    /// by closing the PBL between times t1 and t2.
    /// </summary>
    public static void Findhn(PixelState s, DailyState d, ref float hn, ref float hnnew,
        ref bool converged, ref bool stopiter, ref int ierr)
    {
        s.h1 = hn * s.t1 / s.thrise;
        s.h2 = hn * s.t2 / s.thrise;

        AlexiRad.GetNetRad(s, s.ts1, s.tc1, s.xlwdn1, s.sdn1,
            s.albedo1, s.taubtv1, s.taubtn1,
            out s.rnet1, out s.rnsoil1, out s.rndiv1, out s.swup1, out s.xlwup1);
        AlexiRad.GetNetRad(s, s.ts2, s.tc2, s.xlwdn2, s.sdn2,
            s.albedo2, s.taubtv2, s.taubtn2,
            out s.rnet2, out s.rnsoil2, out s.rndiv2, out s.swup2, out s.xlwup2);

        if (s.rnet2 < 0.0f)
        {
            ierr = 9;
            stopiter = true;
            return; // go to 2000
        }

        AlexiCore.GetSoilHeat(d, s.rnsoil1, s.tloc1, out s.g1);
        AlexiCore.GetSoilHeat(d, s.rnsoil2, s.tloc2, out s.g2);

        // 1500:
        s.zdlamx = 0.4f;
        s.zdlamn = -0.5f;
        FindTa(s, s.trad1, ref s.ta1, ref s.ts1, ref s.tc1,
            s.rnet1, s.rnsoil1, s.rndiv1, s.h1, ref s.hs1, ref s.hc1,
            ref s.xle1, ref s.xles1, ref s.xlec1, ref s.g1,
            ref s.ra1, ref s.rs1, ref s.rx1, s.w1, s.rhocp1, ref s.tac1,
            ref stopiter, ref ierr, 1);
        if (float.IsNaN(s.ta1))
            stopiter = true;
        if (stopiter)
            return; // go to 2000

        s.zdlamx = 0.0f;
        s.zdlamn = -10.0f;
        FindTa(s, s.trad2, ref s.ta2, ref s.ts2, ref s.tc2,
            s.rnet2, s.rnsoil2, s.rndiv2, s.h2, ref s.hs2, ref s.hc2,
            ref s.xle2, ref s.xles2, ref s.xlec2, ref s.g2,
            ref s.ra2, ref s.rs2, ref s.rx2, s.w2, s.rhocp2, ref s.tac2,
            ref stopiter, ref ierr, 2);
        if (float.IsNaN(s.ta2))
            stopiter = true;
        if (stopiter)
            return; // go to 2000

        // ts2 < tc2 with fc > 0.1: body commented out in the source.
        if (s.ts2 < s.tc2 && s.fc > 0.1f)
        {
        }

        if (s.ta2 < s.ta1 + 0.5f)
            s.ta2 = s.ta1 + 0.5f;

        GrowPbl(s, ref hnnew, ref stopiter, ref ierr);

        float h2new = hnnew * s.t2 / s.thrise;
        float h2max = s.rnet2 - s.g2;
        if (h2new > h2max)
            hnnew = h2max * s.thrise / s.t2;
        if (hnnew <= 0.0f)
        {
            stopiter = true;
            ierr = 30;
        }
    }

    /// <summary>
    /// findta(...) — air temperature and component fluxes at one time.
    /// `it` selects time 1 (allow condensation) vs time 2.
    /// </summary>
    public static void FindTa(PixelState s, float trad, ref float ta, ref float ts, ref float tc,
        float rn, float rnsoil, float rndiv, float h, ref float hs, ref float hc,
        ref float xle, ref float xles, ref float xlec, ref float g,
        ref float ra, ref float rs, ref float rx, float wind, float rhocp,
        ref float tac, ref bool stopiter, ref int ierr, int it)
    {
        int niter = 20; // niter>10 keeps LE & H imbalance < 10 W/m2 in FIFE
        float ftheta = s.ftheta;

        ts = (trad - ftheta * tc) / (1.0f - ftheta);

        // Check system energy balance. Does assumed H require negative LE?
        xle = rn - h - g;
        if (xle < 0.0f)
        {
            xle = 0.0f;
            g = rn - h - xle;
            xles = 0.0f;
            xlec = 0.0f;
            hs = rnsoil - g;
            hc = h - hs;
            for (int i = 1; i <= niter; i++)
            {
                AlexiAtmos.GetResistance(s, h, wind, ta, ts, tc, rhocp,
                    out ra, out rs, out rx);
                ta = trad - h * ra / rhocp
                     - ftheta * hc * (rx * 0.5f) / rhocp
                     - (1.0f - ftheta) * (h - hc) * rs / rhocp;
                tac = h * ra / rhocp + ta;
                tc = hc * (rx * 0.5f) / rhocp + tac;
                ts = (trad - ftheta * tc) / (1.0f - ftheta);
            }
            return; // go to 2000 — that's it
        }

        // LE is positive: find canopy and soil temperatures and fluxes.
        AlexiAtmos.GetResistance(s, h, wind, ta, ts, tc, rhocp,
            out ra, out rs, out rx);
        tac = h * ra / rhocp + ta;

        // TC iteration loop
        for (int i = 1; i <= niter; i++)
        {
            float esat = 6.108f * MathF.Pow(10.0f, 7.5f * tac / (237.3f + tac)); // [mb]
            float xlam = (2.501f - 0.00237f * tac) * 1e6f; // latent heat [J/kg]
            float s_ = 2.17e-3f * xlam * esat / ((273.15f + tac) * (273.15f + tac)); // [mb/K]
            xlec = rndiv * 1.3f * s.fg * (s_ / (s_ + 0.66f));
            hc = rndiv - xlec;
            ta = trad - h * ra / rhocp
                 - ftheta * hc * (rx * 0.5f) / rhocp
                 - (1.0f - ftheta) * (h - hc) * rs / rhocp;
            AlexiAtmos.GetResistance(s, h, wind, ta, ts, tc, rhocp,
                out ra, out rs, out rx);
            tac = h * ra / rhocp + ta;
            tc = hc * (rx * 0.5f) / rhocp + tac;
            ts = (trad - ftheta * tc) / (1.0f - ftheta);
        }

        // Soil and canopy energy balances.
        if (rs > 0.0f)
            hs = rhocp * (ts - tac) / rs;
        else
            hs = h - hc;
        xles = rnsoil - hs - g;
        if (it == 1)
            return; // go to 2000 — allow condensation at time t1

        if (xles > -1.0e-3f)
            return; // good solution for soil/canopy fluxes

        // Throttle back LEc: Hs too large because Hc too small.
        xles = 0.0f;
        hs = rnsoil - g;
        hc = h - hs;
        xlec = rndiv - hc;
        if (xlec < -1.0e-3f)
        {
            // Hc too big; it cannot exceed rndiv.
            xlec = 0.0f;
            hc = rndiv;
        }
        ta = trad - h * ra / rhocp
             - ftheta * hc * (rx * 0.5f) / rhocp
             - (1.0f - ftheta) * (h - hc) * rs / rhocp;
        tac = h * ra / rhocp + ta;
        tc = hc * (rx * 0.5f) / rhocp + tac;
        ts = (trad - ftheta * tc) / (1.0f - ftheta);
        if (rs > 0.0f)
            hs = rhocp * (ts - tac) / rs;
        else
            hs = h - hc;

        // 2000: system budget bookkeeping (locals unused in the source) —
        // xlesum, xleresid, hsum, hresid, resid are computed but never used.
    }

    /// <summary>
    /// findta_parallel(...) — parallel-resistance variant (dead code in the
    /// 0.1 build; ported for completeness).
    /// </summary>
    public static void FindTaParallel(PixelState s, float trad, ref float ta, ref float ts, ref float tc,
        float rn, float rnsoil, float rndiv, ref float h, ref float hs, ref float hc,
        ref float xle, ref float xles, ref float xlec, ref float g,
        ref float ra, ref float rs, ref float rx, float wind, float rhocp,
        ref float tac, ref bool stopiter, ref int ierr)
    {
        float ftheta = s.ftheta;

        float tdiff = ts - trad;
        if (tdiff > 15.0f)
            ts = trad + 15.0f;
        if (tdiff < 0.0f)
            ts = trad;
        tc = (trad - ts * (1.0f - ftheta)) / ftheta;
        hc = rhocp * (tc - ta) / (ra + 0.5f * rx);
        hs = (ts - ta) * rhocp / (ra + rs);
        xlec = rndiv - hc;
        xles = rnsoil - hs - g;

        if (xles > -1.0e-3f)
        {
            // good solution for the soil and canopy energy fluxes
        }
        else
        {
            xles = 0.0f;
            hs = rnsoil - g;
            ts = ta + hs * (ra + rs) / rhocp;
            tc = (trad - ts * (1.0f - ftheta)) / ftheta;
            hc = rhocp * (tc - ta) / (ra + 0.5f * rx);
            xlec = rndiv - hc;
            if (xlec < -1.0e-3f)
            {
                xlec = 0.0f;
                hc = rndiv;
                tc = (hc * (ra + 0.5f * rx) / rhocp) + ta;
                ts = (trad - ftheta * tc) / (1.0f - ftheta);
                hs = rhocp * (ts - ta) / (rs + ra);
                g = rnsoil - hs;
            }
        }

        h = hc + hs;
        xle = xlec + xles;
    }

    /// <summary>
    /// growPBL(hnnew,stopiter,ierr) — new HN guess from the precomputed
    /// potential-temperature profile integral.
    /// </summary>
    public static void GrowPbl(PixelState s, ref float hnnew, ref bool stopiter, ref int ierr)
    {
        s.th1 = s.ta1 + 273.15f;
        s.th2 = s.ta2 + 273.15f;
        float dta = s.th2 - s.th1;
        int it = (int)(dta / s.dtheta); // Fortran INT(): truncation toward zero
        if (it > GridDims.Mt)
            return; // go to 2000 (stopiter/ierr assignment is commented out)
        if (it < 1)
            it = 1; // Fortran would read tabtheta(0) out of bounds; clamp (documented)

        float thetasum = s.tabtheta[it - 1];
        s.z2 = s.tabz2[it - 1];
        float thshift = s.thpblz1 - s.th1;
        thetasum = thetasum - thshift * (s.z2 - s.z1);
        hnnew = s.thrise * (s.rhocp1 + s.rhocp2)
              * (s.z2 * s.th2 - s.z1 * s.th1 - thetasum)
              / (s.t2 * s.t2 - s.t1 * s.t1);
    }

    /// <summary>
    /// growPBL2(...) — PBL growth by direct profile integration (dead code in
    /// the 0.1 build; ported for completeness).
    /// </summary>
    public static void GrowPbl2(PixelState s, ref float hnnew, ref bool stopiter, ref int ierr)
    {
        // pressure-adjusted potential temperatures are overwritten below in the source
        s.th1 = s.ta1 + 273.15f;
        s.th2 = s.ta2 + 273.15f;
        if (s.th2 < s.th1)
        {
            hnnew = 5.0f;
            return; // go to 2000
        }
        float tlast = s.th1;
        float zlast = s.z1;
        float thetasum = 0.0f;
        int j = s.jz1;
        if (j < 1)
            j = 1; // thpbl(j-1) with j=1 would read index 0 out of bounds

        bool notdone = true;
        float z2 = 0.0f;
        while (notdone)
        {
            if (j > s.jzmax)
            {
                stopiter = true;
                ierr = 3;
                return; // go to 2000
            }
            float dt = s.thpbl[j - 1] - s.thpbl[j - 2];
            float dz = s.zpbl[j - 1] - s.zpbl[j - 2];
            float t = tlast + dt;
            float z = zlast + dz;
            if (tlast < s.th2 && s.th2 <= t)
            {
                z2 = z;
                notdone = false;
            }
            thetasum = thetasum + t * dz;
            j = j + 1;
            tlast = t;
            zlast = z;
        }

        hnnew = s.thrise * (s.rhocp1 + s.rhocp2)
              * (z2 * s.th2 - s.z1 * s.th1 - thetasum)
              / (s.t2 * s.t2 - s.t1 * s.t1);
    }

    /// <summary>
    /// getsoilheat(rnsoil,g,tloc) — diurnal soil heat flux fraction tied to
    /// solar time; uses the day's sunrise/sunset window (common/sun).
    /// </summary>
    public static void GetSoilHeat(DailyState d, float rnsoil, float tloc, out float g)
    {
        const float a = 0.35f;
        const float b = 100000.0f;
        float tnoon = 0.5f * (d.strt + d.end);
        float t = (tloc - tnoon) * 3600.0f;
        float f = a * MathF.Cos(2.0f * 3.14159f * (t + 10800.0f) / b);
        g = f * rnsoil;
    }

    /// <summary>getsoilheat_old(rnsoil,g) — older height/dispersion form (dead code).</summary>
    public static void GetSoilHeatOld(PixelState s, float rnsoil, out float g)
    {
        const float fmin = 0.15f;
        const float fmax = 0.31f;
        float f;
        if (s.height == 0.0f)
            f = fmin;
        else
        {
            f = s.disp / s.height;
            f = MathF.Max(f, fmin);
            f = MathF.Min(f, fmax);
        }
        g = f * rnsoil;
        g = 0.31f * rnsoil; // TESTING: source overwrites unconditionally
    }

    /// <summary>
    /// checksoln(rnsoil,g,ts,trad,converged,ierr) — sanity-check a converged
    /// solution: F = G/RNSOIL must be >= 0.1 and TS-TRAD <= 15 C.
    /// </summary>
    public static void CheckSoln(float rnsoil, float g, float ts, float trad,
        ref bool converged, ref int ierr)
    {
        float f = g / rnsoil;
        if (f < 0.1f)
        {
            converged = false;
            ierr = 12;
        }
        float tdiff = ts - trad;
        if (tdiff > 15.0f)
        {
            converged = false;
            ierr = 13;
        }
    }
}
