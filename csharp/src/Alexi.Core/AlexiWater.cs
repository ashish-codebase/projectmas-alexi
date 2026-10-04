namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_0.1/src/ALEXI_water.f — open-water variant of the two-source
/// flux model. Field names verbatim.
///
/// Note: like alexi(), the source's `nconv = nconv + 1` targets a dead static
/// local in this file; not reproduced. The fc-retry block is commented out.
/// </summary>
public static class AlexiWater
{
    /// <summary>
    /// alexi_water(ia,ja,ierr,ibad,iter) — solve the water-body model for one
    /// clear-sky pixel. On return s.converged carries the solution flag;
    /// ts2/tc2 are set to trad2 (the water body has no canopy/soil split).
    /// </summary>
    /// <summary>
    /// alexi_water — open-water variant. Renamed AlexiWaterRun: C# forbids a
    /// member with the same name as its enclosing class (Fortran name kept
    /// in this doc comment).
    /// </summary>
    public static void AlexiWaterRun(PixelState s, DailyState d, TextWriter stdout,
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
            s.ts2 = s.trad2;
            s.tc2 = s.trad2;
            return; // go to 3000 (3000 falls through to the ts2/tc2 assignment)
        }

        const float delhmx = 0.1f;
        float fcnew = s.fc;
        s.ftheta = 0f; // fthetanew: uninitialized local (BSS zero)

        s.psima = s.psi0;
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
            // fc=fcnew / updatefc are absent in the water path of the source.
            FindFluxWater(s, ref hn, ref hnnew, ref s.converged, ref s.stopiter, ref ierr);
            if (s.stopiter)
                break; // go to 2000

            float delthn = hnnew - hn;
            if (MathF.Abs(delthn) > delhmx)
            {
                if (iter < 300)
                {
                    float f = 0.25f / (1 + iter / 20);
                    hn = hn + f * delthn;
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

        if (!s.converged)
        {
            // fc-retry block is commented out in the source.
        }
        else
        {
            // nconv++ targets a dead static local — not reproduced.
            ierr = 0;
        }

        // 3000:
        AlexiUtl.AlexiErrorcode(s, ierr, stdout);
        s.xlai = xlai0;
        s.ts2 = s.trad2;
        s.tc2 = s.trad2;
    }

    /// <summary>
    /// findflux_water(hn,hnnew,converged,stopiter,ierr) — water-body flux
    /// closure at times t1 and t2.
    /// </summary>
    public static void FindFluxWater(PixelState s, ref float hn, ref float hnnew,
        ref bool converged, ref bool stopiter, ref int ierr)
    {
        float h1 = hn * s.t1 / s.thrise;
        float h2 = hn * s.t2 / s.thrise;
        // lwup1/lwup2 are implicit locals in the source (not common state)
        float lwup1, lwup2;

        GetNetRadWater(s.trad1, s.xlwdn1, s.sdn1,
            out s.rnet1, out s.swup1, out lwup1);
        GetNetRadWater(s.trad2, s.xlwdn2, s.sdn2,
            out s.rnet2, out s.swup2, out lwup2);
        if (s.rnet2 < 0.0f)
        {
            ierr = 9;
            stopiter = true;
            return; // go to 2000
        }

        GetSoilHeatWater(s.rnet1, out s.g1);
        GetSoilHeatWater(s.rnet2, out s.g2);

        s.zdlamx = 0.4f;
        s.zdlamn = -0.5f;
        FindTaWater(s, s.trad1, ref s.ta1, s.rnet1, h1, ref s.xle1, ref s.g1,
            ref s.ra1, ref s.rs1, ref s.rx1, s.w1, s.rhocp1,
            ref stopiter, ref ierr);
        if (stopiter)
            return; // go to 2000

        s.zdlamx = 0.0f;
        s.zdlamn = -10.0f;
        FindTaWater(s, s.trad2, ref s.ta2, s.rnet2, h2, ref s.xle2, ref s.g2,
            ref s.ra2, ref s.rs2, ref s.rx2, s.w2, s.rhocp2,
            ref stopiter, ref ierr);
        if (stopiter)
            return; // go to 2000

        // ts2 < tc2 with fc > 0.1: body commented out in the source.
        if (s.ts2 < s.tc2 && s.fc > 0.1f)
        {
        }

        if (s.ta2 < s.ta1 + 0.5f)
            s.ta2 = s.ta1 + 0.5f;

        AlexiCore.GrowPbl(s, ref hnnew, ref stopiter, ref ierr);

        float h2new = hnnew * s.t2 / s.thrise;
        float h2max = s.rnet2 - s.g2;
        if (h2new > h2max)
            hnnew = h2max * s.thrise / s.t2;
        // note: the water path has no hnnew<=0 stopiter check in the source
    }

    /// <summary>findta_water(...) — air temperature over a water body.</summary>
    public static void FindTaWater(PixelState s, float trad, ref float ta, float rn, float h,
        ref float xle, ref float g, ref float ra, ref float rs, ref float rx,
        float wind, float rhocp, ref bool stopiter, ref int ierr)
    {
        AlexiAtmos.GetResistance(s, h, wind, ta, ta, ta, rhocp,
            out ra, out rs, out rx);
        ta = trad - h * ra / rhocp;
        xle = rn - h - g;
    }

    /// <summary>getsoilheat_water(rnet,g) — fixed 0.55 fraction of net radiation.</summary>
    public static void GetSoilHeatWater(float rnet, out float g)
    {
        g = 0.55f * rnet;
    }

    /// <summary>
    /// getnetrad_water(trad,xlwdn,sdn,rnet,swup,lwup) — net radiation over
    /// water: sky minus water-body emission plus shortwave minus albedo.
    /// </summary>
    public static void GetNetRadWater(float trad, float xlwdn, float sdn,
        out float rnet, out float swup, out float lwup)
    {
        const float sigma = 5.67e-8f;
        const float albedo = 0.1f;
        // emissivity=0.99 is folded into rwater in the source
        float rsky = xlwdn;
        float rwater = 0.99f * sigma * MathF.Pow(trad + 273.15f, 4.0f);
        float rnetlw = rsky - rwater;
        float rnetsw = sdn * (1.0f - albedo);
        rnet = rnetlw + rnetsw;
        swup = sdn * albedo;
        lwup = rwater;
    }
}
