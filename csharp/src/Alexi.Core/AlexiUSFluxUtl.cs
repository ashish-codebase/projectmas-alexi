namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_0.1/src/USflux_utl.f — potential-flux routines and I/O utilities
/// for the US-flux side. Field names verbatim.
/// Defect note: hourly_flux is called as `call hourly_flux(ibad)` in USflux_cloud.f
/// but takes no arguments; the extra argument is ignored (Fortran arity mismatch).
/// Defect note: `w` in hourly_flux/checkUSinput is an implicit uninitialized real
/// passed as the logical writeme argument; the port passes the pixel's writeme flag.
/// </summary>
public static class AlexiUSFluxUtl
{
    /// <summary>hourly_flux — hourly RN/RNSOIL/G/APAR estimates, integrated to daily totals.</summary>
    public static void HourlyFlux(PixelState s, DailyState d, TextWriter stdout)
    {
        if (s.converged)
        {
            GetTdepart2(ref d.xmc, ref d.xbc, ref d.xcc, s.tc2, s.taobs2, s.tloc2, d.strt, d.end);
            GetTdepart2(ref d.xms, ref d.xbs, ref d.xcs, s.ts2, s.taobs2, s.tloc2, d.strt, d.end);
        }

        bool flag = false;

        for (int ihr = 0; ihr < GridDims.Nohr; ihr++)
        {
            // Initialize to 0 for daily integration step (apar is a local save array, never used)
            d.rnet[ihr] = 0.0f;
            d.rnsoil[ihr] = 0.0f;
            d.g[ihr] = 0.0f;
            d.eref[ihr] = 0.0f;

            if (d.tloc[ihr] < d.strt || d.tloc[ihr] > d.end) d.sdn[ihr] = 0.0f;

            float zen = AlexiRad.GetSunzen(s, d, d.tloc[ihr]);

            if (d.tloc[ihr] >= d.strt && d.tloc[ihr] <= d.end)
            {
                int ibad = 0;
                AlexiUtl.CheckValue("SDN   ", 40, d.sdn[ihr], -1.0f, 1500.0f,
                    ref flag, s.writeme, stdout, ref ibad); // m/s
                AlexiUtl.CheckValue("TA    ", 38, d.ta[ihr], -50.0f, 80.0f,
                    ref flag, s.writeme, stdout, ref ibad); // C
                if (flag)
                {
                    s.badinput = true;
                    ibad = 49;
                    return; // goto 500
                }

                float fclear = AlexiUSFluxRad.GetRadComps(s, zen, ref d.sdn[ihr]);
                // source note: "CHANGE: clump should be clumps(ihr)" — verbatim clump kept
                AlexiRad.GetRadProps(s, zen, out float albedo, out float taubtv,
                    out float taubtn, s.clump);

                float tcproxy, tsproxy;
                if (s.converged || d.xmc != GridDims.Bad)
                {
                    float dtc = d.xmc * d.tloc[ihr] * d.tloc[ihr] + d.xbc * d.tloc[ihr] + d.xcc;
                    float dts = d.xms * d.tloc[ihr] * d.tloc[ihr] + d.xbs * d.tloc[ihr] + d.xcs;
                    if (dtc < 0.0f) dtc = 0.0f;
                    if (dts < 0.0f) dts = 0.0f;
                    tcproxy = d.ta[ihr] + dtc;
                    tsproxy = d.ta[ihr] + dts;
                }
                else
                {
                    tcproxy = d.ta[ihr] + 2.0f;
                    tsproxy = d.ta[ihr] + 4.0f;
                }

                AlexiRad.GetNetRad(s, tsproxy, tcproxy, d.xlwdn[ihr], d.sdn[ihr],
                    albedo, taubtv, taubtn,
                    out d.rnet[ihr], out d.rnsoil[ihr], out _, // rndiv local dummy
                    out d.swup[ihr], out d.xlwup[ihr]);

                AlexiCore.GetSoilHeat(d, d.rnsoil[ihr], d.tloc[ihr], out d.g[ihr]);
                Faopm(s, out d.eref[ihr], d.tloc[ihr], d.sdn[ihr], zen,
                    d.ta[ihr], d.pres[ihr], d.ea[ihr], d.wind[ihr]);
            }
            else
            {
                // Nighttime
                if (d.sdn[ihr] == GridDims.Bad) d.sdn[ihr] = 0.0f;
                int ibad = 0;
                AlexiUtl.CheckValue("SDN   ", 40, d.sdn[ihr], -1.0f, 1500.0f,
                    ref flag, s.writeme, stdout, ref ibad);
                AlexiUtl.CheckValue("TA    ", 38, d.ta[ihr], -50.0f, 80.0f,
                    ref flag, s.writeme, stdout, ref ibad);
                if (flag)
                {
                    s.badinput = true;
                    ibad = 49;
                    return; // goto 500
                }
                float tcproxy = d.ta[ihr];
                float tsproxy = d.ta[ihr];
                AlexiRad.GetRadProps(s, zen, out float albedo, out float taubtv,
                    out float taubtn, s.clump);
                AlexiRad.GetNetRadNight(s, tsproxy, tcproxy, d.xlwdn[ihr], d.sdn[ihr],
                    albedo, taubtv, taubtn,
                    out d.rnet[ihr], out _, out _, // dum,dum locals
                    out d.swup[ihr], out d.xlwup[ihr]);
            }
        }

        // Integrate 24-hr fluxes (tstrt/tend args unused by integrate; kept for parity)
        Integrate(d.rnet, out d.rnday, d.strt, d.end, d.tloc, 24);
        Integrate(d.rnsoil, out d.rnsday, d.strt, d.end, d.tloc, 24);
        Integrate(d.g, out d.gday, d.strt, d.end, d.tloc, 24);
        Integrate(d.sdn, out d.sday, d.strt, d.end, d.tloc, 24);
        Integrate(d.swup, out d.swupday, d.strt, d.end, d.tloc, 24);
        Integrate(d.xlwup, out d.xlwupday, d.strt, d.end, d.tloc, 24);
        Integrate(d.xlwdn, out d.xlwdnday, d.strt, d.end, d.tloc, 24);
        Integrate(d.eref, out d.erefday, d.strt, d.end, d.tloc, 24);
        d.rncday = d.rnday - d.rnsday;
    }

    /// <summary>fao_PM — grass reference ET, FAO Penman-Monteith hourly ETo form.</summary>
    public static void Faopm(PixelState s, out float eref, float tloc, float sdn,
        float zen, float ta, float pres, float ea, float wind)
    {
        const float cp = 1010.0e-6f;   // [MJ/kg-K]
        const float epsilon = 0.622f;
        // sigma parameter in source is unused; omitted (C# warns on unused locals)

        float albedo = 0.23f; // hypothetical grass reference

        float xlam = 2.501f - 0.00237f * ta; // latent heat of vaporization [MJ/kg]
        float psych = 0.1f * pres * cp / (xlam * epsilon);
        float esat = 0.6108f * MathF.Exp(17.2694f * ta / (237.3f + ta));
        float s_ = 17.2694f * 237.3f * esat / MathF.Pow(237.3f + ta, 2.0f);
        float d_ = esat - 0.1f * ea;

        float fclear = AlexiUSFluxRad.GetRadComps(s, zen, ref sdn);
        float rnet = AlexiUSFluxRad.GetNetRadSimple(sdn, ta, fclear, zen, albedo);
        float rnetmj = rnet * 3600.0f / 1.0e6f; // [MJ/m2-hr]
        float gmj = 0.1f * rnetmj;

        float xnum1 = 0.408f * s_ * (rnetmj - gmj);
        float xnum2 = psych * 37.0f * wind * d_ / (ta + 273.15f);
        float xnum = xnum1 + xnum2;
        float xden = s_ + psych * (1.0f + 0.24f * wind);

        float erefmmhr = xnum / xden; // [mm/hr]
        eref = erefmmhr * xlam * 1.0e6f / 3600.0f; // [W/m2]
    }

    /// <summary>integrate — sum hourly fluxes over the day; BAD if any hour is BAD.</summary>
    public static void Integrate(float[] f, out float fint, float tstrt, float tend,
        float[] t, int nhr)
    {
        bool fbad = false;
        fint = 0.0f;
        for (int i = 0; i < nhr; i++)
        {
            if (f[i] == GridDims.Bad) fbad = true;
            fint += f[i];
        }
        float xn = 24.0f / (t[1] - t[0]);
        fint = fbad ? GridDims.Bad : fint / xn * 3600.0f * 24.0f / 1.0e6f; // MJ/day
    }

    /// <summary>getTdepart2 — quadratic surface-temperature departure fit coefficients.</summary>
    public static void GetTdepart2(ref float xm, ref float xb, ref float xc,
        float trad, float ta, float t2, float tr, float ts)
    {
        float a = trad - ta;
        xb = a / (((t2 * t2 - tr * tr) * (tr - ts) / (ts * ts - tr * tr)) + t2 - tr);
        xm = xb * (tr - ts) / (ts * ts - tr * tr);
        xc = -xb * tr - xm * tr * tr;
    }

    /// <summary>checkUSinput — range checks for hourly US-flux inputs.</summary>
    public static void CheckUSInput(PixelState s, DailyState d, TextWriter stdout,
        ref bool flag, ref int ibad)
    {
        flag = false;
        bool allbad = true;

        // Hourly quantities: only check hour slots actually used (sdn>0)
        bool jumped = false; // Fortran `goto 100` skips the rest of the routine
        for (int ihr = 0; ihr < d.nohrin; ihr++)
        {
            if (d.sdn[ihr] != GridDims.Bad) allbad = false;
            if (d.sdn[ihr] > 0.0f)
            {
                AlexiUtl.CheckValue("PRES  ", 37, d.pres[ihr], 0.0f, 1500.0f,
                    ref flag, s.writeme, stdout, ref ibad); // mbar
                if (flag) { jumped = true; break; }
            }
        }
        if (!jumped)
        {
            // Check for missing hours
            for (int ihr = 1; ihr < d.nohrin; ihr++)
            {
                float dtime = d.tloc[ihr] - d.tloc[ihr - 1];
                if (dtime > 1.0f)
                    stdout.Write("CRAS HR GAP BETWEEN TLOCS " + d.tloc[ihr - 1] + " " + d.tloc[ihr] + "\n");
            }

            // Canopy parameters
            float xclass = d.iclass;
            AlexiUtl.CheckValue("ICLASS", 42, xclass, 1.0f, 27.0f, ref flag, s.writeme, stdout, ref ibad);
            AlexiUtl.CheckValue("FG    ", 44, s.fg, 0.0f, 1.0f, ref flag, s.writeme, stdout, ref ibad);
            AlexiUtl.CheckValue("XLAI  ", 22, s.xlai, 0.0f, 10.0f, ref flag, s.writeme, stdout, ref ibad);

            // Insolation data
            if (allbad)
            {
                flag = true;
                ibad = 47;
            }
        }
    }

    /// <summary>setPBLheights — ZPBLI at 200 m spacings from 0 m.</summary>
    public static void SetPBLHeights(PixelState s)
    {
        for (int j = 1; j <= GridDims.Mli; j++)
            s.zpbli[j - 1] = 0.0f + (j - 1.0f) * 200.0f;
    }
}
