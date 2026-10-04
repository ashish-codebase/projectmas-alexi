using System.Globalization;

namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_utl.f. Fortran names preserved. State is passed explicitly
/// (no static mutable fields); all reals are float; stdout is injected.
/// </summary>
public static class AlexiUtl
{
    /// <summary>set_constants — domain-wide constants and clear-sky partitioning.</summary>
    public static void SetConstants(PixelState s)
    {
        s.xk = 0.4f;                       // von Karman's constant
        s.pi = 3.1415926537f;
        s.thrise = 3.0f * 3600.0f;         // flux rise time [s]
        s.cp = 1010.0f;
        s.z1 = 50.0f;                      // initial boundary-layer height [m]
        s.fcbare = 0.0f;                   // fc <= fcbare = bare soil

        s.fvisclr = 0.5f;
        s.fnirclr = 0.5f;
        s.dirvisclr = 0.8f;
        s.difvisclr = 0.2f;
        s.dirnirclr = 1.0f;
        s.difnirclr = 0.0f;

        s.hn0 = 80.0f;
        s.psi0 = 0.0f;
    }

    /// <summary>
    /// getprofile(z1) — interpolate THPBLI/ZPBLI (mli levels) to 1-m spacing
    /// arrays THPBL/ZPBL (ml levels); find jz1 where height == z1 (exact float
    /// compare, as in Fortran) and jzmax = last valid index.
    /// </summary>
    public static void GetProfile(PixelState s)
    {
        int j0 = 1;
        float mxhtpbl = s.zpbli[s.nlev - 1];
        s.zpbl[0] = s.zpbli[0];
        s.thpbl[0] = s.thpbli[0];

        bool exited = false;
        int j = 2;
        for (; j <= GridDims.Ml; j++)
        {
            float ht = s.zpbli[0] + j;
            if (ht > mxhtpbl)
            {
                s.jzmax = j - 1;
                exited = true;
                break;
            }
            // Fortran: do while (ht .gt. zpbli(j0+1)) j0=j0+1  -> 0-based index j0
            while (ht > s.zpbli[j0])
                j0++;
            float f = (ht - s.zpbli[j0 - 1]) / (s.zpbli[j0] - s.zpbli[j0 - 1]);
            float temp = s.thpbli[j0 - 1] + f * (s.thpbli[j0] - s.thpbli[j0 - 1]);
            s.zpbl[j - 1] = ht;
            s.thpbl[j - 1] = temp;
            if (ht == s.z1)
                s.jz1 = j;
        }
        if (!exited)
            s.jzmax = j - 1; // loop ran to ml: j = ml+1, jzmax = ml
    }

    /// <summary>
    /// getpbltable(z1, ierr) — trapezoid-integrate the potential-temperature
    /// profile upward from z1 in dtheta steps; fill tabtheta/tabz2.
    /// `go to 100` (profile exhausted) -> itmax = it-1; full loop -> itmax = mt.
    /// ierr is unused in the Fortran source.
    /// </summary>
    public static void GetPblTable(PixelState s)
    {
        float z1 = s.z1;
        s.dtheta = 0.01f;

        int j = 1;
        while (z1 > s.zpbli[j]) // zpbli(j+1) -> index j
            j++;
        float f = (z1 - s.zpbli[j - 1]) / (s.zpbli[j] - s.zpbli[j - 1]);
        s.thpblz1 = s.thpbli[j - 1] + f * (s.thpbli[j] - s.thpbli[j - 1]);

        float tlast = s.thpblz1;
        float zlast = z1;
        float tint = 0.0f;
        bool exhausted = false;

        for (int it = 1; it <= GridDims.Mt; it++)
        {
            float t = tlast + s.dtheta;
            while (t > s.thpbli[j]) // thpbli(j+1) -> index j
            {
                j++;
                if (j >= GridDims.Mli) // j.ge.mli (1-based): profile exhausted
                {
                    s.itmax = it - 1;
                    exhausted = true;
                    break;
                }
            }
            if (exhausted)
                break;
            f = (t - s.thpbli[j - 1]) / (s.thpbli[j] - s.thpbli[j - 1]);
            float z = s.zpbli[j - 1] + f * (s.zpbli[j] - s.zpbli[j - 1]);
            float dz = z - zlast;
            tint += 0.5f * dz * (t + tlast); // trapezoid rule
            s.tabtheta[it - 1] = tint;
            s.tabz2[it - 1] = z;
            tlast = t;
            zlast = z;
        }
        if (!exhausted)
            s.itmax = GridDims.Mt;
    }

    /// <summary>
    /// runinit — per-pixel initialization. Live path only: the fc-branch for
    /// dtloc is overwritten by dtloc=1.50 in the source, so the branch is dead.
    /// strt/end come from common/sun (DailyState), shared across the day.
    /// </summary>
    public static void RunInit(PixelState s, DailyState d, TextWriter stdout)
    {
        int ibad = 0;
        CheckInput(s, stdout, ref ibad);
        if (s.badinput)
            return; // go to 2000

        // partitioning factors = clear-sky values
        s.fvis = s.fvisclr;
        s.fnir = s.fnirclr;
        s.dirvis = s.dirvisclr;
        s.difvis = s.difvisclr;
        s.dirnir = s.dirnirclr;
        s.difnir = s.difnirclr;

        s.zen1 = AlexiRad.GetSunzen(s, d, s.tloc1);
        s.zen2 = AlexiRad.GetSunzen(s, d, s.tloc2);

        float dtloc = 1.50f; // live path: source hardcodes this after the dead branch
        s.t1 = (s.tloc1 - d.strt - dtloc) * 3600.0f;
        s.t2 = (s.tloc2 - d.strt - dtloc) * 3600.0f;

        // directional thermal emissivity coefficients
        s.eleaf = 0.97f;
        float cosTheta = MathF.Cos(s.theta);
        float b1 = -1.97f * cosTheta + 2.87f;
        float b2 = 1.86f * cosTheta - 2.62f;
        s.bem = b1 * s.eleaf + b2;
        s.aem = s.eleaf + 0.025f - s.bem - s.emsoil;
        UpdateFc(s);

        // scale wind from refhtw to zta over grass
        float z0g = 0.005f;
        float dispg = 0.0f;
        float factw1 = (MathF.Log(s.zta - dispg) - MathF.Log(z0g)) /
                       (MathF.Log(s.refhtw - dispg) - MathF.Log(z0g));
        s.w1 *= factw1;
        s.w2 *= factw1;
        s.refhtw = s.zta;

        const float wmin = 3.0f, wmax = 20.0f;
        if (s.w1 < wmin) s.w1 = wmin;
        if (s.w2 < wmin) s.w2 = wmin;
        if (s.w1 > wmax) s.w1 = wmax;
        if (s.w2 > wmax) s.w2 = wmax;

        // volumetric heat capacity of air [J/deg-m3]
        float rho1 = s.pres1 / (287.04f * (s.taobs1 + 273.15f)) * (1.0f - 0.378f * s.ea1 / s.pres1) * 100.0f;
        float rho2 = s.pres2 / (287.04f * (s.taobs2 + 273.15f)) * (1.0f - 0.378f * s.ea2 / s.pres2) * 100.0f;
        s.rhocp1 = rho1 * s.cp;
        s.rhocp2 = rho2 * s.cp;

        GetPblTable(s);
    }

    /// <summary>updatefc — emissivity and radiation properties from fc/xlai.</summary>
    public static void UpdateFc(PixelState s)
    {
        AlexiAtmos.CanopyArch(s);
        s.esfc = s.aem * s.fc * s.fc + s.bem * s.fc + s.emsoil;
        AlexiRad.GetRadProps(s, s.zen1, out s.albedo1, out s.taubtv1, out s.taubtn1, s.clumps1);
        AlexiRad.GetRadProps(s, s.zen2, out s.albedo2, out s.taubtv2, out s.taubtn2, s.clumps2);
    }

    /// <summary>
    /// checkinput — range checks. The PBL profile loop checks are commented out
    /// in the source, so they are not ported. ibad is local (implicit integer).
    /// </summary>
    public static void CheckInput(PixelState s, TextWriter stdout, ref int ibad)
    {
        s.badinput = false;

        CheckValue("EA1", 3, s.ea1, 0.0f, 80.0f, s, stdout, ref ibad);
        CheckValue("EA2", 4, s.ea2, 0.0f, 80.0f, s, stdout, ref ibad);
        CheckValue("TAOBS1", 5, s.taobs1, -50.0f, 60.0f, s, stdout, ref ibad);
        CheckValue("TAOBS2", 6, s.taobs2, -50.0f, 60.0f, s, stdout, ref ibad);
        CheckValue("TRAD1", 7, s.trad1, -50.0f, 60.0f, s, stdout, ref ibad);
        CheckValue("TRAD2", 8, s.trad2, -50.0f, 60.0f, s, stdout, ref ibad);
        CheckValue("SDN1", 9, s.sdn1, 0.0f, 1500.0f, s, stdout, ref ibad);
        CheckValue("SDN2", 10, s.sdn2, 0.0f, 1500.0f, s, stdout, ref ibad);
        CheckValue("PRES1", 11, s.pres1, 200.0f, 1200.0f, s, stdout, ref ibad);
        CheckValue("PRES2", 12, s.pres2, 200.0f, 1200.0f, s, stdout, ref ibad);
        CheckValue("W1", 13, s.w1, 0.0f, 100.0f, s, stdout, ref ibad);
        CheckValue("W2", 14, s.w2, 0.0f, 100.0f, s, stdout, ref ibad);
        CheckValue("XLWDN1", 15, s.xlwdn1, 0.0f, 1000.0f, s, stdout, ref ibad);
        CheckValue("XLWDN2", 16, s.xlwdn2, 0.0f, 1000.0f, s, stdout, ref ibad);
        CheckValue("ZEN1", 17, s.zen1, 0.0f, 3.14f, s, stdout, ref ibad);
        CheckValue("ZEN2", 18, s.zen2, 0.0f, 3.14f, s, stdout, ref ibad);
        CheckValue("YEAR", 19, (float)s.year, 1900.0f, 2100.0f, s, stdout, ref ibad);
        CheckValue("DOY", 20, (float)s.doy, 0.0f, 366.0f, s, stdout, ref ibad);
        CheckValue("THETA", 21, s.theta, 0.0f, 80.0f, s, stdout, ref ibad);

        CheckValue("XLAI", 22, s.xlai, 0.0f, 10.0f, s, stdout, ref ibad);
        CheckValue("FG", 23, s.fg, 0.0f, 1.0f, s, stdout, ref ibad);
        CheckValue("HEIGHT", 25, s.height, 0.0f, 35.0f, s, stdout, ref ibad);
        CheckValue("CLUMP", 26, s.clump, 0.0f, 1.0f, s, stdout, ref ibad);
        CheckValue("XL", 27, s.xl, 0.0f, 0.10f, s, stdout, ref ibad);

        CheckValue("ALEAFV", 28, s.aleafv, 0.0f, 1.0f, s, stdout, ref ibad);
        CheckValue("ALEAFN", 29, s.aleafn, 0.0f, 1.0f, s, stdout, ref ibad);
        CheckValue("ALEAFL", 30, s.aleafl, 0.0f, 1.0f, s, stdout, ref ibad);
        CheckValue("RSOILV", 31, s.rsoilv, 0.0f, 1.0f, s, stdout, ref ibad);
        CheckValue("RSOILN", 32, s.rsoiln, 0.0f, 1.0f, s, stdout, ref ibad);
        CheckValue("EMSOIL", 33, s.emsoil, 0.0f, 1.0f, s, stdout, ref ibad);

        float dtrad = s.trad2 - s.trad1;
        CheckValue("DTRAD", 34, dtrad, 0.0f, 50.0f, s, stdout, ref ibad);
        CheckValue("TR:TA1", 35, s.trad1, s.taobs1 - 39.0f, 60.0f, s, stdout, ref ibad);
        CheckValue("TR:TA2", 36, s.trad2, s.taobs2 - 29.0f, 60.0f, s, stdout, ref ibad);
    }

    /// <summary>
    /// checkvalue — range check with Fortran FORMAT 5000 output
    /// (a6, ' = ', f9.2, ' outside range ', f9.2, ' to ', f9.2), InvariantCulture.
    /// badinput/writeme are arguments in Fortran (COMMON flags in the pixel context,
    /// locals in the US-flux context), so the core takes them explicitly.
    /// </summary>
    public static void CheckValue(string valname, int itag, float value, float xmin, float xmax,
                                  ref bool badinput, bool writeme, TextWriter stdout, ref int ibad)
    {
        if (value < xmin || value > xmax)
        {
            if (writeme)
            {
                stdout.Write(string.Concat(
                    valname.PadRight(6),
                    " = ", F92(value),
                    " outside range ", F92(xmin),
                    " to ", F92(xmax)));
                stdout.Write("\n");
            }
            badinput = true;
            ibad = itag;
        }
    }

    /// <summary>Pixel-context convenience: badinput/writeme come from common/flags.</summary>
    public static void CheckValue(string valname, int itag, float value, float xmin, float xmax,
                                  PixelState s, TextWriter stdout, ref int ibad)
    {
        CheckValue(valname, itag, value, xmin, xmax, ref s.badinput, s.writeme, stdout, ref ibad);
    }

    private static string F92(float v) =>
        v.ToString("F2", CultureInfo.InvariantCulture).PadLeft(9);

    /// <summary>
    /// ALEXI_errorcode — messages for ierr 0..3; blank lines for 4,5,6,7,8,10
    /// (bare write(6,*)); nothing for 9,11,12,13 (writes commented out in source).
    /// </summary>
    public static void AlexiErrorcode(PixelState s, int ierr, TextWriter stdout)
    {
        if (!s.writeme)
            return;
        switch (ierr)
        {
            case 0: stdout.WriteLine("*** CONVERGED ***"); break;
            case 1: stdout.WriteLine("*** HN did not converge. Bail ***"); break;
            case 2: stdout.WriteLine("*** Bad input.  Bail ***"); break;
            case 3: stdout.WriteLine("*** Exceeded max PBL profile layer. Bail ***"); break;
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
            case 10: stdout.WriteLine(); break;
            // 9, 11, 12, 13: commented out in source — no output.
        }
    }
}
