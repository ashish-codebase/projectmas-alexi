using System.Globalization;

namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_0.1/src/ALEXI_rad.f. Field and routine names kept verbatim.
/// All reals are float (real*4). COMMON-block state lives in PixelState/DailyState.
/// </summary>
public static class AlexiRad
{
    // parameter(sigma=5.67e-8) in getnetrad / getnetradnight / getxlwdn
    const float sigma = 5.67e-8f;

    /// <summary>
    /// getsunzen — sun zenith angle [rad] at local time tloc.
    /// Side effect: sets common/sun strt/end (DailyState), exactly as the Fortran.
    /// </summary>
    public static float GetSunzen(PixelState s, DailyState d, float tloc)
    {
        float pid180 = s.pi / 180.0f;
        float pid2 = s.pi / 2.0f;

        // latitude
        float sinlat = MathF.Sin(s.xlat * pid180);
        float coslat = MathF.Cos(s.xlat * pid180);

        // declination
        float kday = (s.year - 1977) * 365.0f + s.doy + 28123.0f;
        float xm = (-1.0f + 0.9856f * kday) * pid180;
        float delnu = 2.0f * 0.01674f * MathF.Sin(xm) + 1.25f * 0.01674f * 0.01674f * MathF.Sin(2.0f * xm);
        float slong = (-79.8280f + 0.9856479f * kday) * pid180 + delnu;
        float decmax = MathF.Sin(23.44f * pid180);
        float decl = MathF.Asin(decmax * MathF.Sin(slong));
        float sindec = MathF.Sin(decl);
        float cosdec = MathF.Cos(decl);

        float hfday = 12.0f / s.pi * MathF.Acos(-(sinlat * sindec) / (coslat * cosdec));
        float eqtm = 9.4564f * MathF.Sin(2.0f * slong) / cosdec - 4.0f * delnu / pid180;
        eqtm = eqtm / 60.0f;

        // longitude corrections
        float dlong = s.xlong - s.stdlng;
        if (dlong > 180.0f) dlong -= 360.0f;
        if (dlong < -180.0f) dlong += 360.0f;
        dlong = dlong / 15.0f; // positive E of stdlng

        d.strt = 12.0f - dlong - eqtm - hfday; // true sunrise local time
        d.end = 12.0f - dlong - eqtm + hfday;

        // sun zenith angle
        float timsun = tloc + eqtm + dlong;
        float hrang = (timsun - 12.0f) * pid2 / 6.0f;
        return MathF.Acos(sinlat * sindec + coslat * cosdec * MathF.Cos(hrang));
    }

    /// <summary>
    /// getradprops — canopy albedo/transmission from Goudriaan 1988 analytical solutions
    /// (Campbell &amp; Norman 1998). Zen-dependent outputs via out args; the saved COMMON
    /// state (taudv/taudn/taudl, albv/albn, emcpy, rdcpyl) is written into PixelState.
    /// </summary>
    public static void GetRadProps(PixelState s, float zen,
        out float albedo, out float taubtv, out float taubtn, float clumps)
    {
        float coszen = MathF.Cos(zen);

        // weighted live/dead leaf average absorptivities
        float ameanv = s.aleafv * s.fg + s.adeadv * (1.0f - s.fg);
        float ameann = s.aleafn * s.fg + s.adeadn * (1.0f - s.fg);
        float ameanl = s.aleafl * s.fg + s.adeadl * (1.0f - s.fg);

        // ---- diffuse components ----
        float akd = -0.0683f * MathF.Log(clumps * s.xlai) + 0.804f; // fit to Fig 15.4, x=1
        float rcpyn = (1.0f - MathF.Sqrt(ameann)) / (1.0f + MathF.Sqrt(ameann)); // Eq 15.7
        float rcpyv = (1.0f - MathF.Sqrt(ameanv)) / (1.0f + MathF.Sqrt(ameanv));
        float rcpyl = (1.0f - MathF.Sqrt(ameanl)) / (1.0f + MathF.Sqrt(ameanl));
        float rdcpyn = 2.0f * akd * rcpyn / (akd + 1.0f); // Eq 15.8
        float rdcpyv = 2.0f * akd * rcpyv / (akd + 1.0f);
        float rdcpyl = 2.0f * akd * rcpyl / (akd + 1.0f);
        s.rdcpyl = rdcpyl; // common/reflection

        // diffuse canopy transmission, visible (Eq 15.11)
        float expfac = MathF.Sqrt(ameanv) * akd * clumps * s.xlai;
        float xnum = (rdcpyv * rdcpyv - 1.0f) * MathF.Exp(-expfac);
        float xden = (rdcpyv * s.rsoilv - 1.0f) + rdcpyv * (rdcpyv - s.rsoilv) * MathF.Exp(-2.0f * expfac);
        s.taudv = xnum / xden;

        // diffuse canopy transmission, NIR
        expfac = MathF.Sqrt(ameann) * akd * clumps * s.xlai;
        xnum = (rdcpyn * rdcpyn - 1.0f) * MathF.Exp(-expfac);
        xden = (rdcpyn * s.rsoiln - 1.0f) + rdcpyn * (rdcpyn - s.rsoiln) * MathF.Exp(-2.0f * expfac);
        s.taudn = xnum / xden;

        // diffuse canopy transmission, longwave (deep-canopy expression, Eq 15.6)
        s.taudl = MathF.Exp(-MathF.Sqrt(ameanl) * akd * clumps * s.xlai);
        s.emcpy = 0.99f * (1.0f - s.taudl); // Bill K's emcpy

        // diffuse surface albedo (Eq 15.9)
        float fact = ((rdcpyn - s.rsoiln) / (rdcpyn * s.rsoiln - 1.0f)) *
                     MathF.Exp(-2.0f * MathF.Sqrt(ameann) * akd * clumps * s.xlai);
        float albdn = (rdcpyn + fact) / (1.0f + rdcpyn * fact);
        fact = ((rdcpyv - s.rsoilv) / (rdcpyv * s.rsoilv - 1.0f)) *
               MathF.Exp(-2.0f * MathF.Sqrt(ameanv) * akd * clumps * s.xlai);
        float albdv = (rdcpyv + fact) / (1.0f + rdcpyv * fact);

        // ---- beam components ----
        float akb = 0.5f / coszen;
        // rcpyn/rcpyv recomputed in source (identical values); reused for beam
        float rbcpyn = 2.0f * akb * rcpyn / (akb + 1.0f); // Eq 15.8
        float rbcpyv = 2.0f * akb * rcpyv / (akb + 1.0f);

        // beam surface albedo (Eq 15.9)
        fact = ((rbcpyn - s.rsoiln) / (rbcpyn * s.rsoiln - 1.0f)) *
               MathF.Exp(-2.0f * MathF.Sqrt(ameann) * akb * clumps * s.xlai);
        float albbn = (rbcpyn + fact) / (1.0f + rbcpyn * fact);
        fact = ((rbcpyv - s.rsoilv) / (rbcpyv * s.rsoilv - 1.0f)) *
               MathF.Exp(-2.0f * MathF.Sqrt(ameanv) * akb * clumps * s.xlai);
        float albbv = (rbcpyv + fact) / (1.0f + rbcpyv * fact);

        // weighted-average albedo
        albedo = s.fvis * (s.dirvis * albbv + s.difvis * albdv) +
                 s.fnir * (s.dirnir * albbn + s.difnir * albdn);
        s.albv = s.dirvis * albbv + s.difvis * albdv; // common/reflection2
        s.albn = s.dirnir * albbn + s.difnir * albdn;

        // beam+scattered canopy transmission, visible (Eq 15.11)
        expfac = MathF.Sqrt(ameanv) * akb * clumps * s.xlai;
        xnum = (rbcpyv * rbcpyv - 1.0f) * MathF.Exp(-expfac);
        xden = (rbcpyv * s.rsoilv - 1.0f) + rbcpyv * (rbcpyv - s.rsoilv) * MathF.Exp(-2.0f * expfac);
        taubtv = xnum / xden;

        // beam+scattered canopy transmission, NIR
        expfac = MathF.Sqrt(ameann) * akb * clumps * s.xlai;
        xnum = (rbcpyn * rbcpyn - 1.0f) * MathF.Exp(-expfac);
        xden = (rbcpyn * s.rsoiln - 1.0f) + rbcpyn * (rbcpyn - s.rsoiln) * MathF.Exp(-2.0f * expfac);
        taubtn = xnum / xden;
    }

    /// <summary>getnetrad — net radiation above canopy, above soil, and canopy divergence (day).</summary>
    public static void GetNetRad(PixelState s, float ts, float tc, float xlwdn, float sdn,
        float albedo, float taubtv, float taubtn,
        out float rnet, out float rnsoil, out float rndiv, out float swup, out float xlwup)
    {
        // rsoill = 1 - emsoil is computed in source but never used
        float tausolar = s.fvis * (s.difvis * s.taudv + s.dirvis * taubtv) +
                         s.fnir * (s.difnir * s.taudn + s.dirnir * taubtn);
        float tauthermal = s.taudl;
        float albedocanopy = albedo;
        float albedosoil = s.fvis * s.rsoilv + s.fnir * s.rsoiln;

        float rsky = xlwdn;
        float rcpy = 0.99f * sigma * MathF.Pow(tc + 273.15f, 4.0f); // BK's equations
        float rsoil = s.emsoil * sigma * MathF.Pow(ts + 273.15f, 4.0f);

        rnsoil = tauthermal * rsky + (1.0f - tauthermal) * rcpy - rsoil +
                 tausolar * (1.0f - albedosoil) * sdn;
        rndiv = (1.0f - tauthermal) * (rsky + rsoil - 2.0f * rcpy) +
                (1.0f - tausolar) * (1.0f - albedocanopy) * sdn;
        rnet = rndiv + rnsoil;

        // upwelling components
        float rnsoillw = tauthermal * rsky + (1.0f - tauthermal) * rcpy - rsoil;
        float rndivlw = (1.0f - tauthermal) * (rsky + rsoil - 2.0f * rcpy);
        float rnlw = rnsoillw + rndivlw;
        xlwup = xlwdn - rnlw;

        float rnsoilsw = tausolar * (1.0f - albedosoil) * sdn;
        float rndivsw = (1.0f - tausolar) * (1.0f - albedocanopy) * sdn;
        float rnsw = rnsoilsw + rndivsw;
        swup = sdn - rnsw;

        xlwup = sdn - swup + xlwdn - rnet; // source recomputes xlwup; last assignment wins
    }

    /// <summary>getnetradnight — same as getnetrad with all SW terms zeroed.</summary>
    public static void GetNetRadNight(PixelState s, float ts, float tc, float xlwdn, float sdn,
        float albedo, float taubtv, float taubtn,
        out float rnet, out float rnsoil, out float rndiv, out float swup, out float xlwup)
    {
        float tausolar = s.fvis * (s.difvis * s.taudv + s.dirvis * taubtv) +
                         s.fnir * (s.difnir * s.taudn + s.dirnir * taubtn);
        float tauthermal = s.taudl;
        float albedocanopy = albedo;
        float albedosoil = s.fvis * s.rsoilv + s.fnir * s.rsoiln;

        float rsky = xlwdn;
        float rcpy = 0.99f * sigma * MathF.Pow(tc + 273.15f, 4.0f);
        float rsoil = s.emsoil * sigma * MathF.Pow(ts + 273.15f, 4.0f);

        rnsoil = tauthermal * rsky + (1.0f - tauthermal) * rcpy - rsoil;
        rndiv = (1.0f - tauthermal) * (rsky + rsoil - 2.0f * rcpy);
        rnet = rndiv + rnsoil;

        float rnsoillw = tauthermal * rsky + (1.0f - tauthermal) * rcpy - rsoil;
        float rndivlw = (1.0f - tauthermal) * (rsky + rsoil - 2.0f * rcpy);
        float rnlw = rnsoillw + rndivlw;
        xlwup = xlwdn - rnlw;

        float rnsw = 0.0f; // rnsoilsw = 0, rndivsw = 0
        swup = sdn - rnsw;

        xlwup = sdn - swup + xlwdn - rnet;
    }

    /// <summary>getxlwdn — thermal radiation from sky (Brutsaert Eq.), clear-sky/cloud weighting.</summary>
    public static void GetXlwdn(float ea, float ta, float fclear, out float xlwdn)
    {
        float tak = ta + 273.15f;
        float esky = 1.24f * MathF.Pow(ea / tak, 1.0f / 7.0f);
        xlwdn = sigma * MathF.Pow(tak, 4.0f) * (esky * fclear + 1.0f - fclear);
    }

    /// <summary>
    /// find_albedo_soil — iteratively adjust rsoilv/rsoiln so computed albv/albn match
    /// observed albvobs/albnobs (Newton-like nudge, max 20 iterations, BAD on failure).
    /// Side effects: sets partition factors, common/sun strt/end, s.zen2, rsoilv/rsoiln.
    /// </summary>
    public static void FindAlbedoSoil(PixelState s, DailyState d, float tloc2, TextWriter stdout)
    {
        const float errmax = 0.002f;

        // assume clear-sky partitioning
        s.fvis = 0.5f; s.fnir = 0.5f;
        s.dirvis = 0.8f; s.difvis = 0.2f;
        s.dirnir = 1.0f; s.difnir = 0.0f;

        s.zen2 = GetSunzen(s, d, tloc2);

        // rsoilv loop
        bool done = false;
        int iter = 1;
        float alast = -9999.0f;
        float deltr = 0.0f;
        while (!done)
        {
            GetRadProps(s, s.zen2, out _, out _, out _, s.clumps2);
            WriteAlbedoLine(stdout, "RSOILV: ", iter, s.rsoilv, s.albv, s.albvobs);
            float adiff = s.albvobs - s.albv;
            if (MathF.Abs(adiff) > errmax)
            {
                float delta = s.albv - alast;
                if (iter == 1)
                {
                    deltr = 0.02f * MathF.Abs(adiff) / adiff; // change by 0.02 with sign of adiff
                }
                else if (iter < 20)
                {
                    float dadr = delta / deltr;
                    float dr = adiff / dadr;
                    deltr = 0.75f * dr; // nudge by 0.75*DR
                }
                else
                {
                    s.rsoilv = GridDims.Bad;
                    done = true;
                }
                if (!done)
                {
                    s.rsoilv += deltr;
                    alast = s.albv;
                    iter += 1;
                }
            }
            else done = true;
        }

        // rsoiln loop
        done = false;
        iter = 1;
        alast = -9999.0f;
        while (!done)
        {
            GetRadProps(s, s.zen2, out _, out _, out _, s.clumps2);
            WriteAlbedoLine(stdout, "RSOILN: ", iter, s.rsoiln, s.albn, s.albnobs);
            float adiff = s.albnobs - s.albn;
            if (MathF.Abs(adiff) > errmax)
            {
                float delta = s.albn - alast;
                if (iter == 1)
                {
                    deltr = 0.02f * MathF.Abs(adiff) / adiff;
                }
                else if (iter < 20)
                {
                    float dadr = delta / deltr;
                    float dr = adiff / dadr;
                    deltr = 0.75f * dr;
                }
                else
                {
                    s.rsoiln = GridDims.Bad;
                    done = true;
                }
                if (!done)
                {
                    s.rsoiln += deltr;
                    alast = s.albn;
                    iter += 1;
                }
            }
            else done = true;
        }
    }

    /// <summary>
    /// find_albedo_veg — same loop as find_albedo_soil, adjusting aleafv/aleafn
    /// (nudge sign flipped: -0.02 on first iteration).
    /// </summary>
    public static void FindAlbedoVeg(PixelState s, DailyState d, float tloc2, TextWriter stdout)
    {
        const float errmax = 0.002f;

        s.fvis = 0.5f; s.fnir = 0.5f;
        s.dirvis = 0.8f; s.difvis = 0.2f;
        s.dirnir = 1.0f; s.difnir = 0.0f;

        s.zen2 = GetSunzen(s, d, tloc2);

        // aleafv loop
        bool done = false;
        int iter = 1;
        float alast = -9999.0f;
        float deltr = 0.0f;
        while (!done)
        {
            GetRadProps(s, s.zen2, out _, out _, out _, s.clumps2);
            WriteAlbedoLine(stdout, "ALEAFV: ", iter, s.aleafv, s.albv, s.albvobs);
            float adiff = s.albvobs - s.albv;
            if (MathF.Abs(adiff) > errmax)
            {
                float delta = s.albv - alast;
                if (iter == 1)
                {
                    deltr = -0.02f * MathF.Abs(adiff) / adiff; // change by 0.02 with - sign of adiff
                }
                else if (iter < 20)
                {
                    float dadr = delta / deltr;
                    float dr = adiff / dadr;
                    deltr = 0.75f * dr;
                }
                else
                {
                    s.aleafv = GridDims.Bad;
                    done = true;
                }
                if (!done)
                {
                    s.aleafv += deltr;
                    alast = s.albv;
                    iter += 1;
                }
            }
            else done = true;
        }

        // aleafn loop
        done = false;
        iter = 1;
        alast = -9999.0f;
        while (!done)
        {
            GetRadProps(s, s.zen2, out _, out _, out _, s.clumps2);
            WriteAlbedoLine(stdout, "ALEAFN: ", iter, s.aleafn, s.albn, s.albnobs);
            float adiff = s.albnobs - s.albn;
            if (MathF.Abs(adiff) > errmax)
            {
                float delta = s.albn - alast;
                if (iter == 1)
                {
                    deltr = -0.02f * MathF.Abs(adiff) / adiff;
                }
                else if (iter < 20)
                {
                    float dadr = delta / deltr;
                    float dr = adiff / dadr;
                    deltr = 0.75f * dr;
                }
                else
                {
                    s.aleafn = GridDims.Bad;
                    done = true;
                }
                if (!done)
                {
                    s.aleafn += deltr;
                    alast = s.albn;
                    iter += 1;
                }
            }
            else done = true;
        }
    }

    /// <summary>Fortran FORMAT 100: (a10,i5,3f11.5).</summary>
    static void WriteAlbedoLine(TextWriter stdout, string label, int iter,
        float v1, float v2, float v3)
    {
        string line = label.PadRight(10) + iter.ToString().PadLeft(5) +
                      v1.ToString("F5", CultureInfo.InvariantCulture).PadLeft(11) +
                      v2.ToString("F5", CultureInfo.InvariantCulture).PadLeft(11) +
                      v3.ToString("F5", CultureInfo.InvariantCulture).PadLeft(11);
        stdout.Write(line + "\n");
    }
}
