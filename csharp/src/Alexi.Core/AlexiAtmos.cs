namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_0.1/src/ALEXI_atmos.f — aerodynamic resistance III (Massman &amp; Weil 1996).
/// Routine and field names kept verbatim. All reals are float (real*4).
/// </summary>
public static class AlexiAtmos
{
    /// <summary>
    /// getresistance — aerodynamic (RA), soil-surface (RS) and 1-sided leaf boundary-layer
    /// (RX) resistances, s/m on a per-unit-ground-area basis.
    /// Reads/updates common/aerodyn + common/stability state in PixelState.
    /// </summary>
    public static void GetResistance(PixelState s, float h, float wind, float ta, float ts,
        float tc, float rhocp, out float ra, out float rs, out float rx)
    {
        // First, use stability from previous iteration
        float ustarb = 0.4f * wind / (s.xlog1 - s.psima);
        float zdla = -(s.refhtw - s.disp) * h * 9.8f * 0.4f /
                     (rhocp * (ta + 273.0f) * ustarb * ustarb * ustarb); // note: 273., not 273.15 (verbatim)
        zdla = MathF.Min(zdla, s.zdlamx);
        zdla = MathF.Max(zdla, s.zdlamn);
        PsimhnKy(zdla, out s.psima, out s.psih);

        // Now recalculate stability corrected ustar
        float ustara = 0.4f * wind / (s.xlog1 - s.psima);
        zdla = zdla * ustarb * ustarb * ustarb / (ustara * ustara * ustara);
        zdla = MathF.Min(zdla, s.zdlamx);
        zdla = MathF.Max(zdla, s.zdlamn);
        PsimhnKy(zdla, out s.psima, out s.psih);
        float ustar = 0.4f * wind / (s.xlog1 - s.psima);

        // ---- bare soil or water ----
        if (s.iswater ||
            (!s.perennial && (s.fc <= s.fcbare || s.height == 0.0f)))
        {
            float tas = 273.15f + (ta + ts) * 0.5f;
            const float xmuo = 1.8325e-5f;
            float xmu = (296.16f + 393.16f) / (tas + 393.16f) *
                        MathF.Pow(tas / 296.16f, 1.5f) * xmuo;
            float rho = rhocp / s.cp;
            float xnu = xmu / rho;
            float restar = ustar * s.z0 / xnu;
            float z0h = s.z0 * MathF.Exp(-s.xk * (4.0f * MathF.Pow(restar, 0.15f) - 5.0f));
            ra = (MathF.Log((s.refhtw - s.disp) / z0h) - s.psih) / (ustar * s.xk);
            rs = 0.0f;
            rx = 10000.0f;
            return;
        }

        // ---- all other cases ----
        // Ra
        ra = (s.xlog1 - s.psima) * (s.xlog1 - s.psih) / (0.16f * wind);

        // Rs (source uses ts-ta; ts-tc is commented out: "TC gets crazy at low fc")
        float tgrad = ts - ta;
        const float zdlamx2 = 0.4f, zdlamn2 = -0.5f;
        float zdla2 = zdla * (s.height - s.disp) / (s.refhtw - s.disp);
        zdla2 = MathF.Min(zdla2, zdlamx2);
        zdla2 = MathF.Max(zdla2, zdlamn2);
        PsimhnKy(zdla2, out float psima2, out float psih2); // locals, not the common state
        float uc = wind * ((s.xlog2 - psima2) / (s.xlog1 - s.psima));
        float us = uc * s.expuxp;
        // xndvi-dependent val1/val2 branches are commented out in source; live path:
        float val1 = 0.0025f, val2 = 0.012f;
        if (tgrad > 1.0f)
            rs = 1.0f / (val1 * MathF.Pow(tgrad, 0.33f) + val2 * us);
        else
            rs = 1.0f / (val1 + val2 * us);
        // rs = rs*rsmin and the rsmin debug write are commented out in source.

        // Rx (boundary-layer resistance at z0+disp)
        float udz = uc * MathF.Exp(s.uexp2);
        rx = (180.0f * MathF.Sqrt(s.xl / udz)) / s.xlai;
        // rx = rsmin/xlai is commented out in source.
    }

    /// <summary>psimhn_ky — stability functions, Kader &amp; Yaglom (1990) form.</summary>
    public static void PsimhnKy(float zdla, out float psima, out float psih)
    {
        float y = -zdla;
        const float y0 = 0.0f;
        if (y < 0.0059f)
        {
            psih = 0.0f;
            psima = 0.0f;
        }
        else
        {
            psih = 1.2f * MathF.Log((0.33f + MathF.Pow(y, 0.78f)) / 0.33f);
            y = MathF.Min(y, 15.025f);
            psima = 1.47f * MathF.Log((0.28f + MathF.Pow(y, 0.75f)) /
                        (0.28f + MathF.Pow(0.0059f + y0, 0.75f)))
                    - 1.29f * (MathF.Pow(y, 0.33f) - MathF.Pow(0.0059f + y0, 0.33f));
        }
    }

    /// <summary>psimhn — stability functions, simplest Cupid scheme.</summary>
    public static void Psimhn(float zdla, out float psima, out float psih)
    {
        if (zdla >= 0.0f)
        {
            psima = -5.0f * zdla;
            psih = -5.0f * zdla;
        }
        else
        {
            psima = 1.88f + (1.0f / (-0.533f + 0.790f * zdla));
            psih = MathF.Exp(0.598f + 0.390f * MathF.Log(-zdla) -
                             0.09f * MathF.Pow(MathF.Log(-zdla), 2.0f));
        }
    }

    /// <summary>
    /// canopyarch — update canopy architecture factors derived from FC.
    /// disp/z0 come from landcover.f (Massman eqs are commented out here);
    /// this routine computes ftheta, wind-extinction factors, and xlog1/xlog2.
    /// </summary>
    public static void CanopyArch(PixelState s)
    {
        const float cd = 0.20f;
        const float z0s = 0.005f;    // bare-soil roughness
        const float z0w = 0.00035f;  // water roughness

        // cover-fraction variables (xlai comes from landcover, not recomputed here)
        s.ftheta = 1.0f - MathF.Exp(-0.5f * s.xlai * s.clump / MathF.Cos(s.theta));
        if (s.ftheta > 0.8f) s.ftheta = 0.8f;

        if (s.perennial)
        {
            // rigid vegetation: z0/disp held fixed (set by landcover) — no-op branch
        }
        else if (s.iswater)
        {
            s.z0 = z0w;
            s.disp = 0.0f;
            s.xlog1 = MathF.Log(s.refhtw / s.z0);
            return; // go to 500: skip wind-extinction block
        }
        else if (s.fc <= s.fcbare || s.height == 0.0f)
        {
            // bare soil: z0=z0s is commented out in source — z0 stays as set upstream
            s.xlog1 = MathF.Log(s.refhtw / s.z0);
            return; // go to 500
        }
        else
        {
            // Massman disp/z0 equations are commented out (done in landcover.f); live path:
            s.z0 = MathF.Max(s.z0, z0s);
            // for positive xlog2, need height-disp > z0
            s.height = MathF.Max(s.height, s.disp + s.z0 + 0.001f);
        }

        // wind extinction factor, A (Goudriaan '77 p.110)
        float xld = s.xlai / s.height;
        float xlm = MathF.Sqrt(4.0f * s.xl / (s.pi * xld));
        s.a = MathF.Sqrt(cd * s.clump * s.xlai * s.height / xlm); // assumes iw=0.5

        // cached factors for getresistance
        s.xlog1 = MathF.Log((s.refhtw - s.disp) / s.z0);
        s.xlog2 = MathF.Log((s.height - s.disp) / s.z0);
        if (s.height > 0.5f)
            s.uexp1 = -s.a * (1.0f - 0.05f / s.height);
        else
            s.uexp1 = -s.a * 0.90f;
        s.expuxp = MathF.Exp(s.uexp1);
        if (s.expuxp > 0.95f) s.expuxp = 0.95f;
        s.a = MathF.Sqrt(cd * (s.xlai / s.fveg) * s.height / xlm); // assumes iw=0.5
        s.uexp2 = -s.a * (1.0f - (s.z0 + s.disp) / s.height);
    }
}
