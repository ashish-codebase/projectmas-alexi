namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_0.1/src/USflux_rad.f — SDN partitioning (Cupid RADIN4 strategy)
/// and the simple net-radiation estimate.
/// </summary>
public static class AlexiUSFluxRad
{
    /// <summary>
    /// getradcomps — partition SDN into VIS/NIR, diffuse/direct. sdn is passed by
    /// reference in Fortran and is zeroed on the night branch; fclear is returned.
    /// Sets common/partition in PixelState.
    /// </summary>
    public static float GetRadComps(PixelState s, float zen, ref float sdn)
    {
        float pid2 = s.pi / 2.0f;
        float pid180 = s.pi / 180.0f;
        float ration = 1.0f; // assume clear at night absent other info
        float fclear;

        float coszen = MathF.Cos(zen);
        if (coszen < 0.01f)
        {
            // nighttime
            fclear = ration;
            s.fvis = 0.5f; s.fnir = 0.5f;
            s.difvis = 1.0f; s.difnir = 1.0f;
            s.dirvis = 0.0f; s.dirnir = 0.0f;
            if (sdn <= 1.0f) sdn = 0.0f;
            if (sdn != 0.0f)
                sdn = 0.0f; // "The sun is shining at night!" write is commented out
            return fclear;
        }

        // potential (clear-sky) visible and NIR components
        float airmas = (MathF.Sqrt(coszen * coszen + 0.0025f) - coszen) / 0.00125f; // atmospheric curvature
        airmas = airmas - 2.8f / MathF.Pow(90.0f - zen / pid180, 2.0f); // refraction correction

        float potbm1 = 600.0f * MathF.Exp(-0.160f * airmas);
        float potvis = (potbm1 + (600.0f - potbm1) * 0.4f) * coszen;
        float potdif = (600.0f - potbm1) * 0.4f * coszen;
        float u = 1.0f / coszen;
        float axlog = MathF.Log10(u);
        float a = MathF.Pow(10.0f, -1.195f + 0.4459f * axlog - 0.0345f * axlog * axlog);
        float watabs = 1320.0f * a;
        float potbm2 = 720.0f * MathF.Exp(-0.05f * airmas) - watabs;
        if (potbm2 < 0.0f) potbm2 = 0.0f;
        float eval = (720.0f - potbm2 - watabs) * 0.54f * coszen;
        float potnir = eval + potbm2 * coszen;

        fclear = sdn / (potvis + potnir);
        fclear = MathF.Min(1.0f, fclear);


        // partition SDN into VIS and NIR
        s.fvis = potvis / (potvis + potnir);
        s.fnir = potnir / (potvis + potnir);

        // direct-beam and diffuse fractions per waveband
        float fb1 = potbm1 * coszen / potvis;
        float fb2 = potbm2 * coszen / potnir;
        float ratiox = fclear;
        if (fclear > 0.9f) ratiox = 0.9f;
        s.dirvis = fb1 * (1.0f - MathF.Pow((0.9f - ratiox) / 0.7f, 0.6667f));
        if (fclear > 0.88f) ratiox = 0.88f;
        s.dirnir = fb2 * (1.0f - MathF.Pow((0.88f - ratiox) / 0.68f, 0.6667f));
        s.dirvis = MathF.Max(0.0f, s.dirvis);
        s.dirnir = MathF.Max(0.0f, s.dirnir);
        s.dirvis = MathF.Min(fb1, s.dirvis);
        s.dirnir = MathF.Min(fb2, s.dirnir);
        if (s.dirvis < 0.01f && s.dirnir > 0.01f) s.dirvis = 0.011f;
        if (s.dirnir < 0.01f && s.dirvis > 0.01f) s.dirnir = 0.011f;
        s.difvis = 1.0f - s.dirvis;
        s.difnir = 1.0f - s.dirnir;
        return fclear;
    }

    /// <summary>
    /// getnetrad_simple — net radiation from SDN + air temperature only.
    /// (zen is an argument in the source but unused; kept for signature parity.)
    /// </summary>
    public static float GetNetRadSimple(float sdn, float ta, float fclear, float zen, float albedo)
    {
        const float sigma = 5.67e-8f;
        float tak = ta + 273.15f;
        float eskyc = 9.2e-6f * tak * tak; // clear: Swinbank '63
        float esky = (1.0f - 0.84f * fclear) * eskyc + 0.84f * fclear; // cloudy: Monteith & Unsworth '90
        const float esfc = 0.98f;

        float rnetl = esfc * (esky - 1.0f) * sigma * MathF.Pow(tak, 4.0f);
        rnetl = esfc * esky * sigma * MathF.Pow(tak, 4.0f) - esfc * sigma * MathF.Pow(tak + 4.0f, 4.0f);
        float rnets = sdn * (1.0f - albedo);
        return rnetl + rnets;
    }
}
