namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_0.1/src/USflux_clear.f — clear-sky day processing:
/// hourly flux integration plus daily partitioning of system/soil/canopy
/// fluxes. Field names verbatim.
///
/// Note: `fawsfc`/`fawrz` in daily_flux_clear_SDN are implicit locals in the
/// source (not the common/availh2o awfsfc/awfrz), assigned BAD and never
/// used; not reproduced.
/// </summary>
public static class AlexiUSFluxClear
{
    /// <summary>clear_day_proc(ibad) — hourly flux, then daily partitioning.</summary>
    public static void ClearDayProc(PixelState s, DailyState d, TextWriter stdout,
        ref int ibad)
    {
        // hourly_flux(ibad) — arity defect: hourly_flux takes no args in the
        // source; ibad is ignored (see AlexiUSFluxUtl header).
        AlexiUSFluxUtl.HourlyFlux(s, d, stdout);
        if (s.badinput)
            return;
        DailyFluxClearSDN(s, d, stdout);
    }

    /// <summary>
    /// daily_flux_clear_SDN — partition daily fluxes using the SDN fraction
    /// of the day (fsdn).
    /// </summary>
    public static void DailyFluxClearSDN(PixelState s, DailyState d, TextWriter stdout)
    {
        // System fluxes
        d.fsdn = s.sdn2 / d.sday;
        d.eday = s.xle2 / d.fsdn;
        d.gday = 0.0f; // assume soil heat flux integrates to 0
        d.hday = d.rnday - d.eday - d.gday;

        // Soil fluxes
        if (s.xle2 != 0.0)
            d.esday = (s.xles2 / s.xle2) * d.eday;
        else
            d.esday = 0.0f;
        d.hsday = d.rnsday - d.gday - d.esday;

        // Canopy flux as residual
        d.ecday = d.eday - d.esday;
        d.hcday = d.rncday - d.ecday;

        // Total evapotranspiration factor
        AlexiUSFluxUtl.Faopm(s, out s.eref2, s.tloc2, s.sdn2, s.zen2,
            s.taobs2, s.pres2, s.ea2, s.w2orig);
        d.fpet = s.xle2 / s.eref2;
        if (d.fpet < 0.0f)
        {
            if (s.writeme)
                stdout.WriteLine("FPET = " + d.fpet + " - set to 0");
            d.fpet = 0.0f;
        }

        // fawsfc = BAD; fawrz = BAD — unused locals in the source.
    }

    /// <summary>
    /// daily_flux_clear_EF — OLD, UNSUPPORTED variant (kept for parity with
    /// the source; not called by the 0.1 driver).
    ///
    /// Note: ft/fvpd are uninitialized implicit locals in the source (the
    /// getstress call is commented out); the port uses 0f and documents it.
    /// epot2/espot2/ecpot2/fawrz/fawsfc are unused implicit locals.
    /// </summary>
    public static void DailyFluxClearEF(PixelState s, DailyState d, TextWriter stdout)
    {
        // System fluxes
        float fevap = 1.1f * s.xle2 / (s.rnet2 - s.g2);
        d.eday = fevap * (d.rnday - d.gday);
        d.hday = d.rnday - d.gday - d.eday;
        stdout.WriteLine("EDAY HDAY RNDAY GDAY " + d.eday + " " + d.hday + " " + d.rnday + " " + d.gday);

        // Soil fluxes
        float fevaps = 1.1f * s.xles2 / (s.rnsoil2 - s.g2);
        d.esday = fevaps * (d.rnsday - d.gday);
        d.hsday = d.rnsday - d.gday - d.esday;

        // Canopy flux as residual
        d.ecday = d.eday - d.esday;
        d.hcday = d.rncday - d.ecday;

        float xlam = (2.501f - 0.00237f * s.ta2) * 1e6f;
        float esat = 6.108f * MathF.Pow(10.0f, 7.5f * s.ta2 / (237.3f + s.ta2));
        float s_ = 2.17e-3f * xlam * esat / ((273.15f + s.ta2) * (273.15f + s.ta2));

        // Potential canopy transpiration [W/m2] — local ecpot2, unused after
        float ecpot2 = s.rndiv2 * 1.3f * s.fg * (s_ / (s_ + 0.66f));

        // Potential soil evaporation [W/m2] — local espot2, unused after
        float tauc = 0.5f;
        const float alpha = 1.3f;
        float tau = MathF.Exp(-0.45f * s.xlai / MathF.Sqrt(2.0f * MathF.Cos(s.zen2)));
        float pts;
        if (tau <= tauc)
            pts = 1.0f;
        else
            pts = alpha - ((alpha - 1.0f) * (1.0f - tau)) / (1.0f - tauc);
        float espot2 = (s.rnsoil2 - d.gday) * pts * (s_ / (s_ + 0.66f));

        // System ET — local epot2, unused after
        float epot2 = ecpot2 + espot2;

        // Vegetation stress factors — fawrz/fawsfc are unused locals in the
        // source; ft/fvpd are uninitialized (getstress commented out).
        float ft = 0f;
        float fvpd = 0f;
        float fawrz = d.ecday / (d.ecpotday * ft * fvpd);
        if (fawrz > 1.0f)
        {
            if (s.writeme)
                stdout.WriteLine("FAWRZ = " + fawrz + " - NOT set to 1");
        }
        if (fawrz < 0.0f)
        {
            if (s.writeme)
                stdout.WriteLine("FAWRZ = " + fawrz + " - set to 0");
            fawrz = 0.0f;
        }
        float fawsfc = d.esday / d.espotday;
        if (fawsfc > 1.0f)
        {
            if (s.writeme)
                stdout.WriteLine("FAWSFC = " + fawsfc + " - NOT set to 1");
        }
        if (fawsfc < 0.0f)
        {
            if (s.writeme)
                stdout.WriteLine("FAWSFC = " + fawsfc + " - set to 0");
            fawsfc = 0.0f;
        }

        // Total evapotranspiration factor (common/moisture + moisture2)
        d.faw = d.eday / d.epotday;
        d.faw2 = s.xle2 / epot2;
        if (d.faw > 1.0f)
        {
            if (s.writeme)
                stdout.WriteLine("FAW = " + d.faw + " - NOT set to 1");
        }
        if (d.faw < 0.0f)
        {
            if (s.writeme)
                stdout.WriteLine("FAW = " + d.faw + " - set to 0");
            d.faw = 0.0f;
        }
    }
}
