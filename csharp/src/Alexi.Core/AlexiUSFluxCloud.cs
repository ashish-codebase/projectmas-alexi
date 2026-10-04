namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_0.1/src/USflux_cloud.f — cloudy-pixel day processing.
/// </summary>
public static class AlexiUSFluxCloud
{
    /// <summary>
    /// cloudy_day_proc — hourly fluxes then cloudy daily-flux placeholder.
    /// Source calls `call hourly_flux(ibad)` against a zero-arg routine (arity defect);
    /// the argument is ignored.
    /// </summary>
    public static void CloudyDayProc(PixelState s, DailyState d, TextWriter stdout)
    {
        AlexiUSFluxUtl.HourlyFlux(s, d, stdout);
        if (s.badinput) return;
        DailyFluxCloudy(d);
    }

    /// <summary>daily_flux_cloudy — cloudy pixels get BAD daily fluxes (gap-fill downstream).</summary>
    public static void DailyFluxCloudy(DailyState d)
    {
        d.eday = GridDims.Bad;
        d.hday = GridDims.Bad;
        d.ecday = GridDims.Bad;
        d.esday = GridDims.Bad;
        d.hcday = GridDims.Bad;
        d.hsday = GridDims.Bad;
    }
}
