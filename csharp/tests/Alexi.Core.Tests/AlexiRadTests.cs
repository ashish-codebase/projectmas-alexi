using Alexi.Core;
using Xunit;

namespace Alexi.Core.Tests;

public class AlexiRadTests
{
    static PixelState MakeRadPixel()
    {
        var s = new PixelState();
        AlexiUtl.SetConstants(s); // sets pi, cp, xk
        // site: 40N, -105, standard meridian -105 (dlong = 0)
        s.xlat = 40f; s.xlong = -105f; s.stdlng = -105f;
        s.year = 2008; s.doy = 172; // summer solstice
        // grass-like canopy
        s.aleafv = 0.8f; s.aleafn = 0.2f; s.aleafl = 0.9f;
        s.adeadv = 0.5f; s.adeadn = 0.4f; s.adeadl = 0.9f;
        s.fg = 1f; s.xlai = 2f;
        s.clumps1 = 1f; s.clumps2 = 1f;
        s.rsoilv = 0.17f; s.rsoiln = 0.25f;
        s.emsoil = 0.92f;
        // clear-sky partitioning
        s.fvis = 0.5f; s.fnir = 0.5f;
        s.dirvis = 0.8f; s.difvis = 0.2f;
        s.dirnir = 1.0f; s.difnir = 0.0f;
        return s;
    }

    [Fact]
    public void GetSunzen_NoonSummerSolsticeMatchesLatitudeDeclination()
    {
        var s = MakeRadPixel();
        var d = new DailyState();
        float zen = AlexiRad.GetSunzen(s, d, 12.0f);
        // at local solar noon, zen = |lat - decl| ≈ |40 - 23.44| deg = 16.56 deg
        float deg = zen * 180f / MathF.PI;
        Assert.InRange(deg, 16.0, 17.1);
    }

    [Fact]
    public void GetSunzen_SetsSymmetricSunriseSunset()
    {
        var s = MakeRadPixel();
        var d = new DailyState();
        AlexiRad.GetSunzen(s, d, 12.0f);
        // dlong = 0, so (strt+end)/2 = 12 - eqtm; eqtm is minutes-scale (< 0.25 h)
        float noon = (d.strt + d.end) / 2f;
        Assert.InRange(noon, 11.75f, 12.25f);
        // day length at 40N summer solstice ≈ 14.9 h
        Assert.InRange(d.end - d.strt, 14.0f, 15.5f);
    }

    [Fact]
    public void GetRadProps_ProducesPhysicalRanges()
    {
        var s = MakeRadPixel();
        float zen = 0.5f; // ~28.6 deg
        AlexiRad.GetRadProps(s, zen, out float albedo, out float taubtv, out float taubtn, s.clumps1);

        Assert.InRange(albedo, 0.05f, 0.35f);   // grass albedo
        Assert.InRange(s.albv, 0.05f, 0.35f);
        Assert.InRange(s.albn, 0.05f, 0.5f);
        Assert.InRange(s.taudl, 0f, 1f);
        Assert.InRange(taubtv, -1f, 1f);
        Assert.InRange(taubtn, -1f, 1f);
        Assert.InRange(s.emcpy, 0f, 1f);
        // denser canopy -> lower LW transmission
        s.xlai = 4f;
        AlexiRad.GetRadProps(s, zen, out _, out _, out _, s.clumps1);
        Assert.True(s.taudl < 0.5f);
    }

    [Fact]
    public void GetXlwdn_HandCheck()
    {
        // ta is degrees C; Fortran adds 273.15 internally
        float ea = 1000f, ta = 15f, fclear = 1f;
        float tak = ta + 273.15f;
        float esky = 1.24f * MathF.Pow(ea / tak, 1.0f / 7.0f);
        float expected = 5.67e-8f * MathF.Pow(tak, 4.0f) * (esky * fclear + 1f - fclear);
        AlexiRad.GetXlwdn(ea, ta, fclear, out float xlwdn);
        Assert.Equal(expected, xlwdn, 3);
        Assert.InRange(xlwdn, 300f, 1500f);
    }

    [Fact]
    public void GetNetRad_FullCanopyNoTransmission()
    {
        var s = MakeRadPixel();
        s.taudl = 1f; s.taudv = 0f; s.taudn = 0f;
        float taubtv = 0f, taubtn = 0f; // tausolar = 0, tauthermal = 1
        float ts = 300f, tc = 290f, xlwdn = 400f, sdn = 600f, albedo = 0.2f;

        AlexiRad.GetNetRad(s, ts, tc, xlwdn, sdn, albedo, taubtv, taubtn,
            out float rnet, out float rnsoil, out float rndiv, out float swup, out float xlwup);

        float rcpy = 0.99f * 5.67e-8f * MathF.Pow(tc + 273.15f, 4.0f);
        float rsoil = s.emsoil * 5.67e-8f * MathF.Pow(ts + 273.15f, 4.0f);
        // tauthermal=1, tausolar=0: LW part of rndiv vanishes but SW part survives:
        float expRndiv = (1f - 1f) * (xlwdn + rsoil - 2f * rcpy) + (1f - 0f) * (1f - albedo) * sdn;
        float expRnsoil = 1f * xlwdn + 0f * rcpy - rsoil + 0f * sdn;
        Assert.Equal(expRnsoil, rnsoil, 2);
        Assert.Equal(expRndiv, rndiv, 2);
        Assert.Equal(expRndiv + expRnsoil, rnet, 2);
        // rnsw = 0 + (1-0)*(1-albedo)*sdn; swup = sdn - rnsw
        Assert.Equal(albedo * sdn, swup, 2);
        // rnlw = rnsoillw + rndivlw = (rsky - rsoil) + 0; xlwup recomputed: sdn-swup+xlwdn-rnet
        Assert.Equal(sdn - albedo * sdn + xlwdn - (expRndiv + expRnsoil), xlwup, 2);
    }

    [Fact]
    public void GetNetRadNight_SwUpEqualsSdn()
    {
        var s = MakeRadPixel();
        s.taudl = 0.5f; s.taudv = 0.1f; s.taudn = 0.2f;
        AlexiRad.GetNetRadNight(s, 300f, 290f, 400f, 0f, 0.2f, 0.1f, 0.3f,
            out _, out _, out _, out float swup, out _);
        Assert.Equal(0f, swup); // sdn=0, rnsw=0
    }

    [Fact]
    public void FindAlbedoSoil_ConvergesImmediatelyWhenObservedMatchesComputed()
    {
        var s = MakeRadPixel();
        var d = new DailyState();
        // pre-compute albv/albn at the same zen so the first iteration converges
        float zen = AlexiRad.GetSunzen(s, d, 14.0f);
        AlexiRad.GetRadProps(s, zen, out _, out _, out _, s.clumps2);
        s.albvobs = s.albv; s.albnobs = s.albn;
        float rsoilv0 = s.rsoilv, rsoiln0 = s.rsoiln;

        var sw = new StringWriter();
        AlexiRad.FindAlbedoSoil(s, d, 14.0f, sw);

        Assert.Equal(rsoilv0, s.rsoilv);
        Assert.Equal(rsoiln0, s.rsoiln);
        string outp = sw.ToString();
        // FORMAT 100: a10 label + i5 iter -> one line per loop, converged at iter 1
        int lines = 0;
        foreach (var l in outp.Split('\n')) if (l.Length > 0) lines++;
        Assert.Equal(2, lines);
        Assert.Contains("RSOILV:", outp);
        Assert.Contains("RSOILN:", outp);
        // clear-sky partitioning is set
        Assert.Equal(0.5f, s.fvis);
        Assert.Equal(0.8f, s.dirvis);
    }

    [Fact]
    public void FindAlbedoSoil_ConvergesWithinTwentyIterationsWhenObservedDiffers()
    {
        var s = MakeRadPixel();
        var d = new DailyState();
        s.albvobs = 0.20f; s.albnobs = 0.35f;
        var sw = new StringWriter();
        AlexiRad.FindAlbedoSoil(s, d, 14.0f, sw);
        // converged (not BAD) and within errmax of the targets
        Assert.NotEqual(GridDims.Bad, s.rsoilv);
        Assert.NotEqual(GridDims.Bad, s.rsoiln);
        AlexiRad.GetRadProps(s, s.zen2, out _, out _, out _, s.clumps2);
        Assert.True(MathF.Abs(s.albvobs - s.albv) <= 0.002f);
        Assert.True(MathF.Abs(s.albnobs - s.albn) <= 0.002f);
    }

    [Fact]
    public void FindAlbedoVeg_ConvergesImmediatelyWhenObservedMatchesComputed()
    {
        var s = MakeRadPixel();
        var d = new DailyState();
        float zen = AlexiRad.GetSunzen(s, d, 14.0f);
        AlexiRad.GetRadProps(s, zen, out _, out _, out _, s.clumps2);
        s.albvobs = s.albv; s.albnobs = s.albn;
        float aleafv0 = s.aleafv, aleafn0 = s.aleafn;

        var sw = new StringWriter();
        AlexiRad.FindAlbedoVeg(s, d, 14.0f, sw);

        Assert.Equal(aleafv0, s.aleafv);
        Assert.Equal(aleafn0, s.aleafn);
        // converged at iteration 1 in both loops: exactly one line each
        int lines = 0;
        foreach (var l in sw.ToString().Split('\n')) if (l.Length > 0) lines++;
        Assert.Equal(2, lines);
        Assert.Contains("ALEAFV:", sw.ToString());
    }
}
