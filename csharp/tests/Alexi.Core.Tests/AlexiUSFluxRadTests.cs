using Alexi.Core;
using Xunit;

namespace Alexi.Core.Tests;

public class AlexiUSFluxRadTests
{
    [Fact]
    public void GetRadComps_NightZerosSdnAndUsesClearAssumption()
    {
        var s = new PixelState();
        AlexiUtl.SetConstants(s);
        float sdn = 50f;
        float fclear = AlexiUSFluxRad.GetRadComps(s, 1.57f, ref sdn); // cos(1.57) ≈ 0.0008 < 0.01
        Assert.Equal(1f, fclear);
        Assert.Equal(0f, sdn);
        Assert.Equal(0.5f, s.fvis);
        Assert.Equal(1f, s.difvis);
        Assert.Equal(0f, s.dirvis);
    }

    [Fact]
    public void GetRadComps_DayPartitionsConsistently()
    {
        var s = new PixelState();
        AlexiUtl.SetConstants(s);
        float sdn = 600f;
        float zen = 0.7f;
        float fclear = AlexiUSFluxRad.GetRadComps(s, zen, ref sdn);

        Assert.True(fclear > 0f && fclear <= 1f);
        Assert.Equal(1f, s.fvis + s.fnir, 4);
        Assert.Equal(1f, s.dirvis + s.difvis, 4);
        Assert.Equal(1f, s.dirnir + s.difnir, 4);
        Assert.True(s.dirvis >= 0f && s.dirvis <= 1f);
        Assert.True(s.dirnir >= 0f && s.dirnir <= 1f);
        // full clear sky: dirvis close to its beam fraction fb1
        if (fclear >= 0.9f)
            Assert.Equal(0f, 0.9f - fclear, 3); // ratiox capped at .9
    }

    [Fact]
    public void GetNetRadSimple_ClearSkyHandCheck()
    {
        float ta = 20f, fclear = 1f, albedo = 0.2f, sdn = 600f;
        float tak = ta + 273.15f;
        float eskyc = 9.2e-6f * tak * tak;
        float esky = (1f - 0.84f * fclear) * eskyc + 0.84f * fclear;
        float expected = 0.98f * esky * 5.67e-8f * MathF.Pow(tak, 4f)
                         - 0.98f * 5.67e-8f * MathF.Pow(tak + 4f, 4f)
                         + sdn * (1f - albedo);
        float rnet = AlexiUSFluxRad.GetNetRadSimple(sdn, ta, fclear, 0.5f, albedo);
        Assert.Equal(expected, rnet, 2);
    }
}
