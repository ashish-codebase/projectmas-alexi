using Alexi.Core;
using Xunit;

namespace Alexi.Core.Tests;

public class AlexiUSFluxUtlTests
{
    [Fact]
    public void Integrate_SumsAndScalesToMJPerDay()
    {
        var f = new float[GridDims.Nohr];
        var t = new float[GridDims.Nohr];
        for (int i = 0; i < GridDims.Nohr; i++) { f[i] = 10f; t[i] = i; }
        AlexiUSFluxUtl.Integrate(f, out float fint, 0f, 24f, t, 24);
        // sum=240, xn=24/(1-0)=24 -> (240/24)*3600*24/1e6 = 0.864
        Assert.Equal(0.864f, fint, 4);
    }

    [Fact]
    public void Integrate_AnyBadHourPoisonsResult()
    {
        var f = new float[GridDims.Nohr];
        var t = new float[GridDims.Nohr];
        for (int i = 0; i < GridDims.Nohr; i++) { f[i] = 10f; t[i] = i; }
        f[5] = GridDims.Bad;
        AlexiUSFluxUtl.Integrate(f, out float fint, 0f, 24f, t, 24);
        Assert.Equal(GridDims.Bad, fint);
    }

    [Fact]
    public void GetTdepart2_QuadraticPassesThroughPoint()
    {
        // fit: trad at t2; curve a(t)=xm t^2 + xb t + xc should give
        // a(t2) = trad-ta and a at sunrise/sunset (tr,ts) consistent
        float trad = 30f, ta = 20f, t2 = 14f, tr = 6f, ts = 18f;
        float xm = 0f, xb = 0f, xc = 0f;
        AlexiUSFluxUtl.GetTdepart2(ref xm, ref xb, ref xc, trad, ta, t2, tr, ts);
        float atT2 = xm * t2 * t2 + xb * t2 + xc;
        Assert.Equal(trad - ta, atT2, 2);
    }

    [Fact]
    public void SetPBLHeights_200mSpacing()
    {
        var s = new PixelState();
        AlexiUSFluxUtl.SetPBLHeights(s);
        Assert.Equal(0f, s.zpbli[0]);
        Assert.Equal(200f, s.zpbli[1]);
        Assert.Equal(200f * (GridDims.Mli - 1), s.zpbli[GridDims.Mli - 1]);
    }

    [Fact]
    public void DailyFluxCloudy_SetsPlaceholdersBad()
    {
        var d = new DailyState();
        AlexiUSFluxCloud.DailyFluxCloudy(d);
        Assert.Equal(GridDims.Bad, d.eday);
        Assert.Equal(GridDims.Bad, d.hday);
        Assert.Equal(GridDims.Bad, d.ecday);
        Assert.Equal(GridDims.Bad, d.esday);
        Assert.Equal(GridDims.Bad, d.hcday);
        Assert.Equal(GridDims.Bad, d.hsday);
    }
}
