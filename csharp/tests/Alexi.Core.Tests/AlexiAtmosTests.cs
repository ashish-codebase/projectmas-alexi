using Alexi.Core;
using Xunit;

namespace Alexi.Core.Tests;

public class AlexiAtmosTests
{
    static PixelState MakeCanopy()
    {
        var s = new PixelState();
        AlexiUtl.SetConstants(s); // pi, cp, xk
        s.refhtw = 2.0f; s.height = 1.2f; s.disp = 0.8f; s.z0 = 0.1f;
        s.xl = 0.1f; s.xlai = 3.0f; s.fc = 0.8f; s.fg = 1f; s.clump = 1f;
        s.fveg = 1f; s.theta = 0.6f;
        s.zdlamx = 0.05f; s.zdlamn = -0.5f;
        return s;
    }

    [Fact]
    public void PsimhnKy_NeutralStabilityIsZero()
    {
        AlexiAtmos.PsimhnKy(0f, out float psima, out float psih);
        Assert.Equal(0f, psima);
        Assert.Equal(0f, psih);
    }

    [Fact]
    public void PsimhnKy_UnstableGivesPositiveCorrections()
    {
        // zdla=-1 -> y=1: psih=1.2*ln(4.0303)=1.6726; psima=2.1258-1.0517=1.0741
        AlexiAtmos.PsimhnKy(-1f, out float psima, out float psih);
        Assert.Equal(1.0741f, psima, 3);
        Assert.Equal(1.6726f, psih, 3);
    }

    [Fact]
    public void Psimhn_StableBranchIsLinear()
    {
        AlexiAtmos.Psimhn(0.02f, out float psima, out float psih);
        Assert.Equal(-0.1f, psima, 5);
        Assert.Equal(-0.1f, psih, 5);
    }

    [Fact]
    public void CanopyArch_WaterSetsWaterRoughnessAndSkipsExtinction()
    {
        var s = MakeCanopy();
        s.iswater = true;
        s.a = 0f; s.uexp1 = 0f; s.expuxp = 1f; s.uexp2 = 0f;
        AlexiAtmos.CanopyArch(s);
        Assert.Equal(0.00035f, s.z0);
        Assert.Equal(0f, s.disp);
        Assert.Equal(MathF.Log(2.0f / 0.00035f), s.xlog1, 4);
        // wind-extinction block skipped (go to 500)
        Assert.Equal(0f, s.a);
        Assert.InRange(s.ftheta, 0f, 0.8f);
    }

    [Fact]
    public void CanopyArch_BareSoilOnlySetsXlog1()
    {
        var s = MakeCanopy();
        s.perennial = false; s.iswater = false; s.fcbare = 0.9f; // fc=0.8 <= fcbare
        s.a = 0f;
        AlexiAtmos.CanopyArch(s);
        Assert.Equal(MathF.Log(2.0f / 0.1f), s.xlog1, 4);
        Assert.Equal(0f, s.a); // extinction block skipped
    }

    [Fact]
    public void CanopyArch_DenseCanopyComputesExtinctionFactors()
    {
        var s = MakeCanopy();
        s.perennial = false; s.iswater = false; s.fcbare = 0.1f;
        AlexiAtmos.CanopyArch(s);
        // z0 clamped up to soil roughness, height keeps clearance
        Assert.InRange(s.z0, 0.005f, 0.1f);
        Assert.True(s.height >= s.disp + s.z0 + 0.001f);
        Assert.True(s.a > 0f);
        Assert.Equal(MathF.Log((2.0f - 0.8f) / s.z0), s.xlog1, 4);
        Assert.Equal(MathF.Log((s.height - s.disp) / s.z0), s.xlog2, 4);
        Assert.True(s.expuxp <= 0.95f);
        Assert.True(s.uexp2 < 0f);
        Assert.True(s.ftheta <= 0.8f);
    }

    [Fact]
    public void GetResistance_BarePathReturnsZeroSoilResistance()
    {
        var s = MakeCanopy();
        s.perennial = false; s.iswater = false; s.fcbare = 0.9f; // bare
        s.xlog1 = MathF.Log(2.0f / 0.1f); s.xlog2 = MathF.Log(1.2f / 0.1f);
        AlexiAtmos.GetResistance(s, 1.0f, 5.0f, 20f, 25f, 30f, 1300f,
            out float ra, out float rs, out float rx);
        Assert.Equal(0f, rs);
        Assert.Equal(10000f, rx);
        Assert.True(ra > 0f);
    }

    [Fact]
    public void GetResistance_CanopyPathGivesPositiveResistances()
    {
        var s = MakeCanopy();
        s.perennial = true;
        s.xlog1 = MathF.Log((2.0f - 0.8f) / 0.1f);
        s.xlog2 = MathF.Log((1.2f - 0.8f) / 0.1f);
        s.expuxp = 0.5f; s.uexp2 = -1f;
        AlexiAtmos.GetResistance(s, 1.0f, 5.0f, 20f, 25f, 30f, 1300f,
            out float ra, out float rs, out float rx);
        Assert.True(ra > 0f);
        Assert.True(rs > 0f);
        Assert.True(rx > 0f);
    }
}
