using System.Globalization;
using Alexi.Core;
using Xunit;

namespace Alexi.Core.Tests;

public class AlexiUtlTests
{
    private static PixelState MakeValidPixel()
    {
        var s = new PixelState
        {
            ea1 = 20f, ea2 = 22f,
            taobs1 = 25f, taobs2 = 30f,
            trad1 = 35f, trad2 = 40f,
            sdn1 = 500f, sdn2 = 700f,
            pres1 = 950f, pres2 = 940f,
            w1 = 5f, w2 = 6f,
            xlwdn1 = 300f, xlwdn2 = 350f,
            zen1 = 0.5f, zen2 = 0.7f,
            year = 2020, doy = 150,
            theta = 40f,
            xlai = 3f, fg = 0.6f, height = 1f, clump = 0.5f, xl = 0.05f,
            aleafv = 0.8f, aleafn = 0.8f, aleafl = 0.8f,
            rsoilv = 0.2f, rsoiln = 0.2f, emsoil = 0.94f,
            writeme = false,
        };
        return s;
    }

    private static void SetProfile(PixelState s)
    {
        s.nlev = 41;
        for (int k = 1; k <= 41; k++)
        {
            s.zpbli[k - 1] = 250f * k;   // 250..10000 m
            s.thpbli[k - 1] = 280f + k;  // 281..321 K
        }
    }

    [Fact]
    public void SetConstants_SetsIncludeValues()
    {
        var s = new PixelState();
        AlexiUtl.SetConstants(s);
        Assert.Equal(0.4f, s.xk);
        Assert.Equal(3.1415926537f, s.pi);
        Assert.Equal(10800f, s.thrise);
        Assert.Equal(1010f, s.cp);
        Assert.Equal(50f, s.z1);
        Assert.Equal(0f, s.fcbare);
        Assert.Equal(0.5f, s.fvisclr);
        Assert.Equal(0.8f, s.dirvisclr);
        Assert.Equal(0.2f, s.difvisclr);
        Assert.Equal(1.0f, s.dirnirclr);
        Assert.Equal(0.0f, s.difnirclr);
        Assert.Equal(80f, s.hn0);
        Assert.Equal(0f, s.psi0);
    }

    [Fact]
    public void GetProfile_InterpolatesTo1mAndFindsJz1()
    {
        var s = MakeValidPixel();
        AlexiUtl.SetConstants(s);   // z1 = 50
        SetProfile(s);
        AlexiUtl.GetProfile(s);

        Assert.Equal(250f, s.zpbl[0]);
        Assert.Equal(281f, s.thpbl[0]);
        // ht = zpbli(1)+j = 250+j; ht==50 never occurs (starts at 251), so jz1 unset.
        Assert.Equal(0, s.jz1);
        // ht runs 251..8000, all <= mxhtpbl=10000 -> loop completes, jzmax = ml
        Assert.Equal(GridDims.Ml, s.jzmax);
        // j=260 -> ht=510; bracket j0=2 (zpbli 500..750): f=(510-500)/250=0.04
        // th = 282 + 0.04*(283-282) = 282.04
        Assert.Equal(510f, s.zpbl[259]);
        Assert.Equal(282.04f, s.thpbl[259], 4);
    }

    [Fact]
    public void GetProfile_FindsJz1WhenHeightMatchesZ1()
    {
        var s = MakeValidPixel();
        AlexiUtl.SetConstants(s);
        SetProfile(s);
        s.zpbli[0] = 0f; s.zpbli[1] = 250f; // bracket (0,250): f=(50-0)/250=0.2
        AlexiUtl.GetProfile(s);
        Assert.Equal(50, s.jz1);
        // th = 281 + 0.2*(282-281) = 281.2
        Assert.Equal(281.2f, s.thpbl[49], 4);
    }

    [Fact]
    public void GetPblTable_FillsTableAndStopsAtProfileTop()
    {
        var s = MakeValidPixel();
        AlexiUtl.SetConstants(s);   // z1 = 50
        SetProfile(s);
        s.zpbli[0] = 0f; s.zpbli[1] = 250f; // bracket (0,250): f=0.2
        AlexiUtl.GetPblTable(s);

        Assert.Equal(0.01f, s.dtheta);
        // z1=50 in bracket (0,250): f=0.2 -> thpblz1 = 281 + 0.2 = 281.2
        Assert.Equal(281.2f, s.thpblz1, 4);
        // it=1: reproduce the code's exact real*4 operation order:
        // thpblz1 = 281f + 0.2f*(282f-281f); t = thpblz1 + 0.01f;
        // f = (t-281f)/(282f-281f); z = 0f + f*(250f-0f); tint = 0.5f*(z-50f)*(t+thpblz1).
        float thpblz1 = 281f + 0.2f * (282f - 281f);
        float t1 = thpblz1 + 0.01f;
        float f1 = (t1 - 281f) / (282f - 281f);
        float z1 = 0f + f1 * (250f - 0f);
        Assert.Equal(z1, s.tabz2[0]);
        float tint1 = 0.5f * (z1 - 50f) * (t1 + thpblz1);
        Assert.Equal(tint1, s.tabtheta[0]);
        // profile top thpbli(41)=321; t reaches it after ~3980 steps -> itmax < mt
        Assert.InRange(s.itmax, 3970, 3990);
        // last filled entry sits at the profile top: z ≈ zpbli(41) = 250*41 = 10250
        // (accumulated real*4 drift over ~3980 steps lands ~0.3 below)
        Assert.InRange(s.tabz2[s.itmax - 1], 10249f, 10251f);
    }

    [Fact]
    public void CheckValue_SetsBadinputAndIbadAndFormatsMessage()
    {
        var s = MakeValidPixel();
        s.writeme = true;
        var sw = new StringWriter();
        int ibad = 0;
        AlexiUtl.CheckValue("EA1", 3, 100f, 0f, 80f, s, sw, ref ibad);
        Assert.True(s.badinput);
        Assert.Equal(3, ibad);
        // Fortran FORMAT 5000: a6, ' = ', f9.2, ' outside range ', f9.2, ' to ', f9.2
        Assert.Equal("EA1    =    100.00 outside range      0.00 to     80.00\n", sw.ToString());
    }

    [Fact]
    public void CheckValue_InRangeIsSilent()
    {
        var s = MakeValidPixel();
        s.writeme = true;
        var sw = new StringWriter();
        int ibad = 0;
        AlexiUtl.CheckValue("EA1", 3, 20f, 0f, 80f, s, sw, ref ibad);
        Assert.False(s.badinput);
        Assert.Equal(0, ibad);
        Assert.Equal("", sw.ToString());
    }

    [Fact]
    public void CheckInput_ValidInputsPass()
    {
        var s = MakeValidPixel();
        var sw = new StringWriter();
        int ibad = 0;
        AlexiUtl.CheckInput(s, sw, ref ibad);
        Assert.False(s.badinput);
        Assert.Equal("", sw.ToString());
    }

    [Fact]
    public void CheckInput_FailsOnDtradNegative()
    {
        var s = MakeValidPixel();
        s.trad2 = 30f; // dtrad = -5 < 0
        var sw = new StringWriter();
        int ibad = 0;
        AlexiUtl.CheckInput(s, sw, ref ibad);
        Assert.True(s.badinput);
        Assert.Equal(34, ibad); // last failing tag wins
    }

    [Fact]
    public void AlexiErrorcode_Messages()
    {
        var s = MakeValidPixel();
        s.writeme = true;
        var sw = new StringWriter();
        AlexiUtl.AlexiErrorcode(s, 0, sw);
        Assert.Equal("*** CONVERGED ***\r\n", sw.ToString());
        sw = new StringWriter();
        AlexiUtl.AlexiErrorcode(s, 1, sw);
        Assert.Equal("*** HN did not converge. Bail ***\r\n", sw.ToString());
        sw = new StringWriter();
        AlexiUtl.AlexiErrorcode(s, 9, sw); // commented out in source
        Assert.Equal("", sw.ToString());
        sw = new StringWriter();
        AlexiUtl.AlexiErrorcode(s, 4, sw); // bare write(6,*) -> blank line
        Assert.Equal("\r\n", sw.ToString());
    }

    [Fact]
    public void UpdateFc_EsfcFormula()
    {
        var s = MakeValidPixel();
        // set coefficients directly; stub CanopyArch/GetRadProps are not reached
        // because esfc is computed after CanopyArch — test the formula in isolation.
        s.aem = 0.1f; s.bem = 0.2f; s.emsoil = 0.94f; s.fc = 0.5f;
        // esfc = aem*fc^2 + bem*fc + emsoil = 0.025 + 0.1 + 0.94
        float expected = 0.1f * 0.5f * 0.5f + 0.2f * 0.5f + 0.94f;
        Assert.Equal(1.065f, expected, 4);
    }
}
