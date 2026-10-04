using Alexi.Core;
using Xunit;

namespace Alexi.Core.Tests;

// Tests for functions ported in the second pass (Cover, TwoSource, Water).
public class NewCoreTests
{
    [Fact]
    public void GetFveg_SparseCanopy_ReturnsUpperClamp()
    {
        // gap = exp(-0.5*1*2) = 0.3679; rhs(fv) = fv*exp(-2/fv)+1-fv is
        // minimized at fv=1 (rhs=exp(-1)=0.3679), so the grid search lands on 1.
        float fveg;
        Cover.GetFveg(1f, 2f, out fveg);
        Assert.Equal(1f, fveg);
    }

    [Fact]
    public void GetFveg_DenseCanopy_FindsInteriorMinimum()
    {
        // gap = exp(-0.5*0.5*2) = 0.6065; the 100-step search minimizes at
        // fv≈0.44 (rhs=0.6046).
        float fveg;
        Cover.GetFveg(0.5f, 2f, out fveg);
        Assert.InRange(fveg, 0.42f, 0.46f);
    }

    [Fact]
    public void CheckSoln_FluxRatioBelowFloor_FailsWith12()
    {
        bool converged = true;
        int ierr = 0;
        AlexiCore.CheckSoln(100f, 5f, 20f, 10f, ref converged, ref ierr);
        // f = g/rnsoil = 0.05 < 0.1
        Assert.False(converged);
        Assert.Equal(12, ierr);
    }

    [Fact]
    public void CheckSoln_SurfaceAirSpreadTooLarge_FailsWith13()
    {
        bool converged = true;
        int ierr = 0;
        AlexiCore.CheckSoln(100f, 50f, 30f, 10f, ref converged, ref ierr);
        // f = 0.5 ok; tdiff = 20 > 15
        Assert.False(converged);
        Assert.Equal(13, ierr);
    }

    [Fact]
    public void CheckSoln_AcceptsConsistentSolution()
    {
        bool converged = true;
        int ierr = 0;
        AlexiCore.CheckSoln(100f, 50f, 20f, 10f, ref converged, ref ierr);
        Assert.True(converged);
        Assert.Equal(0, ierr);
    }

    [Fact]
    public void GetNetRadWater_OpenWater_UsesSigmaFourthPower()
    {
        // trad=300K: rwater = 0.99*5.67e-8*(573.15)^4 ≈ 6057
        // rnet = (100-6057) + 200*0.9 ≈ -5777 ; swup = 20 ; lwup ≈ 6057
        float rnet, swup, lwup;
        AlexiWater.GetNetRadWater(300f, 100f, 200f, out rnet, out swup, out lwup);
        Assert.InRange(lwup, 6000f, 6120f);
        Assert.Equal(20f, swup);
        Assert.InRange(rnet, -5850f, -5700f);
    }
}
