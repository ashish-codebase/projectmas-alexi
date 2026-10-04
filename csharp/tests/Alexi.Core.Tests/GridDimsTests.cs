using Alexi.Core;
using Xunit;

namespace Alexi.Core.Tests;

public class GridDimsTests
{
    [Fact]
    public void Index2_MatchesFortranColumnMajor()
    {
        // Fortran array a(ilg, jlg): element (1,1) is first, (2,1) second.
        Assert.Equal(0, GridDims.Index2(1, 1));
        Assert.Equal(1, GridDims.Index2(2, 1));
        Assert.Equal(GridDims.Ilg, GridDims.Index2(1, 2));
        Assert.Equal((GridDims.Jlg - 1) * GridDims.Ilg + (GridDims.Ilg - 1),
                     GridDims.Index2(GridDims.Ilg, GridDims.Jlg));
    }

    [Fact]
    public void Index3Kz_KIsFastestVarying()
    {
        // ztan(kz, ilg, jlg): k varies fastest in Fortran.
        Assert.Equal(0, GridDims.Index3Kz(1, 1, 1));
        Assert.Equal(1, GridDims.Index3Kz(2, 1, 1));
        Assert.Equal(GridDims.Kz, GridDims.Index3Kz(1, 2, 1));
        Assert.Equal(GridDims.Kz * GridDims.Ilg, GridDims.Index3Kz(1, 1, 2));
    }

    [Fact]
    public void DirectRecord_PreservesRowFlip()
    {
        // USflux_run.f:920: irec = (jlg - j) * ilg + i, with jj = jlg - j + 1.
        Assert.Equal(1L, GridDims.DirectRecord(1, GridDims.Jlg));          // jj = 1
        Assert.Equal((long)(GridDims.Jlg - 1) * GridDims.Ilg + 1, GridDims.DirectRecord(1, 1));
        Assert.Equal((long)GridDims.Ilg * GridDims.Jlg, GridDims.DirectRecord(GridDims.Ilg, 1)); // last record
    }

    [Fact]
    public void BadSentinel_IsExactFloatConstant()
    {
        Assert.Equal(-9999f, GridDims.Bad);
    }

    [Fact]
    public void GridDims_MatchIncludeFiles()
    {
        Assert.Equal(41, GridDims.Mli);
        Assert.Equal(8000, GridDims.Ml);
        Assert.Equal(8000, GridDims.Mt);
        Assert.Equal(30, GridDims.Kz);
        Assert.Equal(1456, GridDims.Ilg);
        Assert.Equal(625, GridDims.Jlg);
        Assert.Equal(1440, GridDims.Kx);
        Assert.Equal(600, GridDims.Ky);
        Assert.Equal(8, GridDims.Kt);
        Assert.Equal(24, GridDims.Nohr);
        Assert.Equal(8, GridDims.Nclass);
        Assert.Equal(6, GridDims.Nclasm);
        Assert.Equal(20, GridDims.Nbin);
    }

    [Fact]
    public void StateArrays_HaveExpectedSizes()
    {
        var grid = new GridState();
        Assert.Equal(GridDims.Mli * GridDims.Ilg * GridDims.Jlg, grid.htht.Length);
        Assert.Equal(GridDims.Kz * GridDims.Ilg * GridDims.Jlg, grid.ztan.Length);
        Assert.Equal(GridDims.Kx * GridDims.Ky * GridDims.Kt, grid.cta.Length);
        var pixel = new PixelState();
        Assert.Equal(GridDims.Ml, pixel.zpbl.Length);
        Assert.Equal(GridDims.Mt, pixel.tabtheta.Length);
        var daily = new DailyState();
        Assert.Equal(GridDims.Nohr, daily.tloc.Length);
    }

    [Fact]
    public void LandcoverTableIndexes_MatchFortranLayout()
    {
        // tablai(nbin, nclasm): bin varies fastest.
        Assert.Equal(0, LandcoverTables.TableIndex(1, 1));
        Assert.Equal(1, LandcoverTables.TableIndex(2, 1));
        Assert.Equal(GridDims.Nbin, LandcoverTables.TableIndex(1, 2));
        // tabaleaf(nclass, 3): class varies fastest.
        Assert.Equal(0, LandcoverTables.BandIndex(1, 1));
        Assert.Equal(1, LandcoverTables.BandIndex(2, 1));
        Assert.Equal(GridDims.Nclass, LandcoverTables.BandIndex(1, 2));
    }
}
