using Alexi.Core;

namespace Alexi.IO;

/// <summary>
/// Ports of ALEXI_0.1/src/USflux.f90 (main driver), landcover.f (pipeline
/// stage run before USflux), and getdate. Field names verbatim.
/// </summary>
public static class AlexiDriver
{
    /// <summary>getdate(MDATE,year,doy) — integer arithmetic, verbatim.</summary>
    public static void GetDate(int MDATE, out int year, out int doy)
    {
        doy = MDATE % 1000;
        year = MDATE / 1000 - doy / 1000;
    }

    /// <summary>
    /// landcover(MDATE,tile,npoints,part) — reads veg_ (unit 102), writes
    /// veg_us (unit 103) consumed later by cover_props (unit 400).
    /// Fortran defects preserved: `class` is an implicit REAL local (0 in
    /// practice when xlai was BAD); the end=200 label just exits the loop.
    /// </summary>
    public static void Landcover(int MDATE, string t, int npoints, string part,
        AlexiPaths paths, LandcoverTables tab, TextWriter stdout)
    {
        stdout.WriteLine("MDATE =  " + MDATE);
        int year, doy;
        GetDate(MDATE, out year, out doy);
        stdout.WriteLine("YEAR DOY =  " + year + " " + doy);

        // loadclasstable — HOMEDIR//'store/landcover.txt'
        Cover.LoadClassTable(tab, paths.TableDir + "store/landcover.txt", stdout);

        string od = Pad(MDATE.ToString(), 7); // I7
        var vegIn = new ListDirectedReader(paths.InputDir + "veg_" + od + "_" + t + "_" + part + ".input");
        TextWriter vegOut = new System.IO.StreamWriter(paths.InputDir + "veg_us" + od + "_" + t + "_" + part + ".input");

        float cls = 0f; // Fortran name class (C# keyword; renamed)
        float[] freq = new float[GridDims.Nclass + 1];
        for (int ja = 1; ja <= npoints; ja++)
        {
            float[] raw = new float[12];
            if (vegIn.ReadList(12, raw) < 0)
                break; // end=200
            int i = (int)raw[0];
            int j = (int)raw[1];
            float xntot = raw[2];
            for (int k = 1; k <= GridDims.Nclass; k++)
                freq[k] = raw[2 + k];
            float xlai = raw[11] / 10f;
            if (xlai < 0f)
                xlai = GridDims.Bad;
            if (xlai > 10f)
                xlai = GridDims.Bad;
            float xndvi = 100f; // unused in output; kept for parity
            _ = xndvi;

            float fc = 0f, height = 0f, xleaf = 0f, aleafv = 0f, aleafn = 0f, aleafl = 0f, adeadv = 0f, adeadn = 0f, adeadl = 0f, z0eff = 0f, dispeff = 0f, rsmin = 0f;
            if (xlai != GridDims.Bad)
            {
                fc = 1f - MathF.Exp(-0.5f * xlai);
                if (fc < 0f)
                    fc = 0f;
                if (fc > 1f)
                    fc = 1f;
                WeightedAvg(tab, freq, ref height, ref xleaf, ref aleafv, ref aleafn, ref aleafl,
                    ref adeadv, ref adeadn, ref adeadl, ref z0eff, ref dispeff, fc, ref cls, xntot, ref rsmin);
            }
            else
            {
                xndvi = GridDims.Bad;
                fc = GridDims.Bad;
                aleafv = GridDims.Bad;
                aleafn = GridDims.Bad;
                aleafl = GridDims.Bad;
                adeadv = GridDims.Bad;
                adeadn = GridDims.Bad;
                adeadl = GridDims.Bad;
                height = GridDims.Bad;
                xleaf = GridDims.Bad;
                z0eff = GridDims.Bad;
                dispeff = GridDims.Bad;
                rsmin = GridDims.Bad;
            }
            if (cls == 1f)
            {
                xndvi = GridDims.Bad;
                fc = GridDims.Bad;
                aleafv = GridDims.Bad;
                aleafn = GridDims.Bad;
                aleafl = GridDims.Bad;
                adeadv = GridDims.Bad;
                adeadn = GridDims.Bad;
                adeadl = GridDims.Bad;
                height = GridDims.Bad;
                xleaf = GridDims.Bad;
                z0eff = GridDims.Bad;
                dispeff = GridDims.Bad;
                rsmin = GridDims.Bad;
            }

            // write(103,1030) — format (2i5,8f13.5,3f13.5)
            vegOut.WriteLine(
                $"{(int)i,5}{(int)j,5}"
                + F(aleafv, 13, 5) + F(aleafn, 13, 5) + F(aleafl, 13, 5) + F(adeadv, 13, 5)
                + F(adeadn, 13, 5) + F(adeadl, 13, 5) + F(height, 13, 5) + F(xleaf, 13, 5)
                + F(z0eff, 13, 5) + F(dispeff, 13, 5) + F(rsmin, 13, 5));
        }
        vegOut.Close();
    }

    /// <summary>weighted_avg(...) — class-weighted landcover aggregation.</summary>
    private static void WeightedAvg(LandcoverTables tab, float[] freq,
        ref float height, ref float xleaf, ref float aleafv, ref float aleafn, ref float aleafl,
        ref float adeadv, ref float adeadn, ref float adeadl, ref float z0eff, ref float dispeff,
        float fc, ref float cls, float xntot, ref float rsmin)
    {
        const float z0w = 0.00035f; // roughness length for water

        height = 0f;
        xleaf = 0f;
        aleafv = 0f;
        aleafn = 0f;
        aleafl = 0f;
        adeadv = 0f;
        adeadn = 0f;
        adeadl = 0f;
        float zsum = 0f;
        float dsum = 0f;
        float xn = 0f;
        float xnv = 0f;
        rsmin = 0f;

        if (cls == 1f)
        { // pixel is predominantly water
            height = GridDims.Bad;
            xleaf = GridDims.Bad;
            aleafv = GridDims.Bad;
            aleafn = GridDims.Bad;
            aleafl = GridDims.Bad;
            adeadv = GridDims.Bad;
            adeadn = GridDims.Bad;
            adeadl = GridDims.Bad;
            z0eff = GridDims.Bad;
            dispeff = GridDims.Bad;
            rsmin = GridDims.Bad;
        }
        else
        {
            for (int ic = 1; ic <= GridDims.Nclass; ic++)
            {
                if (freq[ic] > 0f)
                {
                    float disp, z0;
                    if (ic == 1)
                    {
                        disp = 0f;
                        z0 = z0w;
                    }
                    else
                    {
                        float hc = tab.hmin[ic - 1] + fc * (tab.hmax[ic - 1] - tab.hmin[ic - 1]);
                        height += freq[ic] * hc;
                        xleaf += freq[ic] * tab.xl[ic - 1];
                        aleafv += freq[ic] * tab.alv[ic - 1];
                        aleafn += freq[ic] * tab.aln[ic - 1];
                        aleafl += freq[ic] * tab.all[ic - 1];
                        adeadv += freq[ic] * tab.adv[ic - 1];
                        adeadn += freq[ic] * tab.adn[ic - 1];
                        adeadl += freq[ic] * tab.adl[ic - 1];
                        rsmin += freq[ic] * tab.rs[ic - 1];
                        CanopyArchLai(fc, hc, out z0, out disp);
                        xnv += freq[ic];
                    }
                    float lg = MathF.Log((50f - disp) / z0);
                    zsum += freq[ic] / (lg * lg);
                    dsum += freq[ic] * disp;
                    xn += freq[ic];
                }
            }

            xleaf /= xnv;
            aleafv /= xnv;
            aleafn /= xnv;
            aleafl /= xnv;
            adeadv /= xnv;
            adeadn /= xnv;
            adeadl /= xnv;
            height /= xn;
            rsmin /= xn;
            dispeff = dsum / xn;
            zsum /= xn;
            z0eff = (50f - dispeff) / MathF.Exp(1f / MathF.Sqrt(zsum));
        }
    }

    /// <summary>canopyarch_lai(fc,hc,z0,disp) — simplified Massman canopy arch.</summary>
    private static void CanopyArchLai(float fc, float hc, out float z0, out float disp)
    {
        const float z0s = 0.005f; // roughness for bare soil
        const float fcbare = 0f;

        if (fc <= fcbare || hc == 0f)
        {
            z0 = z0s;
            disp = 0f;
        }
        else
        {
            float dispdh = 2f / 3f;
            float z0dh = 1f / 8f;
            disp = dispdh * hc;
            z0 = z0dh * hc;
            if (z0 < z0s)
                z0 = z0s;
        }
    }

    /// <summary>
    /// USflux(MDATE,t,npoints,part) — main driver. Counters are saved locals
    /// in the Fortran source (not COMMON); returned as RunCounters.
    /// store_output_bin is commented out in the shipped driver; enabled via
    /// storeBinary. iflag(ia,ja) is a real (ia,ja)-indexed array; ia/ja come
    /// from the profile READ (defect: REAL indices), guarded here.
    /// </summary>
    public static RunCounters USflux(int MDATE, string t, int npoints, string part,
        AlexiPaths paths, PixelState s, DailyState d, LandcoverTables tab, PixelStreams io, bool storeBinary)
    {
        RunCounters c = new RunCounters();

        io.stdout!.WriteLine("USflux MDATE =  " + MDATE);
        int doy = MDATE % 1000;
        int year = MDATE / 1000 - doy / 1000; // inline in USflux.f90 (same as getdate)

        // logfile = './ALEXI_LOG/'//'/'//MDATE//'_'//t//'_'//part//'.LOG'
        io.log99 = new System.IO.StreamWriter(paths.LogDir + "/" + Pad(MDATE.ToString(), 7) + "_" + t + "_" + part + ".LOG");

        AlexiUtl.SetConstants(s);

        string od = Pad(MDATE.ToString(), 7);
        io.met = new ListDirectedReader(paths.InputDir + "met_" + od + "_" + t + "_" + part + ".input");
        io.sat = new ListDirectedReader(paths.InputDir + "sat_" + od + "_" + t + "_" + part + ".input");
        io.nparm = new ListDirectedReader(paths.InputDir + "nparm_" + od + "_" + t + "_" + part + ".input");
        io.profile = new ListDirectedReader(paths.InputDir + "profile_" + od + "_" + t + "_" + part + ".input");
        io.veg = new ListDirectedReader(paths.InputDir + "veg_us" + od + "_" + t + "_" + part + ".input");
        io.output103 = new System.IO.StreamWriter(paths.OutputDir + "output_" + od + "_" + t + "_" + part + ".input");

        int[] iflag = new int[GridDims.Ilg * GridDims.Jlg]; // dimension iflag(ilg,jlg)

        for (int m = 1; m <= npoints; m++)
        {
            c.ntot++;
            s.converged = false;
            int ibad = (int)GridDims.Bad;
            int iflagLocal = 0;
            int ierr = 0;
            int iter = 0;
            AlexiInput.ExtractInput(m, io, s, d, tab, out iflagLocal);
            s.fc0 = s.fc;

            if (s.badinput)
            {
                c.nbad++;
                continue; // go to 1000
            }
            else if (!s.iswater_inland &&
                (s.iclass == GridDims.Water || s.iclass == GridDims.Ice || iflagAt(iflag, (int)s.ai, (int)s.aj) == 3))
            {
                s.iswater = true;
                c.nwater++;
                continue; // go to 1000
            }
            // elseif(writeme) then — empty branch

            if (d.clear)
            {
                c.nclear++;
                s.fc = s.fc0;
                if (s.iswater_inland)
                    AlexiWater.AlexiWaterRun(s, d, io.stdout!, ref ierr, ref ibad, out iter);
                else
                    AlexiCore.Alexi(s, d, io.stdout!, ref ierr, ref ibad, out iter);

                if (s.badinput)
                {
                    c.nbad++;
                    setIflag(iflag, (int)s.ai, (int)s.aj, 1);
                    continue; // go to 1000
                }
                else if (!s.converged)
                {
                    // go to 900 — at 900, converged is false here, so:
                    c.nfail++;
                    continue; // go to 1000 (shipped code skips cloudy_day_proc here)
                }
                c.nconv++;
                AlexiUSFluxClear.ClearDayProc(s, d, io.stdout!, ref ibad);
            }
            else
            {
                c.ncloud++;
                AlexiUSFluxCloud.CloudyDayProc(s, d, io.stdout!);
            }

            // 1001: screen_output
            AlexiInput.ScreenOutput(io, s, d);

            // 1000: store_output_bin — commented out in the shipped driver
            if (storeBinary)
                AlexiInput.StoreOutputBin((int)s.ai, (int)s.aj, ref ierr, ref ibad, iflagLocal, iter, io, s, d);
        }

        io.log99!.WriteLine("          Bad input:  " + c.nbad);
        io.log99!.WriteLine("       Water or ice:  " + c.nwater);
        io.log99!.WriteLine("       Cloudy pixel:  " + c.ncloud);
        io.log99!.WriteLine("          Converged:  " + c.nconv);
        io.log99!.WriteLine("   Did not converge:  " + c.nfail);
        io.log99!.WriteLine("       Total points:  " + c.ntot);
        io.log99!.Close();
        io.output103!.Close();

        return c;
    }

    /// <summary>iflag(ia,ja) read — REAL indices from common; guarded (Fortran UB).</summary>
    private static int iflagAt(int[] iflag, int i, int j)
    {
        return inGrid(i, j) ? iflag[GridDims.Index2(i, j)] : 0;
    }

    private static void setIflag(int[] iflag, int i, int j, int v)
    {
        if (inGrid(i, j))
            iflag[GridDims.Index2(i, j)] = v;
    }

    private static bool inGrid(int i, int j)
    {
        return i >= 1 && i <= GridDims.Ilg && j >= 1 && j <= GridDims.Jlg;
    }

    /// <summary>Fortran f{width}.{decimals} fixed format, locale-invariant.</summary>
    private static string F(float v, int width, int decimals)
    {
        string t = v.ToString("F" + decimals, System.Globalization.CultureInfo.InvariantCulture);
        return t.Length >= width ? t : Spaces(width - t.Length) + t;
    }

    private static string Pad(string t, int width) =>
        t.Length >= width ? t : Spaces(width - t.Length) + t;

    private static string Spaces(int n) =>
        n > 0 ? "                                                                                                    ".Substring(0, n) : "";
}

public sealed class RunCounters
{
    public int ntot, nbad, nconv, nfail, ncloud, nclear, nwater;
}
