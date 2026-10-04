using Alexi.Core;
using System.Buffers.Binary;


namespace Alexi.IO;

/// <summary>
/// Ports of ALEXI_0.1/src/pbl_read.f, sfc_read.f (offline data-preparation
/// utilities) and writepoint from USflux_run.f. The source shells out to
/// gunzip/gzip around the .bin reads; that step is environment-specific and
/// is NOT reproduced here — the C# port reads the raw .bin files directly.
/// </summary>
public static class AlexiOffline
{
    /// <summary>
    /// pbl_read(MDATE) — read temperature/pressure profiles, compute ztan
    /// (potential temperature at sigma levels) and interpolate htht.
    /// The source picks kbelow=km-1 (off-by-one bracketing defect), preserved
    /// verbatim with a defensive bounds guard.
    /// </summary>
    public static void PblRead(int MDATE, string dir, GridState gs)
    {
        float[] aboveHgt = FillHgtDomain();

        int iy = MDATE / 1000 - MDATE % 1000 / 1000;
        string f1 = dir + iy + "/" + MDATE + "/temp_profile.bin";
        string f2 = dir + iy + "/" + MDATE + "/pressure_profile.bin";

        // open(36/37, recl=kx*ky*kz*4); read(36,rec=1) theta; read(37,rec=1) pres
        float[] theta = ReadGrid(f1, GridDims.Kx * GridDims.Ky * GridDims.Kz);
        float[] pres = ReadGrid(f2, GridDims.Kx * GridDims.Ky * GridDims.Kz);

        const float p0 = 100000f;
        const float e0 = 0.286f;
        for (int ja = 1; ja <= GridDims.Jlg; ja++)
        {
            for (int ia = 1; ia <= GridDims.Ilg; ia++)
            {
                if (gs.lookup_i[GridDims.Index2(ia, ja)] != -9999f &&
                    gs.lookup_j[GridDims.Index2(ia, ja)] != -9999f)
                {
                    int ip = (int)gs.lookup_i[GridDims.Index2(ia, ja)];
                    int jp = (int)gs.lookup_j[GridDims.Index2(ia, ja)];
                    for (int kc = 1; kc <= GridDims.Kz; kc++)
                    {
                        int src = GridDims.Index3KxKyKt(ip, jp, kc);
                        gs.ztan[GridDims.Index3Kz(kc, ia, ja)] =
                            theta[src] * Pow(p0 / pres[src], e0);
                    }
                }
            }
        }

        const float DZZ = 200f;
        for (int kc = 1; kc <= 41; kc++)
        {
            for (int ja = 1; ja <= GridDims.Jlg; ja++)
            {
                for (int ia = 1; ia <= GridDims.Ilg; ia++)
                {
                    if (gs.lookup_i[GridDims.Index2(ia, ja)] != -9999f &&
                        gs.lookup_j[GridDims.Index2(ia, ja)] != -9999f)
                    {
                        int kbelow = 0;
                        int kabove = 0;
                        float hzht;
                        if (kc == 1)
                        {
                            hzht = 0f;
                            gs.htht[GridDims.Index3Mli(kc, ia, ja)] = gs.ztan[GridDims.Index3Kz(GridDims.Kz, ia, ja)];
                            if (gs.htht[GridDims.Index3Mli(kc, ia, ja)] == 0f)
                                gs.htht[GridDims.Index3Mli(kc, ia, ja)] = GridDims.Bad;
                        }
                        else
                        {
                            hzht = DZZ * (kc - 1);
                        }

                        for (int km = GridDims.Kz - 1; km >= 1; km--)
                        {
                            float hcz = aboveHgt[km];
                            if (hcz > hzht && kbelow == 0)
                                kbelow = km - 1; // source defect: picks one level low
                        }
                        kabove = kbelow + 1; // "we go UP one level"
                        if (kabove < 1 || kbelow < 1 || kabove > GridDims.Kz || kbelow > GridDims.Kz)
                            continue; // guard: source reads index 0 (UB)
                        float zzab = aboveHgt[kabove];
                        float zzbl = aboveHgt[kbelow];
                        if (zzab == zzbl)
                        {
                            gs.htht[GridDims.Index3Mli(kc, ia, ja)] = GridDims.Bad;
                        }
                        else
                        {
                            gs.htht[GridDims.Index3Mli(kc, ia, ja)] =
                                ((hzht - zzbl) / (zzab - zzbl)) *
                                (gs.ztan[GridDims.Index3Kz(kabove, ia, ja)] - gs.ztan[GridDims.Index3Kz(kbelow, ia, ja)])
                                + gs.ztan[GridDims.Index3Kz(kbelow, ia, ja)];
                        }
                    }
                }
            }
        }
    }

    /// <summary>fill_hgt_domain(above_hgt) — sigma-level height domain (1-based).</summary>
    public static float[] FillHgtDomain()
    {
        float[] a = new float[GridDims.Kz + 1];
        float[] v = [15000f, 14000f, 13000f, 12000f, 11000f, 10000f, 9500f, 9000f, 8500f, 8000f,
            7500f, 7000f, 6500f, 6000f, 5500f, 5000f, 4500f, 4000f, 3500f, 3000f,
            2600f, 2200f, 1800f, 1400f, 1000f, 700f, 500f, 300f, 100f, 0f];
        for (int k = 1; k <= GridDims.Kz; k++)
            a[k] = v[k - 1];
        return a;
    }

    /// <summary>
    /// sfc_read(MDATE) — read NLDAS surface series into the hrdata_met cube
    /// arrays. CLASS.dat is read into a local (unused afterwards in the
    /// source), preserved as a local here.
    /// </summary>
    public static void SfcRead(int MDATE, string dir, GridState gs)
    {
        int iy = MDATE / 1000 - MDATE % 1000 / 1000;
        int n = GridDims.Kx * GridDims.Ky * GridDims.Kt;
        CopyInto(gs.cta, ReadGrid(dir + iy + "/" + MDATE + "/t2_series.bin", n));
        CopyInto(gs.cea, ReadGrid(dir + iy + "/" + MDATE + "/q2_series.bin", n));
        CopyInto(gs.cpres, ReadGrid(dir + iy + "/" + MDATE + "/psfc_series.bin", n));
        CopyInto(gs.cwind, ReadGrid(dir + iy + "/" + MDATE + "/wind_surface.bin", n));
        CopyInto(gs.cxlwdn, ReadGrid(dir + iy + "/" + MDATE + "/lwdn.bin", n));
        float[] data42 = ReadGrid("/data/data123/chain/4KM/GBIM/inputs/CLASS.dat", GridDims.Ilg * GridDims.Jlg);
        _ = data42;
    }

    /// <summary>arraynav — domain lat/long arrays around (37.30, -95.90).</summary>
    public static void ArrayNav(GridState gs)
    {
        const float CLAT = 37.30f;
        const float CLON = -95.90f;
        const float dlat = 0.04f;
        const float PI180f = 0.017453292519943295f;
        float hid = GridDims.Ilg / 2f;
        float hjd = GridDims.Jlg / 2f;
        gs.clat = CLAT;
        gs.clon = CLON;
        gs.dx = dlat;
        gs.dy = dlat;
        gs.minlat = CLAT - hjd * dlat;
        gs.maxlat = CLAT + hjd * dlat;
        gs.minlon = CLON - hid * (dlat * Cos(dlat * (PI180f)));
        gs.maxlon = CLON + hid * (dlat * Cos(dlat * (PI180f)));
        for (int ja = 1; ja <= GridDims.Jlg; ja++)
            gs.navlat[ja - 1] = gs.minlat + (ja - 1) * dlat;
        for (int ia = 1; ia <= GridDims.Ilg; ia++)
            gs.navlon[ia - 1] = gs.minlon + (ia - 1) * dlat;
    }

    /// <summary>
    /// writepoint(ia,ja) — debug dump of all state at one pixel, plus the
    /// single-point us.input file. Output layout is list-directed in the
    /// source; reproduced as label+value lines (debug utility).
    /// </summary>
    public static void WritePoint(int ia, int ja, PixelState s, DailyState d, GridState gs, TextWriter out70)
    {
        out70.WriteLine($"Input fields for ALEXI for point ({ia},{ja})");
        out70.WriteLine("-----------------------------------------------------");
        W(out70, "navlat(ja)", gs.navlat[ja - 1]);
        W(out70, "navlon(ia)", gs.navlon[ia - 1]);
        W(out70, "radini", gs.radini[GridDims.Index2(ia, ja)]);
        W(out70, "radfin", gs.radfin[GridDims.Index2(ia, ja)]);
        W(out70, "swbeg", gs.swbeg[GridDims.Index2(ia, ja)]);
        W(out70, "swend", gs.swend[GridDims.Index2(ia, ja)]);
        W(out70, "xlwbeg", gs.xlwbeg[GridDims.Index2(ia, ja)]);
        W(out70, "xlwend", gs.xlwend[GridDims.Index2(ia, ja)]);
        W(out70, "satang", gs.satang[GridDims.Index2(ia, ja)]);
        W(out70, "psr15", gs.psr15[GridDims.Index2(ia, ja)]);
        W(out70, "psr55", gs.psr55[GridDims.Index2(ia, ja)]);
        W(out70, "wsr15", gs.wsr15[GridDims.Index2(ia, ja)]);
        W(out70, "wsr55", gs.wsr55[GridDims.Index2(ia, ja)]);
        W(out70, "tsr15", gs.tsr15[GridDims.Index2(ia, ja)]);
        W(out70, "tsr55", gs.tsr55[GridDims.Index2(ia, ja)]);
        W(out70, "vsr15", gs.vsr15[GridDims.Index2(ia, ja)]);
        W(out70, "vsr55", gs.vsr55[GridDims.Index2(ia, ja)]);
        W(out70, "pran", gs.pran[GridDims.Index2(ia, ja)]);
        for (int k = 1; k <= 41; k++)
            W(out70, $"htht({k})", gs.htht[GridDims.Index3Mli(k, ia, ja)]);
        W(out70, "lscls", gs.lscls[GridDims.Index2(ia, ja)]);

        // Write input file for single point run (us.input)
        float[] thpblk = new float[GridDims.Mli];
        for (int jz = 0; jz < GridDims.Mli; jz++)
            thpblk[jz] = s.thpbli[jz] + 273.15f;
        int iclear = d.clear ? 1 : 0;
        out70.WriteLine($"{s.xlat} {s.xlong} {s.stdlng}");
        out70.WriteLine(s.refhtw.ToString());
        out70.WriteLine($"{s.xl} {s.clump} {s.fg}");
        out70.WriteLine($"{s.xndvi} {s.iclass}");
        out70.WriteLine($"{s.xlai} {s.height}");
        out70.WriteLine($"{s.rsoilv} {s.rsoiln} {s.emsoil}");
        out70.WriteLine(s.nlev.ToString());
        out70.WriteLine(Join(s.zpbli, 1, s.nlev));
        out70.WriteLine(Join(thpblk, 0, s.nlev));
        out70.WriteLine($"{s.year} {s.doy} {s.theta} {iclear} {d.nohrin}");
        out70.WriteLine($"{s.tloc1} {s.taobs1} {s.ea1} {s.w1} {s.pres1} {s.trad1} {s.sdn1} {s.xlwdn1} {s.rnobs1} {s.hobs1} {s.xleobs1} {s.gobs1}");
        out70.WriteLine($"{s.tloc2} {s.taobs2} {s.ea2} {s.w2} {s.pres2} {s.trad2} {s.sdn2} {s.xlwdn2} {s.rnobs2} {s.hobs2} {s.xleobs2} {s.gobs2}");
        for (int ihr = 1; ihr <= d.nohrin; ihr++)
            out70.WriteLine($"{d.tloc[ihr - 1]} {d.wind[ihr - 1]} {d.ta[ihr - 1]} {d.ea[ihr - 1]} {d.sdn[ihr - 1]} {d.pres[ihr - 1]} {d.rnobs[ihr - 1]} {d.hobs[ihr - 1]} {d.xleobs[ihr - 1]} {d.gobs[ihr - 1]} {d.aobs[ihr - 1]} {d.precp[ihr - 1]}");
    }

    private static void W(TextWriter w, string label, float v)
    {
        w.WriteLine($"       {label}= {v}");
    }

    private static string Join(float[] arr, int from, int count)
    {
        string[] parts = new string[count];
        for (int k = 0; k < count; k++)
            parts[k] = arr[from + k].ToString();
        return String.Join(" ", parts);
    }

    private static float Pow(float b, float e) => MathF.Pow(b, e);

    private static float Cos(float x) => MathF.Cos(x);

    private static float[] ReadGrid(string path, int count)
    {
        byte[] raw = System.IO.File.ReadAllBytes(path);
        float[] v = new float[count];
        for (int k = 0; k < count; k++)
            v[k] = System.Buffers.Binary.BinaryPrimitives.ReadSingleLittleEndian(raw.AsSpan(k * 4, 4));
        return v;
    }

    private static void CopyInto(float[] dst, float[] src)
    {
        for (int k = 0; k < Math.Min(dst.Length, src.Length); k++)
            dst[k] = src[k];
    }
}
