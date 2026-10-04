using Alexi.Core;
using System.Globalization;


namespace Alexi.IO;

/// <summary>
/// Port of ALEXI_0.1/src/USflux_run.f — input extraction, screen output,
/// binary output, and file opening. Field names verbatim.
/// </summary>
public static class AlexiInput
{
    /// <summary>
    /// extract_input(m,ibad,MDATE) — extract input at grid cell IA,JA.
    /// Fortran defects preserved/documented:
    ///  - MDATE and ibad arguments are unused in the source body.
    ///  - iclass is hardcoded to 5 (pvars(5) commented out).
    ///  - iflag is a local scalar here; the driver's iflag(ia,ja) array is
    ///    never written by this routine (returned via iflagOut for testability).
    ///  - checkUSinput call is commented out in the source.
    ///  - svars(9) is never assigned (the READ fills 1..8); the fg check is
    ///    always false and not reproduced.
    ///  - on EOF of the sat/nparm/profile READs, `test` keeps its saved value
    ///    (only the met READ sets iostat) — modelled via io.test.
    ///  - iclear is an unused implicit local — not ported.
    /// </summary>
    public static void ExtractInput(int m, PixelStreams io, PixelState s, DailyState d,
        LandcoverTables tab, out int iflagOut)
    {
        // read(102,*,iostat=test,end=200) i1,i2,i3,i4,(mvars(k),k=1,10)
        float[] raw = new float[14];
        io.test = io.met!.ReadList(14, raw);
        if (io.test < 0)
        {
            goto200(io, s);
            iflagOut = 0;
            return;
        }
        float[] mvars = OneBased(10, raw, 4);

        // read(110,*,end=200) i1,i2,(svars(k),k=1,8)
        float[] sraw = new float[10];
        if (io.sat!.ReadList(10, sraw) < 0)
        {
            goto200(io, s);
            iflagOut = 0;
            return;
        }
        float[] svars = OneBased(8, sraw, 2);

        // read(120,*,end=200) i1,i2,(pvars(k),k=1,6)
        float[] praw = new float[8];
        if (io.nparm!.ReadList(8, praw) < 0)
        {
            goto200(io, s);
            iflagOut = 0;
            return;
        }
        float[] pvars = OneBased(6, praw, 2);

        // read(130,*,end=200) i1,i2,(rvars(k),k=1,14)
        float[] rraw = new float[16];
        if (io.profile!.ReadList(16, rraw) < 0)
        {
            goto200(io, s);
            iflagOut = 0;
            return;
        }
        float[] rvars = OneBased(14, rraw, 2);

        s.ai = (int)rraw[0]; // i1,i2 reassigned by every READ in the source
        s.aj = (int)rraw[1];

        int iflag = 0;
        s.writeme = false;
        s.converged = false;
        s.badinput = false;
        // investigate is a static local never set true — debug block is dead.

        //                             P L A N T S
        s.iclass = 5; // hardcoded in source: iclass=5 !pvars(5)
        s.iswater = false;
        if (s.iclass == 0 || s.iclass == GridDims.Water || s.iclass == GridDims.Ice)
        {
            iflag = 3; // iflag=3 denotes water/ice
            s.iswater = true;
        }
        s.iswater_inland = false;

        s.xlai = svars[8] / 10f;
        if (s.xlai < 0f)
        { // LSA SAF bad value is -999; 0's typically indicate a bad solution
            s.xlai = GridDims.Bad;
            s.fc = GridDims.Bad;
            s.height = GridDims.Bad;
            s.z0 = GridDims.Bad;
            s.disp = GridDims.Bad;
        }
        else
        {
            s.fc = 1f - MathF.Exp(-0.5f * s.xlai);
            if (s.fc > 1f)
                s.fc = 1f;
            if (s.fc <= 0f)
                s.fc = 0.01f;
            Cover.CoverProps(s, tab, io.veg!);
            if (s.height == 0f)
            {
                s.height = 0.01f;
                s.disp = 2f / 3f * s.height;
            }
        }
        s.fg = 1f; // svars(9) is never assigned — the -9999 check is always false

        //                               S O I L
        s.rsoilv = pvars[1];
        s.rsoiln = pvars[2];
        s.emsoil = 0.94f;
        s.albedo2 = GridDims.Bad;

        //                              T I M E
        s.xlat = pvars[3];
        s.xlong = pvars[4];
        float r15 = svars[5];
        float r55 = svars[6];
        GetDgmt(s.xlong, ref io.dgmt, ref s.stdlng); // io.dgmt models the saved local
        s.tloc1 = r15 + io.dgmt; // local standard time
        s.tloc2 = r55 + io.dgmt;
        if (s.tloc1 > 24f)
            s.tloc1 -= 24f;
        if (s.tloc2 > 24f)
            s.tloc2 -= 24f;

        //                    S U R F A C E   W E A T H E R
        s.taobs1 = mvars[1] - 273.15f;
        s.taobs2 = mvars[2] - 273.15f;
        if (s.taobs1 > s.taobs2)
            s.badinput = true;
        s.w1 = (mvars[9] + mvars[10]) / 2f;
        s.w2 = s.w1;
        s.w1 = MathF.Abs(s.w1);
        s.w2 = MathF.Abs(s.w2);
        s.w2orig = MathF.Abs(s.w2);
        s.ea1 = mvars[3];
        s.ea2 = mvars[4];
        s.pres1 = mvars[5];
        s.pres2 = mvars[6];
        s.xlwdn1 = mvars[7];
        s.xlwdn2 = mvars[8];
        s.sdn1 = svars[3];
        s.sdn2 = svars[4];
        AlexiRad.GetXlwdn(s.ea1, s.taobs1, 1f, out s.xlwdn1);
        AlexiRad.GetXlwdn(s.ea2, s.taobs2, 1f, out s.xlwdn2);
        if (r55 - r15 <= 0f)
            iflag = 2;

        s.refhtw = 30f; // height of input wind data from WRF [m]
        s.zta = 50f; // model blending height [m]

        //            S U R F A C E   T E M P E R A T U R E
        s.theta = svars[7];
        s.trad1 = svars[1] - 273.15f;
        s.trad2 = svars[2] - 273.15f;
        float dthr = (s.trad2 - s.trad1) / (r55 - r15);
        float diff = dthr * s.rsmin - dthr; // rsmin set by cover_props (garbage if xlai was bad)
        float offset = diff * (r55 - r15);
        s.trad1 -= offset;
        if (s.trad1 > s.trad2)
            s.badinput = true;

        //                    B O U N D A R Y   L A Y E R
        float[] zprof = [0f, 100f, 300f, 500f, 700f, 1000f, 1400f,
            1800f, 2200f, 2600f, 3000f, 3500f, 4000f, 4500f];
        for (int k = 1; k <= 14; k++)
        {
            s.zpbli[k - 1] = zprof[k - 1];
            s.thpbli[k - 1] = rvars[k];
        }
        // nlev is never assigned anywhere in the shipped source (stays 0).

        //          SOIL AND LEAF RADIATIVE PROPERTIES
        s.eleaf = 0.97f; // leaf emissivity
        float cosTheta = MathF.Cos(s.theta);
        float b1 = -1.97f * cosTheta + 2.87f;
        float b2 = 1.86f * cosTheta - 2.62f;
        s.bem = b1 * s.eleaf + b2;
        s.aem = s.eleaf + 0.025f - s.bem - s.emsoil;
        s.esfc = s.aem * s.fc * s.fc + s.bem * s.fc + s.emsoil;

        //                              F L A G S
        if (iflag == 2)
        {
            d.clear = false; // cloudy
        }
        else
        {
            d.clear = true; // clear
        }
        // checkUSinput is commented out in the source.
        if (s.badinput)
            iflag = 1;

        // 200: if(test.lt.0) badinput=.TRUE.
        goto200(io, s);
        iflagOut = iflag;
    }

    /// <summary>Label 200 tail of extract_input.</summary>
    private static void goto200(PixelStreams io, PixelState s)
    {
        if (io.test < 0)
            s.badinput = true;
    }

    /// <summary>1-based view of a 0-based token buffer, as in Fortran arrays.</summary>
    private static float[] OneBased(int count, float[] src, int srcOffset)
    {
        float[] v = new float[count + 1];
        for (int k = 1; k <= count; k++)
            v[k] = src[srcOffset + k - 1];
        return v;
    }

    /// <summary>
    /// getdgmt(xnlon,dgmt,stdlng) — standard longitude and local-time offset.
    /// The source table uses negative (western-hemisphere) longitudes although
    /// the header says 0-360; ported verbatim. No match leaves both arguments
    /// unchanged (they are saved/common state).
    /// </summary>
    public static void GetDgmt(float xnlon, ref float dgmt, ref float stdlng)
    {
        if (xnlon >= -172.5f && xnlon < -157.5f)
        {
            stdlng = -165f; dgmt = -11f;
        }
        else if (xnlon >= -157.5f && xnlon < -142.5f)
        {
            stdlng = -150f; dgmt = -10f;
        }
        else if (xnlon >= -142.5f && xnlon < -127.5f)
        {
            stdlng = -135f; dgmt = -9f;
        }
        else if (xnlon >= -127.5f && xnlon < -112.5f)
        {
            stdlng = -120f; dgmt = -8f;
        }
        else if (xnlon >= -112.5f && xnlon < -97.5f)
        {
            stdlng = -105f; dgmt = -7f;
        }
        else if (xnlon >= -97.5f && xnlon < -82.5f)
        {
            stdlng = -90f; dgmt = -6f;
        }
        else if (xnlon >= -82.5f && xnlon < -67.5f)
        {
            stdlng = -75f; dgmt = -5f;
        }
        else if (xnlon >= -67.5f && xnlon < -52.5f)
        {
            stdlng = -60f; dgmt = -4f;
        }
        else if (xnlon >= -52.5f && xnlon < -37.5f)
        {
            stdlng = -45f; dgmt = -3f;
        }
        else if (xnlon >= -37.5f && xnlon < -22.5f)
            {
                stdlng = -30f; dgmt = -2f;
            }
            else if (xnlon >= -22.5f && xnlon < -7.5f)
            {
                stdlng = -15f; dgmt = -1f;
            }
            else if (xnlon >= -7.5f && xnlon < 7.5f)
            {
                stdlng = 0f; dgmt = 0f;
            }
            else if (xnlon >= 7.5f && xnlon < 22.5f)
            {
                stdlng = 15f; dgmt = 1f;
            }
            else if (xnlon >= 22.5f && xnlon < 37.5f)
            {
                stdlng = 30f; dgmt = 2f;
            }
            else if (xnlon >= 37.5f && xnlon < 52.5f)
            {
                stdlng = 45f; dgmt = 3f;
            }
            else if (xnlon >= 52.5f && xnlon < 67.5f)
            {
                stdlng = 60f; dgmt = 4f;
            }
            else if (xnlon >= 67.5f && xnlon < 82.5f)
            {
                stdlng = 75f; dgmt = 5f;
            }
            else if (xnlon >= 82.5f && xnlon < 97.5f)
            {
                stdlng = 90f; dgmt = 6f;
            }
            else if (xnlon >= 97.5f && xnlon < 112.5f)
            {
                stdlng = 105f; dgmt = 7f;
            }
            else if (xnlon >= 112.5f && xnlon < 127.5f)
            {
                stdlng = 120f; dgmt = 8f;
            }
            else if (xnlon >= 127.5f && xnlon < 142.5f)
            {
                stdlng = 135f; dgmt = 9f;
            }
            else if (xnlon >= 142.5f && xnlon < 157.5f)
            {
                stdlng = 150f; dgmt = 10f;
            }
            else if (xnlon >= 157.5f && xnlon < 172.5f)
            {
                stdlng = 165f; dgmt = 11f;
            }
    }

    /// <summary>
    /// extract_hourly_input_nldas(ia,ja,dgmt) — build hourly series from the
    /// NLDAS/WRF cube arrays. Dead in the shipped pipeline (the call in
    /// extract_input is commented out), kept for parity.
    /// </summary>
    public static void ExtractHourlyInputNLDAS(PixelState s, DailyState d, GridState gs, float dgmt)
    {
        for (int ihr = 1; ihr <= GridDims.Nohr; ihr++)
        {
            d.sdn[ihr - 1] = GridDims.Bad;
            d.ta[ihr - 1] = GridDims.Bad;
            d.xlwdn[ihr - 1] = GridDims.Bad;
            d.rnet[ihr - 1] = GridDims.Bad;
            d.g[ihr - 1] = GridDims.Bad;
            d.rnsoil[ihr - 1] = GridDims.Bad;
            d.pres[ihr - 1] = GridDims.Bad;
            d.tloc[ihr - 1] = GridDims.Bad;
            d.ea[ihr - 1] = GridDims.Bad;
            d.wind[ihr - 1] = GridDims.Bad;
            d.xlst[ihr - 1] = GridDims.Bad;
        }
        for (int ihr = 1; ihr <= GridDims.Nohr; ihr++)
        {
            d.tloc[ihr - 1] = gs.ctloc[GridDims.Index3KxKyKt(s.ai, s.aj, ihr)] + dgmt;
            d.ea[ihr - 1] = gs.cea[GridDims.Index3KxKyKt(s.ai, s.aj, ihr)] * 1000f;
            d.ta[ihr - 1] = gs.cta[GridDims.Index3KxKyKt(s.ai, s.aj, ihr)] - 273.15f;
            d.wind[ihr - 1] = gs.cwind[GridDims.Index3KxKyKt(s.ai, s.aj, ihr)];
            d.pres[ihr - 1] = gs.cpres[GridDims.Index3KxKyKt(s.ai, s.aj, ihr)] / 100f;
            d.sdn[ihr - 1] = GridDims.Bad;
            d.xlst[ihr - 1] = GridDims.Bad;
            d.xlwdn[ihr - 1] = gs.cxlwdn[GridDims.Index3KxKyKt(s.ai, s.aj, ihr)];
        }
        d.nohrin = GridDims.Nohr;
    }

    /// <summary>
    /// screen_output(ia,ja) — write pixel diagnostics to screen (unit 6) and
    /// the per-pixel summary to unit 103. The verbose block only runs when
    /// writeme AND clear AND converged; writeme is never set true in the
    /// shipped extract_input, so the block is dead in practice.
    /// </summary>
    public static void ScreenOutput(PixelStreams io, PixelState s, DailyState d)
    {
        if (s.writeme && d.clear && s.converged)
        {
            io.stdout!.WriteLine($"Lat: {F(s.xlat, 6, 2)}  Lon: {F(s.xlong, 7, 2)}");
            io.stdout!.WriteLine(s.iclass.ToString());
            io.log99!.WriteLine(" LANDCOVER CLASS           ");
            W(io.stdout, "FCORIG (frac cover in)    ", s.fc0);
            W(io.stdout, "FC     (frac cover out)   ", s.fc);
            W(io.stdout, "XLAI   (LAI out)          ", s.xlai);
            W(io.stdout, "HEIGHT (canopy height)    ", s.height);
            W(io.log99, " XL     (leaf size)", s.xl);
            W(io.stdout, "Z0     (roughness)        ", s.z0);
            W(io.stdout, "TLOC1  (GOES obs1)        ", s.tloc1);
            W(io.stdout, "TLOC2  (GOES obs2)        ", s.tloc2);
            W(io.stdout, "TGMT1  (GOES obs1)        ", 0f); // r15: implicit local, never set here
            W(io.stdout, "TGMT2  (GOES obs2)        ", 0f); // r55: implicit local, never set here
            W(io.stdout, "W1     (wind speed1)      ", s.w1);
            W(io.stdout, "W2     (wind speed2)      ", s.w2);
            W(io.log99, " TC1    (canopy temp1)", s.tc1);
            W(io.log99, " TS1    (soil temp1)", s.ts1);
            W(io.log99, " TC2    (canopy temp2)", s.tc2);
            W(io.log99, " TS2    (soil temp2)", s.ts2);
            W(io.stdout, "XLEC2  (canopy latent)    ", s.xlec2);
            W(io.stdout, "XLES2  (soil latent)      ", s.xles2);
            W(io.stdout, "XLE2   (total latent)     ", s.xle2);
            W(io.log99, " HC2    (canopy sensible)", s.hc2);
            W(io.log99, " HS2    (soil sensible)", s.hs2);
            W(io.stdout, "H2     (total sensible)   ", s.h2);
            W(io.stdout, "G2     (ground cond)      ", s.g2);
            W(io.log99, " ZEN1   (sun zenith T1)", s.zen1);
            W(io.log99, " ZEN2   (sun zenith T2)", s.zen2);
            W(io.log99, " RNET1  (net rad T1)", s.rnet1);
            W(io.stdout, "RNET2  (net rad T2)       ", s.rnet2);
            W(io.log99, " TA1    (modelled air temp)", s.ta1);
            W(io.log99, " TA2    (modelled air temp)", s.ta2);
            W(io.log99, " TAOBS1 (observed air temp)", s.taobs1);
            W(io.log99, " TAOBS2 (observed air temp)", s.taobs2);
            W(io.stdout, "TRAD1  (radiometric T1)   ", s.trad1);
            W(io.stdout, "TRAD2  (radiometric T2)   ", s.trad2);
            W(io.stdout, "Z2     (boundary layer ht)", s.z2);
            W(io.log99, " TH1    (potential temp T1)", s.th1 - 273.15f);
            W(io.log99, " TH2    (potential temp T2)", s.th2 - 273.15f);
            W(io.stdout, "RA2    (aerodynamic res)  ", s.ra2);
            W(io.stdout, "RS2    (soil res)         ", s.rs2);
            W(io.stdout, "RX2    (b.l. res)         ", s.rx2);
            W(io.stdout, "SWUP                      ", s.swup2);
            W(io.stdout, "SDN2                      ", s.sdn2);
            W(io.stdout, "LWUP                      ", s.xlwup2);
            W(io.stdout, "LWDN                      ", s.xlwdn2);
            W(io.log99, " RNDAY", d.rnday);
            W(io.log99, " HDAY", d.hday);
            W(io.log99, " GDAY", d.gday);
            W(io.log99, " EDAY", d.eday);
            W(io.log99, " ESDAY", d.esday);
            W(io.log99, " ECDAY", d.ecday);
            W(io.log99, " EREFDAY", d.erefday);
            W(io.log99, " SDAY", d.sday);
        }

        // write(103,1004) int(ai),int(aj),rnet2,xle2,h2,g2,xles2 — format (2I6,5f10.1)
        io.output103!.WriteLine(
            $"{(int)s.ai,6}{(int)s.aj,6}{F(s.rnet2, 10, 1)}{F(s.xle2, 10, 1)}{F(s.h2, 10, 1)}{F(s.g2, 10, 1)}{F(s.xles2, 10, 1)}");
    }

    /// <summary>format(a27,f12.3) line.</summary>
    private static void W(TextWriter w, string label, float v)
    {
        w.WriteLine(Pad(label, 27) + F(v, 12, 3));
    }

    /// <summary>Fortran f{width}.{decimals} fixed format, locale-invariant.</summary>
    private static string F(float v, int width, int decimals)
    {
        string t = v.ToString("F" + decimals, CultureInfo.InvariantCulture);
        return Pad(t, width);
    }

    private static string Pad(string t, int width) =>
        t.Length >= width ? t : Spaces(width - t.Length) + t;

    private static string PadEnd(string t, int width) =>
        t.Length >= width ? t : t + Spaces(width - t.Length);

    private static string Spaces(int n) =>
        n > 0 ? "                                                                                                    ".Substring(0, n) : "";

    /// <summary>
    /// store_output_bin(i,j,ierr,ibad,iflag,iter) — store gridded output in
    /// binary files. Commented out in the shipped driver; enabled via the
    /// --store-binary CLI flag. xiter/tloc1o/tloc2o/t1gmt/t2gmt/class/cflag/
    /// cierr/cibad are unused implicit locals — not reproduced.
    /// </summary>
    public static void StoreOutputBin(int i, int j, ref int ierr, ref int ibad,
        int iflagIn, int iter, PixelStreams io, PixelState s, DailyState d)
    {
        if (!s.converged)
        { // didn't converge — set all instantaneous fluxes & daily fluxes to bad
            d.fsdn = GridDims.Bad;
            s.h2 = GridDims.Bad;
            s.xle2 = GridDims.Bad;
            s.g2 = GridDims.Bad;
            s.rnet2 = GridDims.Bad;
            s.ts2 = GridDims.Bad;
            s.tc2 = GridDims.Bad;
            s.t2 = GridDims.Bad;
            s.ta1 = GridDims.Bad;
            s.ta2 = GridDims.Bad;
            d.fpet = GridDims.Bad;
            s.z2 = GridDims.Bad;
            s.xlec2 = GridDims.Bad;
            s.xles2 = GridDims.Bad;
            s.hc2 = GridDims.Bad;
            s.hs2 = GridDims.Bad;
            s.rs2 = GridDims.Bad;
            s.rx2 = GridDims.Bad;
            s.ra2 = GridDims.Bad;
            s.xlwup2 = GridDims.Bad;
            s.swup2 = GridDims.Bad;
            d.eday = GridDims.Bad;
            d.ecday = GridDims.Bad;
            d.esday = GridDims.Bad;
            d.hday = GridDims.Bad;
            d.hcday = GridDims.Bad;
            d.hsday = GridDims.Bad;
        }

        float dtrad;
        if (d.clear && !s.badinput && !s.iswater)
            dtrad = s.trad2 - s.trad1;
        else
            dtrad = GridDims.Bad;

        if (s.badinput || s.iswater)
        {
            d.sday = GridDims.Bad;
            s.albedo2 = GridDims.Bad;
            d.erefday = GridDims.Bad;
            d.xlwupday = GridDims.Bad;
            d.xlwdnday = GridDims.Bad;
            d.swupday = GridDims.Bad;
            d.rnday = GridDims.Bad;
            d.gday = GridDims.Bad;
            d.hday = GridDims.Bad;
            d.eday = GridDims.Bad;
            d.ecday = GridDims.Bad;
            d.esday = GridDims.Bad;
            s.xndvi = GridDims.Bad;
        }

        if (s.iswater)
        {
            s.sdn1 = GridDims.Bad;
            s.sdn2 = GridDims.Bad;
            ierr = (int)GridDims.Bad; // ierr=BAD (integer receives -9999)
            s.trad2 = GridDims.Bad;
            s.xlai = GridDims.Bad;
            s.height = GridDims.Bad;
            s.z0 = GridDims.Bad;
            s.fc = GridDims.Bad;
            s.theta = GridDims.Bad;
            s.taobs1 = GridDims.Bad;
            s.taobs2 = GridDims.Bad;
            s.xlat = GridDims.Bad;
            s.xlong = GridDims.Bad;
            ibad = (int)GridDims.Bad;
            // iflag=3 — local copy in the source; iflagIn is the caller's value
            s.esfc = GridDims.Bad;
        }

        if (s.badinput || !d.clear || s.iswater)
        {
            s.w1 = GridDims.Bad; // "These didn't get scaled ..."
            s.w2 = GridDims.Bad;
        }

        // Compute albedo2 assuming clear-sky partitioning
        if (s.albvobs != GridDims.Bad && s.albnobs != GridDims.Bad &&
            s.rsoilv != GridDims.Bad && s.rsoiln != GridDims.Bad &&
            s.aleafv != GridDims.Bad && s.aleafn != GridDims.Bad)
        {
            s.fvis = 0.5f;
            s.fnir = 0.5f;
            s.dirvis = 0.8f;
            s.difvis = 0.2f;
            s.dirnir = 1.0f;
            s.difnir = 0.0f;
            s.zen2 = AlexiRad.GetSunzen(s, d, s.tloc2);
            float taubtv, taubtn; // implicit locals in the source (not common taubtv2)
            AlexiRad.GetRadProps(s, s.zen2, out s.albedo2, out taubtv, out taubtn, s.clumps2);
        }
        else
        {
            s.albedo2 = GridDims.Bad;
            s.albv = GridDims.Bad; // albv/albn are common/reflection2 fields
            s.albn = GridDims.Bad;
        }

        // Write data — live binwrite/binopen units only (commented calls omitted)
        BinaryGridWriter b;
        if (io.binary.ContainsKey(136)) { b = io.binary[136]; b.Write(i, j, s.fc); }
        if (io.binary.ContainsKey(137)) { b = io.binary[137]; b.Write(i, j, s.xlai); }
        if (io.binary.ContainsKey(154)) { b = io.binary[154]; b.Write(i, j, s.xlwdn2); }
        if (io.binary.ContainsKey(156)) { b = io.binary[156]; b.Write(i, j, s.trad1); }
        if (io.binary.ContainsKey(157)) { b = io.binary[157]; b.Write(i, j, s.trad2); }
        if (io.binary.ContainsKey(158)) { b = io.binary[158]; b.Write(i, j, dtrad); }
        if (io.binary.ContainsKey(160)) { b = io.binary[160]; b.Write(i, j, s.h2); }
        if (io.binary.ContainsKey(161)) { b = io.binary[161]; b.Write(i, j, s.xle2); }
        if (io.binary.ContainsKey(162)) { b = io.binary[162]; b.Write(i, j, s.g2); }
        if (io.binary.ContainsKey(163)) { b = io.binary[163]; b.Write(i, j, s.rnet2); }
        if (io.binary.ContainsKey(168)) { b = io.binary[168]; b.Write(i, j, s.z2); }
        if (io.binary.ContainsKey(169)) { b = io.binary[169]; b.Write(i, j, s.xlec2); }
        if (io.binary.ContainsKey(170)) { b = io.binary[170]; b.Write(i, j, s.xles2); }
        if (io.binary.ContainsKey(173)) { b = io.binary[173]; b.Write(i, j, s.xlwup2); }
        if (io.binary.ContainsKey(174)) { b = io.binary[174]; b.Write(i, j, s.swup2); }
        if (io.binary.ContainsKey(175)) { b = io.binary[175]; b.Write(i, j, s.ra2); }
        if (io.binary.ContainsKey(176)) { b = io.binary[176]; b.Write(i, j, s.rs2); }
        if (io.binary.ContainsKey(177)) { b = io.binary[177]; b.Write(i, j, s.rx2); }
        if (io.binary.ContainsKey(191)) { b = io.binary[191]; b.Write(i, j, s.rsoilv); }
        if (io.binary.ContainsKey(192)) { b = io.binary[192]; b.Write(i, j, s.rsoiln); }
    }

    /// <summary>
    /// open_output(year,doy) — open the direct-access binary grids. Live
    /// binopen/binopen2 calls only. binopen: dir//cvar//cyyyyddd//'.dat';
    /// binopen2 uses TABDIR with hardcoded year=2012, doy=0.
    /// </summary>
    public static void OpenOutput(AlexiPaths paths, int year, int doy, PixelStreams io)
    {
        int iyyyyddd = year * 1000 + doy;
        string tag = Pad(iyyyyddd.ToString(), 7); // write(cyyyyddd,'(i7)')

        string[] liveName = ["FCOV", "XLAI", "LWD2", "TRD1", "TRD2", "DTRD", "SENS", "LATN", "GSOL", "RNET",
            "ZPBL", "XLEC", "XLES", "LWU2", "SWU2", "RA2_", "RS2_", "RX2_", "RSLV", "RSLN"];
        int[] liveUnit = [136, 137, 154, 156, 157, 158, 160, 161, 162, 163,
            168, 169, 170, 173, 174, 175, 176, 177, 191, 192];
        for (int k = 0; k < liveName.Length; k++)
            io.binary[liveUnit[k]] = new BinaryGridWriter(paths.OutputDir + liveName[k] + tag + ".dat");

        // call binopen("MLAI",laidir,year,doy,301)
        io.binary[301] = new BinaryGridWriter(paths.LaiDir + "MLAI" + tag + ".dat");

        // binopen2 units live in TABDIR with hardcoded (2012., 0)
        string tag2 = Pad((2012 * 1000).ToString(), 7);
        string[] tabName = ["ELEV", "FMAX", "CRUF", "RSVA", "RSNA", "ALVA", "ALNA", "DIFF"];
        int[] tabUnit = [302, 303, 304, 305, 306, 307, 308, 309];
        for (int k = 0; k < tabName.Length; k++)
            io.binary[tabUnit[k]] = new BinaryGridWriter(paths.TableDir + tabName[k] + tag2 + ".dat");
    }
}
