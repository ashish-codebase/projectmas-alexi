using System.Globalization;

namespace Alexi.Core;

/// <summary>
/// Port of ALEXI_0.1/src/USflux_cover.f — landcover table loading and
/// per-pixel canopy-cover properties. Field names verbatim.
///
/// Note: cover_props reads one list-directed record (13 values) per call
/// from the veg-file stream (Fortran unit 400). The source's `do while`
/// loop over the veg file is commented out — one line per pixel.
///
/// Note: `perennial = tabperen(iclass)` is immediately overwritten with
/// `.FALSE.` in the source ("not robustly implemented") — ported verbatim,
/// so all pixels take the non-perennial branch.
/// </summary>
public static class Cover
{
    /// <summary>
    /// load_tables — read the 14-column landcover.txt into LandcoverTables
    /// (the second `common /lookup/` layout, USflux.inc:31). The first
    /// record (header line) is skipped, as in the source.
    /// </summary>
    public static void LoadTables(LandcoverTables tab, string landcoverFile, TextWriter stdout)
    {
        using var r = new StreamReader(landcoverFile);
        r.ReadLine(); // read(50,*) — skip first record (header)

        string? line;
        while ((line = r.ReadLine()) != null)
        {
            var tok = line.Split(new[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
            if (tok.Length < 14)
            {
                stdout.WriteLine($"load_tables: short record in {Path.GetFileName(landcoverFile)}: {line}");
                continue;
            }
            int iclass = ParseInt(tok[0]);
            if (iclass < 1 || iclass > GridDims.Nclass)
                continue;
            int c = iclass - 1;
            tab.itabclass[c] = iclass;
            float a = ParseF(tok[1]);
            float b = ParseF(tok[2]);
            float l = ParseF(tok[3]);
            tab.tabaleaf[LandcoverTables.BandIndex(iclass, 1)] = a;
            tab.tabaleaf[LandcoverTables.BandIndex(iclass, 2)] = b;
            tab.tabaleaf[LandcoverTables.BandIndex(iclass, 3)] = l;
            float dv = ParseF(tok[4]);
            float dn = ParseF(tok[5]);
            float dl = ParseF(tok[6]);
            tab.tabadead[LandcoverTables.BandIndex(iclass, 1)] = dv;
            tab.tabadead[LandcoverTables.BandIndex(iclass, 2)] = dn;
            tab.tabadead[LandcoverTables.BandIndex(iclass, 3)] = dl;
            tab.tabbeta[c] = ParseF(tok[7]);
            tab.tabhmin[c] = ParseF(tok[8]);
            tab.tabhmax[c] = ParseF(tok[9]);
            tab.tabxl[c] = ParseF(tok[10]);
            tab.itabclass[c] = ParseInt(tok[11]);
            tab.tabperen[c] = ParseF(tok[12]) != 0f;
            tab.tabfcmin[c] = ParseF(tok[13]);
            tab.tabdesc[c] = tok.Length > 14 ? string.Join(' ', tok.Skip(14)) : tok[13];
        }
    }

    /// <summary>
    /// loadclasstable — loads the FIRST common/lookup layout (alv..xl) from
    /// the same store/landcover.txt. The source reads 11 items per record
    /// (dum,alv,aln,all,adv,adn,adl,rs,hmin,hmax,xl); list-directed READ
    /// discards the remaining tokens on the record, so the shipped 14-column
    /// file parses correctly.
    /// </summary>
    public static void LoadClassTable(LandcoverTables tab, string landcoverFile, TextWriter stdout)
    {
        using var r = new StreamReader(landcoverFile);
        r.ReadLine(); // read(30,*) — skip header record
        for (int ic = 1; ic <= GridDims.Nclass; ic++)
        {
            string? line = r.ReadLine();
            if (line == null)
                break;
            var tok = line.Split(new[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
            if (tok.Length < 11)
                continue;
            int c = ic - 1;
            tab.alv[c] = ParseF(tok[1]);
            tab.aln[c] = ParseF(tok[2]);
            tab.all[c] = ParseF(tok[3]);
            tab.adv[c] = ParseF(tok[4]);
            tab.adn[c] = ParseF(tok[5]);
            tab.adl[c] = ParseF(tok[6]);
            tab.rs[c] = ParseF(tok[7]);
            tab.hmin[c] = ParseF(tok[8]);
            tab.hmax[c] = ParseF(tok[9]);
            tab.xl[c] = ParseF(tok[10]);
        }
        stdout.WriteLine("RS = " + string.Join(" ", tab.rs)); // write(6,*)"RS = ",rs
    }

    /// <summary>
    /// cover_props(ia,ja) — per-pixel cover properties from the lookup
    /// tables, then overwritten with one line of the veg file.
    /// </summary>
    public static void CoverProps(PixelState s, LandcoverTables tab, ListDirectedReader veg)
    {
        // Table-driven initial values (second lookup layout)
        s.aleafv = tab.tabaleaf[LandcoverTables.BandIndex(s.iclass, 1)];
        s.aleafn = tab.tabaleaf[LandcoverTables.BandIndex(s.iclass, 2)];
        s.aleafl = tab.tabaleaf[LandcoverTables.BandIndex(s.iclass, 3)];
        s.adeadv = tab.tabadead[LandcoverTables.BandIndex(s.iclass, 1)];
        s.adeadn = tab.tabadead[LandcoverTables.BandIndex(s.iclass, 2)];
        s.adeadl = tab.tabadead[LandcoverTables.BandIndex(s.iclass, 3)];
        s.beta = tab.tabbeta[s.iclass - 1];
        s.perennial = tab.tabperen[s.iclass - 1];
        // not robustly implemented:
        s.perennial = false;

        // First lookup layout (nclass arrays) — see LandcoverTables header.
        s.xl = tab.xl[s.iclass - 1];
        s.height = tab.hmin[s.iclass - 1]
                  + (tab.hmax[s.iclass - 1] - tab.hmin[s.iclass - 1]) * s.fc;

        if (s.iswater)
        {
            s.fg = 1.0f;
            s.fc = s.fcbare;
        }
        else if (!s.perennial)
        {
            s.fg = 1.0f;
            if (s.fc < s.fcbare) s.fc = s.fcbare;
            s.xlai = -2.0f * MathF.Log(1.0f - s.fc);
            s.clump0 = 1.0f;
        }
        else
        {
            s.fg = 1.0f;
            float fcmin = tab.tabfcmin[s.iclass - 1];
            if (s.fc < fcmin)
            {
                s.fg = s.fc / fcmin;
                s.fc = fcmin;
            }
            if (s.fc < s.fcbare) s.fc = s.fcbare;
        }

        s.clump = s.clump0;
        s.clumps1 = s.clump;
        s.clumps2 = s.clump;
        s.clump0 = 1.0f;
        if (s.iclass == 12)
        {
            s.clumps1 = 0.5f;
            s.clumps2 = 0.5f;
            s.clump0 = 0.5f;
            s.clump = 0.5f;
        }

        GetFveg(s.clump0, s.xlai, out s.fveg);

        // Overwrite with one veg-file record (Fortran unit 400):
        // read(400,*) iin,jin,aleafv,aleafn,aleafl,adeadv,adeadn,adeadl,height,xl,z0,disp,rsmin2
        var v = new float[13];
        if (veg.ReadList(13, v) < 0)
        {
            s.aleafv = GridDims.Bad;
            return; // Fortran would error on short record; caller marks badinput via stream
        }
        int iin = (int)v[0]; // unused in source
        int jin = (int)v[1]; // unused in source
        s.aleafv = v[2];
        s.aleafn = v[3];
        s.aleafl = v[4];
        s.adeadv = v[5];
        s.adeadn = v[6];
        s.adeadl = v[7];
        s.height = v[8];
        s.xl = v[9];
        s.z0 = v[10];
        s.disp = v[11];
        // rsmin2 (v[12]) is an unused local in the source — not ported.

        AlexiAtmos.CanopyArch(s);
    }

    /// <summary>
    /// getfveg — solve fv*exp(-xlai/2/fv)+(1-fv) = exp(-xlai*clump0/2)
    /// by a 100-step grid search, clamped to [0.1, 1.0].
    /// </summary>
    public static void GetFveg(float clump0, float xlai, out float fveg)
    {
        float gap = MathF.Exp(-0.5f * clump0 * xlai);
        float xmin = 9999f;
        float fvmin = 0f;
        for (int i = 1; i <= 100; i++)
        {
            float fv = i / 100f;
            float rhs = fv * MathF.Exp(-0.5f * xlai / fv) + (1.0f - fv);
            float diff = MathF.Abs(gap - rhs);
            if (diff < xmin)
            {
                fvmin = fv;
                xmin = diff;
            }
        }
        fveg = fvmin;
        fveg = MathF.Min(fveg, 1.0f);
        fveg = MathF.Max(fveg, 0.1f);
    }

    private static float ParseF(string t) =>
        float.Parse(t, NumberStyles.Float, CultureInfo.InvariantCulture);

    private static int ParseInt(string t) =>
        int.Parse(t, NumberStyles.Integer, CultureInfo.InvariantCulture);
}
