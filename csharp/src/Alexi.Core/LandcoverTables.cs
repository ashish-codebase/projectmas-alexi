namespace Alexi.Core;

/// <summary>
/// Landcover lookup tables — the SECOND, incompatible `common /lookup/` layout
/// declared at USflux.inc:31 (itabclass/tablai/.../tabdesc), used by
/// USflux_cover.f and USflux_run.f.
///
/// Phase 1 #1 of the conversion plan: Fortran aliases this block with
/// lookup_i/lookup_j (USflux.inc:27), which is undefined behaviour.
/// In C# the two layouts are separate structures; the aliasing is NOT
/// reproduced. Confirm intended semantics with the model author before
/// trusting reference outputs that depend on the alias.
///
/// Field names verbatim from source. tabperen is logical in Fortran
/// (used as `perennial = tabperen(iclass)`), so bool here.
/// </summary>
public sealed class LandcoverTables
{
    public int[] itabclass = new int[GridDims.Nclass];
    public float[] tablai = new float[GridDims.Nbin * GridDims.Nclasm];   // (nbin, nclasm)
    public float[] tabfpar = new float[GridDims.Nbin * GridDims.Nclasm];   // (nbin, nclasm)
    public float[] tabbeta = new float[GridDims.Nclass];
    public float[] tabaleaf = new float[GridDims.Nclass * 3];              // (nclass, 3)
    public float[] tabadead = new float[GridDims.Nclass * 3];              // (nclass, 3)
    public float[] tabhmin = new float[GridDims.Nclass];
    public float[] tabhmax = new float[GridDims.Nclass];
    public float[] tabxl = new float[GridDims.Nclass];
    public float[] tabgvs = new float[GridDims.Nclass];
    public float[] tabdstom = new float[GridDims.Nclass];
    public bool[] tabperen = new bool[GridDims.Nclass];
    public float[] tabfcmin = new float[GridDims.Nclass];
    public string[] tabdesc = new string[GridDims.Nclass];                // CHARACTER(LEN=50)

    // ---- First `common /lookup/` layout (USflux_cover.f / landcover.f) ----
    // 10 arrays of length nclass. In Fortran this block aliases the second
    // layout (undefined behaviour); in C# it is stored separately and is
    // loaded by LoadClassTable (11-column record format, landcover.f).
    public float[] alv = new float[GridDims.Nclass];
    public float[] aln = new float[GridDims.Nclass];
    public float[] all = new float[GridDims.Nclass];
    public float[] adv = new float[GridDims.Nclass];
    public float[] adn = new float[GridDims.Nclass];
    public float[] adl = new float[GridDims.Nclass];
    public float[] hmin = new float[GridDims.Nclass];
    public float[] hmax = new float[GridDims.Nclass];
    public float[] xl = new float[GridDims.Nclass];
    public float[] rs = new float[GridDims.Nclass];

    /// <summary>Index into tablai/tabfpar for 1-based (bin, class).</summary>
    public static int TableIndex(int bin, int clasm) => (clasm - 1) * GridDims.Nbin + (bin - 1);

    /// <summary>Index into tabaleaf/tabadead for 1-based (class, band).</summary>
    public static int BandIndex(int cl, int band) => (band - 1) * GridDims.Nclass + (cl - 1);
}
