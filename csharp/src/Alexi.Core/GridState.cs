namespace Alexi.Core;

/// <summary>
/// Grid-scale arrays from USflux_grids.inc, flag.inc, and the met arrays of USflux.inc.
/// Flattened, column-major (Fortran layout preserved): see GridDims.Index* helpers.
/// Total heap ≈ 390 MB; allocated once per run, not per pixel.
/// </summary>
public sealed class GridState
{
    // common/canopy1 — satellite-derived surface inputs (ilg, jlg)
    public float[] tsr15 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] tsr55 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] vsr15 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] vsr55 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] wsr15 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] wsr55 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] psr15 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] psr55 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] radini = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] radfin = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] ztan = new float[GridDims.Kz * GridDims.Ilg * GridDims.Jlg];
    public float[] zzan = new float[GridDims.Kz * GridDims.Ilg * GridDims.Jlg];
    public float[] htht = new float[GridDims.Mli * GridDims.Ilg * GridDims.Jlg];
    public float[] pran = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] lscls = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] swbeg = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] swend = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] xlwbeg = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] xlwend = new float[GridDims.Ilg * GridDims.Jlg];

    // common/canopy2
    public float[] satang = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] rise15 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] rise55 = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] navlat = new float[GridDims.Jlg];
    public float[] navlon = new float[GridDims.Ilg];
    public float[] mlai = new float[GridDims.Ilg * GridDims.Jlg];

    // common/canopy3 — domain scalars
    public float clat, clon, dx, dy, minlat, minlon, maxlat, maxlon;

    // common/fflag (flag.inc; implicit f -> real*4)
    public float[] flag = new float[GridDims.Filg * GridDims.Fjlg];

    // common/hrdata_met (GETD met grids, kx, ky, kt)
    public float[] ctloc = new float[GridDims.Kx * GridDims.Ky * GridDims.Kt];
    public float[] cta = new float[GridDims.Kx * GridDims.Ky * GridDims.Kt];
    public float[] cea = new float[GridDims.Kx * GridDims.Ky * GridDims.Kt];
    public float[] cwind = new float[GridDims.Kx * GridDims.Ky * GridDims.Kt];
    public float[] csdn = new float[GridDims.Kx * GridDims.Ky * GridDims.Kt];
    public float[] cxlwdn = new float[GridDims.Kx * GridDims.Ky * GridDims.Kt];
    public float[] cpres = new float[GridDims.Kx * GridDims.Ky * GridDims.Kt];
    public float[] clst = new float[GridDims.Kx * GridDims.Ky * GridDims.Kt];
    public float[] clapse = new float[GridDims.Ilg * GridDims.Jlg];

    // common/slookup
    public float[] insol_i = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] insol_j = new float[GridDims.Ilg * GridDims.Jlg];

    // common/lookup, layout A (pbl_read.f, sfc_read.f, USflux_run.f):
    // lookup_i(ilg, jlg), lookup_j(ilg, jlg) — pixel-to-met-grid lookup.
    public float[] lookup_i = new float[GridDims.Ilg * GridDims.Jlg];
    public float[] lookup_j = new float[GridDims.Ilg * GridDims.Jlg];
}
