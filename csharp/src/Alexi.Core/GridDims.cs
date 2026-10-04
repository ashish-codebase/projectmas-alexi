namespace Alexi.Core;

/// <summary>
/// Grid and model parameters from ALEXI_parm.inc, USflux_parm.inc, flag.inc.
/// All Fortran `parameter` statements; the model is entirely real*4.
/// </summary>
public static class GridDims
{
    // ALEXI_parm.inc — PBL parameters
    public const int Mli = 41;      // number of input PBL levels
    public const int Ml = 8000;     // number of interpolated PBL levels
    public const int Mt = 8000;     // size of PBL lookup tables

    // ALEXI IERR codes
    public const float Bad = -9999f;
    public const int Success = 0;

    // USflux_parm.inc — MENAE domain grid
    public const int Kz = 30;
    public const int Ilg = 1456;
    public const int Jlg = 625;

    // GETD input grid
    public const int Kx = 1440;
    public const int Ky = 600;
    public const int Kt = 8;

    // Clear-sky INSOL grid
    public const int Sx = 1440;
    public const int Sy = 720;

    // CRAS hours with positive solar radiation, and ALEXI time-series hours
    public const int Nohrc = 24;
    public const int Nohr = 24;

    // SIB2 landcover classes (note: ICE == WATER == RIPARIAN == 1 in the source)
    public const int Nclass = 8;
    public const int Ice = 1;
    public const int Water = 1;
    public const int Riparian = 1;

    // MODIS lookup table
    public const int Nclasm = 6;
    public const int Nbin = 20;

    // NCAR graphics (unused in the port; kept for completeness)
    public const int MXPL = 50;
    public const int MNCLR = 43;
    public const int MXCLR = 250;
    public const int MXPLTS = 28;
    public const int NOUTGRAPHICS = 8;

    // flag.inc: filg/fjlg are integer*8 in Fortran; 1456/625 fit in int.
    public const int Filg = Ilg;
    public const int Fjlg = Jlg;

    /// <summary>
    /// Flattened index for a Fortran (i, j) array of shape (ilg, jlg),
    /// column-major: element (i, j) at offset (j-1)*ilg + (i-1).
    /// i, j are 1-based as in Fortran.
    /// </summary>
    public static int Index2(int i, int j) => (j - 1) * Ilg + (i - 1);

    /// <summary>
    /// Flattened index for ztan/zzan shape (kz, ilg, jlg): k is fastest-varying.
    /// </summary>
    public static int Index3Kz(int k, int i, int j) => Index2(i, j) * Kz + (k - 1);

    /// <summary>
    /// Flattened index for htht shape (mli, ilg, jlg).
    /// </summary>
    public static int Index3Mli(int k, int i, int j) => Index2(i, j) * Mli + (k - 1);

    /// <summary>
    /// Flattened index for ctloc shape (kx, ky, kt): i fastest, then j, then k.
    /// </summary>
    public static int Index3KxKyKt(int i, int j, int k) => (k - 1) * Kx * Ky + (j - 1) * Kx + (i - 1);

    /// <summary>
    /// Direct-access record number used by USflux_run.f:920, including the row flip jj = jlg - j + 1.
    /// Record numbers are 1-based; byte offset = (record - 1) * 4.
    /// </summary>
    public static long DirectRecord(int i, int j) => ((long)(Jlg - j) * Ilg) + i;
}
