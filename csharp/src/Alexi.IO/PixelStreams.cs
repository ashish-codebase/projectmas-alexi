using Alexi.Core;

namespace Alexi.IO;

/// <summary>
/// Per-run text/binary streams, mirroring the Fortran unit numbers used by
/// USflux_run.f and USflux.f90. Fields are assigned when the run opens them.
/// </summary>
public sealed class PixelStreams
{
    public ListDirectedReader? met = null; // unit 102
    public ListDirectedReader? sat = null; // unit 110
    public ListDirectedReader? nparm = null; // unit 120
    public ListDirectedReader? profile = null; // unit 130
    public ListDirectedReader? veg = null; // unit 400 (veg_us written by the landcover stage)

    public TextWriter? stdout = null; // unit 6 (console)
    public TextWriter? log99 = null; // unit 99 (per-run LOG file)
    public TextWriter? output103 = null; // unit 103 (per-pixel summary)

    public readonly Dictionary<int, BinaryGridWriter> binary = new Dictionary<int, BinaryGridWriter>();

    /// Saved iostat of extract_input's first READ (Fortran SAVE on locals);
    /// the sat/nparm/profile READs jump to the EOF label without setting it.
    public int test;

    /// Saved implicit local dgmt in extract_input (Fortran SAVE).
    public float dgmt;
}
