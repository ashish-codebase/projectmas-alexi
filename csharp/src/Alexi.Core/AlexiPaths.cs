namespace Alexi.Core;

/// <summary>
/// Injected directory configuration replacing the hardcoded paths in
/// USflux_dir.inc (HOMEDIR/INDIR/TABDIR/OUTDIR) and the /data/... paths
/// in USflux_run.f, sfc_read.f, landcover.f. Never resolved from CWD.
/// </summary>
public sealed record AlexiPaths(
    string InputDir,
    string OutputDir,
    string TableDir,
    string LogDir,
    string AlbedoDir,
    string LaiDir);

/// <summary>Run options for one model invocation (pyalexi CLI contract).</summary>
public sealed record AlexiRunOptions(
    int Mdate,
    string TimeTag,
    string Part,
    int NPoints,
    AlexiPaths Paths,
    int GridI = GridDims.Ilg,
    int GridJ = GridDims.Jlg);
