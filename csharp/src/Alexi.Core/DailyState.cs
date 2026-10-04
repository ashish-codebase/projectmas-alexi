namespace Alexi.Core;

/// <summary>
/// Daily / hourly state — the Fortran COMMON blocks of USflux.inc.
/// Field names verbatim from source for parity mapping.
/// </summary>
public sealed class DailyState
{
    // common/availh2o
    public float brz, awfrz, awcrz, bsfc, awfsfc, awcsfc;

    // common/cover
    public int iclass;
    public float beta, gvs, xndvi, fpar, dstom;

    // common/dayflux
    public float aparday, acday, rnday, rnsday, rncday, gday, eday, esday, ecday;
    public float hday, hsday, hcday, sday;

    // common/dayflux2
    public float erefday;

    // common/dayflux3
    public float xlwupday, xlwdnday, swupday;

    // common/dayobsflux
    public float rnobsday, hobsday, eobsday, gobsday, aobsday, precpday;

    // common/daypotflux
    public float epotday, espotday, ecpotday;

    // common/departure
    public float xmc, xbc, xcc, xms, xbs, xcs;

    // common/hrdata
    public float[] tloc = new float[GridDims.Nohr];
    public float[] ta = new float[GridDims.Nohr];
    public float[] ea = new float[GridDims.Nohr];
    public float[] wind = new float[GridDims.Nohr];
    public float[] sdn = new float[GridDims.Nohr];
    public float[] xlwdn = new float[GridDims.Nohr];
    public float[] rnet = new float[GridDims.Nohr];
    public float[] rnsoil = new float[GridDims.Nohr];
    public float[] g = new float[GridDims.Nohr];
    public float[] pres = new float[GridDims.Nohr];
    public int nohrin;

    // common/hrdata2
    public float[] epot = new float[GridDims.Nohr];
    public float[] ecpot = new float[GridDims.Nohr];
    public float[] espot = new float[GridDims.Nohr];

    // common/hrdata3
    public float[] xlst = new float[GridDims.Nohr];

    // common/hrdata4
    public float[] eref = new float[GridDims.Nohr];

    // common/hrdata5
    public float[] xlwup = new float[GridDims.Nohr];
    public float[] swup = new float[GridDims.Nohr];

    // common/hrobs
    public float[] rnobs = new float[GridDims.Nohr];
    public float[] hobs = new float[GridDims.Nohr];
    public float[] xleobs = new float[GridDims.Nohr];
    public float[] gobs = new float[GridDims.Nohr];
    public float[] aobs = new float[GridDims.Nohr];
    public float[] precp = new float[GridDims.Nohr];

    // common/moisture
    public float faw, awf;

    // common/moisture2
    public float faw2;

    // common/stress
    public float fpet, fsdn;

    // common/sun
    public float strt, end;

    // common/usflags
    public bool clear;
}
