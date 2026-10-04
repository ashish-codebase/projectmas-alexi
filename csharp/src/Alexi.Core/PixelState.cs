namespace Alexi.Core;

/// <summary>
/// Per-pixel model state — the Fortran COMMON blocks of ALEXI_f90.inc and date.inc.
/// Field names are kept verbatim from the Fortran source for 1:1 parity mapping.
/// All reals are float (the model is entirely real*4); implicit i-n names are int.
/// </summary>
public sealed class PixelState
{
    // common/absorption
    public float aleafn, aleafv, aleafl, adeadv, adeadn, adeadl;

    // common/aerodyn
    public float xlog1, xlog2, a, uexp1, uexp2, expuxp, psima, psih;

    // common/canopy
    public float refhtw, disp, z0, height, xl, xlai, fc, fg, clump, rsmin;

    // common/clumping
    public float clumps1, clumps2, clump0, fveg;

    // common/cover (USflux.inc:10)
    public int iclass;
    public float beta, gvs, xndvi, fpar, dstom;

    // common/cover2
    public bool perennial, iswater;
    public float fcbare;

    // common/cover3
    public bool iswater_inland;
    public int ai, aj;

    // common/constants (set by set_constants in AlexiUtl)
    public float pi, cp, xk;

    // common/data1
    public float t1, tloc1, trad1, taobs1, ea1, w1, pres1, ta1, ts1, tc1;
    public float h1, hs1, hc1, xle1, xles1, xlec1, g1, rnet1, rnsoil1, rndiv1;
    public float sdn1, par1, xlwdn1, ra1, rs1, rx1, th1, z1, rhocp1, tac1;
    public float albedo1, taubtv1, taubtn1, zen1, tb1, xlwup1, swup1;

    // common/data2
    public float t2, tloc2, trad2, taobs2, ea2, w2, pres2, ta2, ts2, tc2;
    public float h2, hs2, hc2, xle2, xles2, xlec2, g2, rnet2, rnsoil2, rndiv2;
    public float sdn2, par2, xlwdn2, ra2, rs2, rx2, th2, z2, rhocp2, tac2;
    public float albedo2, taubtv2, taubtn2, zen2, tb2, xlwup2, swup2;

    // common/data22
    public float w2orig, eref2;

    // common/emission
    public float emleaf, emdead, emsoil, emcpy;

    // common/emission2
    public float esfc, aem, bem, eleaf, esoil;

    // common/flags
    public bool writeme, badinput, converged, stopiter;

    // common/initial
    public float hn0, psi0, fc0;

    // common/model
    public float zta;

    // common/obs
    public float hobs1, hobs2, xleobs1, xleobs2, gobs1, gobs2, rnobs1, rnobs2;

    // common/partition
    public float difvis, difnir, dirvis, dirnir, fvis, fnir;

    // common/partition_clear
    public float difvisclr, difnirclr, dirvisclr, dirnirclr, fvisclr, fnirclr;

    // common/pbl
    public float thrise;
    public float[] zpbli = new float[GridDims.Mli];
    public float[] thpbli = new float[GridDims.Mli];
    public float[] zpbl = new float[GridDims.Ml];
    public float[] thpbl = new float[GridDims.Ml];
    public int jzmax, jz1, nlev;

    // common/pbl2
    public float[] tabtheta = new float[GridDims.Mt];
    public float[] tabz2 = new float[GridDims.Mt];
    public float dtheta;
    public int itmax;
    public float thpblz1;

    // common/reflection
    public float rsoiln, rsoilv, rdcpyl;

    // common/reflection2
    public float albv, albn, albvobs, albnobs, albobs;

    // common/site
    public float xlat, xlong, stdlng;

    // common/stability
    public float zdlamx, zdlamn;

    // common/timestamp (date.inc: year, doy; JDATE/MDATE are integer, CDATE/YDATE character)
    public int year, doy, jdate, mdate;
    public string cdate = new string(' ', 8);
    public string ydate = new string(' ', 17);

    // common/transmission
    public float taudn, taudv, taudl;

    // common/view
    public float theta, ftheta;
}
