using System.Globalization;
using Alexi.Core;

namespace Alexi.IO;

/// <summary>
/// Port of ALEXI_0.1/src/interp_LAI.f — standalone LAI grid interpolation tool.
/// Reads <inputDir>/interp.in ("file1, file2, outfile, f, ilg, jlg"), then
/// linearly blends two direct-access float32 grids (recl=4) into an output grid.
/// </summary>
public static class InterpLAI
{
    public static void Run(string inputDir)
    {
        const int ilgmx = 2000, jlgmx = 2000;
        const float Bad = GridDims.Bad;

        // INPUT contains no embedded blanks, so the index(INPUT,' ')-1 truncation is a no-op
        string filei = inputDir + "interp.in";
        string line = System.IO.File.ReadAllText(filei);

        // read (20,*) file1,file2,outfile,f,ilg,jlg — list-directed: comma/whitespace separated
        string[] parts = line.Trim().Split([',', ' ', '\t', '\r', '\n'], StringSplitOptions.RemoveEmptyEntries);
        string file1 = parts[0], file2 = parts[1], outfile = parts[2];
        float f = float.Parse(parts[3], CultureInfo.InvariantCulture);
        int ilg = int.Parse(parts[4], CultureInfo.InvariantCulture);
        int jlg = int.Parse(parts[5], CultureInfo.InvariantCulture);

        if (ilg > ilgmx || jlg > jlgmx)
        {
            System.Console.WriteLine("Need to increase ILGMX JLGMX");
            return; // stop
        }

        // direct-access unformatted recl=4: record n lives at byte offset (n-1)*4
        byte[] b1 = System.IO.File.ReadAllBytes(file1);
        byte[] b2 = System.IO.File.ReadAllBytes(file2);
        var outBuf = new byte[ilg * jlg * 4];

        for (int ja = 1; ja <= jlg; ja++)
        {
            for (int ia = 1; ia <= ilg; ia++)
            {
                int irec = (ja - 1) * ilg + ia;
                int off = (irec - 1) * 4;
                float v1 = System.Buffers.Binary.BinaryPrimitives.ReadSingleLittleEndian(b1.AsSpan(off, 4));
                float v2 = System.Buffers.Binary.BinaryPrimitives.ReadSingleLittleEndian(b2.AsSpan(off, 4));
                float val;
                if (v1 >= -0.1f && v2 >= -0.1f)
                {
                    val = v1 + f * (v2 - v1);
                    if (val <= 0.0f) val = 0.0f;
                }
                else
                {
                    val = Bad;
                }
                System.Buffers.Binary.BinaryPrimitives.WriteSingleLittleEndian(outBuf.AsSpan(off, 4), val);
            }
        }

        System.IO.File.WriteAllBytes(outfile, outBuf);
        // close 20/21/22 + stop: implicit on return
    }
}
