using Alexi.Core;


namespace Alexi.IO;

/// <summary>
/// Fortran direct-access unformatted binary grid (recl=4, access='direct').
/// Record n maps to byte offset (n-1)*4. binwrite flips rows:
/// irec=(jlg-j)*ilg+i (GridDims.DirectRecord). binread does NOT flip
/// (irec=(j-1)*ilg+i) — a defect in the source, preserved here.
/// Byte order is native little-endian (x86), matching gfortran unformatted
/// output; the source's byteswap helpers (binread_swap/laiswapr4) are no-ops
/// on this platform. The grid is buffered in memory and flushed on close,
/// which reproduces Fortran's zero-filled direct-access file semantics.
/// </summary>
public sealed class BinaryGridWriter
{
    private readonly string path;
    private readonly byte[] data;

    public BinaryGridWriter(string path)
    {
        this.path = path;
        this.data = new byte[GridDims.Ilg * GridDims.Jlg * 4]; // zero-filled like Fortran
    }

    /// <summary>binwrite(i,j,...,val): jj=jlg-j+1; irec=(jj-1)*ilg+i.</summary>
    public void Write(int i, int j, float val)
    {
        long off = (GridDims.DirectRecord(i, j) - 1) * 4;
        System.Buffers.Binary.BinaryPrimitives.WriteSingleLittleEndian(data.AsSpan((int)off, 4), val);
    }

    ~ BinaryGridWriter()
    {
        System.IO.File.WriteAllBytes(path, data);
    }
}

public sealed class BinaryGridReader
{
    private readonly byte[] data;

    public BinaryGridReader(string path)
    {
        data = System.IO.File.ReadAllBytes(path);
    }

    /// <summary>binread(i,j,...,val): irec=(j-1)*ilg+i — no row flip (source defect).</summary>
    public float ReadUnflipped(int i, int j)
    {
        long rec = ((long)j - 1) * GridDims.Ilg + i;
        return ReadAt(rec);
    }

    /// <summary>binread_swap layout: irec=(jlg-j)*ilg+i; byteswap is a no-op on little-endian.</summary>
    public float ReadFlipped(int i, int j)
    {
        return ReadAt(GridDims.DirectRecord(i, j));
    }

    private float ReadAt(long rec)
    {
        int off = (int)((rec - 1) * 4);
        return System.Buffers.Binary.BinaryPrimitives.ReadSingleLittleEndian(data.AsSpan(off, 4));
    }
}
