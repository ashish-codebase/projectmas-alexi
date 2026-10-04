namespace Alexi.IO;

/// <summary>
/// Fortran unit 6 (screen). A minimal System.IO.TextWriter that routes
/// output to the console, for the ported routines' stdout parameters.
/// </summary>
public sealed class StdoutWriter : System.IO.TextWriter
{
    public StdoutWriter()
    {
    }

    public override System.Text.Encoding Encoding => System.Text.Encoding.UTF8;

    public override void WriteLine(string? s = null)
    {
        System.Console.WriteLine(s);
    }

    public override void Write(string? s = null)
    {
        System.Console.Write(s);
    }

    public override void Flush()
    {
    }

    public override void Close()
    {
    }
}
