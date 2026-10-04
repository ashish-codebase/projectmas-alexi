using System.Globalization;

namespace Alexi.Core;

/// <summary>
/// Fortran list-directed sequential reader — READ(unit,*) semantics.
/// A single READ consumes whole records (lines) until all items are
/// satisfied; any tokens left over on the last consumed record are
/// discarded. ReadList returns a Fortran-style iostat: 0 = OK, &lt; 0 = EOF.
/// </summary>
public sealed class ListDirectedReader
{
    private readonly StreamReader stream;

    public ListDirectedReader(string path)
    {
        stream = new StreamReader(path);
    }

    public int ReadList(int count, float[] vals)
    {
        int got = 0;
        while (got < count)
        {
            string? line = stream.ReadLine();
            if (line == null)
                return -1;
            foreach (string tok in Tokens(line))
            {
                if (got >= count)
                    break; // list-directed READ discards the rest of the record
                vals[got++] = ParseF(tok);
            }
        }
        return 0;
    }

    private static IEnumerable<string> Tokens(string line) =>
        line.Split([' ', '\t', ','], StringSplitOptions.RemoveEmptyEntries);

    private static float ParseF(string t) =>
        float.Parse(t, NumberStyles.Float, CultureInfo.InvariantCulture);
}
