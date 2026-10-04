using Alexi.Core;
using Alexi.IO;

// CLI front end. Argument contract from USflux_main.f:
//   getarg(1)=MDATE (YYYYDDD), getarg(2)=t (grid tag), getarg(3)=npoints, getarg(4)=part
// Optional path overrides (replacing the hardcoded HOMEDIR/INDIR/OUTDIR paths):
//   --input-dir D  --output-dir D  --table-dir D  --log-dir D  --albedo-dir D  --lai-dir D
//   --store-binary   (enable store_output_bin; commented out in the shipped driver)
// Subcommand: interp-lai <inputDir>  (port of interp_LAI.f)
if (args.Length >= 1 && args[0] == "interp-lai")
{
    string dir = args.Length >= 2 ? args[1] : "./";
    InterpLAI.Run(dir);
    return 0;
}

if (args.Length < 4)
{
    Console.Error.WriteLine("usage: alexi <date_YYYYDDD> <grid_tag> <npoints> <part> [--input-dir D ...] [--store-binary]");
    Console.Error.WriteLine("       alexi interp-lai <inputDir>");
    return 2;
}

int mdate = int.Parse(args[0]);
string t = args[1];
int npoints = int.Parse(args[2]);
string part = args[3];

// Defaults mirror USflux_dir.inc / the shipped driver's hardcoded paths.
string inputDir = "./INPUTS/ALEXI_INPUT/";
string outputDir = "./OUTPUTS/ALEXI_OUTPUT/";
string tableDir = "info/recipe/ALEXI_0.1/";
string logDir = "./ALEXI_LOG";
string albedoDir = "./";
string laiDir = "./";
bool storeBinary = false;

int k = 4;
while (k < args.Length)
{
    switch (args[k])
    {
        case "--input-dir":
            inputDir = args[++k];
            break;
        case "--output-dir":
            outputDir = args[++k];
            break;
        case "--table-dir":
            tableDir = args[++k];
            break;
        case "--log-dir":
            logDir = args[++k];
            break;
        case "--albedo-dir":
            albedoDir = args[++k];
            break;
        case "--lai-dir":
            laiDir = args[++k];
            break;
        case "--store-binary":
            storeBinary = true;
            break;
        default:
            Console.Error.WriteLine($"unknown option: {args[k]}");
            return 2;
    }
    k++;
}

AlexiPaths paths = new AlexiPaths(inputDir, outputDir, tableDir, logDir, albedoDir, laiDir);

PixelState s = new PixelState();
DailyState d = new DailyState();
LandcoverTables tab = new LandcoverTables();
PixelStreams io = new PixelStreams();
StdoutWriter stdout = new StdoutWriter();

// main (USflux_main.f): landcover stage, then USflux
AlexiDriver.Landcover(mdate, t, npoints, part, paths, tab, stdout);
io.stdout = stdout;
RunCounters c = AlexiDriver.USflux(mdate, t, npoints, part, paths, s, d, tab, io, storeBinary);

Console.WriteLine($"done: ntot={c.ntot} nbad={c.nbad} nwater={c.nwater} ncloud={c.ncloud} nconv={c.nconv} nfail={c.nfail}");
return 0;
