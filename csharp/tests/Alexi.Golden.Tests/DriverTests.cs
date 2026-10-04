using Alexi.IO;
using Xunit;

namespace Alexi.Golden.Tests;

// Tests for Alexi.IO routines (driver/input layer).
public class DriverTests
{
    [Fact]
    public void GetDate_ParsesYYYYDDD()
    {
        int year, doy;
        AlexiDriver.GetDate(2008172, out year, out doy);
        Assert.Equal(2008, year);
        Assert.Equal(172, doy);
    }

    [Fact]
    public void GetDate_ZeroDoyKeepsYear()
    {
        int year, doy;
        AlexiDriver.GetDate(2012000, out year, out doy);
        Assert.Equal(2012, year);
        Assert.Equal(0, doy);
    }

    [Fact]
    public void GetDgmt_MountainMeridian()
    {
        float dgmt = 0f, stdlng = 0f;
        AlexiInput.GetDgmt(-105f, ref dgmt, ref stdlng);
        Assert.Equal(-7f, dgmt);
        Assert.Equal(-105f, stdlng);
    }

    [Fact]
    public void GetDgmt_CentralMeridian()
    {
        float dgmt = 0f, stdlng = 0f;
        AlexiInput.GetDgmt(-96f, ref dgmt, ref stdlng);
        Assert.Equal(-6f, dgmt);
        Assert.Equal(-90f, stdlng);
    }
}
