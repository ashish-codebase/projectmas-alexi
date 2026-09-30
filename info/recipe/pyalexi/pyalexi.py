import subprocess
import os
import argparse
import sys
import datetime
import logging
import multiprocessing


def alexi(year, doy, grid_name, npoints, part):
    date = "{}{:03d}".format(year, doy)
    print('this is how part looks in pyalexi {}'.format(part))
    # alexi_path = os.path.join(os.getcwd(), "alexi_proc")  # first time
    subprocess.call(['alexi_proc', '{}'.format(date), '{}'.format(grid_name), '{}'.format(npoints),'{}'.format(part)])


def main():
    args = arg_parse()
    if args.point is None and args.region is None:
        num_parts = multiprocessing.cpu_count()
        parts = range(num_parts)
    else:
        num_parts = 1
        parts = [0]

    doys = range(args.start_doy, args.end_doy + 1)
    years = range(args.start_year, args.end_year + 1)
    print(years)
    npoints = None
    part = 0
    alexi(args.year, args.doy, args.grid_name, npoints, part)
    # alexi(year, doy, grid_name, npoints, part)


def arg_parse():
    """"""
    parser = argparse.ArgumentParser(
        description='run ALEXI',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    default_year = datetime.date.today().year
    # default_year = 2017 # FOR TESTING
    default_today_doy = datetime.date.today().toordinal() - datetime.date(default_year - 1, 12, 31).toordinal()
    default_yesterday_doy = default_today_doy - 1
    # print(default_yesterday_doy)
    parser.add_argument(
        '--year', type=int, default=default_year,
        help='year of dataset')
    parser.add_argument(
        '--doy', type=int, default=default_yesterday_doy,
        help='start day of processing. *Note: leave blank for Real-time')
    # tiles = [60, 61, 62, 63, 64, 83, 84, 85, 86, 87, 88, 107, 108, 109, 110, 111, 112]
    parser.add_argument("--grid_name", nargs='*', type=str, default=None,
                        help="tile from 15x15 deg tile grid system")
    parser.add_argument("-c", "--corr", nargs='*', type=bool, default=False,
                        help="apply dthr correction")
    parser.add_argument(
        '-d', '--debug', default=logging.INFO, const=logging.DEBUG,
        help='Debug level logging', action='store_const', dest='loglevel')
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = arg_parse()

    logging.basicConfig(level=args.loglevel, format='%(message)s')
    logging.info('\n{0}'.format('#' * 80))
    logging.info('{0:<20s} {1}'.format(
        'Run Time Stamp:', datetime.datetime.now().isoformat(' ')))
    logging.info('{0:<20s} {1}'.format('Current Directory:', os.getcwd()))
    logging.info('{0:<20s} {1}'.format(
        'Script:', os.path.basename(sys.argv[0])))
    main()
