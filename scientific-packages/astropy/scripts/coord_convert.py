#!/usr/bin/env python3
"""
Coordinate conversion utility for astronomical coordinates.

This script provides batch coordinate transformations between different
astronomical coordinate systems using astropy.
"""

import sys
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u


def convert_coordinates(coords_input, input_frame='icrs', output_frame='galactic',
                       input_format='decimal', output_format='decimal'):
    """
    Convert astronomical coordinates between different frames.

    Parameters
    ----------
    coords_input : list of tuples or str
        Input coordinates as (lon, lat) pairs or strings
    input_frame : str
        Input coordinate frame (icrs, fk5, galactic, etc.)
    output_frame : str
        Output coordinate frame
    input_format : str
        Format of input coordinates ('decimal', 'sexagesimal', 'hourangle')
    output_format : str
        Format for output display ('decimal', 'sexagesimal', 'hourangle')

    Returns
    -------
    list
        Converted coordinates
    """
    results = []

    for coord in coords_input:
        try:
            # Parse input coordinate
            if input_format == 'decimal':
                if isinstance(coord, str):
                    parts = coord.split()
                    lon, lat = float(parts[0]), float(parts[1])
                else:
                    lon, lat = coord
                c = SkyCoord(lon*u.degree, lat*u.degree, frame=input_frame)

            elif input_format == 'sexagesimal':
                c = SkyCoord(coord, frame=input_frame, unit=(u.hourangle, u.deg))

            elif input_format == 'hourangle':
                if isinstance(coord, str):
                    parts = coord.split()
                    lon, lat = parts[0], parts[1]
                else:
                    lon, lat = coord
                c = SkyCoord(lon, lat, frame=input_frame, unit=(u.hourangle, u.deg))

            # Transform to output frame
            if output_frame == 'icrs':
                c_out = c.icrs
            elif output_frame == 'fk5':
                c_out = c.fk5
            elif output_frame == 'fk4':
                c_out = c.fk4
            elif output_frame == 'galactic':
                c_out = c.galactic
            elif output_frame == 'supergalactic':
                c_out = c.supergalactic
            else:
                c_out = c.transform_to(output_frame)

            results.append(c_out)

        except Exception as e:
            print(f"Error converting coordinate {coord}: {e}", file=sys.stderr)
            results.append(None)

    return results


def format_output(coords, frame, output_format='decimal'):
    """Format coordinates for display."""
    output = []

    for c in coords:
        if c is None:
            output.append("ERROR")
            continue

        if frame in ['icrs', 'fk5', 'fk4']:
            lon_name, lat_name = 'RA', 'Dec'
            lon = c.ra
            lat = c.dec
        elif frame == 'galactic':
            lon_name, lat_name = 'l', 'b'
            lon = c.l
            lat = c.b
        elif frame == 'supergalactic':
            lon_name, lat_name = 'sgl', 'sgb'
            lon = c.sgl
            lat = c.sgb
        else:
            lon_name, lat_name = 'lon', 'lat'
            lon = c.spherical.lon
            lat = c.spherical.lat

        if output_format == 'decimal':
            out_str = f"{lon.degree:12.8f} {lat.degree:+12.8f}"
        elif output_format == 'sexagesimal':
            if frame in ['icrs', 'fk5', 'fk4']:
                out_str = f"{lon.to_string(unit=u.hourangle, sep=':', pad=True)} "
                out_str += f"{lat.to_string(unit=u.degree, sep=':', pad=True)}"
            else:
                out_str = f"{lon.to_string(unit=u.degree, sep=':', pad=True)} "
                out_str += f"{lat.to_string(unit=u.degree, sep=':', pad=True)}"
        elif output_format == 'hourangle':
            out_str = f"{lon.to_string(unit=u.hourangle, sep=' ', pad=True)} "
            out_str += f"{lat.to_string(unit=u.degree, sep=' ', pad=True)}"

        output.append(out_str)

    return output


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Convert astronomical coordinates between different frames',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Supported frames: icrs, fk5, fk4, galactic, supergalactic

Input formats:
  decimal     : Degrees (e.g., "10.68 41.27")
  sexagesimal : HMS/DMS (e.g., "00:42:44.3 +41:16:09")
  hourangle   : Hours and degrees (e.g., "10.5h 41.5d")

Examples:
  %(prog)s --from icrs --to galactic "10.68 41.27"
  %(prog)s --from icrs --to galactic --input decimal --output sexagesimal "150.5 -30.2"
  %(prog)s --from galactic --to icrs "120.5 45.3"
  %(prog)s --file coords.txt --from icrs --to galactic
        """
    )

    parser.add_argument('coordinates', nargs='*',
                       help='Coordinates to convert (lon lat pairs)')
    parser.add_argument('-f', '--from', dest='input_frame', default='icrs',
                       help='Input coordinate frame (default: icrs)')
    parser.add_argument('-t', '--to', dest='output_frame', default='galactic',
                       help='Output coordinate frame (default: galactic)')
    parser.add_argument('-i', '--input', dest='input_format', default='decimal',
                       choices=['decimal', 'sexagesimal', 'hourangle'],
                       help='Input format (default: decimal)')
    parser.add_argument('-o', '--output', dest='output_format', default='decimal',
                       choices=['decimal', 'sexagesimal', 'hourangle'],
                       help='Output format (default: decimal)')
    parser.add_argument('--file', dest='input_file',
                       help='Read coordinates from file (one per line)')
    parser.add_argument('--header', action='store_true',
                       help='Print header line with coordinate names')

    args = parser.parse_args()

    # Get coordinates from file or command line
    if args.input_file:
        try:
            with open(args.input_file, 'r') as f:
                coords = [line.strip() for line in f if line.strip()]
        except FileNotFoundError:
            print(f"Error: File '{args.input_file}' not found.", file=sys.stderr)
            sys.exit(1)
    else:
        if not args.coordinates:
            print("Error: No coordinates provided.", file=sys.stderr)
            parser.print_help()
            sys.exit(1)

        # Combine pairs of arguments
        if args.input_format == 'decimal':
            coords = []
            i = 0
            while i < len(args.coordinates):
                if i + 1 < len(args.coordinates):
                    coords.append(f"{args.coordinates[i]} {args.coordinates[i+1]}")
                    i += 2
                else:
                    print(f"Warning: Odd number of coordinates, skipping last value",
                          file=sys.stderr)
                    break
        else:
            coords = args.coordinates

    # Convert coordinates
    converted = convert_coordinates(coords,
                                   input_frame=args.input_frame,
                                   output_frame=args.output_frame,
                                   input_format=args.input_format,
                                   output_format=args.output_format)

    # Format and print output
    formatted = format_output(converted, args.output_frame, args.output_format)

    # Print header if requested
    if args.header:
        if args.output_frame in ['icrs', 'fk5', 'fk4']:
            if args.output_format == 'decimal':
                print(f"{'RA (deg)':>12s} {'Dec (deg)':>13s}")
            else:
                print(f"{'RA':>25s} {'Dec':>26s}")
        elif args.output_frame == 'galactic':
            if args.output_format == 'decimal':
                print(f"{'l (deg)':>12s} {'b (deg)':>13s}")
            else:
                print(f"{'l':>25s} {'b':>26s}")

    for line in formatted:
        print(line)


if __name__ == '__main__':
    main()
