#!/usr/bin/env python3
"""
Cosmological calculator using astropy.cosmology.

This script provides quick calculations of cosmological distances,
ages, and other quantities for given redshifts.
"""

import sys
import argparse
import numpy as np
from astropy.cosmology import FlatLambdaCDM, Planck18, Planck15, WMAP9
import astropy.units as u


def calculate_cosmology(redshifts, cosmology='Planck18', H0=None, Om0=None):
    """
    Calculate cosmological quantities for given redshifts.

    Parameters
    ----------
    redshifts : array-like
        Redshift values
    cosmology : str
        Cosmology to use ('Planck18', 'Planck15', 'WMAP9', 'custom')
    H0 : float, optional
        Hubble constant for custom cosmology (km/s/Mpc)
    Om0 : float, optional
        Matter density parameter for custom cosmology

    Returns
    -------
    dict
        Dictionary containing calculated quantities
    """
    # Select cosmology
    if cosmology == 'Planck18':
        cosmo = Planck18
    elif cosmology == 'Planck15':
        cosmo = Planck15
    elif cosmology == 'WMAP9':
        cosmo = WMAP9
    elif cosmology == 'custom':
        if H0 is None or Om0 is None:
            raise ValueError("Must provide H0 and Om0 for custom cosmology")
        cosmo = FlatLambdaCDM(H0=H0 * u.km/u.s/u.Mpc, Om0=Om0)
    else:
        raise ValueError(f"Unknown cosmology: {cosmology}")

    z = np.atleast_1d(redshifts)

    results = {
        'redshift': z,
        'cosmology': str(cosmo),
        'luminosity_distance': cosmo.luminosity_distance(z),
        'angular_diameter_distance': cosmo.angular_diameter_distance(z),
        'comoving_distance': cosmo.comoving_distance(z),
        'comoving_volume': cosmo.comoving_volume(z),
        'age': cosmo.age(z),
        'lookback_time': cosmo.lookback_time(z),
        'H': cosmo.H(z),
        'scale_factor': 1.0 / (1.0 + z)
    }

    return results, cosmo


def print_results(results, verbose=False, csv=False):
    """Print calculation results."""

    z = results['redshift']

    if csv:
        # CSV output
        print("z,D_L(Mpc),D_A(Mpc),D_C(Mpc),Age(Gyr),t_lookback(Gyr),H(km/s/Mpc)")
        for i in range(len(z)):
            print(f"{z[i]:.6f},"
                  f"{results['luminosity_distance'][i].value:.6f},"
                  f"{results['angular_diameter_distance'][i].value:.6f},"
                  f"{results['comoving_distance'][i].value:.6f},"
                  f"{results['age'][i].value:.6f},"
                  f"{results['lookback_time'][i].value:.6f},"
                  f"{results['H'][i].value:.6f}")
    else:
        # Formatted table output
        if verbose:
            print(f"\nCosmology: {results['cosmology']}")
            print("-" * 80)

        print(f"\n{'z':>8s} {'D_L':>12s} {'D_A':>12s} {'D_C':>12s} "
              f"{'Age':>10s} {'t_lb':>10s} {'H(z)':>10s}")
        print(f"{'':>8s} {'(Mpc)':>12s} {'(Mpc)':>12s} {'(Mpc)':>12s} "
              f"{'(Gyr)':>10s} {'(Gyr)':>10s} {'(km/s/Mpc)':>10s}")
        print("-" * 80)

        for i in range(len(z)):
            print(f"{z[i]:8.4f} "
                  f"{results['luminosity_distance'][i].value:12.3f} "
                  f"{results['angular_diameter_distance'][i].value:12.3f} "
                  f"{results['comoving_distance'][i].value:12.3f} "
                  f"{results['age'][i].value:10.4f} "
                  f"{results['lookback_time'][i].value:10.4f} "
                  f"{results['H'][i].value:10.4f}")

        if verbose:
            print("\nLegend:")
            print("  z    : Redshift")
            print("  D_L  : Luminosity distance")
            print("  D_A  : Angular diameter distance")
            print("  D_C  : Comoving distance")
            print("  Age  : Age of universe at z")
            print("  t_lb : Lookback time to z")
            print("  H(z) : Hubble parameter at z")


def convert_quantity(value, quantity_type, cosmo, to_redshift=False):
    """
    Convert between redshift and cosmological quantity.

    Parameters
    ----------
    value : float
        Value to convert
    quantity_type : str
        Type of quantity ('luminosity_distance', 'age', etc.)
    cosmo : Cosmology
        Cosmology object
    to_redshift : bool
        If True, convert quantity to redshift; else convert z to quantity
    """
    from astropy.cosmology import z_at_value

    if to_redshift:
        # Convert quantity to redshift
        if quantity_type == 'luminosity_distance':
            z = z_at_value(cosmo.luminosity_distance, value * u.Mpc)
        elif quantity_type == 'age':
            z = z_at_value(cosmo.age, value * u.Gyr)
        elif quantity_type == 'lookback_time':
            z = z_at_value(cosmo.lookback_time, value * u.Gyr)
        elif quantity_type == 'comoving_distance':
            z = z_at_value(cosmo.comoving_distance, value * u.Mpc)
        else:
            raise ValueError(f"Unknown quantity type: {quantity_type}")
        return z
    else:
        # Convert redshift to quantity
        if quantity_type == 'luminosity_distance':
            return cosmo.luminosity_distance(value)
        elif quantity_type == 'age':
            return cosmo.age(value)
        elif quantity_type == 'lookback_time':
            return cosmo.lookback_time(value)
        elif quantity_type == 'comoving_distance':
            return cosmo.comoving_distance(value)
        else:
            raise ValueError(f"Unknown quantity type: {quantity_type}")


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Calculate cosmological quantities for given redshifts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available cosmologies: Planck18, Planck15, WMAP9, custom

Examples:
  %(prog)s 0.5 1.0 1.5
  %(prog)s 0.5 --cosmology Planck15
  %(prog)s 0.5 --cosmology custom --H0 70 --Om0 0.3
  %(prog)s --range 0 3 0.5
  %(prog)s 0.5 --verbose
  %(prog)s 0.5 1.0 --csv
  %(prog)s --convert 1000 --from luminosity_distance --cosmology Planck18
        """
    )

    parser.add_argument('redshifts', nargs='*', type=float,
                       help='Redshift values to calculate')
    parser.add_argument('-c', '--cosmology', default='Planck18',
                       choices=['Planck18', 'Planck15', 'WMAP9', 'custom'],
                       help='Cosmology to use (default: Planck18)')
    parser.add_argument('--H0', type=float,
                       help='Hubble constant for custom cosmology (km/s/Mpc)')
    parser.add_argument('--Om0', type=float,
                       help='Matter density parameter for custom cosmology')
    parser.add_argument('-r', '--range', nargs=3, type=float, metavar=('START', 'STOP', 'STEP'),
                       help='Generate redshift range (start stop step)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Print verbose output with cosmology details')
    parser.add_argument('--csv', action='store_true',
                       help='Output in CSV format')
    parser.add_argument('--convert', type=float,
                       help='Convert a quantity to redshift')
    parser.add_argument('--from', dest='from_quantity',
                       choices=['luminosity_distance', 'age', 'lookback_time', 'comoving_distance'],
                       help='Type of quantity to convert from')

    args = parser.parse_args()

    # Handle conversion mode
    if args.convert is not None:
        if args.from_quantity is None:
            print("Error: Must specify --from when using --convert", file=sys.stderr)
            sys.exit(1)

        # Get cosmology
        if args.cosmology == 'Planck18':
            cosmo = Planck18
        elif args.cosmology == 'Planck15':
            cosmo = Planck15
        elif args.cosmology == 'WMAP9':
            cosmo = WMAP9
        elif args.cosmology == 'custom':
            if args.H0 is None or args.Om0 is None:
                print("Error: Must provide --H0 and --Om0 for custom cosmology",
                      file=sys.stderr)
                sys.exit(1)
            cosmo = FlatLambdaCDM(H0=args.H0 * u.km/u.s/u.Mpc, Om0=args.Om0)

        z = convert_quantity(args.convert, args.from_quantity, cosmo, to_redshift=True)
        print(f"\n{args.from_quantity.replace('_', ' ').title()} = {args.convert}")
        print(f"Redshift z = {z:.6f}")
        print(f"(using {args.cosmology} cosmology)")
        return

    # Get redshifts
    if args.range:
        start, stop, step = args.range
        redshifts = np.arange(start, stop + step/2, step)
    elif args.redshifts:
        redshifts = np.array(args.redshifts)
    else:
        print("Error: No redshifts provided.", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    # Calculate
    try:
        results, cosmo = calculate_cosmology(redshifts, args.cosmology,
                                            H0=args.H0, Om0=args.Om0)
        print_results(results, verbose=args.verbose, csv=args.csv)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
