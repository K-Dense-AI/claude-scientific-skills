#!/usr/bin/env python3
"""
Quick FITS file inspection tool.

This script provides a convenient way to inspect FITS file structure,
headers, and basic statistics without writing custom code each time.
"""

import sys
from pathlib import Path
from astropy.io import fits
import numpy as np


def print_fits_info(filename, detailed=False, ext=None):
    """
    Print comprehensive information about a FITS file.

    Parameters
    ----------
    filename : str
        Path to FITS file
    detailed : bool
        If True, print detailed statistics for each HDU
    ext : int or str, optional
        Specific extension to examine in detail
    """
    print(f"\n{'='*70}")
    print(f"FITS File: {filename}")
    print(f"{'='*70}\n")

    try:
        with fits.open(filename) as hdul:
            # Print file structure
            print("File Structure:")
            print("-" * 70)
            hdul.info()
            print()

            # If specific extension requested
            if ext is not None:
                print(f"\nDetailed view of extension: {ext}")
                print("-" * 70)
                hdu = hdul[ext]
                print_hdu_details(hdu, detailed=True)
                return

            # Print header and data info for each HDU
            for i, hdu in enumerate(hdul):
                print(f"\n{'='*70}")
                print(f"HDU {i}: {hdu.name}")
                print(f"{'='*70}")
                print_hdu_details(hdu, detailed=detailed)

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FITS file: {e}")
        sys.exit(1)


def print_hdu_details(hdu, detailed=False):
    """Print details for a single HDU."""

    # Header information
    print("\nHeader Information:")
    print("-" * 70)

    # Key header keywords
    important_keywords = ['SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND',
                         'OBJECT', 'TELESCOP', 'INSTRUME', 'OBSERVER',
                         'DATE-OBS', 'EXPTIME', 'FILTER', 'AIRMASS',
                         'RA', 'DEC', 'EQUINOX', 'CTYPE1', 'CTYPE2']

    header = hdu.header
    for key in important_keywords:
        if key in header:
            value = header[key]
            comment = header.comments[key]
            print(f"  {key:12s} = {str(value):20s} / {comment}")

    # NAXIS keywords
    if 'NAXIS' in header:
        naxis = header['NAXIS']
        for i in range(1, naxis + 1):
            key = f'NAXIS{i}'
            if key in header:
                print(f"  {key:12s} = {str(header[key]):20s} / {header.comments[key]}")

    # Data information
    if hdu.data is not None:
        print("\nData Information:")
        print("-" * 70)

        data = hdu.data
        print(f"  Data type: {data.dtype}")
        print(f"  Shape: {data.shape}")

        # For image data
        if hasattr(data, 'ndim') and data.ndim >= 1:
            try:
                # Calculate statistics
                finite_data = data[np.isfinite(data)]
                if len(finite_data) > 0:
                    print(f"  Min: {np.min(finite_data):.6g}")
                    print(f"  Max: {np.max(finite_data):.6g}")
                    print(f"  Mean: {np.mean(finite_data):.6g}")
                    print(f"  Median: {np.median(finite_data):.6g}")
                    print(f"  Std: {np.std(finite_data):.6g}")

                    # Count special values
                    n_nan = np.sum(np.isnan(data))
                    n_inf = np.sum(np.isinf(data))
                    if n_nan > 0:
                        print(f"  NaN values: {n_nan}")
                    if n_inf > 0:
                        print(f"  Inf values: {n_inf}")
            except Exception as e:
                print(f"  Could not calculate statistics: {e}")

        # For table data
        if hasattr(data, 'columns'):
            print(f"\n  Table Columns ({len(data.columns)}):")
            for col in data.columns:
                print(f"    {col.name:20s} {col.format:10s} {col.unit or ''}")

            if detailed:
                print(f"\n  First few rows:")
                print(data[:min(5, len(data))])
    else:
        print("\n  No data in this HDU")

    # WCS information if present
    try:
        from astropy.wcs import WCS
        wcs = WCS(hdu.header)
        if wcs.has_celestial:
            print("\nWCS Information:")
            print("-" * 70)
            print(f"  Has celestial WCS: Yes")
            print(f"  CTYPE: {wcs.wcs.ctype}")
            if wcs.wcs.crval is not None:
                print(f"  CRVAL: {wcs.wcs.crval}")
            if wcs.wcs.crpix is not None:
                print(f"  CRPIX: {wcs.wcs.crpix}")
            if wcs.wcs.cdelt is not None:
                print(f"  CDELT: {wcs.wcs.cdelt}")
    except Exception:
        pass


def main():
    """Main function for command-line usage."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Inspect FITS file structure and contents',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s image.fits
  %(prog)s image.fits --detailed
  %(prog)s image.fits --ext 1
  %(prog)s image.fits --ext SCI
        """
    )

    parser.add_argument('filename', help='FITS file to inspect')
    parser.add_argument('-d', '--detailed', action='store_true',
                       help='Show detailed statistics for each HDU')
    parser.add_argument('-e', '--ext', type=str, default=None,
                       help='Show details for specific extension only (number or name)')

    args = parser.parse_args()

    # Convert extension to int if numeric
    ext = args.ext
    if ext is not None:
        try:
            ext = int(ext)
        except ValueError:
            pass  # Keep as string for extension name

    print_fits_info(args.filename, detailed=args.detailed, ext=ext)


if __name__ == '__main__':
    main()
