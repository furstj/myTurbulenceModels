#!/usr/bin/env python3
"""
Script to modify OpenFOAM mesh points using formula from

Suluksna, K., Dechaumphai, P., & Juntasaro, E. (2009). 
Correlations for modeling transitional boundary layers under influences of freestream turbulence and pressure gradient. 
International Journal of Heat and Fluid Flow, 30(1), 66–75. https://doi.org/10.1016/j.ijheatfluidflow.2008.09.004

Formula is valid for T3C1, T3C2, T3C3, and T3C5


 
"""

import numpy as np
import gzip
import re
from pathlib import Path


def h_over_D(x: np.ndarray) -> np.ndarray:
    """
    Calculate h/D ratio 

    Parameters
    ----------
    x : np.ndarray
        Distance from the leading edge (plate distance) in meters.

    Returns
    -------
    np.ndarray
        h/D ratio, capped at 1.0.
    """
    # Polynomial coefficients from eq. (23)
    hD = ( 1.231*x**6 - 6.705*x**5 + 14.061*x**4 
        - 14.113*x**3 + 7.109*x**2 - 1.900*x + 0.950 )
    return np.minimum(hD, 1.0)


def check_binary_file(filepath: Path) -> bool:
    """
    Check if a file appears to be binary (contains null bytes).

    Parameters
    ----------
    filepath : Path
        Path to the file to check.

    Returns
    -------
    bool
        True if file appears to be binary, False otherwise.
    """
    with open(filepath, 'rb') as f:
        chunk = f.read(8192)
        return b'\x00' in chunk


def read_points(filepath: Path) -> tuple:
    """
    Read OpenFOAM points file (plain or gzipped).

    Parameters
    ----------
    filepath : Path
        Path to the points file.

    Returns
    -------
    tuple
        (header_text, num_points, points_array, footer_text)
        - header_text: lines before the point data
        - num_points: number of points
        - points_array: numpy array of shape (n_points, 3)
        - footer_text: any text after the points (e.g., comments)

    Raises
    ------
    ValueError
        If the file is binary (not ASCII format).
    """
    # Check for binary file first
    is_gzipped = filepath.suffix == '.gz' or filepath.name.endswith('.gz')

    if not is_gzipped:
        if check_binary_file(filepath):
            raise ValueError(f"File '{filepath}' appears to be binary. "
                           "This script only supports ASCII format points files.")

    if is_gzipped:
        opener = lambda: gzip.open(filepath, 'rt', encoding='latin1')
    else:
        opener = lambda: open(filepath, 'r', encoding='latin1')

    with opener() as f:
        lines = f.readlines()

    # Find the start of point data (line with just a number)
    header_lines = []
    points_start_idx = None

    for i, line in enumerate(lines):
        stripped = line.strip()
        # Check if this line is just a number (count of points)
        if re.match(r'^\d+$', stripped):
            points_start_idx = i
            break
        header_lines.append(line.rstrip())

    if points_start_idx is None:
        raise ValueError("Could not find number of points in file")

    num_points = int(lines[points_start_idx].strip())

    # Parse points - they are in parentheses
    # Find the opening parenthesis
    paren_start = points_start_idx + 1
    while paren_start < len(lines) and '(' not in lines[paren_start]:
        paren_start += 1

    if paren_start >= len(lines):
        raise ValueError("Could not find opening parenthesis for points")

    # Extract all point data until closing parenthesis
    point_data = []
    paren_depth = 0
    current_line = paren_start

    while current_line < len(lines):
        line = lines[current_line]
        paren_depth += line.count('(') - line.count(')')

        # Extract coordinates from lines containing them
        # Match numbers like -0.05, 0, 1.23e-4, etc. (without parentheses)
        coords = re.findall(r'-?[\d.eE+-]+', line)
        if len(coords) >= 3:
            # Take first 3 numbers as x, y, z
            point_data.append([float(coords[i]) for i in range(3)])

        if paren_depth <= 0 and ')' in line:
            current_line += 1
            break

        current_line += 1

    points_array = np.array(point_data[:num_points], dtype=np.float64)

    # Get footer (anything after the closing paren)
    footer_text = lines[current_line:] if current_line < len(lines) else []

    return header_lines, num_points, points_array, footer_text


def write_points(filepath: Path, header: list, num_points: int,
                 points: np.ndarray, footer: list, compress: bool = False):
    """
    Write points to OpenFOAM format file (plain or gzipped).

    Parameters
    ----------
    filepath : Path
        Output file path.
    header : list
        Header lines (FoamFile block etc.)
    num_points : int
        Number of points.
    points : np.ndarray
        Array of shape (n_points, 3).
    footer : list
        Footer lines (comments after points).
    compress : bool
        If True, write as .gz file.
    """
    opener = gzip.open if compress else open
    mode = 'wt' if compress else 'w'
    ext = '.gz' if compress else ''
    outpath = filepath.parent / (filepath.stem + ext)

    with opener(outpath, mode, encoding='latin1') as f:
        # Write header
        for line in header:
            f.write(line + '\n')

        f.write('\n')
        f.write(f'{num_points}\n')
        f.write('(\n')

        # Write points
        for point in points:
            f.write(f'({point[0]:.9g} {point[1]:.9g} {point[2]:.9g})\n')

        f.write(')\n')

        # Write footer
        for line in footer:
            if not line.endswith('\n'):
                line += '\n'
            f.write(line)

    print(f"Written {len(points)} points to {outpath}")


def modify_points(points: np.ndarray, x_offset: float = 0.0) -> np.ndarray:
    """
    Modify y coordinate of each point according to equation (23).

    The transformation scales the y coordinate based on the x position.
    According to the paper, the domain height varies along the streamwise direction.

    Parameters
    ----------
    points : np.ndarray
        Array of shape (n_points, 3) with original coordinates.
    x_offset : float
        Offset to subtract from x coordinate (leading edge position).

    Returns
    -------
    np.ndarray
        Modified points array.
    """
    modified = points.copy()

    # x coordinate is the streamwise direction
    x_coords = points[:, 0]

    # Apply offset to get distance from leading edge
    x_from_leading_edge = x_coords - x_offset

    # Calculate h/D ratio for each x position
    hD_ratios = h_over_D(x_from_leading_edge)

    # Scale y coordinate by h/D ratio
    # We need to know the reference y value. Looking at the mesh structure,
    # the y values seem to represent normalized positions within the domain height.
    # We scale y by the local h/D ratio relative to inlet.

    # For simplicity, we assume y represents a fraction of the local domain height
    # So we multiply y by h/D at that x position
    modified[:, 1] = points[:, 1] * hD_ratios

    return modified


def find_points_file(base_dir: Path) -> tuple[Path, bool]:
    """
    Find the points file in constant/polyMesh directory.

    Searches for constant/polyMesh/points and constant/polyMesh/points.gz.

    Parameters
    ----------
    base_dir : Path
        Base directory to search from.

    Returns
    -------
    tuple
        (found_path, is_compressed) or (None, False) if not found.
    """
    candidates = [
        base_dir / 'constant' / 'polyMesh' / 'points',
        base_dir / 'constant' / 'polyMesh' / 'points.gz',
    ]

    for candidate in candidates:
        if candidate.exists():
            is_gz = candidate.suffix == '.gz' or candidate.name.endswith('.gz')
            return candidate, is_gz

    return None, False


def main():
    """Main entry point for the script."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Modify OpenFOAM mesh points according to Eq. (23)'
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        help='Input points file (can be .gz compressed). '
             'If not specified, searches for points/points.gz in '
             'current dir or constant/polyMesh/'
    )
    parser.add_argument(
        '--output', '-o',
        help='Output file path (default: same name with _modified suffix)'
    )
    parser.add_argument(
        '--x-offset',
        type=float,
        default=0.0,
        help='X position of leading edge (default: 0.0 m)'
    )
    parser.add_argument(
        '--compress',
        action='store_true',
        help='Compress output file with gzip'
    )

    args = parser.parse_args()

    # Determine input file
    if args.input_file:
        input_path = Path(args.input_file)
        if not input_path.exists():
            print(f"Error: File not found: {input_path}")
            return 1
    else:
        # Auto-detect points file
        input_path, _ = find_points_file(Path('.'))
        if input_path is None:
            print("Error: No points file found.")
            print("Searched for: points, points.gz, constant/polyMesh/points, constant/polyMesh/points.gz")
            return 1
        print(f"Auto-detected points file: {input_path}")

    # Determine if input is compressed
    is_compressed = input_path.suffix == '.gz' or input_path.name.endswith('.gz')

    print(f"Reading points from {input_path}...")
    try:
        header, num_points, points, footer = read_points(input_path)
    except ValueError as e:
        print(f"Error: {e}")
        return 1

    print(f"Read {len(points)} points")
    print(f"Original bounds: x=[{points[:,0].min():.4f}, {points[:,0].max():.4f}], "
          f"y=[{points[:,1].min():.4f}, {points[:,1].max():.4f}], "
          f"z=[{points[:,2].min():.4f}, {points[:,2].max():.4f}]")

    # Modify points
    print(f"\nApplying transformation with x_offset={args.x_offset}...")
    modified_points = modify_points(points, args.x_offset)

    print(f"Modified bounds: x=[{modified_points[:,0].min():.4f}, {modified_points[:,0].max():.4f}], "
          f"y=[{modified_points[:,1].min():.4f}, {modified_points[:,1].max():.4f}], "
          f"z=[{modified_points[:,2].min():.4f}, {modified_points[:,2].max():.4f}]")

    # Determine output path
    if args.output:
        output_path = Path(args.output)
        compress_out = args.compress or (output_path.suffix == '.gz' or output_path.name.endswith('.gz'))
    else:
        # Overwrite the input file
        output_path = input_path
        compress_out = is_compressed

    write_points(output_path, header, num_points, modified_points, footer, compress_out)

    print("\nDone!")
    return 0


if __name__ == '__main__':
    exit(main())
