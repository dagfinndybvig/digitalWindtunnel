#!/usr/bin/env python3
"""
wing_to_stl.py  —  Build an extruded, lofted 3-D wing from Eppler profile.out files.

Usage
-----
    python wing_to_stl.py root_profile.out [tip_profile.out] [options]

When only a root profile is given the same airfoil shape is used from root to
tip (only chord, twist, sweep and dihedral vary).  Supply an optional tip
profile to blend smoothly between two different airfoil sections.
"""

import argparse
import math
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
Point2D = tuple[float, float]
Point3D = tuple[float, float, float]
Triangle3D = tuple[Point3D, Point3D, Point3D]


# ---------------------------------------------------------------------------
# Profile parsing  (reads the coordinate table from profile.out)
# ---------------------------------------------------------------------------

def parse_profile_out(path: Path) -> tuple[list[Point2D], str]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    airfoil_name = path.stem
    for line in lines:
        if "AIRFOIL" in line and "THICKNESS" in line:
            parts = line.split("AIRFOIL", 1)[1].split("%", 1)[0].strip().split()
            if parts:
                airfoil_name = parts[0]
            break

    header_idx = None
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("N") and "X" in s and "Y" in s:
            header_idx = i
            break

    if header_idx is None:
        raise ValueError(f"{path}: cannot find coordinate table header 'N X Y …'")

    points: list[Point2D] = []
    for line in lines[header_idx + 1:]:
        parts = line.strip().split()
        if len(parts) < 3:
            break
        try:
            int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
        except ValueError:
            break
        points.append((x, y))

    if len(points) < 3:
        raise ValueError(f"{path}: not enough coordinate points found.")

    # Remove duplicate closing point if present
    if _almost_same(points[0], points[-1]):
        points.pop()

    return points, airfoil_name


def _almost_same(a: Point2D, b: Point2D, tol: float = 1e-9) -> bool:
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol


# ---------------------------------------------------------------------------
# 2-D polygon helpers
# ---------------------------------------------------------------------------

def _polygon_area(pts: list[Point2D]) -> float:
    """Signed shoelace area; positive = CCW."""
    n = len(pts)
    s = 0.0
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return 0.5 * s


def _cross2(a: Point2D, b: Point2D, c: Point2D) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _point_in_triangle(p: Point2D, a: Point2D, b: Point2D, c: Point2D) -> bool:
    eps = 1e-12
    c1, c2, c3 = _cross2(a, b, p), _cross2(b, c, p), _cross2(c, a, p)
    return not ((c1 < -eps or c2 < -eps or c3 < -eps) and
                (c1 > eps or c2 > eps or c3 > eps))


def _ear_clip(pts_ccw: list[Point2D]) -> list[tuple[int, int, int]]:
    """Ear-clipping triangulation of a simple CCW polygon."""
    n = len(pts_ccw)
    if n < 3:
        raise ValueError("Need at least 3 points.")
    indices = list(range(n))
    tris: list[tuple[int, int, int]] = []
    eps = 1e-12
    max_loops = n * n
    loops = 0
    while len(indices) > 3:
        m = len(indices)
        ear_found = False
        for i in range(m):
            ip = indices[(i - 1) % m]
            ic = indices[i]
            ix = indices[(i + 1) % m]
            a, b, c = pts_ccw[ip], pts_ccw[ic], pts_ccw[ix]
            if _cross2(a, b, c) <= eps:
                continue
            if any(_point_in_triangle(pts_ccw[j], a, b, c)
                   for j in indices if j not in (ip, ic, ix)):
                continue
            tris.append((ip, ic, ix))
            del indices[i]
            ear_found = True
            break
        loops += 1
        if not ear_found or loops > max_loops:
            raise ValueError("Ear-clip triangulation failed; polygon may be degenerate.")
    tris.append((indices[0], indices[1], indices[2]))
    return tris


# ---------------------------------------------------------------------------
# 3-D math helpers
# ---------------------------------------------------------------------------

def _normalize(v: Point3D) -> Point3D:
    n = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return (v[0] / n, v[1] / n, v[2] / n) if n else (0.0, 0.0, 0.0)


def _tri_normal(tri: Triangle3D) -> Point3D:
    (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = tri
    u = (x2 - x1, y2 - y1, z2 - z1)
    v = (x3 - x1, y3 - y1, z3 - z1)
    return _normalize((u[1]*v[2] - u[2]*v[1],
                       u[2]*v[0] - u[0]*v[2],
                       u[0]*v[1] - u[1]*v[0]))


# ---------------------------------------------------------------------------
# Profile resampling  (ensures both root and tip have the same point count)
# ---------------------------------------------------------------------------

def _arc_lengths(pts: list[Point2D]) -> list[float]:
    s = [0.0]
    for i in range(1, len(pts)):
        dx = pts[i][0] - pts[i-1][0]
        dy = pts[i][1] - pts[i-1][1]
        s.append(s[-1] + math.hypot(dx, dy))
    return s


def _resample(pts: list[Point2D], n: int) -> list[Point2D]:
    """Resample pts to exactly n evenly-spaced points by arc length."""
    s = _arc_lengths(pts)
    total = s[-1]
    out: list[Point2D] = []
    j = 0
    for k in range(n):
        target = total * k / n
        while j < len(s) - 2 and s[j + 1] < target:
            j += 1
        seg = s[j + 1] - s[j]
        t = (target - s[j]) / seg if seg > 1e-12 else 0.0
        x = pts[j][0] + t * (pts[j+1][0] - pts[j][0])
        y = pts[j][1] + t * (pts[j+1][1] - pts[j][1])
        out.append((x, y))
    return out


# ---------------------------------------------------------------------------
# Wing lofting
# ---------------------------------------------------------------------------

def _lerp2d(a: list[Point2D], b: list[Point2D], t: float) -> list[Point2D]:
    """Linearly interpolate between two same-length profiles."""
    return [(a[i][0] * (1 - t) + b[i][0] * t,
             a[i][1] * (1 - t) + b[i][1] * t)
            for i in range(len(a))]


def _twist_profile(pts: list[Point2D], twist_deg: float) -> list[Point2D]:
    """Rotate profile about the quarter-chord point (x=0.25, y=0)."""
    cx, cy = 0.25, 0.0
    r = math.radians(twist_deg)
    cos_r, sin_r = math.cos(r), math.sin(r)
    result = []
    for x, y in pts:
        dx, dy = x - cx, y - cy
        result.append((cx + dx * cos_r - dy * sin_r,
                       cy + dx * sin_r + dy * cos_r))
    return result


def _place_section(
    unit_profile: list[Point2D],   # normalised chord (0..1)
    chord: float,                   # mm
    le_x: float,                    # leading-edge X offset (sweep)
    le_y: float,                    # leading-edge Y offset (dihedral)
    span_z: float,                  # spanwise position (mm)
) -> list[Point3D]:
    """Scale a unit profile and position it in 3-D space."""
    return [(le_x + p[0] * chord,
             le_y + p[1] * chord,
             span_z)
            for p in unit_profile]


def build_wing_mesh(
    root_pts: list[Point2D],
    tip_pts: list[Point2D],
    span_mm: float,
    root_chord_mm: float,
    taper: float,
    sweep_deg: float,
    dihedral_deg: float,
    washout_deg: float,
    sections: int,
    tip_style: str,
) -> list[Triangle3D]:

    n_pts = max(len(root_pts), len(tip_pts))
    # Ensure both profiles have the same point count
    root_rs = _resample(root_pts, n_pts)
    tip_rs  = _resample(tip_pts,  n_pts)

    # Ensure both are CCW (positive shoelace area)
    if _polygon_area(root_rs) < 0:
        root_rs.reverse()
    if _polygon_area(tip_rs) < 0:
        tip_rs.reverse()

    sweep_r    = math.radians(sweep_deg)
    dihedral_r = math.radians(dihedral_deg)

    # Build all spanwise cross-sections as lists of Point3D
    station_rings: list[list[Point3D]] = []
    for k in range(sections + 1):
        t = k / sections                          # 0 = root, 1 = tip
        chord  = root_chord_mm * (1.0 - t * (1.0 - taper))
        span_z = span_mm * t
        le_x   = span_z * math.tan(sweep_r)
        le_y   = span_z * math.tan(dihedral_r)
        twist  = washout_deg * t                  # linear washout

        # Interpolate profile shape, then twist
        unit = _lerp2d(root_rs, tip_rs, t)
        unit = _twist_profile(unit, twist)

        ring = _place_section(unit, chord, le_x, le_y, span_z)
        station_rings.append(ring)

    tris: list[Triangle3D] = []

    # ------------------------------------------------------------------
    # Root cap (k=0), outward normal toward -Z
    # ------------------------------------------------------------------
    ring0_2d = [(p[0], p[1]) for p in station_rings[0]]
    cap_tris = _ear_clip(ring0_2d)
    for i0, i1, i2 in cap_tris:
        p0, p1, p2 = station_rings[0][i0], station_rings[0][i1], station_rings[0][i2]
        tris.append((p0, p2, p1))   # reversed winding for -Z normal

    # ------------------------------------------------------------------
    # Side walls (lofted quads between adjacent stations)
    # ------------------------------------------------------------------
    for k in range(sections):
        ring_a = station_rings[k]
        ring_b = station_rings[k + 1]
        n = len(ring_a)
        for i in range(n):
            j = (i + 1) % n
            a0, a1 = ring_a[i], ring_a[j]
            b0, b1 = ring_b[i], ring_b[j]
            tris.append((a0, a1, b1))
            tris.append((a0, b1, b0))

    # ------------------------------------------------------------------
    # Tip cap
    # ------------------------------------------------------------------
    tip_ring = station_rings[-1]
    tip_2d   = [(p[0], p[1]) for p in tip_ring]

    if tip_style == "round":
        _add_rounded_tip(tris, tip_ring)
    else:
        # Flat cap
        cap_tris = _ear_clip(tip_2d)
        for i0, i1, i2 in cap_tris:
            p0, p1, p2 = tip_ring[i0], tip_ring[i1], tip_ring[i2]
            tris.append((p0, p1, p2))

    return tris


def _add_rounded_tip(tris: list[Triangle3D], ring: list[Point3D]) -> None:
    """
    Close the tip with a simple dome: fan-triangulate through a mid-latitude
    ring and an apex point, shrinking toward the ring centroid.
    """
    n = len(ring)
    cx = sum(p[0] for p in ring) / n
    cy = sum(p[1] for p in ring) / n
    cz = ring[0][2]

    rx = (max(p[0] for p in ring) - min(p[0] for p in ring)) * 0.5
    ry = (max(p[1] for p in ring) - min(p[1] for p in ring)) * 0.5
    dome_h = min(rx, ry) * 0.5   # low-profile dome

    # Intermediate latitude ring at 45°
    cos45 = math.cos(math.radians(45))
    sin45 = math.sin(math.radians(45))
    mid_ring: list[Point3D] = []
    for px, py, pz in ring:
        dx, dy = px - cx, py - cy
        mid_ring.append((cx + dx * cos45,
                         cy + dy * cos45,
                         cz + dome_h * sin45))

    apex: Point3D = (cx, cy, cz + dome_h)

    # Band: base ring → mid ring
    for i in range(n):
        j = (i + 1) % n
        a0, a1 = ring[i],     ring[j]
        b0, b1 = mid_ring[i], mid_ring[j]
        tris.append((a0, a1, b1))
        tris.append((a0, b1, b0))

    # Fan: mid ring → apex
    for i in range(n):
        j = (i + 1) % n
        tris.append((mid_ring[i], mid_ring[j], apex))


# ---------------------------------------------------------------------------
# STL output
# ---------------------------------------------------------------------------

def write_ascii_stl(path: Path, name: str, triangles: list[Triangle3D]) -> None:
    safe = "".join(c if c.isalnum() or c in "_-" else "_" for c in name) or "wing"
    with path.open("w", encoding="ascii") as f:
        f.write(f"solid {safe}\n")
        for tri in triangles:
            n = _tri_normal(tri)
            if n == (0.0, 0.0, 0.0):
                continue
            f.write(f"  facet normal {n[0]:.8e} {n[1]:.8e} {n[2]:.8e}\n")
            f.write("    outer loop\n")
            for vx, vy, vz in tri:
                f.write(f"      vertex {vx:.8e} {vy:.8e} {vz:.8e}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {safe}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a lofted 3-D wing STL from one or two Eppler profile.out files.\n"
            "All linear dimensions are in millimetres; angles in degrees."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Positional: root profile is required; tip profile is optional
    parser.add_argument(
        "root_profile",
        help="Path to profile.out for the root section (required).",
    )
    parser.add_argument(
        "tip_profile",
        nargs="?",
        default=None,
        help="Path to profile.out for the tip section (optional; defaults to same as root).",
    )

    # Output
    parser.add_argument(
        "-o", "--output",
        default="wing.stl",
        help="Output STL file path (default: wing.stl).",
    )

    # Wing geometry
    parser.add_argument("--span-mm",       type=float, default=300.0,
                        help="Half-span length in mm (default: 300).")
    parser.add_argument("--root-chord-mm", type=float, default=100.0,
                        help="Root chord length in mm (default: 100).")
    parser.add_argument("--taper",         type=float, default=0.6,
                        help="Taper ratio = tip_chord / root_chord (default: 0.6).")
    parser.add_argument("--sweep-deg",     type=float, default=15.0,
                        help="Leading-edge sweep angle in degrees (default: 15).")
    parser.add_argument("--dihedral-deg",  type=float, default=3.0,
                        help="Dihedral angle in degrees (default: 3).")
    parser.add_argument("--washout-deg",   type=float, default=-2.0,
                        help="Geometric twist at tip in degrees; negative = washout (default: -2).")
    parser.add_argument("--tip",           choices=["flat", "round"], default="flat",
                        help="Tip cap style: flat or round (default: flat).")
    parser.add_argument("--sections",      type=int, default=20,
                        help="Number of spanwise loft sections (default: 20).")

    args = parser.parse_args()

    # Validate
    if args.span_mm <= 0:
        parser.error("--span-mm must be > 0.")
    if args.root_chord_mm <= 0:
        parser.error("--root-chord-mm must be > 0.")
    if not (0.0 < args.taper <= 1.0):
        parser.error("--taper must be in (0, 1].")
    if args.sections < 2:
        parser.error("--sections must be >= 2.")

    root_path = Path(args.root_profile)
    tip_path  = Path(args.tip_profile) if args.tip_profile else root_path

    print(f"Root profile : {root_path}")
    print(f"Tip  profile : {tip_path}")

    root_pts, root_name = parse_profile_out(root_path)
    tip_pts,  tip_name  = parse_profile_out(tip_path)

    wing_name = (f"wing_{root_name}" if root_name == tip_name
                 else f"wing_{root_name}_to_{tip_name}")

    print(f"\nBuilding wing '{wing_name}':")
    print(f"  Span        : {args.span_mm:.1f} mm")
    print(f"  Root chord  : {args.root_chord_mm:.1f} mm")
    print(f"  Tip chord   : {args.root_chord_mm * args.taper:.1f} mm  (taper {args.taper})")
    print(f"  Sweep       : {args.sweep_deg:.1f}°")
    print(f"  Dihedral    : {args.dihedral_deg:.1f}°")
    print(f"  Washout     : {args.washout_deg:.1f}°")
    print(f"  Tip style   : {args.tip}")
    print(f"  Sections    : {args.sections}")

    tris = build_wing_mesh(
        root_pts      = root_pts,
        tip_pts       = tip_pts,
        span_mm       = args.span_mm,
        root_chord_mm = args.root_chord_mm,
        taper         = args.taper,
        sweep_deg     = args.sweep_deg,
        dihedral_deg  = args.dihedral_deg,
        washout_deg   = args.washout_deg,
        sections      = args.sections,
        tip_style     = args.tip,
    )

    out_path = Path(args.output)
    write_ascii_stl(out_path, wing_name, tris)

    print(f"\nWrote {out_path}  ({len(tris)} triangles)")


if __name__ == "__main__":
    main()
