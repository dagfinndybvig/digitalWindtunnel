#!/usr/bin/env python3
import argparse
import math
from pathlib import Path


Point2D = tuple[float, float]
Point3D = tuple[float, float, float]
Triangle3D = tuple[Point3D, Point3D, Point3D]


def parse_profile_out(path: Path) -> tuple[list[Point2D], str]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    airfoil_name = "airfoil"
    for line in lines:
        if "AIRFOIL" in line and "THICKNESS" in line:
            parts = line.split("AIRFOIL", 1)[1].split("%", 1)[0].strip().split()
            if parts:
                airfoil_name = parts[0]
            break

    header_idx = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("N") and "X" in stripped and "Y" in stripped:
            header_idx = i
            break

    if header_idx is None:
        raise ValueError("Could not find coordinate table header ('N X Y ...').")

    points: list[Point2D] = []
    for line in lines[header_idx + 1 :]:
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
        raise ValueError("Not enough coordinate points found in profile.out.")

    if almost_same(points[0], points[-1]):
        points.pop()

    return points, airfoil_name


def almost_same(a: Point2D, b: Point2D, tol: float = 1e-9) -> bool:
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol


def polygon_area(points: list[Point2D]) -> float:
    area2 = 0.0
    n = len(points)
    for i in range(n):
        x1, y1 = points[i]
        x2, y2 = points[(i + 1) % n]
        area2 += x1 * y2 - x2 * y1
    return 0.5 * area2


def cross2(a: Point2D, b: Point2D, c: Point2D) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def point_in_triangle(p: Point2D, a: Point2D, b: Point2D, c: Point2D) -> bool:
    eps = 1e-12
    c1 = cross2(a, b, p)
    c2 = cross2(b, c, p)
    c3 = cross2(c, a, p)
    has_neg = (c1 < -eps) or (c2 < -eps) or (c3 < -eps)
    has_pos = (c1 > eps) or (c2 > eps) or (c3 > eps)
    return not (has_neg and has_pos)


def triangulate_ear_clip(points_ccw: list[Point2D]) -> list[tuple[int, int, int]]:
    n = len(points_ccw)
    if n < 3:
        raise ValueError("Need at least 3 points to triangulate.")

    indices = list(range(n))
    triangles: list[tuple[int, int, int]] = []
    eps = 1e-12
    max_loops = n * n
    loops = 0

    while len(indices) > 3:
        ear_found = False
        m = len(indices)
        for i in range(m):
            i_prev = indices[(i - 1) % m]
            i_curr = indices[i]
            i_next = indices[(i + 1) % m]

            a = points_ccw[i_prev]
            b = points_ccw[i_curr]
            c = points_ccw[i_next]

            if cross2(a, b, c) <= eps:
                continue

            contains_point = False
            for j in indices:
                if j in (i_prev, i_curr, i_next):
                    continue
                if point_in_triangle(points_ccw[j], a, b, c):
                    contains_point = True
                    break

            if contains_point:
                continue

            triangles.append((i_prev, i_curr, i_next))
            del indices[i]
            ear_found = True
            break

        loops += 1
        if not ear_found or loops > max_loops:
            raise ValueError("Triangulation failed; polygon may be invalid.")

    triangles.append((indices[0], indices[1], indices[2]))
    return triangles


def normalize(v: Point3D) -> Point3D:
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n == 0.0:
        return (0.0, 0.0, 0.0)
    return (v[0] / n, v[1] / n, v[2] / n)


def triangle_normal(tri: Triangle3D) -> Point3D:
    (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = tri
    ux, uy, uz = x2 - x1, y2 - y1, z2 - z1
    vx, vy, vz = x3 - x1, y3 - y1, z3 - z1
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx
    return normalize((nx, ny, nz))


def build_extruded_mesh(points: list[Point2D], span_mm: float) -> list[Triangle3D]:
    if span_mm <= 0:
        raise ValueError("span-mm must be > 0.")

    points_ccw = points[:]
    if polygon_area(points_ccw) < 0:
        points_ccw.reverse()

    cap_tris = triangulate_ear_clip(points_ccw)
    tris: list[Triangle3D] = []

    # Bottom cap (z=0), outward normal toward -z.
    for i0, i1, i2 in cap_tris:
        p0 = (points_ccw[i0][0], points_ccw[i0][1], 0.0)
        p1 = (points_ccw[i1][0], points_ccw[i1][1], 0.0)
        p2 = (points_ccw[i2][0], points_ccw[i2][1], 0.0)
        tris.append((p0, p2, p1))

    # Top cap (z=span), outward normal toward +z.
    for i0, i1, i2 in cap_tris:
        p0 = (points_ccw[i0][0], points_ccw[i0][1], span_mm)
        p1 = (points_ccw[i1][0], points_ccw[i1][1], span_mm)
        p2 = (points_ccw[i2][0], points_ccw[i2][1], span_mm)
        tris.append((p0, p1, p2))

    # Side walls.
    n = len(points_ccw)
    for i in range(n):
        j = (i + 1) % n
        x0, y0 = points_ccw[i]
        x1, y1 = points_ccw[j]
        b0 = (x0, y0, 0.0)
        b1 = (x1, y1, 0.0)
        t1 = (x1, y1, span_mm)
        t0 = (x0, y0, span_mm)
        tris.append((b0, b1, t1))
        tris.append((b0, t1, t0))

    return tris


def write_ascii_stl(path: Path, name: str, triangles: list[Triangle3D]) -> None:
    safe_name = "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in name) or "airfoil"
    with path.open("w", encoding="ascii", errors="strict") as f:
        f.write(f"solid {safe_name}\n")
        for tri in triangles:
            n = triangle_normal(tri)
            if n == (0.0, 0.0, 0.0):
                continue
            f.write(f"  facet normal {n[0]:.8e} {n[1]:.8e} {n[2]:.8e}\n")
            f.write("    outer loop\n")
            for vx, vy, vz in tri:
                f.write(f"      vertex {vx:.8e} {vy:.8e} {vz:.8e}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {safe_name}\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert airfoil coordinates in profile.out into an extruded ASCII STL."
    )
    parser.add_argument("input", nargs="?", default="profile.out", help="Path to profile.out")
    parser.add_argument("output", nargs="?", default="airfoil_section.stl", help="Output STL path")
    parser.add_argument("--chord-mm", type=float, default=100.0, help="Chord length scaling in millimeters")
    parser.add_argument("--span-mm", type=float, default=20.0, help="Extrusion span in millimeters")
    args = parser.parse_args()

    if args.chord_mm <= 0:
        raise ValueError("chord-mm must be > 0.")

    input_path = Path(args.input)
    output_path = Path(args.output)

    points, airfoil_name = parse_profile_out(input_path)
    scaled = [(x * args.chord_mm, y * args.chord_mm) for x, y in points]
    triangles = build_extruded_mesh(scaled, args.span_mm)
    write_ascii_stl(output_path, f"airfoil_{airfoil_name}", triangles)

    print(f"Wrote {output_path}")
    print(f"Airfoil: {airfoil_name}")
    print(f"Points: {len(points)}")
    print(f"Chord: {args.chord_mm:.2f} mm")
    print(f"Span:  {args.span_mm:.2f} mm")
    print(f"Triangles: {len(triangles)}")


if __name__ == "__main__":
    main()
