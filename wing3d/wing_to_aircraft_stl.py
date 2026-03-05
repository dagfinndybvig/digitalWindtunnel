#!/usr/bin/env python3
"""
wing_to_aircraft_stl.py

Build a quick aircraft preview STL by:
1) loading one half-wing STL,
2) mirroring it into left/right wings,
3) adding a proportional cigar-shaped fuselage.
"""

import argparse
import math
import struct
from pathlib import Path


Point3D = tuple[float, float, float]
Triangle3D = tuple[Point3D, Point3D, Point3D]


def parse_ascii_stl(text: str) -> list[Triangle3D]:
    triangles: list[Triangle3D] = []
    vertices: list[Point3D] = []
    for line in text.splitlines():
        s = line.strip()
        if not s.lower().startswith("vertex"):
            continue
        parts = s.split()
        if len(parts) < 4:
            continue
        try:
            v = (float(parts[1]), float(parts[2]), float(parts[3]))
        except ValueError:
            continue
        vertices.append(v)
        if len(vertices) == 3:
            triangles.append((vertices[0], vertices[1], vertices[2]))
            vertices.clear()
    return triangles


def parse_binary_stl(data: bytes) -> list[Triangle3D]:
    if len(data) < 84:
        raise ValueError("File is too small to be a valid binary STL.")
    tri_count = struct.unpack("<I", data[80:84])[0]
    expected_size = 84 + tri_count * 50
    if expected_size > len(data):
        raise ValueError("Binary STL appears truncated.")

    triangles: list[Triangle3D] = []
    offset = 84
    for _ in range(tri_count):
        block = data[offset : offset + 50]
        values = struct.unpack("<12f", block[:48])  # normal(3) + vertices(9)
        p1 = (values[3], values[4], values[5])
        p2 = (values[6], values[7], values[8])
        p3 = (values[9], values[10], values[11])
        triangles.append((p1, p2, p3))
        offset += 50
    return triangles


def read_stl(path: Path) -> list[Triangle3D]:
    data = path.read_bytes()
    try:
        text = data.decode("utf-8")
        triangles = parse_ascii_stl(text)
        if triangles:
            return triangles
    except UnicodeDecodeError:
        pass
    triangles = parse_binary_stl(data)
    if not triangles:
        raise ValueError(f"{path}: no triangles found.")
    return triangles


def mesh_bounds(tris: list[Triangle3D]) -> tuple[float, float, float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for tri in tris:
        for x, y, z in tri:
            xs.append(x)
            ys.append(y)
            zs.append(z)
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


def triangle_normal(tri: Triangle3D) -> Point3D:
    (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = tri
    ux, uy, uz = x2 - x1, y2 - y1, z2 - z1
    vx, vy, vz = x3 - x1, y3 - y1, z3 - z1
    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx
    n = math.sqrt(nx * nx + ny * ny + nz * nz)
    if n == 0.0:
        return (0.0, 0.0, 0.0)
    return (nx / n, ny / n, nz / n)


def orient_outward(tri: Triangle3D, center: Point3D) -> Triangle3D:
    n = triangle_normal(tri)
    if n == (0.0, 0.0, 0.0):
        return tri
    cx = (tri[0][0] + tri[1][0] + tri[2][0]) / 3.0
    cy = (tri[0][1] + tri[1][1] + tri[2][1]) / 3.0
    cz = (tri[0][2] + tri[1][2] + tri[2][2]) / 3.0
    vx = cx - center[0]
    vy = cy - center[1]
    vz = cz - center[2]
    if n[0] * vx + n[1] * vy + n[2] * vz < 0.0:
        return (tri[0], tri[2], tri[1])
    return tri


def translate_mesh(tris: list[Triangle3D], dx: float, dy: float, dz: float) -> list[Triangle3D]:
    return [
        (
            (a[0] + dx, a[1] + dy, a[2] + dz),
            (b[0] + dx, b[1] + dy, b[2] + dz),
            (c[0] + dx, c[1] + dy, c[2] + dz),
        )
        for a, b, c in tris
    ]


def mirror_mesh_across_centerline(tris: list[Triangle3D]) -> list[Triangle3D]:
    mirrored: list[Triangle3D] = []
    for a, b, c in tris:
        ma = (a[0], a[1], -a[2])
        mb = (b[0], b[1], -b[2])
        mc = (c[0], c[1], -c[2])
        mirrored.append((ma, mc, mb))  # reverse winding after reflection
    return mirrored


def root_section_stats(tris: list[Triangle3D], root_z: float, fallback_bounds: tuple[float, float, float, float]) -> tuple[float, float]:
    x_min_fallback, x_max_fallback, y_min_fallback, y_max_fallback = fallback_bounds
    xs: list[float] = []
    ys: list[float] = []
    for tri in tris:
        for x, y, z in tri:
            if abs(z - root_z) <= 1e-6:
                xs.append(x)
                ys.append(y)
    if len(xs) < 3:
        return (0.5 * (x_min_fallback + x_max_fallback), 0.5 * (y_min_fallback + y_max_fallback))
    quarter_chord_x = min(xs) + 0.25 * (max(xs) - min(xs))
    y_center = sum(ys) / len(ys)
    return quarter_chord_x, y_center


def build_cigar_fuselage(
    center: Point3D,
    length: float,
    radius: float,
    axial_segments: int,
    radial_segments: int,
) -> list[Triangle3D]:
    if length <= 0 or radius <= 0:
        raise ValueError("Fuselage length and radius must be > 0.")
    if axial_segments < 4:
        raise ValueError("axial_segments must be >= 4.")
    if radial_segments < 8:
        raise ValueError("radial_segments must be >= 8.")

    cx, cy, cz = center
    x_positions = [(-0.5 + i / axial_segments) * length for i in range(axial_segments + 1)]

    rings: list[list[Point3D]] = []
    for x_local in x_positions:
        profile = 1.0 - (2.0 * x_local / length) ** 2
        local_radius = radius * math.sqrt(max(0.0, profile))
        x = cx + x_local
        if local_radius <= 1e-10:
            rings.append([(x, cy, cz)])
            continue
        ring: list[Point3D] = []
        for j in range(radial_segments):
            theta = 2.0 * math.pi * j / radial_segments
            ring.append((x, cy + local_radius * math.cos(theta), cz + local_radius * math.sin(theta)))
        rings.append(ring)

    tris: list[Triangle3D] = []
    for i in range(len(rings) - 1):
        a = rings[i]
        b = rings[i + 1]
        if len(a) == 1:
            pole = a[0]
            for j in range(radial_segments):
                k = (j + 1) % radial_segments
                tri = (pole, b[k], b[j])
                tris.append(orient_outward(tri, center))
            continue
        if len(b) == 1:
            pole = b[0]
            for j in range(radial_segments):
                k = (j + 1) % radial_segments
                tri = (a[j], a[k], pole)
                tris.append(orient_outward(tri, center))
            continue
        for j in range(radial_segments):
            k = (j + 1) % radial_segments
            tri1 = (a[j], b[j], b[k])
            tri2 = (a[j], b[k], a[k])
            tris.append(orient_outward(tri1, center))
            tris.append(orient_outward(tri2, center))
    return tris


def write_ascii_stl(path: Path, name: str, triangles: list[Triangle3D]) -> None:
    safe_name = "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in name) or "aircraft"
    with path.open("w", encoding="ascii") as f:
        f.write(f"solid {safe_name}\n")
        for tri in triangles:
            n = triangle_normal(tri)
            if n == (0.0, 0.0, 0.0):
                continue
            f.write(f"  facet normal {n[0]:.8e} {n[1]:.8e} {n[2]:.8e}\n")
            f.write("    outer loop\n")
            for x, y, z in tri:
                f.write(f"      vertex {x:.8e} {y:.8e} {z:.8e}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {safe_name}\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mirror a half-wing STL and combine it with a proportional cigar-shaped fuselage."
    )
    parser.add_argument("wing_stl", help="Input half-wing STL path (ASCII or binary).")
    parser.add_argument("-o", "--output", default="aircraft_preview.stl", help="Output STL path.")
    parser.add_argument(
        "--fuselage-length-ratio",
        type=float,
        default=0.65,
        help="Fuselage length as ratio of full wingspan (default: 0.65).",
    )
    parser.add_argument(
        "--fuselage-radius-ratio",
        type=float,
        default=0.12,
        help="Fuselage radius as ratio of wing chord extent (default: 0.12).",
    )
    parser.add_argument(
        "--wing-offset-mm",
        type=float,
        default=None,
        help="Centerline-to-root offset for each wing (default: 0.9 * fuselage radius).",
    )
    parser.add_argument(
        "--axial-segments",
        type=int,
        default=28,
        help="Fuselage longitudinal segments (default: 28).",
    )
    parser.add_argument(
        "--radial-segments",
        type=int,
        default=36,
        help="Fuselage radial segments (default: 36).",
    )
    args = parser.parse_args()

    if args.fuselage_length_ratio <= 0:
        parser.error("--fuselage-length-ratio must be > 0.")
    if args.fuselage_radius_ratio <= 0:
        parser.error("--fuselage-radius-ratio must be > 0.")

    wing_path = Path(args.wing_stl)
    out_path = Path(args.output)

    base_wing = read_stl(wing_path)
    x_min, x_max, y_min, y_max, z_min, z_max = mesh_bounds(base_wing)
    if z_max <= z_min:
        raise ValueError("Wing STL has invalid span extent.")

    half_span = max(abs(z_min), abs(z_max))
    full_span = 2.0 * half_span
    chord_extent = max(1e-6, x_max - x_min)
    thickness_extent = max(1e-6, y_max - y_min)

    fuselage_length = full_span * args.fuselage_length_ratio
    fuselage_radius = max(chord_extent * args.fuselage_radius_ratio, 0.55 * thickness_extent)
    wing_offset = args.wing_offset_mm if args.wing_offset_mm is not None else fuselage_radius * 0.9

    right_wing = translate_mesh(base_wing, 0.0, 0.0, wing_offset)
    left_wing = mirror_mesh_across_centerline(right_wing)

    attach_x, attach_y = root_section_stats(
        base_wing,
        z_min,
        (x_min, x_max, y_min, y_max),
    )
    fuselage_center = (attach_x, attach_y, 0.0)
    fuselage = build_cigar_fuselage(
        center=fuselage_center,
        length=fuselage_length,
        radius=fuselage_radius,
        axial_segments=args.axial_segments,
        radial_segments=args.radial_segments,
    )

    combined = right_wing + left_wing + fuselage
    write_ascii_stl(out_path, f"{wing_path.stem}_aircraft_preview", combined)

    print(f"Input wing triangles : {len(base_wing)}")
    print(f"Mirrored pair        : {len(right_wing) + len(left_wing)}")
    print(f"Fuselage triangles   : {len(fuselage)}")
    print(f"Fuselage length      : {fuselage_length:.2f} mm")
    print(f"Fuselage radius      : {fuselage_radius:.2f} mm")
    print(f"Wing root offset     : {wing_offset:.2f} mm")
    print(f"Wrote {out_path}  ({len(combined)} triangles)")


if __name__ == "__main__":
    main()
