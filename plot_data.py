#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys


def parse_profile_out(filename='profile.out'):
    """Parse profile.out to extract airfoil geometry and velocity distributions."""
    x, y = [], []
    velocities = {}
    airfoil_name = "Unknown"
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find the airfoil name (e.g., " AIRFOIL 1098           18.97% THICKNESS")
    for line in lines:
        if 'AIRFOIL' in line and 'THICKNESS' in line:
            # Extract the airfoil ID (between "AIRFOIL" and "%")
            parts = line.split('AIRFOIL')
            if len(parts) > 1:
                name_part = parts[1].split('%')[0].strip()
                # Get just the ID (first word)
                airfoil_name = name_part.split()[0] if name_part else "Unknown"
            break
    
    # Find the header line that contains the angle values
    # Format: "  N     X        Y      2.00   8.00  10.00  12.00"
    header_idx = None
    angle_values = []
    
    for i, line in enumerate(lines):
        if line.strip().startswith('N') and 'X' in line and 'Y' in line:
            header_idx = i
            # Parse angle values from header
            parts = line.strip().split()
            # angles start at index 3 (after N, X, Y)
            for part in parts[3:]:
                try:
                    angle_values.append(float(part))
                except ValueError:
                    continue
            break
    
    if header_idx is None:
        raise ValueError("Could not find velocity distribution header in output")
    
    # Initialize velocity arrays for each angle
    for angle in angle_values:
        velocities[angle] = []
    
    # Read the velocity data table
    for line in lines[header_idx + 1:]:
        parts = line.strip().split()
        if len(parts) < 3:
            break
        try:
            xval = float(parts[1])
            yval = float(parts[2])
        except (ValueError, IndexError):
            break
        
        x.append(xval)
        y.append(yval)
        
        # Parse velocity values
        for j, angle in enumerate(angle_values):
            idx = 3 + j
            if idx < len(parts):
                try:
                    velocities[angle].append(float(parts[idx]))
                except (ValueError, IndexError):
                    velocities[angle].append(np.nan)
            else:
                velocities[angle].append(np.nan)
    
    return np.array(x), np.array(y), velocities, angle_values, airfoil_name


def plot_combined(x, y, velocities, angles, airfoil_name, output_file='combined_panel.png', variant_label=None):
    """Create a combined plot with airfoil profile and velocity distributions."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    
    # Add overall title with airfoil name
    title = f'Eppler Airfoil Analysis: {airfoil_name}'
    if variant_label:
        title += f' ({variant_label})'
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
    
    # --- Top plot: Velocity distribution ---
    colors = plt.cm.viridis(np.linspace(0, 1, len(angles)))
    
    for i, angle in enumerate(sorted(angles)):
        v = velocities[angle]
        ax1.plot(x, v, color=colors[i], label=f'α = {angle}°', linewidth=1.5)
    
    ax1.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='V/U∞ = 1')
    ax1.set_xlabel('x/c', fontsize=12)
    ax1.set_ylabel('V/U∞', fontsize=12)
    ax1.set_title('Velocity Distribution', fontsize=14)
    ax1.legend(loc='best', fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    
    # --- Bottom plot: Airfoil profile ---
    # Separate upper and lower surfaces
    # Find the leading edge (minimum x, typically around x=0)
    idx_min_x = np.argmin(np.abs(x - 0.0))
    
    # Upper surface: from leading edge to trailing edge (first half)
    # Lower surface: from trailing edge back to leading edge (second half)
    # The airfoil goes: trailing edge -> upper surface -> leading edge -> lower surface -> trailing edge
    
    # Find the index where y transitions from positive to negative (roughly at leading edge)
    # Actually, looking at the data, x goes from 1.0 down to ~0 and back to 1.0
    # Let's find the LE by looking at where x is minimum
    
    # First half is upper surface (descending x), second half is lower surface (ascending x)
    n_points = len(x)
    mid_idx = n_points // 2
    
    x_upper = x[:mid_idx+1]
    y_upper = y[:mid_idx+1]
    x_lower = x[mid_idx:]
    y_lower = y[mid_idx:]
    
    ax2.fill(x_upper, y_upper, color='gray', alpha=0.3)
    ax2.fill(x_lower, y_lower, color='gray', alpha=0.3)
    ax2.plot(x, y, 'k-', linewidth=1.5)
    ax2.plot(x, y, 'k.', markersize=2, alpha=0.5)  # Show points
    
    ax2.set_xlabel('x/c', fontsize=12)
    ax2.set_ylabel('y/c', fontsize=12)
    ax2.set_title('Airfoil Profile', fontsize=14)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(-0.05, 1.05)
    
    # Calculate and display thickness
    thickness = np.max(y_upper) - np.min(y_lower)
    max_thickness_x = x_upper[np.argmax(y_upper)]
    ax2.annotate(f't/c = {thickness*100:.1f}%\n@ x/c = {max_thickness_x:.2f}',
                 xy=(max_thickness_x, np.max(y_upper)),
                 xytext=(0.7, 0.02),
                 fontsize=10,
                 arrowprops=dict(arrowstyle='->', color='red', alpha=0.7),
                 )
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {output_file}")
    
    # Also display info
    print(f"\nAirfoil Analysis Summary:")
    print(f"  Airfoil: {airfoil_name}")
    print(f"  Number of points: {len(x)}")
    print(f"  Thickness: {thickness*100:.2f}%")
    print(f"  Angles analyzed: {sorted(angles)}")


def main():
    # Allow filename as command line argument
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = 'profile.out'
    variant_label = sys.argv[2] if len(sys.argv) > 2 else None
    output_file = sys.argv[3] if len(sys.argv) > 3 else 'combined_panel.png'
    
    try:
        x, y, velocities, angles, airfoil_name = parse_profile_out(filename)
        plot_combined(x, y, velocities, angles, airfoil_name, output_file=output_file, variant_label=variant_label)
    except FileNotFoundError:
        print(f"Error: Could not find {filename}")
        print("Run the profile program first: ./profile < input.dat")
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing {filename}: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
