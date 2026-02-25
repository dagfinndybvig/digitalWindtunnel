import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

def parse_profile_out(filename):
    x, y = [], []
    velocities = {2.0: [], 8.0: [], 10.0: [], 12.0: []}
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Find the start of the velocity distribution table
    for i, line in enumerate(lines):
        if line.strip().startswith('N') and 'X' in line and 'Y' in line:
            start = i + 1
            break
    # Read the table
    for line in lines[start:]:
        parts = line.strip().split()
        if len(parts) < 7:
            break
        try:
            xval = float(parts[1])
            yval = float(parts[2])
            v2 = float(parts[3])
            v8 = float(parts[4])
            v10 = float(parts[5])
            v12 = float(parts[6])
        except ValueError:
            continue  # skip lines that can't be parsed
        x.append(xval)
        y.append(yval)
        velocities[2.0].append(v2)
        velocities[8.0].append(v8)
        velocities[10.0].append(v10)
        velocities[12.0].append(v12)
    return np.array(x), np.array(y), velocities

profile_x, profile_y, velocities = parse_profile_out('profile.out')

plt.figure(figsize=(8, 6))

# Plot velocity distributions
for angle, v in velocities.items():
    plt.plot(profile_x, v, label=f'Velocity {angle}°')

# Overlay airfoil profile, scaled to fit velocity axis
profile_y_scaled = 0.4 + (profile_y - np.min(profile_y)) / (np.max(profile_y) - np.min(profile_y)) * 0.3
plt.plot(profile_x, profile_y_scaled, 'k-', label='Airfoil profile (scaled)')

plt.xlabel('x/c')
plt.ylabel('V/U_inf or scaled y/c')
plt.title('Airfoil Profile and Velocity Distributions')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('combined_panel.png')
plt.savefig('combined_panel.png', dpi=300)
plt.show()
