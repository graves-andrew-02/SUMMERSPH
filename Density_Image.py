""" 
This is made to be run on a google colab session, as such the file path may need changing for running on different jupyter notebook style python sessions.

As for the content of this file, it takes the smoothing length and a file, then produces an image of the density of the image (and gives just the a perpective along the -z axis)
"""

import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from tqdm import tqdm
from numba import jit
import pandas as pd

# Load data
# File path
number = 275
file_path = f'/content/save{number}.txt'


# Print first few lines of the file
print(f"First 5 lines of {file_path}:")
with open(file_path, 'r') as file:
    for i in range(5):
        try:
            print(next(file).strip())
        except StopIteration:
            break
print("-" * 20)


# Lists to store each column
x, y, z = [], [], []
vx, vy, vz = [], [], []
energy, density, mass = [], [], []
alpha = []

# Read and parse the file
with open(file_path, 'r') as file:
    next(file)  # Skip the header line
    for line in file:
        values = line.split()
        if len(values) == 9:  # Ensure line has all 9 fields
            x.append(float(values[0]))
            y.append(float(values[1]))
            z.append(float(values[2]))
            vx.append(float(values[3]))
            vy.append(float(values[4]))
            vz.append(float(values[5]))
            energy.append(float(values[6]))
            mass.append(float(values[7]))
            alpha.append(float(values[8]))

# Convert lists to numpy arrays (optional, for numerical operations)
x = np.array(x)
y = np.array(y)
z = np.array(z)
vx = np.array(vx)
vy = np.array(vy)
vz = np.array(vz)
energy = np.array(energy)
mass = np.array(mass)
alpha = np.array(alpha)

size_to_image = 100 # Adjusted mask range
mask = (x< size_to_image) & (x>-size_to_image) & (y< size_to_image) & (y >-size_to_image) & (z < size_to_image) & (z >-size_to_image)
x = x[mask]
y = y[mask]
z = z[mask]
vx = vx[mask]
vy = vy[mask]
vz = vz[mask]
energy = energy[mask]
mass = mass[mask]
alpha = alpha[mask]

x_sun = x[-1]
y_sun = y[-1]
z_sun = z[-1]

x = x[:-1]
y = y[:-1]
z = z[:-1]
vx = vx[:-1]
vy = vy[:-1]
vz = vz[:-1]
energy = energy[:-1]
mass = mass[:-1]
alpha = alpha[:-1]

# Cubic spline kernel
@jit
def cubic_spline_kernel(r, h):
    q = r / h
    sigma = 1 / (np.pi * h**3)
    w = np.zeros_like(r)
    mask1 = q <= 1
    mask2 = (q > 1) & (q <= 2)
    w[mask1] = sigma * (1 - 1.5*q[mask1]**2 + 0.75*q[mask1]**3)
    w[mask2] = sigma * 0.25 * (2 - q[mask2])**3
    return w

print(f"Number of data points after masking: {len(x)}")

# Constants
h = 1.25
kernel_radius =  2*h

# Bounds
xmin, xmax = x.min(), x.max()
ymin, ymax = y.min(), y.max()
zmin, zmax = z.min(), z.max()

# Grid resolution
grid_resolution = 120
xi = np.linspace(xmin, xmax, grid_resolution)
yi = np.linspace(ymin, ymax, grid_resolution)
zi = np.linspace(zmin, zmax, grid_resolution)
dx = xi[1] - xi[0]
dz = zi[1] - zi[0]  # for integration scaling

# Create 3D grid
X, Y, Z = np.meshgrid(xi, yi, zi, indexing='ij')
grid_density = np.zeros_like(X)

# KD-Tree for efficient neighbor search
positions = np.vstack((x, y, z)).T
tree = cKDTree(positions)

# Compute density grid
for idx in tqdm(range(grid_density.size)):
    i, j, k = np.unravel_index(idx, grid_density.shape)
    gx, gy, gz = X[i, j, k], Y[i, j, k], Z[i, j, k]
    neighbors_idx = tree.query_ball_point([gx, gy, gz], r=kernel_radius)
    if not neighbors_idx:
        continue
    r_vecs = positions[neighbors_idx] - np.array([gx, gy, gz])
    r_mags = np.linalg.norm(r_vecs, axis=1)
    w_vals = cubic_spline_kernel(r_mags, h)
    grid_density[i, j, k] = np.sum(mass[neighbors_idx] * w_vals)


# Integrate along z-axis (axis=2), multiplied by dz to scale correctly
projected_density = np.sum(grid_density, axis=2)# * dz

# Plot the 2D integrated density
plt.imshow(projected_density.T, origin='lower',
           extent=[xmin, xmax, ymin, ymax],
           cmap='inferno')
plt.colorbar(label='Integrated Density')
plt.title('Integrated SPH Density (Projection along Z)')
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x_sun,y_sun,'.', color = 'red', markersize = 1.5)
#plt.tight_layout()
plt.show()
