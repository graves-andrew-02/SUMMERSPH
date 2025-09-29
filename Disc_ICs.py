# Constants (in appropriate units)
G = 4 * np.pi**2 #39.47841760435743  # AU^3 / (Msun * yr^2)
M_star = 5.0           # Msun
N_particles = 12000
R_max = 100; R_min = 0   # AU
mass_particle = 5 / N_particles   # placeholder mass per particle

# Generate radial positions with uniform surface density
n = 0
positions = []
while n < N_particles:
  test_pos = 2*R_max*(np.random.rand(3)-0.5)
  #rejection test
  if np.sum(test_pos**2)**0.5 <= R_max and np.sum(test_pos**2)**0.5 >= R_min:
    positions.append(test_pos)
    n += 1
x = np.array([pos[0] for pos in positions])
y = np.array([pos[1] for pos in positions])
z = np.array([pos[2] for pos in positions])

v_circ = 1
vx = -v_circ * y/r
vy = v_circ * x/r
vz = np.zeros(N_particles)

energy = 0.25 * np.ones(N_particles)
mass = np.full(N_particles, mass_particle)

# Create DataFrame and save
df = pd.DataFrame({
    'x': x, 'y': y, 'z': z,
    'vx': vx, 'vy': vy, 'vz': vz,
    'energy': energy,
    #'sph_mass': sph_mass,
    'mass': mass
})

# Save to file
df.to_csv("Collapse.txt", sep=' ', index=False, float_format='%.15e')
print(np.sqrt(R_max**3/N_particles))
