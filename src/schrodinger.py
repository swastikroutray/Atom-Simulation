import numpy as np
import json
import os
from math import factorial, sqrt
from scipy.special import genlaguerre, sph_harm

a0 = 1.0

# Radial wave function R_{n,l}(r)
def R_nl(n, l, r):
    rho = 2.0 * r / (n * a0)
    norm = sqrt((2.0 / (n * a0))**3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l)))
    # scipy.special.genlaguerre(degree, alpha)
    L = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    return norm * np.exp(-rho / 2.0) * (rho**l) * L

# Full complex orbital psi(r, theta, phi)
def psi(n, l, m, r, theta, phi):
    # sph_harm(m, l, phi, theta) in standard scipy (azimuth=phi, colatitude=theta)
    Y = sph_harm(m, l, phi, theta)
    return R_nl(n, l, r) * Y

def sample_orbital_points(n, l, m, target_count=50000, batch_size=200000):
    """
    Vectorized rejection sampling with the proper 3D volume element dV = r^2 sin(theta) dr dtheta dphi.
    Target spatial density: P(r, theta, phi) = |psi|^2 * r^2 * sin(theta)
    """
    accepted_points = []
    
    # 1. Estimate maximum envelope value across a grid
    r_max = max(30.0, float(2.5 * n**2 * a0))
    r_grid = np.linspace(0.01, r_max, 300)
    th_grid = np.linspace(0, np.pi, 150)
    ph_grid = np.linspace(0, 2 * np.pi, 50)
    
    R_mesh, TH_mesh, PH_mesh = np.meshgrid(r_grid, th_grid, ph_grid, indexing='ij')
    psi_mesh = psi(n, l, m, R_mesh, TH_mesh, PH_mesh)
    # Volume element weight included:
    target_density = (np.abs(psi_mesh)**2) * (R_mesh**2) * np.sin(TH_mesh)
    max_density = np.max(target_density) * 1.15  # 15% safety buffer

    # 2. Fast batched rejection sampling
    while len(accepted_points) < target_count:
        # Uniform proposal inside the bounding sphere/box
        r_prop = np.random.uniform(0.0, r_max, batch_size)
        th_prop = np.random.uniform(0.0, np.pi, batch_size)
        ph_prop = np.random.uniform(0.0, 2.0 * np.pi, batch_size)

        psi_vals = psi(n, l, m, r_prop, th_prop, ph_prop)
        prop_density = (np.abs(psi_vals)**2) * (r_prop**2) * np.sin(th_prop)

        # Acceptance criterion
        rand_thresh = np.random.uniform(0.0, max_density, batch_size)
        mask = rand_thresh < prop_density

        r_acc = r_prop[mask]
        th_acc = th_prop[mask]
        ph_acc = ph_prop[mask]
        psi_acc = psi_vals[mask]

        # Convert to Cartesian (x, y, z)
        x = r_acc * np.sin(th_acc) * np.cos(ph_acc)
        y = r_acc * np.sin(th_acc) * np.sin(ph_acc)
        z = r_acc * np.cos(th_acc)

        for i in range(len(x)):
            if len(accepted_points) >= target_count:
                break
            accepted_points.append({
                "x": float(x[i]),
                "y": float(y[i]),
                "z": float(z[i]),
                "psi_re": float(psi_acc[i].real),
                "psi_im": float(psi_acc[i].imag),
                "phase": float(np.angle(psi_acc[i]))
            })

    return accepted_points

def main():
    print("\n=== Hydrogen Orbital Generator (Vectorized & Corrected) ===\n")
    n = int(input("Enter n (1..7): "))
    l = int(input(f"Enter l (0..{n-1}): "))
    m = int(input(f"Enter m (-{l}..{l}): "))
    target_samples = int(input("How many particle samples? (ex: 50000): "))

    print("\nGenerating probability cloud...")
    pts = sample_orbital_points(n, l, m, target_count=target_samples)

    os.makedirs("orbitals", exist_ok=True)
    filename = f"orbitals/orbital_n{n}_l{l}_m{m}.json"
    with open(filename, "w") as f:
        json.dump({
            "n": n,
            "l": l,
            "m": m,
            "sample_count": len(pts),
            "points": pts
        }, f)

    print(f"\nSaved {len(pts)} samples to: {filename}")

if __name__ == "__main__":
    main()