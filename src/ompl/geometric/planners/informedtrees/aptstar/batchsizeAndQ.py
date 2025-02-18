import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters for M and decay factor range
M_vals = np.linspace(1, 199, 100)      # Batch size M values
decay_factors = np.linspace(0, 1, 100) # Decay factor range (ratio)
lambda_ = 500

# Define functions for each decay function
def q_brachistochrone(M, decay_factor):
    theta = np.pi / 2 * decay_factor
    brachistochrone_factor = np.sqrt(np.sin(theta))
    return M * brachistochrone_factor

def q_iterative(M, decay_factor):
    return np.where(
        (0.1 <= decay_factor) & (decay_factor < 0.3), M * 0.4,
        np.where(
            (0.3 <= decay_factor) & (decay_factor < 0.5), M * 0.5,
            np.where(
                (0.5 <= decay_factor) & (decay_factor < 0.8), M * 0.65,
                np.where(
                    (0.8 <= decay_factor) & (decay_factor <= 1.0), M * 1.0,
                    M * 0.1
                )
            )
        )
    )

def q_linear(M, decay_factor):
    return M * decay_factor

def q_parabolic(M, decay_factor):
    return M * np.sqrt(decay_factor)

def q_logarithmic(M, decay_factor):
    return M * np.log(1 + lambda_ * decay_factor) / np.log(1 + lambda_)

# Generate mesh grid for M and decay_factor
M_mesh, decay_factor_mesh = np.meshgrid(M_vals, decay_factors)

# Calculate q values for each decay function
q_brachistochrone_vals = q_brachistochrone(M_mesh, decay_factor_mesh)
q_iterative_vals = q_iterative(M_mesh, decay_factor_mesh)
q_linear_vals = q_linear(M_mesh, decay_factor_mesh)
q_parabolic_vals = q_parabolic(M_mesh, decay_factor_mesh)
q_logarithmic_vals = q_logarithmic(M_mesh, decay_factor_mesh)

# Set up 3D plot
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot each decay function with terrain-like color maps
ax.plot_surface(M_mesh, decay_factor_mesh, q_brachistochrone_vals, cmap='Blues', alpha=0.6, edgecolor='none')
ax.plot_surface(M_mesh, decay_factor_mesh, q_iterative_vals, cmap='Greens', alpha=0.6, edgecolor='none')
ax.plot_surface(M_mesh, decay_factor_mesh, q_linear_vals, cmap='Oranges', alpha=0.6, edgecolor='none')
ax.plot_surface(M_mesh, decay_factor_mesh, q_parabolic_vals, cmap='Purples', alpha=0.6, edgecolor='none')
ax.plot_surface(M_mesh, decay_factor_mesh, q_logarithmic_vals, cmap='Reds', alpha=0.6, edgecolor='none')

# Add contour (terrain map) effects
ax.contour3D(M_mesh, decay_factor_mesh, q_brachistochrone_vals, 20, cmap='Blues', linestyles="solid")
ax.contour3D(M_mesh, decay_factor_mesh, q_iterative_vals, 20, cmap='Greens', linestyles="solid")
ax.contour3D(M_mesh, decay_factor_mesh, q_linear_vals, 20, cmap='Oranges', linestyles="solid")
ax.contour3D(M_mesh, decay_factor_mesh, q_parabolic_vals, 20, cmap='Purples', linestyles="solid")
ax.contour3D(M_mesh, decay_factor_mesh, q_logarithmic_vals, 20, cmap='Reds', linestyles="solid")

# Set axis labels
ax.set_xlabel('Batch Size M')
ax.set_ylabel('Decay Factor')
ax.set_zlabel('q Value')

plt.show()
