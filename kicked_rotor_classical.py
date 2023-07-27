import numpy as np
import matplotlib.pyplot as plt

# Define the parameter and number of iterations
k = 0.5
N = 1000

# Function to calculate the trajectory for a given initial condition
def calculate_trajectory(theta0, p0, k, N):
    theta = np.zeros(N)
    p = np.zeros(N)
    theta[0] = theta0
    p[0] = p0

    for n in range(N - 1):
        theta[n + 1] = (theta[n] + p[n]) % (2 * np.pi)
        p[n + 1] = p[n] + k * np.sin(theta[n + 1])

    return theta, p

# initial points uniformly spaced in phase space
theta0_vals = np.linspace(0, 2 * np.pi, 8)
p0_vals = np.linspace(0, 2 * np.pi, 20)

# Plotting for different initial conditions
for theta0 in theta0_vals:
    for p0 in p0_vals:
        theta, p = calculate_trajectory(theta0, p0, k, N)
        plt.plot(theta, p, ".", markersize=1)

plt.title(f"Standard map for k = {k}")
plt.xlabel("θ")
plt.ylabel("P")
plt.ylim(0, 2 * np.pi + 0.1)
plt.show()

