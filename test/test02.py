import matplotlib.pyplot as plt
import numpy as np

# Set problem parameters/constants
L = 10  # length of domain
dx = 0.1  # spatial discretization size
nx = int(L / dx)  # number of discretization points
dt = 0.001  # time step
nt = int(60 / 0.001)  # number of time steps


# Define a function for temperature dependent and location dependent thermal diffusivity
def k(x, T):
    # This is just a dummy function for illustration.
    # You should replace this with your actual function.
    # return 0.5 + 0.01 * x + 0.001 * T
    return 1


# Initialize solution: the grid of u
u = np.empty(nx)

# Initial condition everywhere inside the grid
u_initial = 293.15
u.fill(u_initial)

# Boundary conditions
q_left = 100
u[0] = u[1] - q_left * dx / k(0, u[1])
# u[-1] = 293.15

# Time-stepping loop
for t in range(nt - 1):
    u_old = u.copy()
    for i in range(1, nx - 1):
        u[i] = u_old[i] + k(i * dx, u_old[i]) * dt / dx ** 2 * (u_old[i + 1] - 2 * u_old[i] + u_old[i - 1])
    # Heat flux boundary condition at x = 0
    u[0] = u_old[1] + q_left * dx / k(0, u_old[1])
    # Temperature boundary condition at x = L
    # u[-1] = Thot
    if (t * dt) % 5 == 0:
        plt.plot(np.linspace(0, L, nx), u - 293.15, label=f'{t * dt}')

# All done! Plot the solution
plt.plot(np.linspace(0, L, nx), u - 293.15, label=f'{nt * dt}')
# plt.plot(np.linspace(0, L, nx), u)
plt.legend().set_visible(True)
plt.xlabel("Position in rod")
plt.ylabel("Temperature")
plt.title("Heat conduction with variable diffusivity")
plt.show()
