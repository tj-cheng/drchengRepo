import numpy as np
import qutip as qt

# Define the number of spins
N = 4

# Define Pauli matrices
sigma_x = qt.sigmax()
sigma_y = qt.sigmay()
sigma_z = qt.sigmaz()

# Define the Heisenberg Hamiltonian parameters
J = 1.0  # Coupling strength

# Create the unperturbed Hamiltonian (Heisenberg model)
H_0 = sum([J * (qt.tensor(sigma_x, sigma_x, *[qt.qeye(2) for _ in range(N - 4 - i)], sigma_y, sigma_y, *[qt.qeye(2) for _ in range(i)]) +
                sigma_y, sigma_y, *[qt.qeye(2) for _ in range(i)], sigma_x, sigma_x, *[qt.qeye(2) for _ in range(N - 4 - i)]) for i in range(N - 3)])

# Define the perturbing Hamiltonian (e.g., a magnetic field term)
B = 0.5  # Strength of the magnetic field
H_perturb = B * sum([sigma_z, *[qt.qeye(2) for _ in range(N - 1)]])

# Calculate the ground state and eigenvalues of the unperturbed Hamiltonian
eigenstates, eigenenergies = H_0.eigenstates()

# Ground state energy
E0 = eigenenergies[0]

# Calculate the second-order perturbation correction
delta_E_2nd_order = 0.0

for i in range(1, len(eigenstates)):
    H_perturb_element = qt.expect(H_perturb, eigenstates[i])
    delta_E_2nd_order += (H_perturb_element ** 2) / (E0 - eigenenergies[i])

# Print the result
print(f"Second-Order Perturbation Correction: {delta_E_2nd_order:.4f}")

