import pickle
with open("final_results.pkl", "rb") as f:
    final_results = pickle.load(f)

import numpy as np
import matplotlib.pyplot as plt

# Assuming final_results is already defined and has the following structure:
# final_results[pdb_id][chain] = list of tuples (res_name, res_number, phi, psi)

# Collect phi and psi values for asparagine (ASN) residues
phi_list = []
psi_list = []

for pdb_id, chains in final_results.items():
    for chain_id, residues in chains.items():
        for res_name, res_number, phi, psi in residues:
            # Filter for asparagine (ASN) only
            if res_name.upper() == "ASN":
                # Make sure angles are available (they are numbers, not "N/A")
                if phi != "N/A" and psi != "N/A":
                    phi_list.append(phi)
                    psi_list.append(psi)

# Convert lists to numpy arrays
phi_array = np.array(phi_list)
psi_array = np.array(psi_list)

# Define bin edges: from -180 to 180 degrees in 3° increments
x_bins = np.arange(-180, 180 + 3, 3)
y_bins = np.arange(-180, 180 + 3, 3)

# Create a 2D histogram of the phi and psi values
hist, xedges, yedges = np.histogram2d(phi_array, psi_array, bins=[x_bins, y_bins])

# Plot the heatmap of the histogram (Ramachandran plot)
plt.figure(figsize=(8, 6))
plt.imshow(hist.T, origin='lower', extent=[-180, 180, -180, 180],
           cmap='viridis', aspect='auto')
plt.xlabel('Phi (°)')
plt.ylabel('Psi (°)')
plt.title('Ramachandran Plot for Asparagine (ASN)')
plt.colorbar(label='Frequency')
plt.show()
