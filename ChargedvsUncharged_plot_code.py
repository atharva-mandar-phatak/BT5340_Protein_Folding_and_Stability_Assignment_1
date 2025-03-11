import pickle
import matplotlib.pyplot as plt
import numpy as np

# Define the set of charged amino acids (three-letter codes)
CHARGED = {"LYS", "ARG", "HIS", "ASP", "GLU"}

# Load precomputed final_results from the pickle file
with open("/Users/hiteshkandarpa/Desktop/Acads/Protein folding/Assignment-1/final_results.pkl", "rb") as f:
    final_results = pickle.load(f)

# Containers for the central residue phi/psi angles
charged_context = []    # (phi, psi) for ASN where left and right windows have at least one charged residue
uncharged_context = []  # (phi, psi) for ASN where that condition is not met

# Process each PDB and chain in final_results
for pdb_id, chains in final_results.items():
    for chain_id, residues in chains.items():
        sorted_residues = sorted(residues, key=lambda x: x[1])  # Sort by residue number
        if len(sorted_residues) < 7:
            continue
        for i in range(3, len(sorted_residues) - 3):
            window = sorted_residues[i - 3: i + 4]  # a, b, c, X, d, e, f
            center = window[3]  # Central residue X
            res_name, res_number, phi, psi = center
            if res_name != "ASN" or phi == "N/A" or psi == "N/A":
                continue
            left_charged = any(r[0] in CHARGED for r in window[:3])
            right_charged = any(r[0] in CHARGED for r in window[4:])
            if left_charged and right_charged:
                charged_context.append((phi, psi))
            else:
                uncharged_context.append((phi, psi))

# Convert lists to numpy arrays
charged_context = np.array(charged_context)
uncharged_context = np.array(uncharged_context)

# Define bin edges (5° bins from -180° to 180°)
bin_edges = np.arange(-180, 185, 5)

# Function to compute frequency-based coloring
def compute_frequency(x, y, bins):
    hist, xedges, yedges = np.histogram2d(x, y, bins=[bins, bins])
    hist = hist.T  # Transpose to match coordinate system
    return hist, xedges, yedges

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

# Plot Charged Context Scatter Plot
if charged_context.size > 0:
    phi_charged, psi_charged = charged_context[:, 0], charged_context[:, 1]
    hist_charged, xedges, yedges = compute_frequency(phi_charged, psi_charged, bin_edges)
    im1 = axes[0].imshow(hist_charged, cmap="plasma", origin="lower",
                          extent=[-180, 180, -180, 180], aspect="auto")
    axes[0].set_title("Charged Context")
    axes[0].set_xlabel("Phi (°)")
    axes[0].set_ylabel("Psi (°)")
    cbar1 = plt.colorbar(im1, ax=axes[0])
    cbar1.set_label("Residue Count")

# Plot Uncharged Context Scatter Plot
if uncharged_context.size > 0:
    phi_uncharged, psi_uncharged = uncharged_context[:, 0], uncharged_context[:, 1]
    hist_uncharged, _, _ = compute_frequency(phi_uncharged, psi_uncharged, bin_edges)
    im2 = axes[1].imshow(hist_uncharged, cmap="plasma", origin="lower",
                          extent=[-180, 180, -180, 180], aspect="auto")
    axes[1].set_title("Uncharged Context")
    axes[1].set_xlabel("Phi (°)")
    cbar2 = plt.colorbar(im2, ax=axes[1])
    cbar2.set_label("Residue Count")

plt.suptitle("Ramachandran Scatter Plots (Colored by Residue Count)")
plt.tight_layout()
plt.show()
