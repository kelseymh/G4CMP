from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

file_path = r"/your/file/path/mesh.yaml"

def extract_phonon_data_with_modes(file_path):
    frequencies = []
    q_point_weights = []
    mode_numbers = []  # To store mode numbers
    
    current_q_weight = None
    current_mode = None
    skip_eigenvector_block = False  # Flag to skip the eigenvector block

    with open(file_path, 'r') as f:
        for line in tqdm(f, desc="Reading file line by line"):
            line = line.strip()

            # Detect the start of the eigenvector block and skip it
            if "eigenvector:" in line:
                skip_eigenvector_block = True
                continue  # Skip the eigenvector line itself
            
            # Stop skipping lines when we reach the end of the eigenvector block
            if skip_eigenvector_block:
                if line.startswith("- # atom"):  # Check if we encounter a new mode/section
                    skip_eigenvector_block = True
                elif line.startswith("- #"):
                    current_mode = int(line.split("#")[1].strip())
                    #mode_numbers.append(current_mode)

                    skip_eigenvector_block = False
                continue  # Continue skipping until we find the end of the block

            # Extract q-point weight
            if "weight:" in line:
                current_q_weight = float(line.split(":")[1].strip())
              #  q_point_weights.append(current_q_weight)
            
            # Extract mode number
            if line.startswith("- #") and not line.startswith("- # atom"):  # Mode number line
               current_mode = int(line.split("#")[1].strip())

                
            # Extract frequency
            if "frequency:" in line:
                frequency = float(line.split(":")[1].strip())
                frequencies.append(frequency)
                q_point_weights.append(current_q_weight)
                mode_numbers.append(current_mode)  # Associate mode with frequency
    
    return np.array(frequencies), np.array(q_point_weights), np.array(mode_numbers)

def compute_dos(frequencies, q_point_weights, num_bins=100):
    # Create a histogram of frequencies to calculate DOS
    dos, bin_edges = np.histogram(frequencies, bins=num_bins, weights=q_point_weights, density=False)
    
    
    # Get bin centers for plotting
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    
    return bin_centers, dos

def plot_dos(bin_centers, total_dos, mode_1_frequencies, mode_1_dos, mode_2_frequencies, mode_2_dos, mode_3_frequencies, mode_3_dos):
    plt.plot(bin_centers, total_dos, label='Total DOS', color='blue')
    plt.plot(mode_1_frequencies, mode_1_dos, label='Mode 1 Contribution', linestyle='--', color='red')
    plt.plot(mode_2_frequencies, mode_2_dos, label='Mode 2 Contribution', linestyle='--', color='green')
    plt.plot(mode_3_frequencies, mode_3_dos, label='Mode 3 Contribution', linestyle='--', color='orange')
    #plt.plot(ff, dd*300, color='black', linewidth=1.5, label='Density of States (DOS)')
    plt.xlabel('Frequency (THz)')
    plt.xlim(0,10)
    plt.ylim(0,1500)
    plt.ylabel('DOS (states/THz)')
    plt.title('Phonon DOS and Mode Contributions')
    plt.legend()
    plt.show()

# Main function to extract data, compute DOS, and plot mode-specific contributions
    
print("Extracting frequencies, weights, and mode numbers manually...")
frequencies, q_point_weights, mode_numbers = extract_phonon_data_with_modes(file_path)

print("Computing total DOS...")
bin_centers, total_dos = compute_dos(frequencies, q_point_weights, num_bins=math.ceil((max(frequencies) - min(frequencies)) * 50))
plt.plot(bin_centers, total_dos, label='$Total$ $DOS$', color='black')

# Plot individual modes 
for i in range(1, 4):
    mode_frequencies = frequencies[mode_numbers == i]
    mode_weights = q_point_weights[mode_numbers == i]
    
    bins = math.ceil((max(mode_frequencies) - min(mode_frequencies)) * 50)
    m, mode_dos = compute_dos(mode_frequencies, mode_weights, num_bins=bins)
    
    labels = {1: r"$Longitudinal$", 2: r"$Transverse$ $Fast$", 3: r"$Transverse$ $Slow$"}
    linestyles = {1: '-', 2: '-', 3: '-'}
    
    plt.plot(m, mode_dos, linestyle=linestyles[i], label=labels[i])

# Add reference lines
a=max(frequencies[mode_numbers == 3])
plt.axvline(a, label=r"$\omega_a$", color="red", linewidth=0.8, linestyle='--')

# Axis labels and limits
plt.xlabel('$Frequency$ $(THz)$')
plt.ylabel('$DOS$ $(states/THz)$')
plt.xlim(0, 5)
plt.ylim(0, 5000)
plt.title('Phonon DOS', fontsize=12)

# Tight layout for compact spacing
plt.legend(loc='upper left', frameon=False)
plt.tight_layout()

# Save for publication
plt.show()

