import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 14,  # Base font size
    "axes.labelsize": 16,  # Label size
    "xtick.labelsize": 14,  # X-axis tick labels
    "ytick.labelsize": 14,  # Y-axis tick labels
    "legend.fontsize": 14,  # Legend text size
    "lines.linewidth": 2.5  # Thicker lines for clarity
})

# Load Particle Trajectory Data
filename = "./particles.txt"
df = pd.read_csv(filename, delim_whitespace=True)

particle_0 = df[df["ID"] == 0]
particle_1 = df[df["ID"] == 1]

# Particle Trajectory Plot
plt.figure(figsize=(9, 7))  
plt.plot(particle_0["X"], particle_0["Y"], label="Particle 1", linestyle="-", color="r", marker="o", markersize=1)
plt.plot(particle_1["X"], particle_1["Y"], label="Particle 2", linestyle="-", color="b", marker="s", markersize=1)

plt.xlim([0, 20])
plt.ylim([0, 20])

plt.xlabel("X Position (Angstroms)", fontsize=16, fontweight='bold')
plt.ylabel("Y Position (Angstroms)", fontsize=16, fontweight='bold')
plt.legend(loc="best", frameon=True, fancybox=True, shadow=True)
plt.grid(visible=True, linestyle="--", linewidth=0.7, alpha=0.7)

plt.show()

# Load Kinetic Energy Data
file_path = "./kinetic_energy.txt"  
data = pd.read_csv(file_path, delim_whitespace=True)

# Kinetic Energy Plot
plt.figure(figsize=(10, 6))
plt.plot(data["TIME"], data["KE"], color="r", linestyle="-", marker="o", markersize=1)

plt.xlim([0, 50])
plt.ylim([0, 30])

plt.xlabel("Time", fontsize=16, fontweight='bold')
plt.ylabel("Kinetic Energy", fontsize=16, fontweight='bold')
plt.grid(visible=True, linestyle="--", linewidth=0.7, alpha=0.7)

plt.show()

plt.rcParams.update({
    "font.size": 14,  # Base font size
    "axes.labelsize": 16,  # Label size
    "xtick.labelsize": 14,  # X-axis tick labels
    "ytick.labelsize": 14,  # Y-axis tick labels
    "legend.fontsize": 14,  # Legend text size
    "lines.linewidth": 2.5  # Thicker lines for clarity
})

# Data extracted from your screenshot
num_cores = [1, 2, 4, 8, 16, 24, 32, 40, 48]

# Raw timings for each run
#runtime_N15000_notemp = [2211.56, 780.90, 780.47, 779.32, 1916.35, 1995.36, 1775.90, 1618.44, 1619.37]
runtime_N15000_temp07617_T1 = [2226.41, 1467.01, 1459.99, 1454.64, 1878.45, 851.22, 1311.91, 1311.46, 756.68]
#runtime_N15000_temp07617_T015 = [3654.71, 2423.89, 2423.81, 2424.38, 2078.80, 2073.01, 1415.69, 2204.70, 1264.10]
runtime_N50000_temp01904_T015 = [3654.71, 2137.66, 2426.75, 2426.18, 2425.00, 2206.63, 1409.49, 2620.69, 1261.39]

# Normalizing runtimes by single-core runtime for each run
#speedup_N15000_notemp = [runtime_N15000_notemp[0]/t for t in runtime_N15000_notemp]
speedup_N15000_temp07617_T1 = [runtime_N15000_temp07617_T1[0]/t for t in runtime_N15000_temp07617_T1]
#speedup_N15000_temp07617_T015 = [runtime_N15000_temp07617_T015[0]/t for t in runtime_N15000_temp07617_T015]
speedup_N50000_temp01904_T015 = [runtime_N50000_temp01904_T015[0]/t for t in runtime_N50000_temp01904_T015]

# Plotting
plt.figure(figsize=(9, 7))

#plt.plot(num_cores, speedup_N15000_notemp, '-o', label='N=15000, no temp set')
plt.plot(num_cores, speedup_N15000_temp07617_T1, '-o', color = 'b', label='N=15000, temp=0.7617')
#plt.plot(num_cores, speedup_N15000_temp07617_T015, '-o', label='N=15000, temp=0.7617 (T=0.15)')
plt.plot(num_cores, speedup_N50000_temp01904_T015, '-x', color = 'r', label='N=50000, temp=0.1904')

# Labels and title
plt.xlabel("Number of Cores", fontsize=16, fontweight='bold')
plt.ylabel("Relative Speed-Up 1/($t / t_0$)", fontsize=16, fontweight='bold')
plt.xticks(num_cores)
plt.grid(visible=True, linestyle="--", linewidth=0.7, alpha=0.7)
plt.legend(loc="best", frameon=True, fancybox=True, shadow=True)

# Show plot
plt.tight_layout()
plt.show()