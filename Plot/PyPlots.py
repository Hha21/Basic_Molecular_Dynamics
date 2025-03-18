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