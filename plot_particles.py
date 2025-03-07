import pandas as pd
import matplotlib.pyplot as plt

filename = "particles.txt"
df = pd.read_csv(filename, delim_whitespace=True)

particle_0 = df[df["ID"] == 0]
particle_1 = df[df["ID"] == 1]

plt.figure(figsize=(8, 6))
plt.plot(particle_0["X"], particle_0["Y"], label="Particle 0", marker="o", linestyle="-")
plt.plot(particle_1["X"], particle_1["Y"], label="Particle 1", marker="s", linestyle="-")

plt.xlabel("X Position (Angstroms)")
plt.ylabel("Y Position (Angstroms)")
plt.title("Particle Motion in X-Y Plane")
plt.legend()
plt.grid()

plt.show()