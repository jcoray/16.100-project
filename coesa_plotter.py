from skaero.atmosphere import coesa
import matplotlib.pyplot as plt
import numpy as np

max_alt = 80e3
altitudes = np.linspace(0, max_alt, num=2000)

ts = []
ps = []
for a in altitudes:
	_, Temp, p, rho = coesa.table(a)
	ts.append(Temp)
	ps.append(p)

# Pressure plot
plt.plot(ps, altitudes)
plt.xlabel("Pressure (Pa)")
plt.ylabel("Altitude (m)")
plt.title("COESA Pressure Model")
plt.grid()
plt.show()

# Temperature plot
plt.plot(ts, altitudes)
plt.xlabel("Temperature (K)")
plt.ylabel("Altitude (m)")
plt.title("COESA Temperature Model")
plt.grid()
plt.show()