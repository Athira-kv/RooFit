import numpy as np
import matplotlib.pyplot as plt

# Load the data, skipping the header
data = np.loadtxt("k_signal_background_results_simfit.txt", skiprows=1)

momentum = data[:, 1]  # Momentum values
signal_values = data[:, 3]  # Signal values
deno = data[:, 2] #denominator 

# Prepare to plot
plt.figure(figsize=(12, 8))

efficiency = signal_values/deno;
error = np.sqrt(efficiency * (1 - efficiency) / deno)

plt.errorbar(momentum, efficiency,yerr=error, marker='o', color='b',linestyle='none')




#plt.title('efficiency (#pi  vs Momentum', fontsize=16)
plt.xlabel('Momentum (GeV/c)', fontsize=14)
plt.ylabel(r'$\epsilon(\pi^+ \rightarrow \pi^+) = \frac{N_{\pi^+ \rightarrow \pi^+}^{Sig}}{N_{{All events}}^{Sig}}$', fontsize=18)
plt.legend()
plt.grid()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(0.7,1.0)
# Style the axes
plt.tick_params(axis='both', which='both', direction='in', length=6, width=2, color='black')
#plt.tight_layout()  # Adjust layout to prevent clipping
#plt.savefig("signal_yield_ratio_category_0_over_4.png")  # Save the figure
plt.show()  # Show the plot