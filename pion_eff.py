import numpy as np
import matplotlib.pyplot as plt

# Load data from the file
data = np.loadtxt('k_signal_background_results_simfit.txt', skiprows=1)

# Extract columns
column2 = data[:, 1]   # Second column (x-axis)
column4 = data[:, 2]   # (common denominator)

numerator_columns = [4,5,6]  
ratios = []
errors = []

particle_labels = {3: 'pi', 4: 'K', 5: 'proton', 6: 'noID'}


for i, num_col in enumerate(numerator_columns):
    numerator = data[:, num_col]  # Current numerator column
    ratio = numerator / column4
    error = np.sqrt(ratio * (1 - ratio) / column4)  # Binomial error
    ratios.append(ratio)
    errors.append(error)

    # Get the label from the particle_labels dictionary
    label = particle_labels.get(num_col, f'Unknown ({num_col})')  # Default label if not found

    # Plot each ratio with error bars
    plt.errorbar(column2, ratio, yerr=error, fmt='o', label=f'$\epsilon$({label})', capsize=5)






plt.xlabel('Momentum (GeV/c)', fontsize=14)
plt.ylabel(r'$\epsilon(\pi^+ \rightarrow j^+) = \frac{N_{\pi^+ \rightarrow j^+}^{Sig}}{N_{{All events}}^{Sig}}$', fontsize=18)
plt.legend()
plt.grid()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim(0.0,0.2)
# Style the axes
plt.tick_params(axis='both', which='both', direction='in', length=6, width=2, color='black')
#plt.tight_layout()  # Adjust layout to prevent clipping
#plt.savefig("signal_yield_ratio_category_0_over_4.png")  # Save the figure
plt.show()  # Show the plot
