import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict


def main():
    # Read the data from the input file
    data_file = 'k_covvalues.txt'  # Specify your input file
    data = np.loadtxt(data_file, skiprows=1)
    
    # Define categories and momentum bins
    categories = data[:,0]
    momenta = data[:, 1]
    n_as = data[:, 2]
    n_pi_s = data[:, 3]
    n_k_s = data[:, 4]
    n_p_s = data[:, 5]
    n_u_s = data[:, 6]

    numerator_columns = [3,4,5,6]
    ratios = []
    errors = []

    particle_labels = {3: 'pi', 4: 'K', 5: 'proton', 6: 'noID'}
    for i, num_col in enumerate(numerator_columns):
        numerator = data[:, num_col]  # Current numerator column
        ratio = numerator / n_as
        error = np.sqrt(ratio * (1 - ratio) / n_as)  # Binomial error
        ratios.append(ratio)
        errors.append(error)

        # Get the label from the particle_labels dictionary
        label = particle_labels.get(num_col, f'Unknown ({num_col})')  # Default label if not found

        # Plot each ratio with error bars
        plt.errorbar(momenta, ratio, yerr=error, fmt='o', label=f'$\epsilon$({label})', capsize=5)

    plt.xlabel('Momentum (GeV/c)', fontsize=14)
    plt.ylabel(r'$\epsilon(\pi^+ \rightarrow j^+) = \frac{N_{\pi^+ \rightarrow j^+}^{Sig}}{N_{{All events}}^{Sig}}$', fontsize=18)
    plt.legend()
    plt.grid()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0.0,1.0)
    # Style the axes
    plt.tick_params(axis='both', which='both', direction='in', length=6, width=2, color='black')
    #plt.tight_layout()  # Adjust layout to prevent clipping
    #plt.savefig("signal_yield_ratio_category_0_over_4.png")  # Save the figure
    plt.show()  # Show the plot

        
if __name__ == '__main__':
    main()