import numpy as np
import matplotlib.pyplot as plt

# Function representing the logarithmic decay
def decay(ratio, lambda_, maxSamples, minSamples):
    decay_factor = np.log(1 + lambda_ * ratio) / np.log(1 + lambda_)
    return maxSamples - (maxSamples - minSamples) * decay_factor

# Parameters
minSamples = 70
maxSamples = 300
lambda_ = 1000

# Generate ratio values
ratios = np.linspace(0, 1, 500)
batch_sizes = [decay(r, lambda_, maxSamples, minSamples) for r in ratios]

# Plotting
plt.plot(ratios, batch_sizes, label=f"Lambda={lambda_}")
plt.title("Logarithmic Decay of batchSize vs. ratio")
plt.xlabel("Ratio")
plt.ylabel("BatchSize")
plt.legend()
plt.grid(True)
plt.show()
