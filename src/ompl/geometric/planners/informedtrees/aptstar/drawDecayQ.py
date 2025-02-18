import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# Function to align the range
def remap(x, old_min, old_max, new_min, new_max):
    return new_min + (x - old_min) * (new_max - new_min) / (old_max - old_min)
# Function representing an exponential decay
def exponential_function(ratio, maxCharges, minCharges):
    decay_factor = 1 - math.exp(-5 * ratio)  # Modify decay speed with constant (e.g., 5)
    return maxCharges - (maxCharges - minCharges) * decay_factor
# Function representing a polynomial decay
def polynomial_function(ratio, maxCharges, minCharges, power=2):
    decay_factor = ratio ** power  # Polynomial decay with adjustable power
    return maxCharges - (maxCharges - minCharges) * decay_factor
# Function representing a logarithmic decay
def logarithmic_function(ratio, maxCharges, minCharges):
    decay_factor = math.log(1 + 9 * ratio) / math.log(10)  # Adjust base and multiplier
    return maxCharges - (maxCharges - minCharges) * decay_factor
# Function representing a hyperbolic tangent (tanh) decay
def tanh_function(ratio, maxCharges, minCharges):
    decay_factor = (math.tanh(6 * (ratio - 0.5)) + 1) / 2  # Shift and scale tanh
    return maxCharges - (maxCharges - minCharges) * decay_factor
def iterat_function(iteration):
    if 0.1 <= iteration < 0.3:
        return 1.9
    elif 0.3 <= iteration < 0.5:
        return 1.5
    elif 0.5 <= iteration < 0.8:
        return 0.8
    elif 0.8 <= iteration <= 1.0:
        return 0.1
    else:
        return 1.9
# Function representing a sigmoid decay
def sigmoid_function(ratio, maxCharges, minCharges):
    decay_factor = 1 / (1 + math.exp(-10 * (ratio - 0.5)))  # Shift and steepness
    return maxCharges - (maxCharges - minCharges) * decay_factor
# Generate ratio values from 0 to 1
m = np.linspace(1, 199, 10000)
ratios = (m) / (199)
# Assuming maxCharges and minCharges values for demonstration
maxCharges = 1.9
minCharges = 0.1
# Compute decay values for each ratio
values_exponential = [exponential_function(r, maxCharges, minCharges) for r in ratios]
values_polynomial = [polynomial_function(r, maxCharges, minCharges) for r in ratios]
values_logarithmic = [logarithmic_function(r, maxCharges, minCharges) for r in ratios]
values_inter = [iterat_function(r) for r in ratios]
values_tanh = [tanh_function(r, maxCharges, minCharges) for r in ratios]
values_sigmoid = [sigmoid_function(r, maxCharges, minCharges) for r in ratios]
font_props = {'family': 'times new roman', 'size': 14}
# Remap to align
values_logarithmic = [remap(y,0.1,1.87,0.1,1.9) for y in values_logarithmic]
values_exponential = [remap(y,0.1,1.87,0.1,1.9) for y in values_exponential]
plt.figure(figsize=(5, 5))
# Plotting the decay curve
# plt.plot(ratios, values_exponential, label="APT*-E", color='#9DD0C7', linewidth=2.5)
# plt.plot(ratios, values_polynomial, label="APT*-P", color='#8AB1D2', linewidth=2.5)
# # plt.plot(ratios, values_reciprocal, label="APT*-R", color='#D9BDD8', linewidth=2.5)
# plt.plot(ratios, values_logarithmic, label="APT*-L", color='#9180AC', linewidth=2.5)
# plt.plot(ratios, values_inter, label="APT*-I", color='#D9BDD8', linewidth=2.5)
# # plt.plot(ratios, values_sigmoid, label="APT*-S", color='#F4A261', linewidth=2.5)
# # plt.plot(m, values_sigmoid, label="APT*-S", color='#E58579', linewidth=2.5)
# plt.plot(ratios, values_tanh, label="APT*-T", color='#E58579', linewidth=2.5)
plt.plot(m, values_exponential, label="APT*-E", color='#9AC9DB', linewidth=2.5)
plt.plot(m, values_polynomial, label="APT*-P", color='#2878B5', linewidth=2.5)
plt.plot(m, values_logarithmic, label="APT*-L", color='#F8AC8C', linewidth=2.5)
plt.plot(m, values_inter, label="APT*-I", color='#FF8884', linewidth=2.5)
plt.plot(m, values_tanh, label="APT*-T", color='#C82423', linewidth=2.5)
# Set up range of y axis
plt.ylim(bottom=0,top=2.0)
plt.yticks(fontproperties='Times New Roman', size=10)
plt.xticks(fontproperties='Times New Roman', size=10)
plt.legend(prop=font_props)
plt.savefig('/home/liding/Documents/charge_method.pdf', dpi=800, bbox_inches='tight')  # Save as PDF
plt.show()