import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Function representing the brachistochrone curve decay
def decayBrachistochrone(ratio, maxSamples, minSamples):
    theta = math.pi / 2 * ratio
    brachistochrone_factor = math.sqrt(math.sin(theta))
    return minSamples + (maxSamples - minSamples) * brachistochrone_factor

def decayLinear(ratio, maxSamples, minSamples):
    return minSamples + (maxSamples - minSamples) * ratio

def decayParabola(ratio, maxSamples, minSamples):
    parabola_factor = math.sqrt(ratio)
    return minSamples + (maxSamples - minSamples) * parabola_factor

def decayLog(ratio, maxSamples, minSamples):
    decay_factor = np.log(1 + lambda_ * ratio) / np.log(1 + lambda_)
    return minSamples + (maxSamples - minSamples) * decay_factor

def decayIterat(iteration):
    if 0.1 <= iteration < 0.3:
        return 80
    elif 0.3 <= iteration < 0.5:
        return 100
    elif 0.5 <= iteration < 0.8:
        return 130
    elif 0.8 <= iteration <= 1.0:
        return 199
    else:
        return 0
# Generate ratio values from 0 to 1
ratios = np.linspace(0, 1, 10000)

# Assuming maxSamples and minSamples values for demonstration
maxSamples = 199
minSamples = 1
lambda_= 500

# Compute decay values for each ratio
valuesBrachistochrone = [decayBrachistochrone(r, maxSamples, minSamples) for r in ratios]
valuesLinear = [decayLinear(r, maxSamples, minSamples) for r in ratios]
valuesParabola = [decayParabola(r, maxSamples, minSamples) for r in ratios]
valuesLog = [decayLog(r, maxSamples, minSamples) for r in ratios]
valuesIterat = [decayIterat(r) for r in ratios]

font_props = {'family': 'times new roman', 'size': 14}

plt.figure(figsize=(5,5))
# Plotting the decay curve
plt.plot(ratios,valuesBrachistochrone, label="FIT* - B", color='#9DD0C7',linewidth=2.5)
plt.plot(ratios,valuesIterat, label="FIT* - I", color='#8AB1D2',linewidth=2.5)
plt.plot(ratios,valuesLinear,label="FIT* - L", color='#D9BDD8',linewidth=2.5)
plt.plot(ratios,valuesParabola, label="FIT* - P", color='#9180AC',linewidth=2.5)
plt.plot(ratios,valuesLog, label="FIT* - SL", color='#E58579',linewidth=2.5)
# plt.xlabel('1',fontproperties=font_props)
# plt.ylabel('1',fontproperties=font_props)
plt.yticks(fontproperties = 'Times New Roman', size = 10)
plt.xticks(fontproperties = 'Times New Roman', size = 10)
# plt.title("Decay me")
plt.legend(prop=font_props)
# plt.grid(True)
plt.savefig('/home/liding/Documents/decay_method.pdf',dpi=800 ,bbox_inches='tight') # 保存成PDF放大后不失真（默认保存在了当前文件夹下）
plt.show()
