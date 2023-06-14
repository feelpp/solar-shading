import pandas as pd
import matplotlib.pyplot as plt

# read the csv file
data = pd.read_csv('/home/u4/csmi/2022/devora/solar-shading/benchmark/distributions/comparison.csv')

# plot the method
plt.figure(figsize=(20,12))

methods = []
avg_times = []

# unique colors for each unique method
colors = ['orange', 'blue', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

for (flag, method), time in data.groupby(['Flag', 'Method']):
    AvgTime = time['Duration'].mean()
    methods.append(f"{flag}_{method}")
    avg_times.append(AvgTime)

# Sort methods and times based on the times
methods_sorted_with_times = sorted(zip(methods, avg_times), key=lambda x: x[1])
methods_sorted, avg_times_sorted = zip(*methods_sorted_with_times)

# color mapping for each unique method
color_mapping = {method.split('_')[-1]: colors[i % len(colors)] for i, method in enumerate(set(methods))}

# assign colors
colors_assigned = [color_mapping[method.split('_')[-1]] for method in methods_sorted]

bars = plt.bar(methods_sorted, avg_times_sorted, width=0.3, color=colors_assigned, edgecolor='white', linewidth=1.2)

plt.xlabel('Method', fontsize=14)
plt.ylabel('Average Time (milliseconds)', fontsize=14)
plt.title('Comparison of different methods and flags', fontsize=16)

# Create legend
for i, method in enumerate(set(methods)):
    plt.bar(0, 0, color=color_mapping[method.split('_')[-1]], label=method.split('_')[-1], edgecolor='black', linewidth=1.2)

# Rotate x-labels to 45 degrees and decrease the font size
plt.xticks(rotation=45, ha='right', fontsize=8)

plt.tight_layout()
plt.savefig('comparison.png')
plt.show()
