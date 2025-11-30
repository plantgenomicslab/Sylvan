import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

data = {'OOB': [2.6700843746674252e-05, 0.05903881785117915, 0.06016962714034246],
 'discarded': [7004, 9322, 9335],
 'discarded_buscos': [0, 19, 19],
 'kept': [30448, 34395, 34408],
 'kept_buscos': [1507, 1544, 1544]}

data["Discarded Genes"] = data["discarded"]
data["Kept Genes"] = data["kept"]
data["Kept Buscos "] = data["kept_buscos"]
data["Discarded Buscos"] = data["discarded_buscos"]
data.pop("discarded")
data.pop("kept")
data.pop("kept_buscos")
data.pop("discarded_buscos")

def format_with_commas(value, pos):
	if value >= 1:
		return "{:,.0f}".format(value)
	else:
		return "{:.10f}".format(value).rstrip('0')

fig, axs = plt.subplots(3, 2, figsize=(8, 7))

# Flatten the 2D array of subplots for easier indexing
axs = axs.flatten()

# Plot each item in data
for i, (key, values) in enumerate(data.items()):
	ax = axs[i]
	x_values = range(1, len(values) + 1)
	ax.plot(x_values, values, marker='o')
	ax.set_xlabel('Iteration')
	ax.yaxis.set_major_formatter(FuncFormatter(format_with_commas))
	ax.set_title(f'{key}')

fig.delaxes(axs[5])
# Adjust layout
plt.tight_layout()

# Show the plots
# plt.savefig('data_plots.png')
plt.show()