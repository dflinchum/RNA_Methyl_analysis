import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('JN_intersected_loci_enhancers2.csv')
loci_per_enh = df.groupby('enhancer_name')['CpG_ID'].count()

# Convert the series to a dataframe
loci_per_enh_df = loci_per_enh.reset_index()
loci_per_enh_df.columns = ['Enhancer', 'Count']

plt.figure(figsize=(10, 6), dpi=300)

# Use seaborn's countplot
sns.countplot(x='Count', data=loci_per_enh_df, palette=['red'])

plt.xlabel('Number of Differentially Methylated Loci')
plt.ylabel('Number of Enhancers')
plt.title('JN: Distribution of Loci per Enhancer')

# Set the desired ticks using matplotlib.pyplot.xticks
desired_ticks = [0,24,49,99,124]
plt.xticks(desired_ticks, rotation=45, ha='right')

# Adjust the x-axis limits  <-- This is the line that extends the limits
plt.xlim(left=-5)  # Extend the left limit by 5 units

# Set the y-axis limit to 400
plt.ylim(top=400)

plt.tight_layout()
plt.show()
plt.savefig('JN_loci_per_enhancer.png')
print("the thick plottens")
