#而總比數是504895
#（那先畫前20筆 和other類別（504895-加總前20筆））100%
#（然後再畫地21筆-40筆 和other類別（504895-加總前20筆-21筆-40筆加總））100%


import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
file_path = '/mnt/data/glycan_category_counts.csv'
data = pd.read_csv(file_path)

# Calculate total count
total_count = data['Count'].sum()

# Number of entries per pie chart
entries_per_chart = 20

# Iterate through the data in chunks of 20
for i in range(0, len(data), entries_per_chart):
    # Get the subset of data for this pie chart
    subset = data.iloc[i:i + entries_per_chart]
    
    # Calculate the "Others" count
    others_count = total_count - subset['Count'].sum()
    
    # Append the "Others" category to the subset
    subset = subset.append({'glycan': 'Others', 'Count': others_count}, ignore_index=True)
    
    # Plot the pie chart
    plt.figure(figsize=(10, 6))
    plt.pie(subset['Count'], labels=subset['glycan'], autopct='%1.1f%%', startangle=140)
    plt.title(f'Glycan Distribution (Entries {i + 1} to {i + entries_per_chart})')
    plt.axis('equal')
    plt.show()
