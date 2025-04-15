#統計peak list幾種peak
import csv
import ast

# Load the dataset
file_path = 'full_datasethead100.csv'  # 請將這裡替換為您的原始資料檔案路徑

# Extract all keys from the 'peak_d' dictionaries in each row
all_keys = set()
with open(file_path, 'r') as file:
    reader = csv.DictReader(file)
    for row in reader:
        peak_d = row['peak_d']
        if peak_d:
            peak_dict = ast.literal_eval(peak_d)
            all_keys.update(peak_dict.keys())

# Sort the keys
sorted_keys = sorted(all_keys)
print('len():',len(sorted_keys))# len(): 19690
# Display the sorted keys
print(sorted_keys)