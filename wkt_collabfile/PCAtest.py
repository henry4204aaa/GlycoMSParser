
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA

# Read data
mode = 1
if mode == 1:
    file_path = 'NGdata1014.csv'
    column_to_drop = ['Structure', 'ID', 'Source', 'IUPACname(optional)', 'Glycanannotation2', 'unique_ID']
    column_of_class = 'Structure'
elif mode == 2:
    file_path = 'smallCandyset.csv'
    column_to_drop = ['RT', 'filename', 'GlycoPost_ID', 'unique_ID']
    column_of_class = 'glycan'

# Load data
data = pd.read_csv(file_path)

# Filter classes with less than 5 samples
def filter_data_for_training(data, column_of_class, min_samples=5):
    class_counts = data[column_of_class].value_counts()
    valid_classes = class_counts[class_counts >= min_samples].index
    filtered_data = data[data[column_of_class].isin(valid_classes)].copy()
    return filtered_data

# Filter data
training_data = filter_data_for_training(data, column_of_class)

# Encode class labels
label_encoder = LabelEncoder()
training_data['encoded_class'] = label_encoder.fit_transform(training_data[column_of_class])

# Split features and labels
X = training_data.drop(columns=['encoded_class'] + column_to_drop)  # Features, remove unnecessary columns

# Perform PCA (3 components)
pca = PCA(n_components=3)
pca_result = pca.fit_transform(X)

# Create a DataFrame with PCA results
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2', 'PC3'])
pca_df['unique_ID'] = training_data['unique_ID']  # Retain unique_ID
pca_df['Structure'] = training_data['Structure']  # Retain original Structure values

# Output the PCA result DataFrame
print(pca_df)