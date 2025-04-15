import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
from sklearn.metrics import classification_report
import plotly.express as px
import plotly.graph_objs as go
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression

# Read data
mode = 1
if mode == 1:
    file_path = 'NGdata1014.csv'
    column_to_drop = ['Structure', 'ID', 'Source', 'IUPACname(optional)', 'Glycanannotation2', 'unique_ID']
    column_of_class = 'Structure'

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
y = training_data['encoded_class']  # Labels
unique_ids = training_data['unique_ID']  # Sample unique_ID
structures = training_data['Structure']  # Original structures

# Perform PCA (3 components)
pca = PCA(n_components=3)
pca_result = pca.fit_transform(X)

# Train and evaluate a simple classifier on PCA-transformed data
X_train_pca, X_test_pca, y_train_pca, y_test_pca = train_test_split(pca_result, y, test_size=0.3, random_state=42)
classifier_pca = LogisticRegression(max_iter=1000)
classifier_pca.fit(X_train_pca, y_train_pca)
y_pred_pca = classifier_pca.predict(X_test_pca)
print("Classification Report for PCA-transformed data:")
print(classification_report(y_test_pca, y_pred_pca))

# Create a DataFrame with PCA results and Structure labels
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2', 'PC3'])
pca_df['unique_ID'] = unique_ids  # Retain unique_ID
pca_df['Structure'] = data['Structure']  # Retain original Structure values

# Plotly 3D scatter plot, color and legend by Structure
fig = go.Figure()
for structure_label in pca_df['Structure'].unique():
    structure_data = pca_df[pca_df['Structure'] == structure_label]
    fig.add_trace(go.Scatter3d(
        x=structure_data['PC1'],
        y=structure_data['PC2'],
        z=structure_data['PC3'],
        mode='markers',
        marker=dict(size=2),
        name=f'Structure {structure_label}',
        text=structure_data[['unique_ID', 'Structure']].apply(lambda row: f"ID: {row['unique_ID']}, Structure: {row['Structure']}", axis=1),
        hoverinfo='text'
    ))

# Update layout
fig.update_layout(
    scene=dict(
        xaxis_title='PC1',
        yaxis_title='PC2',
        zaxis_title='PC3'
    ),
    margin=dict(l=0, r=0, b=0, t=0),
    legend=dict(x=-0.1, y=0.5)  # Move legend to the left side
)

# Show plot
fig.show()
