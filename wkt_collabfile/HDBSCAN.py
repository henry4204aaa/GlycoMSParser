import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
import hdbscan
import plotly.express as px
from sklearn.metrics import classification_report
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
import plotly.colors as pc
import plotly.graph_objs as go

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
y = training_data['encoded_class']  # Labels
unique_ids = training_data['unique_ID']  # Sample unique_ID
structures = training_data['Structure']  # Original structures

# Perform PCA (3 components)
pca = PCA(n_components=3)
pca_result = pca.fit_transform(X)

# Perform HDBSCAN clustering
hdb = hdbscan.HDBSCAN(min_cluster_size=5, min_samples=5)
cluster_labels = hdb.fit_predict(X)

# Create a DataFrame with PCA results and cluster labels
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2', 'PC3'])
pca_df['Cluster'] = cluster_labels
pca_df['unique_ID'] = unique_ids  # Retain unique_ID
pca_df['Structure'] = data['Structure']  # Retain original Structure values

# Generate a color palette with more distinct and colorful colors
num_clusters = len(pca_df['Cluster'].unique())
colors = px.colors.qualitative.Pastel1 if num_clusters <= 9 else px.colors.qualitative.Set3 if num_clusters <= 12 else px.colors.qualitative.Dark24

# Plotly 3D scatter plot, color by cluster labels, legend by Structure
fig = go.Figure()
for structure_label in pca_df['Structure'].unique():
    structure_data = pca_df[pca_df['Structure'] == structure_label]
    fig.add_trace(go.Scatter3d(
        x=structure_data['PC1'],
        y=structure_data['PC2'],
        z=structure_data['PC3'],
        mode='markers',
        marker=dict(size=2, color=structure_data['Cluster'].apply(lambda x: colors[x % len(colors)])),
        name=f'Structure {structure_label}',
        text=structure_data[['unique_ID', 'Structure', 'Cluster']].apply(lambda row: f"ID: {row['unique_ID']}, Structure: {row['Structure']}, Cluster: {row['Cluster']}", axis=1),
        hoverinfo='text'
    ))

# Update layout
fig.update_layout(
    scene=dict(
        xaxis_title='PC1',
        yaxis_title='PC2',
        zaxis_title='PC3'
    ),
    margin=dict(l=0, r=0, b=0, t=0)
)

# Show plot
fig.show()

