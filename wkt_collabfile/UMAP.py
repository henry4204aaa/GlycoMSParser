
#重要參數umap_model = umap.UMAP(n_components=2, random_state=42)
#N'NGdata1014.csv' f1 0.9  'Odata1014.csv' 0.7
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report
import plotly.express as px
import plotly.graph_objs as go
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
import umap

# Read data
mode = 1
if mode == 1:
    file_path = 'NGdata1014.csv'
    #file_path = 'Odata1014.csv'
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

# Perform UMAP (2 components)
umap_model = umap.UMAP(n_components=2, random_state=42)
umap_result = umap_model.fit_transform(X)

# Train and evaluate a simple classifier on UMAP-transformed data
X_train_umap, X_test_umap, y_train_umap, y_test_umap = train_test_split(umap_result, y, test_size=0.3, random_state=42)
classifier_umap = LogisticRegression(max_iter=1000)
classifier_umap.fit(X_train_umap, y_train_umap)
y_pred_umap = classifier_umap.predict(X_test_umap)
print("Classification Report for UMAP-transformed data:")
print(classification_report(y_test_umap, y_pred_umap))

# Create a DataFrame with UMAP results and Structure labels
umap_df = pd.DataFrame(umap_result, columns=['UMAP1', 'UMAP2'])
umap_df['unique_ID'] = unique_ids  # Retain unique_ID
umap_df['Structure'] = data['Structure']  # Retain original Structure values

# Plotly 2D scatter plot, color and legend by Structure
fig = go.Figure()
for structure_label in umap_df['Structure'].unique():
    structure_data = umap_df[umap_df['Structure'] == structure_label]
    fig.add_trace(go.Scatter(
        x=structure_data['UMAP1'],
        y=structure_data['UMAP2'],
        mode='markers',
        marker=dict(size=6),
        name=f'Structure {structure_label}',
        text=structure_data[['unique_ID', 'Structure']].apply(lambda row: f"ID: {row['unique_ID']}, Structure: {row['Structure']}", axis=1),
        hoverinfo='text'
    ))

# Update layout
fig.update_layout(
    xaxis_title='UMAP1',
    yaxis_title='UMAP2',
    margin=dict(l=0, r=0, b=0, t=0),
    legend=dict(x=-0.1, y=0.5)  # Move legend to the left side
)

# Show plot
fig.show()

'''              precision    recall  f1-score   support

           0       0.50      1.00      0.67         1
           1       1.00      1.00      1.00        12
           2       1.00      1.00      1.00         5
           3       1.00      1.00      1.00        11
           4       0.78      1.00      0.88         7
           5       1.00      1.00      1.00         9
           6       1.00      1.00      1.00         7
           7       1.00      1.00      1.00         5
           8       0.00      0.00      0.00         1
           9       0.91      1.00      0.95        10
          10       1.00      1.00      1.00         3
          11       0.93      1.00      0.97        14
          12       0.88      0.64      0.74        11
          13       1.00      1.00      1.00         1
          14       0.00      0.00      0.00         1
          15       1.00      1.00      1.00         5
          16       1.00      1.00      1.00         3
          17       0.00      0.00      0.00         1
          18       0.75      1.00      0.86         6
          19       0.50      1.00      0.67         2
          20       0.00      0.00      0.00         2
          21       0.50      1.00      0.67         1
          22       0.73      1.00      0.85        11
          23       0.00      0.00      0.00         1
          24       1.00      1.00      1.00        13
          26       1.00      1.00      1.00         4
          27       1.00      1.00      1.00         6
          28       1.00      1.00      1.00         9
          29       1.00      1.00      1.00         3
          30       1.00      1.00      1.00         3
          31       0.00      0.00      0.00         1
          32       1.00      0.75      0.86         4
          33       0.00      0.00      0.00         2
          34       0.00      0.00      0.00         1
          35       1.00      1.00      1.00         4
          36       1.00      1.00      1.00         4
          37       1.00      0.75      0.86         4
          38       0.86      1.00      0.92         6
          39       1.00      1.00      1.00         2
          40       1.00      1.00      1.00         5
          41       1.00      1.00      1.00         2
          42       1.00      1.00      1.00         2

    accuracy                           0.92       205
   macro avg       0.75      0.79      0.76       205
weighted avg       0.89      0.92      0.90       205

'''