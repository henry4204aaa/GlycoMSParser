import plotly.express as px
from sklearn.decomposition import PCA
import pandas as pd

# 加載數據
file_path = 'All_data.csv'
data = pd.read_csv(file_path)

# 刪除非數值型欄位 'IUPACname(optional)' 和 'Glycanannotation2'
data_cleaned = data.drop(columns=['IUPACname(optional)', 'Glycanannotation2'])

# 提取數值型特徵，並使用均值填補 NaN 值
numeric_columns = data_cleaned.select_dtypes(include=[float, int]).columns
data_filled = data_cleaned[numeric_columns].fillna(data_cleaned[numeric_columns].mean())

# 使用數值型特徵進行 PCA（3 個主成分）
pca = PCA(n_components=3)
pca_result = pca.fit_transform(data_filled)

# 創建一個 DataFrame，包含 PCA 結果和原始索引及類別
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2', 'PC3'])
pca_df['Index'] = pca_df.index  # 加入每個樣本的索引
pca_df['Structure'] = data['Structure']  # 假設你根據 'Structure' 列進行上色

# 使用 Plotly 繪製 3D 散點圖，懸停顯示索引，點的大小設置較小
fig = px.scatter_3d(pca_df, x='PC1', y='PC2', z='PC3', color='Structure', hover_data=['Index'], size_max=5)
fig.update_traces(marker=dict(size=2))  # 調整點的大小

# 顯示圖形
fig.show()
