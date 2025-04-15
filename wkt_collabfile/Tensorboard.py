import tensorflow as tf
from tensorboard.plugins import projector
import pandas as pd
import os

# 加載數據
#file_path = 'All_data.csv'
file_path = 'testtraindata.csv'
data = pd.read_csv(file_path)

# 創建一個 log 目錄
log_dir = "logs/feature_projector"
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

# 提取數值型特徵列作為嵌入向量
numeric_columns = data.select_dtypes(include=[float, int]).columns
features_np = data[numeric_columns].values  # 將數值特徵轉換為 numpy 格式

# 創建一個 TensorFlow 變量來保存特徵
feature_var = tf.Variable(features_np, name='feature_embedding')

# 保存變量
checkpoint = tf.train.Checkpoint(feature_embedding=feature_var)
checkpoint.save(os.path.join(log_dir, "feature_embedding.ckpt"))

# 創建 Projector 配置文件
config = projector.ProjectorConfig()
embedding = config.embeddings.add()
embedding.tensor_name = "feature_embedding/.ATTRIBUTES/VARIABLE_VALUE"

# 使用 'Structure' 列作為元數據標籤
metadata_path = os.path.join(log_dir, 'metadata.tsv')
with open(metadata_path, 'w') as f:
    for structure in data['Structure']:
        f.write(f"{structure}\n")  # 將 'Structure' 信息寫入元數據文件

# 指定 metadata 文件的路徑
embedding.metadata_path = 'metadata.tsv'

# 保存配置
projector.visualize_embeddings(log_dir, config)

print(f"Embedding Projector 文件已保存至: {log_dir}")
