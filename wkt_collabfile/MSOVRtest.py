import pandas as pd
import joblib
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score

# 加載CSV數據
def load_data(file_path):
    return pd.read_csv(file_path)

# 加載測試數據
input_data = 'testtestdata.csv'  # 測試集檔案路徑
test_data = load_data(input_data)

# 使用與訓練相同的 LabelEncoder 來處理測試集
label_encoder = LabelEncoder()

# 假設測試集中的結構標籤也是文字類型
test_data['Structure'] = label_encoder.fit_transform(test_data['Structure'])

# 分離測試集的特徵和標籤
X_test = test_data.drop(columns=['Structure'])  # 特徵
y_test = test_data['Structure']  # 標籤

# 加載隨機森林模型
clf_rf = joblib.load('OVRrandom_forest_model.joblib')

# 使用隨機森林模型進行預測
y_pred_rf = clf_rf.predict(X_test)

# 評估隨機森林模型的性能
accuracy_rf = accuracy_score(y_test, y_pred_rf)
print("隨機森林模型的準確率:", accuracy_rf)

# 加載 XGBoost 模型
clf_xgb = joblib.load('OVRxgboost_model.joblib')

# 使用 XGBoost 模型進行預測
y_pred_xgb = clf_xgb.predict(X_test)

# 評估 XGBoost 模型的性能
accuracy_xgb = accuracy_score(y_test, y_pred_xgb)
print("XGBoost 模型的準確率:", accuracy_xgb)


'''
import pandas as pd
import joblib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
import umap  # 如果需要降維處理

# 加載測試數據
def load_data(file_path):
    return pd.read_csv(file_path)

# 加載未見過的測試數據
test_data = load_data('new_test_data.csv')  # 假設這是未來的測試數據

# 確保有 ID 列來跟踪樣本
if 'ID' not in test_data.columns:
    test_data['ID'] = test_data.index

# 將特徵與 ID 分開
X_test = test_data.drop(columns=['ID'])  # 特徵
ids = test_data['ID']  # 样本 ID

# 加載已訓練的模型
clf_rf = joblib.load('OVRrandom_forest_model.joblib')

# 使用模型進行預測和獲取預測概率
y_pred_rf = clf_rf.predict(X_test)
y_proba_rf = clf_rf.predict_proba(X_test)  # 獲取每個樣本的預測概率

# 構建結果表格，將每個樣本的 ID、預測類別及其概率顯示出來
result_df = pd.DataFrame({
    'ID': ids,
    'Predicted_Class': y_pred_rf,
    'Confidence_Score': np.max(y_proba_rf, axis=1)  # 最大概率作為可信度
})

# 保存結果
result_df.to_csv('prediction_results_with_confidence_test.csv', index=False)
print("測試結果與可信度已保存至 'prediction_results_with_confidence_test.csv'")

# 使用 UMAP 進行降維並可視化樣本的相似性
reducer = umap.UMAP()
embedding = reducer.fit_transform(y_proba_rf)  # 使用預測的概率作為樣本間的距離依據

# 可視化
plt.figure(figsize=(10, 7))
sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=y_pred_rf, palette='Set2', legend='full')
plt.title("UMAP: 样本的二維分群圖")
plt.show()

'''