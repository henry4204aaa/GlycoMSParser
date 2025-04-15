import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, classification_report
from xgboost import XGBClassifier
import joblib
from tqdm import tqdm
import numpy as np
import warnings

# 抑制所有警告訊息
warnings.filterwarnings('ignore')

# 加載數據
file_path = 'All_data.csv'
data = pd.read_csv(file_path)

# 為每個樣本自動生成唯一的 ID
data['ID'] = data.index + 1  # 使用行索引來生成唯一 ID，從1開始

# 過濾樣本數少於5的類別
def filter_data_for_training(data, min_samples=5):
    class_counts = data['Structure'].value_counts()
    valid_classes = class_counts[class_counts >= min_samples].index
    filtered_data = data[data['Structure'].isin(valid_classes)].copy()
    return filtered_data

# 過濾資料
training_data = filter_data_for_training(data)

# 將類別標籤轉換為數字
label_encoder = LabelEncoder()
training_data['Structure'] = label_encoder.fit_transform(training_data['Structure'])

# 分離特徵和標籤
X = training_data.drop(columns=['Structure', 'ID'])  # 特徵，去掉 ID 和標籤列
y = training_data['Structure']  # 標籤
ids = training_data['ID']  # 样本 ID

# 訓練過程的進度條
print("開始訓練 XGBoost 模型...")

# 定義 XGBoost 模型與 OneVsRestClassifier
clf_xgb = OneVsRestClassifier(XGBClassifier(n_estimators=15, eval_metric='mlogloss'))

# 定義交叉驗證
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 準備記錄錯誤樣本和所有樣本的列表
errors = []
results = []

# 使用交叉驗證預測，並顯示進度條
y_pred = cross_val_predict(clf_xgb, X, y, cv=skf, method='predict', verbose=1)
y_proba = cross_val_predict(clf_xgb, X, y, cv=skf, method='predict_proba', verbose=1)  # 獲取預測概率

# 記錄所有樣本的信息，包括信心分數
for i, (actual, predicted) in enumerate(zip(y, y_pred)):
    confidence_score = np.max(y_proba[i])
    sample_info = {
        'ID': ids.iloc[i],
        'Actual': label_encoder.inverse_transform([actual])[0],
        'Predicted': label_encoder.inverse_transform([predicted])[0],
        'Confidence_Score': confidence_score
    }
    if actual != predicted:
        errors.append(sample_info)  # 記錄錯誤樣本
    results.append(sample_info)  # 記錄所有樣本

# 計算混淆矩陣
conf_matrix = confusion_matrix(y, y_pred)
print("混淆矩陣:\n", conf_matrix)

# 計算分類報告，包括準確率和F1分數（單獨和加權）
classification_report_str = classification_report(y, y_pred, target_names=label_encoder.classes_, output_dict=False)
classification_report_dict = classification_report(y, y_pred, target_names=label_encoder.classes_, output_dict=True)
print("分類報告:\n", classification_report_str)

# 將錯誤樣本和所有樣本分別保存為 CSV 文件
errors_df = pd.DataFrame(errors)
results_df = pd.DataFrame(results)

# 儲存錯誤樣本並將其置於文件的開頭
errors_df.to_csv('XGbOVRCV_model_errors.csv', index=False)
results_df.to_csv('XGbOVRCV_model_all_samples.csv', index=False)

print("所有樣本與錯誤樣本的結果已保存至 'XGbOVRCV_model_all_samples.csv' 和 'XGbOVRCV_model_errors.csv'")

# 保存模型
clf_xgb.fit(X, y)
joblib.dump(clf_xgb, 'XGbOVRCV.joblib')
print("模型已保存至 'XGbOVRCV.joblib'")
'''使用 tqdm 或其他工具來顯示訓練過程的進度
刪除最小樣本5
OneVsRestClassifier＆cross-validation& RF訓練n_estimators=15
RF對於每一筆資料的信心程度ac &f1 score(每個類別單獨計算分數和加權平均的 分數)需要一格混淆矩陣呈現
模型保存
在每次交叉驗證後，記錄哪些樣本被模型錯誤分類，以及這些錯誤分類的具體情況。這可以通過以下幾步實現：

記錄錯誤樣本：在每次 fit 和 predict 之後，對比預測結果和實際標籤，並記錄那些被錯誤分類的樣本。這些錯誤樣本包括它們的 ID、實際標籤、預測結果以及預測的置信度。

輸出錯誤樣本信息：將這些錯誤樣本的信息保存到一個 CSV 文件中，這樣你可以在訓練後檢查這些樣本，並分析為何模型無法正確分類它們。
'''

'''
On 2024/10/05
% pwd
/Users/henry/henry_old/Desktop/python/massAI/2024project
(MS) henry@HenrydeMacBook-Pro 2024project % python XGbOVRCV.py
開始訓練 XGBoost 模型...
混淆矩陣:
 [[10  0  0  0  0  0  0  0  0  0  0  0  0  1]
 [ 1  9  0  0  0  0  0  0  0  0  0  0  0  0]
 [ 0  1  7  0  0  0  0  0  0  0  0  0  0  0]
 [ 0  0  0  5  0  0  0  0  1  0  0  0  0  0]
 [ 0  0  0  0 12  2  0  0  0  0  0  0  0  0]
 [ 0  0  1  0  0 18  0  0  0  0  0  0  0  0]
 [ 0  0  0  2  0  0  6  0  0  0  0  0  0  0]
 [ 0  0  0  0  1  0  0  9  1  0  0  0  0  0]
 [ 0  0  0  0  0  0  0  2  5  0  0  0  0  0]
 [ 0  0  0  0  0  0  0  0  0 10  2  0  0  0]
 [ 0  0  0  0  0  0  1  0  0  0  9  0  0  0]
 [ 0  0  0  0  0  0  0  0  0  1  0 19  0  0]
 [ 0  0  0  0  0  0  0  0  0  0  0  0  9  0]
 [ 0  0  0  0  0  0  0  0  0  0  0  0  0  8]]
分類報告:
               precision    recall  f1-score   support

      F1H3N5       0.91      0.91      0.91        11
      F1H4N5       0.90      0.90      0.90        10
    F1H5N4S1       0.88      0.88      0.88         8
    F1H5N5K2       0.71      0.83      0.77         6
    F1H5N5S1       0.92      0.86      0.89        14
    F1H6N4S1       0.90      0.95      0.92        19
    F1H6N4S2       0.86      0.75      0.80         8
    F2H5N5K1       0.82      0.82      0.82        11
    F2H6N4S1       0.71      0.71      0.71         7
    F2H6N5S1       0.91      0.83      0.87        12
    F2H7N4S1       0.82      0.90      0.86        10
    F2H7N4S2       1.00      0.95      0.97        20
        H3N4       1.00      1.00      1.00         9
        H3N5       0.89      1.00      0.94         8

    accuracy                           0.89       153
   macro avg       0.87      0.88      0.87       153
weighted avg       0.89      0.89      0.89       153

所有樣本與錯誤樣本的結果已保存至 'XGbOVRCV_model_all_samples.csv' 和 'XGbOVRCV_model_errors.csv'
模型已保存至 'XGbOVRCV.joblib'
'''