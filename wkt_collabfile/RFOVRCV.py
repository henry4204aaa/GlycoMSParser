import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, classification_report
import joblib
from tqdm import tqdm
import numpy as np



mode=1
if mode==1:
    #file_path = 'NGdata1014.csv'
    file_path = 'Odata1014.csv'
    print('this is result from :',file_path)
    column_to_drop=['Structure', 'ID','Source','IUPACname(optional)','Glycanannotation2','unique_ID']
    column_of_class='Structure'
elif mode==2:
    file_path = 'smallCandyset.csv'
    column_to_drop=['RT', 'filename', 'GlycoPost_ID']
    column_of_class='glycan'
# 加載數據

data = pd.read_csv(file_path)





# 過濾樣本數少於5的類別
def filter_data_for_training(data, column_of_class, min_samples=5):
    # 計算每個類別的樣本數
    class_counts = data[column_of_class].value_counts()
    print("過濾前的類別及其樣本數：")
    for cls, count in class_counts[class_counts >= 0].items():
        print(f"類別: {cls}, 樣本數: {count}")
    print('---------------------------------')
    # 過濾樣本數少於指定數目的類別
    valid_classes = class_counts[class_counts >= min_samples].index
    
    # 過濾後的資料
    filtered_data = data[data[column_of_class].isin(valid_classes)].copy()

    # 印出每個類別的名稱和樣本數
    print("過濾後的類別及其樣本數：")
    for cls, count in class_counts[class_counts >= min_samples].items():
        print(f"類別: {cls}, 樣本數: {count}")
    
    return filtered_data

# 使用資料進行過濾，並指定類別欄位
training_data = filter_data_for_training(data, column_of_class)
print('--------a-------')



# 將類別標籤轉換為數字
label_encoder = LabelEncoder()
training_data[column_of_class] = label_encoder.fit_transform(training_data[column_of_class])

# 分離特徵和標籤
X = training_data.drop(columns=[column_of_class,]+column_to_drop)  # 特徵，去掉 ID 和標籤列
y = training_data[column_of_class]  # 標籤
ids = training_data['ID']  # 样本 ID

# 訓練過程的進度條
print("開始訓練隨機森林模型...")

# 定義隨機森林模型與 OneVsRestClassifier
clf_rf = OneVsRestClassifier(RandomForestClassifier(n_estimators=25, random_state=42))
#clf_rf = RandomForestClassifier(n_estimators=25, random_state=42)
# 定義交叉驗證
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 準備記錄錯誤樣本和所有樣本的列表
errors = []
results = []

# 使用交叉驗證預測，並顯示進度條
y_pred = cross_val_predict(clf_rf, X, y, cv=skf, method='predict', verbose=1)
y_proba = cross_val_predict(clf_rf, X, y, cv=skf, method='predict_proba', verbose=1)  # 獲取預測概率

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

if mode==1:
    # 儲存錯誤樣本並將其置於文件的開頭
    errors_df.to_csv('RFOVRCV_model_errors.csv', index=False)
    results_df.to_csv('RFOVRCV_model_all_samples.csv', index=False)

    print("所有樣本與錯誤樣本的結果已保存至 'RFOVRCV_model_all_samples.csv' 和 'RFOVRCV_model_errors.csv'")

    # 保存模型
    clf_rf.fit(X, y)
    joblib.dump(clf_rf, 'RFOVRCV.joblib')
    print("模型已保存至 'RFOVRCV.joblib'")

elif mode==2:
    # 儲存錯誤樣本並將其置於文件的開頭
    errors_df.to_csv('RFOVRCV_model_errors_CC.csv', index=False)
    results_df.to_csv('RFOVRCV_model_CC_all_samples.csv', index=False)

    print("所有樣本與錯誤樣本的結果已保存至 'RFOVRCV_model_CC_all_samples.csv' 和 'RFOVRCV_model_CC_errors.csv'")

    # 保存模型
    clf_rf.fit(X, y)
    joblib.dump(clf_rf, 'RFOVRCV_CC.joblib')
    print("模型已保存至 'RFOVRCV_CC.joblib'")

'''使用 tqdm 或其他工具來顯示訓練過程的進度
刪除最小樣本5
OneVsRestClassifier＆cross-validation& RF訓練n_estimators=15
RF對於每一筆資料的信心程度ac &f1 score(每個類別單獨計算分數和加權平均的 分數)需要一格混淆矩陣呈現
模型保存
在每次交叉驗證後，記錄哪些樣本被模型錯誤分類，以及這些錯誤分類的具體情況。這可以通過以下幾步實現：

記錄錯誤樣本：在每次 fit 和 predict 之後，對比預測結果和實際標籤，並記錄那些被錯誤分類的樣本。這些錯誤樣本包括它們的 ID、實際標籤、預測結果以及預測的置信度。

輸出錯誤樣本信息：將這些錯誤樣本的信息保存到一個 CSV 文件中，這樣你可以在訓練後檢查這些樣本，並分析為何模型無法正確分類它們。
'''
#要加一個drop class 功能 
'''
On 2024/10/05
% pwd
/Users/henry/henry_old/Desktop/python/massAI/2024project
(MS) henry@HenrydeMacBook-Pro 2024project % python RFOVRCV.py 
開始訓練隨機森林模型...
混淆矩陣:
 [[11  0  0  0  0  0  0  0  0  0  0  0  0  0]
 [ 0 10  0  0  0  0  0  0  0  0  0  0  0  0]
 [ 0  0  8  0  0  0  0  0  0  0  0  0  0  0]
 [ 0  0  0  6  0  0  0  0  0  0  0  0  0  0]
 [ 0  0  0  0 14  0  0  0  0  0  0  0  0  0]
 [ 0  0  0  0  0 19  0  0  0  0  0  0  0  0]
 [ 0  0  0  0  0  0  8  0  0  0  0  0  0  0]
 [ 0  0  0  0  0  0  0 10  1  0  0  0  0  0]
 [ 0  0  0  0  0  0  0  3  4  0  0  0  0  0]
 [ 0  0  0  0  0  0  0  0  0 12  0  0  0  0]
 [ 0  0  0  0  0  0  0  0  0  0 10  0  0  0]
 [ 0  0  0  0  0  0  0  0  0  0  0 20  0  0]
 [ 0  0  0  0  0  0  0  0  0  0  0  0  9  0]
 [ 0  0  0  0  0  0  0  0  0  0  0  0  0  8]]
分類報告:
1
所有樣本與錯誤樣本的結果已保存至 'RFOVRCV_model_all_samples.csv' 和 'RFOVRCV_model_errors.csv'
模型已保存至 'RFOVRCV.joblib'
(MS) henry@HenrydeMBP 2024project % python RFOVRCV.py
過濾前的類別及其樣本數：
類別: H1N2G1F1, 樣本數: 39
類別: H2N2G1, 樣本數: 39
類別: H1N1S2, 樣本數: 39
類別: H1N1S1G1, 樣本數: 33
類別: H1N2S1F1, 樣本數: 31
類別: H1N1S1, 樣本數: 30
類別: H1N2F1, 樣本數: 24
類別: H2N2, 樣本數: 20
類別: H2N1G1, 樣本數: 18
類別: H2N2S1, 樣本數: 17
類別: H1N1, 樣本數: 17
類別: H1N2G1, 樣本數: 16
類別: H1N1G1, 樣本數: 15
類別: H1N1S3, 樣本數: 12
類別: H1N1KDN1, 樣本數: 12
類別: H3N2S1F1, 樣本數: 10
類別: H1N2KDN1, 樣本數: 10
類別: N1S1, 樣本數: 10
類別: H1N1S1F1, 樣本數: 10
類別: H1N1G2, 樣本數: 9
類別: H1N1G1F1, 樣本數: 9
類別: H2N1S2, 樣本數: 9
類別: N2KDN1, 樣本數: 9
類別: H2N3, 樣本數: 8
類別: H1N2S1, 樣本數: 8
類別: N1G1, 樣本數: 8
類別: N2S1, 樣本數: 8
類別: H1N1S4, 樣本數: 7
類別: H3N2F1, 樣本數: 7
類別: H2N2F1, 樣本數: 6
類別: H1N3F1, 樣本數: 5
類別: H1N1KDN1S1, 樣本數: 5
類別: H1N2KDN1F1, 樣本數: 5
類別: N1KDN1, 樣本數: 4
類別: N2, 樣本數: 4
類別: H1N1KDN2, 樣本數: 4
類別: H3N2, 樣本數: 4
類別: H2N2S2, 樣本數: 4
類別: H1N1S2G1, 樣本數: 3
類別: H1N2F1S1G1, 樣本數: 3
類別: H3N2G1F1, 樣本數: 2
類別: H1N1F1, 樣本數: 2
類別: N4F1, 樣本數: 1
類別: H2N1S1, 樣本數: 1
---------------------------------
過濾後的類別及其樣本數：
類別: H1N2G1F1, 樣本數: 39
類別: H2N2G1, 樣本數: 39
類別: H1N1S2, 樣本數: 39
類別: H1N1S1G1, 樣本數: 33
類別: H1N2S1F1, 樣本數: 31
類別: H1N1S1, 樣本數: 30
類別: H1N2F1, 樣本數: 24
類別: H2N2, 樣本數: 20
類別: H2N1G1, 樣本數: 18
類別: H2N2S1, 樣本數: 17
類別: H1N1, 樣本數: 17
類別: H1N2G1, 樣本數: 16
類別: H1N1G1, 樣本數: 15
類別: H1N1S3, 樣本數: 12
類別: H1N1KDN1, 樣本數: 12
類別: H3N2S1F1, 樣本數: 10
類別: H1N2KDN1, 樣本數: 10
類別: N1S1, 樣本數: 10
類別: H1N1S1F1, 樣本數: 10
類別: H1N1G2, 樣本數: 9
類別: H1N1G1F1, 樣本數: 9
類別: H2N1S2, 樣本數: 9
類別: N2KDN1, 樣本數: 9
類別: H2N3, 樣本數: 8
類別: H1N2S1, 樣本數: 8
類別: N1G1, 樣本數: 8
類別: N2S1, 樣本數: 8
類別: H1N1S4, 樣本數: 7
類別: H3N2F1, 樣本數: 7
類別: H2N2F1, 樣本數: 6
類別: H1N3F1, 樣本數: 5
類別: H1N1KDN1S1, 樣本數: 5
類別: H1N2KDN1F1, 樣本數: 5
--------a-------
開始訓練隨機森林模型...
混淆矩陣:
 [[17  0  0 ...  0  0  0]
 [ 0 13  0 ...  0  0  0]
 [ 0  0  6 ...  0  0  0]
 ...
 [ 0  0  0 ... 10  0  0]
 [ 0  0  0 ...  0  8  1]
 [ 0  0  0 ...  0  0  8]]
分類報告:
               precision    recall  f1-score   support

        H1N1       0.89      1.00      0.94        17
      H1N1G1       1.00      0.87      0.93        15
    H1N1G1F1       1.00      0.67      0.80         9
      H1N1G2       1.00      0.78      0.88         9
    H1N1KDN1       1.00      0.92      0.96        12
  H1N1KDN1S1       1.00      0.80      0.89         5
      H1N1S1       0.93      0.93      0.93        30
    H1N1S1F1       1.00      0.90      0.95        10
    H1N1S1G1       0.92      1.00      0.96        33
      H1N1S2       0.95      0.95      0.95        39
      H1N1S3       1.00      0.92      0.96        12
      H1N1S4       0.88      1.00      0.93         7
      H1N2F1       0.96      1.00      0.98        24
      H1N2G1       1.00      1.00      1.00        16
    H1N2G1F1       0.84      0.92      0.88        39
    H1N2KDN1       0.91      1.00      0.95        10
  H1N2KDN1F1       1.00      0.60      0.75         5
      H1N2S1       1.00      0.75      0.86         8
    H1N2S1F1       1.00      1.00      1.00        31
      H1N3F1       1.00      1.00      1.00         5
      H2N1G1       0.85      0.94      0.89        18
      H2N1S2       0.82      1.00      0.90         9
        H2N2       0.82      0.90      0.86        20
      H2N2F1       1.00      0.50      0.67         6
      H2N2G1       0.90      0.95      0.93        39
      H2N2S1       0.79      0.65      0.71        17
        H2N3       0.75      0.75      0.75         8
      H3N2F1       1.00      0.71      0.83         7
    H3N2S1F1       0.90      0.90      0.90        10
        N1G1       1.00      1.00      1.00         8
        N1S1       0.91      1.00      0.95        10
      N2KDN1       1.00      0.89      0.94         9
        N2S1       0.73      1.00      0.84         8

    accuracy                           0.92       505
   macro avg       0.93      0.88      0.90       505
weighted avg       0.92      0.92      0.91       505



過濾前的類別及其樣本數：
類別: F2H7N4S2, 樣本數: 45
類別: F1H6N4S1, 樣本數: 34
類別: F2H7N4S1, 樣本數: 33
類別: F1H5N4S1, 樣本數: 30
類別: F1H3N5, 樣本數: 30
類別: H3N5, 樣本數: 28
類別: F1H5N4S2, 樣本數: 28
類別: F1H6N4S2, 樣本數: 27
類別: F1H4N5, 樣本數: 27
類別: F2H6N4S1, 樣本數: 25
類別: F1H5N5S1, 樣本數: 24
類別: H3N4, 樣本數: 21
類別: H5N4G1, 樣本數: 18
類別: F2H6N5S1, 樣本數: 17
類別: H3N6G2, 樣本數: 17
類別: Man5, 樣本數: 17
類別: H3N6S1G1, 樣本數: 14
類別: F1H5N4, 樣本數: 14
類別: H5N4S1, 樣本數: 13
類別: F3H7N5, 樣本數: 13
類別: F1H5N5, 樣本數: 13
類別: H4N4G1, 樣本數: 13
類別: F2H7N4G2, 樣本數: 13
類別: F1H6N4, 樣本數: 13
類別: F2H5N5K1, 樣本數: 11
類別: F1H7N4S2, 樣本數: 11
類別: F1H3N4, 樣本數: 11
類別: Man6, 樣本數: 10
類別: F1H4N4, 樣本數: 10
類別: Man8, 樣本數: 10
類別: H5N4, 樣本數: 10
類別: Man7, 樣本數: 10
類別: Man9, 樣本數: 10
類別: F2H5N5, 樣本數: 9
類別: F1H6N5S1, 樣本數: 7
類別: F1H7N4S1, 樣本數: 6
類別: H5N3, 樣本數: 6
類別: F2H7N4S1G1, 樣本數: 6
類別: H3N6S2, 樣本數: 6
類別: F1H5N5K2, 樣本數: 6
類別: H4N4S1, 樣本數: 6
類別: F2H7N4G1, 樣本數: 5
類別: F3H7N4S1, 樣本數: 5
類別: F1H6N4S1G1, 樣本數: 4
類別: H6N4S1, 樣本數: 4
類別: F1H5N5S1K1, 樣本數: 4
類別: F1H4N5S1, 樣本數: 4
類別: F1H4N5K1, 樣本數: 4
類別: H4N4, 樣本數: 4
類別: F3H8N4S1, 樣本數: 3
類別: F2H6N5, 樣本數: 3
類別: H6N3, 樣本數: 3
類別: F2H6N5K1, 樣本數: 2
類別: F1H5N4S1G1, 樣本數: 2
類別: F1H5N4G1, 樣本數: 2
類別: F3H9N5S3, 樣本數: 2
類別: F2H5N4S1G1, 樣本數: 2
類別: H5N4 (Normal), 樣本數: 1
類別: H5N4 (Gal-Gal), 樣本數: 1
類別: F2H6N4G1, 樣本數: 1
類別: F2H7N5, 樣本數: 1
---------------------------------
過濾後的類別及其樣本數：
類別: F2H7N4S2, 樣本數: 45
類別: F1H6N4S1, 樣本數: 34
類別: F2H7N4S1, 樣本數: 33
類別: F1H5N4S1, 樣本數: 30
類別: F1H3N5, 樣本數: 30
類別: H3N5, 樣本數: 28
類別: F1H5N4S2, 樣本數: 28
類別: F1H6N4S2, 樣本數: 27
類別: F1H4N5, 樣本數: 27
類別: F2H6N4S1, 樣本數: 25
類別: F1H5N5S1, 樣本數: 24
類別: H3N4, 樣本數: 21
類別: H5N4G1, 樣本數: 18
類別: F2H6N5S1, 樣本數: 17
類別: H3N6G2, 樣本數: 17
類別: Man5, 樣本數: 17
類別: H3N6S1G1, 樣本數: 14
類別: F1H5N4, 樣本數: 14
類別: H5N4S1, 樣本數: 13
類別: F3H7N5, 樣本數: 13
類別: F1H5N5, 樣本數: 13
類別: H4N4G1, 樣本數: 13
類別: F2H7N4G2, 樣本數: 13
類別: F1H6N4, 樣本數: 13
類別: F2H5N5K1, 樣本數: 11
類別: F1H7N4S2, 樣本數: 11
類別: F1H3N4, 樣本數: 11
類別: Man6, 樣本數: 10
類別: F1H4N4, 樣本數: 10
類別: Man8, 樣本數: 10
類別: H5N4, 樣本數: 10
類別: Man7, 樣本數: 10
類別: Man9, 樣本數: 10
類別: F2H5N5, 樣本數: 9
類別: F1H6N5S1, 樣本數: 7
類別: F1H7N4S1, 樣本數: 6
類別: H5N3, 樣本數: 6
類別: F2H7N4S1G1, 樣本數: 6
類別: H3N6S2, 樣本數: 6
類別: F1H5N5K2, 樣本數: 6
類別: H4N4S1, 樣本數: 6
類別: F2H7N4G1, 樣本數: 5
類別: F3H7N4S1, 樣本數: 5
--------a-------
開始訓練隨機森林模型...
/opt/anaconda3/envs/MS/lib/python3.9/site-packages/sklearn/multiclass.py:551: RuntimeWarning: invalid value encountered in divide
  Y /= np.sum(Y, axis=1)[:, np.newaxis]
混淆矩陣:
 [[11  0  0 ...  0  0  0]
 [ 0 30  0 ...  0  0  0]
 [ 0  0 10 ...  0  0  0]
 ...
 [ 0  0  0 ...  8  1  0]
 [ 0  0  0 ...  0  9  0]
 [ 0  0  0 ...  0  1  9]]
分類報告:
               precision    recall  f1-score   support

      F1H3N4       1.00      1.00      1.00        11
      F1H3N5       1.00      1.00      1.00        30
      F1H4N4       1.00      1.00      1.00        10
      F1H4N5       0.87      0.96      0.91        27
      F1H5N4       1.00      0.93      0.96        14
    F1H5N4S1       0.83      0.97      0.89        30
    F1H5N4S2       0.90      0.96      0.93        28
      F1H5N5       0.89      0.62      0.73        13
    F1H5N5K2       1.00      0.83      0.91         6
    F1H5N5S1       0.96      0.92      0.94        24
      F1H6N4       1.00      0.85      0.92        13
    F1H6N4S1       0.86      0.88      0.87        34
    F1H6N4S2       0.73      0.81      0.77        27
    F1H6N5S1       1.00      0.43      0.60         7
    F1H7N4S1       1.00      0.67      0.80         6
    F1H7N4S2       1.00      0.91      0.95        11
      F2H5N5       0.80      0.89      0.84         9
    F2H5N5K1       0.85      1.00      0.92        11
    F2H6N4S1       0.84      0.84      0.84        25
    F2H6N5S1       0.88      0.88      0.88        17
    F2H7N4G1       0.80      0.80      0.80         5
    F2H7N4G2       1.00      1.00      1.00        13
    F2H7N4S1       0.83      0.88      0.85        33
  F2H7N4S1G1       1.00      0.50      0.67         6
    F2H7N4S2       0.92      1.00      0.96        45
    F3H7N4S1       1.00      0.20      0.33         5
      F3H7N5       0.91      0.77      0.83        13
        H3N4       1.00      0.90      0.95        21
        H3N5       0.90      0.96      0.93        28
      H3N6G2       0.80      0.71      0.75        17
    H3N6S1G1       0.58      0.79      0.67        14
      H3N6S2       0.67      0.67      0.67         6
      H4N4G1       0.92      0.85      0.88        13
      H4N4S1       1.00      0.67      0.80         6
        H5N3       1.00      1.00      1.00         6
        H5N4       0.91      1.00      0.95        10
      H5N4G1       0.89      0.94      0.92        18
      H5N4S1       0.82      0.69      0.75        13
        Man5       0.94      1.00      0.97        17
        Man6       0.91      1.00      0.95        10
        Man7       1.00      0.80      0.89        10
        Man8       0.82      0.90      0.86        10
        Man9       0.82      0.90      0.86        10

    accuracy                           0.89       682
   macro avg       0.90      0.84      0.86       682
weighted avg       0.89      0.89      0.88       682

'''