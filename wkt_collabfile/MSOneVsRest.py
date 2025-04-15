'''
這是RandomForest＆XGBoost基本版
使用 tqdm 或其他工具來顯示訓練過程的進度
刪除最小樣本5
OneVsRestClassifier＆cross-validation& RF訓練n_estimators=15
RF對於每一筆資料的信心程度與ac &f1 score(每個類別單獨計算分數和加權平均的 分數)
模型保存

OneVsRestClassifier＆cross-validation&XGBoost訓練n_estimators=15
XGBoost對於每一筆資料的信心程度與ac &f1 score(每個類別單獨計算分數和加權平均的 分數)
模型保存'''

import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.multiclass import OneVsRestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import make_scorer, f1_score, accuracy_score, classification_report
import joblib
from tqdm import tqdm  # 用來顯示進度條

# 加載CSV數據
def load_data(file_path):
    return pd.read_csv(file_path)

# 過濾數據函數，排除掉樣本數少於5的類別
def filter_data_for_training(data, min_samples=5):
    class_counts = data['Structure'].value_counts()
    valid_classes = class_counts[class_counts >= min_samples].index
    filtered_data = data[data['Structure'].isin(valid_classes)].copy()
    return filtered_data

# 加載資料
input_data = 'testtraindata.csv'
data = load_data(input_data)

# 過濾資料，排除少於5個樣本的類別
training_data = filter_data_for_training(data)

# 將類別標籤轉換為數字
label_encoder = LabelEncoder()
training_data['Structure'] = label_encoder.fit_transform(training_data['Structure'])

# 分離特徵和標籤
X = training_data.drop(columns=['Structure'])  # 特徵
y = training_data['Structure']  # 標籤

# 定義交叉驗證
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 定義評估指標
scoring = {
    'accuracy': make_scorer(accuracy_score),
    'f1_weighted': make_scorer(f1_score, average='weighted')
}

# 訓練過程的進度條
print("開始訓練隨機森林模型...")
clf_rf = OneVsRestClassifier(RandomForestClassifier(n_estimators=15))

# 訓練與交叉驗證
accuracy_rf = []
f1_weighted_rf = []

for train_index, test_index in tqdm(skf.split(X, y), total=5):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]
    
    clf_rf.fit(X_train, y_train)
    y_pred = clf_rf.predict(X_test)
    accuracy_rf.append(accuracy_score(y_test, y_pred))
    f1_weighted_rf.append(f1_score(y_test, y_pred, average='weighted'))

# 計算加權平均的準確率與 F1 分數
print(f"隨機森林模型加權平均準確率: {sum(accuracy_rf) / len(accuracy_rf):.4f}")
print(f"隨機森林模型加權平均 F1 分數: {sum(f1_weighted_rf) / len(f1_weighted_rf):.4f}")

# 保存隨機森林模型
joblib.dump(clf_rf, 'OVRrandom_forest_model.joblib')

# 訓練 XGBoost 模型
print("開始訓練 XGBoost 模型...")
clf_xgb = OneVsRestClassifier(XGBClassifier(n_estimators=15, use_label_encoder=False, eval_metric='mlogloss'))

# 訓練與交叉驗證
accuracy_xgb = []
f1_weighted_xgb = []

for train_index, test_index in tqdm(skf.split(X, y), total=5):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]
    
    clf_xgb.fit(X_train, y_train)
    y_pred = clf_xgb.predict(X_test)
    accuracy_xgb.append(accuracy_score(y_test, y_pred))
    f1_weighted_xgb.append(f1_score(y_test, y_pred, average='weighted'))

# 計算加權平均的準確率與 F1 分數
print(f"XGBoost 模型加權平均準確率: {sum(accuracy_xgb) / len(accuracy_xgb):.4f}")
print(f"XGBoost 模型加權平均 F1 分數: {sum(f1_weighted_xgb) / len(f1_weighted_xgb):.4f}")

# 保存 XGBoost 模型
joblib.dump(clf_xgb, 'OVRxgboost_model.joblib')




