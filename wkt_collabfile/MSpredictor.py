'''
python MSpredictor.py --file-path '/Users/henry/henry_old/Desktop/python/massAI/2024project/zfOGbrain_converted.csv' --drop-columns 'IUPACname(optional)' 'Glycanannotation2' --encode-structure

這個命令的功能如下：
--file-path：指定要讀取的 CSV 文件路徑。
--drop-columns：刪除 IUPACname(optional) 和 Glycanannotation2 這兩列。
--encode-structure：啟用對 Structure 列進行編碼。

此外，如果需要自定義其他參數，例如 min-samples 或 cv-folds，你可以在命令中加入這些參數：

python MSpredictor.py --file-path '/Users/henry/henry_old/Desktop/python/massAI/2024project/zfOGbrain_converted.csv' --drop-columns 'IUPACname(optional)' 'Glycanannotation2' --encode-structure --min-samples 5 --cv-folds 3
這樣可以在過濾類別的時候，只保留有 5 個樣本以上的類別，並使用 3 折交叉驗證進行模型訓練。
python MSpredictor.py --file-path 'zfNGintestine_converted.csv' --drop-columns 'IUPACname(optional)' 'Glycanannotation2' --encode-structure --min-samples 5 --cv-folds 3

'''


import pandas as pd
import argparse
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score
from tqdm import tqdm

# 命令列參數設置
def get_args():
    parser = argparse.ArgumentParser(description="機器學習模型訓練與參數設置")

    # 文件路徑相關
    parser.add_argument('--file-path', type=str, required=True, help="要讀取的 CSV 檔案路徑")

    # 資料處理相關
    parser.add_argument('--drop-columns', nargs='+', default=['IUPACname(optional)', 'Glycanannotation2'], 
                        help="要刪除的列名")
    parser.add_argument('--encode-structure', action='store_true', help="是否對 'Structure' 欄位進行編碼")
    parser.add_argument('--min-samples', type=int, default=10, 
                        help="過濾掉樣本數小於 N 的類別 (預設值: 10)")

    # 訓練與模型選擇相關
    parser.add_argument('--test-size', type=float, default=0.3, 
                        help="訓練集和測試集的比例 (預設值: 0.3)")
    parser.add_argument('--cv-folds', type=int, default=5, 
                        help="交叉驗證的分組數量 (預設值: 5)")
    parser.add_argument('--models', nargs='+', default=['SVM', 'RandomForest', 'KNN', 'XGBoost'],
                        help="要訓練的模型 (預設值: ['SVM', 'RandomForest', 'KNN', 'XGBoost'])")
    
    return parser.parse_args()

def main():
    args = get_args()

    # 讀取資料
    data = pd.read_csv(args.file_path)

    # 1. 刪除指定的列
    if args.drop_columns:
        data_cleaned = data.drop(columns=args.drop_columns)
    else:
        data_cleaned = data.copy()

    # 2. 是否對 'Structure' 欄位進行編碼
    if args.encode_structure:
        label_encoder = LabelEncoder()
        data_cleaned['Structure'] = label_encoder.fit_transform(data_cleaned['Structure'])
    
    def show_class_distribution(data, structure_column='Structure'):#要加上row和column數
        """Display the number of samples for each class in the 'Structure' column."""
        class_counts = data[structure_column].value_counts()
        print("Class Distribution:")
        print(class_counts)
        return class_counts

    def get_user_inputs(default_args_min_samples=args.min_samples,default_args_CV_fold=args.cv_folds):
        """Ask the user for 'min-samples' and 'cv-folds' input."""
        userinput=input(f"type y to set min_samples={args.min_samples} and cv_folds={args.cv_folds} or type n to reset:").lower()
        while True:
            if userinput=='y':
                return default_args_min_samples, default_args_CV_fold
            if userinput=='n':
                min_samples = int(input("Enter the minimum number of samples per class (min-samples): "))
                cv_folds = int(input("Enter the number of cross-validation folds (cv-folds): "))
                return min_samples, cv_folds
            else:
                print('Invalid input. Inputmust be y or n')
                userinput=input(f"type y to set min_samples={args.min_samples} and cv_folds={args.cv_folds} or type n to reset:").lower()

    # Display class distribution in the dataset
    class_distribution = show_class_distribution(data)
    min_samples, cv_folds = get_user_inputs()
    # 3. 過濾樣本數小於 N 的類別
    class_counts = data_cleaned['Structure'].value_counts()
    valid_classes = class_counts[class_counts >= args.min_samples].index
    data_cleaned = data_cleaned[data_cleaned['Structure'].isin(valid_classes)]

    # 4. 分割數據集為特徵和標籤
    X = data_cleaned.drop(columns=['Structure'])  # 特徵
    y = data_cleaned['Structure']  # 標籤

    # 5. 分割為訓練集和測試集
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=args.test_size, random_state=42)

    # 6. 檢查並過濾測試集中訓練集不存在的類別
    valid_train_classes = set(y_train)
    X_test = X_test[y_test.isin(valid_train_classes)]
    y_test = y_test[y_test.isin(valid_train_classes)]

    # 構建四個機器學習模型的 pipeline
    pipelines = {
        'SVM': Pipeline([('svm', SVC())]),
        'RandomForest': Pipeline([('rf', RandomForestClassifier())]),
        'KNN': Pipeline([('knn', KNeighborsClassifier())]),
        'XGBoost': Pipeline([('xgb', XGBClassifier(eval_metric='mlogloss'))])
    }

    # 設置各模型的參數範圍
    ##這裡都要改成F1
    param_grids = {#記得重設每個模型的參數 或導出給用戶設定
        'SVM': {'svm__C': [0.1, 1, 10], 'svm__kernel': ['linear', 'rbf']},
        'RandomForest': {'rf__n_estimators': [50, 100, 200], 'rf__max_depth': [None, 10, 20]},
        'KNN': {'knn__n_neighbors': [3, 5, 7]},
        'XGBoost': {'xgb__n_estimators': [50, 100, 200], 'xgb__max_depth': [3, 5, 7]}
    }

    # 儲存最佳模型及其表現
    best_models = {}
    accuracies = {}

    # 使用 tqdm 展示進度條
    for model_name in tqdm(args.models, desc="Training Models"):
        if model_name in pipelines:
            print(f"\nTraining {model_name}...")
            grid_search = GridSearchCV(pipelines[model_name], param_grids[model_name], cv=args.cv_folds)
            grid_search.fit(X_train, y_train)

            # 儲存最佳模型
            best_models[model_name] = grid_search.best_estimator_

            # 預測測試集並計算精度
            y_pred = grid_search.predict(X_test)
            accuracy = accuracy_score(y_test, y_pred)
            accuracies[model_name] = accuracy
            print(f"{model_name} Accuracy: {accuracy}")

    # 輸出所有模型的準確度
    print("\nModel Performance:")
    for model_name, accuracy in accuracies.items():
        print(f"{model_name}: {accuracy}")


if __name__ == "__main__":
    main()


'''1. Transformer
適用性: Transformer 模型最初是為自然語言處理設計的，但它們的多頭自注意力機制可以很好地處理質譜數據，尤其是具有長範圍依賴關係的數據。
優勢: 它能夠動態地關注不同碎片信號的強度、位置和質量比（m/z），這對於聚醣結構的預測是關鍵。利用 Transformer，你可以同時考慮質譜圖中的所有信號，從而學習到更細緻的特徵關係。
變體: 可以考慮使用 BERT 或 GPT 變體，根據預測目標進行微調。此外，Vision Transformer (ViT) 也可以用於處理質譜圖像形式的輸入。
2. Self-Attention Neural Networks (SAN)
適用性: 自注意力網絡可以在質譜數據中找到與碎片化相關的關鍵信號。這類網絡能夠根據數據的整體結構動態地調整各特徵的重要性。
優勢: 與傳統 CNN 不同，Self-Attention 網絡能更好地應對不同碎片之間的關聯性，並識別出哪些信號更具診斷性。
3. Recurrent Neural Networks (RNN) with Attention Mechanism
適用性: 雖然 RNN 已經被部分新技術取代，但當加入注意力機制時，它們依然能在處理序列數據（例如質譜數據）方面表現良好。
優勢: 注意力機制可以讓 RNN 模型學會在序列中不同的時間步上聚焦，識別出具有代表性的碎片和信號強度。
變體: 可以考慮使用 BiLSTM with Attention，這種模型可以從雙向學習質譜中的依賴關係，從而更好地捕捉結構特徵。
4. Graph Neural Networks (GNN) with Attention
適用性: 如果將質譜數據視為一個圖，其中峰值和碎片之間的關係表示邊，那麼 GNN 是一個強大的工具。引入注意力機制（如 Graph Attention Network, GAT）可以讓模型學會在圖中找到最重要的信號。
優勢: 尤其適用於具有複雜結構的數據，例如聚醣，因為它能夠模擬信號之間的結構化關係。
5. Set Transformer
適用性: Set Transformer 專門為處理無序集合設計，質譜數據中的碎片信號可以視為無序集合。該模型利用注意力機制來學習碎片之間的相互關聯。
優勢: 能處理不同數量的質譜信號，並且在無需排序的情況下有效地捕捉整體特徵。
'''