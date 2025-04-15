import pandas as pd
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="機器學習模型訓練與參數設置")

    # 文件路徑相關
    parser.add_argument('--file-path', type=str, help="要 append 的 CSV 檔案路徑")
    #python MSpredictor.py --file-path '/Users/
    parser.add_argument('--first', type=str, help="第一次創建總資料集")
    
    return parser.parse_args()

def load_data(file_path):
    """
    加載CSV數據
    """
    return pd.read_csv(file_path)

def compare_columns(file_paths):
    """
    比較多個CSV文件的欄位名稱，確保一致
    """
    columns_set = None
    for file_path in file_paths:
        data = load_data(file_path)
        if columns_set is None:
            columns_set = set(data.columns)
            print(f"{file_path} 的欄位: {data.columns.tolist()}")
        else:
            if columns_set != set(data.columns):
                print(f"錯誤: {file_path} 的欄位與其他檔案不一致")
                print(f"{file_path} 的欄位: {data.columns.tolist()}")
                print(f"先前檔案的欄位: {list(columns_set)}")
            else:
                print(f"{file_path} 的欄位一致")

def combine_data(existing_data, new_file_path):
    """
    只將新資料庫與已合併的資料進行合併
    """
    new_data = load_data(new_file_path)
    combined_data = pd.concat([existing_data, new_data], ignore_index=True)
    return combined_data

args = get_args()

if args.first == 'Y':
    # 初次運行：加載並合併初始的資料庫
    file_paths = [
        'zfNGbrain_converted.csv',
        'zfNGovary_converted.csv',
        'zfNGintestine_converted.csv'
    ]

    # 比對欄位名稱
    compare_columns(file_paths)

    # 初始化合併資料
    nulldata = pd.DataFrame()
    for file_path in file_paths:
        all_data = combine_data(nulldata, file_path)

    # 將合併後的所有資料輸出到一個新的 CSV 文件
    all_data.to_csv('All_data.csv', index=False)
    print("初次合併的資料已成功保存到 'All_data.csv'")
else:
    # 當有新的資料庫時，只需要加載新資料並進行合併
    if args.file_path:
        # 加載現有資料
        try:
            All_data = pd.read_csv('All_data.csv')
        except FileNotFoundError:
            print("錯誤：找不到 'All_data.csv'，請先進行初次合併。")
            sys.exit(1)
        compare_columns([All_data,args.file_path])
        # 合併新資料
        All_data = combine_data(All_data, args.file_path)

        # 將合併後的所有資料輸出到一個新的 CSV 文件
        All_data.to_csv('All_data.csv', index=False)
        print("新資料已成功合併並保存到 'All_data.csv'")
    else:
        print('錯誤: 請提供 --file-path 參數')

import pandas as pd
'''
是否要額外加上ＩＤ蘭為？
# 加載數據
data = pd.read_csv('All_data.csv')

# 為每個樣本自動生成唯一的 ID
data['ID'] = data.index + 1  # 使用行索引來生成唯一 ID，從1開始

# 檢查資料集
print(data.head())

# 保存帶有 ID 的資料集
data.to_csv('All_data.csv', index=False)
'''