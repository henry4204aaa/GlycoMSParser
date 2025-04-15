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
    """加載CSV數據
    """
    return pd.read_csv(file_path)

def add_id_column(data, source_name):
    """為每個原始資料集添加 ID 與源來欄位
    """
    data['ID'] = data.index + 1  # 為每行數據添加 ID
    data['Source'] = source_name  # 添加源來欄位
    return data

def compare_columns(file_paths):
    """比較多個CSV文件的欄位名稱，確保一致
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
    """只將新資料庫與已合併的資料進行合併
    """
    new_data = load_data(new_file_path)
    combined_data = pd.concat([existing_data, new_data], ignore_index=True)
    return combined_data

args = get_args()
#outputpath='Ndata1014'##
outputpath='NGdata1014.csv'##
print('outputpath:',outputpath)
if args.first == 'Y':
    # 初次運行：加載並合併初始的資料庫 python DataCombiner.py --first 'Y'
    #file_paths = [
    #    'zfOGbrainfixed20241012_adduid_dropunknown.csv',
    #    'zfOGintestinefixed20241012_adduid_dropunknown.csv',
    #    'zfOGovaryfixed20241012_adduid_dropunknown.csv'
    #]
    file_paths = [
        'zfNGbrainfixed20241012_adduid_dropunknown.csv',
        'zfNGintestinefixed20241012_adduid_dropunknown.csv',
        'zfNGovaryfixed20241012_adduid_dropunknown.csv'
    ]

    # 比對欄位名稱
    compare_columns(file_paths)

    # 初始化合併資料
    all_data = pd.DataFrame()
    for file_path in file_paths:
        data = load_data(file_path)
        source_name = file_path.split('.')[0]  # 使用檔名作為源來標誌
        data_with_id = add_id_column(data, source_name)
        all_data = pd.concat([all_data, data_with_id], ignore_index=True)

    # 將合併後的所有資料輸出到一個新的 CSV 文件
    all_data.to_csv(outputpath, index=False)
    print("新資料已成功合併並保存到",outputpath)
else:
    # 當有新的資料庫時，只需要加載新資料並進行合併
    if args.file_path:
        # 加載現有資料
        try:
            All_data = pd.read_csv(outputpath)
        except FileNotFoundError:
            print("錯誤：找不到 outputpath，請先進行初次合併。")
            sys.exit(1)
        
        # 比對欄位名稱
        compare_columns([args.file_path, outputpath])
        
        # 加載新資料並添加 ID 與 Source 欄位
        new_data = load_data(args.file_path)
        source_name = args.file_path.split('.')[0]  # 使用檔名作為源來標誌
        new_data_with_id = add_id_column(new_data, source_name)
        
        # 合併新資料
        All_data = pd.concat([All_data, new_data_with_id], ignore_index=True)

        # 將合併後的所有資料輸出到一個新的 CSV 文件
        All_data.to_csv(outputpath, index=False)
        print("新資料已成功合併並保存到",outputpath)
    else:
        print('錯誤: 請提供 --file-path 參數')