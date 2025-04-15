import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
import os
import gc

# 步驟 1：逐步將 Excel 轉換為 CSV 文件（結合方案）
def excel_to_csv_in_chunks(file_path, output_csv_path, chunk_size=1000):
    total_rows = pd.read_excel(file_path, sheet_name=0, engine='openpyxl').shape[0]
    pbar = tqdm(total=total_rows, desc="Converting Excel to CSV")
    start_row = 0

    # 嘗試讀取中間結果文件來繼續進行
    if os.path.exists(output_csv_path):
        try:
            processed_data = pd.read_csv(output_csv_path)
            start_row = len(processed_data)
            pbar.update(start_row)
        except pd.errors.ParserError:
            print("中間文件讀取失敗，重新從頭開始處理。")
            os.remove(output_csv_path)
            start_row = 0

    # 分塊讀取 Excel 文件並逐步保存
    with pd.ExcelFile(file_path, engine='openpyxl') as xls:
        for i in range(start_row, total_rows, chunk_size):
            chunk = pd.read_excel(xls, sheet_name=0, skiprows=range(1, i + 1), nrows=chunk_size, header=0 if i == 0 else None)
            if i == 0 and start_row == 0:
                chunk.to_csv(output_csv_path, index=False, mode='w')
            else:
                chunk.to_csv(output_csv_path, index=False, header=False, mode='a')
            pbar.update(len(chunk))
    pbar.close()

# 步驟 2：使用 Pandas 分塊處理 CSV 文件
def process_csv_in_chunks(input_csv_path, output_csv_path, chunk_size=1000):
    # 第一步：收集所有可能的 key
    if 'all_keys.pkl' in os.listdir():
        with open('all_keys.pkl', 'rb') as f:
            all_keys = pickle.load(f)
    else:
        all_keys = set()
        for chunk in pd.read_csv(input_csv_path, chunksize=chunk_size):
            for peak_d in chunk['peak_d']:
                try:
                    peak_dict = eval(peak_d)
                    all_keys.update(peak_dict.keys())
                except:
                    continue
        all_keys = list(all_keys)
        with open('all_keys.pkl', 'wb') as f:
            pickle.dump(all_keys, f)

    # 第二步：分塊處理 'peak_d' 列並保存結果
    with tqdm(total=os.path.getsize(input_csv_path), desc="Processing chunks") as pbar:
        for chunk in pd.read_csv(input_csv_path, chunksize=chunk_size):
            try:
                # 展開 'peak_d' 列
                peak_d_expanded = chunk['peak_d'].apply(eval).apply(pd.Series)
                peak_d_expanded = peak_d_expanded.reindex(columns=all_keys).fillna(1.0)
                peak_d_expanded = peak_d_expanded.applymap(lambda x: 1.0 if x == 1.0 else np.log10(x + 1))
                processed_chunk = pd.concat([chunk.drop(columns=['peak_d']), peak_d_expanded], axis=1)

                # 寫入 CSV 文件（逐步寫入）
                if not os.path.exists(output_csv_path):
                    processed_chunk.to_csv(output_csv_path, index=False, mode='w')
                else:
                    processed_chunk.to_csv(output_csv_path, index=False, header=False, mode='a')
            except Exception as e:
                print(f"Error processing chunk: {e}")
            finally:
                # 顯式刪除對象並進行垃圾回收
                del chunk, peak_d_expanded, processed_chunk
                gc.collect()
                pbar.update(chunk_size)

# 定義路徑與檔案
excel_file_path = 'full_dataset.xlsx'
temp_csv_path = 'full_dataset.csv'

# 如果 full_dataset.csv 存在，跳過步驟 1
if not os.path.exists(temp_csv_path):
    # 執行步驟 1：逐步將 Excel 轉換為 CSV
    excel_to_csv_in_chunks(excel_file_path, temp_csv_path)
final_output_csv_path = 'processed_largeCandyset.csv'

# 執行步驟 1：逐步將 Excel 轉換為 CSV
excel_to_csv_in_chunks(excel_file_path, temp_csv_path)

# 執行步驟 2：使用 Pandas 分塊處理 CSV 文件
process_csv_in_chunks(temp_csv_path, final_output_csv_path)