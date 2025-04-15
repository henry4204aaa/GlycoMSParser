import pandas as pd
from tqdm import tqdm
import numpy as np
import os

# 處理 'peak_d' 列的函數
def process_peak_d_column(df):
    # 展開 'peak_d' 字典
    peak_d_expanded = df['peak_d'].apply(eval).apply(pd.Series)
    
    # 填補 NaN 為 1.0
    peak_d_expanded_filled = peak_d_expanded.fillna(1.0)
    
    # 只對原有數值進行 log10(x + 1) 轉換，補上的 1.0 保持不變
    peak_d_expanded_log = peak_d_expanded_filled.apply(lambda x: np.where(x == 1.0, 1.0, np.log10(x + 1)))
    
    # 合併處理好的數據與其他欄位
    processed_df = pd.concat([df.drop(columns=['peak_d']), peak_d_expanded_log], axis=1)
    
    return processed_df

# 分批讀取 Excel 並保存為 CSV 的函數，並展示進度條
def xlsx_to_csv_in_chunks(file_path, output_csv_path, chunk_size=1000):
    # 獲取文件總行數，僅讀取標題行來避免內存占用
    total_rows = pd.read_excel(file_path, sheet_name=0).shape[0]

    # 初始化進度條
    pbar = tqdm(total=total_rows, desc="Processing rows")

    start_row = 0

    # 嘗試讀取中間結果文件來繼續進行
    if os.path.exists(output_csv_path):
        try:
            processed_data = pd.read_csv(output_csv_path)
            start_row = len(processed_data)  # 繼續從已處理的行數開始
            pbar.update(start_row)
        except pd.errors.ParserError:
            print("中間文件讀取失敗，可能存在結構不一致的問題。將重新從頭開始處理。")
            os.remove(output_csv_path)  # 刪除損壞的中間結果文件
            start_row = 0

    # 分塊讀取 Excel 文件
    for row in range(start_row, total_rows, chunk_size):
        # 使用 skiprows 和 nrows 分批讀取
        chunk = pd.read_excel(file_path, skiprows=range(1, row + 1), nrows=chunk_size)

        # 處理 'peak_d' 列
        chunk_processed = process_peak_d_column(chunk)

        # 保存中間結果
        if row == 0 and start_row == 0:
            chunk_processed.to_csv(output_csv_path, index=False, mode='w')
        else:
            chunk_processed.to_csv(output_csv_path, index=False, header=False, mode='a')

        # 更新進度條
        pbar.update(min(chunk_size, total_rows - row))

    pbar.close()

# 定義路徑與檔案
file_path = 'full_dataset.xlsx'
output_csv_path = 'smallCandyset.csv'

# 分批轉換並處理所有數據
xlsx_to_csv_in_chunks(file_path, output_csv_path)
