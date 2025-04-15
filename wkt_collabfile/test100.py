import pandas as pd
import numpy as np

# 第一步：收集前 100 筆資料中的所有可能的 key
def collect_all_keys_sample(file_path, num_rows=100):
    all_keys = set()
    # 使用 Pandas 讀取前 num_rows 行
    df = pd.read_excel(file_path, nrows=num_rows, engine='openpyxl')

    # 掃描所有的字典 key
    for peak_d in df['peak_d']:
        try:
            peak_dict = eval(peak_d)
            all_keys.update(peak_dict.keys())
        except:
            continue
    return list(all_keys)

# 第二步：處理 'peak_d' 列並寫入 CSV
def process_and_write_csv_sample(file_path, output_csv_path, all_keys, num_rows=100):
    # 使用 Pandas 讀取前 num_rows 行
    df = pd.read_excel(file_path, nrows=num_rows, engine='openpyxl')

    # 展開 'peak_d' 列
    peak_d_expanded = df['peak_d'].apply(eval).apply(lambda x: pd.Series(x, index=all_keys))

    # 填補 NaN 為 1.0
    peak_d_expanded_filled = peak_d_expanded.fillna(1.0)

    # 只對原有數值進行 log10(x + 1) 轉換，補上的 1.0 保持不變
    peak_d_expanded_log = peak_d_expanded_filled.apply(lambda x: np.where(x == 1.0, 1.0, np.log10(x + 1)))

    # 合併處理好的數據與其他欄位
    processed_df = pd.concat([df.drop(columns=['peak_d']), pd.DataFrame(peak_d_expanded_log, columns=all_keys)], axis=1)

    # 寫入 CSV 文件
    processed_df.to_csv(output_csv_path, index=False, mode='w')

# 定義路徑與檔案
file_path = 'full_dataset.xlsx'
output_csv_path = 'smallCandy100set.csv'

# 第一步：收集前 100 行的所有 key
all_keys = collect_all_keys_sample(file_path)

# 第二步：處理前 100 行並寫入 CSV
process_and_write_csv_sample(file_path, output_csv_path, all_keys)
