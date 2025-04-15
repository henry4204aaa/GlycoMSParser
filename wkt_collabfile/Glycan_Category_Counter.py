import pandas as pd

def count_glycan_categories(file_path):
    # 讀取 Excel 文件
    df = pd.read_excel(file_path)
    
    # 確認 'glycan' 欄位存在
    if 'glycan' not in df.columns:
        raise ValueError("Excel 文件中不存在名為 'glycan' 的欄位")
    
    # 使用 value_counts() 統計每個類別的數量
    glycan_counts = df['glycan'].value_counts()
    
    # 輸出結果
    print(glycan_counts)
    
    # 如果需要保存統計結果為 CSV 文件
    glycan_counts.to_csv('glycan_category_counts.csv', header=['Count'])

# 指定 Excel 文件的路徑
file_path = 'full_dataset.xlsx'

# 執行函數進行統計
count_glycan_categories(file_path)
