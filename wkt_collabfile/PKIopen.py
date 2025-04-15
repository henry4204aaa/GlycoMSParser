import pickle

# 文件路徑
file_path = 'all_keys.pkl'

# 打開文件並讀取內容
with open(file_path, 'rb') as file:
    data = pickle.load(file)

# 打印內容
print(data)
