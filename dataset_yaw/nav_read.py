import pandas as pd

# 读取 .nav 文件，使用空格作为分隔符
df = pd.read_csv('KF_GINS_Navresult.nav', delimiter=r'\s+', engine='python')

# 打印数据框的列数和列名以确认
print(f"列数: {df.shape[1]}")
print(f"列名: {df.columns.tolist()}")

# 确认数据框有足够的列
if df.shape[1] >= 4:
    # 提取第三列和第四列数据
    lat_lon_df = df.iloc[:, [1,2,3]]
    lat_lon_df.columns = ['time_stamp','lat', 'lon']

    # 将提取的数据保存到新的 csv 文件中
    lat_lon_df.to_csv('lat_lon.csv', index=False)

    # 打印前几行数据以确认
    print(lat_lon_df.head())
else:
    print("数据框的列数不足，无法提取第三列和第四列数据。")