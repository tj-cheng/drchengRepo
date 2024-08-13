#!/bin/bash

# 假設您的資料文件名為 data.txt
data_file="data.txt"

# 使用 awk 計算機率對應的數量總和
probability_sum=$(awk -F ',' '{ sum[$3] += $2 } END { for (p in sum) print p, sum[p] }' "$data_file")

# 使用 awk 計算密度對應的數量總和
density_sum=$(awk -F ',' '{ sum[$4] += $2 } END { for (d in sum) print d, sum[d] }' "$data_file")

# 輸出結果
echo "機率對應的數量總和："
echo "$probability_sum"

echo "密度對應的數量總和："
echo "$density_sum"

