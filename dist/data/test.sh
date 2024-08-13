#!/bin/bash

# 計算需要的元素數量
num_elements=11  # 由於包含 0 和 1，所以有 11 個元素

# 初始化一維矩陣
rgprob=()
dd=()

# 生成矩陣元素
for ((i=0; i<num_elements; i++)); 
do
    value=$(echo "scale=1; $i * 0.1" | bc)  # 使用 bc 進行浮點運算
    rgprob+=($value)
    value1=$(echo "scale=1; 1+$i * 0.1" | bc)  # 使用 bc 進行浮點運算
    dd+=($value1)
done

# 輸出結果
echo "${rgprob[@]}"
echo "${dd[@]}"

