#!/bin/bash
#SBATCH -N 1
#SBATCH -p c01
#SBATCH -n 6
#SBATCH -o extract.log
# 设置副本索引范围，例如从0到29（总共30个副本）
start_replica=0
end_replica=29

# 遍历所有副本索引
for (( i=$start_replica; i<=$end_replica; i++ ))
do
  # 格式化索引以确保两位数（如果有需要的话）
  replica_index=$(printf "%02d" $i)

  # 提取副本i的轨迹并保存到指定目录和文件名
  extract_trajectory follow_structure --replica $i walker/trajectory.${replica_index}.pdb
  # 打印进度信息
  # 打印进度信息
  echo "follow_structure for replica $i"
done