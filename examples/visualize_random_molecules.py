# 添加用于随机可视化20个分子的示例脚本

创建一个简单的示例脚本，展示如何随机选择并可视化分子：

<function_calls>
<invoke name="create_file">
<parameter name="content">#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
随机可视化分子示例脚本
从生成的分子库中选择20个分子进行可视化
"""

import os
import sys

# 添加主目录到Python路径，确保可以导入utils模块
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.visualize_molecules import MoleculeVisualizer

def main():
    """
    主函数 - 随机选择20个分子并可视化
    """
    # 默认使用top_10000_energy_molecules.csv数据文件
    data_file = os.path.join('data', 'top_10000_energy_molecules.csv')
    
    # 确保文件存在
    if not os.path.exists(data_file):
        print(f"错误：找不到数据文件 {data_file}")
        print("请确保运行目录正确或指定正确的文件路径")
        return
        
    print(f"使用数据文件: {data_file}")
    
    # 创建可视化器实例
    visualizer = MoleculeVisualizer(data_file)
    
    # 运行可视化，随机选择20个分子
    visualizer.run_visualization(n_molecules=20, save_results=True)
    
    print("可视化完成!")

if __name__ == "__main__":
    main()
