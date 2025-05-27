#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
生成高性能分子库示例
使用高性能生成器生成1000个含氟分子
"""

import os
import sys

# 添加主目录到Python路径
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.mol_generator_high_performance import HighPerformanceMoleculeLibraryGenerator

def main():
    """
    主函数 - 生成1000个含氟分子
    """
    print("初始化高性能分子生成器...")
    
    # 创建生成器实例，使用4个工作进程
    generator = HighPerformanceMoleculeLibraryGenerator(max_workers=4)
    
    print(f"开始生成1000个含氟分子 (使用 {generator.max_workers} 个并行进程)")
    
    # 生成分子库
    report = generator.generate_large_scale_library(
        num_molecules=1000,  # 生成1000个分子
        min_stability=30.0,  # 最小稳定性阈值
        max_energy=60.0,     # 最大能量阈值
        export_path="results/generated_molecules_example.csv"  # 导出路径
    )
    
    # 打印生成报告
    print("\n=== 生成报告 ===")
    print(f"总生成分子数: {report['total_molecules']}")
    print(f"成功生成分子数: {report['successful_molecules']}")
    print(f"生成速度: {report['molecules_per_second']:.2f} 分子/秒")
    print(f"成功率: {report['success_rate']*100:.2f}%")
    print(f"平均能量: {report['average_energy']:.2f} kcal/mol")
    print(f"平均稳定性: {report['average_stability']:.2f}")
    print(f"使用内存: {report['memory_usage']:.2f} MB")
    print(f"结果已保存到: {report['export_path']}")
    
    print("\n生成完成!")

if __name__ == "__main__":
    main()
