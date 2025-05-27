# 含氟分子生成器项目

[![Python Tests](https://github.com/yiozeeyoung/fluorinated-molecule-generator/actions/workflows/python-test.yml/badge.svg)](https://github.com/yiozeeyoung/fluorinated-molecule-generator/actions/workflows/python-test.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Python](https://img.shields.io/badge/python-3.9%2B-blue)
[![RDKit](https://img.shields.io/badge/RDKit-2024.03.5-orange)](https://www.rdkit.org/)

这是一个用于生成大规模含氟分子库的高性能分子生成器项目，支持量子化学计算和高级分子分析。通过多进程并行处理，能够高效生成、分析和可视化数万级分子。

## 项目特性

- **高性能生成**：支持多进程并行生成数万级分子
- **量子化学计算**：集成DFT理论计算和分子性质预测
- **智能优化**：基于能量和稳定性的分子优化
- **可视化分析**：丰富的分子结构和性质可视化工具
- **数据库管理**：高效的分子数据存储和检索

## 文件结构

### 主要生成器
- `mol_generator_high_performance.py` - 高性能分子生成器（推荐）
- `mol_generator_quantum_enhanced.py` - 量子增强生成器
- `mol_generator_advanced.py` - 高级生成器
- `mol_genrator_fixed.py` - 修复版生成器
- `mol_genrator.py` - 原始版本

### 可视化工具
- `visualize_molecules.py` - 分子可视化工具
- `batch_processor.py` - 批处理和性能分析工具

### 数据文件
- `top_10000_energy_molecules.csv` - 按能量排序的万级分子数据
- `advanced_fluorinated_molecules_results.csv` - 高级生成结果
- `quantum_enhanced_fluorinated_molecules.csv` - 量子增强结果
- `fluorinated_molecules.db` - 分子数据库

### 配置文件
- `quick_start.py` - 快速启动脚本
- `requirements.txt` - 依赖包列表

## 快速开始

### 1. 安装依赖
```bash
pip install -r requirements.txt
```

### 2. 快速生成分子
```bash
python quick_start.py
```

### 3. 高性能生成
```python
from mol_generator_high_performance import HighPerformanceMoleculeLibraryGenerator

generator = HighPerformanceMoleculeLibraryGenerator(max_workers=8)
report = generator.generate_large_scale_library(
    num_molecules=10000,
    min_stability=30.0,
    max_energy=60.0
)
```

### 4. 可视化分子
```python
from visualize_molecules import MoleculeVisualizer

visualizer = MoleculeVisualizer('top_10000_energy_molecules.csv')
visualizer.run_visualization(n_molecules=20)
```

## 主要功能

### 分子生成
- 基于SMILES的智能分子生成
- 多种氟化策略（单氟、多氟、环状氟化）
- 分子模板驱动生成
- 实时能量优化

### 量子化学计算
- DFT能量计算
- 分子轨道分析
- 电子密度计算
- 振动频率分析

### 性质预测
- 分子量、LogP、极性表面积
- 氢键供体/受体数量
- 旋转键数量
- 稳定性评分

### 可视化分析
- 分子结构网格图
- 性质分布图表
- 能量-稳定性相关性
- 批量分子比较

## 性能特点

- **生成速度**：可达1000+分子/秒
- **并行处理**：支持多核并行计算
- **内存优化**：流式处理大规模数据
- **容错机制**：自动处理生成失败的分子

## 应用领域

- 药物分子设计
- 氟化学研究
- 材料科学
- 分子模拟
- 化学信息学

## 开发团队

含氟分子生成器开发组

## 版本历史

- v1.0 - 基础分子生成功能
- v2.0 - 添加高性能并行处理
- v3.0 - 集成量子化学计算
- v4.0 - 完善可视化和分析工具

## 许可证

MIT License