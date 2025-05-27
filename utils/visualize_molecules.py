#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
含氟分子可视化工具
从生成的分子数据库中随机选择分子进行可视化显示
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
import random
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

class MoleculeVisualizer:
    """分子可视化器"""
    
    def __init__(self, csv_file_path):
        """
        初始化可视化器
        
        Args:
            csv_file_path (str): 分子数据CSV文件路径
        """
        self.csv_file_path = csv_file_path
        self.molecules_df = None
        self.load_data()
    
    def load_data(self):
        """加载分子数据"""
        try:
            self.molecules_df = pd.read_csv(self.csv_file_path)
            print(f"成功加载 {len(self.molecules_df)} 个分子数据")
            print(f"数据列: {list(self.molecules_df.columns)}")
        except Exception as e:
            print(f"加载数据失败: {e}")
            return False
        return True
    
    def validate_smiles(self, smiles):
        """验证SMILES字符串是否有效"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except:
            return False
    
    def get_random_molecules(self, n=20):
        """
        随机选择n个分子
        
        Args:
            n (int): 选择的分子数量
            
        Returns:
            pd.DataFrame: 随机选择的分子数据
        """
        if self.molecules_df is None or len(self.molecules_df) == 0:
            print("没有可用的分子数据")
            return None
        
        # 过滤有效的SMILES
        valid_molecules = []
        for idx, row in self.molecules_df.iterrows():
            if self.validate_smiles(row['smiles']):
                valid_molecules.append(row)
        
        if len(valid_molecules) == 0:
            print("没有找到有效的分子")
            return None
        
        # 随机选择
        n = min(n, len(valid_molecules))
        selected_indices = random.sample(range(len(valid_molecules)), n)
        selected_molecules = [valid_molecules[i] for i in selected_indices]
        
        return pd.DataFrame(selected_molecules)
    
    def draw_molecule_grid(self, molecules_df, save_path=None):
        """
        绘制分子网格图
        
        Args:
            molecules_df (pd.DataFrame): 分子数据
            save_path (str): 保存路径
        """
        if molecules_df is None or len(molecules_df) == 0:
            print("没有分子数据可供可视化")
            return
        
        n_molecules = len(molecules_df)
        cols = 4  # 每行4个分子
        rows = (n_molecules + cols - 1) // cols  # 向上取整
        
        fig, axes = plt.subplots(rows, cols, figsize=(16, 4*rows))
        fig.suptitle(f'随机选择的 {n_molecules} 个含氟分子', fontsize=16, fontweight='bold')
        
        # 确保axes是二维数组
        if rows == 1:
            axes = axes.reshape(1, -1)
        elif cols == 1:
            axes = axes.reshape(-1, 1)
        
        molecules = []
        legends = []
        
        for idx, (_, row) in enumerate(molecules_df.iterrows()):
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is not None:
                molecules.append(mol)
                # 创建分子信息标签
                legend = (f"ID: {row['id']}\n"
                         f"SMILES: {row['smiles']}\n"
                         f"能量: {row['energy']:.2f}\n"
                         f"稳定性: {row['stability_score']}\n"
                         f"分子量: {row['molecular_weight']:.1f}\n"
                         f"F原子数: {row['fluorine_count']}")
                legends.append(legend)
        
        # 使用RDKit绘制分子
        if molecules:
            img = Draw.MolsToGridImage(
                molecules,
                molsPerRow=cols,
                subImgSize=(300, 300),
                legends=legends,
                useSVG=False
            )
            
            # 显示图像
            plt.figure(figsize=(20, 5*rows))
            plt.imshow(img)
            plt.axis('off')
            plt.title(f'随机选择的 {len(molecules)} 个含氟分子结构', fontsize=16, pad=20)
            
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches='tight', 
                           facecolor='white', edgecolor='none')
                print(f"分子结构图已保存到: {save_path}")
            
            plt.show()
    
    def create_detailed_visualization(self, molecules_df, save_dir=None):
        """
        创建详细的分子可视化，包括单独的分子图和属性图
        
        Args:
            molecules_df (pd.DataFrame): 分子数据
            save_dir (str): 保存目录
        """
        if molecules_df is None or len(molecules_df) == 0:
            print("没有分子数据可供可视化")
            return
        
        # 创建保存目录
        if save_dir and not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        # 1. 绘制分子属性分布图
        self.plot_molecular_properties(molecules_df, save_dir)
        
        # 2. 绘制个别分子的详细结构
        self.plot_individual_molecules(molecules_df, save_dir)
        
        # 3. 创建分子属性对比表
        self.create_properties_table(molecules_df, save_dir)
    
    def plot_molecular_properties(self, molecules_df, save_dir=None):
        """绘制分子属性分布图"""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('选定含氟分子的属性分布', fontsize=16, fontweight='bold')
        
        # 能量分布
        axes[0, 0].hist(molecules_df['energy'], bins=10, alpha=0.7, color='skyblue', edgecolor='black')
        axes[0, 0].set_title('能量分布')
        axes[0, 0].set_xlabel('能量')
        axes[0, 0].set_ylabel('频次')
        
        # 稳定性评分分布
        axes[0, 1].hist(molecules_df['stability_score'], bins=10, alpha=0.7, color='lightgreen', edgecolor='black')
        axes[0, 1].set_title('稳定性评分分布')
        axes[0, 1].set_xlabel('稳定性评分')
        axes[0, 1].set_ylabel('频次')
        
        # 分子量分布
        axes[0, 2].hist(molecules_df['molecular_weight'], bins=10, alpha=0.7, color='lightcoral', edgecolor='black')
        axes[0, 2].set_title('分子量分布')
        axes[0, 2].set_xlabel('分子量')
        axes[0, 2].set_ylabel('频次')
        
        # 氟原子数分布
        axes[1, 0].hist(molecules_df['fluorine_count'], bins=range(1, int(molecules_df['fluorine_count'].max())+2), 
                       alpha=0.7, color='gold', edgecolor='black')
        axes[1, 0].set_title('氟原子数分布')
        axes[1, 0].set_xlabel('氟原子数')
        axes[1, 0].set_ylabel('频次')
        
        # logP分布
        if 'logp' in molecules_df.columns:
            axes[1, 1].hist(molecules_df['logp'], bins=10, alpha=0.7, color='plum', edgecolor='black')
            axes[1, 1].set_title('logP分布')
            axes[1, 1].set_xlabel('logP')
            axes[1, 1].set_ylabel('频次')
        
        # 能量vs稳定性散点图
        axes[1, 2].scatter(molecules_df['energy'], molecules_df['stability_score'], 
                          alpha=0.7, c=molecules_df['fluorine_count'], cmap='viridis')
        axes[1, 2].set_title('能量 vs 稳定性')
        axes[1, 2].set_xlabel('能量')
        axes[1, 2].set_ylabel('稳定性评分')
        cbar = plt.colorbar(axes[1, 2].collections[0], ax=axes[1, 2])
        cbar.set_label('氟原子数')
        
        plt.tight_layout()
        
        if save_dir:
            save_path = os.path.join(save_dir, 'molecular_properties_distribution.png')
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"属性分布图已保存到: {save_path}")
        
        plt.show()
    
    def plot_individual_molecules(self, molecules_df, save_dir=None):
        """绘制个别分子的详细结构"""
        # 选择前6个分子进行详细展示
        top_molecules = molecules_df.head(6)
        
        for idx, (_, row) in enumerate(top_molecules.iterrows()):
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is not None:
                # 创建分子图像
                drawer = rdMolDraw2D.MolDraw2DCairo(400, 400)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                
                # 保存单个分子图像
                if save_dir:
                    img_path = os.path.join(save_dir, f"molecule_{row['id']}.png")
                    with open(img_path, 'wb') as f:
                        f.write(drawer.GetDrawingText())
                    print(f"分子 {row['id']} 的结构图已保存到: {img_path}")
    
    def create_properties_table(self, molecules_df, save_dir=None):
        """创建分子属性对比表"""
        # 选择关键属性
        key_columns = ['id', 'smiles', 'energy', 'stability_score', 
                      'molecular_weight', 'fluorine_count']
        
        if 'logp' in molecules_df.columns:
            key_columns.append('logp')
        
        table_data = molecules_df[key_columns].copy()
        
        # 数值格式化
        for col in ['energy', 'molecular_weight']:
            if col in table_data.columns:
                table_data[col] = table_data[col].round(2)
        
        if 'logp' in table_data.columns:
            table_data['logp'] = table_data['logp'].round(3)
        
        # 保存为CSV
        if save_dir:
            table_path = os.path.join(save_dir, 'selected_molecules_properties.csv')
            table_data.to_csv(table_path, index=False, encoding='utf-8-sig')
            print(f"分子属性表已保存到: {table_path}")
        
        # 显示表格
        print("\n=== 选定分子的属性对比 ===")
        print(table_data.to_string(index=False))
        
        return table_data
    
    def run_visualization(self, n_molecules=20, save_results=True):
        """
        运行完整的可视化流程
        
        Args:
            n_molecules (int): 选择的分子数量
            save_results (bool): 是否保存结果
        """
        print(f"开始可视化 {n_molecules} 个随机选择的含氟分子...")
        
        # 随机选择分子
        selected_molecules = self.get_random_molecules(n_molecules)
        
        if selected_molecules is None:
            print("无法获取有效的分子数据")
            return
        
        print(f"成功选择了 {len(selected_molecules)} 个分子")
        
        # 创建保存目录
        save_dir = None
        if save_results:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            save_dir = f"molecule_visualization_{timestamp}"
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            print(f"结果将保存到目录: {save_dir}")
        
        # 绘制分子网格图
        grid_save_path = os.path.join(save_dir, 'molecules_grid.png') if save_dir else None
        self.draw_molecule_grid(selected_molecules, grid_save_path)
        
        # 创建详细可视化
        self.create_detailed_visualization(selected_molecules, save_dir)
        
        print("可视化完成！")
        
        return selected_molecules


def main():
    """主函数"""
    # 设置文件路径
    csv_file_path = "top_10000_energy_molecules.csv"
    
    # 检查文件是否存在
    if not os.path.exists(csv_file_path):
        print(f"错误: 找不到文件 {csv_file_path}")
        print("请确保分子数据文件在当前目录中")
        return
    
    # 创建可视化器
    visualizer = MoleculeVisualizer(csv_file_path)
    
    # 运行可视化
    try:
        selected_molecules = visualizer.run_visualization(n_molecules=20, save_results=True)
        
        if selected_molecules is not None:
            print(f"\n=== 统计信息 ===")
            print(f"总共可视化分子数: {len(selected_molecules)}")
            print(f"平均能量: {selected_molecules['energy'].mean():.2f}")
            print(f"平均稳定性: {selected_molecules['stability_score'].mean():.1f}")
            print(f"平均分子量: {selected_molecules['molecular_weight'].mean():.1f}")
            print(f"氟原子数范围: {selected_molecules['fluorine_count'].min()} - {selected_molecules['fluorine_count'].max()}")
        
    except Exception as e:
        print(f"可视化过程中出现错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
