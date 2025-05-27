"""
高性能分子生成器 - 批处理工具和性能优化脚本
=============================================

提供批量生成、性能监控和结果分析的辅助工具
"""

import os
import sys
import time
import argparse
import json
import multiprocessing as mp
from typing import Dict, List
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mol_generator_high_performance import (
    HighPerformanceMoleculeLibraryGenerator,
    MoleculeDatabase,
    ProgressMonitor
)

class BatchProcessor:
    """批处理器"""
    
    def __init__(self, config_file: str = None):
        self.config = self._load_config(config_file)
    
    def _load_config(self, config_file: str) -> Dict:
        """加载配置文件"""
        default_config = {
            "batch_sizes": [1000, 5000, 10000, 20000],
            "min_stability": 40.0,
            "max_energy": 50.0,
            "max_workers": min(mp.cpu_count(), 32),
            "output_dir": "batch_results"
        }
        
        if config_file and os.path.exists(config_file):
            with open(config_file, 'r', encoding='utf-8') as f:
                user_config = json.load(f)
                default_config.update(user_config)
        
        return default_config
    
    def run_batch_generation(self):
        """运行批量生成"""
        os.makedirs(self.config["output_dir"], exist_ok=True)
        
        results = []
        
        for batch_size in self.config["batch_sizes"]:
            print(f"\n{'='*60}")
            print(f"开始生成 {batch_size:,} 个分子的批次")
            print(f"{'='*60}")
            
            # 创建专用数据库
            db_path = os.path.join(self.config["output_dir"], f"molecules_{batch_size}.db")
            generator = HighPerformanceMoleculeLibraryGenerator(
                max_workers=self.config["max_workers"]
            )
            generator.db = MoleculeDatabase(db_path)
            
            # 生成分子
            start_time = time.time()
            report = generator.generate_large_scale_library(
                num_molecules=batch_size,
                min_stability=self.config["min_stability"],
                max_energy=self.config["max_energy"]
            )
            generation_time = time.time() - start_time
            
            # 分析结果
            analysis = generator.analyze_large_library()
            
            # 导出结果
            export_file = os.path.join(
                self.config["output_dir"], 
                f"molecules_{batch_size}.csv"
            )
            generator.export_top_molecules(min(batch_size, 5000), 'energy')
            
            # 记录结果
            batch_result = {
                'batch_size': batch_size,
                'generation_time': generation_time,
                'molecules_per_second': report['molecules_per_second'],
                'success_rate': report['success_rate'],
                'total_molecules': analysis['total_molecules'],
                'avg_energy': analysis['energy_stats']['mean'],
                'avg_stability': analysis['stability_stats']['mean'],
                'convergence_rate': analysis['convergence_rate']
            }
            results.append(batch_result)
            
            print(f"批次完成：{generation_time:.1f}秒，"
                  f"速度：{report['molecules_per_second']:.1f} mol/s")
        
        # 保存批处理报告
        self._save_batch_report(results)
        self._plot_performance_charts(results)
        
        return results
    
    def _save_batch_report(self, results: List[Dict]):
        """保存批处理报告"""
        report_file = os.path.join(self.config["output_dir"], "batch_report.json")
        
        report = {
            'config': self.config,
            'results': results,
            'summary': {
                'total_batches': len(results),
                'total_molecules': sum(r['total_molecules'] for r in results),
                'avg_speed': np.mean([r['molecules_per_second'] for r in results]),
                'avg_success_rate': np.mean([r['success_rate'] for r in results])
            }
        }
        
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        
        print(f"批处理报告已保存到: {report_file}")
    
    def _plot_performance_charts(self, results: List[Dict]):
        """绘制性能图表"""
        df = pd.DataFrame(results)
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 生成速度 vs 批量大小
        axes[0, 0].plot(df['batch_size'], df['molecules_per_second'], 'b-o')
        axes[0, 0].set_xlabel('批量大小')
        axes[0, 0].set_ylabel('生成速度 (分子/秒)')
        axes[0, 0].set_title('生成速度 vs 批量大小')
        axes[0, 0].grid(True, alpha=0.3)
        
        # 成功率 vs 批量大小
        axes[0, 1].plot(df['batch_size'], df['success_rate']*100, 'g-o')
        axes[0, 1].set_xlabel('批量大小')
        axes[0, 1].set_ylabel('成功率 (%)')
        axes[0, 1].set_title('成功率 vs 批量大小')
        axes[0, 1].grid(True, alpha=0.3)
        
        # 平均能量分布
        axes[1, 0].bar(range(len(df)), df['avg_energy'], color='orange', alpha=0.7)
        axes[1, 0].set_xlabel('批次')
        axes[1, 0].set_ylabel('平均能量 (kcal/mol)')
        axes[1, 0].set_title('各批次平均能量')
        axes[1, 0].set_xticks(range(len(df)))
        axes[1, 0].set_xticklabels([f"{size//1000}k" for size in df['batch_size']])
        
        # 平均稳定性分布
        axes[1, 1].bar(range(len(df)), df['avg_stability'], color='purple', alpha=0.7)
        axes[1, 1].set_xlabel('批次')
        axes[1, 1].set_ylabel('平均稳定性评分')
        axes[1, 1].set_title('各批次平均稳定性')
        axes[1, 1].set_xticks(range(len(df)))
        axes[1, 1].set_xticklabels([f"{size//1000}k" for size in df['batch_size']])
        
        plt.tight_layout()
        chart_file = os.path.join(self.config["output_dir"], "performance_charts.png")
        plt.savefig(chart_file, dpi=300, bbox_inches='tight')
        plt.show()
        
        print(f"性能图表已保存到: {chart_file}")

class PerformanceProfiler:
    """性能分析器"""
    
    def __init__(self):
        self.results = []
    
    def profile_generation_methods(self):
        """分析不同生成方法的性能"""
        print("性能分析：不同生成方法对比")
        print("="*50)
        
        # 测试不同的worker数量
        worker_counts = [1, 2, 4, 8, min(16, mp.cpu_count())]
        test_molecules = 1000
        
        for workers in worker_counts:
            print(f"\n测试 {workers} 个工作进程...")
            
            generator = HighPerformanceMoleculeLibraryGenerator(max_workers=workers)
            
            start_time = time.time()
            report = generator.generate_large_scale_library(
                num_molecules=test_molecules,
                min_stability=30.0,
                max_energy=60.0
            )
            elapsed = time.time() - start_time
            
            result = {
                'workers': workers,
                'molecules': test_molecules,
                'time': elapsed,
                'speed': report['molecules_per_second'],
                'success_rate': report['success_rate'],
                'efficiency': report['molecules_per_second'] / workers
            }
            self.results.append(result)
            
            print(f"  耗时: {elapsed:.1f}s")
            print(f"  速度: {report['molecules_per_second']:.1f} mol/s")
            print(f"  效率: {result['efficiency']:.1f} mol/s/worker")
        
        self._plot_worker_performance()
        return self.results
    
    def _plot_worker_performance(self):
        """绘制工作进程性能图"""
        df = pd.DataFrame(self.results)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # 绝对速度
        ax1.plot(df['workers'], df['speed'], 'b-o', label='生成速度')
        ax1.set_xlabel('工作进程数')
        ax1.set_ylabel('生成速度 (分子/秒)')
        ax1.set_title('生成速度 vs 进程数')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # 效率（每进程速度）
        ax2.plot(df['workers'], df['efficiency'], 'r-o', label='每进程效率')
        ax2.set_xlabel('工作进程数')
        ax2.set_ylabel('效率 (分子/秒/进程)')
        ax2.set_title('进程效率 vs 进程数')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig('worker_performance.png', dpi=300, bbox_inches='tight')
        plt.show()

class MoleculeAnalyzer:
    """分子分析器"""
    
    def __init__(self, db_path: str):
        self.db = MoleculeDatabase(db_path)
    
    def comprehensive_analysis(self, sample_size: int = 5000):
        """全面分析分子库"""
        print(f"正在进行分子库全面分析（采样{sample_size:,}个分子）...")
        
        import sqlite3
        with sqlite3.connect(self.db.db_path) as conn:
            # 获取统计数据
            query = f"""
                SELECT 
                    COUNT(*) as total,
                    AVG(energy) as avg_energy,
                    MIN(energy) as min_energy,
                    MAX(energy) as max_energy,
                    AVG(stability_score) as avg_stability,
                    AVG(molecular_weight) as avg_mw,
                    AVG(fluorine_count) as avg_fluorines,
                    SUM(optimization_converged) * 100.0 / COUNT(*) as convergence_rate
                FROM molecules
            """
            stats = conn.execute(query).fetchone()
            
            # 采样详细分析
            sample_query = f"""
                SELECT * FROM molecules 
                ORDER BY RANDOM() 
                LIMIT {sample_size}
            """
            df = pd.read_sql_query(sample_query, conn)
        
        # 打印统计信息
        print(f"\n{'='*60}")
        print(f"分子库统计分析")
        print(f"{'='*60}")
        print(f"总分子数: {stats[0]:,}")
        print(f"平均能量: {stats[1]:.2f} kcal/mol")
        print(f"能量范围: {stats[2]:.2f} - {stats[3]:.2f} kcal/mol")
        print(f"平均稳定性: {stats[4]:.2f}")
        print(f"平均分子量: {stats[5]:.2f}")
        print(f"平均氟原子数: {stats[6]:.1f}")
        print(f"收敛率: {stats[7]:.1f}%")
        
        # 生成详细图表
        self._plot_detailed_analysis(df)
        
        return {
            'total_molecules': stats[0],
            'energy_stats': {'avg': stats[1], 'min': stats[2], 'max': stats[3]},
            'avg_stability': stats[4],
            'avg_molecular_weight': stats[5],
            'avg_fluorines': stats[6],
            'convergence_rate': stats[7]
        }
    
    def _plot_detailed_analysis(self, df: pd.DataFrame):
        """绘制详细分析图表"""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        # 能量分布
        axes[0, 0].hist(df['energy'], bins=50, alpha=0.7, edgecolor='black')
        axes[0, 0].set_xlabel('能量 (kcal/mol)')
        axes[0, 0].set_ylabel('频次')
        axes[0, 0].set_title('能量分布')
        axes[0, 0].grid(True, alpha=0.3)
        
        # 稳定性分布
        axes[0, 1].hist(df['stability_score'], bins=30, alpha=0.7, edgecolor='black', color='green')
        axes[0, 1].set_xlabel('稳定性评分')
        axes[0, 1].set_ylabel('频次')
        axes[0, 1].set_title('稳定性分布')
        axes[0, 1].grid(True, alpha=0.3)
        
        # 分子量分布
        axes[0, 2].hist(df['molecular_weight'], bins=40, alpha=0.7, edgecolor='black', color='orange')
        axes[0, 2].set_xlabel('分子量')
        axes[0, 2].set_ylabel('频次')
        axes[0, 2].set_title('分子量分布')
        axes[0, 2].grid(True, alpha=0.3)
        
        # 氟原子数分布
        fluorine_counts = df['fluorine_count'].value_counts().sort_index()
        axes[1, 0].bar(fluorine_counts.index, fluorine_counts.values, alpha=0.7, color='purple')
        axes[1, 0].set_xlabel('氟原子数')
        axes[1, 0].set_ylabel('分子数量')
        axes[1, 0].set_title('氟化模式分布')
        axes[1, 0].grid(True, alpha=0.3)
        
        # 能量 vs 稳定性散点图
        axes[1, 1].scatter(df['energy'], df['stability_score'], alpha=0.5, s=10)
        axes[1, 1].set_xlabel('能量 (kcal/mol)')
        axes[1, 1].set_ylabel('稳定性评分')
        axes[1, 1].set_title('能量 vs 稳定性')
        axes[1, 1].grid(True, alpha=0.3)
        
        # LogP 分布
        axes[1, 2].hist(df['logp'], bins=30, alpha=0.7, edgecolor='black', color='red')
        axes[1, 2].set_xlabel('LogP')
        axes[1, 2].set_ylabel('频次')
        axes[1, 2].set_title('LogP 分布')
        axes[1, 2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('detailed_molecular_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()

def main():
    """主程序"""
    parser = argparse.ArgumentParser(description='高性能分子生成器批处理工具')
    parser.add_argument('--mode', choices=['batch', 'profile', 'analyze'], 
                       default='batch', help='运行模式')
    parser.add_argument('--config', help='配置文件路径')
    parser.add_argument('--db', help='分析的数据库路径')
    parser.add_argument('--molecules', type=int, default=10000, help='生成分子数量')
    
    args = parser.parse_args()
    
    if args.mode == 'batch':
        print("运行批处理模式...")
        processor = BatchProcessor(args.config)
        results = processor.run_batch_generation()
        print(f"批处理完成，共处理 {len(results)} 个批次")
        
    elif args.mode == 'profile':
        print("运行性能分析模式...")
        profiler = PerformanceProfiler()
        results = profiler.profile_generation_methods()
        print("性能分析完成")
        
    elif args.mode == 'analyze':
        if not args.db:
            args.db = "fluorinated_molecules.db"
        
        if not os.path.exists(args.db):
            print(f"数据库文件不存在: {args.db}")
            return
        
        print(f"分析数据库: {args.db}")
        analyzer = MoleculeAnalyzer(args.db)
        results = analyzer.comprehensive_analysis()
        print("分析完成")

if __name__ == "__main__":
    main()
