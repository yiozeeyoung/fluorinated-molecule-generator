#!/usr/bin/env python3
"""
高性能含氟分子生成器 - 快速启动脚本
========================================

一键式启动脚本，支持不同规模的分子生成任务
"""

import os
import sys
import json
import argparse
import multiprocessing as mp
from datetime import datetime

def get_system_info():
    """获取系统信息"""
    import psutil
    
    cpu_count = mp.cpu_count()
    memory_gb = psutil.virtual_memory().total / (1024**3)
    
    return {
        'cpu_cores': cpu_count,
        'memory_gb': memory_gb,
        'recommended_workers': min(cpu_count, 32),
        'recommended_batch_size': min(2000, cpu_count * 250)
    }

def create_quick_config(scale: str, output_dir: str) -> dict:
    """根据规模创建快速配置"""
    system_info = get_system_info()
    
    configs = {
        'small': {
            'molecules': 1000,
            'batch_size': 250,
            'max_workers': min(4, system_info['cpu_cores']),
            'description': '小规模测试 (1千个分子)'
        },
        'medium': {
            'molecules': 10000,
            'batch_size': 1000,
            'max_workers': min(8, system_info['cpu_cores']),
            'description': '中等规模 (1万个分子)'
        },
        'large': {
            'molecules': 50000,
            'batch_size': 2000,
            'max_workers': min(16, system_info['cpu_cores']),
            'description': '大规模 (5万个分子)'
        },
        'ultra': {
            'molecules': 100000,
            'batch_size': 5000,
            'max_workers': min(32, system_info['cpu_cores']),
            'description': '超大规模 (10万个分子)'
        }
    }
    
    config = configs.get(scale, configs['medium'])
    
    return {
        'description': f'自动生成配置 - {config["description"]}',
        'generation_settings': {
            'target_molecules': config['molecules'],
            'batch_size': config['batch_size'],
            'max_workers': config['max_workers']
        },
        'quality_criteria': {
            'min_stability': 35.0,
            'max_energy': 60.0
        },
        'output_settings': {
            'output_dir': output_dir,
            'export_formats': ['csv'],
            'auto_analysis': True
        },
        'system_info': system_info
    }

def estimate_time_and_resources(molecules: int, workers: int) -> dict:
    """估算时间和资源需求"""
    # 基于基准测试的经验数据
    base_speed = 50  # 分子/秒/worker (保守估计)
    estimated_speed = workers * base_speed * 0.8  # 考虑并行效率损失
    
    estimated_time_seconds = molecules / estimated_speed
    estimated_memory_mb = workers * 100 + molecules * 0.01  # 保守估计
    
    return {
        'estimated_time_minutes': estimated_time_seconds / 60,
        'estimated_memory_mb': estimated_memory_mb,
        'estimated_speed': estimated_speed
    }

def main():
    """主程序"""
    parser = argparse.ArgumentParser(
        description='高性能含氟分子生成器快速启动',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python quick_start.py --scale small                # 生成1千个分子(测试)
  python quick_start.py --scale medium               # 生成1万个分子
  python quick_start.py --scale large                # 生成5万个分子
  python quick_start.py --scale ultra                # 生成10万个分子
  python quick_start.py --custom 25000 --workers 16  # 自定义规模
  python quick_start.py --benchmark                  # 运行性能基准测试
        """
    )
    
    parser.add_argument('--scale', choices=['small', 'medium', 'large', 'ultra'],
                       help='预设规模')
    parser.add_argument('--custom', type=int, 
                       help='自定义分子数量')
    parser.add_argument('--workers', type=int,
                       help='工作进程数量')
    parser.add_argument('--output-dir', default=None,
                       help='输出目录')
    parser.add_argument('--min-stability', type=float, default=35.0,
                       help='最小稳定性阈值')
    parser.add_argument('--max-energy', type=float, default=60.0,
                       help='最大能量阈值')
    parser.add_argument('--benchmark', action='store_true',
                       help='运行性能基准测试')
    parser.add_argument('--dry-run', action='store_true',
                       help='仅显示配置信息，不实际运行')
    
    args = parser.parse_args()
    
    # 获取系统信息
    system_info = get_system_info()
    
    print("高性能含氟分子生成器")
    print("=" * 50)
    print(f"系统信息:")
    print(f"  CPU核心数: {system_info['cpu_cores']}")
    print(f"  内存: {system_info['memory_gb']:.1f} GB")
    print(f"  推荐工作进程: {system_info['recommended_workers']}")
    print()
    
    # 运行基准测试
    if args.benchmark:
        print("运行性能基准测试...")
        from batch_processor import PerformanceProfiler
        profiler = PerformanceProfiler()
        profiler.profile_generation_methods()
        return
    
    # 确定配置
    if args.custom:
        molecules = args.custom
        workers = args.workers or system_info['recommended_workers']
        output_dir = args.output_dir or f"custom_{molecules}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        
        config = {
            'molecules': molecules,
            'workers': workers,
            'min_stability': args.min_stability,
            'max_energy': args.max_energy,
            'output_dir': output_dir
        }
        scale_description = f"自定义规模 ({molecules:,} 个分子)"
        
    elif args.scale:
        output_dir = args.output_dir or f"{args.scale}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        config_dict = create_quick_config(args.scale, output_dir)
        
        config = {
            'molecules': config_dict['generation_settings']['target_molecules'],
            'workers': config_dict['generation_settings']['max_workers'],
            'min_stability': config_dict['quality_criteria']['min_stability'],
            'max_energy': config_dict['quality_criteria']['max_energy'],
            'output_dir': output_dir
        }
        scale_description = config_dict['description']
        
    else:
        # 默认中等规模
        args.scale = 'medium'
        output_dir = f"medium_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        config_dict = create_quick_config('medium', output_dir)
        
        config = {
            'molecules': config_dict['generation_settings']['target_molecules'],
            'workers': config_dict['generation_settings']['max_workers'],
            'min_stability': config_dict['quality_criteria']['min_stability'],
            'max_energy': config_dict['quality_criteria']['max_energy'],
            'output_dir': output_dir
        }
        scale_description = config_dict['description']
    
    # 估算资源需求
    estimates = estimate_time_and_resources(config['molecules'], config['workers'])
    
    print(f"生成配置:")
    print(f"  规模: {scale_description}")
    print(f"  目标分子数: {config['molecules']:,}")
    print(f"  工作进程数: {config['workers']}")
    print(f"  最小稳定性: {config['min_stability']}")
    print(f"  最大能量: {config['max_energy']} kcal/mol")
    print(f"  输出目录: {config['output_dir']}")
    print()
    
    print(f"资源估算:")
    print(f"  预计耗时: {estimates['estimated_time_minutes']:.1f} 分钟")
    print(f"  预计内存: {estimates['estimated_memory_mb']:.0f} MB")
    print(f"  预计速度: {estimates['estimated_speed']:.0f} 分子/秒")
    print()
    
    # 检查资源
    if estimates['estimated_memory_mb'] > system_info['memory_gb'] * 800:  # 80% 内存使用率警告
        print("⚠️  警告: 预计内存使用量较高，建议减少工作进程数量")
        
    if config['workers'] > system_info['cpu_cores']:
        print("⚠️  警告: 工作进程数超过CPU核心数，可能影响性能")
    
    if args.dry_run:
        print("干运行模式，配置检查完成。")
        return
    
    # 确认运行
    try:
        confirm = input("是否开始生成? (y/N): ").lower().strip()
        if confirm not in ['y', 'yes']:
            print("已取消生成。")
            return
    except KeyboardInterrupt:
        print("\n已取消生成。")
        return
    
    # 创建输出目录
    os.makedirs(config['output_dir'], exist_ok=True)
    
    # 保存配置
    config_file = os.path.join(config['output_dir'], 'generation_config.json')
    with open(config_file, 'w', encoding='utf-8') as f:
        json.dump({
            'config': config,
            'system_info': system_info,
            'estimates': estimates,
            'timestamp': datetime.now().isoformat()
        }, f, indent=2, ensure_ascii=False)
    
    print(f"开始生成 {config['molecules']:,} 个分子...")
    print(f"配置已保存到: {config_file}")
    print("=" * 50)
    
    # 启动生成
    try:
        from mol_generator_high_performance import HighPerformanceMoleculeLibraryGenerator
        
        generator = HighPerformanceMoleculeLibraryGenerator(
            max_workers=config['workers'],
            batch_size=min(config['molecules'] // 10, 2000)
        )
        
        # 生成分子
        report = generator.generate_large_scale_library(
            num_molecules=config['molecules'],
            min_stability=config['min_stability'],
            max_energy=config['max_energy']
        )
        
        # 分析结果
        analysis = generator.analyze_large_library()
        
        # 导出结果
        output_file = os.path.join(config['output_dir'], 'generated_molecules.csv')
        generator.export_top_molecules(min(config['molecules'], 10000), 'energy')
        
        # 生成报告
        final_report = {
            'config': config,
            'generation_report': report,
            'analysis': analysis,
            'files': {
                'database': generator.db.db_path,
                'csv_export': output_file
            },
            'completion_time': datetime.now().isoformat()
        }
        
        report_file = os.path.join(config['output_dir'], 'final_report.json')
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(final_report, f, indent=2, ensure_ascii=False, default=str)
        
        print(f"\n🎉 生成完成!")
        print(f"总分子数: {analysis['total_molecules']:,}")
        print(f"成功率: {report['success_rate']*100:.1f}%")
        print(f"生成速度: {report['molecules_per_second']:.1f} 分子/秒")
        print(f"结果保存在: {config['output_dir']}")
        print(f"详细报告: {report_file}")
        
    except KeyboardInterrupt:
        print("\n生成被用户中断。")
    except Exception as e:
        print(f"\n生成过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
