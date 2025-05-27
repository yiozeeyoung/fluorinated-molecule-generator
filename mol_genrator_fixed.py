"""
高效含氟分子生成器 - 修复版本
============================

专注于能量优化和稳定性评估的含氟分子生成系统
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, Draw
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple, Optional
import logging
from dataclasses import dataclass

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class MoleculeResult:
    """分子结果数据类"""
    molecule: Chem.Mol
    smiles: str
    energy: float
    stability_score: float
    properties: Dict
    optimization_converged: bool

class StabilityGuidedFluorineMoleculeGenerator:
    """基于稳定性引导的含氟分子生成器"""
    
    def __init__(self):
        pass
        
    def generate_simple_fluoroalkane(self, smiles_template: str, fluorine_positions: List[int]) -> Chem.Mol:
        """从SMILES模板生成含氟分子"""
        try:
            # 从SMILES创建分子
            mol = Chem.MolFromSmiles(smiles_template)
            if mol is None:
                return None
            
            # 创建可编辑的分子
            rw_mol = Chem.RWMol(mol)
            
            # 在指定位置替换氢原子为氟原子
            for pos in fluorine_positions:
                if pos < rw_mol.GetNumAtoms():
                    atom = rw_mol.GetAtomWithIdx(pos)
                    if atom.GetAtomicNum() == 6:  # 碳原子
                        # 添加氟原子
                        fluorine_idx = rw_mol.AddAtom(Chem.Atom(9))  # 氟原子
                        rw_mol.AddBond(pos, fluorine_idx, Chem.BondType.SINGLE)
            
            # 转换回普通分子并清理
            result_mol = rw_mol.GetMol()
            Chem.SanitizeMol(result_mol)
            
            return result_mol
            
        except Exception as e:
            logger.error(f"分子生成错误: {e}")
            return None
    
    def generate_fluoroalkane_library(self, num_molecules: int = 20) -> List[Chem.Mol]:
        """生成含氟烷烃库"""
        molecules = []
        
        # 预定义的烷烃SMILES模板
        alkane_templates = [
            "CCCC",      # 丁烷
            "CCCCC",     # 戊烷
            "CCCCCC",    # 己烷
            "CCCCCCC",   # 庚烷
            "CCCCCCCC",  # 辛烷
            "CC(C)CCC",  # 2-甲基戊烷
            "CCCC(C)C",  # 2-甲基戊烷
            "C1CCCCC1",  # 环己烷
            "C1CCCC1",   # 环戊烷
        ]
        
        logger.info(f"开始生成{num_molecules}个含氟分子...")
        
        for i in range(num_molecules):
            try:
                # 随机选择模板
                template = random.choice(alkane_templates)
                base_mol = Chem.MolFromSmiles(template)
                
                if base_mol is None:
                    continue
                
                # 选择氟化位置
                carbon_count = sum(1 for atom in base_mol.GetAtoms() if atom.GetAtomicNum() == 6)
                max_fluorines = min(carbon_count, 3)  # 限制氟原子数量
                num_fluorines = random.randint(1, max_fluorines)
                
                carbon_indices = [i for i, atom in enumerate(base_mol.GetAtoms()) 
                                if atom.GetAtomicNum() == 6]
                fluorine_positions = random.sample(carbon_indices, 
                                                 min(num_fluorines, len(carbon_indices)))
                
                # 生成含氟分子
                fluorinated_mol = self.generate_simple_fluoroalkane(template, fluorine_positions)
                
                if fluorinated_mol is not None:
                    molecules.append(fluorinated_mol)
                    
            except Exception as e:
                logger.error(f"生成分子{i+1}时出错: {e}")
                continue
        
        logger.info(f"成功生成{len(molecules)}个含氟分子")
        return molecules
    
    def generate_diverse_fluorinated_molecules(self, num_molecules: int = 50) -> List[Chem.Mol]:
        """生成多样化的含氟分子"""
        molecules = []
        
        # 扩展的分子模板库
        molecular_templates = {
            'linear_alkanes': [
                "CCCC", "CCCCC", "CCCCCC", "CCCCCCC", "CCCCCCCC"
            ],
            'branched_alkanes': [
                "CC(C)CC", "CCC(C)C", "CC(C)CCC", "CCCC(C)C",
                "CC(C)(C)CC", "CCC(C)(C)C"
            ],
            'cyclic_compounds': [
                "C1CCCC1", "C1CCCCC1", "C1CCCCCCC1",
                "C1CCC(C)CC1", "C1CCCC(C)C1"
            ],
            'aromatics': [
                "c1ccccc1", "c1ccc(C)cc1", "c1ccc(CC)cc1"
            ],
            'hetero_cycles': [
                "C1CCOC1", "C1CCNC1", "C1CCSC1"
            ]
        }
        
        logger.info(f"开始生成{num_molecules}个多样化含氟分子...")
        
        for i in range(num_molecules):
            try:
                # 随机选择分子类型和模板
                mol_type = random.choice(list(molecular_templates.keys()))
                template = random.choice(molecular_templates[mol_type])
                
                # 生成含氟衍生物
                fluorinated_mol = self._generate_fluorinated_derivative(template)
                
                if fluorinated_mol is not None:
                    molecules.append(fluorinated_mol)
                    
            except Exception as e:
                logger.error(f"生成分子{i+1}时出错: {e}")
                continue
        
        logger.info(f"成功生成{len(molecules)}个多样化含氟分子")
        return molecules
    
    def _generate_fluorinated_derivative(self, template_smiles: str) -> Chem.Mol:
        """从模板生成含氟衍生物"""
        try:
            base_mol = Chem.MolFromSmiles(template_smiles)
            if base_mol is None:
                return None
            
            # 创建可编辑分子
            rw_mol = Chem.RWMol(base_mol)
            
            # 智能选择氟化位置
            fluorination_sites = self._select_fluorination_sites(base_mol)
            
            # 应用氟化
            for site_idx, fluorine_count in fluorination_sites.items():
                for _ in range(fluorine_count):
                    try:
                        # 添加氟原子
                        fluorine_idx = rw_mol.AddAtom(Chem.Atom(9))
                        rw_mol.AddBond(site_idx, fluorine_idx, Chem.BondType.SINGLE)
                    except Exception as e:
                        logger.debug(f"添加氟原子失败: {e}")
                        continue
            
            # 转换并净化分子
            result_mol = rw_mol.GetMol()
            Chem.SanitizeMol(result_mol)
            
            return result_mol
            
        except Exception as e:
            logger.error(f"氟化衍生物生成错误: {e}")
            return None
    
    def _select_fluorination_sites(self, mol: Chem.Mol) -> Dict[int, int]:
        """智能选择氟化位置"""
        sites = {}
        
        # 获取可氟化的碳原子
        carbon_atoms = [(i, atom) for i, atom in enumerate(mol.GetAtoms()) 
                       if atom.GetAtomicNum() == 6]
        
        if not carbon_atoms:
            return sites
        
        # 限制氟化程度
        max_fluorination = min(len(carbon_atoms), 4)
        num_sites = random.randint(1, max_fluorination)
        
        # 随机选择位点
        selected_carbons = random.sample(carbon_atoms, num_sites)
        
        for idx, atom in selected_carbons:
            # 每个碳原子最多添加1-2个氟原子
            max_f_per_carbon = min(2, 4 - atom.GetTotalDegree())
            if max_f_per_carbon > 0:
                sites[idx] = random.randint(1, max_f_per_carbon)
        
        return sites
    
    def calculate_energy(self, mol: Chem.Mol) -> float:
        """计算分子能量"""
        try:
            mol_with_h = Chem.AddHs(mol)
            
            # 嵌入3D坐标
            if AllChem.EmbedMolecule(mol_with_h, randomSeed=42) != 0:
                # 尝试多次嵌入
                for i in range(5):
                    if AllChem.EmbedMolecule(mol_with_h, randomSeed=42+i) == 0:
                        break
                else:
                    return float('inf')
            
            # 首先尝试MMFF94
            try:
                # 检查MMFF参数是否可用
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol_with_h)
                if mmff_props is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(mol_with_h, mmff_props)
                    if ff is not None:
                        energy = ff.CalcEnergy()
                        return energy
            except Exception as mmff_error:
                logger.debug(f"MMFF94计算失败: {mmff_error}")
            
            # 如果MMFF失败，尝试UFF
            try:
                ff = AllChem.UFFGetMoleculeForceField(mol_with_h)
                if ff is not None:
                    energy = ff.CalcEnergy()
                    return energy
            except Exception as uff_error:
                logger.debug(f"UFF计算失败: {uff_error}")
            
            # 如果两种力场都失败，返回基于分子描述符的估算能量
            return self._estimate_energy_from_descriptors(mol)
            
        except Exception as e:
            logger.error(f"能量计算错误: {e}")
            return self._estimate_energy_from_descriptors(mol)
    def optimize_geometry(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """优化分子几何结构"""
        try:
            mol_with_h = Chem.AddHs(mol)
            
            # 嵌入3D坐标
            embed_result = AllChem.EmbedMolecule(mol_with_h, randomSeed=42)
            if embed_result != 0:
                # 尝试多次嵌入
                for i in range(5):
                    if AllChem.EmbedMolecule(mol_with_h, randomSeed=42+i) == 0:
                        embed_result = 0
                        break
                if embed_result != 0:
                    return mol, False
            
            # 首先尝试MMFF优化
            converged = False
            try:
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol_with_h)
                if mmff_props is not None:
                    result = AllChem.MMFFOptimizeMolecule(mol_with_h, mmffVariant='MMFF94', maxIters=500)
                    converged = (result == 0)
            except Exception as mmff_error:
                logger.debug(f"MMFF优化失败: {mmff_error}")
            
            # 如果MMFF优化失败，尝试UFF优化
            if not converged:
                try:
                    result = AllChem.UFFOptimizeMolecule(mol_with_h, maxIters=500)
                    converged = (result == 0)
                except Exception as uff_error:
                    logger.debug(f"UFF优化失败: {uff_error}")
            
            optimized_mol = Chem.RemoveHs(mol_with_h)
            return optimized_mol, converged
            
        except Exception as e:
            logger.error(f"几何优化错误: {e}")
            return mol, False
    
    def calculate_stability_score(self, mol: Chem.Mol) -> float:
        """计算稳定性评分"""
        try:
            score = 100.0  # 基础分数
            
            # 检查分子结构合理性
            try:
                Chem.SanitizeMol(mol)
            except:
                score -= 50
            
            # 氟原子数量检查
            fluorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
            if fluorine_count > 0:
                score += 10  # 氟原子存在加分
            
            # 检查过度氟化
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            if carbon_count > 0 and fluorine_count / carbon_count > 0.5:
                score -= 20  # 过度氟化扣分
            
            # 分子量合理性
            mw = Descriptors.ExactMolWt(mol)
            if 50 <= mw <= 300:
                score += 5
            else:
                score -= 10
            
            return max(0, min(100, score))
            
        except Exception as e:
            logger.error(f"稳定性评分错误: {e}")
            return 0.0
    
    def calculate_molecular_properties(self, mol: Chem.Mol) -> Dict:
        """计算分子性质"""
        try:
            properties = {
                'molecular_weight': Descriptors.ExactMolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'tpsa': Descriptors.TPSA(mol),
                'hbd': Lipinski.NumHDonors(mol),
                'hba': Lipinski.NumHAcceptors(mol),
                'fluorine_count': sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9),
                'carbon_count': sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6),
                'heavy_atoms': mol.GetNumHeavyAtoms()
            }
            return properties
            
        except Exception as e:
            logger.error(f"分子性质计算错误: {e}")
            return {}
    
    def evaluate_molecule(self, mol: Chem.Mol) -> MoleculeResult:
        """全面评估分子"""
        try:
            # 几何优化
            optimized_mol, converged = self.optimize_geometry(mol)
            
            # 能量计算
            energy = self.calculate_energy(optimized_mol)
            
            # 稳定性评估
            stability_score = self.calculate_stability_score(optimized_mol)
            
            # 分子性质
            properties = self.calculate_molecular_properties(optimized_mol)
            
            # SMILES
            smiles = Chem.MolToSmiles(optimized_mol)
            
            return MoleculeResult(
                molecule=optimized_mol,
                smiles=smiles,
                energy=energy,
                stability_score=stability_score,
                properties=properties,
                optimization_converged=converged
            )
            
        except Exception as e:
            logger.error(f"分子评估错误: {e}")
            return None

    def _estimate_energy_from_descriptors(self, mol: Chem.Mol) -> float:
        """基于分子描述符估算能量"""
        try:
            # 基于分子性质的简单能量估算
            mw = Descriptors.ExactMolWt(mol)
            tpsa = Descriptors.TPSA(mol)
            logp = Descriptors.MolLogP(mol)
            heavy_atoms = mol.GetNumHeavyAtoms()
            
            # 简单的能量估算公式（基于经验关系）
            # 较大分子通常具有更高的绝对能量
            base_energy = heavy_atoms * 2.5  # 每个重原子约2.5 kcal/mol
            
            # 氟原子的贡献
            fluorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
            fluorine_contribution = fluorine_count * 1.2  # 氟原子稳定化效应
            
            # 分子紧密度影响
            compactness = heavy_atoms / (tpsa + 1)  # 避免除零
            compactness_effect = compactness * 0.5
            
            estimated_energy = base_energy - fluorine_contribution + compactness_effect
            
            return estimated_energy
            
        except Exception as e:
            logger.error(f"能量估算错误: {e}")
            return 50.0  # 默认值

class MoleculeAnalyzer:
    """分子分析工具"""
    
    @staticmethod
    def analyze_library(results: List[MoleculeResult]):
        """分析分子库"""
        if not results:
            print("没有分子可分析")
            return
        
        print(f"\n=== 分子库分析报告 ===")
        print(f"总分子数: {len(results)}")
        
        # 能量统计
        energies = [r.energy for r in results if r.energy != float('inf')]
        if energies:
            print(f"\n能量统计:")
            print(f"  平均能量: {np.mean(energies):.2f} kcal/mol")
            print(f"  最低能量: {np.min(energies):.2f} kcal/mol")
            print(f"  最高能量: {np.max(energies):.2f} kcal/mol")
        
        # 稳定性统计
        stabilities = [r.stability_score for r in results]
        print(f"\n稳定性统计:")
        print(f"  平均稳定性: {np.mean(stabilities):.2f}")
        print(f"  最高稳定性: {np.max(stabilities):.2f}")
        
        # 分子性质统计
        mw_list = [r.properties.get('molecular_weight', 0) for r in results]
        f_count_list = [r.properties.get('fluorine_count', 0) for r in results]
        
        print(f"\n分子性质统计:")
        print(f"  平均分子量: {np.mean(mw_list):.2f}")
        print(f"  平均氟原子数: {np.mean(f_count_list):.1f}")
        
        # 优化收敛统计
        converged_count = sum(1 for r in results if r.optimization_converged)
        print(f"\n优化统计:")
        print(f"  收敛分子数: {converged_count}/{len(results)} ({converged_count/len(results)*100:.1f}%)")
    
    @staticmethod
    def show_best_molecules(results: List[MoleculeResult], top_n: int = 5):
        """显示最优分子"""
        if not results:
            return
        
        # 按稳定性排序
        sorted_by_stability = sorted(results, key=lambda x: x.stability_score, reverse=True)
        
        print(f"\n=== 最稳定的分子 (Top {top_n}) ===")
        for i, result in enumerate(sorted_by_stability[:top_n], 1):
            print(f"{i}. SMILES: {result.smiles}")
            print(f"   稳定性: {result.stability_score:.2f}")
            print(f"   能量: {result.energy:.2f} kcal/mol")
            print(f"   分子量: {result.properties.get('molecular_weight', 0):.2f}")
            print(f"   氟原子数: {result.properties.get('fluorine_count', 0)}")
            print()
        
        # 按能量排序
        valid_energy_results = [r for r in results if r.energy != float('inf')]
        if valid_energy_results:
            sorted_by_energy = sorted(valid_energy_results, key=lambda x: x.energy)
            
            print(f"\n=== 能量最低的分子 (Top {top_n}) ===")
            for i, result in enumerate(sorted_by_energy[:top_n], 1):
                print(f"{i}. SMILES: {result.smiles}")
                print(f"   能量: {result.energy:.2f} kcal/mol")
                print(f"   稳定性: {result.stability_score:.2f}")
                print(f"   分子量: {result.properties.get('molecular_weight', 0):.2f}")
                print(f"   氟原子数: {result.properties.get('fluorine_count', 0)}")
                print()
    
    @staticmethod
    def visualize_molecules(results: List[MoleculeResult], max_display: int = 12):
        """可视化分子"""
        if not results:
            return
        
        try:
            molecules = [r.molecule for r in results[:max_display]]
            legends = [f"S:{r.stability_score:.1f}\nE:{r.energy:.1f}" 
                      for r in results[:max_display]]
            
            img = Draw.MolsToGridImage(molecules, 
                                     molsPerRow=4,
                                     subImgSize=(200, 200),
                                     legends=legends)
            
            plt.figure(figsize=(12, 8))
            plt.imshow(img)
            plt.title("Generated Fluorinated Molecules")
            plt.axis('off')
            plt.show()
            
        except Exception as e:
            logger.error(f"分子可视化错误: {e}")
    
    @staticmethod
    def plot_energy_distribution(results: List[MoleculeResult]):
        """绘制能量分布"""
        try:
            energies = [r.energy for r in results if r.energy != float('inf')]
            
            if not energies:
                print("没有有效的能量数据")
                return
            
            plt.figure(figsize=(10, 6))
            plt.hist(energies, bins=20, alpha=0.7, edgecolor='black')
            plt.xlabel('Energy (kcal/mol)')
            plt.ylabel('Frequency')
            plt.title('Energy Distribution of Generated Molecules')
            plt.grid(True, alpha=0.3)
            plt.show()
            
        except Exception as e:
            logger.error(f"能量分布图绘制错误: {e}")

def main():
    """主程序"""
    print("含氟分子生成器 - 修复版本")
    print("=" * 40)
    
    # 初始化生成器
    generator = StabilityGuidedFluorineMoleculeGenerator()
    analyzer = MoleculeAnalyzer()
    
    # 生成分子库
    molecules = generator.generate_fluoroalkane_library(num_molecules=30)
    
    if not molecules:
        print("未生成任何分子")
        return
    
    # 评估所有分子
    print("正在评估分子...")
    results = []
    for i, mol in enumerate(molecules):
        result = generator.evaluate_molecule(mol)
        if result:
            results.append(result)
        if (i + 1) % 10 == 0:
            print(f"已评估 {i + 1}/{len(molecules)} 个分子")
    
    if not results:
        print("没有成功评估的分子")
        return
    
    # 分析结果
    analyzer.analyze_library(results)
    analyzer.show_best_molecules(results)
    
    # 可视化
    print("\n生成可视化...")
    analyzer.visualize_molecules(results)
    analyzer.plot_energy_distribution(results)
    
    # 导出数据
    print("\n导出数据...")
    export_data = []
    for result in results:
        row = {
            'SMILES': result.smiles,
            'Energy': result.energy,
            'StabilityScore': result.stability_score,
            'OptimizationConverged': result.optimization_converged
        }
        row.update(result.properties)
        export_data.append(row)
    
    df = pd.DataFrame(export_data)
    df.to_csv("fluorinated_molecules_results.csv", index=False)
    print(f"已导出{len(results)}个分子数据到 fluorinated_molecules_results.csv")
    
    print("\n程序执行完毕！")

if __name__ == "__main__":
    main()
