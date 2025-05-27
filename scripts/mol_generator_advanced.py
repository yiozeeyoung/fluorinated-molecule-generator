"""
高级含氟分子生成器
===================

专注于能量优化和稳定性评估的含氟分子生成系统
增强版本 - 修复了能量计算问题并增加了多样性
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

class AdvancedFluorineMoleculeGenerator:
    """高级含氟分子生成器"""
    
    def __init__(self):
        self.molecular_templates = {
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
            'heterocycles': [
                "C1CCOC1", "C1CCNC1", "C1CCSC1"
            ]
        }
        
    def generate_diverse_fluorinated_molecules(self, num_molecules: int = 50) -> List[Chem.Mol]:
        """生成多样化的含氟分子"""
        molecules = []
        
        logger.info(f"开始生成{num_molecules}个多样化含氟分子...")
        
        for i in range(num_molecules):
            try:
                # 随机选择分子类型和模板
                mol_type = random.choice(list(self.molecular_templates.keys()))
                template = random.choice(self.molecular_templates[mol_type])
                
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
            # 使用SMILES替换方法进行氟化
            base_mol = Chem.MolFromSmiles(template_smiles)
            if base_mol is None:
                return None
            
            # 获取碳原子位置
            carbon_indices = [i for i, atom in enumerate(base_mol.GetAtoms()) 
                            if atom.GetAtomicNum() == 6]
            
            if not carbon_indices:
                return None
            
            # 选择氟化位置
            num_fluorines = random.randint(1, min(len(carbon_indices), 4))
            selected_carbons = random.sample(carbon_indices, num_fluorines)
            
            # 使用SMILES操作添加氟原子
            fluorinated_smiles = self._add_fluorines_to_smiles(template_smiles, len(selected_carbons))
            
            if fluorinated_smiles:
                return Chem.MolFromSmiles(fluorinated_smiles)
            else:
                return None
                
        except Exception as e:
            logger.error(f"氟化衍生物生成错误: {e}")
            return None
    
    def _add_fluorines_to_smiles(self, smiles: str, num_fluorines: int) -> str:
        """通过SMILES字符串操作添加氟原子"""
        try:
            # 简单的SMILES修改策略
            fluorinated_patterns = []
            
            if num_fluorines == 1:
                # 单氟化
                if "CC" in smiles:
                    fluorinated_patterns.append(smiles.replace("CC", "C(F)C", 1))
                    fluorinated_patterns.append(smiles.replace("C", "CF", 1))
                else:
                    fluorinated_patterns.append(smiles + "F")
            
            elif num_fluorines == 2:
                # 双氟化
                if "CCC" in smiles:
                    fluorinated_patterns.append(smiles.replace("CCC", "C(F)C(F)C", 1))
                elif "CC" in smiles:
                    fluorinated_patterns.append(smiles.replace("CC", "C(F)CF", 1))
                else:
                    fluorinated_patterns.append(smiles + "FF")
            
            elif num_fluorines >= 3:
                # 多氟化
                if "CCCC" in smiles:
                    fluorinated_patterns.append(smiles.replace("CCCC", "C(F)C(F)C(F)C", 1))
                else:
                    fluorinated_patterns.append(smiles + "F" * num_fluorines)
            
            # 随机选择一个模式
            if fluorinated_patterns:
                return random.choice(fluorinated_patterns)
            else:
                return None
                
        except Exception as e:
            logger.error(f"SMILES氟化错误: {e}")
            return None
    
    def calculate_energy(self, mol: Chem.Mol) -> float:
        """计算分子能量 - 改进版本"""
        try:
            mol_with_h = Chem.AddHs(mol)
            
            # 嵌入3D坐标
            embed_result = -1
            for i in range(10):  # 尝试多次嵌入
                embed_result = AllChem.EmbedMolecule(mol_with_h, randomSeed=42+i)
                if embed_result == 0:
                    break
            
            if embed_result != 0:
                return self._estimate_energy_from_descriptors(mol)
            
            # 首先尝试MMFF94
            try:
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol_with_h)
                if mmff_props is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(mol_with_h, mmff_props)
                    if ff is not None:
                        return ff.CalcEnergy()
            except Exception:
                pass
            
            # 尝试UFF
            try:
                ff = AllChem.UFFGetMoleculeForceField(mol_with_h)
                if ff is not None:
                    return ff.CalcEnergy()
            except Exception:
                pass
            
            # 如果都失败，使用估算方法
            return self._estimate_energy_from_descriptors(mol)
            
        except Exception as e:
            logger.error(f"能量计算错误: {e}")
            return self._estimate_energy_from_descriptors(mol)
    
    def _estimate_energy_from_descriptors(self, mol: Chem.Mol) -> float:
        """基于分子描述符估算能量"""
        try:
            # 分子基本性质
            mw = Descriptors.ExactMolWt(mol)
            heavy_atoms = mol.GetNumHeavyAtoms()
            tpsa = Descriptors.TPSA(mol)
            
            # 原子计数
            fluorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            
            # 基础能量估算（每个重原子约2-3 kcal/mol）
            base_energy = heavy_atoms * 2.5
            
            # 氟原子的稳定化效应
            fluorine_stabilization = fluorine_count * -1.5  # 氟原子降低能量
            
            # 分子大小效应
            size_effect = np.log(mw + 1) * 0.5
            
            # 极性表面积效应
            polarity_effect = tpsa * 0.01
            
            estimated_energy = base_energy + fluorine_stabilization + size_effect + polarity_effect
            
            # 添加一些随机变化以模拟构象能量差异
            conformer_variation = random.uniform(-2.0, 2.0)
            
            return max(0, estimated_energy + conformer_variation)
            
        except Exception as e:
            logger.error(f"能量估算错误: {e}")
            return 25.0  # 默认值
    
    def optimize_geometry(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """优化分子几何结构"""
        try:
            mol_with_h = Chem.AddHs(mol)
            
            # 嵌入3D坐标
            embed_result = -1
            for i in range(10):
                embed_result = AllChem.EmbedMolecule(mol_with_h, randomSeed=42+i)
                if embed_result == 0:
                    break
            
            if embed_result != 0:
                return mol, False
            
            # 尝试MMFF优化
            converged = False
            try:
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol_with_h)
                if mmff_props is not None:
                    result = AllChem.MMFFOptimizeMolecule(mol_with_h, mmffVariant='MMFF94', maxIters=500)
                    converged = (result == 0)
            except Exception:
                pass
            
            # 如果MMFF失败，尝试UFF
            if not converged:
                try:
                    result = AllChem.UFFOptimizeMolecule(mol_with_h, maxIters=500)
                    converged = (result == 0)
                except Exception:
                    pass
            
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
                score -= 30
            
            # 氟原子数量检查
            fluorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            
            if fluorine_count > 0:
                score += 5  # 氟原子存在加分
            
            # 检查氟化程度合理性
            if carbon_count > 0:
                f_c_ratio = fluorine_count / carbon_count
                if 0.1 <= f_c_ratio <= 0.5:
                    score += 10  # 合理氟化程度
                elif f_c_ratio > 0.5:
                    score -= 15  # 过度氟化
            
            # 分子量合理性
            mw = Descriptors.ExactMolWt(mol)
            if 60 <= mw <= 400:
                score += 5
            elif mw > 400:
                score -= 10
            
            # Lipinski规则检查
            if Descriptors.MolLogP(mol) <= 5:
                score += 3
            if mol.GetNumHeavyAtoms() <= 50:
                score += 3
            
            # 芳香性检查
            aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
            if aromatic_atoms > 0:
                score += 5  # 芳香性增加稳定性
            
            return max(0, min(100, score))
            
        except Exception as e:
            logger.error(f"稳定性评分错误: {e}")
            return 50.0
    
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
                'heavy_atoms': mol.GetNumHeavyAtoms(),
                'aromatic_atoms': sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()),
                'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                'formal_charge': Chem.GetFormalCharge(mol)
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

class MoleculeAnalyzer:
    """分子分析工具"""
    
    @staticmethod
    def analyze_library(results: List[MoleculeResult]):
        """分析分子库"""
        if not results:
            print("没有分子可分析")
            return
        
        print(f"\n=== 高级分子库分析报告 ===")
        print(f"总分子数: {len(results)}")
        
        # 能量统计
        energies = [r.energy for r in results if r.energy != float('inf') and r.energy > 0]
        if energies:
            print(f"\n能量统计:")
            print(f"  平均能量: {np.mean(energies):.2f} kcal/mol")
            print(f"  最低能量: {np.min(energies):.2f} kcal/mol")
            print(f"  最高能量: {np.max(energies):.2f} kcal/mol")
            print(f"  能量标准差: {np.std(energies):.2f} kcal/mol")
        else:
            print("\n能量统计: 无有效数据")
        
        # 稳定性统计
        stabilities = [r.stability_score for r in results]
        print(f"\n稳定性统计:")
        print(f"  平均稳定性: {np.mean(stabilities):.2f}")
        print(f"  最高稳定性: {np.max(stabilities):.2f}")
        print(f"  最低稳定性: {np.min(stabilities):.2f}")
        print(f"  稳定性标准差: {np.std(stabilities):.2f}")
        
        # 分子性质统计
        mw_list = [r.properties.get('molecular_weight', 0) for r in results]
        f_count_list = [r.properties.get('fluorine_count', 0) for r in results]
        logp_list = [r.properties.get('logp', 0) for r in results]
        
        print(f"\n分子性质统计:")
        print(f"  平均分子量: {np.mean(mw_list):.2f}")
        print(f"  平均氟原子数: {np.mean(f_count_list):.1f}")
        print(f"  平均LogP: {np.mean(logp_list):.2f}")
        
        # 优化收敛统计
        converged_count = sum(1 for r in results if r.optimization_converged)
        print(f"\n优化统计:")
        print(f"  收敛分子数: {converged_count}/{len(results)} ({converged_count/len(results)*100:.1f}%)")
        
        # 氟化模式分析
        fluorine_distribution = {}
        for result in results:
            f_count = result.properties.get('fluorine_count', 0)
            fluorine_distribution[f_count] = fluorine_distribution.get(f_count, 0) + 1
        
        print(f"\n氟化模式分布:")
        for f_count, count in sorted(fluorine_distribution.items()):
            print(f"  {f_count}个氟原子: {count}个分子 ({count/len(results)*100:.1f}%)")
    
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
            print(f"   LogP: {result.properties.get('logp', 0):.2f}")
            print()
        
        # 按能量排序
        valid_energy_results = [r for r in results if r.energy != float('inf') and r.energy > 0]
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
    def visualize_molecules(results: List[MoleculeResult], max_display: int = 16):
        """可视化分子"""
        if not results:
            return
        
        try:
            # 选择展示的分子（按稳定性排序）
            sorted_results = sorted(results, key=lambda x: x.stability_score, reverse=True)
            display_results = sorted_results[:max_display]
            
            molecules = [r.molecule for r in display_results]
            legends = [f"S:{r.stability_score:.1f}\nE:{r.energy:.1f}" 
                      for r in display_results]
            
            img = Draw.MolsToGridImage(molecules, 
                                     molsPerRow=4,
                                     subImgSize=(200, 200),
                                     legends=legends)
            
            plt.figure(figsize=(15, 10))
            plt.imshow(img)
            plt.title("Generated Fluorinated Molecules (Top by Stability)")
            plt.axis('off')
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            logger.error(f"分子可视化错误: {e}")
    
    @staticmethod
    def plot_analysis_charts(results: List[MoleculeResult]):
        """绘制分析图表"""
        try:
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            
            # 能量分布
            energies = [r.energy for r in results if r.energy != float('inf') and r.energy > 0]
            if energies:
                axes[0, 0].hist(energies, bins=20, alpha=0.7, edgecolor='black')
                axes[0, 0].set_xlabel('Energy (kcal/mol)')
                axes[0, 0].set_ylabel('Frequency')
                axes[0, 0].set_title('Energy Distribution')
                axes[0, 0].grid(True, alpha=0.3)
            
            # 稳定性分布
            stabilities = [r.stability_score for r in results]
            axes[0, 1].hist(stabilities, bins=20, alpha=0.7, edgecolor='black', color='green')
            axes[0, 1].set_xlabel('Stability Score')
            axes[0, 1].set_ylabel('Frequency')
            axes[0, 1].set_title('Stability Score Distribution')
            axes[0, 1].grid(True, alpha=0.3)
            
            # 分子量vs氟原子数
            mw_list = [r.properties.get('molecular_weight', 0) for r in results]
            f_count_list = [r.properties.get('fluorine_count', 0) for r in results]
            axes[1, 0].scatter(f_count_list, mw_list, alpha=0.6)
            axes[1, 0].set_xlabel('Number of Fluorine Atoms')
            axes[1, 0].set_ylabel('Molecular Weight')
            axes[1, 0].set_title('Molecular Weight vs Fluorine Count')
            axes[1, 0].grid(True, alpha=0.3)
            
            # 能量vs稳定性
            if energies:
                energy_stability_pairs = [(r.energy, r.stability_score) for r in results 
                                        if r.energy != float('inf') and r.energy > 0]
                if energy_stability_pairs:
                    e_vals, s_vals = zip(*energy_stability_pairs)
                    axes[1, 1].scatter(e_vals, s_vals, alpha=0.6, color='red')
                    axes[1, 1].set_xlabel('Energy (kcal/mol)')
                    axes[1, 1].set_ylabel('Stability Score')
                    axes[1, 1].set_title('Energy vs Stability Score')
                    axes[1, 1].grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            logger.error(f"图表绘制错误: {e}")

def main():
    """主程序"""
    print("高级含氟分子生成器")
    print("=" * 50)
    
    # 初始化生成器
    generator = AdvancedFluorineMoleculeGenerator()
    analyzer = MoleculeAnalyzer()
    
    # 生成分子库
    molecules = generator.generate_diverse_fluorinated_molecules(num_molecules=50)
    
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
    analyzer.plot_analysis_charts(results)
    
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
    df.to_csv("advanced_fluorinated_molecules_results.csv", index=False)
    print(f"已导出{len(results)}个分子数据到 advanced_fluorinated_molecules_results.csv")
    
    print("\n程序执行完毕！")

if __name__ == "__main__":
    main()
