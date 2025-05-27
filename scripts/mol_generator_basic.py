"""
高效含氟分子生成器 - 专注于能量优化和稳定性评估
=====================================

本模块实现了一个专门用于生成稳定、低能量含氟分子的系统。
主要特点：
- 基于分子力学的能量计算和优化
- 多重稳定性评估机制
- 智能的氟原子位置优化
- 批量生成和筛选功能
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdForceFieldHelpers
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import Draw
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

class MolecularEnergyCalculator:
    """分子能量计算器"""
    
    def __init__(self):
        self.available_forcefields = ['MMFF94', 'UFF']
        
    def calculate_mmff_energy(self, mol: Chem.Mol) -> float:
        """使用MMFF94力场计算能量"""
        try:
            mol_copy = Chem.Mol(mol)
            # 添加氢原子
            mol_with_h = Chem.AddHs(mol_copy)
            
            # 嵌入3D坐标
            if AllChem.EmbedMolecule(mol_with_h, randomSeed=42) != 0:
                logger.warning("3D坐标嵌入失败，使用备用方法")
                return float('inf')
            
            # 创建MMFF力场
            ff = AllChem.MMFFGetMoleculeForceField(mol_with_h)
            if ff is None:
                logger.warning("MMFF力场创建失败")
                return float('inf')
                
            energy = ff.CalcEnergy()
            return energy
            
        except Exception as e:
            logger.error(f"MMFF能量计算错误: {e}")
            return float('inf')
    
    def calculate_uff_energy(self, mol: Chem.Mol) -> float:
        """使用UFF力场计算能量"""
        try:
            mol_copy = Chem.Mol(mol)
            mol_with_h = Chem.AddHs(mol_copy)
            
            if AllChem.EmbedMolecule(mol_with_h, randomSeed=42) != 0:
                return float('inf')
            
            ff = AllChem.UFFGetMoleculeForceField(mol_with_h)
            if ff is None:
                return float('inf')
                
            energy = ff.CalcEnergy()
            return energy
            
        except Exception as e:
            logger.error(f"UFF能量计算错误: {e}")
            return float('inf')
    
    def calculate_energy(self, mol: Chem.Mol, method: str = 'MMFF94') -> float:
        """计算分子能量"""
        if method == 'MMFF94':
            return self.calculate_mmff_energy(mol)
        elif method == 'UFF':
            return self.calculate_uff_energy(mol)
        else:
            raise ValueError(f"不支持的力场方法: {method}")

class GeometryOptimizer:
    """几何结构优化器"""
    
    def __init__(self):
        self.max_iterations = 1000
        self.energy_tolerance = 1e-6
        
    def optimize_geometry_mmff(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """使用MMFF94优化几何结构"""
        try:
            mol_copy = Chem.Mol(mol)
            mol_with_h = Chem.AddHs(mol_copy)
            
            # 嵌入3D坐标
            if AllChem.EmbedMolecule(mol_with_h, randomSeed=42) != 0:
                return mol, False
            
            # MMFF优化
            result = AllChem.MMFFOptimizeMolecule(mol_with_h, maxIters=self.max_iterations)
            
            # 移除氢原子，返回重原子结构
            mol_optimized = Chem.RemoveHs(mol_with_h)
            converged = (result == 0)
            
            return mol_optimized, converged
            
        except Exception as e:
            logger.error(f"MMFF几何优化错误: {e}")
            return mol, False
    
    def optimize_geometry_uff(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """使用UFF优化几何结构"""
        try:
            mol_copy = Chem.Mol(mol)
            mol_with_h = Chem.AddHs(mol_copy)
            
            if AllChem.EmbedMolecule(mol_with_h, randomSeed=42) != 0:
                return mol, False
            
            result = AllChem.UFFOptimizeMolecule(mol_with_h, maxIters=self.max_iterations)
            mol_optimized = Chem.RemoveHs(mol_with_h)
            converged = (result == 0)
            
            return mol_optimized, converged
            
        except Exception as e:
            logger.error(f"UFF几何优化错误: {e}")
            return mol, False
    
    def generate_conformers(self, mol: Chem.Mol, num_confs: int = 10) -> List[Chem.Mol]:
        """生成多个构象"""
        try:
            mol_copy = Chem.Mol(mol)
            mol_with_h = Chem.AddHs(mol_copy)
            
            # 生成多个构象
            conf_ids = AllChem.EmbedMultipleConfs(mol_with_h, numConfs=num_confs, randomSeed=42)
            
            if not conf_ids:
                return [mol]
            
            conformers = []
            for conf_id in conf_ids:
                conf_mol = Chem.Mol(mol_with_h, confId=conf_id)
                conf_mol = Chem.RemoveHs(conf_mol)
                conformers.append(conf_mol)
            
            return conformers
            
        except Exception as e:
            logger.error(f"构象生成错误: {e}")
            return [mol]

class StabilityAnalyzer:
    """分子稳定性分析器"""
    
    def __init__(self):
        self.energy_calculator = MolecularEnergyCalculator()
    
    def assess_thermodynamic_stability(self, mol: Chem.Mol) -> float:
        """评估热力学稳定性"""
        try:
            # 基于分子能量和结构特征的稳定性评分
            energy = self.energy_calculator.calculate_energy(mol)
            
            if energy == float('inf'):
                return 0.0
            
            # 计算分子量归一化的能量
            mol_weight = Descriptors.ExactMolWt(mol)
            normalized_energy = energy / mol_weight if mol_weight > 0 else float('inf')
            
            # 转换为稳定性评分 (能量越低，稳定性越高)
            stability_score = max(0, 100 - abs(normalized_energy))
            
            return min(100, stability_score)
            
        except Exception as e:
            logger.error(f"热力学稳定性评估错误: {e}")
            return 0.0
    
    def check_structural_strain(self, mol: Chem.Mol) -> float:
        """检查结构应变"""
        try:
            strain_score = 100.0  # 基础分数
            
            # 检查小环应变
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if len(ring) < 5:  # 小于5元环有应变
                    strain_score -= 20
                elif len(ring) == 5:  # 5元环有一定应变
                    strain_score -= 5
            
            # 检查四价氟原子（不合理）
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 9:  # 氟原子
                    if atom.GetDegree() > 1:
                        strain_score -= 50  # 氟原子不应该有多个键
            
            # 检查碳原子的化合价
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 6:  # 碳原子
                    if atom.GetDegree() > 4:
                        strain_score -= 30  # 碳原子化合价不应超过4
            
            return max(0, strain_score)
            
        except Exception as e:
            logger.error(f"结构应变检查错误: {e}")
            return 0.0
    
    def evaluate_fluorine_effects(self, mol: Chem.Mol) -> float:
        """评估氟原子的稳定性贡献"""
        try:
            stability_bonus = 0.0
            
            # 统计氟原子
            fluorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
            
            if fluorine_count == 0:
                return 0.0
            
            # 氟原子位置分析
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 9:  # 氟原子
                    neighbors = [mol.GetAtomWithIdx(n.GetIdx()) for n in atom.GetNeighbors()]
                    
                    # 氟原子连接到碳原子是稳定的
                    if len(neighbors) == 1 and neighbors[0].GetAtomicNum() == 6:
                        stability_bonus += 10
                    
                    # 检查是否在链末端（更稳定）
                    carbon_neighbor = neighbors[0] if neighbors else None
                    if carbon_neighbor and carbon_neighbor.GetDegree() <= 2:
                        stability_bonus += 5
            
            return min(50, stability_bonus)  # 最大加分50
            
        except Exception as e:
            logger.error(f"氟原子效应评估错误: {e}")
            return 0.0
    
    def assess_molecular_stability(self, mol: Chem.Mol) -> float:
        """综合评估分子稳定性"""
        try:
            # 各项稳定性评分
            thermo_score = self.assess_thermodynamic_stability(mol)
            strain_score = self.check_structural_strain(mol)
            fluorine_score = self.evaluate_fluorine_effects(mol)
            
            # 加权综合评分
            total_score = (thermo_score * 0.5 + strain_score * 0.3 + fluorine_score * 0.2)
            
            return min(100, max(0, total_score))
            
        except Exception as e:
            logger.error(f"综合稳定性评估错误: {e}")
            return 0.0

class StabilityGuidedFluorineMoleculeGenerator:
    """基于稳定性引导的含氟分子生成器"""
    
    def __init__(self):
        self.energy_calculator = MolecularEnergyCalculator()
        self.geometry_optimizer = GeometryOptimizer()
        self.stability_analyzer = StabilityAnalyzer()
        
    def create_carbon_skeleton(self, carbon_count: int, branching: bool = False) -> Chem.RWMol:
        """创建碳骨架"""
        mol = Chem.RWMol()
        
        # 添加碳原子
        carbon_indices = []
        for i in range(carbon_count):
            carbon = Chem.Atom(6)  # 碳原子
            idx = mol.AddAtom(carbon)
            carbon_indices.append(idx)
        
        # 连接主链
        for i in range(carbon_count - 1):
            mol.AddBond(carbon_indices[i], carbon_indices[i + 1], Chem.BondType.SINGLE)
        
        # 添加分支（如果启用）
        if branching and carbon_count >= 4:
            # 在中间碳原子添加甲基分支
            branch_position = random.randint(1, carbon_count - 2)
            methyl_carbon = Chem.Atom(6)
            methyl_idx = mol.AddAtom(methyl_carbon)
            mol.AddBond(carbon_indices[branch_position], methyl_idx, Chem.BondType.SINGLE)
        
        return mol
    
    def add_fluorine_atoms(self, mol: Chem.RWMol, fluorine_positions: List[int]) -> Chem.RWMol:
        """在指定位置添加氟原子"""
        for pos in fluorine_positions:
            if pos < mol.GetNumAtoms():
                carbon_atom = mol.GetAtomWithIdx(pos)
                if carbon_atom.GetAtomicNum() == 6:  # 确保是碳原子
                    # 检查碳原子的可用化合价
                    current_bonds = carbon_atom.GetDegree()
                    if current_bonds < 4:  # 碳原子最多4个键
                        fluorine = Chem.Atom(9)  # 氟原子
                        f_idx = mol.AddAtom(fluorine)
                        mol.AddBond(pos, f_idx, Chem.BondType.SINGLE)
        
        return mol
    
    def generate_linear_fluoroalkane(self, carbon_count: int, fluorine_ratio: float = 0.3) -> Chem.Mol:
        """生成线性氟代烷烃"""
        try:
            # 创建碳骨架
            mol = self.create_carbon_skeleton(carbon_count, branching=False)
            
            # 计算氟原子数量
            max_fluorine = min(carbon_count * 2, int(carbon_count * 4 * fluorine_ratio))
            fluorine_count = random.randint(1, max(1, max_fluorine))
            
            # 随机选择氟原子位置（偏向链末端）
            carbon_indices = list(range(carbon_count))
            weights = [3 if i in [0, carbon_count-1] else 1 for i in carbon_indices]
            
            fluorine_positions = []
            for _ in range(fluorine_count):
                if carbon_indices:
                    pos = random.choices(carbon_indices, weights=weights[:len(carbon_indices)])[0]
                    fluorine_positions.append(pos)
                    # 移除已选择的位置，避免过度氟化
                    idx = carbon_indices.index(pos)
                    carbon_indices.pop(idx)
                    weights.pop(idx)
            
            # 添加氟原子
            mol = self.add_fluorine_atoms(mol, fluorine_positions)
              # 获取分子并进行清理
            final_mol = mol.GetMol()
            
            # 清理分子
            try:
                Chem.SanitizeMol(final_mol)
            except:
                logger.warning("分子清理失败")
                return None
            
            return final_mol
            
        except Exception as e:
            logger.error(f"线性氟代烷烃生成错误: {e}")
            return None
    
    def generate_branched_fluoroalkane(self, main_chain_length: int, fluorine_ratio: float = 0.3) -> Chem.Mol:
        """生成分支氟代烷烃"""
        try:
            mol = self.create_carbon_skeleton(main_chain_length, branching=True)
            
            total_carbons = mol.GetNumAtoms()
            max_fluorine = int(total_carbons * 2 * fluorine_ratio)
            fluorine_count = random.randint(1, max(1, max_fluorine))
            
            # 选择氟原子位置
            carbon_indices = [i for i in range(total_carbons) if mol.GetAtomWithIdx(i).GetAtomicNum() == 6]
            fluorine_positions = random.sample(carbon_indices, min(fluorine_count, len(carbon_indices)))
            mol = self.add_fluorine_atoms(mol, fluorine_positions)
            
            # 获取分子并进行清理
            final_mol = mol.GetMol()
            try:
                Chem.SanitizeMol(final_mol)
            except:
                logger.warning("分子清理失败")
                return None
            
            return final_mol
            
        except Exception as e:
            logger.error(f"分支氟代烷烃生成错误: {e}")
            return None
    
    def generate_cyclic_fluoroalkane(self, ring_size: int = 6, fluorine_count: int = 2) -> Chem.Mol:
        """生成环状氟代烷烃"""
        try:
            mol = Chem.RWMol()
            
            # 创建环状碳骨架
            carbon_indices = []
            for i in range(ring_size):
                carbon = Chem.Atom(6)
                idx = mol.AddAtom(carbon)
                carbon_indices.append(idx)
            
            # 形成环
            for i in range(ring_size):
                next_i = (i + 1) % ring_size
                mol.AddBond(carbon_indices[i], carbon_indices[next_i], Chem.BondType.SINGLE)
            
            # 添加氟原子
            fluorine_positions = random.sample(carbon_indices, min(fluorine_count, ring_size))            
            mol = self.add_fluorine_atoms(mol, fluorine_positions)
            
            # 获取分子并进行清理
            final_mol = mol.GetMol()
            try:
                Chem.SanitizeMol(final_mol)
            except:
                logger.warning("分子清理失败")
                return None
            
            return final_mol
            
        except Exception as e:
            logger.error(f"环状氟代烷烃生成错误: {e}")
            return None
    
    def optimize_molecule(self, mol: Chem.Mol) -> MoleculeResult:
        """优化分子并返回结果"""
        try:
            # 几何优化
            optimized_mol, converged = self.geometry_optimizer.optimize_geometry_mmff(mol)
            
            # 能量计算
            energy = self.energy_calculator.calculate_energy(optimized_mol)
            
            # 稳定性评估
            stability_score = self.stability_analyzer.assess_molecular_stability(optimized_mol)
            
            # 计算分子性质
            properties = self.calculate_molecular_properties(optimized_mol)
            
            # 生成SMILES
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
            logger.error(f"分子优化错误: {e}")
            return None
    
    def calculate_molecular_properties(self, mol: Chem.Mol) -> Dict:
        """计算分子性质"""
        try:
            properties = {
                'molecular_weight': Descriptors.ExactMolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'tpsa': Descriptors.TPSA(mol),
                'hbd': Lipinski.NumHDonors(mol),
                'hba': Lipinski.NumHAcceptors(mol),
                'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
                'aromatic_rings': Descriptors.NumAromaticRings(mol),
                'fluorine_count': sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9),
                'carbon_count': sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6),
                'heavy_atoms': mol.GetNumHeavyAtoms()
            }
            return properties
            
        except Exception as e:
            logger.error(f"分子性质计算错误: {e}")
            return {}

class MoleculeLibraryGenerator:
    """分子库生成器"""
    
    def __init__(self):
        self.generator = StabilityGuidedFluorineMoleculeGenerator()
        
    def generate_molecule_library(self, 
                                num_molecules: int = 100,
                                carbon_range: Tuple[int, int] = (3, 15),
                                fluorine_ratio_range: Tuple[float, float] = (0.1, 0.5),
                                include_branched: bool = True,
                                include_cyclic: bool = True,
                                min_stability_threshold: float = 50.0) -> List[MoleculeResult]:
        """生成分子库"""
        
        molecules = []
        generation_methods = ['linear']
        
        if include_branched:
            generation_methods.append('branched')
        if include_cyclic:
            generation_methods.append('cyclic')
        
        logger.info(f"开始生成{num_molecules}个含氟分子...")
        
        for i in range(num_molecules):
            try:
                # 随机选择生成方法
                method = random.choice(generation_methods)
                carbon_count = random.randint(*carbon_range)
                fluorine_ratio = random.uniform(*fluorine_ratio_range)
                
                # 生成分子
                if method == 'linear':
                    mol = self.generator.generate_linear_fluoroalkane(carbon_count, fluorine_ratio)
                elif method == 'branched':
                    mol = self.generator.generate_branched_fluoroalkane(carbon_count, fluorine_ratio)
                elif method == 'cyclic':
                    ring_size = min(carbon_count, 8)  # 限制环大小
                    fluorine_count = max(1, int(ring_size * fluorine_ratio))
                    mol = self.generator.generate_cyclic_fluoroalkane(ring_size, fluorine_count)
                
                if mol is None:
                    continue
                
                # 验证分子
                try:
                    Chem.SanitizeMol(mol)
                except:
                    logger.warning(f"分子{i+1}验证失败，跳过")
                    continue
                
                # 优化和评估
                result = self.generator.optimize_molecule(mol)
                
                if result and result.stability_score >= min_stability_threshold:
                    molecules.append(result)
                    if len(molecules) % 10 == 0:
                        logger.info(f"已生成{len(molecules)}个合格分子...")
                
            except Exception as e:
                logger.error(f"生成分子{i+1}时出错: {e}")
                continue
        
        logger.info(f"完成！共生成{len(molecules)}个合格的含氟分子")
        return molecules
    
    def filter_molecules_by_criteria(self, 
                                   molecules: List[MoleculeResult],
                                   max_energy: float = None,
                                   min_stability: float = None,
                                   molecular_weight_range: Tuple[float, float] = None) -> List[MoleculeResult]:
        """根据条件筛选分子"""
        
        filtered = []
        
        for mol_result in molecules:
            # 能量筛选
            if max_energy is not None and mol_result.energy > max_energy:
                continue
            
            # 稳定性筛选
            if min_stability is not None and mol_result.stability_score < min_stability:
                continue
            
            # 分子量筛选
            if molecular_weight_range is not None:
                mw = mol_result.properties.get('molecular_weight', 0)
                if not (molecular_weight_range[0] <= mw <= molecular_weight_range[1]):
                    continue
            
            filtered.append(mol_result)
        
        return filtered

class MoleculeVisualizer:
    """分子可视化工具"""
    
    @staticmethod
    def show_molecule(mol: Chem.Mol, title: str = "Molecule", size: Tuple[int, int] = (400, 400)):
        """显示单个分子"""
        try:
            img = Draw.MolToImage(mol, size=size)
            plt.figure(figsize=(6, 4))
            plt.imshow(img)
            plt.title(title)
            plt.axis('off')
            plt.show()
        except Exception as e:
            logger.error(f"分子可视化错误: {e}")
    
    @staticmethod
    def show_molecule_grid(molecules: List[MoleculeResult], 
                          max_molecules: int = 20,
                          mols_per_row: int = 4,
                          img_size: Tuple[int, int] = (200, 200)):
        """网格显示多个分子"""
        try:
            mol_list = [result.molecule for result in molecules[:max_molecules]]
            legends = [f"E:{result.energy:.1f}\nS:{result.stability_score:.1f}" 
                      for result in molecules[:max_molecules]]
            
            img = Draw.MolsToGridImage(mol_list, 
                                     molsPerRow=mols_per_row,
                                     subImgSize=img_size,
                                     legends=legends)
            
            plt.figure(figsize=(12, 8))
            plt.imshow(img)
            plt.title("Generated Fluorinated Molecules")
            plt.axis('off')
            plt.show()
            
        except Exception as e:
            logger.error(f"分子网格可视化错误: {e}")
    
    @staticmethod
    def plot_energy_distribution(molecules: List[MoleculeResult]):
        """绘制能量分布图"""
        try:
            energies = [mol.energy for mol in molecules if mol.energy != float('inf')]
            
            plt.figure(figsize=(10, 6))
            plt.hist(energies, bins=30, alpha=0.7, edgecolor='black')
            plt.xlabel('Energy (kcal/mol)')
            plt.ylabel('Frequency')
            plt.title('Energy Distribution of Generated Molecules')
            plt.grid(True, alpha=0.3)
            plt.show()
            
        except Exception as e:
            logger.error(f"能量分布图绘制错误: {e}")
    
    @staticmethod
    def plot_stability_vs_energy(molecules: List[MoleculeResult]):
        """绘制稳定性vs能量散点图"""
        try:
            energies = []
            stabilities = []
            
            for mol in molecules:
                if mol.energy != float('inf'):
                    energies.append(mol.energy)
                    stabilities.append(mol.stability_score)
            
            plt.figure(figsize=(10, 6))
            plt.scatter(energies, stabilities, alpha=0.6)
            plt.xlabel('Energy (kcal/mol)')
            plt.ylabel('Stability Score')
            plt.title('Stability Score vs Energy')
            plt.grid(True, alpha=0.3)
            plt.show()
            
        except Exception as e:
            logger.error(f"稳定性-能量图绘制错误: {e}")

class MoleculeExporter:
    """分子数据导出工具"""
    
    @staticmethod
    def export_to_sdf(molecules: List[MoleculeResult], filename: str):
        """导出为SDF格式"""
        try:
            writer = Chem.SDWriter(filename)
            
            for i, mol_result in enumerate(molecules):
                mol = mol_result.molecule
                
                # 添加属性
                mol.SetProp("Energy", str(mol_result.energy))
                mol.SetProp("StabilityScore", str(mol_result.stability_score))
                mol.SetProp("SMILES", mol_result.smiles)
                
                # 添加分子性质
                for prop, value in mol_result.properties.items():
                    mol.SetProp(prop, str(value))
                
                writer.write(mol)
            
            writer.close()
            logger.info(f"已导出{len(molecules)}个分子到 {filename}")
            
        except Exception as e:
            logger.error(f"SDF导出错误: {e}")
    
    @staticmethod
    def export_to_csv(molecules: List[MoleculeResult], filename: str):
        """导出为CSV格式"""
        try:
            data = []
            
            for mol_result in molecules:
                row = {
                    'SMILES': mol_result.smiles,
                    'Energy': mol_result.energy,
                    'StabilityScore': mol_result.stability_score,
                    'OptimizationConverged': mol_result.optimization_converged
                }
                
                # 添加分子性质
                row.update(mol_result.properties)
                data.append(row)
            
            df = pd.DataFrame(data)
            df.to_csv(filename, index=False)
            logger.info(f"已导出{len(molecules)}个分子到 {filename}")
            
        except Exception as e:
            logger.error(f"CSV导出错误: {e}")

def analyze_molecule_library(molecules: List[MoleculeResult]):
    """分析分子库统计信息"""
    if not molecules:
        print("没有分子可分析")
        return
    
    print(f"\n=== 分子库分析报告 ===")
    print(f"总分子数: {len(molecules)}")
    
    # 能量统计
    energies = [mol.energy for mol in molecules if mol.energy != float('inf')]
    if energies:
        print(f"\n能量统计:")
        print(f"  平均能量: {np.mean(energies):.2f} kcal/mol")
        print(f"  最低能量: {np.min(energies):.2f} kcal/mol")
        print(f"  最高能量: {np.max(energies):.2f} kcal/mol")
        print(f"  能量标准差: {np.std(energies):.2f} kcal/mol")
    
    # 稳定性统计
    stabilities = [mol.stability_score for mol in molecules]
    print(f"\n稳定性统计:")
    print(f"  平均稳定性: {np.mean(stabilities):.2f}")
    print(f"  最高稳定性: {np.max(stabilities):.2f}")
    print(f"  最低稳定性: {np.min(stabilities):.2f}")
    
    # 分子性质统计
    mw_list = [mol.properties.get('molecular_weight', 0) for mol in molecules]
    f_count_list = [mol.properties.get('fluorine_count', 0) for mol in molecules]
    
    print(f"\n分子性质统计:")
    print(f"  平均分子量: {np.mean(mw_list):.2f}")
    print(f"  平均氟原子数: {np.mean(f_count_list):.1f}")
    
    # 优化收敛统计
    converged_count = sum(1 for mol in molecules if mol.optimization_converged)
    print(f"\n优化统计:")
    print(f"  收敛分子数: {converged_count}/{len(molecules)} ({converged_count/len(molecules)*100:.1f}%)")

def find_most_stable_molecules(molecules: List[MoleculeResult], top_n: int = 10) -> List[MoleculeResult]:
    """找到最稳定的分子"""
    sorted_molecules = sorted(molecules, key=lambda x: x.stability_score, reverse=True)
    return sorted_molecules[:top_n]

def find_lowest_energy_molecules(molecules: List[MoleculeResult], top_n: int = 10) -> List[MoleculeResult]:
    """找到能量最低的分子"""
    valid_molecules = [mol for mol in molecules if mol.energy != float('inf')]
    sorted_molecules = sorted(valid_molecules, key=lambda x: x.energy)
    return sorted_molecules[:top_n]

def main_example():
    """主要使用示例"""
    print("含氟分子生成器 - 能量和稳定性优化版本")
    print("=" * 50)
    
    # 初始化生成器
    library_generator = MoleculeLibraryGenerator()
    visualizer = MoleculeVisualizer()
    exporter = MoleculeExporter()
    
    # 配置生成参数
    generation_config = {
        'num_molecules': 50,  # 生成分子数量
        'carbon_range': (4, 12),  # 碳原子数范围
        'fluorine_ratio_range': (0.1, 0.4),  # 氟化程度范围
        'include_branched': True,  # 包含分支结构
        'include_cyclic': True,   # 包含环状结构
        'min_stability_threshold': 60.0  # 最低稳定性阈值
    }
    
    print(f"生成参数:")
    for key, value in generation_config.items():
        print(f"  {key}: {value}")
    print()
    
    # 生成分子库
    molecules = library_generator.generate_molecule_library(**generation_config)
    
    if not molecules:
        print("未生成任何合格分子，请检查参数设置")
        return
    
    # 分析结果
    analyze_molecule_library(molecules)
    
    # 找到最优分子
    print(f"\n=== 最稳定的分子 (Top 5) ===")
    most_stable = find_most_stable_molecules(molecules, 5)
    for i, mol in enumerate(most_stable, 1):
        print(f"{i}. SMILES: {mol.smiles}")
        print(f"   稳定性: {mol.stability_score:.2f}, 能量: {mol.energy:.2f} kcal/mol")
        print(f"   分子量: {mol.properties.get('molecular_weight', 0):.2f}")
        print(f"   氟原子数: {mol.properties.get('fluorine_count', 0)}")
        print()
    
    print(f"\n=== 能量最低的分子 (Top 5) ===")
    lowest_energy = find_lowest_energy_molecules(molecules, 5)
    for i, mol in enumerate(lowest_energy, 1):
        print(f"{i}. SMILES: {mol.smiles}")
        print(f"   能量: {mol.energy:.2f} kcal/mol, 稳定性: {mol.stability_score:.2f}")
        print(f"   分子量: {mol.properties.get('molecular_weight', 0):.2f}")
        print(f"   氟原子数: {mol.properties.get('fluorine_count', 0)}")
        print()
    
    # 可视化
    print("生成可视化图表...")
    
    # 显示最优分子结构
    if len(most_stable) > 0:
        print("显示最稳定的分子结构:")
        visualizer.show_molecule(most_stable[0].molecule, 
                               f"Most Stable Molecule (Score: {most_stable[0].stability_score:.1f})")
    
    # 显示分子网格
    print("显示分子库网格:")
    visualizer.show_molecule_grid(molecules[:16], max_molecules=16)
    
    # 显示统计图表
    print("显示能量分布:")
    visualizer.plot_energy_distribution(molecules)
    
    print("显示稳定性vs能量关系:")
    visualizer.plot_stability_vs_energy(molecules)
    
    # 导出数据
    print("\n导出数据...")
    
    # 过滤高质量分子
    high_quality = library_generator.filter_molecules_by_criteria(
        molecules,
        min_stability=70.0,
        molecular_weight_range=(50, 300)
    )
    
    print(f"高质量分子数量: {len(high_quality)}")
    
    # 导出为不同格式
    if high_quality:
        exporter.export_to_csv(high_quality, "high_quality_fluorinated_molecules.csv")
        exporter.export_to_sdf(high_quality[:10], "top_fluorinated_molecules.sdf")
    
    print("\n程序执行完毕！")

def generate_custom_molecule_example():
    """自定义分子生成示例"""
    print("\n=== 自定义分子生成示例 ===")
    
    generator = StabilityGuidedFluorineMoleculeGenerator()
    visualizer = MoleculeVisualizer()
    
    # 生成特定类型的分子
    print("生成线性氟代烷烃...")
    linear_mol = generator.generate_linear_fluoroalkane(carbon_count=8, fluorine_ratio=0.25)
    if linear_mol:
        result = generator.optimize_molecule(linear_mol)
        if result:
            print(f"线性分子 - SMILES: {result.smiles}")
            print(f"稳定性: {result.stability_score:.2f}, 能量: {result.energy:.2f}")
            visualizer.show_molecule(result.molecule, "Linear Fluoroalkane")
    
    print("生成分支氟代烷烃...")
    branched_mol = generator.generate_branched_fluoroalkane(main_chain_length=6, fluorine_ratio=0.3)
    if branched_mol:
        result = generator.optimize_molecule(branched_mol)
        if result:
            print(f"分支分子 - SMILES: {result.smiles}")
            print(f"稳定性: {result.stability_score:.2f}, 能量: {result.energy:.2f}")
            visualizer.show_molecule(result.molecule, "Branched Fluoroalkane")
    
    print("生成环状氟代烷烃...")
    cyclic_mol = generator.generate_cyclic_fluoroalkane(ring_size=6, fluorine_count=2)
    if cyclic_mol:
        result = generator.optimize_molecule(cyclic_mol)
        if result:
            print(f"环状分子 - SMILES: {result.smiles}")
            print(f"稳定性: {result.stability_score:.2f}, 能量: {result.energy:.2f}")
            visualizer.show_molecule(result.molecule, "Cyclic Fluoroalkane")

if __name__ == "__main__":
    # 运行主要示例
    main_example()
    
    # 运行自定义示例
    generate_custom_molecule_example()
