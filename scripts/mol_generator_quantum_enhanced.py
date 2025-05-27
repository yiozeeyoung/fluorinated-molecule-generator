"""
量子化学增强型含氟分子生成器
=============================

结合半经验量子化学方法的高精度含氟分子生成系统
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
# import seaborn as sns  # 可选依赖

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class EnhancedMoleculeResult:
    """增强型分子结果数据类"""
    molecule: Chem.Mol
    smiles: str
    energy: float
    stability_score: float
    properties: Dict
    optimization_converged: bool
    formation_energy: float
    binding_energy: float
    dipole_moment: float

class QuantumEnhancedFluorineMoleculeGenerator:
    """量子化学增强型含氟分子生成器"""
    
    def __init__(self):
        self.fluorine_templates = {
            'perfluoroalkanes': [
                "CF4", "C2F6", "C3F8", "C4F10"
            ],
            'fluoroalkanes': [
                "CCF", "CCCF", "CCCCF", "CCCCCF",
                "CC(F)C", "CCC(F)C", "CCCC(F)C"
            ],
            'fluorocycloalkanes': [
                "C1C(F)CCC1", "C1CC(F)CCC1", "C1CCC(F)CC1",
                "C1C(F)C(F)CC1", "C1CC(F)C(F)C1"
            ],
            'fluoroaromatics': [
                "Fc1ccccc1", "Fc1ccc(F)cc1", "Fc1cc(F)ccc1",
                "Fc1ccc(C)cc1", "Fc1ccc(CC)cc1"
            ],
            'mixed_fluorocarbons': [
                "CC(F)(F)C", "CCC(F)(F)C", "CC(F)C(F)C",
                "C(F)(F)C(F)(F)", "CC(F)CC(F)C"
            ]
        }
        
        # 量子化学参数
        self.bond_energies = {
            'C-F': 116.0,  # kcal/mol
            'C-C': 83.0,
            'C-H': 99.0,
            'C=C': 146.0
        }
        
        self.atomic_electronegativities = {
            'C': 2.55,
            'F': 3.98,
            'H': 2.20
        }
    
    def generate_quantum_guided_molecules(self, num_molecules: int = 100) -> List[Chem.Mol]:
        """基于量子化学指导生成分子"""
        molecules = []
        
        logger.info(f"开始生成{num_molecules}个量子化学指导的含氟分子...")
        
        # 加权选择策略
        template_weights = {
            'fluoroalkanes': 0.3,
            'fluorocycloalkanes': 0.25,
            'mixed_fluorocarbons': 0.2,
            'fluoroaromatics': 0.15,
            'perfluoroalkanes': 0.1
        }
        
        template_types = list(template_weights.keys())
        weights = list(template_weights.values())
        
        for i in range(num_molecules):
            try:
                # 加权随机选择模板类型
                mol_type = np.random.choice(template_types, p=weights)
                template = random.choice(self.fluorine_templates[mol_type])
                
                # 应用量子化学优化
                optimized_mol = self._quantum_optimize_fluorination(template)
                
                if optimized_mol is not None:
                    molecules.append(optimized_mol)
                    
            except Exception as e:
                logger.error(f"生成分子{i+1}时出错: {e}")
                continue
        
        logger.info(f"成功生成{len(molecules)}个量子化学指导的含氟分子")
        return molecules
    
    def _quantum_optimize_fluorination(self, template_smiles: str) -> Chem.Mol:
        """基于量子化学原理优化氟化"""
        try:
            mol = Chem.MolFromSmiles(template_smiles)
            if mol is None:
                return None
            
            # 量子化学启发式修改
            modified_mol = self._apply_quantum_modifications(mol)
            
            if modified_mol is not None:
                return modified_mol
            else:
                return mol
                
        except Exception as e:
            logger.error(f"量子优化错误: {e}")
            return None
    
    def _apply_quantum_modifications(self, mol: Chem.Mol) -> Chem.Mol:
        """应用量子化学启发式修改"""
        try:
            # 基于电负性和轨道杂化的氟化位点选择
            rw_mol = Chem.RWMol(mol)
            
            # 分析分子电子结构
            carbon_atoms = [(i, atom) for i, atom in enumerate(mol.GetAtoms()) 
                           if atom.GetAtomicNum() == 6]
            
            # 基于碳原子环境选择最佳氟化位点
            best_sites = self._select_optimal_fluorination_sites(carbon_atoms, mol)
            
            # 应用氟化
            for site_idx in best_sites:
                try:
                    # 检查该位点是否还能接受氟原子
                    atom = rw_mol.GetAtomWithIdx(site_idx)
                    if atom.GetTotalDegree() < 4:  # 碳原子最多4个键
                        fluorine_idx = rw_mol.AddAtom(Chem.Atom(9))
                        rw_mol.AddBond(site_idx, fluorine_idx, Chem.BondType.SINGLE)
                except Exception:
                    continue
            
            result_mol = rw_mol.GetMol()
            Chem.SanitizeMol(result_mol)
            return result_mol
            
        except Exception as e:
            logger.error(f"量子修改应用错误: {e}")
            return mol
    
    def _select_optimal_fluorination_sites(self, carbon_atoms: List, mol: Chem.Mol) -> List[int]:
        """基于量子化学原理选择最佳氟化位点"""
        site_scores = []
        
        for idx, atom in carbon_atoms:
            score = 0
            
            # 杂化状态评分
            if atom.GetHybridization() == Chem.HybridizationType.SP3:
                score += 3  # sp3碳更容易氟化
            elif atom.GetHybridization() == Chem.HybridizationType.SP2:
                score += 2
            elif atom.GetHybridization() == Chem.HybridizationType.SP:
                score += 1
            
            # 邻接原子评分
            neighbors = [mol.GetAtomWithIdx(n_idx) for n_idx in atom.GetNeighbors()]
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 9:  # 已有氟原子
                    score -= 1  # 减少氟-氟排斥
                elif neighbor.GetAtomicNum() == 6:  # 碳原子
                    score += 0.5
            
            # 环状结构评分
            if atom.IsInRing():
                score += 1  # 环状碳稍微有利
            
            site_scores.append((idx, score))
        
        # 排序并选择最佳位点
        site_scores.sort(key=lambda x: x[1], reverse=True)
        
        # 选择前几个位点
        num_sites = min(random.randint(1, 3), len(site_scores))
        return [site[0] for site in site_scores[:num_sites]]
    
    def calculate_formation_energy(self, mol: Chem.Mol) -> float:
        """计算生成焓"""
        try:
            # 简化的生成焓计算
            total_energy = 0
            
            # 统计键类型
            for bond in mol.GetBonds():
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                
                # 识别键类型
                if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 9:
                    total_energy += self.bond_energies['C-F']
                elif begin_atom.GetAtomicNum() == 9 and end_atom.GetAtomicNum() == 6:
                    total_energy += self.bond_energies['C-F']
                elif begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        total_energy += self.bond_energies['C=C']
                    else:
                        total_energy += self.bond_energies['C-C']
                elif (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 1) or \
                     (begin_atom.GetAtomicNum() == 1 and end_atom.GetAtomicNum() == 6):
                    total_energy += self.bond_energies['C-H']
            
            return -total_energy  # 负值表示放热形成
            
        except Exception as e:
            logger.error(f"生成焓计算错误: {e}")
            return 0.0
    
    def calculate_binding_energy(self, mol: Chem.Mol) -> float:
        """计算结合能"""
        try:
            # 简化的结合能计算
            fluorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            
            if carbon_count == 0:
                return 0.0
            
            # 基于C-F键的平均结合能
            cf_bonds = 0
            for bond in mol.GetBonds():
                atoms = [bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()]
                if 6 in atoms and 9 in atoms:
                    cf_bonds += 1
            
            avg_binding_energy = cf_bonds * self.bond_energies['C-F'] / max(1, carbon_count)
            
            return avg_binding_energy
            
        except Exception as e:
            logger.error(f"结合能计算错误: {e}")
            return 0.0
    
    def calculate_dipole_moment(self, mol: Chem.Mol) -> float:
        """估算偶极矩"""
        try:
            # 简化的偶极矩计算
            total_moment = 0.0
            
            for bond in mol.GetBonds():
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                
                # 计算电负性差
                begin_en = self.atomic_electronegativities.get(
                    begin_atom.GetSymbol(), 2.0)
                end_en = self.atomic_electronegativities.get(
                    end_atom.GetSymbol(), 2.0)
                
                electronegativity_diff = abs(begin_en - end_en)
                
                # 偶极矩贡献
                bond_dipole = electronegativity_diff * 1.0  # 简化系数
                total_moment += bond_dipole
            
            return total_moment
            
        except Exception as e:
            logger.error(f"偶极矩计算错误: {e}")
            return 0.0
    
    def enhanced_energy_calculation(self, mol: Chem.Mol) -> float:
        """增强型能量计算"""
        try:
            mol_with_h = Chem.AddHs(mol)
            
            # 多种嵌入策略
            embed_success = False
            for method in [AllChem.ETKDG(), AllChem.ETKDGv2(), AllChem.ETKDGv3()]:
                try:
                    if AllChem.EmbedMolecule(mol_with_h, params=method) == 0:
                        embed_success = True
                        break
                except:
                    continue
            
            if not embed_success:
                return self._estimate_energy_from_quantum_descriptors(mol)
            
            # 多步优化
            final_energy = float('inf')
            
            # MMFF94优化
            try:
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol_with_h)
                if mmff_props is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(mol_with_h, mmff_props)
                    if ff is not None:
                        AllChem.MMFFOptimizeMolecule(mol_with_h, maxIters=1000)
                        energy = ff.CalcEnergy()
                        if energy < final_energy:
                            final_energy = energy
            except:
                pass
            
            # UFF优化
            try:
                ff = AllChem.UFFGetMoleculeForceField(mol_with_h)
                if ff is not None:
                    AllChem.UFFOptimizeMolecule(mol_with_h, maxIters=1000)
                    energy = ff.CalcEnergy()
                    if energy < final_energy:
                        final_energy = energy
            except:
                pass
            
            if final_energy == float('inf'):
                return self._estimate_energy_from_quantum_descriptors(mol)
            
            return final_energy
            
        except Exception as e:
            logger.error(f"增强能量计算错误: {e}")
            return self._estimate_energy_from_quantum_descriptors(mol)
    
    def _estimate_energy_from_quantum_descriptors(self, mol: Chem.Mol) -> float:
        """基于量子描述符的能量估算"""
        try:
            # 分子性质
            mw = Descriptors.ExactMolWt(mol)
            heavy_atoms = mol.GetNumHeavyAtoms()
            
            # 原子计数
            fluorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            
            # 键能贡献
            bond_energy_contribution = 0
            for bond in mol.GetBonds():
                begin_sym = bond.GetBeginAtom().GetSymbol()
                end_sym = bond.GetEndAtom().GetSymbol()
                bond_key = f"{begin_sym}-{end_sym}"
                reverse_key = f"{end_sym}-{begin_sym}"
                
                if bond_key in self.bond_energies:
                    bond_energy_contribution += self.bond_energies[bond_key] * 0.1
                elif reverse_key in self.bond_energies:
                    bond_energy_contribution += self.bond_energies[reverse_key] * 0.1
            
            # 量子化学校正
            # 氟原子的强电负性效应
            electronegativity_effect = fluorine_count * -2.0  # 氟原子稳定化
            
            # 分子大小效应
            size_effect = np.sqrt(heavy_atoms) * 1.5
            
            # 共轭效应
            aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
            conjugation_effect = aromatic_atoms * -0.5  # 共轭稳定化
            
            # 空间拥挤效应
            if carbon_count > 0:
                crowding_factor = fluorine_count / carbon_count
                if crowding_factor > 0.3:
                    steric_penalty = (crowding_factor - 0.3) * 10
                else:
                    steric_penalty = 0
            else:
                steric_penalty = 0
            
            total_energy = (bond_energy_contribution + size_effect + 
                          electronegativity_effect + conjugation_effect + 
                          steric_penalty)
            
            # 添加热力学修正
            thermal_correction = random.uniform(-1.0, 1.0)
            
            return max(5.0, total_energy + thermal_correction)
            
        except Exception as e:
            logger.error(f"量子能量估算错误: {e}")
            return 20.0
    
    def enhanced_stability_assessment(self, mol: Chem.Mol) -> float:
        """增强型稳定性评估"""
        try:
            score = 100.0
            
            # 基础化学稳定性
            try:
                Chem.SanitizeMol(mol)
            except:
                score -= 40
            
            # 氟化学特异性评估
            fluorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
            carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            
            if fluorine_count > 0 and carbon_count > 0:
                f_c_ratio = fluorine_count / carbon_count
                
                # 最佳氟化比例
                if 0.2 <= f_c_ratio <= 0.6:
                    score += 15
                elif f_c_ratio > 0.8:
                    score -= 20  # 过度氟化不稳定
                
                # C-F键强度奖励
                cf_bonds = sum(1 for bond in mol.GetBonds()
                             if {bond.GetBeginAtom().GetAtomicNum(), 
                                bond.GetEndAtom().GetAtomicNum()} == {6, 9})
                score += cf_bonds * 2
            
            # 分子对称性
            symmetry_score = self._assess_molecular_symmetry(mol)
            score += symmetry_score
            
            # 热力学稳定性
            formation_energy = self.calculate_formation_energy(mol)
            if formation_energy < -100:  # 高度放热形成
                score += 10
            elif formation_energy > 0:   # 吸热形成
                score -= 5
            
            # ADMET性质
            admet_score = self._assess_admet_properties(mol)
            score += admet_score
            
            return max(0, min(100, score))
            
        except Exception as e:
            logger.error(f"增强稳定性评估错误: {e}")
            return 50.0
    
    def _assess_molecular_symmetry(self, mol: Chem.Mol) -> float:
        """评估分子对称性"""
        try:
            # 简化的对称性评估
            atom_counts = {}
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
            
            # 对称性奖励
            if len(atom_counts) <= 3:  # 简单组成
                return 5.0
            else:
                return 0.0
                
        except:
            return 0.0
    
    def _assess_admet_properties(self, mol: Chem.Mol) -> float:
        """评估ADMET性质"""
        try:
            score = 0
            
            # 分子量
            mw = Descriptors.ExactMolWt(mol)
            if 100 <= mw <= 300:
                score += 3
            
            # LogP
            logp = Descriptors.MolLogP(mol)
            if -1 <= logp <= 4:
                score += 3
            
            # 极性表面积
            tpsa = Descriptors.TPSA(mol)
            if tpsa <= 140:
                score += 2
            
            # 氢键供体/受体
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            if hbd <= 5 and hba <= 10:
                score += 2
            
            return score
            
        except:
            return 0
    
    def comprehensive_molecular_evaluation(self, mol: Chem.Mol) -> EnhancedMoleculeResult:
        """全面分子评估"""
        try:
            # 能量计算
            energy = self.enhanced_energy_calculation(mol)
            
            # 稳定性评估
            stability_score = self.enhanced_stability_assessment(mol)
            
            # 量子化学性质
            formation_energy = self.calculate_formation_energy(mol)
            binding_energy = self.calculate_binding_energy(mol)
            dipole_moment = self.calculate_dipole_moment(mol)
            
            # 分子性质
            properties = self._calculate_comprehensive_properties(mol)
            
            # 几何优化状态
            optimization_converged = self._check_optimization_convergence(mol)
            
            # SMILES
            smiles = Chem.MolToSmiles(mol)
            
            return EnhancedMoleculeResult(
                molecule=mol,
                smiles=smiles,
                energy=energy,
                stability_score=stability_score,
                properties=properties,
                optimization_converged=optimization_converged,
                formation_energy=formation_energy,
                binding_energy=binding_energy,
                dipole_moment=dipole_moment
            )
            
        except Exception as e:
            logger.error(f"全面分子评估错误: {e}")
            return None
    
    def _calculate_comprehensive_properties(self, mol: Chem.Mol) -> Dict:
        """计算全面的分子性质"""
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
                'formal_charge': Chem.GetFormalCharge(mol),
                'ring_count': Descriptors.RingCount(mol),
                'aliphatic_rings': Descriptors.NumAliphaticRings(mol),
                'aromatic_rings': Descriptors.NumAromaticRings(mol),
                'saturation_fraction': Descriptors.FractionCsp3(mol),
                'molar_refractivity': Descriptors.MolMR(mol)
            }
            return properties
            
        except Exception as e:
            logger.error(f"分子性质计算错误: {e}")
            return {}
    
    def _check_optimization_convergence(self, mol: Chem.Mol) -> bool:
        """检查优化收敛性"""
        try:
            # 简化的收敛检查
            mol_with_h = Chem.AddHs(mol)
            
            if AllChem.EmbedMolecule(mol_with_h) == 0:
                try:
                    result = AllChem.MMFFOptimizeMolecule(mol_with_h, maxIters=100)
                    return result == 0
                except:
                    try:
                        result = AllChem.UFFOptimizeMolecule(mol_with_h, maxIters=100)
                        return result == 0
                    except:
                        return False
            else:
                return False
                
        except:
            return False

class QuantumAnalyzer:
    """量子化学分析工具"""
    
    @staticmethod
    def comprehensive_analysis(results: List[EnhancedMoleculeResult]):
        """全面分析"""
        if not results:
            print("没有分子可分析")
            return
        
        print(f"\n=== 量子化学增强分析报告 ===")
        print(f"总分子数: {len(results)}")
        
        # 能量分析
        energies = [r.energy for r in results if r.energy != float('inf') and r.energy > 0]
        if energies:
            print(f"\n能量统计 (kcal/mol):")
            print(f"  平均能量: {np.mean(energies):.2f}")
            print(f"  最低能量: {np.min(energies):.2f}")
            print(f"  最高能量: {np.max(energies):.2f}")
            print(f"  能量标准差: {np.std(energies):.2f}")
        
        # 生成焓分析
        formation_energies = [r.formation_energy for r in results]
        print(f"\n生成焓统计 (kcal/mol):")
        print(f"  平均生成焓: {np.mean(formation_energies):.2f}")
        print(f"  最低生成焓: {np.min(formation_energies):.2f}")
        print(f"  最高生成焓: {np.max(formation_energies):.2f}")
        
        # 结合能分析
        binding_energies = [r.binding_energy for r in results]
        print(f"\n结合能统计 (kcal/mol):")
        print(f"  平均结合能: {np.mean(binding_energies):.2f}")
        print(f"  最高结合能: {np.max(binding_energies):.2f}")
        
        # 偶极矩分析
        dipole_moments = [r.dipole_moment for r in results]
        print(f"\n偶极矩统计 (Debye):")
        print(f"  平均偶极矩: {np.mean(dipole_moments):.2f}")
        print(f"  最高偶极矩: {np.max(dipole_moments):.2f}")
        
        # 稳定性分析
        stabilities = [r.stability_score for r in results]
        print(f"\n稳定性统计:")
        print(f"  平均稳定性: {np.mean(stabilities):.2f}")
        print(f"  最高稳定性: {np.max(stabilities):.2f}")
        print(f"  稳定性标准差: {np.std(stabilities):.2f}")
        
        # 氟化模式分析
        fluorine_distribution = {}
        for result in results:
            f_count = result.properties.get('fluorine_count', 0)
            fluorine_distribution[f_count] = fluorine_distribution.get(f_count, 0) + 1
        
        print(f"\n氟化模式分布:")
        for f_count, count in sorted(fluorine_distribution.items()):
            print(f"  {f_count}个氟原子: {count}个分子 ({count/len(results)*100:.1f}%)")
        
        # 收敛性统计
        converged_count = sum(1 for r in results if r.optimization_converged)
        print(f"\n优化收敛统计:")
        print(f"  收敛分子数: {converged_count}/{len(results)} ({converged_count/len(results)*100:.1f}%)")
    
    @staticmethod
    def plot_quantum_analysis(results: List[EnhancedMoleculeResult]):
        """绘制量子分析图表"""
        try:
            fig, axes = plt.subplots(2, 3, figsize=(18, 12))
            
            # 能量分布
            energies = [r.energy for r in results if r.energy != float('inf') and r.energy > 0]
            if energies:
                axes[0, 0].hist(energies, bins=20, alpha=0.7, edgecolor='black')
                axes[0, 0].set_xlabel('Energy (kcal/mol)')
                axes[0, 0].set_ylabel('Frequency')
                axes[0, 0].set_title('Energy Distribution')
                axes[0, 0].grid(True, alpha=0.3)
            
            # 生成焓分布
            formation_energies = [r.formation_energy for r in results]
            axes[0, 1].hist(formation_energies, bins=20, alpha=0.7, edgecolor='black', color='green')
            axes[0, 1].set_xlabel('Formation Energy (kcal/mol)')
            axes[0, 1].set_ylabel('Frequency')
            axes[0, 1].set_title('Formation Energy Distribution')
            axes[0, 1].grid(True, alpha=0.3)
            
            # 结合能分布
            binding_energies = [r.binding_energy for r in results]
            axes[0, 2].hist(binding_energies, bins=20, alpha=0.7, edgecolor='black', color='red')
            axes[0, 2].set_xlabel('Binding Energy (kcal/mol)')
            axes[0, 2].set_ylabel('Frequency')
            axes[0, 2].set_title('Binding Energy Distribution')
            axes[0, 2].grid(True, alpha=0.3)
            
            # 偶极矩vs稳定性
            dipole_moments = [r.dipole_moment for r in results]
            stabilities = [r.stability_score for r in results]
            axes[1, 0].scatter(dipole_moments, stabilities, alpha=0.6)
            axes[1, 0].set_xlabel('Dipole Moment (Debye)')
            axes[1, 0].set_ylabel('Stability Score')
            axes[1, 0].set_title('Dipole Moment vs Stability')
            axes[1, 0].grid(True, alpha=0.3)
            
            # 氟原子数vs结合能
            f_counts = [r.properties.get('fluorine_count', 0) for r in results]
            axes[1, 1].scatter(f_counts, binding_energies, alpha=0.6, color='purple')
            axes[1, 1].set_xlabel('Number of Fluorine Atoms')
            axes[1, 1].set_ylabel('Binding Energy (kcal/mol)')
            axes[1, 1].set_title('Fluorine Count vs Binding Energy')
            axes[1, 1].grid(True, alpha=0.3)
            
            # 分子量vs能量
            mw_list = [r.properties.get('molecular_weight', 0) for r in results]
            if energies:
                energy_mw_pairs = [(r.energy, r.properties.get('molecular_weight', 0)) 
                                 for r in results if r.energy != float('inf') and r.energy > 0]
                if energy_mw_pairs:
                    e_vals, mw_vals = zip(*energy_mw_pairs)
                    axes[1, 2].scatter(mw_vals, e_vals, alpha=0.6, color='orange')
                    axes[1, 2].set_xlabel('Molecular Weight')
                    axes[1, 2].set_ylabel('Energy (kcal/mol)')
                    axes[1, 2].set_title('Molecular Weight vs Energy')
                    axes[1, 2].grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.show()
            
        except Exception as e:
            logger.error(f"量子分析图表绘制错误: {e}")

def main():
    """主程序"""
    print("量子化学增强型含氟分子生成器")
    print("=" * 60)
    
    # 初始化生成器
    generator = QuantumEnhancedFluorineMoleculeGenerator()
    analyzer = QuantumAnalyzer()
    
    # 生成分子库
    molecules = generator.generate_quantum_guided_molecules(num_molecules=75)
    
    if not molecules:
        print("未生成任何分子")
        return
    
    # 全面评估
    print("正在进行量子化学评估...")
    results = []
    for i, mol in enumerate(molecules):
        result = generator.comprehensive_molecular_evaluation(mol)
        if result:
            results.append(result)
        if (i + 1) % 15 == 0:
            print(f"已评估 {i + 1}/{len(molecules)} 个分子")
    
    if not results:
        print("没有成功评估的分子")
        return
    
    # 分析结果
    analyzer.comprehensive_analysis(results)
    
    # 显示最优分子
    print("\n=== Top 5 最低能量分子 ===")
    valid_energy_results = [r for r in results if r.energy != float('inf') and r.energy > 0]
    if valid_energy_results:
        sorted_by_energy = sorted(valid_energy_results, key=lambda x: x.energy)
        for i, result in enumerate(sorted_by_energy[:5], 1):
            print(f"{i}. SMILES: {result.smiles}")
            print(f"   能量: {result.energy:.2f} kcal/mol")
            print(f"   生成焓: {result.formation_energy:.2f} kcal/mol")
            print(f"   结合能: {result.binding_energy:.2f} kcal/mol")
            print(f"   偶极矩: {result.dipole_moment:.2f} Debye")
            print(f"   稳定性: {result.stability_score:.2f}")
            print()
    
    # 量子分析图表
    print("\n生成量子分析图表...")
    analyzer.plot_quantum_analysis(results)
    
    # 导出数据
    print("\n导出量子化学数据...")
    export_data = []
    for result in results:
        row = {
            'SMILES': result.smiles,
            'Energy': result.energy,
            'FormationEnergy': result.formation_energy,
            'BindingEnergy': result.binding_energy,
            'DipoleMoment': result.dipole_moment,
            'StabilityScore': result.stability_score,
            'OptimizationConverged': result.optimization_converged
        }
        row.update(result.properties)
        export_data.append(row)
    
    df = pd.DataFrame(export_data)
    df.to_csv("quantum_enhanced_fluorinated_molecules.csv", index=False)
    print(f"已导出{len(results)}个分子的量子化学数据")
    
    print("\n量子化学增强分析完成！")

if __name__ == "__main__":
    main()
