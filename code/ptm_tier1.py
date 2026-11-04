#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PTM Tier-1 建模: C端酰胺化 + 焦谷氨酸 + 羟脯氨酸 (PyRosetta)
用法: rosetta_env python ptm_tier1.py <序列名> <序列> [修饰说明]
修饰说明格式: "amide" | "pgq" | "hyp@4,10,21" | 组合如 "amide,hyp@6"
输出: ptm_out/<name>_ptm.pdb (单构象, 能量最小化后)
"""
import sys, os
import pyrosetta
from pyrosetta import Pose
from pyrosetta.rosetta.core.chemical import ChemicalManager
from pyrosetta.rosetta.core.pose import make_pose_from_sequence
from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory

def build(name, seq, mods):
    pyrosetta.init("-mute all -run:preserve_header")
    chem = ChemicalManager.get_instance()
    rset = chem.residue_type_set("fa_standard")
    pose = Pose()
    make_pose_from_sequence(pose, seq, "fa_standard")
    # 1. Hyp: 替换指定位置的 Pro 为 HYP
    if "hyp" in mods:
        hyp_type = rset.name_map("HYP")
        for m in mods.split(","):
            if m.startswith("hyp@"):
                for pos in [int(x) for x in m[4:].split(";")]:
                    if pose.residue(pos).name3() == "PRO":
                        newres = ResidueFactory.create_residue(hyp_type)
                        pose.replace_residue(pos, newres, True)
    # 2. 焦谷氨酸 (N端 Gln->PCA): 用 patch 或替换为 PCA 类型
    if "pgq" in mods and rset.has_name("PCA"):
        pca = rset.name_map("PCA")
        newres = ResidueFactory.create_residue(pca)
        pose.replace_residue(1, newres, True)
    # 3. C端酰胺化 patch
    if "amide" in mods:
        # Rosetta patch: Cterm_amidation
        from pyrosetta.rosetta.core.pose import add_variant_type_to_pose_residue
        from pyrosetta.rosetta.core.chemical import VariantType
        add_variant_type_to_pose_residue(pose, VariantType.CTERM_AMIDATION, pose.total_residue())
    # 4. 最小化 (骨架固定, 侧链+末端自由)
    sf = ScoreFunctionFactory.create_score_function("ref2015")
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(True)
    minm = MinMover(mm, sf, "lbfgs_armijo_nonmonotone", 0.01, True)
    minm.apply(pose)
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ptm_out")
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{name}_ptm.pdb")
    pose.dump_pdb(out)
    print(f"OK {name} -> {out}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("用法: ptm_tier1.py <名称> <序列> <修饰>"); sys.exit(1)
    mods = sys.argv[3] if len(sys.argv) > 3 else ""
    build(sys.argv[1], sys.argv[2].replace("X","").upper(), mods)
