import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import gseapy as gp
from src.methods_utils import (
    BASE_DIR,
    RANDOM_STATE, TOP_N, STRICT_FDR, STRICT_LOG2FC,
    RELAX_FDR, RELAX_LOG2FC,
    MIRNA_TARGET_FILE, LNCRNA_MIRNA_FILE,
    bh_fdr,
    load_expression_and_clusters,
    compute_deg,
    save_deg_tables,
    run_enrichr,
    smallrna_to_targets,
    build_smallrna_rna_overlap,
    build_ceRNA3_triplets,
    lncrna_to_mirna,
    postprocess_mirna_mrna_triplets
)


# =====================================================
# main：一键跑 RNA + smallRNA + lncRNA + overlap + ceRNA3 + 2 层 triplets 后处理
# =====================================================
def main():
    # ------------- RNA 分支：DEG + GO-BP -------------
    X_rna, lab_rna, deg_dir_rna, path_dir_rna, _ = load_expression_and_clusters(
        "RNA", RANDOM_STATE
    )
    deg_rna = compute_deg(X_rna, lab_rna, group1=0, group2=1)
    rna_up, rna_down, rna_relaxed = save_deg_tables(
        deg_rna, deg_dir_rna, "RNA", RANDOM_STATE
    )

    # RNA 富集
    run_enrichr(rna_up,
                label="Up-regulated genes (RNA, relaxed)",
                out_prefix="Pathway_up",
                pathway_dir=path_dir_rna,
                data_type="RNA",
                seed=RANDOM_STATE)
    run_enrichr(rna_down,
                label="Down-regulated genes (RNA, relaxed)",
                out_prefix="Pathway_down",
                pathway_dir=path_dir_rna,
                data_type="RNA",
                seed=RANDOM_STATE)

    # ------------- smallRNA 分支：DEG + 映射 + 富集 -------------
    X_s, lab_s, deg_dir_s, path_dir_s, targets_dir_s = load_expression_and_clusters(
        "smallRNA", RANDOM_STATE
    )
    deg_s = compute_deg(X_s, lab_s, group1=0, group2=1)
    small_up, small_down, small_relaxed = save_deg_tables(
        deg_s, deg_dir_s, "smallRNA", RANDOM_STATE
    )

    # smallRNA -> mRNA 映射
    targets_up, targets_down, mirna_db = smallrna_to_targets(
        small_up, small_down, targets_dir_s, RANDOM_STATE
    )

    # 富集的是靶 mRNA
    run_enrichr(targets_up,
                label="Targets of up-regulated smallRNA (relaxed)",
                out_prefix="Pathway_targets_up",
                pathway_dir=path_dir_s,
                data_type="smallRNA",
                seed=RANDOM_STATE)
    run_enrichr(targets_down,
                label="Targets of down-regulated smallRNA (relaxed)",
                out_prefix="Pathway_targets_down",
                pathway_dir=path_dir_s,
                data_type="smallRNA",
                seed=RANDOM_STATE)

    # ------------- smallRNA–target–RNA overlap + triplets -------------
    integrated_dir = os.path.join(BASE_DIR, "data", "integrated_results")
    triplet_dir = os.path.join(integrated_dir, "triplets")
    os.makedirs(triplet_dir, exist_ok=True)

    build_smallrna_rna_overlap(
        rna_up=rna_up,
        rna_down=rna_down,
        targets_up=targets_up,
        targets_down=targets_down,
        small_up=small_up,
        small_down=small_down,
        mirna_db=mirna_db,
        seed=RANDOM_STATE
    )

    # ---------- lncRNA ----------
    X_lnc, lab_lnc, deg_dir_lnc, path_dir_lnc, targets_dir_lnc = load_expression_and_clusters(
        "lncRNA", RANDOM_STATE
    )
    deg_lnc = compute_deg(X_lnc, lab_lnc, group1=0, group2=1)
    lnc_up, lnc_down, _ = save_deg_tables(
        deg_lnc, deg_dir_lnc, "lncRNA", RANDOM_STATE
    )

    # lncRNA GO 富集（可选）
    run_enrichr(lnc_up,
                label="Up-regulated lncRNA (relaxed)",
                out_prefix="Pathway_up",
                pathway_dir=path_dir_lnc,
                data_type="lncRNA",
                seed=RANDOM_STATE)
    run_enrichr(lnc_down,
                label="Down-regulated lncRNA (relaxed)",
                out_prefix="Pathway_down",
                pathway_dir=path_dir_lnc,
                data_type="lncRNA",
                seed=RANDOM_STATE)

    # lncRNA -> miRNA 映射
    mir_up, mir_down, lncrna_db = lncrna_to_mirna(
        lnc_up, lnc_down, targets_dir_lnc, RANDOM_STATE
    )

    # ---------- 3 层 ceRNA ----------
    mirna_rna_trip_file = os.path.join(
        triplet_dir,
        f"miRNA_target_RNA_triplets_seed{RANDOM_STATE}.csv"
    )
    build_ceRNA3_triplets(
        lnc_up=lnc_up,
        lnc_down=lnc_down,
        mirna_up=mir_up,
        mirna_down=mir_down,
        mirna_rna_triplets_file=mirna_rna_trip_file,
        lncrna_mirna_db=lncrna_db,
        seed=RANDOM_STATE
    )

    # ---------- 2 层 miRNA–mRNA triplets 的后处理 + GO ----------
    postprocess_mirna_mrna_triplets(seed=RANDOM_STATE)

    print("\n[DONE] RNA + smallRNA + lncRNA + smallRNA→mRNA + lncRNA→miRNA "
          "+ 2-layer & 3-layer ceRNA + GO-BP 全流程完成")


if __name__ == "__main__":
    main()
