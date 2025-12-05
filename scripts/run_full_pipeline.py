#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

# ====== 把项目根目录加到 sys.path ======
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from src.methods_utils import (
    BASE_DIR,
    RANDOM_STATE,
    load_expression_and_snf_clusters,
    compute_deg,
    save_deg_tables,
    run_enrichr,
    smallrna_to_targets,
    build_smallrna_rna_overlap,
    lncrna_to_mirna,
    build_ceRNA3_triplets,
    postprocess_mirna_mrna_triplets,
    load_deg_full_from_deseq2,
    run_r_deseq2_by_snf
)
from src.multiomics_snf_v2 import run_snf_v2_pipeline



def main():
    # =========================
    # STEP 1: 跑 SNF + EarlyFusion 聚类
    # =========================
    print("===== [STEP 1] SNF + Early-fusion multi-omics clustering =====")
    run_snf_v2_pipeline()  # 会写出 cluster_results_SNF_*.csv 等

    # =========================
    # STEP 1.5: 调 R 跑 DESeq2（基于 SNF cluster 0 vs 1）
    # =========================
    run_r_deseq2_by_snf()

    # =========================
    # STEP 2: 用 SNF 亚型标签跑 DEG + pathway + ceRNA
    # =========================
    print("\n===== [STEP 2] DEG + pathway + ceRNA driven by SNF clusters =====")

    # ---------- 2.1 RNA：DEG + GO ----------
    # 这里只是为了拿 deg_dir / path_dir 等目录信息
    _, _, deg_dir_rna, path_dir_rna, _ = load_expression_and_snf_clusters(
        "RNA", RANDOM_STATE, cluster_source="SNF"
    )

    # 用 R-DESeq2 的结果替代 Python 版 compute_deg
    deg_rna = load_deg_full_from_deseq2("RNA_SNF", RANDOM_STATE)

    # 仍然用你自己的 save_deg_tables 来做 strict / relaxed + up/down
    rna_up, rna_down, rna_relaxed = save_deg_tables(
        deg_rna, deg_dir_rna, "RNA_SNF", RANDOM_STATE
    )

    run_enrichr(
        rna_up,
        label="Up-regulated genes (RNA, SNF, relaxed)",
        out_prefix="Pathway_up",
        pathway_dir=path_dir_rna,
        data_type="RNA_SNF",
        seed=RANDOM_STATE,
    )
    run_enrichr(
        rna_down,
        label="Down-regulated genes (RNA, SNF, relaxed)",
        out_prefix="Pathway_down",
        pathway_dir=path_dir_rna,
        data_type="RNA_SNF",
        seed=RANDOM_STATE,
    )

    # ---------- 2.2 smallRNA：DEG + target + GO ----------
    # 同样只用来拿目录信息 + SNF 标签（下游可能会用到）
    _, _, deg_dir_s, path_dir_s, targets_dir_s = load_expression_and_snf_clusters(
        "smallRNA", RANDOM_STATE, cluster_source="SNF"
    )

    # 从 R-DESeq2 读取 smallRNA DEG
    deg_s = load_deg_full_from_deseq2("smallRNA_SNF", RANDOM_STATE)

    small_up, small_down, small_relaxed = save_deg_tables(
        deg_s, deg_dir_s, "smallRNA_SNF", RANDOM_STATE
    )

    # smallRNA -> mRNA
    targets_up, targets_down, mirna_db = smallrna_to_targets(
        small_up, small_down, targets_dir_s, RANDOM_STATE
    )

    run_enrichr(
        targets_up,
        label="Targets of up-regulated smallRNA (SNF, relaxed)",
        out_prefix="Pathway_targets_up",
        pathway_dir=path_dir_s,
        data_type="smallRNA_SNF",
        seed=RANDOM_STATE,
    )
    run_enrichr(
        targets_down,
        label="Targets of down-regulated smallRNA (SNF, relaxed)",
        out_prefix="Pathway_targets_down",
        pathway_dir=path_dir_s,
        data_type="smallRNA_SNF",
        seed=RANDOM_STATE,
    )

    # ---------- 2.3 smallRNA–target–RNA overlap + triplets ----------
    build_smallrna_rna_overlap(
        rna_up=rna_up,
        rna_down=rna_down,
        targets_up=targets_up,
        targets_down=targets_down,
        small_up=small_up,
        small_down=small_down,
        mirna_db=mirna_db,
        seed=RANDOM_STATE,
    )

    # ---------- 2.4 lncRNA：DEG + GO ----------
    _, _, deg_dir_lnc, path_dir_lnc, targets_dir_lnc = load_expression_and_snf_clusters(
        "lncRNA", RANDOM_STATE, cluster_source="SNF"
    )

    # 从 R-DESeq2 读取 lncRNA DEG
    deg_lnc = load_deg_full_from_deseq2("lncRNA_SNF", RANDOM_STATE)

    lnc_up, lnc_down, _ = save_deg_tables(
        deg_lnc, deg_dir_lnc, "lncRNA_SNF", RANDOM_STATE
    )

    run_enrichr(
        lnc_up,
        label="Up-regulated lncRNA (SNF, relaxed)",
        out_prefix="Pathway_up",
        pathway_dir=path_dir_lnc,
        data_type="lncRNA_SNF",
        seed=RANDOM_STATE,
    )
    run_enrichr(
        lnc_down,
        label="Down-regulated lncRNA (SNF, relaxed)",
        out_prefix="Pathway_down",
        pathway_dir=path_dir_lnc,
        data_type="lncRNA_SNF",
        seed=RANDOM_STATE,
    )

    # ---------- 2.5 lncRNA -> miRNA ----------
    mir_up, mir_down, lncrna_db = lncrna_to_mirna(
        lnc_up, lnc_down, targets_dir_lnc, RANDOM_STATE
    )

    # ---------- 2.6 三层 ceRNA：lncRNA–miRNA–mRNA ----------
    integrated_dir = os.path.join(BASE_DIR, "data", "integrated_results")
    triplet_dir = os.path.join(integrated_dir, "triplets")
    mirna_rna_trip_file = os.path.join(
        triplet_dir, f"miRNA_target_RNA_triplets_seed{RANDOM_STATE}.csv"
    )

    build_ceRNA3_triplets(
        lnc_up=lnc_up,
        lnc_down=lnc_down,
        mirna_up=mir_up,
        mirna_down=mir_down,
        mirna_rna_triplets_file=mirna_rna_trip_file,
        lncrna_mirna_db=lncrna_db,
        seed=RANDOM_STATE,
    )

    # ---------- 2.7 2 层 miRNA–mRNA triplets 后处理 + GO ----------
    postprocess_mirna_mrna_triplets(seed=RANDOM_STATE)

    print("\n[DONE] SNF-driven multi-omics DEG + pathway + ceRNA pipeline 完成。")


if __name__ == "__main__":
    main()


# def main():
#     # =========================
#     # STEP 1: 跑 SNF + EarlyFusion 聚类
#     # =========================
#     print("===== [STEP 1] SNF + Early-fusion multi-omics clustering =====")
#     run_snf_v2_pipeline()  # 会写出 cluster_results_SNF_*.csv 等

#     # =========================
#     # STEP 2: 用 SNF 亚型标签跑 DEG + pathway + ceRNA
#     # =========================
#     print("\n===== [STEP 2] DEG + pathway + ceRNA driven by SNF clusters =====")

#     # ---------- 2.1 RNA：DEG + GO ----------
#     X_rna, lab_snf, deg_dir_rna, path_dir_rna, _ = load_expression_and_snf_clusters(
#         "RNA", RANDOM_STATE, cluster_source="SNF"
#     )
#     deg_rna = compute_deg(X_rna, lab_snf, group1=0, group2=1)
#     rna_up, rna_down, rna_relaxed = save_deg_tables(
#         deg_rna, deg_dir_rna, "RNA_SNF", RANDOM_STATE
#     )

#     run_enrichr(
#         rna_up,
#         label="Up-regulated genes (RNA, SNF, relaxed)",
#         out_prefix="Pathway_up",
#         pathway_dir=path_dir_rna,
#         data_type="RNA_SNF",
#         seed=RANDOM_STATE,
#     )
#     run_enrichr(
#         rna_down,
#         label="Down-regulated genes (RNA, SNF, relaxed)",
#         out_prefix="Pathway_down",
#         pathway_dir=path_dir_rna,
#         data_type="RNA_SNF",
#         seed=RANDOM_STATE,
#     )

#     # ---------- 2.2 smallRNA：DEG + target + GO ----------
#     X_s, lab_snf_s, deg_dir_s, path_dir_s, targets_dir_s = load_expression_and_snf_clusters(
#         "smallRNA", RANDOM_STATE, cluster_source="SNF"
#     )
#     deg_s = compute_deg(X_s, lab_snf_s, group1=0, group2=1)
#     small_up, small_down, small_relaxed = save_deg_tables(
#         deg_s, deg_dir_s, "smallRNA_SNF", RANDOM_STATE
#     )

#     # smallRNA -> mRNA
#     targets_up, targets_down, mirna_db = smallrna_to_targets(
#         small_up, small_down, targets_dir_s, RANDOM_STATE
#     )

#     run_enrichr(
#         targets_up,
#         label="Targets of up-regulated smallRNA (SNF, relaxed)",
#         out_prefix="Pathway_targets_up",
#         pathway_dir=path_dir_s,
#         data_type="smallRNA_SNF",
#         seed=RANDOM_STATE,
#     )
#     run_enrichr(
#         targets_down,
#         label="Targets of down-regulated smallRNA (SNF, relaxed)",
#         out_prefix="Pathway_targets_down",
#         pathway_dir=path_dir_s,
#         data_type="smallRNA_SNF",
#         seed=RANDOM_STATE,
#     )

#     # ---------- 2.3 smallRNA–target–RNA overlap + triplets ----------
#     build_smallrna_rna_overlap(
#         rna_up=rna_up,
#         rna_down=rna_down,
#         targets_up=targets_up,
#         targets_down=targets_down,
#         small_up=small_up,
#         small_down=small_down,
#         mirna_db=mirna_db,
#         seed=RANDOM_STATE,
#     )

#     # ---------- 2.4 lncRNA：DEG + GO ----------
#     X_lnc, lab_snf_lnc, deg_dir_lnc, path_dir_lnc, targets_dir_lnc = load_expression_and_snf_clusters(
#         "lncRNA", RANDOM_STATE, cluster_source="SNF"
#     )
#     deg_lnc = compute_deg(X_lnc, lab_snf_lnc, group1=0, group2=1)
#     lnc_up, lnc_down, _ = save_deg_tables(
#         deg_lnc, deg_dir_lnc, "lncRNA_SNF", RANDOM_STATE
#     )

#     run_enrichr(
#         lnc_up,
#         label="Up-regulated lncRNA (SNF, relaxed)",
#         out_prefix="Pathway_up",
#         pathway_dir=path_dir_lnc,
#         data_type="lncRNA_SNF",
#         seed=RANDOM_STATE,
#     )
#     run_enrichr(
#         lnc_down,
#         label="Down-regulated lncRNA (SNF, relaxed)",
#         out_prefix="Pathway_down",
#         pathway_dir=path_dir_lnc,
#         data_type="lncRNA_SNF",
#         seed=RANDOM_STATE,
#     )

#     # ---------- 2.5 lncRNA -> miRNA ----------
#     mir_up, mir_down, lncrna_db = lncrna_to_mirna(
#         lnc_up, lnc_down, targets_dir_lnc, RANDOM_STATE
#     )

#     # ---------- 2.6 三层 ceRNA：lncRNA–miRNA–mRNA ----------
#     integrated_dir = os.path.join(BASE_DIR, "data", "integrated_results")
#     triplet_dir = os.path.join(integrated_dir, "triplets")
#     mirna_rna_trip_file = os.path.join(
#         triplet_dir, f"miRNA_target_RNA_triplets_seed{RANDOM_STATE}.csv"
#     )

#     build_ceRNA3_triplets(
#         lnc_up=lnc_up,
#         lnc_down=lnc_down,
#         mirna_up=mir_up,
#         mirna_down=mir_down,
#         mirna_rna_triplets_file=mirna_rna_trip_file,
#         lncrna_mirna_db=lncrna_db,
#         seed=RANDOM_STATE,
#     )

#     # ---------- 2.7 2 层 miRNA–mRNA triplets 后处理 + GO ----------
#     postprocess_mirna_mrna_triplets(seed=RANDOM_STATE)

#     print("\n[DONE] SNF-driven multi-omics DEG + pathway + ceRNA pipeline 完成。")


# if __name__ == "__main__":
#     main()

