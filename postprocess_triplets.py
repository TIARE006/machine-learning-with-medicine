import os
import pandas as pd
import gseapy as gp

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SEED = 42

# === 目录结构 ===
INTEGRATED_DIR = os.path.join(BASE_DIR, "data", "integrated_results")
TRIPLET_DIR    = os.path.join(INTEGRATED_DIR, "triplets")
ENRICH_DIR     = os.path.join(INTEGRATED_DIR, "enrichment")
os.makedirs(TRIPLET_DIR, exist_ok=True)
os.makedirs(ENRICH_DIR, exist_ok=True)

# 原始 triplets 文件现在放在 triplets/ 下面
TRIPLET_FILE = os.path.join(
    TRIPLET_DIR,
    f"miRNA_target_RNA_triplets_seed{SEED}.csv"
)


def run_enrichr(gene_list, label, out_prefix):
    """对一组 mRNA 做 GO-BP 富集"""
    if not gene_list:
        print(f"[WARN] No genes for {label}, skip enrichment.")
        return

    print(f"[INFO] Enrichr for {label}, #genes = {len(gene_list)}")

    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=["GO_Biological_Process_2021"],
        cutoff=0.05
    )
    res = enr.results

    out_path = os.path.join(ENRICH_DIR, f"{out_prefix}_GO_BP_seed{SEED}.csv")
    res.to_csv(out_path, index=False)
    print(f"[SAVE] GO-BP result -> {out_path}")


def main():
    print(f"[INFO] Loading triplets: {TRIPLET_FILE}")
    df = pd.read_csv(TRIPLET_FILE)

    # 只保留真正参与调控的（两边都有方向信息）
    df = df[
        (df["miRNA_regulation"].isin(["up", "down"])) &
        (df["RNA_regulation"].isin(["up", "down"]))
    ].copy()

    print(f"[INFO] Total triplets with direction info: {df.shape[0]}")

    # 1) miRNA 上调，mRNA 下调（典型抑制型）
    up_miR_down_RNA = df[
        (df["miRNA_regulation"] == "up") &
        (df["RNA_regulation"] == "down")
    ].copy()

    # 2) miRNA 下调，mRNA 上调（释放抑制）
    down_miR_up_RNA = df[
        (df["miRNA_regulation"] == "down") &
        (df["RNA_regulation"] == "up")
    ].copy()

    # 保存详细表 → triplets/ 目录
    up_down_file = os.path.join(
        TRIPLET_DIR,
        f"triplets_miRup_RNAdown_seed{SEED}.csv"
    )
    down_up_file = os.path.join(
        TRIPLET_DIR,
        f"triplets_miRdown_RNAup_seed{SEED}.csv"
    )

    up_miR_down_RNA.to_csv(up_down_file, index=False)
    down_miR_up_RNA.to_csv(down_up_file, index=False)

    print(f"[SAVE] miR↑–RNA↓ triplets -> {up_down_file} (n={up_miR_down_RNA.shape[0]})")
    print(f"[SAVE] miR↓–RNA↑ triplets -> {down_up_file} (n={down_miR_up_RNA.shape[0]})")

    # 取各自的 mRNA 做 GO-BP 富集 → enrichment/ 目录
    genes_up_miR_down_RNA = sorted(up_miR_down_RNA["TargetGene"].unique())
    genes_down_miR_up_RNA = sorted(down_miR_up_RNA["TargetGene"].unique())

    run_enrichr(
        genes_up_miR_down_RNA,
        label="miR_up_RNA_down_targets",
        out_prefix="GO_miRup_RNAdown_targets"
    )

    run_enrichr(
        genes_down_miR_up_RNA,
        label="miR_down_RNA_up_targets",
        out_prefix="GO_miRdown_RNAup_targets"
    )

    print("[DONE] Triplets post-processing finished.")


if __name__ == "__main__":
    main()

