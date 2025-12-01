# Before Start
1. **Original Aritcle analysis:** Since the original article has already been published in Nature Medicine, if our research merely "reproduces" its subtypes (Subtype 1 vs Subtype 2), it will not be able to be published in a top journal. We must "dig deeper," "subdivide more finely," or "find the mechanism more accurately" based on it.
2. **Key points:** Introduce the "non-coding RNA regulatory network" and "transcription factors as an upstream driving force", and combine them with **"spatial location"** to elevate the research topic.
# 第一部分：干实验——多组学精准分型与靶点发现 (Discovery)
1. **目标：** 利用公共数据，建立新的分型标准，构建调控网络，找到“主其事”的关键分子。
2. **Step 1: 多组学整合聚类 (Integrative Clustering)**
    - **方法：** 不仅使用mRNA，还要整合sRNA (miRNA) 数据。使用 **SNF (Similarity Network Fusion)** 算法==（*我不太懂，我随便写了个算法*）==进行无监督聚类。
    - **预期结果：** 定义 **4个** 稳定的分子亚型：==（举例）==
        1. **Cluster 1: 免疫-焦亡型 (Immuno-Pyroptotic):** 特征为 NF-κB/STAT3 激活，炎症小体 (NLRP3/GSDMD) 高表达。（_对应之前的Subtype 1，最严重_）
        2. **Cluster 2: 神经-突触剥离型 (Neuro-Synaptic):** 特征为神经肌肉接头 (NMJ) 破坏，MuSK 通路抑制，AchR 分解。（_这是新颖点_）
        3. **Cluster 3: 代谢-线粒体缺陷型 (Mito-Metabolic):** 特征为 OXPHOS 下调，线粒体自噬受阻，能量“断供”。（_对应之前的Subtype 2_）
        4. **Cluster 4: 泛素-蛋白酶体型 (Ubiquitin-Proteasomal):** 特征为 Atrogin-1/MuRF1 极端高表达，自噬过度激活，蛋白质“纯降解”。
        5. 任何可能的新结果，我们可以按照上面来解释
3. **Step 2: 构建亚型特异性 ceRNA 调控网络**
    - **方法：** 针对每个 Cluster，利用预测算法构建 **lncRNA - miRNA - mRNA** 互作网络。
    - **目的：** 解释为什么关键转录因子（如STAT3或FoxO）在特定亚型中被异常激活（例如，发现某个 lncRNA 可以降解抑制 STAT3 的 miRNA）。
4. **Step 3: 筛选核心驱动因子 (Master Regulators)**
    - **方法：** 使用 **SCENIC** 或 **DoRothEA** 算法分析转录因子活性。
    - **预期结果：** 锁定每个 Cluster 的“带头大哥”：==（举例）==
        - Cluster 1 Driver: **p-STAT3 / IRF1**
        - Cluster 2 Driver: **ATF4 / CHOP**
        - Cluster 3 Driver: **PGC-1α (缺失)**
        - Cluster 4 Driver: **FoxO3**
---
# **第二部分：湿实验——体内外验证与精准治疗 (Validation & Therapy)**
1. **目标：** 在临床样本和小鼠模型中验证这4种亚型的存在，并证明针对特定亚型的治疗有效。
2. **Step 1: 临床样本验证 (Human Validation)**
    - **方法：** 收集少量恶病质患者肌肉切片（或组织芯片 TMA）。
    - **操作：** 进行 **多重免疫荧光 (Multiplex IF)**，同时染色 4种亚型的标志物（例如：GSDMD, MuSK, PGC-1α, MuRF1, 这是我们上一波生信分析找到的结果）。
    - **逻辑：** 证明这4种病理状态在真实病人中确实存在，且不同病人以某种状态为主。这里就是整个研究最困难的地方，需要把前后衔接起来，我们如何把代码结果和实际的病人结果结合。
3. **Step 2: 小鼠模型异质性验证 (Mouse Heterogeneity)**
    - **模型：** 使用**UNKC6141 胰腺癌原位模型**。
    - **操作：** 在终点处死一批小鼠，**单独**分析每一只小鼠的肌肉转录组或蛋白组。
    - **逻辑：** 通过散点图证明，即使是同一种肿瘤，不同小鼠的肌肉也可能分布在 Cluster 1 到 Cluster 4 的不同状态，或者随病程进展发生状态转化（例如从 Cluster 3 代谢型 恶化为 Cluster 1 焦亡型）。
4. **Step 3: 功能性干预与药物筛选 (Intervention)**
    - **策略：** 选择最严重、临床最棘手的 **Cluster 1 (免疫-焦亡型)** 进行攻克。
    - **方法A (基因层面):** 使用现有的 **GSDMD-KO 小鼠** 构建模型。
        - _预期：_ Cluster 1 特征消失，肌肉萎缩逆转。
    - **方法B (药物层面):** 针对 Cluster 1 的上游驱动 TF（如 STAT3）或下游执行者（如 Caspase-1）给药。
    - **方法C (分子对接):** 针对第一部分发现的某个特定靶点（如某个 lncRNA 或 蛋白结构），利用计算机模拟筛选小分子药物，并在小鼠上测试。
