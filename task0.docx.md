**What Single-Cell Data Is Teaching Us About Cancer Evolution**

**Introduction**

Tumour evolution refers to the acquisition of genetic and epigenetic mutations by tumour cells, leading to the development of distinct subpopulations, or subclones. Each subclone contributes to the tumour's heterogeneity. (Andor et al., 2015) The evolution of cell populations into subclones is facilitated by genetic instability, epigenetic reprogramming and interaction of tumour cells with the microenvironment. (Ciriello et al., 2023)The evaluation of intra-tumour heterogeneity helps in disease prognosis. For example, pan-cancer analyses of whole-exome sequencing data revealed a poor prognosis for tumours with more than two detectable subclones. (Cosgrove et al., 2024)

Next-generation sequencing, including bulk DNA and RNA sequencing, has helped understand differential gene expression and identify biomarkers. However, these studies utilise bulk samples, such as biopsies, containing our populations mixed with host cells, such as oblasts and macrophages, which limits the analysis. Single-cell RNA sequencing allows evaluation of tumour heterogeneity at the cellular level. Over the years, several scRNA-Seq techniques have been developed for sequencing full-length mRNA (SMART-Seq2), 3’prime end mRNA (CEL-Seq2), 5’prime end mRNA (STRT-Seq) and small RNA (Holo-Seq). (Andor et al., 2015)

Although scRNA sequencing is highly useful in understanding tumour evolution and drug resistance, its clinical application remains a challenge due to inadequate sampling and complex computational analyses. (Mossner, Baker and Graham, 2021) Therefore, robust scRNA sequencing pipelines should be developed.  Nonetheless, scRNA studies have provided numerous insights into tumour evolution and response in breast, ovarian, brain, lung and colon cancer.

**Conceptual Framework: Linking Heterogeneity, Resistance and Evolution**

Tumour heterogeneity, a hallmark of tumour evolution, provides a fundamental basis for reconstructing the evolutionary history of tumours and for understanding their biological behaviour. It encompasses the diversity of cellular populations observed both between tumours of the same type in different individuals (intertumor heterogeneity) and within a single tumour (intratumor heterogeneity). This heterogeneity manifests through variations at multiple molecular levels, including genetic mutations, transcriptional alterations, protein expression changes, and epigenetic modifications (Sun and Yu, 2015). Tumour heterogeneity exists ubiquitously across all cancers. Chromosome mutational evolution, together with ecosystem pressures, leads to temporal and spatial segregation of mutation clones and subclones.Dentro et al. (2021)

Therapy resistance poses a significant obstacle to effective cancer treatment. Recent insights into cell plasticity as a new paradigm for understanding resistance to treatment: as cancer progresses, cancer cells experience phenotypic and molecular alterations, collectively known as cell plasticity. These alterations are caused by microenvironment factors, stochastic genetic and epigenetic changes, and/or selective pressure engendered by treatment, resulting in tumour heterogeneity and therapy resistance. Increasing evidence suggests that cancer cells display remarkable intrinsic plasticity and reversibly adapt to dynamic microenvironment conditions. Dynamic interactions between cell states and with the surrounding microenvironment form a flexible tumour ecosystem, which is able to quickly adapt to external pressure, especially treatment. Niu et al. (2024)

Mounting evidence suggests critical roles played by the tumour microenvironment (TME) in multiple aspects of cancer progression, particularly therapeutic resistance. The TME decreases drug penetration, confers proliferative and antiapoptotic advantages to surviving cells, facilitates resistance without causing genetic mutations and epigenetic changes, collectively modifying disease modality and distorting clinical indices. Sun (2015)

As an emerging sequencing technology, single-cell RNA sequencing (scRNA-Seq) has become a powerful tool for describing cell subpopulation classification and cell heterogeneity by achieving high-throughput and multidimensional analysis of individual cells and circumventing the shortcomings of traditional sequencing for detecting the average transcript level of cell populations. Wang et al. (2023)

For example, cytokines secreted by tumour cells, such as interleukins, can recruit macrophages, and the secretion of VEGF can stimulate the migration and proliferation of endothelial cells and promote tumour angiogenesis. Zhou JX et al. systematically compared the expression patterns of 2,558 ligand receptor pairs in seven cell types isolated from melanoma (melanoma cells, T cells, B cells, macrophages, NK cells, CAF cells, and endothelial cells) by scRNA-seq. Based on the specific ligand-receptor connection, the interaction network between different cell types was established. A large number of detected cell communication signals, such as growth factors, chemokines and matrix proteins, were reported to be closely related to tumour proliferation, metastasis, cell adhesion, angiogenesis, and immunoregulatory processes, thereby revealing that the interaction network was of great significance for tumour development and prognosis (Watson et al., 2013).

**Technological Evolution of scRNA-seq**

Single-cell RNA sequencing (scRNA-seq) is a transformative technology that profiles transcriptomes at the resolution of individual cells, enabling the identification of cell-specific gene expression patterns, epigenetic changes, and rare or previously unrecognised cell types. The technique addresses key challenges inherent to single-cell analysis, including the isolation of individual cells and the amplification of minute amounts of RNA for library construction without introducing significant bias *(Wang, S. et al. 2023)*.

Early scRNA-seq methods, such as Smart-seq and Smart-seq2, provided full-length transcript coverage with high sensitivity, allowing for detailed studies of rare cells and splice variants. However, these protocols were labour-intensive, low-throughput, and costly. Subsequent innovations, including Smart-seq3 with 5′-end counting and unique molecular identifiers (UMIs), improved quantification and reduced amplification bias *(Tieng FYF et al, 2025)*.

The introduction of droplet-based microfluidic platforms, such as 10x Genomics' Chromium, has enabled the high-throughput profiling of thousands to tens of thousands of cells simultaneously. These methods combine single-cell encapsulation, barcoding, and cDNA amplification in a scalable workflow, dramatically increasing efficiency and reducing cost per cell *(Wang, S. et al. 2023)*.

Following this, multi-omics single-cell approaches emerged, integrating transcriptomic data with genomic, epigenomic, and proteomic information to study multiple molecular layers within the same cell, providing deeper insights into cellular heterogeneity and regulatory programs *(Catalano et al,.2025)*. Meanwhile, spatial transcriptomics methods, such as Slide-seq and MERFISH, restored tissue context lost in dissociative workflows, revealing how transcriptional programs are organised within tissues and cellular microenvironments *(Williams et al,.2022)*.

Together, these innovations, high-sensitivity protocols, high-throughput droplet methods, multi-omics integration, and spatial resolution have transformed single-cell biology and are now widely applied in cancer research to dissect tumour heterogeneity, trace clonal evolution, and identify rare cell populations that drive disease progression.

**Key Computational Advances**

Alongside experimental progress, powerful computational methods have greatly advanced the interpretation of scRNA-seq data.The first critical step is **unsupervised clustering analysis**, which groups cells based on transcriptional similarity to uncover distinct cell types and tumour subpopulations *(Xiang et al, 2022)*. Numerous clustering algorithms, such as Louvain and Leiden, have been optimised for high-dimensional single-cell data, improving scalability and robustness *(Yu, L.2022)*. Once clusters are identified, marker gene expression is used to annotate and characterise cell identities, providing a foundation for downstream biological interpretation *(Xiang et al,.2022)*.

Single-cell RNA sequencing (scRNA-seq) has revolutionised cancer research by enabling transcriptome profiling at single-cell resolution, revealing tumour heterogeneity and rare malignant clones. Yet, it provides only a static view of gene expression. Computational tools such as pseudotime analysis and RNA velocity now reconstruct dynamic cellular trajectories, predicting how cells transition between states during tumour progression or treatment response. These methods extend scRNA-seq beyond static profiling, capturing the temporal and evolutionary dynamics that drive cancer development *(Teppei 2025)*.

Ligand–receptor inference represents another key computational advance in scRNA-seq analysis. By quantifying ligand and receptor gene expression across cell types, these models reconstruct intercellular communication networks within complex tissues. In cancer research, this approach reveals how tumour, immune, and stromal cells coordinate invasion, angiogenesis, and immune evasion through signalling crosstalk *(Vahid et al,.2023)*. Tools such as CellPhoneDB and NicheNet have enabled systematic mapping of these interactions, uncovering potential therapeutic targets that disrupt pro-tumour signalling pathways and reshape the tumour microenvironment *(Yang et al,.2024)*.

Machine learning and AI have become essential for analysing the high-dimensional and complex data generated by scRNA-seq. These approaches enable efficient clustering, dimensionality reduction, trajectory inference, and cell type annotation, often combining traditional methods like support vector machines with advanced architectures such as graph neural networks and transformers *(Liu et al,.2025)*. In cancer research, ML-driven analyses identify rare malignant or immune cell populations, reconstruct dynamic cell-state transitions, and predict responses to therapy. By automating pattern recognition and refining predictions of cellular behaviour, AI enhances both the accuracy and scalability of single-cell studies, bridging computational innovation with biological discovery *(Liu et al,.2025)*.

**Table 1: Compare Main Scrna-Seq Platforms And Integration Methods (10x Genomics, Slide-Seq, Merfish, etc).**

|  |  |  |  |  |  |
| --- | --- | --- | --- | --- | --- |
| **Platform/ Method** | **Types** | **Key Features** | **Strengths** | **Limitations** | **References** |
| 10x Genomics Chromium | Droplet-based scRNA-seq | Microfluidic droplet encapsulation, barcoding, UMIs | High throughput (up to 10,000 cells per run), cost-effective, well-supported | Fresh samples needed  Detect only 10% mRNA | *(Xiliang et al,. 2021;Gao et al,.2020 )* |
| Smart-seq2 / Smart-seq3 | Plate-based scRNA-seq | Full-length transcript coverage, UMIs in Smart-seq3 | High sensitivity, captures splice variants, ideal for rare cells | Low throughput, labour-intensive, expensive per cell | *(Tieng FYF et al,.2025)* |
| Merfish | Spatial transcriptomics | Multiplexed error-robust FISH, single-molecule resolution | High spatial resolution, measures hundreds to thousands of transcripts simultaneously with subcellular resolution1 | Low signal-to-noise ratio of single-molecule signals and big  data interpretation challenge. | *(Liu J et al,. 2022 ; Ghani et al,.2019)* |
| Seq-Well / inDrop | Spatial transcriptomic | Microwell/droplet encapsulation, barcoding | Scalable, sensitive, and cost-efficient | Lower sensitivity than plate-based methods, some cell loss | *(Simonas et al,.2025; De et al,. 2023)* |
| Slide-seq | Droplet / microwell-based scRNA-seq | Bead-based spatial barcoding | Preserves tissue architecture, maps gene expression spatially | Lower sensitivity per gene, limited transcript coverage | **(***De et al,. 2023;Rodriques et al,.2019 )* |

**Mapping Tumour Heterogeneity through scRNA-seq**

As single-cell RNA sequencing technologies are maturing, researchers are now able to observe cellular heterogeneity within tumours and leukocytes with as much detail as possible, linking distinct cell states, microenvironment interactions and evolutionary trajectories. With bulk RNA sequencing, the diversity of cell populations within each tumour was not visible. In contrast, scRNA-seq dissects the tumours into individual cells, enabling detailed mapping of tumour and non-tumour compartments, and the recognition of sub-clonal, plastic and rare cell populations that drive progression, resistance to therapy and recurrence.

**Understanding the dimensions of heterogeneity**

Tumour heterogeneity comprises multiple orthogonal layers.

* Genetic (clonal) heterogeneity: different subclones within the tumour carry distinct somatic mutations or copy-numbered alterations, reflecting an evolutionary tree of diverging lineages.
* Transcriptomic (cell state) heterogeneity: cells with the same genetic background may adopt different gene expressions, for example, proliferative vs quiescent, epithelial vs mesenchymal, immunologically cold vs immunologically hot.
* Phenotypic & microenvironmental heterogeneity: the different cell types and the way they interact create spatial, temporal and functional diversity within the tumour ecosystem. scRNA provides a microscope on both the cellular composition and the state dynamics of the tumour, enabling links to therapy resistance, relapse and evolutionary adaptation.

**How scRNA-seq maps heterogeneity**

When using scRNA-seq in tumour studies, the tumour tissue is dissociated into single cells, their transcriptomes are then captured and analysed by dimensionality reduction (UMAP/t-SNE), clustering, and trajectory inference. The result we get is a map of cell populations: malignant epithelial clusters, immune and stromal clusters, and rare subsets. Computational analysis may help us see pseudotime of cell states, detect copy-number variation signatures in scRNA data to infer subclones, and compute cell-cell interaction networks via ligand-receptor pairing. Through this, we can identify:

* distinct subclones of cancer cells
* plastic transitions
* rare or emergent populations
* diverse microenvironment niches (immune, stromal, vascular) that correlate with tumour behaviour.

**Key findings across tumour types**

Here we summarise representative studies (see Table 1) and highlight selected findings:

* In breast cancer, Xu et al. identified six distinct natural-killer (NK) cell subtypes in the tumour micro-environment, and malignant epithelial cells that exhibited variation by gene expression, intrinsic subtype and intratumoral programs. This shows that what might be considered a uniform NK population is, in fact, composed of multiple subsets differentiated by functionality and transcriptomics. Similarly, malignant epithelial cells are not a homogeneous block within the tumour, but differ in state, which underlies transcriptional heterogeneity and potentially differential response to therapy.
* In ovarian cancer, studies such as Zhang et al. revealed that high stromal tumours were enriched for CXCL12-expressing cells and immune-excluded NK and CD8 T cells, while Chai et al. identified early surges of immunosuppressive Tregs and exhausted CD8 T cells and found antigen-presenting CAFs interacting with tumour cells via the NECTIN2-TIGIT axis.
* In pancreatic ductal andecarcinoma (PDAC), Elyada et al. found that although the tumour epithelial subtype remained mixed even after chemotherapy, the TME was profoundly remodelled, including changes in TIGIT-NECTIN2 signalling. Zhou et al. integrated spatial and single-cell data and found chemotherapy-resistant tumours enriched for metallothionein-positive CAFs and distinct malignant subclusters in particular niches.
* In lung cancer (NSCLC and lung adenocarcinoma), Wu et al. generated a single-cell atlas of advanced tumours showing that a high transcriptional heterogeneity correlates with tumour-associated neutrophil infiltration, while Bischoff et al. showed two distinct immune-stromal phenotypes that categorised patient outcome.
* In leukaemia, Mumme et al. in paediatric AML tracked blast clusters and immune niche signatures, predictive of relapse, and Bhasin et al. in T-ALL showed that minimal residual disease blasts had up-regulated PD-1 signalling and distinct cell-cell communication networks.

**Biological and clinical implications**

The mapping of tumour heterogeneity via scRNA-seq has several major implications:

* Predicting progression and metastasis: Identifying transcriptionally distinct sub-clones allows prediction of which tumours may recur or metastasise.
* Understanding therapy escape and resistance: Rare or plastic tumour sub-states can be identified before treatment; similarly, distinct immune/stromal niches may foster resistance.
* Enabling precision medicine: Single-cell maps can identify cell-type-specific biomarkers and stratify patients by tumour ecosystem.
* Targeting the microenvironment: heterogeneous stromal/immune compartments revealed by scRNA can serve as a therapeutic target to remodel or normalise the tumour.

By resolving the complexity of tumour and micro-environmental cell populations, scRNA-seq reveals not only the map of heterogeneity, but it also reveals the path by which tumours evolve and resist therapy.

**Mechanism of Therapy Resistance revealed by scRNA-seq**

Several types of therapies are used in the treatment of cancer. The choice of therapy is dependent on the type, location, and stage of the disease. Some patients receive one type of treatment, while most receive combination therapies (*Types of Cancer Treatment*, n.d.). The current types of therapies include surgery, immunotherapy, hyperthermia, stem cell transplant, chemotherapy, photodynamic therapy, radiation therapy, hormone therapy and targeted therapy (Kaur et al., 2023). Despite the advancements that have been made in treating cancer using these therapies, scientists are faced with a limitation of therapy resistance, which tends to be either intrinsic or acquired. Intrinsic resistance occurs when cancer cells fail to respond to therapy from the beginning. In the case of acquired resistance, however, the cells respond positively to treatment but later become unresponsive. The diverse and heterogeneous nature of the disease is the primary trait that makes it resistant to therapy, fuels its progression and recurrence (Damodaran et al., 2025).

Increased drug efflux by ATP-binding cassette (ABC) transporters, drug target mutations, increased DNA repair capacity, resistance to apoptosis, and dysregulated signalling pathways are some examples of mechanisms of drug resistance. These mechanisms, whether intrinsic or acquired, allow cancer cells to endure and proliferate under the selective pressures of chemotherapy, targeted therapy, or even immunotherapy, ultimately recapitulating the efficacy of anticancer drugs. Some of these mechanisms are discussed in detail:

**Epithelial-Mesenchymal Transition (EMT)**

Epithelial-mesenchymal transition is an important process which occurs in tissue regeneration and normal embryonic development. However, the abnormal reactivation of EMT causes epithelial cells to take on mesenchymal properties and become more motile, resulting in the emergence of an invading phenotype, enhanced resistance to immunotherapy and chemotherapy(Huang et al., 2022; Liaghat et al., 2024).

**Metabolic reprogramming**

Metabolic reprogramming, a characteristic feature of cancer, alters the glucose, lipid and amino acid metabolism, enabling tumour cells to adapt to their environment (Liu et al., 2024). Its orchestration is caused by a complex interaction between genetic and protein-level modifications. Metabolites play a crucial role in the tumour microenvironment (TME) by fuelling energy production, biosynthetic precursors, or waste products, which alter cell function. Normal cells and malignant cells have a significant metabolic difference, because tumour cells tend to sustain growth and acquire invasive capabilities simultaneously. For instance, in the Warburg effect, tumours show enhanced aerobic glycolysis, converting sugars into lactic acid even in the presence of oxygen, at higher levels compared to those in normal tissues (DeBerardinis & Chandel, 2020).

**Epigenetic plasticity**

DNA methylation and histone alterations are two of the ways that epigenetics alters genes. In several processes, including the induction and maintenance of senescence, which are critical stages for slow-cycling cancer cells, epigenetic regulation plays a crucial role in phenotypic plasticity and influences cell destiny. A stress-induced premature senescence (SIPS) phenotype of non-cancerous cells is similar to the stress-induced, slow-cycling subset of cancer cells. While replicative senescence is governed by heritable DNA methylation alterations, including global hypomethylation of DNA and focal DNA hypermethylation, SIPS is regulated by temporal histone modifications rather than DNA methylation mechanisms (Menon et al., 2020).

**Drug-tolerant persister cells**

A small population of cancer cells can escape cell death resulting from targeted therapy and chemotherapy by entering a drug-tolerant persistent state. Tumour cells can stay in this state long enough to develop additional mechanisms of acquired drug resistance (Li et al., 2025).

**Immune evasion**

Immune evasion allows tumours to escape the surveillance and destruction by immune cells, leading to poor patient outcomes and complicating therapeutic interventions. Tumour-induced immune suppression is one way of immune evasion. This entails creating an immunosuppressive milieu that impedes the growth and function of immune cells, enabling tumour cells to endure and proliferate unchecked. Tumours can specifically affect the immune system by the expression of checkpoint molecules that impede immune responses, the secretion of immunosuppressive substances, and the recruitment of regulatory immune cells. Immune checkpoint pathways, which are typically involved in self-tolerance and autoimmune regulation, can also be used by tumour cells to avoid immune detection (Tufail et al., 2025).

**Tumour Microenvironment (TME) and Cell - Cell Interactions**

A tumour does not exist in isolation; it grows within a complex neighbourhood of non-cancerous cells collectively called the tumour microenvironment (TME), an intrinsic milieu within which tumour cells thrive *(Li, Z., et al., 2025)*.

![](data:image/png;base64...)

*Figure 1. Composition of the tumour microenvironment (Li, Z., et al., 2025)*

This diverse cellular ecosystem is made up of stromal cells such as cancer-associated fibroblasts (CAFs) and pericytes, that remodel tissue; immune cells like tumour-associated macrophages (TAMs), T cells, B cells, and natural killer (NK) cells that defend or at times protect the tumour; and vascular cells, like endothelial cells that build new blood vessels *(Wang et al., 2021)*. This shows that the tumour cannot survive or adapt on its own; it needs and uses the surrounding cells. Together, these cell types create a very dynamic microenvironment where they constantly exchange molecular messages, growth factors, cytokines, and chemokines, much like residents negotiating in a noisy city square. Some defend the healthy tissues; others unintentionally help the tumour to survive and resist treatment *(Perez-Moreno, 2009)*. This is why understanding how these cells interact with each other and the tumour cells is crucial in developing new and improved cancer treatments.

The big question is *‘How can we eavesdrop on the crosstalk, the conversation that keeps cancer alive, to better understand the TME?’* The answer is Single-cell RNA sequencing. scRNA-seq enables high-resolution profiling of individual cells, which allows us to map which genes are active in each cell type *(Li et al., 2022)*. It exposes who’s talking to whom by pairing ligand receptor expression patterns. This reveals crucial communication circuits such as:

* CXCL12-CXCR4: CAFs luring immune cells away from the core of the tumour *(Peng et al., 2022)*.
* TGFB-TGFBR: stromal signalling that triggers EMT and fibrosis *(Timperi & Romano, 2023)*.
* PD-L1-PD-1: a tumour-immune crosstalk that silences T-cell attack.

These findings show that the tumour somehow understands that it takes a community to survive and thrive. However, thanks to scRNA-seq, we can target these crosstalk pathways to disrupt the supportive interactions between tumour cells and their microenvironment. *(Li et al., 2022).*

**Cancer-associated fibroblasts**

Among the TME diverse cells, CAFs are what we would call the expert builders of the tumour’s protective architecture. They are originally responsible for wound repair; however, they are co-opted by cancer cells to remodel the ECM through MMP‐mediated degradation and collagen crosslinking, thereby creating a fibrotic niche  *(Miao et al., 2024)*, secrete growth factors, and manipulate immune cells, thereby turning tissue repair into tumour maintenance *(Li et al., 2025)*. The perfect puppet for cancer cells. The secreted factors include TGFβ, PDGF, FGF, hepatocyte growth factor (HGF), VEGF, tumour necrosis factor α (TNFα), interferon-γ (IFNγ), CXCL12, IL-6, connective tissue growth factor (CTGFβ), EGF, growth arrest-specific protein 6 (GAS6), galectin-1, secreted frizzled-related protein 1 (SFRP1), sonic hedgehog protein (SHH), and bone morphogenetic protein (BMP), which are tumour-promoting *(see Figure 2)* *(Valkenburg et al,.2018)*.

![](data:image/png;base64...)

*Figure 2. Schematic representation of selected pro-tumorigenic functions of CAFs. CAFs induce (****1****) angiogenesis and tumour growth, (****2****) invasion and metastasis of cancer cells, (****3****) modulation of the immune system, including recruitment and activation of immune suppressors and inhibition of anti-tumour effector cells, and (****4****) therapy-resistance through ECM production and remodelling (Glabman et al,.2022).*

By using scRNA-seq, scientists have uncovered that CAFs are not a uniform population, but an entire team of specialised subtypes:

* CXCL12+ CAFs: these function as molecular gatekeepers, secreting chemokines that bind CXCR4 receptors on immune cells, thereby creating immune-excluded zones where cytotoxicity T cells cannot penetrate (Shah et al., 2022).
* LRRC15+ CAFs: upregulate TGF-β signalling, leading to dense fibrotic barriers that block drug diffusion and encourage epithelial-mesenchymal transition (EMT) in tumour cells  (Cords et al., 2024).

By using scRNA-seq, we can track how these complex CAF populations emerge under therapeutic stress by influencing drug resistance and tumour progression. An eavesdrop into how these cellular construction workers rebuild the tumour’s defences every time we try to demolish it (Ma et al., 2023).

**Immune response exhaustion**

While CAFs are busy remodelling the tumour environment, the immune cells are most times left exhausted and unarmed. This exhaustion from chronic antigen exposure, expressing markers like PD-1 and TOX that signal burnout, hinders the T cells' ability to effectively combat the tumour cells (“Deciphering the Tumour Immune Microenvironment Cell by Cell,” 2023). Moreover, macrophages often switch to the M2 state, which makes them operate more like peacekeepers than soldiers, suppressing inflammation and secreting pro-growth signals. Even dendritic cells lose their ability to alter T cells effectively  (Pattabiram et al., 2025) and contribute to the immune evasion of the tumour. Thanks to scRNA-seq, these interactions can be decoded at single-cell resolution, hence identifying the exact ligand-receptor dialogue that influences immune suppression.

As scientists, scRNA-seq allows us to hear the dialogues of how the cells in the TME communicate with and influence each other and thereby influencing cancer progression and resistance to treatment. Stromal cells, such as CAFs and pericytes, remodel the physical and chemical environment. In contrast, immune cells, including TAMs, T cells, and dendritic cells, negotiate ceasefires through checkpoint signalling. Endothelial cells maintain the nutrient supply for the tumour, and the tumour cells adapt to these cues from their supportive neighbours, thereby evolving resistant phenotypes. So, resistance is an adaptation emerging from collective responses of different cell types communicating through molecular channels. Even better, when scRNA-seq is paired with spatial transcriptomics, not only can scientists eavesdrop on the crosstalk, but they can see where it is taking place.

**Integrating with Multi-Omics and Spatial Data**

Tissue dissociation is a crucial step in scRNA analysis; however, it leads to the loss of spatial information of interactions between different cell types. (Yan et al., 2024) Spatial omics is based on counting transcripts of a gene at a distinct tissue location. In other words, it utilises transcriptomics analysis to analyse expression patterns of genes and cells while maintaining tissue integrity. (Williams et al., 2022)  However, spatial transcriptomics does not have the high resolution of single-cell sequencing. To understand the role of different cell types within a tissue structure, several integration methods have been developed to combine scRNA sequencing and spatial transcriptomics.

In their study on the *‘Spatial single cell analysis of tumour microenvironment remodelling pattern in primary central nervous system lymphoma’*, Xia Yuan et al. performed spatial transcriptomics and matched the corresponding single-cell sequencing data of PCNSL patients.

![](data:image/png;base64...)

*Figure 3. Landscape of the primary central nervous system lymphoma (PCNSL) microenvironment.*

It was found that tumour cells may achieve a “TME remodelling pattern” through an “immune pressure-sensing model”, in which they could choose to reshape the TME into a barrier environment or a cold environment according to the immune pressure. From the integration of spatial transcriptomics with scRNA-seq, they were able to discover the spatial and temporal distribution of and characteristics of immune checkpoint molecules and CAR-T target molecules in immunotherapy through the spatial communication analysis. The data from their study clarified the tumour microenvironment remodelling pattern of PCNSL and provided a reference for immunotherapy treatment options (Xia, Y.,2025).

Integration methods of scRNA sequencing and spatial omics are divided into two categories: deconvolution methods and mapping methods. Deconvolution is applied to spatial barcoding data with scRNA seq data as background, whereas mapping uses high-plex RNA imaging data to localise scRNA-seq subpopulations. Mapping is more flexible as it doesn’t require previously developed cell subtype models. Spatial barcoding is based on transcriptome analysis of mRNA tagged with a barcode in arrays of beads on a tissue slide, whereas high plex RNA imaging is based on labelled imaging methods like FISH, which provides higher resolution up to the subcellular level. The successful integration of scRNA analysis and spatial transcriptomics requires standardisation across sample collection, sequencing technique and data analysis. Since there are numerous techniques for integration, the selection should be dependent on a case-by-case basis.

In terms of integration with other single cell omics technologies for better understanding of cell functions and disease resistance mechanisms, researchers have developed AI models for integration using feature links, for example scMODEL developed by Wang et al(2025), utilised neural networks and deep learning algorithms to integrate single cell omics data with common features, in tonsil tissues, removing variation and preserving significant biological data. (Wang, 2025)

**Translation and Clinical Implications**

The integration of single-cell RNA sequencing (scRNA-seq) into clinical oncology is beginning to redefine how we discover biomarkers and design therapies. By resolving the tumour microenvironment (TME) at cellular resolution, scRNA-seq reveals hidden subpopulations that predict treatment outcomes and disease progression. It also allows researchers to dissect the cellular and molecular heterogeneity of tumours at single-cell resolution, revealing unique transcriptional states that bulk sequencing often masks. These fine-grained profiles have led to the identification of predictive biomarkers—cell-type–specific gene signatures that correlate with therapy outcomes. For example, Chen et al. (2025) used single-cell and bulk RNA sequencing to identify a PRRX2-driven fibroblast signature in colorectal cancer that predicts poor response to immune checkpoint inhibitors and may guide combination therapies targeting TGF-β signalling *(Chen et al., 2025)*. Similarly, single-cell profiling of bladder cancer uncovered a CD6/ALCAM-mediated immune-suppressive pathway associated with BCG therapy resistance, offering a target for improving immunotherapy efficacy *(Jurić et al., 2025)*. Also, single-cell signatures have been linked to response or resistance to immunotherapy and chemotherapy (Sinha et al., 2022). Such biomarkers can help anticipate treatment efficacy and guide patient-specific therapeutic strategies. These examples highlight how scRNA-seq-derived biomarkers can support both prognostic and predictive applications in precision oncology.

Patient stratification represents another major translational opportunity. By characterising immune and stromal heterogeneity, scRNA-seq enables clinicians to distinguish responders from non-responders based on TME composition. By characterising the tumour microenvironment (TME) at single-cell depth, scRNA-seq helps classify patients based on immune landscape features—such as T-cell exhaustion states, myeloid infiltration, or stromal interactions. In breast cancer, Peng et al. (2025) integrated single-cell and bulk RNA data to develop a cuproptosis-related gene signature that accurately stratifies patients by risk and immune landscape, offering a framework for personalised immunotherapy planning *(Peng et al., 2025)*. Integrating these single-cell insights into clinical workflows could enable more accurate selection of patients likely to benefit from immune checkpoint inhibitors or combination therapies. Similar multi-omic analyses in colorectal and cervical cancers have revealed cell-type-specific markers, such as neutrophil-linked PRDX1 and FTH1 genes or CCR7+ CD8⁺ T cells that could refine patient selection for targeted and immune-based therapies *(Yan et al., 2025; Yuan et al., 2025)*.

Bridging discovery and therapy requires functional validation and computational translation. Despite its promise, clinical implementation of scRNA-seq faces major barriers. The technology remains costly, data processing pipelines lack standardisation, and regulatory guidelines for clinical-grade single-cell analyses are still emerging. Patient-derived organoids and xenografts now serve as living testbeds to confirm scRNA-seq predictions and screen for drug sensitivity, while AI-driven platforms use transcriptomic signatures to forecast optimal drug combinations *(Rafique et al,.2021)*. Yet, significant barriers remain. High costs, data standardisation issues, and the lack of unified regulatory frameworks continue to slow clinical adoption. Ensuring data reproducibility, establishing validated biomarkers, and developing regulatory frameworks for clinical use are essential next steps. Overcoming these hurdles will demand interdisciplinary collaboration among computational scientists, clinicians, and regulatory agencies. If achieved, single-cell-informed precision oncology could move beyond describing cellular diversity to harnessing it, transforming every patient’s tumour into a blueprint for their personalised treatment strategy.

**Future Directions and Conclusion**

Single-cell technologies have changed how we think about cancer. We now see tumours as a complex ecosystem with a dynamic community of tumour, immune and stromal cells locked in constant communication. scRNA-seq has given us the ability to listen in on this conversation, showing us how cancers evolve, adapt, and resist treatment at a newer resolution. However, future work needs to move from describing what we see in tumours to predicting and influencing how they change, and in order to do this, there are a few key challenges to be dealt with.

Firstly, there is a need for standardised and reproducible workflows from tissue handling to computational analysis to ensure data consistency across studies and cancer types. Presently, different labs use slightly different ways to collect and process cells, which makes it harder to compare results. Shared protocols and global databases will help build a comprehensive, reliable picture of how different cancers behave across populations and treatment contexts.

Secondly, another challenge lies in watching tumours evolve. Cancer isn’t static; it adapts and changes under the pressure of treatment, think of it as a city reshaping after a storm. By combining scRNA-seq with other multi-omics methods, such as spatial mapping, proteomics and epigenomics, scientists can track not just what genes are active but where and when key changes occur inside a tumour. This type of ‘4D view is beginning to show us how resistant cell populations emerge, migrate and remodel their surroundings.

Thirdly, with AI entering the scene redefining how we interpret single-cell data. New deep learning models like scGen and MultiVI can integrate multiple layers of molecular information to predict how individual tumour cells might respond to therapy, offering a peek into a future where clinicians can model treatment outcomes before giving a drug *(Lotfollahi et al, 2023; Ashuach et al, 2023)*.

Finally, *‘How can we translate these discoveries into patient care?’* For this to happen, there is a need for functional validation through organoids, xenografts and clinical trials. Single-cell insights should guide biomarker discovery, patient stratification and drug combination design, ultimately informing real-time treatment decisions in the clinic.

In conclusion, single-cell RNA sequencing has changed the game by exposing the molecular dialogues driving adaptation and resistance. scRNA-seq has opened new boundaries in precision oncology. As technologies intersect, integrating spatial, temporal and functional elements, single-cell approaches are set to change oncology from static observation to predictive and intervention-driven medicine.

**References**

Andor, N., Graham, T.A., Jansen, M., Xia, L.C., Aktipis, C.A., Petritsch, C., Ji, H.P. and Maley, C.C. (2015). Pan-cancer analysis of the extent and consequences of intratumor heterogeneity. *Nature Medicine*, 22(1), pp.105–113. doi:<https://doi.org/10.1038/nm.3984>.

Ashuach, T., Lopez, R., et al. (2023). MultiVI: Integrative deep probabilistic modelling for single-cell multi-omics. Cell Systems, 15(3), 240–256. DOI:10.1038/s41592-023-01909-9

Bhasin, S.S., Thomas, B.E., Summers, R.J., Sarkar, D., Mumme, H., Pilcher, W., Emam, M., Raikar, S.S., Park, S.I., Castellino, S.M., Graham, D.K., DeRyckere, D. and Bhasin, M.K. (2023). Pediatric T-cell acute lymphoblastic leukaemia blast signature and MRD-associated immune environment changes defined by single-cell transcriptomics analysis. *Scientific Reports*, 13(1), 12556. doi:10.1038/s41598-023-39152-z.

Catalano M, D'Angelo A, De Logu F, Nassini R, Generali D, Roviello G. Navigating Cancer Complexity: Integrative Multi-Omics Methodologies for Clinical Insights. Clin Med Insights Oncol. 2025 Oct 21;19:11795549251384582. Doi: 10.1177/11795549251384582. PMID: 41147019; PMCID: PMC12553891.

Chen, M. et al. (2025). Single-cell and bulk RNA-sequencing reveal PRRX2-driven cancer-associated fibroblast-mediated perineural invasion for predicting the immunotherapy outcome in colorectal cancer. Frontiers in Cell and Developmental Biology. https://doi.org/10.3389/fcell.2025.1620388

Ciriello, G., Magnani, L., Aitken, S.J., Akkari, L., Behjati, S., Hanahan, D., Landau, D.A., Lopez-Bigas, N., Lupiáñez, D.G., Marine, J.-C., Martin-Villalba, A., Natoli, G., Obenauf, A.C., Oricchio, E., Scaffidi, P., Sottoriva, A., Swarbrick, A., Tonon, G., Vanharanta, S. and Zuber, J. (2023). Cancer Evolution: A Multifaceted Affair. *Cancer Discovery*, 14(1), pp.OF1–OF13. doi:<https://doi.org/10.1158/2159-8290.CD-23-0530>.

Cords, L., de Souza, N. and Bodenmiller, B. (2024). Classifying cancer-associated fibroblasts—The good, the bad, and the target. *Cancer Cell*. Doi:<https://doi.org/10.1016/j.ccell.2024.08.011>.

Cosgrove, P.A., Bild, A.H., Dellinger, T.H., Badie, B., Portnow, J. and Nath, A. (2024). Single-Cell Transcriptomics Sheds Light on Tumour Evolution: Perspectives from City of Hope’s Clinical Trial Teams. *Journal of Clinical Medicine*, 13(24), 7507. doi:<https://doi.org/10.3390/jcm13247507>.

De Jonghe, J., Kaminski, T.S., Morse, D.B. et al. spinDrop: a droplet microfluidic platform to maximise single-cell sequencing information content. Nat Commun 14, 4788 (2023). <https://doi.org/10.1038/s41467-023-40322-w>

Deciphering the tumour immune microenvironment cell by cell. (2023). *Immuno-Oncology Technology.* Doi:<https://doi.org/10.1016/j.iotech.2023.100383>.

Dentro, S.C. et al. (2021). Characterising genetic intra-tumour heterogeneity across 2,658 human cancer genomes. *Cell*, 184(8), pp.2239–2254.e39. doi:10.1016/j.cell.2021.03.009.

Gao C, Zhang M, Chen L. The Comparison of Two Single-cell Sequencing Platforms: BD Rhapsody and 10x Genomics Chromium. Curr Genomics. 2020 Dec;21(8):602-609. doi: 10.2174/1389202921999200625220812. PMID: 33414681; PMCID: PMC7770630.

Ghani Naureen and Yoh Isogai, Optimisation of MERFISH Data Analysis and Visualisation, May 21, 2019. <https://www.ucl.ac.uk/~ucbpngh/MiniProject3.pdf>.

Glabman, R. A., Choyke, P. L., & Sato, N. (2022). Cancer-Associated Fibroblasts: Tumorigenicity and Targeting for Cancer Therapy. Cancers, 14(16), 3906. https://doi.org/10.3390/cancers14163906

Jurić, I. et al. (2025). Single-cell RNA-sequencing of BCG naïve and recurrent non-muscle invasive bladder cancer reveals a CD6/ALCAM-mediated immune-suppressive pathway. npj Precision Oncology. https://doi.org/10.1038/s41698-025-01093-3

Leader, A.M., Desai, J., Huang, H. et al. (2021). Single-cell analysis of human non-small cell lung cancer lesions refines tumour classification and patient stratification. *Cancer Cell*, 39, pp.1594–1609.e12. doi:10.1016/j.ccell.2021.10.009.

Li, Z., Li, J., Bai, X. et al.Tumour microenvironment as a complex milieu driving cancer progression: a mini review. Clin Transl Oncol 27, 1943–1952 (2025). <https://doi.org/10.1007/s12094-024-03697-w>

Li, P., Kong, X., He, Y., Liu, Y., Peng, X., Li, Z.-H., Xu, H., Luo, H. and Park, J. (2022). Recent developments in the application of single-cell RNA sequencing in the tumour immune microenvironment and cancer therapy. *Military Medical Research.* Doi:<https://doi.org/10.1186/s40779-022-00414-y>.

Li, Y., Liu, Q., Jing, X., Wang, Y., Jia, X., Yang, X. and Chen, K. (2025). Cancer-Associated Fibroblasts: Heterogeneity, Cancer Pathogenesis, and Therapeutic Targets. *MedComm*, 6(7), e70292. doi:<https://doi.org/10.1002/mco2.70292>.

Liu J, Tran V, Vemuri VNP, Byrne A, Borja M, Kim YJ, Agarwal S, Wang R, Awayan K, Murti A, Taychameekiatchai A, Wang B, Emanuel G, He J, Haliburton J, Oliveira Pisco A, Neff NF. Concordance of MERFISH spatial transcriptomics with bulk and single-cell RNA sequencing. Life Sci Alliance. 2022 Dec 16;6(1):e202201701. doi: 10.26508/lsa.. 202201701. PMID: 36526371; PMCID: PMC9760489.

Liu, X., Zhang, Z., Tan, C. et al. Global trends in machine learning applications for single-cell transcriptomics research. Hereditas 162, 164 (2025). https://doi.org/10.1186/s41065-025-00528-y

Valkenburg KC, de Groot AE, Pienta KJ. Targeting the tumour stroma to improve cancer therapy. Nat Rev Clin Oncol. 2018 Jun;15(6):366-381. doi: 10.1038/s41571-018-0007-1. PMID: 29651130; PMCID: PMC5960434.

Lotfollahi M, Klimovskaia Susmelj A, De Donno C, Hetzel L, Ji Y, Ibarra IL, Srivatsan SR, Naghipourfar M, Daza RM, Martin B, Shendure J, McFaline-Figueroa JL, Boyeau P, Wolf FA, Yakubova N, Günnemann S, Trapnell C, Lopez-Paz D, Theis FJ. Predicting cellular responses to complex perturbations in high-throughput screens. Mol Syst Biol. 2023 Jun 12;19(6):e11517. Doi: 10.15252/msb 202211517. Epub 2023 May 8. PMID: 37154091; PMCID: PMC10258562.

Ma, C., Yang, C., Peng, A., Sun, T., Ji, X., Mi, J., Wei, L., Shen, S. and Feng, Q. (2023). Pan-cancer spatially resolved single-cell analysis reveals the crosstalk between cancer-associated fibroblasts and tumour microenvironment. *Molecular Cancer.* Doi:<https://doi.org/10.1186/s12943-023-01876-x>.

Maia, A., Schöllhorn, A., Schuhmacher, J. and Gouttefangeas, C. (2022). CAF-immune cell crosstalk and its impact in immunotherapy. *Seminars in Immunopathology.* Doi:<https://doi.org/10.1007/s00281-022-00977-x>.

Miao, C., Liu, L. and Cao, Y. (2024). OSCC‐derived EVs Educate Fibroblasts and Remodel Collagen Landscape. *Matrix Biology*, 134, pp.132–143.

Mossner, M., Baker, A.-M.C. and Graham, T.A. (2021). The role of single-cell sequencing in studying tumour evolution. *Faculty Reviews*, 10. doi:<https://doi.org/10.12703/r/10-49>.

Mumme, H., Thomas, T., Pilcher, W., Bhasin, S. et al. (2023). Single-cell analysis reveals altered tumour microenvironments of relapse- and remission-associated pediatric acute myeloid leukaemia. *Nature Communications*, 14. doi:10.1038/s41467-023-41994-0.

Niu, X. et al. (2024). Cancer plasticity in therapy resistance: Mechanisms and novel strategies. *Drug Resistance Updates*, 76, 101114. doi:10.1016/j.drup.2024.101114.

Pattabiram, S., Gangadaran, P., Dhayalan, S., Chatterjee, G., Reyaz, D., Kodali, P., Renganathan, A., Rajendran, R.L., Ahn, B. and ArulJothi, K.N. (2025). Decoding the Tumour Microenvironment: Insights and New Targets from Single-Cell Sequencing and Spatial Transcriptomics. *Current Issues in Molecular Biology.* Doi:<https://doi.org/10.3390/cimb47090730>.

Peng, L., Wang, F., Wang, Z., Tan, J., Huang, L., Tian, X., Liu, G. and Zhou, L. (2022). Cell-cell communication inference and analysis in the tumour microenvironments from single-cell transcriptomics: data resources and computational strategies. *Briefings in Bioinformatics.* Doi:<https://doi.org/10.1093/bib/bbac234>.

Peng, R., Yu, L., & Xu, B. (2025). Integrating single-cell and bulk RNA sequencing data establishes a cuproptosis-related gene predictive signature in breast cancer. Discover Oncology. <https://doi.org/10.1007/s12672-025-03525-9>

Pérez-Moreno, M. (2009). When neighbourhood matters: tumour microenvironment. *Clinical & Translational Oncology.* Doi:<https://doi.org/10.1007/S12094-009-0316-Z>.

Rafique R, Islam SMR, Kazi JU. Machine learning in the prediction of cancer therapy. Comput Struct Biotechnol J. 2021 Jul 8;19:4003-4017. Doi: 10.1016/j.csbj.2021.07.003. PMID: 34377366; PMCID: PMC8321893.

Raghavan, S. et al. (2021). Microenvironment drives cell state, plasticity, and drug response in pancreatic cancer. *Cell*, 184(25), pp.6119–6137.e26. doi:10.1016/j.cell.2021.11.017.

Rodriques SG, Stickels RR, Goeva A, Martin CA, Murray E, Vanderburg CR, Welch J, Chen LM, Chen F, Macosko EZ. Slide-seq: A scalable technology for measuring genome-wide expression at high spatial resolution. Science. 2019 Mar 29;363(6434):1463-1467. doi: 10.1126/science.aaw1219. Epub 2019 Mar 28. PMID: 30923225; PMCID: PMC6927209.

Shah, K., Basu Mallik, S., Gupta, P. and Iyer, A. (2022). Targeting Tumour-Associated Fibroblasts in Cancers. *Frontiers in Oncology.* Doi:<https://doi.org/10.3389/fonc.2022.908156>.

Simonas Juzenas, Karolis Goda, Vaidotas Kiseliovas, Justina Zvirblyte, Alvaro Quintinal-Villalonga, Juozas Siurkus, Juozas Nainys, Linas Mazutis, inDrops-2: a flexible, versatile and cost-efficient droplet microfluidic approach for high-throughput scRNA-seq of fresh and preserved clinical samples, Nucleic Acids Research, Volume 53, Issue 2, 27 January 2025, gkae1312, <https://doi.org/10.1093/nar/gkae1312>

Song, H. et al. (2022). Single-cell analysis of human primary prostate cancer reveals the heterogeneity of tumour-associated epithelial cell states. *Nature Communications*, 13, 141. doi:10.1038/s41467-021-27867-5.

Sun, X. and Yu, Q. (2015). Intra-tumour heterogeneity of cancer cells and its implications for cancer treatment. *Acta Pharmacologica Sinica*, 36(10), pp.1219–1227. doi:10.1038/aps.. 2015.92.

Sun, Y. (2015). Tumour microenvironment and cancer therapy resistance. *Cancer Letters*, 380(1), pp.205–215. doi:10.1016/j.canlet.2015.07.044.

Taavitsainen, S. et al. (2021). Single-cell ATAC and RNA sequencing reveal pre-existing and persistent cells associated with prostate cancer relapse. *Nature Communications*, 12, 5307. doi:10.1038/s41467-021-25624-1.

Teppei Shimamura, RNA velocity and beyond: Current advances in modelling single-cell transcriptional dynamics, Allergology International, Volume 74, Issue 4, 2025, Pages 525-533, ISSN 1323-8930, <https://doi.org/10.1016/j.alit.2025.08.005>.

Tieng FYF, Lee LH, Ab Mutalib NS. Single-cell RNA-sequencing of circulating tumour cells: A practical guide to workflow and translational applications. Cancer Metastasis Rev. 2025 Oct 6;44(4):75. doi: 10.1007/s10555-025-10293-z. PMID: 41053409; PMCID: PMC12500777.

Timperi, E. and Romano, E. (2023). Stromal circuits involving tumour-associated macrophages and cancer-associated fibroblasts. *Frontiers in Immunology.* Doi:<https://doi.org/10.3389/fimmu.2023.1194642>.

Vahid MR, Kurlovs AH, Andreani T, Augé F, Olfati-SDoie R, de Rinaldis E, Rapaport F, Savova V. DiSiR: fast and robust method to identify ligand-receptor interactions at the subunit level from single-cell RNA-sequencing data. NAR Genom Bioinform. 2023 Mar 23;5(1):lqad030. doi: 10.1093/nargab/lqad030. PMID: 36968431; PMCID: PMC10034587.

Wang, G., Zhao, J., Lin, Y., Liu, T., Zhao, Y. and Zhao, H. (2025). scMODAL: a general deep learning framework for comprehensive single-cell multi-omics data alignment with feature links. *Nature Communications*, 16(1). doi:<https://doi.org/10.1038/s41467-025-60333-z>.

Wang, S. et al. (2023). The evolution of Single-Cell RNA Sequencing Technology and Application: Progress and Perspectives. *International Journal of Molecular Sciences*, 24(3), 2943. doi:10.3390/ijms24032943.

Wang, T. et al. (2024). Effect of fibroblast heterogeneity on prognosis and drug response in ovarian cancer (identification of CXCL12-positive fibroblasts associated with chemoresistance). [Open access]. DOI available via PubMed Central.

Wang, W., Wang, L., Wang, L., She, J. and Zhu, J. (2021). Examining heterogeneity of stromal cells in tumour microenvironment based on pan-cancer single-cell RNA sequencing data. *Cancer Biology and Medicine.* Doi:<https://doi.org/10.20892/J.ISSN.2095-3941.2020.0762>.

Watson, I.R., Takahashi, K., Futreal, P.A. and Chin, L. (2013). Emerging patterns of somatic mutations in cancer. *Nature Reviews Genetics*, 14, pp.703–718. doi:10.1038/nrg3539.

Werba, G. et al. (2023). Single-cell RNA sequencing reveals the effects of chemotherapy on human pancreatic adenocarcinoma and its tumour microenvironment. *Nature Communications*, 14. doi:10.1038/s41467-023-36296-4.

Williams, C.G., Lee, H.J., Asatsuma, T., Vento-Tormo, R. and Haque, A. (2022). An introduction to spatial transcriptomics for biomedical research. *Genome Medicine*, 14(1). doi:<https://doi.org/10.1186/s13073-022-01075-1>.

Winkler, J. et al. (2024). Single-cell analysis of breast cancer metastasis reveals epithelial-mesenchymal plasticity signatures associated with poor outcomes. *Journal of Clinical Investigation*, 134(17), e164227. doi:10.1172/JCI164227.

Wu, F. et al. (2021). Single-cell profiling of tumour heterogeneity and the microenvironment in advanced non-small cell lung cancer. *Nature Communications*, 12, 2540. doi:10.1038/s41467-021-22801-0.

Xia, Y., Sun, T., Li, G. et al. Spatial single-cell analysis of tumour microenvironment remodelling pattern in primary central nervous system lymphoma. Leukaemia 37, 1499–1510 (2023). <https://doi.org/10.1038/s41375-023-01908-x>

Xiang Lin, Haoran Liu, Zhi Wei, Senjuti Basu Roy, Nan Gao, An active learning approach for clustering single-cell RNA-seq data, Laboratory Investigation, Volume 102, Issue 3,2022, Pages 227-235, ISSN 0023-6837,<https://doi.org/10.1038/s41374-021-00639-w>.

Xiliang Wang, Yao He, Qiming Zhang, Xianwen Ren, Zemin Zhang, Direct Comparative Analyses of 10X Genomics Chromium and Smart-Seq2, Genomics, Proteomics & Bioinformatics, Volume 19, Issue 2, April 2021, Pages 253–266, https://doi.org/10.1016/j.gpb.2020.02.005

Xu, L. et al. (2024). A comprehensive single-cell breast tumour atlas defines epithelial and immune heterogeneity and interactions predicting anti-PD-1 therapy response. *Cell Reports Medicine*, 5(5), 101511. doi:10.1016/j.xcrm.2024.101511.

Yan, C., Zhu, Y., Chen, M., Yang, K., Cui, F., Zou, Q. and Zhang, Z. (2024). Integration tools for scRNA-seq data and spatial transcriptomics sequencing data. *Briefings in Functional Genomics*, 23(4), pp.295–302. doi:<https://doi.org/10.1093/bfgp/elae002>.

Yan, H. et al. (2025). Integrated single-cell transcriptomics and Mendelian randomisation identify neutrophil-associated biomarker genes in colorectal cancer. Discover Oncology. <https://doi.org/10.1007/s12672-025-03452-9>

Yan, H. et al. (2022). Technique integration of single-cell RNA sequencing with spatially resolved transcriptomics in the tumour microenvironment. *Cancer Cell International.* Doi:<https://doi.org/10.1186/s12935-022-02580-4>.

Yang, W., Wang, P., Xu, S. et al. Deciphering cell–cell communication at single-cell resolution for spatial transcriptomics with a subgraph-based graph attention network. Nat Commun 15, 7101 (2024). https://doi.org/10.1038/s41467-024-51329-2

Yu, L., Cao, Y., Yang, J.Y.H. et al. Benchmarking clustering algorithms for estimating the number of cell types from single-cell RNA-sequencing data. Genome Biol 23, 49 (2022). https://doi.org/10.1186/s13059-022-02622-0

Yuan, Y. et al. (2025). Profiling of immune cell subsets and functional characteristics of cervical cancer based on single-cell RNA sequencing. Frontiers in Immunology. <https://doi.org/10.3389/fimmu.2025.1658705>

Zhang, L. et al. (2024). Single-cell analysis reveals stromal dynamics and tumour-immune infiltration classes in ovarian cancer. *Communications Biology.* (Open access).