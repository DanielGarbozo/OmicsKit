# TCGA-BRCA Pathway Analysis and Clustering Workflow

## Overview

This vignette demonstrates a downstream workflow for Pathway Analysis
(e.g., GSEA) using **OmicsKit**, applied to the TCGA-BRCA dataset. The
workflow covers:

1.  **Exploring enriched pathways** using single and multi-comparison
    visualizations (`splot_PA`, `multiplot_PA`).
2.  **Extracting and annotating genes** from pathway results
    (`getgenesPA`, `addgenesPA`).
3.  **Visualizing gene-level expression** inside selected pathways via
    heatmaps (`heatmap_PA`).
4.  **Reducing redundancy** by computing pathway similarity and
    clustering pathway terms (`geneset_similarity`, `do_clust`).
5.  **Detecting network communities** and identifying biological
    super-terms (`get_network_communities`, `network_clust_gg`).

Several computationally intensive visualizations are shown as
pre-rendered static figures to keep vignette compilation fast and
reproducible.

## 1. Data Loading

First, we load the package and the pre-computed TCGA-BRCA objects
provided in the `data/` folder.

``` r

library(OmicsKit)

# Load pathway analysis results, gene set lists, and ranked genes
data("brca_pa_merged")
data("brca_geneset_list")
data("brca_ranked_genes")

# Load pre-computed clustering and similarity objects
data("brca_pa_similarity")
data("brca_pa_clustering")
```

Let’s take a quick look at the merged pathway analysis results:

``` r

table(brca_pa_merged$COMPARISON, brca_pa_merged$COLLECTION)
#>                           
#>                            GO_BP Hallmark
#>   ERpositive_vs_ERnegative  3898       50
#>   Tumor_vs_Normal           3902       50
```

## 2. Visualization of Enriched Pathways

OmicsKit provides specialized functions to visualize GSEA results
easily.

### Single Comparison (`splot_PA`)

We can visualize the significant pathways for a specific comparison, for
example, *ER positive vs ER negative* within the Hallmark collection.

The code below shows how the plot was generated. The vignette displays a
pre-rendered figure to avoid re-rendering a large pathway plot during
package checks or website builds.

``` r

pa_single <- subset(
  brca_pa_merged,
  COMPARISON == "ERpositive_vs_ERnegative" & COLLECTION == "Hallmark"
)

splot_PA(
  data           = pa_single,
  geneset_col    = "NAME",
  collection_col = "COLLECTION",
  nes_col        = "NES",
  fdr_col        = "FDR",
  order          = "desc",
  fill_limits    = c(0, 5),
  top_n_per_direction = 15,
  clean_labels   = TRUE,
  label_wrap     = 32,
  theme_params   = list(
    axis_title_size   = 12,
    axis_text_size_x  = 10,
    axis_text_size_y  = 7,
    strip_text_size   = 10,
    legend_title_size = 10,
    legend_text_size  = 9,
    bar_size          = 0.3,
    bar_width         = 0.72,
    hline_size        = 0.45
  )
)
```

![Top Hallmark pathway enrichments for the ER positive vs ER negative
comparison.](figures/pathway_splot_er_status.png)

Top Hallmark pathway enrichments for the ER positive vs ER negative
comparison.

### Multiple Comparisons (`multiplot_PA`)

If we want to compare enrichment patterns across multiple conditions, we
can use
[`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md).
To keep this example compact, we focus on two comparisons and two
estrogen-response Hallmark gene sets.

``` r

pa_multi_subset <- subset(
  brca_pa_merged,
  COMPARISON %in% c("Tumor_vs_Normal", "ERpositive_vs_ERnegative") &
    NAME %in% c(
      "HALLMARK_ESTROGEN_RESPONSE_EARLY",
      "HALLMARK_ESTROGEN_RESPONSE_LATE"
    )
)

# Clean labels for display
pa_multi_subset$PATHWAY_LABEL <- pa_multi_subset$NAME
pa_multi_subset$PATHWAY_LABEL[
  pa_multi_subset$NAME == "HALLMARK_ESTROGEN_RESPONSE_EARLY"
] <- "Estrogen response\nearly"
pa_multi_subset$PATHWAY_LABEL[
  pa_multi_subset$NAME == "HALLMARK_ESTROGEN_RESPONSE_LATE"
] <- "Estrogen response\nlate"

comparison_labels <- c(
  "Tumor_vs_Normal" = "Tumor vs\nNormal",
  "ERpositive_vs_ERnegative" = "ER+ vs\nER-"
)

multiplot_PA(
  data             = pa_multi_subset,
  comparison_col   = "COMPARISON",
  facet_col        = "PATHWAY_LABEL",
  axis_y           = "NES",
  fdr_col          = "FDR",
  comparison_order = c("Tumor_vs_Normal", "ERpositive_vs_ERnegative"),
  custom_labels    = comparison_labels,
  ncol_wrap        = 2,
  free_y           = FALSE,
  fill_limits      = c(0, 5),
  theme_params     = list(
    axis_title_size     = 13,
    axis_text_size_x    = 10,
    axis_text_size_y    = 10,
    strip_text_size     = 11,
    hline_size          = 0.6,
    tick_size           = 0.4,
    tick_length         = 0.12,
    panel_spacing_multi = 0.8
  )
)
```

![Multi-comparison plot for selected Hallmark estrogen-response
pathways.](figures/pathway_multiplot_estrogen_response.png)

Multi-comparison plot for selected Hallmark estrogen-response pathways.

## 3. Gene Extraction and Annotation

Often, we need to know exactly *which* genes are associated with an
enriched pathway result. Here, we extract genes for the same ER positive
vs ER negative Hallmark subset.

``` r

pa_single <- subset(
  brca_pa_merged,
  COMPARISON == "ERpositive_vs_ERnegative" & COLLECTION == "Hallmark"
)

# Extract all pathway genes and leading-edge genes.
# The leading-edge option is lightweight and avoids requiring additional top-rank columns.
gene_lists <- plot_or_message(
  getgenesPA(
    pa_data      = pa_single,
    geneset_list = brca_geneset_list,
    ranked_genes = brca_ranked_genes$ERpositive_vs_ERnegative,
    genes        = c("all", "le")
  )
)
#> $all
#> $all$HALLMARK_E2F_TARGETS
#>   [1] "RAD50"    "WDR90"    "PAN2"     "UBR7"     "GSPT1"    "MLH1"    
#>   [7] "PDS5B"    "CDKN1B"   "RFC1"     "LUC7L3"   "BRCA1"    "CDKN1A"  
#>  [13] "PPM1D"    "RBBP7"    "PMS2"     "CTCF"     "TP53"     "IPO7"    
#>  [19] "RAD51C"   "DCTPP1"   "CBX5"     "RPA3"     "NBN"      "WEE1"    
#>  [25] "DUT"      "SHMT1"    "NAP1L1"   "TUBG1"    "TIMELESS" "STAG1"   
#>  [31] "DCK"      "BARD1"    "PA2G4"    "KIF22"    "LIG1"     "RAD1"    
#>  [37] "PRPS1"    "PNN"      "NME1"     "NUP107"   "SMC3"     "POLE4"   
#>  [43] "EIF2S1"   "POLE"     "SMC1A"    "PSMC3IP"  "RAD21"    "NAA38"   
#>  [49] "HUS1"     "RPA2"     "CSE1L"    "POLD2"    "ING3"     "MYC"     
#>  [55] "TOP2A"    "SNRPB"    "PRKDC"    "NOLC1"    "NUDT21"   "RACGAP1" 
#>  [61] "PAICS"    "MXD3"     "BRCA2"    "RPA1"     "ESPL1"    "SLBP"    
#>  [67] "ATAD2"    "HMGB2"    "UNG"      "NUP153"   "HELLS"    "TFRC"    
#>  [73] "SPAG5"    "NOP56"    "BRMS1L"   "PSIP1"    "CIT"      "POLA2"   
#>  [79] "POLD3"    "RNASEH2A" "TRA2B"    "CENPM"    "PPP1R8"   "CDK4"    
#>  [85] "ASF1B"    "XPO1"     "TIPIN"    "TK1"      "E2F8"     "SPC25"   
#>  [91] "PCNA"     "SPC24"    "EXOSC8"   "TMPO"     "CDKN2C"   "HMMR"    
#>  [97] "HNRNPD"   "LMNB1"    "DNMT1"    "PLK4"     "ILF3"     "SMC4"    
#> [103] "CDK1"     "RAN"      "POP7"     "HMGB3"    "RFC3"     "UBE2T"   
#> [109] "XRCC6"    "GINS1"    "RRM2"     "CDKN3"    "TACC3"    "BUB1B"   
#> [115] "POLD1"    "ZW10"     "AURKA"    "MCM4"     "KPNA2"    "MCM2"    
#> [121] "CENPE"    "GINS4"    "MAD2L1"   "TBRG4"    "CKS2"     "NUP205"  
#> [127] "SUV39H1"  "GINS3"    "KIF4A"    "DCLRE1B"  "EED"      "UBE2S"   
#> [133] "TCF19"    "AK2"      "ASF1A"    "KIF18B"   "MKI67"    "DSCC1"   
#> [139] "PTTG1"    "MTHFD2"   "BIRC5"    "CDC25B"   "SMC6"     "MCM3"    
#> [145] "SSRP1"    "DIAPH3"   "MYBL2"    "DLGAP5"   "USP1"     "TRIP13"  
#> [151] "ANP32E"   "MCM6"     "RAD51AP1" "MSH2"     "EZH2"     "CHEK2"   
#> [157] "PLK1"     "SYNCRIP"  "CDCA3"    "DEPDC1"   "RFC2"     "CCNB2"   
#> [163] "CKS1B"    "MELK"     "CDKN2A"   "DONSON"   "PHF5A"    "PRIM2"   
#> [169] "KIF2C"    "AURKB"    "LYAR"     "PRDX4"    "CDC25A"   "MCM7"    
#> [175] "MCM5"     "NCAPD2"   "CHEK1"    "STMN1"    "DEK"      "RANBP1"  
#> [181] "NASP"     "CDCA8"    "LBR"      "TUBB"     "CDC20"    "HMGA1"   
#> [187] "CCNE1"    "CCP110"   "CNOT9"    "CTPS1"    "DDX39A"   "H2AX"    
#> [193] "H2AZ1"    "JPT1"     "MMS22L"   "MRE11"    "ORC2"     "ORC6"    
#> [199] "SRSF1"    "SRSF2"   
#> 
#> $all$HALLMARK_ESTROGEN_RESPONSE_EARLY
#>   [1] "CA12"     "THSD4"    "MLPH"     "XBP1"     "ANXA9"    "TFF1"    
#>   [7] "GREB1"    "TFF3"     "MYB"      "SLC22A5"  "ABAT"     "CELSR1"  
#>  [13] "GFRA1"    "SLC39A6"  "MAPT"     "KCNK15"   "KDM4B"    "SYBU"    
#>  [19] "ADCY9"    "WFS1"     "PGR"      "AR"       "SLC7A2"   "BCL2"    
#>  [25] "UGCG"     "RARA"     "IL6ST"    "TTC39A"   "IGF1R"    "BHLHE40" 
#>  [31] "SCNN1A"   "MAST4"    "MED13L"   "IGFBP4"   "RAB17"    "LRIG1"   
#>  [37] "KRT18"    "MREG"     "SEMA3B"   "TPBG"     "RET"      "SLC19A2" 
#>  [43] "ARL3"     "CISH"     "ELOVL5"   "CCND1"    "STC2"     "REEP1"   
#>  [49] "TJP3"     "SIAH2"    "PDZK1"    "ABCA3"    "ELOVL2"   "SLC27A2" 
#>  [55] "ADCY1"    "NRIP1"    "AFF1"     "PEX11A"   "AMFR"     "MYOF"    
#>  [61] "ABHD2"    "DHRS2"    "ITPK1"    "TOB1"     "PRSS23"   "SULT2B1" 
#>  [67] "SLC1A4"   "ASB13"    "CELSR2"   "CBFA2T3"  "FLNB"     "HSPB8"   
#>  [73] "SEC14L2"  "KRT8"     "MUC1"     "ELF1"     "SNX24"    "GJA1"    
#>  [79] "OLFM1"    "CANT1"    "DYNLT3"   "TIPARP"   "EGR3"     "RAB31"   
#>  [85] "NPY1R"    "BLVRB"    "FKBP4"    "BAG1"     "KRT19"    "OVOL2"   
#>  [91] "FASN"     "P2RY2"    "PMAIP1"   "SLC1A1"   "MSMB"     "CALCR"   
#>  [97] "AREG"     "OLFML3"   "FOS"      "RAPGEFL1" "NADSYN1"  "UNC119"  
#> [103] "NBL1"     "CXCL12"   "PTGES"    "MPPED2"   "SYNGR1"   "RHOBTB3" 
#> [109] "NCOR2"    "DLC1"     "CHPT1"    "RBBP8"    "GLA"      "KLF4"    
#> [115] "IL17RB"   "SOX3"     "RHOD"     "KLF10"    "RPS6KA2"  "TSKU"    
#> [121] "HES1"     "MICB"     "ESRP2"    "INHBB"    "CLDN7"    "FHL2"    
#> [127] "SLC24A3"  "SYT12"    "GAB2"     "TMEM164"  "WWC1"     "JAK2"    
#> [133] "KRT13"    "FDFT1"    "MED24"    "RASGRP1"  "SH3BP5"   "TMPRSS3" 
#> [139] "INPP5F"   "SLC37A1"  "B4GALT1"  "FRK"      "ALDH3B1"  "CD44"    
#> [145] "ZNF185"   "AQP3"     "PDLIM3"   "AKAP1"    "PAPSS2"   "TIAM1"   
#> [151] "TGM2"     "FARP1"    "NAV2"     "ENDOD1"   "ELF3"     "FKBP5"   
#> [157] "ADD3"     "RRP12"    "SCARB1"   "DHRS3"    "MYBL1"    "ISG20L2" 
#> [163] "CLIC3"    "TPD52L1"  "SLC26A2"  "MYC"      "TGIF2"    "SVIL"    
#> [169] "TUBB2B"   "SFN"      "ABLIM1"   "KRT15"    "MYBBP1A"  "PPIF"    
#> [175] "CYP26B1"  "PODXL"    "TFAP2C"   "SLC2A1"   "NXT1"     "BCL11B"  
#> [181] "KLK10"    "DHCR7"    "CALB2"    "SLC16A1"  "OPN3"     "HR"      
#> [187] "LAD1"     "SLC7A5"   "KCNK5"    "FOXC1"    "CCN5"     "DEPTOR"  
#> [193] "EEIG1"    "FCMR"     "KAZN"     "MINDY1"   "NHERF1"   "PLAAT3"  
#> [199] "RETREG1"  "TBC1D30" 
#> 
#> $all$HALLMARK_G2M_CHECKPOINT
#>   [1] "TLE3"     "CCND1"    "NUMA1"    "ARID4A"   "SMAD3"    "CUL3"    
#>   [7] "PURA"     "RPS6KA5"  "GSPT1"    "MNAT1"    "SLC38A1"  "BUB3"    
#>  [13] "YTHDC1"   "SLC12A2"  "CCNT1"    "PDS5B"    "CDKN1B"   "ATRX"    
#>  [19] "TGFB1"    "EGF"      "PRMT5"    "LIG3"     "CDC27"    "ABL1"    
#>  [25] "CTCF"     "PAFAH1B1" "ODF2"     "TOP1"     "NUP98"    "BCL3"    
#>  [31] "RAD23B"   "G3BP1"    "SS18"     "CUL5"     "HMGN2"    "STAG1"   
#>  [37] "HOXC10"   "HNRNPU"   "BARD1"    "UPF1"     "KIF22"    "NOTCH2"  
#>  [43] "HSPA8"    "KPNB1"    "CUL4A"    "KIF5B"    "WRN"      "MEIS1"   
#>  [49] "GINS2"    "CBX1"     "FOXN3"    "POLE"     "SMC1A"    "MEIS2"   
#>  [55] "SQLE"     "SMARCC1"  "RAD21"    "MAPK14"   "HUS1"     "RBM14"   
#>  [61] "RPA2"     "TNPO2"    "CCNF"     "MT2A"     "MYC"      "SMC2"    
#>  [67] "DTYMK"    "TOP2A"    "EWSR1"    "KIF20B"   "NOLC1"    "DR1"     
#>  [73] "RACGAP1"  "TRAIP"    "BRCA2"    "ESPL1"    "RBL1"     "PBK"     
#>  [79] "CASP8AP2" "NEK2"     "ATF5"     "POLA2"    "RASAL2"   "MARCKS"  
#>  [85] "CDC6"     "TRA2B"    "CHMP1A"   "CDK4"     "XPO1"     "NUSAP1"  
#>  [91] "CHAF1A"   "NUP50"    "TMPO"     "CDKN2C"   "HMMR"     "HNRNPD"  
#>  [97] "LMNB1"    "PTTG3P"   "HIF1A"    "SFPQ"     "PLK4"     "ILF3"    
#> [103] "SMC4"     "SAP30"    "CDK1"     "POLQ"     "CUL1"     "FANCC"   
#> [109] "SLC7A1"   "E2F1"     "MTF2"     "HMGB3"    "PRC1"     "HIRA"    
#> [115] "TROAP"    "KIF11"    "CDKN3"    "TACC3"    "NCL"      "FBXO5"   
#> [121] "KIF15"    "AURKA"    "KPNA2"    "MCM2"     "CENPE"    "MAD2L1"  
#> [127] "CKS2"     "SUV39H1"  "KIF4A"    "KIF23"    "EFNA5"    "UBE2S"   
#> [133] "CENPF"    "MKI67"    "PTTG1"    "DMD"      "UBE2C"    "PML"     
#> [139] "CDC7"     "BIRC5"    "SNRPD1"   "CDC25B"   "INCENP"   "EXO1"    
#> [145] "CDC45"    "MCM3"     "MYBL2"    "E2F2"     "MCM6"     "TFDP1"   
#> [151] "TPX2"     "E2F4"     "EZH2"     "PLK1"     "BUB1"     "CCNA2"   
#> [157] "SYNCRIP"  "RAD54L"   "CCNB2"    "CKS1B"    "KATNA1"   "NDC80"   
#> [163] "DBF4"     "STIL"     "PRIM2"    "KIF2C"    "AURKB"    "CDC25A"  
#> [169] "MCM5"     "CHEK1"    "SLC7A5"   "STMN1"    "TTK"      "ODC1"    
#> [175] "NASP"     "LBR"      "AMD1"     "UCK2"     "CENPA"    "DKC1"    
#> [181] "CDC20"    "HMGA1"    "E2F3"     "DDX39A"   "H2AX"     "H2AZ1"   
#> [187] "H2AZ2"    "H2BC12"   "JPT1"     "KMT5A"    "KNL1"     "MAP3K20" 
#> [193] "NSD2"     "ORC5"     "ORC6"     "PRP4K"    "SRSF1"    "SRSF10"  
#> [199] "SRSF2"    "TENT4A"  
#> 
#> $all$HALLMARK_MYC_TARGETS_V1
#>   [1] "FAM120A"   "PRDX3"     "GSPT1"     "IMPDH2"    "BUB3"      "PTGES3"   
#>   [7] "GOT2"      "CANX"      "PSMC6"     "DHX15"     "EIF4E"     "ETF1"     
#>  [13] "CCT2"      "RAD23B"    "AP3S1"     "RSL1D1"    "HNRNPC"    "EIF4G2"   
#>  [19] "EIF1AX"    "NHP2"      "G3BP1"     "RPS5"      "TUFM"      "VDAC1"    
#>  [25] "DUT"       "RPL14"     "MRPL23"    "EIF3J"     "GLO1"      "TRIM28"   
#>  [31] "HNRNPA1"   "NAP1L1"    "RRM1"      "RPL22"     "PSMA6"     "MRPS18B"  
#>  [37] "APEX1"     "SF3A1"     "RPL18"     "RUVBL2"    "HNRNPU"    "PA2G4"    
#>  [43] "ERH"       "NCBP2"     "RPL34"     "KPNB1"     "NDUFAB1"   "CSTF2"    
#>  [49] "COPS5"     "NPM1"      "PSMD3"     "PSMB3"     "HNRNPR"    "RPLP0"    
#>  [55] "NME1"      "PSMA1"     "CLNS1A"    "RPL6"      "RPS2"      "PRPF31"   
#>  [61] "EIF2S1"    "CYC1"      "YWHAE"     "EXOSC7"    "GNL3"      "SLC25A3"  
#>  [67] "UBE2E1"    "CBX3"      "VDAC3"     "LDHA"      "SMARCC1"   "HNRNPA2B1"
#>  [73] "RPS3"      "PSMD8"     "RNPS1"     "C1QBP"     "TXNL4A"    "NCBP1"    
#>  [79] "POLE3"     "POLD2"     "PWP1"      "CNBP"      "HSP90AB1"  "EIF4A1"   
#>  [85] "EIF4H"     "CDK2"      "RRP9"      "MYC"       "PCBP1"     "PSMD1"    
#>  [91] "SNRPB2"    "SNRPD2"    "RPS6"      "NOLC1"     "ACP1"      "DDX21"    
#>  [97] "UBE2L3"    "SET"       "PSMA2"     "PSMC4"     "EEF1B2"    "VBP1"     
#> [103] "PPIA"      "AIMP2"     "HNRNPA3"   "HPRT1"     "NOP56"     "TARDBP"   
#> [109] "HSPE1"     "EIF3D"     "TRA2B"     "PRPS2"     "CCT3"      "PHB2"     
#> [115] "CDK4"      "XPO1"      "RPS10"     "SNRPD3"    "XPOT"      "PCNA"     
#> [121] "HSPD1"     "HNRNPD"    "CCT7"      "PSMA4"     "SSB"       "LSM7"     
#> [127] "COX5A"     "PGK1"      "ABCE1"     "NOP16"     "PSMA7"     "RAN"      
#> [133] "CUL1"      "SSBP1"     "HDDC2"     "EIF2S2"    "XRCC6"     "FBL"      
#> [139] "PSMD14"    "STARD7"    "PABPC1"    "MCM4"      "KPNA2"     "MCM2"     
#> [145] "MAD2L1"    "CAD"       "SRM"       "EIF3B"     "DDX18"     "UBA2"     
#> [151] "TCP1"      "SNRPA"     "TYMS"      "CCT5"      "U2AF1"     "MRPL9"    
#> [157] "SNRPD1"    "PABPC4"    "PSMD7"     "SF3B3"     "SNRPA1"    "CDC45"    
#> [163] "YWHAQ"     "USP1"      "PPM1G"     "MCM6"      "TFDP1"     "CCNA2"    
#> [169] "SYNCRIP"   "SNRPG"     "CCT4"      "HDGF"      "RFC4"      "PRDX4"    
#> [175] "MCM7"      "MCM5"      "HDAC2"     "LSM2"      "DEK"       "RANBP1"   
#> [181] "ODC1"      "ILF2"      "SERBP1"    "PSMB2"     "CDC20"     "SRPK1"    
#> [187] "IFRD1"     "CTPS1"     "EPRS1"     "H2AZ1"     "IARS1"     "KARS1"    
#> [193] "ORC2"      "PHB1"      "RACK1"     "SRSF1"     "SRSF2"     "SRSF3"    
#> [199] "SRSF7"     "TOMM70"   
#> 
#> $all$HALLMARK_ALLOGRAFT_REJECTION
#>   [1] "IKBKB"    "EIF4G3"   "APBB1"    "F2R"      "TLR3"     "TPD52"   
#>   [7] "TGFB1"    "ELANE"    "CARTPT"   "IRF7"     "BRCA1"    "KRT1"    
#>  [13] "LY86"     "THY1"     "AKT1"     "INHBB"    "INHBA"    "BCL3"    
#>  [19] "IL16"     "ITGAL"    "STAB1"    "RPS3A"    "RPL9"     "TIMP1"   
#>  [25] "NLRP3"    "CSF1"     "PTPN6"    "EIF3J"    "UBE2N"    "DYRK3"   
#>  [31] "JAK2"     "MBL2"     "CCND2"    "GPR65"    "EIF3A"    "IGSF6"   
#>  [37] "IL4R"     "GALNT1"   "CCL11"    "RPS9"     "PRKCG"    "TRAF2"   
#>  [43] "NPM1"     "HLA-DOA"  "CCL22"    "IL2"      "CD74"     "IL11"    
#>  [49] "BCL10"    "NME1"     "GCNT1"    "HLA-DMB"  "CRTAM"    "CCL19"   
#>  [55] "ST8SIA4"  "ACVR2A"   "CFP"      "C2"       "DEGS1"    "CD4"     
#>  [61] "FCGR2B"   "IRF8"     "HLA-DRA"  "B2M"      "CD28"     "CD1D"    
#>  [67] "HLA-DQA1" "HLA-DMA"  "GBP2"     "PTPRC"    "SPI1"     "PSMB10"  
#>  [73] "EREG"     "PRKCB"    "IL10"     "LCP2"     "CD8A"     "HCLS1"   
#>  [79] "CAPG"     "TRAT1"    "STAT1"    "CD40LG"   "IL12B"    "CCR5"    
#>  [85] "IL1B"     "CD86"     "MMP9"     "CD3G"     "TGFB2"    "IL18"    
#>  [91] "STAT4"    "WAS"      "FGR"      "CXCR3"    "CCR2"     "IL2RB"   
#>  [97] "NOS2"     "ITK"      "SOCS5"    "HLA-G"    "NCF4"     "CD3E"    
#> [103] "GZMA"     "CD80"     "TAPBP"    "SRGN"     "ITGB2"    "CD2"     
#> [109] "CD3D"     "CD96"     "FASLG"    "IL15"     "CCND3"    "ZAP70"   
#> [115] "IL6"      "SIT1"     "CD8B"     "IL13"     "LTB"      "FAS"     
#> [121] "HLA-E"    "IL12RB1"  "MAP4K1"   "CCR1"     "CD247"    "ACHE"    
#> [127] "IL7"      "IFNGR2"   "CTSS"     "CXCL9"    "BCAT1"    "CXCL13"  
#> [133] "NCR1"     "LY75"     "LIF"      "ETS1"     "RPL39"    "CD40"    
#> [139] "MAP3K7"   "TLR1"     "ELF4"     "HLA-A"    "CCL4"     "IL2RG"   
#> [145] "GLMN"     "KLRD1"    "EIF3D"    "HDAC9"    "CCL2"     "TLR6"    
#> [151] "CD79A"    "CD47"     "LCK"      "ICOSLG"   "FLNA"     "TNF"     
#> [157] "TLR2"     "IL18RAP"  "CCL5"     "CD7"      "MRPL3"    "PRF1"    
#> [163] "ABCE1"    "TAP1"     "IRF4"     "HIF1A"    "IFNG"     "IL12A"   
#> [169] "CCL13"    "IL27RA"   "RPS19"    "ICAM1"    "UBE2D1"   "SOCS1"   
#> [175] "EIF5A"    "IFNGR1"   "MTIF2"    "IL2RA"    "ABI1"     "GZMB"    
#> [181] "TAP2"     "CCL7"     "CSK"      "NCK1"     "HLA-DOB"  "LYN"     
#> [187] "CDKN2A"   "IFNAR2"   "RIPK2"    "EGFR"     "AARS1"    "DARS1"   
#> [193] "F2"       "FYB1"     "IL4"      "IL9"      "PF4"      "RARS1"   
#> [199] "RPL3L"    "WARS1"   
#> 
#> $all$HALLMARK_MYC_TARGETS_V2
#>  [1] "SORD"      "RABEPK"    "IPO4"      "DCTPP1"    "HK2"       "PRMT3"    
#>  [7] "GRWD1"     "TMEM97"    "NOC4L"     "PA2G4"     "NPM1"      "SLC29A2"  
#> [13] "TFB2M"     "FARSA"     "GNL3"      "RRP12"     "LAS1L"     "CBX3"     
#> [19] "UTP20"     "EXOSC5"    "RRP9"      "MYC"       "NOLC1"     "IMP4"     
#> [25] "MYBBP1A"   "AIMP2"     "SLC19A1"   "UNG"       "RCL1"      "NOP56"    
#> [31] "HSPE1"     "MAP3K6"    "WDR74"     "CDK4"      "HSPD1"     "PUS1"     
#> [37] "PPRC1"     "PPAN"      "NOP16"     "DUSP2"     "PLK4"      "MRTO4"    
#> [43] "PES1"      "MCM4"      "TBRG4"     "SRM"       "SUPV3L1"   "DDX18"    
#> [49] "TCOF1"     "MPHOSPH10" "PLK1"      "NOP2"      "BYSL"      "WDR43"    
#> [55] "NIP7"      "MCM5"      "NDUFAF4"   "PHB1"     
#> 
#> $all$HALLMARK_ESTROGEN_RESPONSE_LATE
#>   [1] "CA12"       "XBP1"       "ANXA9"      "TFF1"       "AGR2"      
#>   [6] "SCUBE2"     "TFF3"       "MYB"        "SLC22A5"    "RABEP1"    
#>  [11] "CACNA2D2"   "MAPT"       "DNAJC12"    "WFS1"       "PGR"       
#>  [16] "BCL2"       "IL6ST"      "TSPAN13"    "SCNN1A"     "IGFBP4"    
#>  [21] "SEMA3B"     "TPBG"       "RET"        "ARL3"       "CISH"      
#>  [26] "ELOVL5"     "CCND1"      "SERPINA5"   "ACOX2"      "TJP3"      
#>  [31] "SIAH2"      "PDZK1"      "PTGER3"     "ABCA3"      "PRLR"      
#>  [36] "SLC27A2"    "NRIP1"      "UGDH"       "AFF1"       "TPSAB1"    
#>  [41] "AMFR"       "COX6C"      "MYOF"       "ABHD2"      "DHRS2"     
#>  [46] "DLG5"       "ITPK1"      "TOB1"       "SORD"       "PRSS23"    
#>  [51] "SULT2B1"    "SLC1A4"     "HMGCS2"     "CELSR2"     "CHST8"     
#>  [56] "FGFR3"      "CYP4F11"    "FLNB"       "EMP2"       "HSPB8"     
#>  [61] "ST6GALNAC2" "MOCS2"      "OLFM1"      "CXCL14"     "DNAJC1"    
#>  [66] "ALDH3A2"    "CPE"        "DYNLT3"     "EGR3"       "RAB31"     
#>  [71] "NPY1R"      "SERPINA3"   "BLVRB"      "FKBP4"      "BAG1"      
#>  [76] "PDCD4"      "PLAC1"      "KRT19"      "OVOL2"      "ASCL1"     
#>  [81] "DCXR"       "BATF"       "PLXNB1"     "CALCR"      "AREG"      
#>  [86] "FOS"        "RAPGEFL1"   "NBL1"       "CXCL12"     "SLC2A8"    
#>  [91] "PTGES"      "PRKAR2B"    "NCOR2"      "UNC13B"     "CHPT1"     
#>  [96] "METTL3"     "RBBP8"      "GLA"        "KLF4"       "IGSF1"     
#> [101] "IL17RB"     "SOX3"       "SERPINA1"   "JAK1"       "CA2"       
#> [106] "RPS6KA2"    "CAV1"       "ZFP36"      "SNX10"      "MICB"      
#> [111] "CD9"        "LLGL2"      "TFPI2"      "NAB2"       "GALE"      
#> [116] "TNNC1"      "SLC24A3"    "CDH1"       "SLC29A1"    "MEST"      
#> [121] "TST"        "PTPN6"      "ATP2B4"     "JAK2"       "HOMER2"    
#> [126] "KRT13"      "CKB"        "ETFB"       "FDFT1"      "TMPRSS3"   
#> [131] "FRK"        "ALDH3B1"    "CD44"       "MAPK13"     "PDLIM3"    
#> [136] "KLK11"      "TH"         "PAPSS2"     "TIAM1"      "PKP3"      
#> [141] "GINS2"      "LTF"        "FARP1"      "FKBP5"      "ADD3"      
#> [146] "HSPA4L"     "SCARB1"     "CLIC3"      "TPD52L1"    "XRCC3"     
#> [151] "ID2"        "SLC26A2"    "PCP4"       "SGK1"       "MDK"       
#> [156] "TOP2A"      "SFN"        "HPRT1"      "PPIF"       "CYP26B1"   
#> [161] "LAMC2"      "ISG20"      "RNASEH2A"   "CDC6"       "TFAP2C"    
#> [166] "CCNA1"      "IDH2"       "DUSP2"      "PLK4"       "NXT1"      
#> [171] "KLK10"      "DHCR7"      "FABP5"      "SLC16A1"    "OPN3"      
#> [176] "KIF20A"     "ASS1"       "PERP"       "TRIM29"     "NMU"       
#> [181] "HR"         "S100A9"     "LSR"        "STIL"       "ST14"      
#> [186] "IMPA2"      "SLC7A5"     "GJB3"       "KCNK5"      "BTG3"      
#> [191] "GAL"        "CDC20"      "FOXC1"      "CCN5"       "EEIG1"     
#> [196] "GFUS"       "GPER1"      "LARGE1"     "NHERF1"     "PLAAT3"    
#> 
#> $all$HALLMARK_MTORC1_SIGNALING
#>   [1] "XBP1"     "BHLHE40"  "SYTL2"    "BTG2"     "ELOVL5"   "SORD"    
#>   [7] "QDPR"     "SLC1A4"   "TM7SF2"   "DHCR24"   "FDXR"     "TBK1"    
#>  [13] "IGFBP5"   "EDEM1"    "ADIPOR2"  "CCNG1"    "ALDOA"    "HSPA9"   
#>  [19] "HMGCR"    "LTA4H"    "ACACA"    "ATP6V1D"  "SQSTM1"   "HSPA4"   
#>  [25] "STC1"     "NMT1"     "CDKN1A"   "SCD"      "TRIB3"    "CANX"    
#>  [31] "GSR"      "ACLY"     "UFM1"     "PSMC6"    "GLA"      "SLC6A6"  
#>  [37] "CYP51A1"  "GGA2"     "ACSL3"    "UBE2D3"   "RDH11"    "CD9"     
#>  [43] "ELOVL6"   "PIK3R3"   "GCLC"     "ETF1"     "CTH"      "HK2"     
#>  [49] "INSIG1"   "TCEA1"    "SERP1"    "DHFR"     "FGL2"     "GSK3B"   
#>  [55] "GTF2H1"   "FKBP2"    "NUPR1"    "USO1"     "PSME3"    "GBE1"    
#>  [61] "TMEM97"   "TUBG1"    "LDLR"     "EGLN3"    "COPS5"    "RIT1"    
#>  [67] "ATP2A2"   "GLRX"     "SEC11A"   "P4HA1"    "SKAP2"    "M6PR"    
#>  [73] "DDIT3"    "CFP"      "MTHFD2L"  "SDF2L1"   "SQLE"     "SLA"     
#>  [79] "ADD3"     "LDHA"     "HSPA5"    "SLC37A4"  "IDH1"     "HMGCS1"  
#>  [85] "FADS1"    "SLC1A5"   "PSMA3"    "PSMD13"   "CCNF"     "PPP1R15A"
#>  [91] "CYB5B"    "SLC2A3"   "PSMD12"   "RRP9"     "SLC7A11"  "GOT1"    
#>  [97] "ITGB2"    "NUFIP1"   "PSMB5"    "EEF1E1"   "IMMT"     "G6PD"    
#> [103] "FADS2"    "CORO1A"   "RPA1"     "PSMC4"    "PRDX1"    "PPIA"    
#> [109] "TXNRD1"   "BCAT1"    "UNG"      "IDI1"     "PLOD2"    "TFRC"    
#> [115] "IFI30"    "HPRT1"    "STIP1"    "PFKL"     "HSPE1"    "SSR1"    
#> [121] "RAB1A"    "CXCR4"    "YKT6"     "NFKBIB"   "CACYBP"   "LGMN"    
#> [127] "MLLT11"   "HSP90B1"  "HSPD1"    "NFYC"     "SLC2A1"   "DAPP1"   
#> [133] "PSMA4"    "STARD4"   "ME1"      "ABCF2"    "PGK1"     "EBP"     
#> [139] "PPA1"     "PSPH"     "HMBS"     "EIF2S2"   "ACTR2"    "DDIT4"   
#> [145] "ARPC5L"   "RRM2"     "PSMD14"   "PSMC2"    "DHCR7"    "AURKA"   
#> [151] "MCM4"     "SERPINH1" "GPI"      "MCM2"     "NAMPT"    "NUP205"  
#> [157] "RPN1"     "MAP2K3"   "PDAP1"    "PITPNB"   "SHMT2"    "UCHL5"   
#> [163] "TOMM40"   "CCT6A"    "PNO1"     "MTHFD2"   "GAPDH"    "TUBA4A"  
#> [169] "TPI1"     "CALR"     "PNP"      "PGM1"     "CTSC"     "POLR3G"  
#> [175] "VLDLR"    "TES"      "PLK1"     "BUB1"     "NFIL3"    "ASNS"    
#> [181] "CDC25A"   "ACTR3"    "PDK1"     "PSMG1"    "SLC7A5"   "GMPS"    
#> [187] "PHGDH"    "ENO1"     "SRD5A1"   "IFRD1"    "PSAT1"    "AK4"     
#> [193] "ATP5MC1"  "DDX39A"   "EPRS1"    "ERO1A"    "NHERF1"   "NIBAN1"  
#> [199] "SC5D"     "WARS1"   
#> 
#> $all$HALLMARK_INTERFERON_GAMMA_RESPONSE
#>   [1] "IFITM2"   "EIF4E3"   "RAPGEF6"  "CFB"      "PSME1"    "MVP"     
#>   [7] "PTPN1"    "DHX58"    "ISOC1"    "CASP7"    "IFI35"    "TRAFD1"  
#>  [13] "TOR1B"    "RBCK1"    "TDRD7"    "ZNFX1"    "TXNIP"    "BST2"    
#>  [19] "SP110"    "IRF7"     "RNF31"    "IRF9"     "RNF213"   "CDKN1A"  
#>  [25] "SPPL2A"   "DDX60"    "IRF2"     "NCOA3"    "TNFSF10"  "SELP"    
#>  [31] "LGALS3BP" "RIPK1"    "ST3GAL5"  "TNFAIP6"  "IFIT1"    "STAT3"   
#>  [37] "SRI"      "UBE2L6"   "OGFR"     "CFH"      "OAS3"     "FGL2"    
#>  [43] "PTPN6"    "SSPN"     "PARP14"   "NFKB1"    "ADAR"     "AUTS2"   
#>  [49] "CASP8"    "JAK2"     "CMKLR1"   "IFITM3"   "IL4R"     "NFKBIA"  
#>  [55] "TRIM14"   "STAT2"    "TRIM21"   "ISG15"    "IFI27"    "PSME2"   
#>  [61] "OAS2"     "CD69"     "ARID5B"   "SAMD9L"   "VAMP8"    "IFIT3"   
#>  [67] "SAMHD1"   "SLC25A28" "IFIT2"    "LAP3"     "LATS2"    "BTG1"    
#>  [73] "CD74"     "TRIM25"   "IRF5"     "P2RY14"   "FPR1"     "ST8SIA4" 
#>  [79] "HLA-DRB1" "PSMB8"    "XAF1"     "VAMP5"    "HERC6"    "MX1"     
#>  [85] "IRF8"     "B2M"      "IRF1"     "OASL"     "BPGM"     "SECTM1"  
#>  [91] "HLA-DQA1" "SERPING1" "HLA-DMA"  "IFI44L"   "RSAD2"    "FCGR1A"  
#>  [97] "IL10RA"   "PSMB10"   "KLRK1"    "GBP4"     "LCP2"     "PARP12"  
#> [103] "ARL4A"    "BANK1"    "STAT1"    "PSMA3"    "CD86"     "SOCS3"   
#> [109] "CMPK2"    "MT2A"     "STAT4"    "GCH1"     "IL2RB"    "TNFAIP2" 
#> [115] "CIITA"    "GBP6"     "HLA-G"    "GZMA"     "C1S"      "APOL6"   
#> [121] "TAPBP"    "RTP4"     "CD274"    "LYSMD2"   "CASP3"    "XCL1"    
#> [127] "EPSTI1"   "USP18"    "GPR18"    "IL15"     "BATF2"    "IL6"     
#> [133] "CASP4"    "CASP1"    "C1R"      "PSMA2"    "CSF2RB"   "FAS"     
#> [139] "EIF2AK2"  "PSMB9"    "NOD1"     "ITGB7"    "VCAM1"    "IL7"     
#> [145] "HLA-B"    "CXCL9"    "ZBP1"     "IFI44"    "IFI30"    "IFIH1"   
#> [151] "CD40"     "LY6E"     "ISG20"    "SLAMF7"   "HLA-A"    "PDE4B"   
#> [157] "MYD88"    "CCL2"     "IL18BP"   "MX2"      "NLRC5"    "TNFAIP3" 
#> [163] "CXCL11"   "CCL5"     "PLSCR1"   "TRIM26"   "TAP1"     "IRF4"    
#> [169] "HIF1A"    "PTGS2"    "IL15RA"   "ICAM1"    "PNPT1"    "CXCL10"  
#> [175] "SOCS1"    "CD38"     "NAMPT"    "CCL7"     "IDO1"     "NMI"     
#> [181] "PML"      "MTHFD2"   "PNP"      "PLA2G4A"  "PIM1"     "SOD2"    
#> [187] "PTPN2"    "IFNAR2"   "PFKP"     "RIPK2"    "UPP1"     "PELI1"   
#> [193] "NUP93"    "PSMB2"    "CMTR1"    "HELZ2"    "MARCHF1"  "RIGI"    
#> [199] "TMT1B"    "WARS1"   
#> 
#> $all$HALLMARK_INFLAMMATORY_RESPONSE
#>   [1] "SLC7A2"   "BTG2"     "TPBG"     "ACVR1B"   "P2RX4"    "HPN"     
#>   [7] "SLC1A2"   "SCN1B"    "FFAR2"    "PSEN1"    "P2RY2"    "NPFFR2"  
#>  [13] "LPAR1"    "TLR3"     "GPC3"     "CSF3R"    "BST2"     "AXL"     
#>  [19] "IRF7"     "TACR3"    "NDP"      "PCDH7"    "CALCRL"   "ATP2B1"  
#>  [25] "CDKN1A"   "PTPRE"    "PTGER2"   "TNFSF10"  "IFITM1"   "SERPINE1"
#>  [31] "ABCA1"    "BDKRB1"   "VIP"      "MSR1"     "TNFAIP6"  "ITGB3"   
#>  [37] "SRI"      "APLNR"    "INHBA"    "AHR"      "OLR1"     "IFNAR1"  
#>  [43] "STAB1"    "C3AR1"    "SLC4A4"   "TIMP1"    "P2RX7"    "GABBR1"  
#>  [49] "NLRP3"    "CSF1"     "F3"       "NFKB1"    "BEST1"    "SLC11A2" 
#>  [55] "PDPN"     "SCARF1"   "RGS16"    "CMKLR1"   "IL1R1"    "SLC31A2" 
#>  [61] "RGS1"     "IL4R"     "NFKBIA"   "LDLR"     "RASGRP1"  "ITGA5"   
#>  [67] "CSF3"     "CCRL2"    "PIK3R5"   "CD69"     "OSMR"     "MXD1"    
#>  [73] "HAS2"     "C5AR1"    "MMP14"    "RNF144B"  "ATP2A2"   "KCNMB2"  
#>  [79] "EMP3"     "CCL22"    "CLEC5A"   "GPR132"   "PTGIR"    "CD55"    
#>  [85] "PTAFR"    "EBI3"     "FPR1"     "ACVR2A"   "GP1BA"    "KLF6"    
#>  [91] "SGMS2"    "IRF1"     "NOD2"     "CYBB"     "SELL"     "SELE"    
#>  [97] "IL10RA"   "EREG"     "RELA"     "PTGER4"   "IL10"     "MEP1A"   
#> [103] "LCP2"     "HRH1"     "CD48"     "CCR7"     "IL7R"     "IL12B"   
#> [109] "IL1B"     "KCNA3"    "MEFV"     "IL18"     "GCH1"     "IL2RB"   
#> [115] "RAF1"     "GNA15"    "GPR183"   "MYC"      "HBEGF"    "IRAK2"   
#> [121] "TAPBP"    "RTP4"     "FZD5"     "TACR1"    "TNFRSF1B" "IL15"    
#> [127] "IL6"      "OSM"      "TNFSF15"  "NMUR1"    "CCL17"    "EIF2AK2" 
#> [133] "RHOG"     "IFNGR2"   "ADRM1"    "CXCL9"    "SLAMF1"   "CD14"    
#> [139] "CD70"     "LIF"      "TNFSF9"   "KCNJ2"    "CXCR6"    "CD40"    
#> [145] "SLC31A1"  "LY6E"     "DCBLD2"   "TLR1"     "PDE4B"    "CCL2"    
#> [151] "GNAI3"    "IL18R1"   "LCK"      "TNFRSF9"  "ICOSLG"   "LTA"     
#> [157] "TLR2"     "IL18RAP"  "CXCL11"   "CCL5"     "EDN1"     "SEMA4D"  
#> [163] "PROK2"    "HIF1A"    "IL1A"     "ICAM4"    "KIF1B"    "AQP9"    
#> [169] "ITGB8"    "SLC7A1"   "IL15RA"   "ICAM1"    "SPHK1"    "CXCL10"  
#> [175] "ADM"      "NAMPT"    "ATP2C1"   "ABI1"     "PLAUR"    "CXCL6"   
#> [181] "CCL7"     "ROS1"     "NMI"      "CD82"     "ADORA2B"  "MET"     
#> [187] "CX3CL1"   "LYN"      "LAMP3"    "CCL20"    "MARCO"    "RIPK2"   
#> [193] "CHST2"    "PVR"      "OPRK1"    "ADGRE1"   "CCL24"    "CXCL8"   
#> [199] "SELENOS"  "SLC28A2" 
#> 
#> $all$HALLMARK_IL6_JAK_STAT3_SIGNALING
#>  [1] "IL6ST"     "ACVR1B"    "STAM2"     "PTPN1"     "IL13RA1"   "GRB2"     
#>  [7] "CSF3R"     "TGFB1"     "IRF9"      "CD36"      "LEPR"      "PDGFC"    
#> [13] "IL17RB"    "STAT3"     "ITGB3"     "CD9"       "IL9R"      "IFNAR1"   
#> [19] "CSF1"      "ACVRL1"    "JUN"       "IL1R1"     "ITGA4"     "IL4R"     
#> [25] "STAT2"     "PIK3R5"    "IL3RA"     "TYK2"      "OSMR"      "HMOX1"    
#> [31] "PTPN11"    "CD44"      "A2M"       "EBI3"      "CNTFR"     "PLA2G2A"  
#> [37] "TNFRSF1A"  "IRF1"      "MAP3K8"    "STAT1"     "IL1B"      "HAX1"     
#> [43] "SOCS3"     "CSF2RA"    "IL17RA"    "INHBE"     "CBL"       "TNFRSF1B" 
#> [49] "IL6"       "CSF2RB"    "LTB"       "FAS"       "IL12RB1"   "TNFRSF12A"
#> [55] "CCR1"      "IL7"       "IFNGR2"    "CXCL9"     "CXCL13"    "CD14"     
#> [61] "IL2RG"     "MYD88"     "IL18R1"    "TNF"       "TLR2"      "CXCL11"   
#> [67] "IL10RB"    "LTBR"      "IL15RA"    "CXCL10"    "SOCS1"     "CSF2"     
#> [73] "IFNGR1"    "CD38"      "IL2RA"     "CCL7"      "CXCL3"     "BAK1"     
#> [79] "IL1R2"     "PIM1"      "PTPN2"     "CXCL1"     "TNFRSF21"  "CRLF2"    
#> [85] "DNTT"      "PF4"       "REG1A"    
#> 
#> $all$HALLMARK_TNFA_SIGNALING_VIA_NFKB
#>   [1] "SLC16A6"  "RHOB"     "IL6ST"    "BHLHE40"  "BTG2"     "CCND1"   
#>   [7] "NINJ1"    "DUSP5"    "SMAD3"    "TIPARP"   "EGR3"     "NR4A2"   
#>  [13] "DUSP4"    "KLF2"     "IER3"     "SDC4"     "DUSP1"    "AREG"    
#>  [19] "KLF9"     "FOS"      "REL"      "PLK2"     "FOSB"     "SQSTM1"  
#>  [25] "DENND5A"  "PFKFB3"   "PER1"     "ATP2B1"   "CDKN1A"   "PTPRE"   
#>  [31] "GEM"      "TSC22D1"  "KLF4"     "EIF1"     "SERPINE1" "ABCA1"   
#>  [37] "KDM6B"    "KLF10"    "EGR1"     "JAG1"     "TGIF1"    "ZFP36"   
#>  [43] "GADD45B"  "PDLIM5"   "HES1"     "TNFAIP6"  "IRS2"     "JUNB"    
#>  [49] "BCL6"     "INHBA"    "OLR1"     "BCL3"     "CSF1"     "JUN"     
#>  [55] "F3"       "FOSL2"    "NFKB1"    "CLCF1"    "TRIB1"    "TNIP2"   
#>  [61] "NFE2L2"   "NR4A1"    "NFKBIA"   "EGR2"     "CEBPD"    "F2RL1"   
#>  [67] "LDLR"     "CCRL2"    "TNIP1"    "CD69"     "B4GALT1"  "MXD1"    
#>  [73] "NFAT5"    "ZC3H12A"  "CD44"     "SIK1"     "ZBTB10"   "IFIT2"   
#>  [79] "GADD45A"  "IER2"     "G0S2"     "BTG1"     "PMEPA1"   "STAT5A"  
#>  [85] "PLAU"     "MSC"      "EFNA1"    "KLF6"     "IRF1"     "FJX1"    
#>  [91] "TNC"      "GFPT2"    "LITAF"    "TUBB2A"   "DNAJB4"   "DRAM1"   
#>  [97] "NR4A3"    "RELA"     "PTGER4"   "IER5"     "MAP3K8"   "IL7R"    
#> [103] "IL12B"    "IL1B"     "PLEK"     "CCNL1"    "ATF3"     "SNN"     
#> [109] "SOCS3"    "MCL1"     "PPP1R15A" "IL18"     "GCH1"     "PHLDA1"  
#> [115] "TNFAIP2"  "ID2"      "SLC2A3"   "GPR183"   "MYC"      "HBEGF"   
#> [121] "SGK1"     "CD80"     "PHLDA2"   "TNFAIP8"  "IL6"      "LAMB3"   
#> [127] "IFNGR2"   "SERPINB2" "LIF"      "TNFSF9"   "BIRC3"    "IFIH1"   
#> [133] "CCL4"     "PDE4B"    "BIRC2"    "CXCL2"    "TRAF1"    "CCL2"    
#> [139] "MARCKS"   "SERPINB8" "NFKB2"    "TNFRSF9"  "ICOSLG"   "RNF19B"  
#> [145] "TNF"      "MAFF"     "TNFAIP3"  "TLR2"     "CXCL11"   "ETS2"    
#> [151] "CCL5"     "EHD1"     "PANX1"    "EDN1"     "SAT1"     "KYNU"    
#> [157] "SPSB1"    "TAP1"     "IL1A"     "DUSP2"    "CD83"     "SLC2A6"  
#> [163] "IL23A"    "TRIP10"   "TANK"     "BMP2"     "YRDC"     "B4GALT5" 
#> [169] "PTGS2"    "RELB"     "IL15RA"   "VEGFA"    "ICAM1"    "SPHK1"   
#> [175] "CXCL10"   "CSF2"     "NAMPT"    "MAP2K3"   "PLAUR"    "FUT4"    
#> [181] "CXCL6"    "CFLAR"    "CXCL3"    "RCAN1"    "NFKBIE"   "BCL2A1"  
#> [187] "PNRC1"    "PTX3"     "SOD2"     "CCL20"    "NFIL3"    "RIPK2"   
#> [193] "CXCL1"    "FOSL1"    "CEBPB"    "BTG3"     "ACKR3"    "CCN1"    
#> [199] "PLPP3"    "RIGI"    
#> 
#> $all$HALLMARK_MITOTIC_SPINDLE
#>   [1] "PREX1"    "NUMA1"    "ARFIP2"   "TAOK2"    "KIF3B"    "RASA1"   
#>   [7] "SHROOM2"  "FLNB"     "SHROOM1"  "RAPGEF6"  "ARHGEF12" "APC"     
#>  [13] "PCM1"     "ARHGEF3"  "RHOT2"    "NF1"      "CYTH2"    "ARFGEF1" 
#>  [19] "WASL"     "BCL2L11"  "TSC1"     "ABR"      "RICTOR"   "FGD6"    
#>  [25] "LATS1"    "EZR"      "RFC1"     "ITSN1"    "DYNLL2"   "TRIO"    
#>  [31] "CDC42"    "TUBGCP2"  "SPTAN1"   "ARHGAP29" "CLIP1"    "MYH10"   
#>  [37] "CD2AP"    "TUBGCP6"  "CSNK1D"   "TLK1"     "ARF6"     "PALLD"   
#>  [43] "CDC27"    "ABL1"     "CTTN"     "PXN"      "PPP4R2"   "PAFAH1B1"
#>  [49] "PCGF5"    "TUBGCP5"  "CNTROB"   "ARHGAP5"  "RABGAP1"  "PKD2"    
#>  [55] "PDLIM5"   "ALMS1"    "MID1IP1"  "ROCK1"    "ALS2"     "DST"     
#>  [61] "EPB41"    "NIN"      "CDC42EP2" "KIFAP3"   "ARL8A"    "MARK4"   
#>  [67] "STAU1"    "RAB3GAP1" "SSH2"     "CDK5RAP2" "HDAC6"    "TBCD"    
#>  [73] "SYNPO"    "KPTN"     "AKAP13"   "NEDD9"    "CEP250"   "RANBP9"  
#>  [79] "ARHGDIA"  "CDC42BPA" "KIF22"    "DOCK4"    "NOTCH2"   "RAPGEF5" 
#>  [85] "KIF5B"    "CDC42EP4" "BCAR1"    "CAPZB"    "OPHN1"    "ARHGEF11"
#>  [91] "TIAM1"    "LLGL1"    "KLC1"     "VCL"      "TUBD1"    "NET1"    
#>  [97] "DLG1"     "SMC3"     "DOCK2"    "FARP1"    "ATG4B"    "DYNC1H1" 
#> [103] "SMC1A"    "YWHAE"    "SORBS2"   "CEP192"   "ARHGAP4"  "MAP1S"   
#> [109] "HOOK3"    "ARAP3"    "KNTC1"    "GSN"      "WASF2"    "SPTBN1"  
#> [115] "RALBP1"   "ARHGEF7"  "PCNT"     "FGD4"     "UXT"      "MAP3K11" 
#> [121] "RASA2"    "TOP2A"    "CEP72"    "MAPRE1"   "MYH9"     "CKAP5"   
#> [127] "KIF20B"   "CLASP1"   "ARHGEF2"  "LRPPRC"   "RACGAP1"  "KATNB1"  
#> [133] "BRCA2"    "EPB41L2"  "ESPL1"    "ECT2"     "GEMIN4"   "SAC3D1"  
#> [139] "ARHGAP10" "KIF3C"    "RHOF"     "NEK2"     "TUBGCP3"  "WASF1"   
#> [145] "CLIP2"    "RASAL2"   "MARCKS"   "SASS6"    "FLNA"     "NUSAP1"  
#> [151] "PLEKHG2"  "LMNB1"    "BCR"      "CCDC88A"  "ARHGAP27" "SMC4"    
#> [157] "KIF1B"    "CDK1"     "BIN1"     "MYO9B"    "PRC1"     "ACTN4"   
#> [163] "STK38L"   "CEP57"    "KIF11"    "FBXO5"    "MYO1E"    "SUN2"    
#> [169] "KIF15"    "AURKA"    "CENPE"    "ABI1"     "KIF4A"    "KIF23"   
#> [175] "SOS1"     "CENPF"    "BIRC5"    "TUBA4A"   "INCENP"   "DLGAP5"  
#> [181] "NCK1"     "TPX2"     "PIF1"     "PLK1"     "BUB1"     "CCNB2"   
#> [187] "KATNA1"   "NDC80"    "ANLN"     "KIF2C"    "FSCN1"    "CDC42EP1"
#> [193] "MID1"     "TTK"      "NCK2"     "CEP131"   "CNTRL"    "CPAP"    
#> [199] "SEPTIN9" 
#> 
#> $all$HALLMARK_SPERMATOGENESIS
#>   [1] "PCSK4"   "IFT88"   "RAD17"   "STAM2"   "GSTM3"   "CLGN"    "NPHP1"  
#>   [8] "PHF7"    "TUBA3C"  "SPATA6"  "HSPA1L"  "HSPA2"   "TEKT2"   "JAM3"   
#>  [15] "ZC3H14"  "SIRT1"   "SLC12A2" "IP6K1"   "HOXB1"   "CHRM4"   "STRBP"  
#>  [22] "PHKG2"   "GRM8"    "SHE"     "IDE"     "TNNI3"   "COIL"    "PRKAR2A"
#>  [29] "GAD1"    "CCT6B"   "ACRBP"   "PGK2"    "PAPOLB"  "MAP7"    "GSG1"   
#>  [36] "PEBP1"   "ALOX15"  "GPR182"  "GMCL1"   "ACE"     "NPY5R"   "NEFH"   
#>  [43] "OAZ3"    "SCG5"    "PACRG"   "DPEP3"   "CAMK4"   "DDX25"   "PGS1"   
#>  [50] "BRAF"    "MTOR"    "YBX2"    "DMC1"    "PIAS2"   "CLVS1"   "ADCYAP1"
#>  [57] "NOS1"    "IL13RA2" "VDAC3"   "MAST2"   "MLF1"    "HSPA4L"  "SYCP1"  
#>  [64] "ARL4A"   "SCG3"    "TSN"     "PARP2"   "ELOVL3"  "LDHC"    "MLLT10" 
#>  [71] "TULP2"   "TCP11"   "TOPBP1"  "GAPDHS"  "TKTL1"   "DCC"     "NEK2"   
#>  [78] "CFTR"    "CLPB"    "CRISP2"  "CNIH2"   "GFI1"    "SNAP91"  "TALDO1" 
#>  [85] "CHFR"    "CCNA1"   "NF2"     "MTNR1A"  "AGFG1"   "CDK1"    "PCSK1N" 
#>  [92] "TLE4"    "CDKN3"   "AURKA"   "SLC2A5"  "RPL39L"  "ACRV1"   "POMC"   
#>  [99] "NCAPH"   "CSNK2A2" "EZH2"    "BUB1"    "CCNB2"   "DBF4"    "RFC4"   
#> [106] "KIF2C"   "PSMG1"   "TTK"     "DMRT1"   "IL12RB2" "LPIN1"   "ART3"   
#> [113] "ACTL7B"  "ADAD1"   "ADAM2"   "AKAP4"   "CST8"    "DDX4"    "DNAJB8" 
#> [120] "H1-6"    "HBZ"     "HTR5A"   "MEP1B"   "NAA11"   "ODF1"    "PDHA2"  
#> [127] "PRM2"    "SEPTIN4" "SPMAP2"  "TNP1"    "TNP2"    "TSSK2"   "ZC2HC1C"
#> [134] "ZNRF4"   "ZPBP"   
#> 
#> $all$HALLMARK_GLYCOLYSIS
#>   [1] "TFF3"     "FUT8"     "TPBG"     "STC2"     "COG2"     "CACNA1H" 
#>   [7] "CYB5A"    "ALG1"     "XYLT2"    "IDUA"     "SLC35A3"  "PFKFB1"  
#>  [13] "ANG"      "ARPP19"   "GLCE"     "POLR3K"   "MPI"      "FKBP4"   
#>  [19] "BIK"      "IER3"     "GPC1"     "GPC4"     "ALDOA"    "IL13RA1" 
#>  [25] "CITED2"   "RBCK1"    "ECD"      "GPC3"     "AGL"      "ENO2"    
#>  [31] "PAXIP1"   "PGAM2"    "B4GALT7"  "DCN"      "CHST1"    "GNPDA1"  
#>  [37] "STC1"     "COL5A1"   "GOT2"     "PYGL"     "NOL3"     "GMPPB"   
#>  [43] "VCAN"     "B3GALT6"  "QSOX1"    "ARTN"     "ALDH9A1"  "GALE"    
#>  [49] "IRS2"     "FAM162A"  "GCLC"     "PHKA2"    "PPP2CB"   "CTH"     
#>  [55] "CASP6"    "HK2"      "KIF2A"    "FBP2"     "AGRN"     "GUSB"    
#>  [61] "SDHC"     "CLDN3"    "MED24"    "MXI1"     "SOD1"     "EGLN3"   
#>  [67] "GFPT1"    "PGLS"     "COPB2"    "SRD5A3"   "AK3"      "B4GALT1" 
#>  [73] "SPAG4"    "TPST1"    "CD44"     "SDC2"     "GMPPA"    "ALDOB"   
#>  [79] "P4HA2"    "MERTK"    "GLRX"     "BPNT1"    "AKR1A1"   "P4HA1"   
#>  [85] "LHX9"     "CHPF"     "NDUFV3"   "PMM2"     "PRPS1"    "SLC25A10"
#>  [91] "DPYSL4"   "B3GAT3"   "ALDH7A1"  "ABCB6"    "CHPF2"    "PAM"     
#>  [97] "PGM2"     "NT5E"     "PC"       "LHPP"     "HDLBP"    "PKP2"    
#> [103] "ELF3"     "GNE"      "MIF"      "LDHA"     "GYS1"     "NDST3"   
#> [109] "HSPA5"    "ZNF292"   "SLC37A4"  "IDH1"     "HOMER1"   "HAX1"    
#> [115] "TGFBI"    "GYS2"     "CAPN5"    "GALK2"    "PPFIA4"   "RPE"     
#> [121] "SDC3"     "HS6ST2"   "B4GALT4"  "CHST12"   "GALK1"    "LDHC"    
#> [127] "SLC16A3"  "DLD"      "GOT1"     "EFNA3"    "ANKZF1"   "GAL3ST1" 
#> [133] "GPR87"    "SDC1"     "G6PD"     "MIOX"     "PSMC4"    "GAPDHS"  
#> [139] "KDELR3"   "CLN6"     "TKTL1"    "PPIA"     "TXN"      "PGAM1"   
#> [145] "PYGB"     "PLOD2"    "ISG20"    "IGFBP3"   "HS2ST1"   "PDK3"    
#> [151] "ANGPTL4"  "CXCR4"    "TALDO1"   "HMMR"     "NSDHL"    "ME1"     
#> [157] "PGK1"     "NANP"     "CLDN9"    "EXT2"     "MDH2"     "SAP30"   
#> [163] "CDK1"     "VEGFA"    "DDIT4"    "SOX9"     "CHST6"    "AURKA"   
#> [169] "B3GAT1"   "B4GALT2"  "EXT1"     "KIF20A"   "ME2"      "TPI1"    
#> [175] "ADORA2B"  "B3GNT3"   "MET"      "RRAGD"    "TGFA"     "VLDLR"   
#> [181] "SLC25A13" "DEPDC1"   "MDH1"     "CHST4"    "PFKP"     "CHST2"   
#> [187] "EGFR"     "STMN1"    "NASP"     "CENPA"    "ENO1"     "UGP2"    
#> [193] "LCT"      "PLOD1"    "DSC2"     "AK4"      "ERO1A"    "GFUS"    
#> [199] "PKM"      "RARS1"   
#> 
#> $all$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
#>  [1] "HHEX"    "EGLN2"   "NQO1"    "GPX4"    "CAT"     "SRXN1"   "ATOX1"  
#>  [8] "GSR"     "PRDX2"   "TXNRD2"  "NDUFB4"  "JUNB"    "PDLIM1"  "GCLC"   
#> [15] "ERCC2"   "HMOX2"   "SOD1"    "GPX3"    "NDUFS2"  "LSP1"    "SBNO2"  
#> [22] "GLRX"    "OXSR1"   "MSRA"    "FES"     "GCLM"    "CDKN2D"  "MBP"    
#> [29] "IPCEF1"  "FTL"     "MGST1"   "G6PD"    "PRDX6"   "PRDX1"   "ABCC1"  
#> [36] "TXN"     "TXNRD1"  "NDUFA6"  "STK25"   "MPO"     "GLRX2"   "PRNP"   
#> [43] "SOD2"    "PFKP"    "PRDX4"   "LAMTOR5" "PTPA"    "SCAF4"   "SELENOS"
#> 
#> $all$HALLMARK_UNFOLDED_PROTEIN_RESPONSE
#>   [1] "XBP1"     "WFS1"     "STC2"     "SLC30A5"  "SLC1A4"   "SPCS1"   
#>   [7] "DNAJA4"   "EXOC2"    "EDEM1"    "CNOT2"    "BAG3"     "DCP2"    
#>  [13] "HSPA9"    "ERN1"     "PARN"     "DCP1A"    "EEF2"     "TSPYL2"  
#>  [19] "ATF6"     "CNOT6"    "DCTN1"    "RPS14"    "GOSR2"    "SPCS3"   
#>  [25] "DNAJB9"   "EIF4A2"   "IMP3"     "EIF4E"    "EIF2AK3"  "IFIT1"   
#>  [31] "SHC1"     "IGFBP1"   "NOP14"    "SERP1"    "NHP2"     "CNOT4"   
#>  [37] "DNAJC3"   "WIPI1"    "SEC31A"   "ARFGAP1"  "LSM1"     "ATP6V0D1"
#>  [43] "NFYB"     "HERPUD1"  "KIF5B"    "EXOSC10"  "EIF4A3"   "EXOSC1"  
#>  [49] "NPM1"     "SEC11A"   "FUS"      "ZBTB17"   "SDAD1"    "CXXC1"   
#>  [55] "EIF2S1"   "YIF1A"    "EXOSC9"   "TUBB2A"   "FKBP14"   "PAIP1"   
#>  [61] "HSPA5"    "EXOSC5"   "ATF3"     "EXOSC2"   "EIF4A1"   "EDC4"    
#>  [67] "RRP9"     "EXOSC4"   "NOLC1"    "ALDH18A1" "BANF1"    "LSM4"    
#>  [73] "TATDN2"   "KDELR3"   "GEMIN4"   "YWHAZ"    "POP4"     "NOP56"   
#>  [79] "SSR1"     "CCL2"     "SRPRB"    "PREB"     "EIF4G1"   "XPOT"    
#>  [85] "HSP90B1"  "DDX10"    "NFYA"     "VEGFA"    "DDIT4"    "HYOU1"   
#>  [91] "KHSRP"    "EIF4EBP1" "PDIA5"    "ATF4"     "MTHFD2"   "CALR"    
#>  [97] "CEBPG"    "CKS1B"    "ASNS"     "CHAC1"    "SLC7A5"   "CEBPB"   
#> [103] "PDIA6"    "DKC1"     "PSAT1"    "ERO1A"    "H2AX"     "IARS1"   
#> [109] "MTREX"    "NABP1"    "SKIC3"    "SRPRA"    "TARS1"   
#> 
#> $all$HALLMARK_BILE_ACID_METABOLISM
#>   [1] "NUDT12"   "AR"       "HSD17B4"  "PEX12"    "EFHC1"    "ABCA3"   
#>   [7] "SLC27A2"  "ABCA2"    "PEX11A"   "PEX19"    "ABCD3"    "SLC23A1" 
#>  [13] "CROT"     "PEX11G"   "PEX7"     "SERPINA6" "SULT2B1"  "DIO1"    
#>  [19] "RXRA"     "DHCR24"   "FDXR"     "NEDD4"    "HSD3B7"   "PECR"    
#>  [25] "RETSAT"   "GNMT"     "LIPE"     "AMACR"    "LONP2"    "EPHX2"   
#>  [31] "PAOX"     "ISOC1"    "BCAR3"    "ABCA8"    "ATXN1"    "HACL1"   
#>  [37] "ABCA9"    "PFKM"     "PEX1"     "APOA1"    "ABCA6"    "PXMP2"   
#>  [43] "SCP2"     "CYP46A1"  "CAT"      "RXRG"     "ABCA5"    "HSD17B6" 
#>  [49] "PEX16"    "ABCA1"    "HSD17B11" "ALDH9A1"  "ALDH1A1"  "TTR"     
#>  [55] "AGXT"     "ALDH8A1"  "SLC27A5"  "MLYCD"    "SLC29A1"  "PNPLA8"  
#>  [61] "HSD3B1"   "SOAT2"    "CH25H"    "ABCA4"    "NR3C2"    "NR0B2"   
#>  [67] "SOD1"     "DIO2"     "AKR1D1"   "GSTK1"    "PEX6"     "GNPAT"   
#>  [73] "PEX26"    "SLC23A2"  "PHYH"     "PIPOX"    "CYP7A1"   "ABCG8"   
#>  [79] "SULT1B1"  "ABCG4"    "IDH1"     "PRDX5"    "FADS1"    "ABCD2"   
#>  [85] "GCLM"     "ACSL1"    "ACSL5"    "BMP6"     "FADS2"    "CYP27A1" 
#>  [91] "CYP8B1"   "IDI1"     "RBP1"     "TFCP2L1"  "NR1I2"    "KLF1"    
#>  [97] "LCK"      "IDH2"     "SLCO1A2"  "PEX13"    "SLC35B2"  "ABCD1"   
#> [103] "AQP9"     "OPTN"     "NPC1"     "CYP7B1"   "CYP39A1"  "BBOX1"   
#> [109] "HAO1"     "GC"       "NR1H4"    "SLC67A1" 
#> 
#> $all$HALLMARK_WNT_BETA_CATENIN_SIGNALING
#>  [1] "FRAT1"  "HDAC11" "AXIN2"  "LEF1"   "PSEN2"  "NOTCH4" "NCOR2"  "MAML1" 
#>  [9] "NUMB"   "JAG1"   "HEY1"   "TP53"   "DLL1"   "KAT2A"  "HDAC5"  "CCND2" 
#> [17] "HEY2"   "FZD8"   "RBPJ"   "DKK4"   "NCSTN"  "AXIN1"  "CTNNB1" "JAG2"  
#> [25] "NKD1"   "FZD1"   "WNT1"   "DVL2"   "MYC"    "WNT5B"  "CSNK1E" "TCF7"  
#> [33] "CUL1"   "PTCH1"  "GNAI1"  "ADAM17" "DKK1"   "NOTCH1" "PPARD"  "HDAC2" 
#> [41] "SKP2"   "WNT6"  
#> 
#> $all$HALLMARK_HYPOXIA
#>   [1] "CA12"     "FBP1"     "BCL2"     "BHLHE40"  "TPBG"     "STC2"    
#>   [7] "SIAH2"    "CCNG2"    "TGFB3"    "PPP1R3C"  "SELENBP1" "MAP3K1"  
#>  [13] "SULT2B1"  "NEDD4L"   "STBD1"    "TIPARP"   "BNIP3L"   "ATP7A"   
#>  [19] "IER3"     "GPC1"     "AKAP12"   "SDC4"     "GPC4"     "DUSP1"   
#>  [25] "ALDOA"    "NDST2"    "CITED2"   "GPC3"     "TPST2"    "HK1"     
#>  [31] "CDKN1B"   "IDS"      "ENO2"     "FOS"      "TPD52"    "PGAM2"   
#>  [37] "NDST1"    "DCN"      "PCK1"     "KIF5A"    "PFKFB3"   "STC1"    
#>  [43] "COL5A1"   "CDKN1A"   "PDGFB"    "SLC6A6"   "GCK"      "B3GALT6" 
#>  [49] "RORA"     "SERPINE1" "ENO3"     "NR3C1"    "CAV1"     "ZFP36"   
#>  [55] "VHL"      "ILVBL"    "IGFBP1"   "IRS2"     "FAM162A"  "PYGM"    
#>  [61] "GAA"      "CASP6"    "HK2"      "AMPD3"    "BCAN"     "PKLR"    
#>  [67] "GRHPR"    "JUN"      "F3"       "BGN"      "FOSL2"    "WSB1"    
#>  [73] "INHA"     "DTNA"     "KLF7"     "SRPX"     "PGF"      "MXI1"    
#>  [79] "GBE1"     "HEXA"     "RBPJ"     "HMOX1"    "LOX"      "SDC2"    
#>  [85] "FOXO3"    "ALDOB"    "P4HA2"    "LXN"      "GLRX"     "B4GALNT2"
#>  [91] "P4HA1"    "SLC25A1"  "CDKN1C"   "MT1E"     "BTG1"     "TGM2"    
#>  [97] "DPYSL4"   "DDIT3"    "HAS1"     "PAM"      "EFNA1"    "KLF6"    
#> [103] "PGM2"     "HDLBP"    "MIF"      "ALDOC"    "LDHA"     "GYS1"    
#> [109] "HSPA5"    "ZNF292"   "NAGK"     "SLC37A4"  "SCARB1"   "KDM3A"   
#> [115] "PRDX5"    "TGFBI"    "ATF3"     "PPFIA4"   "ANXA2"    "PPP1R15A"
#> [121] "HOXB9"    "MT2A"     "SDC3"     "SLC2A3"   "GALK1"    "LDHC"    
#> [127] "EFNA3"    "MYH9"     "EDN2"     "ANKZF1"   "IL6"      "GAPDHS"  
#> [133] "KDELR3"   "TKTL1"    "PLAC8"    "S100A4"   "ETS1"     "PFKL"    
#> [139] "ISG20"    "IGFBP3"   "PDK3"     "ANGPTL4"  "ERRFI1"   "PHKG1"   
#> [145] "CXCR4"    "MAFF"     "TNFAIP3"  "JMJD6"    "SLC2A1"   "PPARGC1A"
#> [151] "HS3ST1"   "PGK1"     "KLHL24"   "SAP30"    "TMEM45A"  "VEGFA"   
#> [157] "CP"       "DDIT4"    "PRKCA"    "ADM"      "SLC2A5"   "GPI"     
#> [163] "PKP1"     "PLAUR"    "EXT1"     "NCAN"     "NDRG1"    "GAPDH"   
#> [169] "TPI1"     "ADORA2B"  "PLIN2"    "PGM1"     "RRAGD"    "PNRC1"   
#> [175] "VLDLR"    "CHST3"    "XPNPEP1"  "PIM1"     "GCNT2"    "TES"     
#> [181] "NFIL3"    "PFKP"     "CHST2"    "EGFR"     "PDK1"     "CSRP2"   
#> [187] "ENO1"     "UGP2"     "ACKR3"    "AK4"      "BRS3"     "CAVIN1"  
#> [193] "CAVIN3"   "CCN1"     "CCN2"     "CCN5"     "ERO1A"    "LALBA"   
#> [199] "LARGE1"   "NOCT"    
#> 
#> $all$HALLMARK_UV_RESPONSE_UP
#>   [1] "NAT1"     "RHOB"     "TMBIM6"   "IL6ST"    "BTG2"     "CYB5R1"  
#>   [7] "RET"      "SULT1A1"  "IGFBP2"   "SPR"      "PRKCD"    "DLG4"    
#>  [13] "OLFM1"    "HSPA2"    "PLCL1"    "FKBP4"    "MAOA"     "MAPK8IP2"
#>  [19] "PPT1"     "BCL2L11"  "FGF18"    "ALDOA"    "HTR7"     "ENO2"    
#>  [25] "FOS"      "EIF5"     "GRINA"    "MSX1"     "FOSB"     "ATP6V1C1"
#>  [31] "SQSTM1"   "TACR3"    "EPHX1"    "PPP1R2"   "SLC25A4"  "ABCB1"   
#>  [37] "POLG2"    "CDO1"     "APOM"     "TGFBRAP1" "ACAA1"    "CA2"     
#>  [43] "CYP1A1"   "JUNB"     "PTPRD"    "DGAT1"    "C4BPB"    "NPTXR"   
#>  [49] "HYAL2"    "TST"      "MRPL23"   "NR4A1"    "COL2A1"   "SLC6A12" 
#>  [55] "PRKACA"   "GRPEL1"   "UROD"     "NFKBIA"   "GPX3"     "STARD3"  
#>  [61] "NXF1"     "RASGRP1"  "HSPA13"   "HNRNPU"   "CLCN2"    "RRAD"    
#>  [67] "NTRK3"    "DNAJB1"   "BSG"      "HMOX1"    "MMP14"    "TCHH"    
#>  [73] "EIF2S3"   "AQP3"     "PDLIM3"   "CDKN1C"   "BTG1"     "FURIN"   
#>  [79] "RAB27A"   "KCNH2"    "SHOX2"    "IRF1"     "CDKN2B"   "NPTX2"   
#>  [85] "RXRB"     "MARK2"    "CLTB"     "AP2S1"    "FMO1"     "PPAT"    
#>  [91] "PARP2"    "ATF3"     "CREG1"    "GCH1"     "CYB5B"    "POLE3"   
#>  [97] "WIZ"      "CDK2"     "MGAT1"    "CASP3"    "DDX21"    "POLR2H"  
#> [103] "PSMC3"    "ATP6V1F"  "ALAS1"    "CCND3"    "IL6"      "E2F5"    
#> [109] "DNAJA1"   "ARRB2"    "TFRC"     "CHKA"     "PPIF"     "STK25"   
#> [115] "STIP1"    "CNP"      "CXCL2"    "LHX2"     "SLC6A8"   "YKT6"    
#> [121] "CDC34"    "HLA-F"    "SIGMAR1"  "TAP1"     "CCK"      "CDC5L"   
#> [127] "BMP2"     "ICAM1"    "PRPF3"    "NKX2-5"   "FEN1"     "RPN1"    
#> [133] "PDAP1"    "GGH"      "TUBA4A"   "TYRO3"    "EPCAM"    "BID"     
#> [139] "BAK1"     "CEBPG"    "GLS"      "SOD2"     "CHRNA5"   "LYN"     
#> [145] "RFC4"     "ASNS"     "KLHDC3"   "AMD1"     "BTG3"     "GAL"     
#> [151] "CCNE1"    "AGO2"     "CTSV"     "H2AX"     "NUP58"    "ONECUT1" 
#> [157] "SELENOW"  "TARS1"   
#> 
#> $all$HALLMARK_CHOLESTEROL_HOMEOSTASIS
#>  [1] "CPEB2"     "TP53INP1"  "SEMA3B"    "ABCA2"     "ALCAM"     "TM7SF2"   
#>  [7] "FASN"      "ATXN2"     "GSTM2"     "HSD17B7"   "HMGCR"     "SCD"      
#> [13] "TRIB3"     "PMVK"      "CYP51A1"   "CLU"       "JAG1"      "PPARG"    
#> [19] "AVPR1A"    "CD9"       "GPX8"      "ANXA13"    "ECH1"      "MAL2"     
#> [25] "GUSB"      "ANTXR2"    "ETHE1"     "FDFT1"     "TMEM97"    "LDLR"     
#> [31] "STX5"      "FBXO6"     "LPL"       "ACSS2"     "CTNNB1"    "ACTG1"    
#> [37] "MVK"       "LGALS3"    "ADH4"      "SREBF2"    "SQLE"      "ALDOC"    
#> [43] "HMGCS1"    "ATF3"      "LSS"       "GLDC"      "FADS2"     "TNFRSF12A"
#> [49] "IDI1"      "ACAT2"     "CHKA"      "PCYT2"     "ATF5"      "ANXA5"    
#> [55] "PDK3"      "ERRFI1"    "LGMN"      "PLSCR1"    "NSDHL"     "STARD4"   
#> [61] "EBP"       "GNAI1"     "DHCR7"     "MVD"       "CXCL16"    "PLAUR"    
#> [67] "FABP5"     "FDPS"      "S100A11"   "CBS"       "PNRC1"     "NFIL3"    
#> [73] "NIBAN1"    "SC5D"     
#> 
#> $all$HALLMARK_KRAS_SIGNALING_UP
#>   [1] "TSPAN1"   "TSPAN13"  "SEMA3B"   "PTCD2"    "KIF5C"    "CROT"    
#>   [7] "FUCA1"    "CAB39L"   "SCN1B"    "DNMBP"    "PLAT"     "CBR4"    
#>  [13] "MAP3K1"   "MTMR10"   "FBXO4"    "CPE"      "SPARCL1"  "GADD45G" 
#>  [19] "CFB"      "ANO1"     "SERPINA3" "ITGBL1"   "IGF2"     "AKAP12"  
#>  [25] "F13A1"    "CDADC1"   "ITGA2"    "NAP1L2"   "EPB41L3"  "MMP11"   
#>  [31] "GNG11"    "MMP10"    "VWA5A"    "TPH1"     "EVI5"     "ABCB1"   
#>  [37] "SPON1"    "JUP"      "CIDEA"    "TSPAN7"   "KLF4"     "ANXA10"  
#>  [43] "USP12"    "APOD"     "ALDH1A2"  "MAP7"     "CA2"      "AKT2"    
#>  [49] "FLT4"     "AVL9"     "PRKG2"    "SNAP25"   "PIGR"     "PRRX1"   
#>  [55] "INHBA"    "CBX8"     "PECAM1"   "PLVAP"    "MAFB"     "CFH"     
#>  [61] "ETV1"     "ACE"      "NIN"      "C3AR1"    "SDCCAG8"  "PTPRR"   
#>  [67] "RBP4"     "RBM4"     "SCG5"     "TOR1AIP2" "IL33"     "TRIB1"   
#>  [73] "DUSP6"    "RGS16"    "NRP1"     "CMKLR1"   "CCND2"    "NR0B2"   
#>  [79] "PEG3"     "ATG10"    "F2RL1"    "LAT2"     "ARG1"     "MMD"     
#>  [85] "ENG"      "TFPI"     "RELN"     "ZNF277"   "TMEM100"  "NGF"     
#>  [91] "RABGAP1L" "GLRX"     "G0S2"     "STRN"     "PLAU"     "PLEK2"   
#>  [97] "PSMB8"    "DOCK2"    "SPP1"     "ADAM8"    "TMEM176A" "EMP1"    
#> [103] "IRF8"     "LAPTM5"   "BPGM"     "CD37"     "GFPT2"    "TMEM176B"
#> [109] "CLEC4A"   "IL10RA"   "EREG"     "LY96"     "TLR8"     "HSD11B1" 
#> [115] "SCG3"     "IL7R"     "IL1B"     "PRDM1"    "FCER1G"   "BTC"     
#> [121] "MMP9"     "PPP1R15A" "CSF2RA"   "ID2"      "RETN"     "AMMECR1" 
#> [127] "GYPC"     "GPNMB"    "PCP4"     "HBEGF"    "LCP1"     "CBL"     
#> [133] "WNT7A"    "ITGB2"    "SATB1"    "PPBP"     "ANKH"     "TNFRSF1B"
#> [139] "IKZF1"    "TRIB2"    "MAP4K1"   "CTSS"     "SPRY2"    "EPHB2"   
#> [145] "HOXD11"   "LIF"      "GABRA3"   "BIRC3"    "ETS1"     "WDR33"   
#> [151] "HKDC1"    "PTBP2"    "IGFBP3"   "DCBLD2"   "IL2RG"    "TRAF1"   
#> [157] "USH1C"    "HDAC9"    "ANGPTL4"  "PDCD1LG2" "ZNF639"   "CXCR4"   
#> [163] "SNAP91"   "TNFAIP3"  "MYCN"     "ST6GAL1"  "ADAMDEC1" "FGF9"    
#> [169] "BMP2"     "YRDC"     "PCSK1N"   "PTGS2"    "IL1RL2"   "SOX9"    
#> [175] "GALNT3"   "CXCL10"   "ALDH1A3"  "CSF2"     "ETV5"     "ETV4"    
#> [181] "ADAM17"   "PLAUR"    "TNNT2"    "BTBD3"    "TMEM158"  "GPRC5B"  
#> [187] "CCL20"    "SLPI"     "MPZL2"    "MALL"     "KCNN4"    "ADGRA2"  
#> [193] "ADGRL4"   "CCSER2"   "CFHR2"    "ERO1A"    "GUCY1A1"  "H2BC3"   
#> [199] "NR1H4"    "PRELID3B"
#> 
#> $all$HALLMARK_NOTCH_SIGNALING
#>  [1] "SKP1"    "CCND1"   "LFNG"    "ARRB1"   "FBXW11"  "PSEN2"   "HEYL"   
#>  [8] "WNT5A"   "JAG1"    "HES1"    "DLL1"    "PSENEN"  "KAT2A"   "NOTCH2" 
#> [15] "NOTCH3"  "FZD1"    "WNT2"    "DTX1"    "APH1A"   "FZD5"    "RBX1"   
#> [22] "DTX2"    "SAP30"   "CUL1"    "FZD7"    "DTX4"    "TCF7L2"  "PRKCA"  
#> [29] "NOTCH1"  "MAML2"   "PPARD"   "ST3GAL6"
#> 
#> $all$HALLMARK_COMPLEMENT
#>   [1] "GATA3"    "F7"       "TMPRSS6"  "CD46"     "DUSP5"    "PLAT"    
#>   [7] "PRKCD"    "AKAP10"   "CTSO"     "ANG"      "USP8"     "CFB"     
#>  [13] "PSEN1"    "ZEB1"     "PCLO"     "CASP7"    "PRDM4"    "SERPINC1"
#>  [19] "GRB2"     "GMFB"     "LIPA"     "ZFPM2"    "CASP9"    "RNF4"    
#>  [25] "LTA4H"    "GPD2"     "LRP1"     "NOTCH4"   "DOCK10"   "IRF7"    
#>  [31] "BRPF3"    "F10"      "CD36"     "PDGFB"    "ATOX1"    "F8"      
#>  [37] "IRF2"     "HSPA1A"   "GNAI2"    "VCPIP1"   "PRCP"     "S100A13" 
#>  [43] "KCNIP2"   "SERPINE1" "DOCK9"    "FN1"      "SERPINA1" "CDH13"   
#>  [49] "CLU"      "STX4"     "MMP13"    "CA2"      "TFPI2"    "USP14"   
#>  [55] "PPP2CB"   "OLR1"     "TIMP2"    "CTSD"     "C4BPB"    "CFH"     
#>  [61] "PRSS36"   "TIMP1"    "SH2B3"    "KIF2A"    "F3"       "HNF4A"   
#>  [67] "CALM3"    "ITGAM"    "DUSP6"    "JAK2"     "ERAP2"    "GNG2"    
#>  [73] "PPP4C"    "CDA"      "USP16"    "RASGRP1"  "CALM1"    "PIK3CG"  
#>  [79] "KCNIP3"   "PIK3R5"   "DOCK4"    "SIRT6"    "GNB2"     "USP15"   
#>  [85] "MMP14"    "APOC1"    "KLKB1"    "DPP4"     "APOBEC3F" "LAP3"    
#>  [91] "RCE1"     "ADAM9"    "CTSH"     "ITIH1"    "CD55"     "CSRP1"   
#>  [97] "LGALS3"   "GZMK"     "RABIF"    "PDP1"     "LTF"      "DGKH"    
#> [103] "GP1BA"    "C2"       "APOBEC3G" "IRF1"     "CPM"      "SERPING1"
#> [109] "DYRK2"    "PHEX"     "C1QC"     "CDK5R1"   "CTSB"     "HSPA5"   
#> [115] "MT3"      "GNGT2"    "HPCAL4"   "LCP2"     "C3"       "ACTN2"   
#> [121] "SCG3"     "PIK3CA"   "C1QA"     "CD40LG"   "MMP8"     "PLEK"    
#> [127] "FCER1G"   "GCA"      "MMP15"    "WAS"      "RAF1"     "CASP10"  
#> [133] "CBLB"     "KLK1"     "GZMA"     "C1S"      "CASP3"    "SRC"     
#> [139] "COL4A2"   "IL6"      "FCN1"     "CASP4"    "CR1"      "CASP1"   
#> [145] "C1R"      "S100A12"  "PSMB9"    "RHOG"     "CTSS"     "SERPINB2"
#> [151] "PCSK9"    "CD59"     "ADRA2B"   "ANXA5"    "F5"       "GNAI3"   
#> [157] "FYN"      "LAMP2"    "LCK"      "MAFF"     "LGMN"     "TNFAIP3" 
#> [163] "CCL5"     "PLSCR1"   "EHD1"     "ME1"      "KYNU"     "CR2"     
#> [169] "PLA2G7"   "GNB4"     "SPOCK2"   "FDX1"     "CP"       "CASP5"   
#> [175] "PLAUR"    "GZMB"     "PFN1"     "DGKG"     "MMP12"    "PREP"    
#> [181] "PRSS3"    "PLA2G4A"  "CTSC"     "XPNPEP1"  "PIM1"     "S100A9"  
#> [187] "LYN"      "L3MBTL4"  "CXCL1"    "CEBPB"    "APOA4"    "C9"      
#> [193] "CPQ"      "CTSL"     "CTSV"     "F2"       "GP9"      "MSRB1"   
#> [199] "PLG"      "RBSN"    
#> 
#> $all$HALLMARK_IL2_STAT5_SIGNALING
#>   [1] "XBP1"      "RHOB"      "BCL2"      "ENPP1"     "IGF1R"     "BHLHE40"  
#>   [7] "LRIG1"     "CISH"      "AHNAK"     "P2RX4"     "SERPINB6"  "FAH"      
#>  [13] "RHOH"      "ALCAM"     "BCL2L1"    "SNX9"      "MUC1"      "IKZF4"    
#>  [19] "SPRED2"    "ECM1"      "F2RL2"     "GATA1"     "AMACR"     "SNX14"    
#>  [25] "BATF"      "SERPINC1"  "PTH1R"     "IL1RL1"    "RNH1"      "TWSG1"    
#>  [31] "SOCS2"     "GPX4"      "MYO1C"     "DENND5A"   "MAP6"      "SYNGR2"   
#>  [37] "SHE"       "BMPR2"     "MAPKAPK2"  "TNFRSF18"  "SH3BGRL2"  "NCOA3"    
#>  [43] "ABCB1"     "PTGER2"    "TNFSF10"   "ITIH5"     "SELP"      "RORA"     
#>  [49] "WLS"       "ITGAV"     "ENO3"      "PRAF2"     "COL6A1"    "CA2"      
#>  [55] "GADD45B"   "PLEC"      "NFKBIZ"    "TLR7"      "POU2F1"    "AHR"      
#>  [61] "HK2"       "IKZF2"     "SCN9A"     "HIPK2"     "FLT3LG"    "HOPX"     
#>  [67] "FGL2"      "CSF1"      "SPRY4"     "PRKCH"     "PENK"      "SYT11"    
#>  [73] "CYFIP1"    "SMPDL3A"   "IRF6"      "RGS16"     "NRP1"      "CCND2"    
#>  [79] "PHTF2"     "IFITM3"    "GPR65"     "IL4R"      "SWAP70"    "IL3RA"    
#>  [85] "CCR4"      "GALM"      "TNFSF11"   "CDC42SE2"  "MXD1"      "PTRH2"    
#>  [91] "CD44"      "RABGAP1L"  "GSTO1"     "P4HA1"     "CDKN1C"    "TIAM1"    
#>  [97] "TGM2"      "FURIN"     "CD81"      "SLC29A2"   "APLP1"     "SPP1"     
#> [103] "KLF6"      "EMP1"      "NT5E"      "IRF8"      "SELL"      "SLC39A8"  
#> [109] "IL10RA"    "AGER"      "GBP4"      "CAPN3"     "IL10"      "MAP3K8"   
#> [115] "ARL4A"     "DHRS3"     "CAPG"      "CD48"      "SLC1A5"    "EOMES"    
#> [121] "TTC39B"    "CD86"      "IGF2R"     "GABARAPL1" "GPR83"     "PHLDA1"   
#> [127] "IL2RB"     "AHCY"      "SLC2A3"    "HUWE1"     "MYC"       "CST7"     
#> [133] "CTSZ"      "CASP3"     "UMPS"      "TNFRSF1B"  "CCND3"     "TNFRSF4"  
#> [139] "CDCP1"     "IL13"      "LTB"       "ITGAE"     "LIF"       "LCLAT1"   
#> [145] "ADAM19"    "ST3GAL4"   "ANXA4"     "S100A1"    "TRAF1"     "CDC6"     
#> [151] "IL18R1"    "TNFRSF8"   "TNFRSF9"   "MAFF"      "CKAP4"     "PLSCR1"   
#> [157] "ICOS"      "PUS1"      "LRRC8C"    "COCH"      "IRF4"      "CTLA4"    
#> [163] "CD83"      "CD79B"     "ITGA6"     "BMP2"      "BATF3"     "PDCD2L"   
#> [169] "PTCH1"     "CXCL10"    "SOCS1"     "CSF2"      "MYO1E"     "IFNGR1"   
#> [175] "ETV4"      "IL2RA"     "PLAGL1"    "PRNP"      "NDRG1"     "PLIN2"    
#> [181] "PNP"       "RRAGD"     "IL1R2"     "PIM1"      "DCPS"      "NOP2"     
#> [187] "NFIL3"     "TNFRSF21"  "GLIPR2"    "NCS1"      "ODC1"      "UCK2"     
#> [193] "CCNE1"     "DRC1"      "EEF1AKMT1" "ETFBKMT"   "GUCY1B1"   "HYCC2"    
#> [199] "PLPP1"    
#> 
#> $all$HALLMARK_INTERFERON_ALPHA_RESPONSE
#>  [1] "IFITM2"   "ELF1"     "PSME1"    "DHX58"    "LPAR6"    "IFI35"   
#>  [7] "TRAFD1"   "TDRD7"    "TXNIP"    "UBA7"     "BST2"     "SP110"   
#> [13] "IRF7"     "RNF31"    "IRF9"     "DDX60"    "IRF2"     "PARP9"   
#> [19] "PROCR"    "IFITM1"   "LGALS3BP" "UBE2L6"   "OGFR"     "SAMD9"   
#> [25] "OAS1"     "CSF1"     "PARP14"   "ADAR"     "CASP8"    "IFITM3"  
#> [31] "TMEM140"  "IL4R"     "TRIM14"   "STAT2"    "TRIM21"   "ISG15"   
#> [37] "IFI27"    "CCRL2"    "PSME2"    "TRIM5"    "MOV10"    "SAMD9L"  
#> [43] "IFIT3"    "SLC25A28" "IFIT2"    "LAP3"     "CD74"     "NUB1"    
#> [49] "TRIM25"   "GMPR"     "PSMB8"    "HERC6"    "MX1"      "B2M"     
#> [55] "IRF1"     "OASL"     "SELL"     "IFI44L"   "GBP2"     "RSAD2"   
#> [61] "HLA-C"    "GBP4"     "PARP12"   "PSMA3"    "CMPK2"    "C1S"     
#> [67] "RTP4"     "EPSTI1"   "USP18"    "IL15"     "BATF2"    "CASP1"   
#> [73] "EIF2AK2"  "PSMB9"    "IL7"      "IFI44"    "IFI30"    "IFIH1"   
#> [79] "LY6E"     "ISG20"    "CNP"      "CD47"     "CXCL11"   "PLSCR1"  
#> [85] "TRIM26"   "TAP1"     "PNPT1"    "CXCL10"   "NMI"      "LAMP3"   
#> [91] "RIPK2"    "NCOA7"    "CMTR1"    "HELZ2"    "MVB12A"   "TENT5A"  
#> [97] "WARS1"   
#> 
#> $all$HALLMARK_PI3K_AKT_MTOR_SIGNALING
#>   [1] "VAV3"     "PRKAG1"   "PLA2G12A" "CAB39L"   "MAPK9"    "PTEN"    
#>   [7] "TBK1"     "GNA14"    "MKNK2"    "TSC2"     "ATF1"     "IRAK4"   
#>  [13] "DUSP3"    "GRB2"     "CDKN1B"   "CLTC"     "YWHAB"    "ACACA"   
#>  [19] "RALB"     "SQSTM1"   "ITPR2"    "CDKN1A"   "PRKAR2A"  "TRIB3"   
#>  [25] "PRKAA2"   "PIN1"     "RPTOR"    "RIPK1"    "THEM4"    "PIKFYVE" 
#>  [31] "EIF4E"    "CAB39"    "UBE2D3"   "AKT1"     "PIK3R3"   "AKT1S1"  
#>  [37] "ECSIT"    "MAPK10"   "PAK4"     "HRAS"     "UBE2N"    "GSK3B"   
#>  [43] "ARPC3"    "CAMK4"    "PLCB1"    "STAT2"    "ARHGDIA"  "MAPK8"   
#>  [49] "PITX2"    "PTPN11"   "GNGT1"    "RIT1"     "NGF"      "TRAF2"   
#>  [55] "PPP1CA"   "SMAD2"    "RAC1"     "MAPK1"    "TIAM1"    "DDIT3"   
#>  [61] "TNFRSF1A" "SLA"      "PRKCB"    "ADCY2"    "PLCG1"    "ARF1"    
#>  [67] "RAF1"     "MKNK1"    "CDK2"     "MAPKAP1"  "MAP2K6"   "SFN"     
#>  [73] "FGF17"    "FASLG"    "RPS6KA1"  "NOD1"     "PPP2R1B"  "MAP3K7"  
#>  [79] "CFL1"     "IL2RG"    "MYD88"    "CXCR4"    "NFKBIB"   "LCK"     
#>  [85] "CDK4"     "HSP90B1"  "SLC2A1"   "DAPP1"    "AP2M1"    "CDK1"    
#>  [91] "CSNK2B"   "E2F1"     "ACTR2"    "RPS6KA3"  "MAP2K3"   "PFN1"    
#>  [97] "CALR"     "NCK1"     "ACTR3"    "EGFR"     "PDK1"     "FGF22"   
#> [103] "FGF6"     "GRK2"     "IL4"     
#> 
#> $all$HALLMARK_UV_RESPONSE_DN
#>   [1] "APBB2"    "INPP4B"   "IGF1R"    "BHLHE40"  "AGGF1"    "IRS1"    
#>   [7] "SPOP"     "CDON"     "SMAD3"    "KCNMA1"   "RXRA"     "PTEN"    
#>  [13] "GJA1"     "AMPH"     "ZMIZ1"    "DBP"      "IGFBP5"   "MAGI2"   
#>  [19] "CAP2"     "ICA1"     "PDGFRB"   "SYNJ2"    "YTHDC1"   "MAP2K5"  
#>  [25] "DUSP1"    "PMP22"    "MMP16"    "MGMT"     "CITED2"   "LPAR1"   
#>  [31] "ATXN1"    "CDKN1B"   "RUNX1"    "ATRX"     "PIAS3"    "COL1A2"  
#>  [37] "RBPMS"    "COL1A1"   "PTPRM"    "DDAH1"    "COL3A1"   "SFMBT1"  
#>  [43] "MIOS"     "ATRN"     "MGLL"     "NIPBL"    "COL5A2"   "SYNE1"   
#>  [49] "PRKAR2B"  "TGFBR3"   "DLC1"     "BDNF"     "PRDM2"    "ATP2B1"  
#>  [55] "RGS4"     "PRKCE"    "SMAD7"    "TJP1"     "PHF3"     "SERPINE1"
#>  [61] "FZD2"     "VAV2"     "NR3C1"    "BCKDHB"   "NR1D2"    "PPARG"   
#>  [67] "CAV1"     "NEK7"     "PDLIM5"   "CDK13"    "ERBB2"    "ITGB3"   
#>  [73] "COL11A1"  "TGFBR2"   "DAB2"     "SRI"      "EFEMP1"   "PIK3R3"  
#>  [79] "FHL2"     "SCN8A"    "INSIG1"   "F3"       "ATP2B4"   "NFKB1"   
#>  [85] "PTPN21"   "NRP1"     "MRPS31"   "PEX14"    "LDLR"     "SNAI2"   
#>  [91] "CDC42BPA" "NOTCH2"   "SIPA1L1"  "TFPI"     "HAS2"     "SDC2"    
#>  [97] "DYRK1A"   "MT1E"     "GRK5"     "BMPR1A"   "GCNT1"    "DLG1"    
#> [103] "ACVR2A"   "LAMC1"    "MAP1B"    "WDR37"    "ID1"      "ADD3"    
#> [109] "AKT3"     "MAPK14"   "RASA2"    "ANXA2"    "CELF2"    "FBLN5"   
#> [115] "MYC"      "MTA1"     "ABCC1"    "KIT"      "KALRN"    "ANXA4"   
#> [121] "PTGFR"    "FYN"      "PIK3CD"   "CACNA1A"  "LTBP1"    "SLC7A1"  
#> [127] "PLCB4"    "PRKCA"    "RND3"     "ARHGEF9"  "ATP2C1"   "ADORA2B" 
#> [133] "MET"      "VLDLR"    "SCHIP1"   "NFIB"     "ADGRL2"   "CCN1"    
#> [139] "DMAC2L"   "PLPP3"    "SCAF8"    "SLC67A1"  "TENT4A"   "TOGARAM1"
#> 
#> $all$HALLMARK_PROTEIN_SECRETION
#>  [1] "KRT18"    "SCAMP1"   "COG2"     "ARFIP1"   "AP2B1"    "MON2"    
#>  [7] "TOM1L1"   "TSG101"   "AP3B1"    "SNX2"     "SNAP23"   "NAPA"    
#> [13] "GALC"     "ICA1"     "RAB14"    "ARFGEF1"  "ARFGEF2"  "GBF1"    
#> [19] "ATP7A"    "PPT1"     "COPB1"    "VAMP4"    "NAPG"     "CLN5"    
#> [25] "TMED10"   "CD63"     "VAMP3"    "TMX1"     "ATP6V1H"  "TPD52"   
#> [31] "CLTC"     "SEC24D"   "CLCN3"    "YIPF6"    "GOSR2"    "VPS4B"   
#> [37] "GLA"      "ATP6V1B1" "ABCA1"    "ADAM10"   "STX16"    "STX12"   
#> [43] "GOLGA4"   "RER1"     "AP1G1"    "RAB22A"   "SGMS1"    "OCRL"    
#> [49] "AP3S1"    "DST"      "RAB5A"    "TMED2"    "ERGIC3"   "SSPN"    
#> [55] "RAB2A"    "SEC31A"   "LMAN1"    "USO1"     "SOD1"     "SCRN1"   
#> [61] "SCAMP3"   "COPB2"    "SEC22B"   "VPS45"    "BET1"     "COPE"    
#> [67] "MAPK1"    "STX7"     "ARCN1"    "M6PR"     "PAM"      "VAMP7"   
#> [73] "BNIP3"    "SH3GL2"   "TSPAN8"   "AP2S1"    "IGF2R"    "ARF1"    
#> [79] "CAV2"     "DNM1L"    "STAM"     "GNAS"     "LAMP2"    "YKT6"    
#> [85] "RAB9A"    "AP2M1"    "ARFGAP3"  "KIF1B"    "CLTA"     "RPS6KA3" 
#> [91] "ATP1A1"   "ZW10"     "CTSC"     "ANP32E"   "EGFR"     "DOP1A"   
#> 
#> $all$HALLMARK_OXIDATIVE_PHOSPHORYLATION
#>   [1] "ACADSB"   "GLUD1"    "ALDH6A1"  "COX6C"    "CYB5A"    "MRPS30"  
#>   [7] "PRDX3"    "SURF1"    "NDUFS4"   "PDHB"     "NDUFA2"   "ETFDH"   
#>  [13] "COX15"    "RHOT2"    "CPT1A"    "ATP6AP1"  "RETSAT"   "UQCRQ"   
#>  [19] "RHOT1"    "ATP6V1G1" "NDUFA5"   "SLC25A20" "ISCU"     "CASP7"   
#>  [25] "PDK4"     "COX7C"    "HSPA9"    "SUCLA2"   "ATP6V1H"  "NDUFA7"  
#>  [31] "COX11"    "COX17"    "BCKDHA"   "GPX4"     "ATP6V1D"  "ATP6V1C1"
#>  [37] "ATP6V0C"  "GOT2"     "ATP6V0E1" "HADHB"    "MAOB"     "UQCR11"  
#>  [43] "SLC25A4"  "NDUFV2"   "ACADVL"   "MRPL35"   "UQCRC2"   "OXA1L"   
#>  [49] "COX6A1"   "ISCA1"    "NDUFC1"   "ACAT1"    "ACAA1"    "MTRR"    
#>  [55] "NNT"      "NDUFB1"   "ECHS1"    "NDUFS7"   "NDUFB4"   "NDUFC2"  
#>  [61] "TIMM9"    "SLC25A12" "CS"       "NDUFA3"   "MFN2"     "TIMM13"  
#>  [67] "ATP6V1E1" "IDH3A"    "TCIRG1"   "ECH1"     "ABCB7"    "DLST"    
#>  [73] "VDAC1"    "PDHX"     "ATP1B1"   "NDUFS1"   "NDUFS8"   "MTRF1"   
#>  [79] "BDH2"     "COX8A"    "ACADM"    "SDHC"     "ETFB"     "GRPEL1"  
#>  [85] "MTX2"     "SLC25A11" "NDUFS2"   "OGDH"     "ATP6V0B"  "NDUFB8"  
#>  [91] "PMPCA"    "SLC25A6"  "NDUFA8"   "NDUFAB1"  "MRPS11"   "IDH3B"   
#>  [97] "MRPL34"   "UQCRC1"   "NDUFB7"   "HADHA"    "NDUFB6"   "COX7A2"  
#> [103] "ACO2"     "PDP1"     "PHYH"     "AIFM1"    "NDUFB5"   "NQO2"    
#> [109] "CYC1"     "COX7A2L"  "NDUFB2"   "TIMM10"   "SLC25A3"  "NDUFS6"  
#> [115] "TIMM17A"  "VDAC3"    "LDHA"     "DECR1"    "OAT"      "VDAC2"   
#> [121] "NDUFA4"   "SDHA"     "IDH1"     "UQCRB"    "NDUFS3"   "UQCR10"  
#> [127] "NDUFA1"   "SDHB"     "MGST3"    "NDUFB3"   "DLD"      "OPA1"    
#> [133] "ETFA"     "ATP6V1F"  "ALAS1"    "LRPPRC"   "COX6B1"   "COX5B"   
#> [139] "IMMT"     "BAX"      "NDUFV1"   "ACAA2"    "HSD17B10" "COX10"   
#> [145] "HTRA2"    "NDUFA6"   "SDHD"     "SUCLG1"   "CYB5R3"   "FXN"     
#> [151] "FH"       "COX7B"    "IDH3G"    "PHB2"     "DLAT"     "MRPS15"  
#> [157] "IDH2"     "MRPS12"   "AFG3L2"   "COX5A"    "COX4I1"   "TIMM8B"  
#> [163] "UQCRFS1"  "HCCS"     "TIMM50"   "MDH2"     "MRPS22"   "FDX1"    
#> [169] "CYCS"     "POR"      "GPI"      "SUPV3L1"  "SLC25A5"  "MRPL11"  
#> [175] "MRPL15"   "TOMM22"   "NDUFA9"   "POLR2F"   "MDH1"     "LDHB"    
#> [181] "PDHA1"    "UQCRH"    "ATP5F1A"  "ATP5F1B"  "ATP5F1C"  "ATP5F1D" 
#> [187] "ATP5F1E"  "ATP5MC1"  "ATP5MC2"  "ATP5MC3"  "ATP5ME"   "ATP5MF"  
#> [193] "ATP5MG"   "ATP5PB"   "ATP5PD"   "ATP5PF"   "ATP5PO"   "ECI1"    
#> [199] "MPC1"     "TOMM70"  
#> 
#> $all$HALLMARK_DNA_REPAIR
#>   [1] "DCTN4"   "BCAM"    "GMPR2"   "ADCY6"   "XPC"     "NME3"    "POLD4"  
#>   [8] "SMAD5"   "TAF9"    "TSG101"  "SURF1"   "AAAS"    "POLL"    "CANT1"  
#>  [15] "VPS37D"  "RRM2B"   "TK2"     "DDB2"    "IMPDH2"  "CCNO"    "POLR2A" 
#>  [22] "GTF2F1"  "ERCC4"   "MPG"     "NUDT9"   "COX17"   "TARBP2"  "SUPT4H1"
#>  [29] "GPX4"    "GTF2H3"  "VPS28"   "NME4"    "ERCC8"   "CETN2"   "NT5C"   
#>  [36] "DAD1"    "POLR2E"  "EIF1B"   "AK1"     "ARL6IP1" "EDF1"    "TP53"   
#>  [43] "POLR3GL" "BRF2"    "POLB"    "GTF2H5"  "RPA3"    "RAD52"   "ERCC2"  
#>  [50] "PRIM1"   "SNAPC5"  "MRPL40"  "ERCC1"   "TAF10"   "USP11"   "RAE1"   
#>  [57] "DUT"     "TMED2"   "NFX1"    "ERCC5"   "NPR2"    "DDB1"    "GTF2H1" 
#>  [64] "SNAPC4"  "ZNF707"  "CDA"     "POLR2G"  "NCBP2"   "AK3"     "ELL"    
#>  [71] "RFC5"    "TAF12"   "POLA1"   "LIG1"    "REV3L"   "NME1"    "POLR1D" 
#>  [78] "POLR2K"  "POLE4"   "RNMT"    "BOLA2"   "ERCC3"   "POLR2I"  "PDE6G"  
#>  [85] "STX3"    "BCAP31"  "HCLS1"   "ITPA"    "POLR2J"  "RALA"    "RPA2"   
#>  [92] "APRT"    "POLH"    "GTF2B"   "CMPK2"   "DGCR8"   "CSTF3"   "GTF2A2" 
#>  [99] "UMPS"    "POLR2H"  "GUK1"    "NUDT21"  "POLR2C"  "TAF13"   "ZWINT"  
#> [106] "POLR3C"  "POM121"  "RBX1"    "ADRM1"   "SAC3D1"  "SDCBP"   "TAF6"   
#> [113] "TAF1C"   "HPRT1"   "PDE4B"   "POLA2"   "POLD3"   "PCNA"    "CLP1"   
#> [120] "SEC61A1" "RFC3"    "RAD51"   "ADA"     "POLD1"   "GTF3C5"  "SUPT5H" 
#> [127] "DGUOK"   "FEN1"    "TYMS"    "PNP"     "POLR1C"  "SSRP1"   "POLR2D" 
#> [134] "VPS37B"  "SF3A3"   "UPF3B"   "POLR2F"  "RFC2"    "RFC4"    "AGO4"   
#> [141] "ALYREF"  "ELOA"    "GSDME"   "MPC2"    "NELFB"   "NELFCD"  "NELFE"  
#> [148] "NT5C3A"  "POLR1H"  "SRSF6"  
#> 
#> $all$HALLMARK_KRAS_SIGNALING_DN
#>   [1] "GAMT"      "SIDT1"     "GP2"       "BTG2"      "RGS11"     "CPEB3"    
#>   [7] "CCDC106"   "CACNG1"    "SLC29A3"   "LFNG"      "CPB1"      "BMPR1B"   
#>  [13] "CAPN9"     "IGFBP2"    "GPRC5C"    "CACNA1F"   "THRB"      "CELSR2"   
#>  [19] "PDK2"      "FGFR3"     "IDUA"      "C5"        "MYO15A"    "PDE6B"    
#>  [25] "BRDT"      "ZBTB16"    "NR4A2"     "EFHD1"     "HNF1A"     "MTHFR"    
#>  [31] "NRIP2"     "ADRA2C"    "SPHK2"     "TAS2R4"    "FGF16"     "P2RY4"    
#>  [37] "COPZ2"     "MAST3"     "EGF"       "EPHA5"     "TNNI3"     "ARPP21"   
#>  [43] "SERPINA10" "ITIH3"     "LYPD3"     "THNSL2"    "TFAP2B"    "SLC25A23" 
#>  [49] "KRT1"      "ASB7"      "ATP6V1B1"  "FGGY"      "SKIL"      "MYOT"     
#>  [55] "HTR1D"     "P2RX6"     "PTPRJ"     "MFSD6"     "GDNF"      "PCDHB1"   
#>  [61] "PNMT"      "MAGIX"     "ABCB11"    "PRODH"     "NTF3"      "NR6A1"    
#>  [67] "SYNPO"     "KRT13"     "OXT"       "HSD11B2"   "COL2A1"    "HTR1B"    
#>  [73] "ENTPD7"    "ARHGDIG"   "RIBC2"     "NPHS1"     "KCNQ2"     "SLC6A3"   
#>  [79] "BARD1"     "SLC5A5"    "ACTC1"     "TG"        "PAX3"      "SPRR3"    
#>  [85] "NGB"       "CLSTN3"    "YBX2"      "SCGB1A1"   "CD207"     "SLC16A7"  
#>  [91] "DLK2"      "AKR1B10"   "ALOX12B"   "KCND1"     "CNTFR"     "GP1BA"    
#>  [97] "IRS4"      "NOS1"      "SHOX2"     "MX1"       "SNCB"      "MSH5"     
#> [103] "UPK3B"     "IFI44L"    "RSAD2"     "ABCG4"     "WNT16"     "YPEL1"    
#> [109] "RYR2"      "CD40LG"    "IL12B"     "PLAG1"     "KCNE2"     "MYH7"     
#> [115] "SNN"       "MEFV"      "TGFB2"     "KLHDC8A"   "SLC12A3"   "CCR8"     
#> [121] "LGALS7"    "NUDT11"    "SGK1"      "CD80"      "EDN2"      "TFF2"     
#> [127] "KCNN1"     "KRT15"     "ITGB1BP2"  "KCNMB1"    "CAMK1D"    "CLDN8"    
#> [133] "DCC"       "SERPINB2"  "TCL1A"     "CLPS"      "TFCP2L1"   "PDCD1"    
#> [139] "PTGFR"     "KRT5"      "CKM"       "CCNA1"     "SPTBN2"    "EDN1"     
#> [145] "IFNG"      "CALML5"    "SLC38A3"   "EDAR"      "RYR1"      "KRT4"     
#> [151] "SLC30A3"   "GTF3C5"    "PKP1"      "SOX10"     "STAG3"     "KLK7"     
#> [157] "CALCB"     "CYP39A1"   "CPA2"      "SMPX"      "GPR19"     "TEX15"    
#> [163] "DTNB"      "KLK8"      "CDKAL1"    "TGM1"      "TLX1"      "GPR3"     
#> [169] "CHST2"     "CLDN16"    "TCF7L1"    "SLC6A14"   "AMBN"      "ATP4A"    
#> [175] "CDH16"     "CHRNG"     "COQ8A"     "CYP11B2"   "FGF22"     "FSHB"     
#> [181] "GRID2"     "IL5"       "INSL5"     "KMT2D"     "MACROH2A2" "NPY4R"    
#> [187] "PAX4"      "PRKN"      "PROP1"     "SCN10A"    "SELENOP"   "SSTR4"    
#> [193] "TENM2"     "TENT5C"    "TSHB"      "UGT2B17"   "VPREB1"    "VPS50"    
#> [199] "ZC2HC1C"   "ZNF112"   
#> 
#> $all$HALLMARK_PEROXISOME
#>   [1] "ABCC8"    "HSD17B4"  "ELOVL5"   "SEMA3C"   "HMGCL"    "SLC27A2" 
#>   [7] "PEX11A"   "ABCD3"    "SERPINA6" "SULT2B1"  "DIO1"     "DLG4"    
#>  [13] "DHCR24"   "CDK7"     "CRAT"     "MVP"      "ABCC5"    "HSD3B7"  
#>  [19] "RETSAT"   "LONP2"    "CRABP2"   "EPHX2"    "ABCB4"    "ISOC1"   
#>  [25] "ATXN1"    "SCP2"     "FIS1"     "PEX2"     "CAT"      "RXRG"    
#>  [31] "IDE"      "STS"      "SLC25A4"  "ABCB1"    "ALB"      "VPS4B"   
#>  [37] "PEX11B"   "ACAA1"    "HSD17B11" "ALDH9A1"  "RDH11"    "ALDH1A1" 
#>  [43] "TTR"      "HAO2"     "MLYCD"    "ABCB9"    "ERCC1"    "ECH1"    
#>  [49] "CLN8"     "ACOT8"    "HRAS"     "SLC25A17" "HSD11B2"  "CTBP1"   
#>  [55] "SOD1"     "PEX14"    "EHHADH"   "GSTK1"    "PEX6"     "ACOX1"   
#>  [61] "GNPAT"    "SCGB1A1"  "PEX5"     "BCL10"    "CADM1"    "SIAH1"   
#>  [67] "SLC23A2"  "CEL"      "ERCC3"    "SMARCC1"  "IDH1"     "PRDX5"   
#>  [73] "DHRS3"    "FADS1"    "ABCD2"    "CNBP"     "TOP2A"    "ACSL1"   
#>  [79] "FABP6"    "ACSL5"    "PRDX1"    "CLN6"     "CACNA1B"  "IDI1"    
#>  [85] "NR1I2"    "NUDT19"   "YWHAH"    "ESR2"     "ITGB1BP1" "ACSL4"   
#>  [91] "IDH2"     "PEX13"    "SLC35B2"  "ABCD1"    "PABPC1"   "SLC25A19"
#>  [97] "FDPS"     "TSPO"     "CRABP1"   "SOD2"     "MSH2"     "CTPS1"   
#> [103] "ECI2"     "UGT2B17" 
#> 
#> $all$HALLMARK_APICAL_JUNCTION
#>   [1] "GAMT"      "EVL"       "WNK4"      "IRS1"      "LIMA1"     "TAOK2"    
#>   [7] "RASA1"     "SHROOM2"   "PCDH1"     "CADM2"     "PTEN"      "AMIGO2"   
#>  [13] "CTNNA1"    "ITGA3"     "NEGR1"     "CRAT"      "NF1"       "JAM3"     
#>  [19] "BAIAP2"    "LAMA3"     "GTF2F1"    "WASL"      "TMEM8B"    "KRT31"    
#>  [25] "TSC1"      "CERCAM"    "CDH11"     "ITGA2"     "FBN1"      "CLDN5"    
#>  [31] "CDH6"      "SLIT2"     "AMIGO1"    "PKD1"      "CRB3"      "VWF"      
#>  [37] "CD34"      "ADAM23"    "CLDN19"    "SORBS3"    "MYH10"     "TRO"      
#>  [43] "GNAI2"     "JUP"       "ARHGEF6"   "NRXN2"     "INPPL1"    "TJP1"     
#>  [49] "LAYN"      "VCAN"      "TSPAN4"    "PPP2R2C"   "THY1"      "CDSN"     
#>  [55] "ALOX15B"   "NLGN3"     "EXOC4"     "VAV2"      "PARVA"     "COL16A1"  
#>  [61] "NLGN2"     "MMP2"      "STX4"      "COL17A1"   "DMP1"      "AKT2"     
#>  [67] "CTNND1"    "CNTN1"     "FLNC"      "CDH8"      "CLDN11"    "SHC1"     
#>  [73] "CLDN7"     "SYMPK"     "AMH"       "PIK3R3"    "ATP1A3"    "PECAM1"   
#>  [79] "LDLRAP1"   "INSIG1"    "CLDN18"    "CDH1"      "MADCAM1"   "HRAS"     
#>  [85] "TIAL1"     "HADH"      "ITGA10"    "TNFRSF11B" "CADM3"     "NRAP"     
#>  [91] "RRAS"      "PIK3CB"    "MYL9"      "TUBG1"     "ADRA1B"    "MAPK11"   
#>  [97] "ACTC1"     "PARD6G"    "B4GALT1"   "ACTA1"     "MAPK13"    "NEXN"     
#> [103] "ACTG1"     "PTK2"      "CD99"      "ADAM9"     "SKAP2"     "CLDN15"   
#> [109] "NFASC"     "ACTN3"     "VCL"       "DLG1"      "KCNH2"     "CDH15"    
#> [115] "DSC1"      "ADAMTS5"   "GRB7"      "ITGA9"     "CD276"     "THBS3"    
#> [121] "RAC2"      "IKBKG"     "PTPRC"     "SGCE"      "CD209"     "AKT3"     
#> [127] "BMP1"      "ACTN2"     "MAPK14"    "MYL12B"    "CD86"      "TGFBI"    
#> [133] "PLCG1"     "MMP9"      "ITGB1"     "CNN2"      "SDC3"      "PBX2"     
#> [139] "CAP1"      "MDK"       "CD274"     "MYH9"      "ACTG2"     "SRC"      
#> [145] "SYK"       "DHX16"     "MAP4K2"    "LAMB3"     "VASP"      "EPB41L2"  
#> [151] "ADAM15"    "VCAM1"     "CLDN8"     "CDH4"      "ICAM2"     "ACTN1"    
#> [157] "RHOF"      "LAMC2"     "CLDN14"    "ITGB4"     "TRAF1"     "YWHAH"    
#> [163] "ICAM5"     "CLDN4"     "NF2"       "ZYX"       "ARPC2"     "DSC3"     
#> [169] "CLDN9"     "ICAM4"     "SPEG"      "ACTN4"     "ICAM1"     "SLC30A3"  
#> [175] "GNAI1"     "MVD"       "ACTB"      "CALB2"     "PFN1"      "COL9A1"   
#> [181] "CDK8"      "SIRPA"     "NRTN"      "MPZL1"     "CX3CL1"    "CLDN6"    
#> [187] "MPZL2"     "MSN"       "EGFR"      "CDH3"      "FSCN1"     "RSU1"     
#> [193] "FYB1"      "MAP3K20"   "NECTIN1"   "NECTIN2"   "NECTIN3"   "NECTIN4"  
#> [199] "NHERF4"    "PALS1"    
#> 
#> $all$HALLMARK_APICAL_SURFACE
#>  [1] "GATA3"    "ATP8B1"   "RTN4RL1"  "GSTM3"    "SCUBE1"   "SHROOM2" 
#>  [7] "HSPB1"    "SLC34A3"  "MDGA1"    "ADIPOR2"  "TMEM8B"   "AFAP1L2" 
#> [13] "SULF2"    "BRCA1"    "LYPD3"    "CROCC"    "THY1"     "ADAM10"  
#> [19] "NTNG1"    "SLC2A4"   "NCOA6"    "SRPX"     "CD160"    "FLOT2"   
#> [25] "B4GALT1"  "PKHD1"    "GAS1"     "AKAP7"    "IL2RB"    "ATP6V0A4"
#> [31] "PCSK9"    "DCBLD2"   "IL2RG"    "GHRL"     "EPHB4"    "MAL"     
#> [37] "APP"      "PLAUR"    "EFNA5"    "CX3CL1"   "LYN"      "RHCG"    
#> [43] "CRYBG1"   "SLC22A12"
#> 
#> $all$HALLMARK_HEME_METABOLISM
#>   [1] "BCAM"     "BTG2"     "BTRC"     "TMEM9B"   "ALAD"     "PIGQ"    
#>   [7] "ADD1"     "VEZF1"    "EZH1"     "ALDH6A1"  "NUDT4"    "LRP10"   
#>  [13] "SELENBP1" "DCAF10"   "RBM5"     "TRIM58"   "BLVRA"    "EPOR"    
#>  [19] "SLC22A4"  "RNF123"   "HEBP1"    "SLC25A38" "HAGH"     "ABCG2"   
#>  [25] "ARHGEF12" "GDE1"     "SIDT2"    "DCAF11"   "GAPVD1"   "KHNYN"   
#>  [31] "MINPP1"   "GATA1"    "ATP6V0A1" "BLVRB"    "TRAK2"    "SLC30A1" 
#>  [37] "BNIP3L"   "NCOA4"    "FN3K"     "TNS1"     "CLCN3"    "DAAM1"   
#>  [43] "FECH"     "YPEL5"    "CAST"     "TSPAN5"   "HBB"      "ADIPOR1" 
#>  [49] "ACP5"     "KLF3"     "TAL1"     "CAT"      "NFE2L1"   "ELL2"    
#>  [55] "PPOX"     "ERMAP"    "SYNJ1"    "PRDX2"    "ISCA1"    "CCDC28A" 
#>  [61] "CDC27"    "FBXO9"    "RAP1GAP"  "NR3C1"    "FBXO34"   "NNT"     
#>  [67] "CA2"      "KAT2B"    "NFE2"     "NEK7"     "EPB42"    "BMP2K"   
#>  [73] "GYPE"     "TOP1"     "LMO2"     "GCLC"     "PSMD9"    "MKRN1"   
#>  [79] "SNCA"     "CTNS"     "LPIN2"    "ALAS2"    "UROS"     "EPB41"   
#>  [85] "HBD"      "UCP2"     "MYL4"     "TCEA1"    "CIR1"     "EIF2AK1" 
#>  [91] "SLC4A1"   "SLC11A2"  "ATG4A"    "ANK1"     "GLRX5"    "MXI1"    
#>  [97] "UROD"     "OSBP2"    "SPTA1"    "ARL2BP"   "PGLS"     "FOXJ2"   
#> [103] "MBOAT2"   "MARK3"    "CDR2"     "RHCE"     "USP15"    "BSG"     
#> [109] "RHD"      "FOXO3"    "P4HA2"    "AQP3"     "RANBP10"  "ALDH1L1" 
#> [115] "MPP1"     "CLIC2"    "ABCB6"    "PPP2R5B"  "PC"       "RNF19A"  
#> [121] "ENDOD1"   "BPGM"     "FBXO7"    "IGSF3"    "KEL"      "CTSB"    
#> [127] "C3"       "MGST3"    "TNRC6B"   "SPTB"     "RAD23A"   "GCLM"    
#> [133] "GYPC"     "XPO7"     "SLC7A11"  "CTSE"     "UBAC1"    "MOCOS"   
#> [139] "BACH1"    "CCND3"    "SEC14L1"  "SDCBP"    "RCL1"     "TFRC"    
#> [145] "DCUN1D1"  "ACSL6"    "HTRA2"    "MOSPD1"   "PICALM"   "KLF1"    
#> [151] "RIOK3"    "SLC6A9"   "SLC6A8"   "LAMP2"    "HTATIP2"  "TFDP2"   
#> [157] "CPOX"     "SLC2A1"   "HBQ1"     "ICAM4"    "PDZK1IP1" "MFHAS1"  
#> [163] "XK"       "HMBS"     "OPTN"     "AGPAT4"   "MAP2K3"   "ADD2"    
#> [169] "FTCD"     "TRIM10"   "E2F2"     "SLC30A10" "SLC10A3"  "NARF"    
#> [175] "SMOX"     "HDGF"     "ASNS"     "SLC25A37" "TMCC2"    "RBM38"   
#> [181] "GMPS"     "ACKR1"    "AHSP"     "CA1"      "CROCCP2"  "DMTN"    
#> [187] "GYPA"     "GYPB"     "H1-0"     "H4C3"     "HBBP1"    "HBZ"     
#> [193] "KDM7A"    "MARCHF2"  "MARCHF8"  "RHAG"     "SLC66A2"  "TENT5C"  
#> [199] "TSPO2"    "TYR"     
#> 
#> $all$HALLMARK_ANGIOGENESIS
#>  [1] "SERPINA5" "THBD"     "POSTN"    "COL3A1"   "MSX1"     "LUM"     
#>  [7] "FGFR1"    "COL5A2"   "LRPAP1"   "STC1"     "VTN"      "VCAN"    
#> [13] "VAV2"     "ITGAV"    "JAG1"     "SLCO2A1"  "KCNJ8"    "OLR1"    
#> [19] "TIMP1"    "FSTL1"    "NRP1"     "CCND2"    "PDGFA"    "LPL"     
#> [25] "PTK2"     "JAG2"     "SPP1"     "S100A4"   "VEGFA"    "APP"     
#> [31] "CXCL6"    "PRG2"     "TNFRSF21" "APOH"     "PF4"      "PGLYRP1" 
#> 
#> $all$HALLMARK_XENOBIOTIC_METABOLISM
#>   [1] "ESR1"      "FBP1"      "MCCC2"     "TMBIM6"    "IGFBP4"    "CSAD"     
#>   [7] "ELOVL5"    "MPP2"      "ACOX2"     "UGDH"      "NINJ1"     "CROT"     
#>  [13] "CYB5A"     "FAH"       "SERPINA6"  "ENTPD5"    "ETFDH"     "SAR1B"    
#>  [19] "ARPP19"    "CFB"       "SLC35D1"   "ACSM1"     "BLVRB"     "CYFIP2"   
#>  [25] "RETSAT"    "GNMT"      "DHPS"      "MAOA"      "CYP26A1"   "ACOX3"    
#>  [31] "PTGR1"     "DCXR"      "NQO1"      "PDK4"      "GSTM4"     "DHRS7"    
#>  [37] "MAN1A1"    "HACL1"     "CYP2E1"    "NFS1"      "PINK1"     "TAT"      
#>  [43] "IGF1"      "ASL"       "SLC46A3"   "AOX1"      "GSS"       "LEAP2"    
#>  [49] "PTGES"     "CYP4F2"    "CAT"       "PTGES3"    "NMT1"      "F10"      
#>  [55] "EPHX1"     "CD36"      "VTN"       "DHRS1"     "GAD1"      "CES1"     
#>  [61] "JUP"       "GSR"       "CDO1"      "GSTT2"     "SLC6A6"    "SERPINE1" 
#>  [67] "RAP1GAP"   "ADH1C"     "CA2"       "ABHD6"     "ENPEP"     "SERTAD1"  
#>  [73] "PDLIM5"    "GSTA3"     "CYP17A1"   "ALDH9A1"   "DDAH2"     "IGFBP1"   
#>  [79] "CYP1A1"    "ITIH4"     "ATOH8"     "GCLC"      "ALDH2"     "ACP2"     
#>  [85] "CASP6"     "CNDP2"     "ABCC3"     "LPIN2"     "FBLN1"     "ECH1"     
#>  [91] "BPHL"      "RBP4"      "HNF4A"     "CYP2S1"    "SPINT2"    "SLC12A4"  
#>  [97] "MBL2"      "IL1R1"     "SLC6A12"   "PEMT"      "TMEM97"    "CDA"      
#> [103] "ARG1"      "CYP2J2"    "HES6"      "PROS1"     "HMOX1"     "ACOX1"    
#> [109] "BCAR1"     "TPST1"     "HGFAC"     "ATP2A2"    "APOE"      "PMM1"     
#> [115] "GSTO1"     "GCKR"      "ITIH1"     "SLC35B1"   "PAPSS2"    "ACO2"     
#> [121] "PTGDS"     "AKR1C3"    "UPB1"      "SSR3"      "TTPA"      "TNFRSF1A" 
#> [127] "MTHFD1"    "PC"        "IRF8"      "ANGPTL3"   "FMO3"      "TMEM176B" 
#> [133] "LONP1"     "PSMB10"    "ALDH3A1"   "IDH1"      "FMO1"      "HSD11B1"  
#> [139] "AKR1C2"    "SLC1A5"    "CRP"       "ABCD2"     "GABARAPL1" "PYCR1"    
#> [145] "TGFB2"     "MT2A"      "GCH1"      "LCAT"      "AHCY"      "ID2"      
#> [151] "ADH5"      "ACP1"      "ALAS1"     "COMT"      "CYP27A1"   "FAS"      
#> [157] "SLC22A1"   "BCAT1"     "HPRT1"     "ABCC2"     "CYP2C18"   "ARG2"     
#> [163] "AP4B1"     "NDRG2"     "DDT"       "CCL25"     "ETS2"      "CBR1"     
#> [169] "KYNU"      "TDO2"      "PGD"       "XDH"       "PTS"       "VNN1"     
#> [175] "AQP9"      "NPC1"      "PGRMC1"    "POR"       "SHMT2"     "GART"     
#> [181] "EPHA2"     "HSD17B2"   "DDC"       "PPARD"     "GCNT2"     "SMOX"     
#> [187] "FETUB"     "UPP1"      "ADH7"      "CYP1A2"    "F11"       "FABP1"    
#> [193] "G6PC1"     "HRG"       "KARS1"     "MARCHF6"   "PLG"       "REG1A"    
#> [199] "TKFC"      "TYR"      
#> 
#> $all$HALLMARK_FATTY_ACID_METABOLISM
#>   [1] "SLC22A5"   "HSD17B4"   "REEP6"     "ELOVL5"    "GLUL"      "ALAD"     
#>   [7] "HMGCL"     "UGDH"      "HSDL2"     "BMPR1B"    "BLVRA"     "SERINC1"  
#>  [13] "SUCLG2"    "HMGCS2"    "PDHB"      "ACOT2"     "ETFDH"     "LTC4S"    
#>  [19] "DHCR24"    "ALDH3A2"   "CRAT"      "PSME1"     "CPT1A"     "RETSAT"   
#>  [25] "CYP4A11"   "MAOA"      "ADIPOR2"   "PTPRG"     "FASN"      "AUH"      
#>  [31] "D2HGDH"    "GSTZ1"     "ALDOA"     "HSD17B7"   "RDH16"     "ACADS"    
#>  [37] "SUCLA2"    "AOC3"      "ENO2"      "INMT"      "RAP1GDS1"  "AQP7"     
#>  [43] "ACSS1"     "GPD2"      "CRYZ"      "MGLL"      "XIST"      "CYP4A22"  
#>  [49] "CPT2"      "EPHX1"     "CD36"      "GPD1"      "HADHB"     "CA4"      
#>  [55] "ACADVL"    "CIDEA"     "MCEE"      "NTHL1"     "ENO3"      "ADH1C"    
#>  [61] "ACAA1"     "HPGD"      "BCKDHB"    "CA2"       "HSD17B11"  "PCBD1"    
#>  [67] "ECHS1"     "ALDH9A1"   "RDH11"     "ALDH1A1"   "GCDH"      "CYP1A1"   
#>  [73] "UBE2L6"    "NBN"       "HAO2"      "MLYCD"     "UROS"      "ECH1"     
#>  [79] "GRHPR"     "BPHL"      "DLST"      "ACOT8"     "TP53INP2"  "ACADM"    
#>  [85] "SDHC"      "HADH"      "HIBCH"     "UROD"      "ERP29"     "EHHADH"   
#>  [91] "APEX1"     "ACOX1"     "IDH3B"     "G0S2"      "LGALS1"    "ACO2"     
#>  [97] "HSP90AA1"  "CD1D"      "HSPH1"     "CEL"       "MIF"       "ALDH3A1"  
#> [103] "LDHA"      "DECR1"     "SDHA"      "IDH1"      "HMGCS1"    "FMO1"     
#> [109] "ACADL"     "ACSM3"     "GABARAPL1" "DLD"       "SMS"       "ACSL1"    
#> [115] "ACSL5"     "PRDX6"     "GAPDHS"    "ACAA2"     "HSD17B10"  "IDI1"     
#> [121] "ACAT2"     "OSTC"      "SDHD"      "SUCLG1"    "METAP1"    "FH"       
#> [127] "IDH3G"     "YWHAH"     "CPOX"      "NCAPH2"    "CBR1"      "NSDHL"    
#> [133] "ME1"       "ACSL4"     "TDO2"      "HCCS"      "PTS"       "MDH2"     
#> [139] "IL4I1"     "VNN1"      "CBR3"      "S100A10"   "AADAT"     "ADSL"     
#> [145] "MDH1"      "CA6"       "PPARA"     "ODC1"      "PDHA1"     "ADH7"     
#> [151] "ECI1"      "ECI2"      "FABP1"     "FABP2"     "GAD2"      "H2AZ1"    
#> [157] "KMT5A"     "MIX23"    
#> 
#> $all$HALLMARK_P53_PATHWAY
#>   [1] "ABAT"     "ANKRA2"   "WWP1"     "GLS2"     "HEXIM1"   "BTG2"    
#>   [7] "SLC19A2"  "CGRRF1"   "PRKAB1"   "ACVR1B"   "XPC"      "NINJ1"   
#>  [13] "KIF13B"   "SERTAD3"  "FUCA1"    "SP1"      "TOB1"     "RB1"     
#>  [19] "MDM2"     "CTSF"     "RAB40C"   "RXRA"     "RPS27L"   "SESN1"   
#>  [25] "FDXR"     "SLC35D1"  "CSRNP2"   "ZBTB16"   "ABCC5"    "CYFIP2"  
#>  [31] "DDB2"     "RETSAT"   "MKNK2"    "BAIAP2"   "ABHD4"    "CCNG1"   
#>  [37] "ZNF365"   "MXD4"     "DCXR"     "IER3"     "CDKN2AIP" "ISCU"    
#>  [43] "F2R"      "TRAFD1"   "STOM"     "RRP8"     "TXNIP"    "NHLH2"   
#>  [49] "FOS"      "PLK2"     "BLCAP"    "TGFB1"    "VWA5A"    "IP6K2"   
#>  [55] "TSPYL2"   "PVT1"     "PHLDA3"   "EPHX1"    "GPX2"     "CDKN1A"  
#>  [61] "PTPRE"    "TRIB3"    "PPM1D"    "TSC22D1"  "MAPKAPK3" "KLF4"    
#>  [67] "PROCR"    "RCHY1"    "ALOX15B"  "AK1"      "TP63"     "CDH13"   
#>  [73] "TP53"     "RAD51C"   "CEBPA"    "NOL8"     "INHBB"    "HDAC3"   
#>  [79] "FAM162A"  "APAF1"    "CTSD"     "DGKA"     "JUN"      "HRAS"    
#>  [85] "ERCC5"    "HINT1"    "ZMAT3"    "RGS16"    "CCND2"    "OSGIN1"  
#>  [91] "RAD9A"    "NUPR1"    "RPL18"    "PDGFA"    "RRAD"     "PLXNB2"  
#>  [97] "PLK3"     "EPS8L2"   "VAMP8"    "ZFP36L1"  "MXD1"     "RALGDS"  
#> [103] "HMOX1"    "GM2A"     "TCHH"     "FOXO3"    "PMM1"     "GADD45A" 
#> [109] "BTG1"     "JAG2"     "DEF6"     "CD81"     "DDIT3"    "RPL36"   
#> [115] "CLCA2"    "VDR"      "CDKN2B"   "CDK5R1"   "TM7SF3"   "DRAM1"   
#> [121] "TRIAP1"   "HSPA4L"   "IER5"     "FGF13"    "TRAF4"    "CCNK"    
#> [127] "POLH"     "ATF3"     "PPP1R15A" "TPD52L1"  "SLC7A11"  "HBEGF"   
#> [133] "PITPNC1"  "SLC3A2"   "TNNI1"    "SFN"      "SDC1"     "CCND3"   
#> [139] "CASP1"    "BAX"      "AEN"      "FAS"      "TPRKB"    "POM121"  
#> [145] "EI24"     "LIF"      "TNFSF9"   "S100A4"   "IFI30"    "RPS12"   
#> [151] "ITGB4"    "FBXW7"    "RNF19B"   "PCNA"     "SAT1"     "TAX1BP3" 
#> [157] "TCN2"     "PTPN14"   "TAP1"     "KRT17"    "IL1A"     "BMP2"    
#> [163] "PRMT2"    "NUDT15"   "SEC61A1"  "DDIT4"    "SPHK1"    "ADA"     
#> [169] "APP"      "DNTTIP2"  "SOCS1"    "TM4SF1"   "SERPINB5" "S100A10" 
#> [175] "CD82"     "NDRG1"    "NOTCH1"   "BAK1"     "EPHA2"    "TGFA"    
#> [181] "PERP"     "RHBDF2"   "STEAP3"   "KLK8"     "CDKN2A"   "ST14"    
#> [187] "IRAK1"    "UPP1"     "RAP2B"    "LDHB"     "CCP110"   "COQ8A"   
#> [193] "ELP1"     "H1-2"     "H2AC25"   "H2AJ"     "IRAG2"    "PIDD1"   
#> [199] "RACK1"    "WRAP73"  
#> 
#> $all$HALLMARK_HEDGEHOG_SIGNALING
#>  [1] "CELSR1" "TLE3"   "RASA1"  "NF1"    "UNC5C"  "RTN1"   "GLI1"   "NRCAM" 
#>  [9] "SCG2"   "LDB1"   "SLIT1"  "THY1"   "HEY1"   "CRMP1"  "NKX6-1" "NRP1"  
#> [17] "HEY2"   "SHH"    "OPHN1"  "CNTFR"  "NRP2"   "CDK5R1" "AMOT"   "MYH9"  
#> [25] "ACHE"   "DPYSL2" "ETS2"   "VEGFA"  "L1CAM"  "PTCH1"  "PML"    "TLE1"  
#> [33] "VLDLR"  "CDK6"   "ADGRG1" "PLG"   
#> 
#> $all$HALLMARK_ANDROGEN_RESPONSE
#>   [1] "SPDEF"    "INPP4B"   "GPD1L"    "ELOVL5"   "CCND1"    "STK39"   
#>   [7] "ABHD2"    "BMPR1B"   "SORD"     "APPBP2"   "KRT8"     "DHCR24"  
#>  [13] "IQGAP2"   "ZMIZ1"    "KRT19"    "AKAP12"   "AZGP1"    "NCOA4"   
#>  [19] "LIFR"     "MAK"      "HMGCR"    "TPD52"    "SEC24D"   "NKX3-1"  
#>  [25] "HERC3"    "STEAP4"   "ELL2"     "SCD"      "GSR"      "TSC22D1" 
#>  [31] "SLC38A2"  "VAPA"     "KLK3"     "CAMKK2"   "ITGAV"    "SPCS3"   
#>  [37] "DNAJB9"   "HPGD"     "MAP7"     "ACSL3"    "HSD17B14" "PDLIM5"  
#>  [43] "TMEM50A"  "AKT1"     "KLK2"     "SRP19"    "INSIG1"   "NGLY1"   
#>  [49] "PTK2B"    "PTPN21"   "HOMER2"   "PIAS1"    "XRCC5"    "TARP"    
#>  [55] "LMAN1"    "MAF"      "PA2G4"    "ARID5B"   "B4GALT1"  "ZBTB10"  
#>  [61] "MERTK"    "PMEPA1"   "UBE2I"    "B2M"      "ELK4"     "FKBP5"   
#>  [67] "RRP12"    "PGM3"     "HMGCS1"   "RAB4A"    "FADS1"    "TMPRSS2" 
#>  [73] "SLC26A2"  "SGK1"     "ADAMTS1"  "ANKH"     "SMS"      "TNFAIP8" 
#>  [79] "CCND3"    "CDC14B"   "ADRM1"    "IDI1"     "ACTN1"    "UAP1"    
#>  [85] "GNAI3"    "SAT1"     "DBI"      "ABCC4"    "MYL12A"   "UBE2J1"  
#>  [91] "XRCC6"    "RPS6KA3"  "ALDH1A3"  "NDRG1"    "SRF"      "CDK6"    
#>  [97] "CENPN"    "GUCY1A1"  "H1-0"     "PLPP1"    "SELENOP" 
#> 
#> $all$HALLMARK_APOPTOSIS
#>   [1] "RHOB"      "RARA"      "ERBB3"     "KRT18"     "BTG2"      "CCND1"    
#>   [7] "ADD1"      "MADD"      "RNASEL"    "PLAT"      "BCL2L1"    "TIMP3"    
#>  [13] "HSPB1"     "EGR3"      "FDXR"      "PSEN1"     "RHOT2"     "CREBBP"   
#>  [19] "HGF"       "RETSAT"    "PDCD4"     "BNIP3L"    "XIAP"      "LEF1"     
#>  [25] "PDGFRB"    "BIK"       "IER3"      "PPT1"      "CASP7"     "BCL2L11"  
#>  [31] "F2R"       "PSEN2"     "PMAIP1"    "MGMT"      "BCL2L2"    "CDKN1B"   
#>  [37] "TXNIP"     "ENO2"      "CASP9"     "IGFBP6"    "DCN"       "SLC20A1"  
#>  [43] "GPX4"      "LUM"       "SQSTM1"    "TGFBR3"    "GUCY2D"    "SPTAN1"   
#>  [49] "BRCA1"     "CDKN1A"    "SMAD7"     "FEZ1"      "GSR"       "AIFM3"    
#>  [55] "TNFSF10"   "CLU"       "MMP2"      "CAV1"      "GSTM1"     "GADD45B"  
#>  [61] "ERBB2"     "AVPR1A"    "ROCK1"     "DPYD"      "ETF1"      "DIABLO"   
#>  [67] "CTH"       "CASP6"     "TIMP2"     "WEE1"      "TIMP1"     "NEFH"     
#>  [73] "JUN"       "BGN"       "GPX1"      "CASP8"     "DNAJC3"    "CCND2"    
#>  [79] "IFITM3"    "NEDD9"     "SOD1"      "GPX3"      "PAK1"      "CD69"     
#>  [85] "HMOX1"     "CD44"      "CYLD"      "CTNNB1"    "PTK2"      "GADD45A"  
#>  [91] "PEA15"     "LGALS3"    "BCL10"     "BMF"       "DDIT3"     "PPP2R5B"  
#>  [97] "EMP1"      "IRF1"      "DAP3"      "PLCB2"     "EREG"      "VDAC2"    
#> [103] "RELA"      "GSN"       "BCAP31"    "IL1B"      "IGF2R"     "ATF3"     
#> [109] "MCL1"      "TGFB2"     "IL18"      "GCH1"      "GNA15"     "CDK2"     
#> [115] "TOP2A"     "CD2"       "SATB1"     "CASP3"     "ANKH"      "FASLG"    
#> [121] "IL6"       "CASP4"     "CASP1"     "BAX"       "FAS"       "TNFRSF12A"
#> [127] "DNM1L"     "IFNB1"     "DNAJA1"    "HMGB2"     "DAP"       "CD14"     
#> [133] "LMNA"      "BIRC3"     "ISG20"     "DFFA"      "TNF"       "CCNA1"    
#> [139] "SAT1"      "PRF1"      "TAP1"      "IL1A"      "EBP"       "BMP2"     
#> [145] "APP"       "IFNGR1"    "CD38"      "TSPO"      "CASP2"     "ANXA1"    
#> [151] "CFLAR"     "PPP3R1"    "CDC25B"    "BID"       "BCL2L10"   "SOD2"     
#> [157] "BTG3"      "F2"        "H1-0"      "PLPPR4"    "SC5D"     
#> 
#> $all$HALLMARK_ADIPOGENESIS
#>   [1] "CMBL"     "REEP6"    "REEP5"    "CYP4B1"   "CCNG2"    "PTGER3"  
#>   [7] "ADCY6"    "BAZ2A"    "COQ5"     "FAH"      "TOB1"     "SULT1A1" 
#>  [13] "PRDX3"    "QDPR"     "OMD"      "SLC27A1"  "HSPB8"    "ELMOD3"  
#>  [19] "LTC4S"    "SPARCL1"  "CRAT"     "RETSAT"   "LIPE"     "PDCD4"   
#>  [25] "UQCRQ"    "ADIPOR2"  "EPHX2"    "GPHN"     "FZD4"     "NDUFA5"  
#>  [31] "ABCB8"    "RREB1"    "MAP4K3"   "ALDOA"    "CD302"    "STOM"    
#>  [37] "DHRS7"    "NKIRAS1"  "LIFR"     "ACADS"    "BCKDHA"   "SCP2"    
#>  [43] "ADIPOQ"   "GPX4"     "GPD2"     "MGLL"     "FABP4"    "ITSN1"   
#>  [49] "ARAF"     "ESYT1"    "PPP1R15B" "DHRS7B"   "CPT2"     "GHITM"   
#>  [55] "CAT"      "GPAM"     "PFKFB3"   "SNCG"     "NMT1"     "CD36"    
#>  [61] "UQCR11"   "CIDEA"    "MRAP"     "ACLY"     "PHLDB1"   "CD151"   
#>  [67] "COX6A1"   "ITIH5"    "CHUK"     "AGPAT3"   "ABCA1"    "DNAJB9"  
#>  [73] "SORBS1"   "LPCAT3"   "PPARG"    "LEP"      "ITGA7"    "ECHS1"   
#>  [79] "LAMA4"    "RAB34"    "VEGFB"    "G3BP2"    "ELOVL6"   "DBT"     
#>  [85] "DGAT1"    "CS"       "BCL6"     "COL15A1"  "UCK1"     "APLP2"   
#>  [91] "ALDH2"    "DRAM2"    "IDH3A"    "UCP2"     "RNF11"    "ECH1"    
#>  [97] "UBC"      "PPM1B"    "TST"      "RTN3"     "SSPN"     "COX8A"   
#> [103] "UBQLN1"   "ACADM"    "SDHC"     "HADH"     "BCL2L13"  "ENPP2"   
#> [109] "HIBCH"    "ETFB"     "PEMT"     "GRPEL1"   "SOD1"     "GBE1"    
#> [115] "GPX3"     "PEX14"    "DNAJC15"  "NDUFAB1"  "LPL"      "ACOX1"   
#> [121] "MYLK"     "APOE"     "UQCRC1"   "GADD45A"  "PIM3"     "SLC25A1" 
#> [127] "NDUFB7"   "STAT5A"   "ACO2"     "SLC25A10" "PHYH"     "AIFM1"   
#> [133] "COL4A1"   "CYC1"     "DECR1"    "C3"       "IDH1"     "SCARB1"  
#> [139] "SAMM50"   "ARL4A"    "ACADL"    "NDUFS3"   "SLC1A5"   "UQCR10"  
#> [145] "SDHB"     "MGST3"    "RETN"     "DLD"      "JAGN1"    "ORM1"    
#> [151] "IMMT"     "YWHAG"    "ACAA2"    "SLC19A1"  "PFKL"     "CHCHD10" 
#> [157] "ANGPT1"   "MTCH2"    "SUCLG1"   "ANGPTL4"  "RIOK3"    "COX7B"   
#> [163] "IDH3G"    "DDT"      "PREB"     "PTCD3"    "TALDO1"   "DLAT"    
#> [169] "CDKN2C"   "ME1"      "COQ9"     "ESRRA"    "MDH2"     "TKT"     
#> [175] "TANK"     "POR"      "DHCR7"    "IFNGR1"   "COQ3"     "AK2"     
#> [181] "MRPL15"   "PLIN2"    "CMPK1"    "PGM1"     "ATP1B3"   "SLC5A6"  
#> [187] "MCCC1"    "ATL2"     "ADIG"     "ATP5PO"   "CAVIN1"   "CAVIN2"  
#> [193] "GPAT4"    "MIGA2"    "MTARC2"   "NABP1"    "RMDN3"    "SLC66A3" 
#> [199] "SOWAHC"   "SQOR"    
#> 
#> $all$HALLMARK_PANCREAS_BETA_CELLS
#>  [1] "ABCC8"   "SPCS1"   "SRP14"   "HNF1A"   "G6PC2"   "STXBP1"  "SYT13"  
#>  [8] "INSM1"   "SRP9"    "GCK"     "PCSK2"   "SCGN"    "ELP4"    "FOXO1"  
#> [15] "LMO2"    "MAFB"    "CHGA"    "NKX6-1"  "PKLR"    "PCSK1"   "NKX2-2" 
#> [22] "DPP4"    "ISL1"    "SEC11A"  "VDR"     "AKT3"    "DCX"     "FOXA2"  
#> [29] "SRPRB"   "PAX6"    "PAK3"    "PDX1"    "GCG"     "IAPP"    "INS"    
#> [36] "NEUROD1" "NEUROG3" "PAX4"    "SLC2A2"  "SST"    
#> 
#> $all$HALLMARK_COAGULATION
#>   [1] "TMPRSS6"  "ACOX2"    "RAPGEF3"  "HPN"      "PLAT"     "PRSS23"  
#>   [7] "HMGCS2"   "TIMP3"    "CTSO"     "CFD"      "MST1"     "ANG"     
#>  [13] "CFB"      "F2RL2"    "HTRA1"    "THBD"     "ISCU"     "CRIP2"   
#>  [19] "SERPINC1" "COMP"     "SPARC"    "ITGA2"    "THBS1"    "FBN1"    
#>  [25] "MASP2"    "CASP9"    "APOA1"    "MMP11"    "MMP10"    "LTA4H"   
#>  [31] "LRP1"     "VWF"      "GNG12"    "PROZ"     "F10"      "CTSK"    
#>  [37] "PROC"     "PDGFB"    "F8"       "S100A13"  "SERPINE1" "FN1"     
#>  [43] "SERPINA1" "CLU"      "MMP2"     "ITGB3"    "CD9"      "TFPI2"   
#>  [49] "PECAM1"   "OLR1"     "CFH"      "PEF1"     "USP11"    "TIMP1"   
#>  [55] "F3"       "HNF4A"    "FGA"      "WDR1"     "DUSP6"    "APOC2"   
#>  [61] "KLF7"     "DUSP14"   "F12"      "MBL2"     "FGG"      "PROS1"   
#>  [67] "GNB2"     "ARF4"     "MSRB2"    "MMP14"    "APOC1"    "KLKB1"   
#>  [73] "DPP4"     "C8B"      "RAC1"     "A2M"      "ADAM9"    "CTSH"    
#>  [79] "ITIH1"    "CSRP1"    "FURIN"    "PLAU"     "RABIF"    "GP1BA"   
#>  [85] "C2"       "MMP3"     "CFI"      "SERPING1" "P2RY1"    "DCT"     
#>  [91] "CTSB"     "MEP1A"    "GSN"      "BMP1"     "C3"       "CPB2"    
#>  [97] "C1QA"     "MMP8"     "CAPN5"    "PLEK"     "MMP9"     "MMP15"   
#> [103] "LEFTY2"   "SIRT2"    "CTSE"     "GDA"      "C1S"      "CAPN2"   
#> [109] "RGN"      "C1R"      "TF"       "SERPINB2" "C8G"      "S100A1"  
#> [115] "FYN"      "LAMP2"    "MAFF"     "LGMN"     "MMP1"     "ANXA1"   
#> [121] "PREP"     "SH2B2"    "MMP7"     "KLK8"     "APOC3"    "C8A"     
#> [127] "C9"       "CPN1"     "CPQ"      "CTSV"     "F11"      "F13B"    
#> [133] "F2"       "F9"       "GP9"      "HRG"      "PF4"      "PLG"     
#> 
#> $all$HALLMARK_MYOGENESIS
#>   [1] "SPDEF"     "ADCY9"     "ERBB3"     "BHLHE40"   "KCNH1"     "OCEL1"    
#>   [7] "STC2"      "REEP1"     "CACNG1"    "PPP1R3C"   "CACNA1H"   "ITGB5"    
#>  [13] "HRC"       "VIPR1"     "RB1"       "MEF2D"     "CASQ1"     "HSPB8"    
#>  [19] "CFD"       "PVALB"     "CAMK2B"    "CRAT"      "TNNT1"     "ATP6AP1"  
#>  [25] "BAG1"      "LAMA2"     "SGCD"      "TSC2"      "NQO1"      "PSEN2"    
#>  [31] "COX7A1"    "CASQ2"     "MYOG"      "SGCG"      "MYH11"     "COL6A3"   
#>  [37] "FKBP1B"    "SPARC"     "AGL"       "PFKM"      "MB"        "IGF1"     
#>  [43] "MAPRE3"    "DMPK"      "BDKRB2"    "PGAM2"     "COL1A1"    "LDB3"     
#>  [49] "COL3A1"    "TGFB1"     "MYL3"      "MYO1C"     "AEBP1"     "SYNGR2"   
#>  [55] "TNNT3"     "SORBS3"    "CD36"      "SPTAN1"    "CDKN1A"    "SCD"      
#>  [61] "FHL1"      "FST"       "ADAM12"    "SGCA"      "MYL7"      "IGFBP7"   
#>  [67] "APOD"      "MYH1"      "MYH8"      "FXYD1"     "AK1"       "ENO3"     
#>  [73] "GABARAPL2" "CDH13"     "CLU"       "MYH4"      "CHRNA1"    "EIF4A2"   
#>  [79] "SORBS1"    "AKT2"      "GADD45B"   "ITGA7"     "FOXO4"     "SH2B1"    
#>  [85] "APLNR"     "TNNC1"     "MYOM2"     "PYGM"      "COL15A1"   "GAA"      
#>  [91] "FABP3"     "GNAO1"     "COX6A2"    "CTF1"      "COL6A2"    "MYH2"     
#>  [97] "MYL4"      "HDAC5"     "ANKRD2"    "MYOM1"     "SH3BGR"    "SSPN"     
#> [103] "AGRN"      "DTNA"      "CKB"       "CHRNB1"    "MEF2C"     "CKMT2"    
#> [109] "MYH3"      "GPX3"      "SLN"       "TAGLN"     "PDLIM7"    "FLII"     
#> [115] "PLXNB2"    "LSP1"      "PKIA"      "ACTC1"     "ACTA1"     "RIT1"     
#> [121] "MYLK"      "MEF2A"     "DAPK2"     "TNNC2"     "MYL6B"     "TCAP"     
#> [127] "FGF2"      "PTGIS"     "ACTN3"     "PDE4DIP"   "GJA5"      "NOS1"     
#> [133] "KCNH2"     "MYOZ1"     "TPM3"      "ATP2A1"    "MYBPH"     "MAPK12"   
#> [139] "PC"        "NAV2"      "DES"       "HSPB2"     "TPM2"      "GSN"      
#> [145] "ACTN2"     "PPFIA4"    "MYH7"      "ITGB1"     "TPD52L1"   "MYBPC3"   
#> [151] "PICK1"     "SIRT2"     "NCAM1"     "HBEGF"     "TNNI1"     "MYH9"     
#> [157] "SVIL"      "ACSL1"     "ABLIM1"    "COL4A2"    "CNN3"      "ACHE"     
#> [163] "PTP4A3"    "IGFBP3"    "SMTN"      "ITGB4"     "SLC6A8"    "CKM"      
#> [169] "SOD3"      "SPEG"      "EFS"       "BIN1"      "TNNI2"     "RYR1"     
#> [175] "SPHK1"     "APP"       "FDPS"      "KIFC3"     "TNNT2"     "DMD"      
#> [181] "PRNP"      "TEAD4"     "WWTR1"     "NOTCH1"    "KLF5"      "EPHB3"    
#> [187] "SCHIP1"    "CRYAB"     "MRAS"      "LPIN1"     "IFRD1"     "CAV3"     
#> [193] "CHRNG"     "CSRP3"     "DENND2B"   "LARGE1"    "MYF6"      "MYL1"     
#> [199] "MYL11"     "MYL2"     
#> 
#> $all$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
#>   [1] "RHOB"      "IGFBP4"    "MATN3"     "FUCA1"     "IGFBP2"    "ITGB5"    
#>   [7] "TIMP3"     "GJA1"      "FSTL3"     "ECM1"      "CAP2"      "HTRA1"    
#>  [13] "NTM"       "LAMA2"     "SGCD"      "LAMA3"     "LOXL1"     "PDGFRB"   
#>  [19] "ECM2"      "GPC1"      "SDC4"      "SLIT3"     "EDIL3"     "PMP22"    
#>  [25] "POSTN"     "MAGEE1"    "COMP"      "SGCG"      "SPOCK1"    "FMOD"     
#>  [31] "AREG"      "COL6A3"    "SPARC"     "CDH11"     "ENO2"      "ITGA2"    
#>  [37] "SFRP4"     "COL12A1"   "THBS1"     "COL1A2"    "FBN1"      "CDH6"     
#>  [43] "DCN"       "COL1A1"    "SLIT2"     "COL3A1"    "TGFB1"     "COL8A2"   
#>  [49] "LRRC15"    "CXCL12"    "MSX1"      "LRP1"      "LUM"       "COL5A2"   
#>  [55] "FAP"       "COL5A1"    "TGFBR3"    "BDNF"      "COL5A3"    "RGS4"     
#>  [61] "SCG2"      "EFEMP2"    "ABI3BP"    "PTHLH"     "MFAP5"     "GEM"      
#>  [67] "ADAM12"    "NID2"      "WNT5A"     "VCAN"      "QSOX1"     "THBS2"    
#>  [73] "SERPINE1"  "THY1"      "OXTR"      "FN1"       "ITGAV"     "COL16A1"  
#>  [79] "MXRA5"     "MMP2"      "COL7A1"    "PCOLCE"    "GADD45B"   "BASP1"    
#>  [85] "ITGB3"     "COL11A1"   "DAB2"      "PRRX1"     "TFPI2"     "ELN"      
#>  [91] "INHBA"     "VEGFC"     "COL6A2"    "DST"       "FBLN1"     "ACTA2"    
#>  [97] "MATN2"     "TIMP1"     "SGCB"      "MEST"      "PDLIM4"    "FSTL1"    
#> [103] "JUN"       "BGN"       "GREM1"     "ANPEP"     "TNFRSF11B" "FERMT2"   
#> [109] "MYL9"      "CTHRC1"    "MGP"       "TAGLN"     "ITGA5"     "SNAI2"    
#> [115] "FZD8"      "NOTCH2"    "DPYSL3"    "LOXL2"     "LOX"       "MMP14"    
#> [121] "CD44"      "WIPF1"     "MYLK"      "FBN2"      "EMP3"      "TPM1"     
#> [127] "GADD45A"   "LGALS1"    "PMEPA1"    "FGF2"      "NNMT"      "TGM2"     
#> [133] "CADM1"     "LAMC1"     "FBLN2"     "CDH2"      "COL4A1"    "APLP1"    
#> [139] "SPP1"      "MMP3"      "NT5E"      "GAS1"      "TNC"       "LAMA1"    
#> [145] "TPM2"      "BMP1"      "CAPG"      "TGFBI"     "COPA"      "ITGB1"    
#> [151] "ID2"       "FBLN5"     "CALD1"     "COL4A2"    "SDC1"      "IL15"     
#> [157] "TPM4"      "GLIPR1"    "IL6"       "CRLF1"     "FAS"       "TNFRSF12A"
#> [163] "VIM"       "VCAM1"     "PLOD2"     "PCOLCE2"   "FOXC2"     "CD59"     
#> [169] "LAMC2"     "IGFBP3"    "SNTB1"     "SLC6A8"    "FLNA"      "TNFAIP3"  
#> [175] "PFN2"      "SAT1"      "PPIB"      "VEGFA"     "MMP1"      "SERPINH1" 
#> [181] "IL32"      "PLAUR"     "SERPINE2"  "GPX7"      "DKK1"      "CXCL6"    
#> [187] "CALU"      "PLOD3"     "PTX3"      "SFRP1"     "MCM7"      "CXCL1"    
#> [193] "PVR"       "PLOD1"     "CCN1"      "CCN2"      "COLGALT1"  "CXCL8"    
#> [199] "P3H1"      "PRSS2"    
#> 
#> $all$HALLMARK_TGF_BETA_SIGNALING
#>  [1] "PPM1A"    "SMAD3"    "APC"      "RAB31"    "XIAP"     "BCAR3"   
#>  [7] "THBS1"    "LTBP2"    "TGFB1"    "SLC20A1"  "RHOA"     "ARID4B"  
#> [13] "NCOR2"    "BMPR2"    "CDK9"     "SMAD7"    "TJP1"     "SMAD1"   
#> [19] "SERPINE1" "KLF10"    "TGIF1"    "SKIL"     "UBE2D3"   "ACVR1"   
#> [25] "JUNB"     "HIPK2"    "CDH1"     "NOG"      "TRIM33"   "ENG"     
#> [31] "SMAD6"    "ID3"      "HDAC1"    "PPP1CA"   "CTNNB1"   "CDKN1C"  
#> [37] "PMEPA1"   "TGFBR1"   "FURIN"    "BMPR1A"   "ID1"      "SPTBN1"  
#> [43] "PPP1R15A" "LEFTY2"   "ID2"      "FNTA"     "IFNGR2"   "MAP3K7"  
#> [49] "SKI"      "SMURF1"   "SMURF2"   "FKBP1A"   "BMP2"     "WWTR1"   
#> 
#> 
#> $le
#> $le$HALLMARK_E2F_TARGETS
#>   [1] "RAD50"    "WDR90"    "PAN2"     "UBR7"     "GSPT1"    "MLH1"    
#>   [7] "PDS5B"    "CDKN1B"   "RFC1"     "LUC7L3"   "BRCA1"    "CDKN1A"  
#>  [13] "PPM1D"    "RBBP7"    "PMS2"     "CTCF"     "TP53"     "IPO7"    
#>  [19] "RAD51C"   "DCTPP1"   "CBX5"     "RPA3"     "NBN"      "WEE1"    
#>  [25] "DUT"      "SHMT1"    "NAP1L1"   "TUBG1"    "TIMELESS" "STAG1"   
#>  [31] "DCK"      "BARD1"    "PA2G4"    "KIF22"    "LIG1"     "RAD1"    
#>  [37] "PRPS1"    "PNN"      "NME1"     "NUP107"   "SMC3"     "POLE4"   
#>  [43] "EIF2S1"   "POLE"     "SMC1A"    "PSMC3IP"  "RAD21"    "NAA38"   
#>  [49] "HUS1"     "RPA2"     "CSE1L"    "POLD2"    "ING3"     "MYC"     
#>  [55] "TOP2A"    "SNRPB"    "PRKDC"    "NOLC1"    "NUDT21"   "RACGAP1" 
#>  [61] "PAICS"    "MXD3"     "BRCA2"    "RPA1"     "ESPL1"    "SLBP"    
#>  [67] "ATAD2"    "HMGB2"    "UNG"      "NUP153"   "HELLS"    "TFRC"    
#>  [73] "SPAG5"    "NOP56"    "BRMS1L"   "PSIP1"    "CIT"      "POLA2"   
#>  [79] "POLD3"    "RNASEH2A" "TRA2B"    "CENPM"    "PPP1R8"   "CDK4"    
#>  [85] "ASF1B"    "XPO1"     "TIPIN"    "TK1"      "E2F8"     "SPC25"   
#>  [91] "PCNA"     "SPC24"    "EXOSC8"   "TMPO"     "CDKN2C"   "HMMR"    
#>  [97] "HNRNPD"   "LMNB1"    "DNMT1"    "PLK4"     "ILF3"     "SMC4"    
#> [103] "CDK1"     "RAN"      "POP7"     "HMGB3"    "RFC3"     "UBE2T"   
#> [109] "XRCC6"    "GINS1"    "RRM2"     "CDKN3"    "TACC3"    "BUB1B"   
#> [115] "POLD1"    "ZW10"     "AURKA"    "MCM4"     "KPNA2"    "MCM2"    
#> [121] "CENPE"    "GINS4"    "MAD2L1"   "TBRG4"    "CKS2"     "NUP205"  
#> [127] "SUV39H1"  "GINS3"    "KIF4A"    "DCLRE1B"  "EED"      "UBE2S"   
#> [133] "TCF19"    "AK2"      "ASF1A"    "KIF18B"   "MKI67"    "DSCC1"   
#> [139] "PTTG1"    "MTHFD2"   "BIRC5"    "CDC25B"   "SMC6"     "MCM3"    
#> [145] "SSRP1"    "DIAPH3"   "MYBL2"    "DLGAP5"   "USP1"     "TRIP13"  
#> [151] "ANP32E"   "MCM6"     "RAD51AP1" "MSH2"     "EZH2"     "CHEK2"   
#> [157] "PLK1"     "SYNCRIP"  "CDCA3"    "DEPDC1"   "RFC2"     "CCNB2"   
#> [163] "CKS1B"    "MELK"     "CDKN2A"   "DONSON"   "PHF5A"    "PRIM2"   
#> [169] "KIF2C"    "AURKB"    "LYAR"     "PRDX4"    "CDC25A"   "MCM7"    
#> [175] "MCM5"     "NCAPD2"   "CHEK1"    "STMN1"    "DEK"      "RANBP1"  
#> [181] "NASP"     "CDCA8"    "LBR"      "TUBB"     "CDC20"    "HMGA1"   
#> [187] "CCNE1"    "CCP110"   "CNOT9"    "CTPS1"    "DDX39A"   "H2AX"    
#> [193] "H2AZ1"    "JPT1"     "MMS22L"   "MRE11"    "ORC2"     "ORC6"    
#> [199] "SRSF1"    "SRSF2"   
#> 
#> $le$HALLMARK_ESTROGEN_RESPONSE_EARLY
#>   [1] "CA12"     "THSD4"    "MLPH"     "XBP1"     "ANXA9"    "TFF1"    
#>   [7] "GREB1"    "TFF3"     "MYB"      "SLC22A5"  "ABAT"     "CELSR1"  
#>  [13] "GFRA1"    "SLC39A6"  "MAPT"     "KCNK15"   "KDM4B"    "SYBU"    
#>  [19] "ADCY9"    "WFS1"     "PGR"      "AR"       "SLC7A2"   "BCL2"    
#>  [25] "UGCG"     "RARA"     "IL6ST"    "TTC39A"   "IGF1R"    "BHLHE40" 
#>  [31] "SCNN1A"   "MAST4"    "MED13L"   "IGFBP4"   "RAB17"    "LRIG1"   
#>  [37] "KRT18"    "MREG"     "SEMA3B"   "TPBG"     "RET"      "SLC19A2" 
#>  [43] "ARL3"     "CISH"     "ELOVL5"   "CCND1"    "STC2"     "REEP1"   
#>  [49] "TJP3"     "SIAH2"    "PDZK1"    "ABCA3"    "ELOVL2"   "SLC27A2" 
#>  [55] "ADCY1"    "NRIP1"    "AFF1"     "PEX11A"   "AMFR"     "MYOF"    
#>  [61] "ABHD2"    "DHRS2"    "ITPK1"    "TOB1"     "PRSS23"   "SULT2B1" 
#>  [67] "SLC1A4"   "ASB13"    "CELSR2"   "CBFA2T3"  "FLNB"     "HSPB8"   
#>  [73] "SEC14L2"  "KRT8"     "MUC1"     "ELF1"     "SNX24"    "GJA1"    
#>  [79] "OLFM1"    "CANT1"    "DYNLT3"   "TIPARP"   "EGR3"     "RAB31"   
#>  [85] "NPY1R"    "BLVRB"    "FKBP4"    "BAG1"     "KRT19"    "OVOL2"   
#>  [91] "FASN"     "P2RY2"    "PMAIP1"   "SLC1A1"   "MSMB"     "CALCR"   
#>  [97] "AREG"     "OLFML3"   "FOS"      "RAPGEFL1" "NADSYN1"  "UNC119"  
#> [103] "NBL1"     "CXCL12"   "PTGES"    "MPPED2"   "SYNGR1"   "RHOBTB3" 
#> [109] "NCOR2"    "DLC1"     "CHPT1"    "RBBP8"    "GLA"      "KLF4"    
#> [115] "IL17RB"   "SOX3"     "RHOD"     "KLF10"    "RPS6KA2"  "TSKU"    
#> [121] "HES1"     "MICB"     "ESRP2"    "INHBB"    "CLDN7"    "FHL2"    
#> [127] "SLC24A3"  "SYT12"    "GAB2"     "TMEM164"  "WWC1"     "JAK2"    
#> [133] "KRT13"    "FDFT1"    "MED24"    "RASGRP1"  "SH3BP5"   "TMPRSS3" 
#> [139] "INPP5F"   "SLC37A1"  "B4GALT1"  "FRK"      "ALDH3B1"  "CD44"    
#> [145] "ZNF185"   "AQP3"     "PDLIM3"   "AKAP1"    "PAPSS2"   "TIAM1"   
#> [151] "TGM2"     "FARP1"    "NAV2"     "ENDOD1"   "ELF3"     "FKBP5"   
#> [157] "ADD3"     "RRP12"    "SCARB1"   "DHRS3"    "MYBL1"    "ISG20L2" 
#> [163] "CLIC3"    "TPD52L1"  "SLC26A2"  "MYC"      "TGIF2"    "SVIL"    
#> [169] "TUBB2B"   "SFN"      "ABLIM1"   "KRT15"    "MYBBP1A"  "PPIF"    
#> [175] "CYP26B1"  "PODXL"    "TFAP2C"   "SLC2A1"   "NXT1"     "BCL11B"  
#> [181] "KLK10"    "DHCR7"    "CALB2"    "SLC16A1"  "OPN3"     "HR"      
#> [187] "LAD1"     "SLC7A5"   "KCNK5"    "FOXC1"    "CCN5"     "DEPTOR"  
#> [193] "EEIG1"    "FCMR"     "KAZN"     "MINDY1"   "NHERF1"   "PLAAT3"  
#> [199] "RETREG1"  "TBC1D30" 
#> 
#> $le$HALLMARK_G2M_CHECKPOINT
#>   [1] "TLE3"     "CCND1"    "NUMA1"    "ARID4A"   "SMAD3"    "CUL3"    
#>   [7] "PURA"     "RPS6KA5"  "GSPT1"    "MNAT1"    "SLC38A1"  "BUB3"    
#>  [13] "YTHDC1"   "SLC12A2"  "CCNT1"    "PDS5B"    "CDKN1B"   "ATRX"    
#>  [19] "TGFB1"    "EGF"      "PRMT5"    "LIG3"     "CDC27"    "ABL1"    
#>  [25] "CTCF"     "PAFAH1B1" "ODF2"     "TOP1"     "NUP98"    "BCL3"    
#>  [31] "RAD23B"   "G3BP1"    "SS18"     "CUL5"     "HMGN2"    "STAG1"   
#>  [37] "HOXC10"   "HNRNPU"   "BARD1"    "UPF1"     "KIF22"    "NOTCH2"  
#>  [43] "HSPA8"    "KPNB1"    "CUL4A"    "KIF5B"    "WRN"      "MEIS1"   
#>  [49] "GINS2"    "CBX1"     "FOXN3"    "POLE"     "SMC1A"    "MEIS2"   
#>  [55] "SQLE"     "SMARCC1"  "RAD21"    "MAPK14"   "HUS1"     "RBM14"   
#>  [61] "RPA2"     "TNPO2"    "CCNF"     "MT2A"     "MYC"      "SMC2"    
#>  [67] "DTYMK"    "TOP2A"    "EWSR1"    "KIF20B"   "NOLC1"    "DR1"     
#>  [73] "RACGAP1"  "TRAIP"    "BRCA2"    "ESPL1"    "RBL1"     "PBK"     
#>  [79] "CASP8AP2" "NEK2"     "ATF5"     "POLA2"    "RASAL2"   "MARCKS"  
#>  [85] "CDC6"     "TRA2B"    "CHMP1A"   "CDK4"     "XPO1"     "NUSAP1"  
#>  [91] "CHAF1A"   "NUP50"    "TMPO"     "CDKN2C"   "HMMR"     "HNRNPD"  
#>  [97] "LMNB1"    "PTTG3P"   "HIF1A"    "SFPQ"     "PLK4"     "ILF3"    
#> [103] "SMC4"     "SAP30"    "CDK1"     "POLQ"     "CUL1"     "FANCC"   
#> [109] "SLC7A1"   "E2F1"     "MTF2"     "HMGB3"    "PRC1"     "HIRA"    
#> [115] "TROAP"    "KIF11"    "CDKN3"    "TACC3"    "NCL"      "FBXO5"   
#> [121] "KIF15"    "AURKA"    "KPNA2"    "MCM2"     "CENPE"    "MAD2L1"  
#> [127] "CKS2"     "SUV39H1"  "KIF4A"    "KIF23"    "EFNA5"    "UBE2S"   
#> [133] "CENPF"    "MKI67"    "PTTG1"    "DMD"      "UBE2C"    "PML"     
#> [139] "CDC7"     "BIRC5"    "SNRPD1"   "CDC25B"   "INCENP"   "EXO1"    
#> [145] "CDC45"    "MCM3"     "MYBL2"    "E2F2"     "MCM6"     "TFDP1"   
#> [151] "TPX2"     "E2F4"     "EZH2"     "PLK1"     "BUB1"     "CCNA2"   
#> [157] "SYNCRIP"  "RAD54L"   "CCNB2"    "CKS1B"    "KATNA1"   "NDC80"   
#> [163] "DBF4"     "STIL"     "PRIM2"    "KIF2C"    "AURKB"    "CDC25A"  
#> [169] "MCM5"     "CHEK1"    "SLC7A5"   "STMN1"    "TTK"      "ODC1"    
#> [175] "NASP"     "LBR"      "AMD1"     "UCK2"     "CENPA"    "DKC1"    
#> [181] "CDC20"    "HMGA1"    "E2F3"     "DDX39A"   "H2AX"     "H2AZ1"   
#> [187] "H2AZ2"    "H2BC12"   "JPT1"     "KMT5A"    "KNL1"     "MAP3K20" 
#> [193] "NSD2"     "ORC5"     "ORC6"     "PRP4K"    "SRSF1"    "SRSF10"  
#> [199] "SRSF2"    "TENT4A"  
#> 
#> $le$HALLMARK_MYC_TARGETS_V1
#>   [1] "FAM120A"   "PRDX3"     "GSPT1"     "IMPDH2"    "BUB3"      "PTGES3"   
#>   [7] "GOT2"      "CANX"      "PSMC6"     "DHX15"     "EIF4E"     "ETF1"     
#>  [13] "CCT2"      "RAD23B"    "AP3S1"     "RSL1D1"    "HNRNPC"    "EIF4G2"   
#>  [19] "EIF1AX"    "NHP2"      "G3BP1"     "RPS5"      "TUFM"      "VDAC1"    
#>  [25] "DUT"       "RPL14"     "MRPL23"    "EIF3J"     "GLO1"      "TRIM28"   
#>  [31] "HNRNPA1"   "NAP1L1"    "RRM1"      "RPL22"     "PSMA6"     "MRPS18B"  
#>  [37] "APEX1"     "SF3A1"     "RPL18"     "RUVBL2"    "HNRNPU"    "PA2G4"    
#>  [43] "ERH"       "NCBP2"     "RPL34"     "KPNB1"     "NDUFAB1"   "CSTF2"    
#>  [49] "COPS5"     "NPM1"      "PSMD3"     "PSMB3"     "HNRNPR"    "RPLP0"    
#>  [55] "NME1"      "PSMA1"     "CLNS1A"    "RPL6"      "RPS2"      "PRPF31"   
#>  [61] "EIF2S1"    "CYC1"      "YWHAE"     "EXOSC7"    "GNL3"      "SLC25A3"  
#>  [67] "UBE2E1"    "CBX3"      "VDAC3"     "LDHA"      "SMARCC1"   "HNRNPA2B1"
#>  [73] "RPS3"      "PSMD8"     "RNPS1"     "C1QBP"     "TXNL4A"    "NCBP1"    
#>  [79] "POLE3"     "POLD2"     "PWP1"      "CNBP"      "HSP90AB1"  "EIF4A1"   
#>  [85] "EIF4H"     "CDK2"      "RRP9"      "MYC"       "PCBP1"     "PSMD1"    
#>  [91] "SNRPB2"    "SNRPD2"    "RPS6"      "NOLC1"     "ACP1"      "DDX21"    
#>  [97] "UBE2L3"    "SET"       "PSMA2"     "PSMC4"     "EEF1B2"    "VBP1"     
#> [103] "PPIA"      "AIMP2"     "HNRNPA3"   "HPRT1"     "NOP56"     "TARDBP"   
#> [109] "HSPE1"     "EIF3D"     "TRA2B"     "PRPS2"     "CCT3"      "PHB2"     
#> [115] "CDK4"      "XPO1"      "RPS10"     "SNRPD3"    "XPOT"      "PCNA"     
#> [121] "HSPD1"     "HNRNPD"    "CCT7"      "PSMA4"     "SSB"       "LSM7"     
#> [127] "COX5A"     "PGK1"      "ABCE1"     "NOP16"     "PSMA7"     "RAN"      
#> [133] "CUL1"      "SSBP1"     "HDDC2"     "EIF2S2"    "XRCC6"     "FBL"      
#> [139] "PSMD14"    "STARD7"    "PABPC1"    "MCM4"      "KPNA2"     "MCM2"     
#> [145] "MAD2L1"    "CAD"       "SRM"       "EIF3B"     "DDX18"     "UBA2"     
#> [151] "TCP1"      "SNRPA"     "TYMS"      "CCT5"      "U2AF1"     "MRPL9"    
#> [157] "SNRPD1"    "PABPC4"    "PSMD7"     "SF3B3"     "SNRPA1"    "CDC45"    
#> [163] "YWHAQ"     "USP1"      "PPM1G"     "MCM6"      "TFDP1"     "CCNA2"    
#> [169] "SYNCRIP"   "SNRPG"     "CCT4"      "HDGF"      "RFC4"      "PRDX4"    
#> [175] "MCM7"      "MCM5"      "HDAC2"     "LSM2"      "DEK"       "RANBP1"   
#> [181] "ODC1"      "ILF2"      "SERBP1"    "PSMB2"     "CDC20"     "SRPK1"    
#> [187] "IFRD1"     "CTPS1"     "EPRS1"     "H2AZ1"     "IARS1"     "KARS1"    
#> [193] "ORC2"      "PHB1"      "RACK1"     "SRSF1"     "SRSF2"     "SRSF3"    
#> [199] "SRSF7"     "TOMM70"   
#> 
#> $le$HALLMARK_ALLOGRAFT_REJECTION
#>   [1] "IKBKB"    "EIF4G3"   "APBB1"    "F2R"      "TLR3"     "TPD52"   
#>   [7] "TGFB1"    "ELANE"    "CARTPT"   "IRF7"     "BRCA1"    "KRT1"    
#>  [13] "LY86"     "THY1"     "AKT1"     "INHBB"    "INHBA"    "BCL3"    
#>  [19] "IL16"     "ITGAL"    "STAB1"    "RPS3A"    "RPL9"     "TIMP1"   
#>  [25] "NLRP3"    "CSF1"     "PTPN6"    "EIF3J"    "UBE2N"    "DYRK3"   
#>  [31] "JAK2"     "MBL2"     "CCND2"    "GPR65"    "EIF3A"    "IGSF6"   
#>  [37] "IL4R"     "GALNT1"   "CCL11"    "RPS9"     "PRKCG"    "TRAF2"   
#>  [43] "NPM1"     "HLA-DOA"  "CCL22"    "IL2"      "CD74"     "IL11"    
#>  [49] "BCL10"    "NME1"     "GCNT1"    "HLA-DMB"  "CRTAM"    "CCL19"   
#>  [55] "ST8SIA4"  "ACVR2A"   "CFP"      "C2"       "DEGS1"    "CD4"     
#>  [61] "FCGR2B"   "IRF8"     "HLA-DRA"  "B2M"      "CD28"     "CD1D"    
#>  [67] "HLA-DQA1" "HLA-DMA"  "GBP2"     "PTPRC"    "SPI1"     "PSMB10"  
#>  [73] "EREG"     "PRKCB"    "IL10"     "LCP2"     "CD8A"     "HCLS1"   
#>  [79] "CAPG"     "TRAT1"    "STAT1"    "CD40LG"   "IL12B"    "CCR5"    
#>  [85] "IL1B"     "CD86"     "MMP9"     "CD3G"     "TGFB2"    "IL18"    
#>  [91] "STAT4"    "WAS"      "FGR"      "CXCR3"    "CCR2"     "IL2RB"   
#>  [97] "NOS2"     "ITK"      "SOCS5"    "HLA-G"    "NCF4"     "CD3E"    
#> [103] "GZMA"     "CD80"     "TAPBP"    "SRGN"     "ITGB2"    "CD2"     
#> [109] "CD3D"     "CD96"     "FASLG"    "IL15"     "CCND3"    "ZAP70"   
#> [115] "IL6"      "SIT1"     "CD8B"     "IL13"     "LTB"      "FAS"     
#> [121] "HLA-E"    "IL12RB1"  "MAP4K1"   "CCR1"     "CD247"    "ACHE"    
#> [127] "IL7"      "IFNGR2"   "CTSS"     "CXCL9"    "BCAT1"    "CXCL13"  
#> [133] "NCR1"     "LY75"     "LIF"      "ETS1"     "RPL39"    "CD40"    
#> [139] "MAP3K7"   "TLR1"     "ELF4"     "HLA-A"    "CCL4"     "IL2RG"   
#> [145] "GLMN"     "KLRD1"    "EIF3D"    "HDAC9"    "CCL2"     "TLR6"    
#> [151] "CD79A"    "CD47"     "LCK"      "ICOSLG"   "FLNA"     "TNF"     
#> [157] "TLR2"     "IL18RAP"  "CCL5"     "CD7"      "MRPL3"    "PRF1"    
#> [163] "ABCE1"    "TAP1"     "IRF4"     "HIF1A"    "IFNG"     "IL12A"   
#> [169] "CCL13"    "IL27RA"   "RPS19"    "ICAM1"    "UBE2D1"   "SOCS1"   
#> [175] "EIF5A"    "IFNGR1"   "MTIF2"    "IL2RA"    "ABI1"     "GZMB"    
#> [181] "TAP2"     "CCL7"     "CSK"      "NCK1"     "HLA-DOB"  "LYN"     
#> [187] "CDKN2A"   "IFNAR2"   "RIPK2"    "EGFR"     "AARS1"    "DARS1"   
#> [193] "F2"       "FYB1"     "IL4"      "IL9"      "PF4"      "RARS1"   
#> [199] "RPL3L"    "WARS1"   
#> 
#> $le$HALLMARK_MYC_TARGETS_V2
#>  [1] "SORD"      "RABEPK"    "IPO4"      "DCTPP1"    "HK2"       "PRMT3"    
#>  [7] "GRWD1"     "TMEM97"    "NOC4L"     "PA2G4"     "NPM1"      "SLC29A2"  
#> [13] "TFB2M"     "FARSA"     "GNL3"      "RRP12"     "LAS1L"     "CBX3"     
#> [19] "UTP20"     "EXOSC5"    "RRP9"      "MYC"       "NOLC1"     "IMP4"     
#> [25] "MYBBP1A"   "AIMP2"     "SLC19A1"   "UNG"       "RCL1"      "NOP56"    
#> [31] "HSPE1"     "MAP3K6"    "WDR74"     "CDK4"      "HSPD1"     "PUS1"     
#> [37] "PPRC1"     "PPAN"      "NOP16"     "DUSP2"     "PLK4"      "MRTO4"    
#> [43] "PES1"      "MCM4"      "TBRG4"     "SRM"       "SUPV3L1"   "DDX18"    
#> [49] "TCOF1"     "MPHOSPH10" "PLK1"      "NOP2"      "BYSL"      "WDR43"    
#> [55] "NIP7"      "MCM5"      "NDUFAF4"   "PHB1"     
#> 
#> $le$HALLMARK_ESTROGEN_RESPONSE_LATE
#>   [1] "CA12"       "XBP1"       "ANXA9"      "TFF1"       "AGR2"      
#>   [6] "SCUBE2"     "TFF3"       "MYB"        "SLC22A5"    "RABEP1"    
#>  [11] "CACNA2D2"   "MAPT"       "DNAJC12"    "WFS1"       "PGR"       
#>  [16] "BCL2"       "IL6ST"      "TSPAN13"    "SCNN1A"     "IGFBP4"    
#>  [21] "SEMA3B"     "TPBG"       "RET"        "ARL3"       "CISH"      
#>  [26] "ELOVL5"     "CCND1"      "SERPINA5"   "ACOX2"      "TJP3"      
#>  [31] "SIAH2"      "PDZK1"      "PTGER3"     "ABCA3"      "PRLR"      
#>  [36] "SLC27A2"    "NRIP1"      "UGDH"       "AFF1"       "TPSAB1"    
#>  [41] "AMFR"       "COX6C"      "MYOF"       "ABHD2"      "DHRS2"     
#>  [46] "DLG5"       "ITPK1"      "TOB1"       "SORD"       "PRSS23"    
#>  [51] "SULT2B1"    "SLC1A4"     "HMGCS2"     "CELSR2"     "CHST8"     
#>  [56] "FGFR3"      "CYP4F11"    "FLNB"       "EMP2"       "HSPB8"     
#>  [61] "ST6GALNAC2" "MOCS2"      "OLFM1"      "CXCL14"     "DNAJC1"    
#>  [66] "ALDH3A2"    "CPE"        "DYNLT3"     "EGR3"       "RAB31"     
#>  [71] "NPY1R"      "SERPINA3"   "BLVRB"      "FKBP4"      "BAG1"      
#>  [76] "PDCD4"      "PLAC1"      "KRT19"      "OVOL2"      "ASCL1"     
#>  [81] "DCXR"       "BATF"       "PLXNB1"     "CALCR"      "AREG"      
#>  [86] "FOS"        "RAPGEFL1"   "NBL1"       "CXCL12"     "SLC2A8"    
#>  [91] "PTGES"      "PRKAR2B"    "NCOR2"      "UNC13B"     "CHPT1"     
#>  [96] "METTL3"     "RBBP8"      "GLA"        "KLF4"       "IGSF1"     
#> [101] "IL17RB"     "SOX3"       "SERPINA1"   "JAK1"       "CA2"       
#> [106] "RPS6KA2"    "CAV1"       "ZFP36"      "SNX10"      "MICB"      
#> [111] "CD9"        "LLGL2"      "TFPI2"      "NAB2"       "GALE"      
#> [116] "TNNC1"      "SLC24A3"    "CDH1"       "SLC29A1"    "MEST"      
#> [121] "TST"        "PTPN6"      "ATP2B4"     "JAK2"       "HOMER2"    
#> [126] "KRT13"      "CKB"        "ETFB"       "FDFT1"      "TMPRSS3"   
#> [131] "FRK"        "ALDH3B1"    "CD44"       "MAPK13"     "PDLIM3"    
#> [136] "KLK11"      "TH"         "PAPSS2"     "TIAM1"      "PKP3"      
#> [141] "GINS2"      "LTF"        "FARP1"      "FKBP5"      "ADD3"      
#> [146] "HSPA4L"     "SCARB1"     "CLIC3"      "TPD52L1"    "XRCC3"     
#> [151] "ID2"        "SLC26A2"    "PCP4"       "SGK1"       "MDK"       
#> [156] "TOP2A"      "SFN"        "HPRT1"      "PPIF"       "CYP26B1"   
#> [161] "LAMC2"      "ISG20"      "RNASEH2A"   "CDC6"       "TFAP2C"    
#> [166] "CCNA1"      "IDH2"       "DUSP2"      "PLK4"       "NXT1"      
#> [171] "KLK10"      "DHCR7"      "FABP5"      "SLC16A1"    "OPN3"      
#> [176] "KIF20A"     "ASS1"       "PERP"       "TRIM29"     "NMU"       
#> [181] "HR"         "S100A9"     "LSR"        "STIL"       "ST14"      
#> [186] "IMPA2"      "SLC7A5"     "GJB3"       "KCNK5"      "BTG3"      
#> [191] "GAL"        "CDC20"      "FOXC1"      "CCN5"       "EEIG1"     
#> [196] "GFUS"       "GPER1"      "LARGE1"     "NHERF1"     "PLAAT3"    
#> 
#> $le$HALLMARK_MTORC1_SIGNALING
#>   [1] "XBP1"     "BHLHE40"  "SYTL2"    "BTG2"     "ELOVL5"   "SORD"    
#>   [7] "QDPR"     "SLC1A4"   "TM7SF2"   "DHCR24"   "FDXR"     "TBK1"    
#>  [13] "IGFBP5"   "EDEM1"    "ADIPOR2"  "CCNG1"    "ALDOA"    "HSPA9"   
#>  [19] "HMGCR"    "LTA4H"    "ACACA"    "ATP6V1D"  "SQSTM1"   "HSPA4"   
#>  [25] "STC1"     "NMT1"     "CDKN1A"   "SCD"      "TRIB3"    "CANX"    
#>  [31] "GSR"      "ACLY"     "UFM1"     "PSMC6"    "GLA"      "SLC6A6"  
#>  [37] "CYP51A1"  "GGA2"     "ACSL3"    "UBE2D3"   "RDH11"    "CD9"     
#>  [43] "ELOVL6"   "PIK3R3"   "GCLC"     "ETF1"     "CTH"      "HK2"     
#>  [49] "INSIG1"   "TCEA1"    "SERP1"    "DHFR"     "FGL2"     "GSK3B"   
#>  [55] "GTF2H1"   "FKBP2"    "NUPR1"    "USO1"     "PSME3"    "GBE1"    
#>  [61] "TMEM97"   "TUBG1"    "LDLR"     "EGLN3"    "COPS5"    "RIT1"    
#>  [67] "ATP2A2"   "GLRX"     "SEC11A"   "P4HA1"    "SKAP2"    "M6PR"    
#>  [73] "DDIT3"    "CFP"      "MTHFD2L"  "SDF2L1"   "SQLE"     "SLA"     
#>  [79] "ADD3"     "LDHA"     "HSPA5"    "SLC37A4"  "IDH1"     "HMGCS1"  
#>  [85] "FADS1"    "SLC1A5"   "PSMA3"    "PSMD13"   "CCNF"     "PPP1R15A"
#>  [91] "CYB5B"    "SLC2A3"   "PSMD12"   "RRP9"     "SLC7A11"  "GOT1"    
#>  [97] "ITGB2"    "NUFIP1"   "PSMB5"    "EEF1E1"   "IMMT"     "G6PD"    
#> [103] "FADS2"    "CORO1A"   "RPA1"     "PSMC4"    "PRDX1"    "PPIA"    
#> [109] "TXNRD1"   "BCAT1"    "UNG"      "IDI1"     "PLOD2"    "TFRC"    
#> [115] "IFI30"    "HPRT1"    "STIP1"    "PFKL"     "HSPE1"    "SSR1"    
#> [121] "RAB1A"    "CXCR4"    "YKT6"     "NFKBIB"   "CACYBP"   "LGMN"    
#> [127] "MLLT11"   "HSP90B1"  "HSPD1"    "NFYC"     "SLC2A1"   "DAPP1"   
#> [133] "PSMA4"    "STARD4"   "ME1"      "ABCF2"    "PGK1"     "EBP"     
#> [139] "PPA1"     "PSPH"     "HMBS"     "EIF2S2"   "ACTR2"    "DDIT4"   
#> [145] "ARPC5L"   "RRM2"     "PSMD14"   "PSMC2"    "DHCR7"    "AURKA"   
#> [151] "MCM4"     "SERPINH1" "GPI"      "MCM2"     "NAMPT"    "NUP205"  
#> [157] "RPN1"     "MAP2K3"   "PDAP1"    "PITPNB"   "SHMT2"    "UCHL5"   
#> [163] "TOMM40"   "CCT6A"    "PNO1"     "MTHFD2"   "GAPDH"    "TUBA4A"  
#> [169] "TPI1"     "CALR"     "PNP"      "PGM1"     "CTSC"     "POLR3G"  
#> [175] "VLDLR"    "TES"      "PLK1"     "BUB1"     "NFIL3"    "ASNS"    
#> [181] "CDC25A"   "ACTR3"    "PDK1"     "PSMG1"    "SLC7A5"   "GMPS"    
#> [187] "PHGDH"    "ENO1"     "SRD5A1"   "IFRD1"    "PSAT1"    "AK4"     
#> [193] "ATP5MC1"  "DDX39A"   "EPRS1"    "ERO1A"    "NHERF1"   "NIBAN1"  
#> [199] "SC5D"     "WARS1"   
#> 
#> $le$HALLMARK_INTERFERON_GAMMA_RESPONSE
#>   [1] "IFITM2"   "EIF4E3"   "RAPGEF6"  "CFB"      "PSME1"    "MVP"     
#>   [7] "PTPN1"    "DHX58"    "ISOC1"    "CASP7"    "IFI35"    "TRAFD1"  
#>  [13] "TOR1B"    "RBCK1"    "TDRD7"    "ZNFX1"    "TXNIP"    "BST2"    
#>  [19] "SP110"    "IRF7"     "RNF31"    "IRF9"     "RNF213"   "CDKN1A"  
#>  [25] "SPPL2A"   "DDX60"    "IRF2"     "NCOA3"    "TNFSF10"  "SELP"    
#>  [31] "LGALS3BP" "RIPK1"    "ST3GAL5"  "TNFAIP6"  "IFIT1"    "STAT3"   
#>  [37] "SRI"      "UBE2L6"   "OGFR"     "CFH"      "OAS3"     "FGL2"    
#>  [43] "PTPN6"    "SSPN"     "PARP14"   "NFKB1"    "ADAR"     "AUTS2"   
#>  [49] "CASP8"    "JAK2"     "CMKLR1"   "IFITM3"   "IL4R"     "NFKBIA"  
#>  [55] "TRIM14"   "STAT2"    "TRIM21"   "ISG15"    "IFI27"    "PSME2"   
#>  [61] "OAS2"     "CD69"     "ARID5B"   "SAMD9L"   "VAMP8"    "IFIT3"   
#>  [67] "SAMHD1"   "SLC25A28" "IFIT2"    "LAP3"     "LATS2"    "BTG1"    
#>  [73] "CD74"     "TRIM25"   "IRF5"     "P2RY14"   "FPR1"     "ST8SIA4" 
#>  [79] "HLA-DRB1" "PSMB8"    "XAF1"     "VAMP5"    "HERC6"    "MX1"     
#>  [85] "IRF8"     "B2M"      "IRF1"     "OASL"     "BPGM"     "SECTM1"  
#>  [91] "HLA-DQA1" "SERPING1" "HLA-DMA"  "IFI44L"   "RSAD2"    "FCGR1A"  
#>  [97] "IL10RA"   "PSMB10"   "KLRK1"    "GBP4"     "LCP2"     "PARP12"  
#> [103] "ARL4A"    "BANK1"    "STAT1"    "PSMA3"    "CD86"     "SOCS3"   
#> [109] "CMPK2"    "MT2A"     "STAT4"    "GCH1"     "IL2RB"    "TNFAIP2" 
#> [115] "CIITA"    "GBP6"     "HLA-G"    "GZMA"     "C1S"      "APOL6"   
#> [121] "TAPBP"    "RTP4"     "CD274"    "LYSMD2"   "CASP3"    "XCL1"    
#> [127] "EPSTI1"   "USP18"    "GPR18"    "IL15"     "BATF2"    "IL6"     
#> [133] "CASP4"    "CASP1"    "C1R"      "PSMA2"    "CSF2RB"   "FAS"     
#> [139] "EIF2AK2"  "PSMB9"    "NOD1"     "ITGB7"    "VCAM1"    "IL7"     
#> [145] "HLA-B"    "CXCL9"    "ZBP1"     "IFI44"    "IFI30"    "IFIH1"   
#> [151] "CD40"     "LY6E"     "ISG20"    "SLAMF7"   "HLA-A"    "PDE4B"   
#> [157] "MYD88"    "CCL2"     "IL18BP"   "MX2"      "NLRC5"    "TNFAIP3" 
#> [163] "CXCL11"   "CCL5"     "PLSCR1"   "TRIM26"   "TAP1"     "IRF4"    
#> [169] "HIF1A"    "PTGS2"    "IL15RA"   "ICAM1"    "PNPT1"    "CXCL10"  
#> [175] "SOCS1"    "CD38"     "NAMPT"    "CCL7"     "IDO1"     "NMI"     
#> [181] "PML"      "MTHFD2"   "PNP"      "PLA2G4A"  "PIM1"     "SOD2"    
#> [187] "PTPN2"    "IFNAR2"   "PFKP"     "RIPK2"    "UPP1"     "PELI1"   
#> [193] "NUP93"    "PSMB2"    "CMTR1"    "HELZ2"    "MARCHF1"  "RIGI"    
#> [199] "TMT1B"    "WARS1"   
#> 
#> $le$HALLMARK_INFLAMMATORY_RESPONSE
#>   [1] "SLC7A2"   "BTG2"     "TPBG"     "ACVR1B"   "P2RX4"    "HPN"     
#>   [7] "SLC1A2"   "SCN1B"    "FFAR2"    "PSEN1"    "P2RY2"    "NPFFR2"  
#>  [13] "LPAR1"    "TLR3"     "GPC3"     "CSF3R"    "BST2"     "AXL"     
#>  [19] "IRF7"     "TACR3"    "NDP"      "PCDH7"    "CALCRL"   "ATP2B1"  
#>  [25] "CDKN1A"   "PTPRE"    "PTGER2"   "TNFSF10"  "IFITM1"   "SERPINE1"
#>  [31] "ABCA1"    "BDKRB1"   "VIP"      "MSR1"     "TNFAIP6"  "ITGB3"   
#>  [37] "SRI"      "APLNR"    "INHBA"    "AHR"      "OLR1"     "IFNAR1"  
#>  [43] "STAB1"    "C3AR1"    "SLC4A4"   "TIMP1"    "P2RX7"    "GABBR1"  
#>  [49] "NLRP3"    "CSF1"     "F3"       "NFKB1"    "BEST1"    "SLC11A2" 
#>  [55] "PDPN"     "SCARF1"   "RGS16"    "CMKLR1"   "IL1R1"    "SLC31A2" 
#>  [61] "RGS1"     "IL4R"     "NFKBIA"   "LDLR"     "RASGRP1"  "ITGA5"   
#>  [67] "CSF3"     "CCRL2"    "PIK3R5"   "CD69"     "OSMR"     "MXD1"    
#>  [73] "HAS2"     "C5AR1"    "MMP14"    "RNF144B"  "ATP2A2"   "KCNMB2"  
#>  [79] "EMP3"     "CCL22"    "CLEC5A"   "GPR132"   "PTGIR"    "CD55"    
#>  [85] "PTAFR"    "EBI3"     "FPR1"     "ACVR2A"   "GP1BA"    "KLF6"    
#>  [91] "SGMS2"    "IRF1"     "NOD2"     "CYBB"     "SELL"     "SELE"    
#>  [97] "IL10RA"   "EREG"     "RELA"     "PTGER4"   "IL10"     "MEP1A"   
#> [103] "LCP2"     "HRH1"     "CD48"     "CCR7"     "IL7R"     "IL12B"   
#> [109] "IL1B"     "KCNA3"    "MEFV"     "IL18"     "GCH1"     "IL2RB"   
#> [115] "RAF1"     "GNA15"    "GPR183"   "MYC"      "HBEGF"    "IRAK2"   
#> [121] "TAPBP"    "RTP4"     "FZD5"     "TACR1"    "TNFRSF1B" "IL15"    
#> [127] "IL6"      "OSM"      "TNFSF15"  "NMUR1"    "CCL17"    "EIF2AK2" 
#> [133] "RHOG"     "IFNGR2"   "ADRM1"    "CXCL9"    "SLAMF1"   "CD14"    
#> [139] "CD70"     "LIF"      "TNFSF9"   "KCNJ2"    "CXCR6"    "CD40"    
#> [145] "SLC31A1"  "LY6E"     "DCBLD2"   "TLR1"     "PDE4B"    "CCL2"    
#> [151] "GNAI3"    "IL18R1"   "LCK"      "TNFRSF9"  "ICOSLG"   "LTA"     
#> [157] "TLR2"     "IL18RAP"  "CXCL11"   "CCL5"     "EDN1"     "SEMA4D"  
#> [163] "PROK2"    "HIF1A"    "IL1A"     "ICAM4"    "KIF1B"    "AQP9"    
#> [169] "ITGB8"    "SLC7A1"   "IL15RA"   "ICAM1"    "SPHK1"    "CXCL10"  
#> [175] "ADM"      "NAMPT"    "ATP2C1"   "ABI1"     "PLAUR"    "CXCL6"   
#> [181] "CCL7"     "ROS1"     "NMI"      "CD82"     "ADORA2B"  "MET"     
#> [187] "CX3CL1"   "LYN"      "LAMP3"    "CCL20"    "MARCO"    "RIPK2"   
#> [193] "CHST2"    "PVR"      "OPRK1"    "ADGRE1"   "CCL24"    "CXCL8"   
#> [199] "SELENOS"  "SLC28A2" 
#> 
#> $le$HALLMARK_IL6_JAK_STAT3_SIGNALING
#>  [1] "IL6ST"     "ACVR1B"    "STAM2"     "PTPN1"     "IL13RA1"   "GRB2"     
#>  [7] "CSF3R"     "TGFB1"     "IRF9"      "CD36"      "LEPR"      "PDGFC"    
#> [13] "IL17RB"    "STAT3"     "ITGB3"     "CD9"       "IL9R"      "IFNAR1"   
#> [19] "CSF1"      "ACVRL1"    "JUN"       "IL1R1"     "ITGA4"     "IL4R"     
#> [25] "STAT2"     "PIK3R5"    "IL3RA"     "TYK2"      "OSMR"      "HMOX1"    
#> [31] "PTPN11"    "CD44"      "A2M"       "EBI3"      "CNTFR"     "PLA2G2A"  
#> [37] "TNFRSF1A"  "IRF1"      "MAP3K8"    "STAT1"     "IL1B"      "HAX1"     
#> [43] "SOCS3"     "CSF2RA"    "IL17RA"    "INHBE"     "CBL"       "TNFRSF1B" 
#> [49] "IL6"       "CSF2RB"    "LTB"       "FAS"       "IL12RB1"   "TNFRSF12A"
#> [55] "CCR1"      "IL7"       "IFNGR2"    "CXCL9"     "CXCL13"    "CD14"     
#> [61] "IL2RG"     "MYD88"     "IL18R1"    "TNF"       "TLR2"      "CXCL11"   
#> [67] "IL10RB"    "LTBR"      "IL15RA"    "CXCL10"    "SOCS1"     "CSF2"     
#> [73] "IFNGR1"    "CD38"      "IL2RA"     "CCL7"      "CXCL3"     "BAK1"     
#> [79] "IL1R2"     "PIM1"      "PTPN2"     "CXCL1"     "TNFRSF21"  "CRLF2"    
#> [85] "DNTT"      "PF4"       "REG1A"    
#> 
#> $le$HALLMARK_TNFA_SIGNALING_VIA_NFKB
#>   [1] "SLC16A6"  "RHOB"     "IL6ST"    "BHLHE40"  "BTG2"     "CCND1"   
#>   [7] "NINJ1"    "DUSP5"    "SMAD3"    "TIPARP"   "EGR3"     "NR4A2"   
#>  [13] "DUSP4"    "KLF2"     "IER3"     "SDC4"     "DUSP1"    "AREG"    
#>  [19] "KLF9"     "FOS"      "REL"      "PLK2"     "FOSB"     "SQSTM1"  
#>  [25] "DENND5A"  "PFKFB3"   "PER1"     "ATP2B1"   "CDKN1A"   "PTPRE"   
#>  [31] "GEM"      "TSC22D1"  "KLF4"     "EIF1"     "SERPINE1" "ABCA1"   
#>  [37] "KDM6B"    "KLF10"    "EGR1"     "JAG1"     "TGIF1"    "ZFP36"   
#>  [43] "GADD45B"  "PDLIM5"   "HES1"     "TNFAIP6"  "IRS2"     "JUNB"    
#>  [49] "BCL6"     "INHBA"    "OLR1"     "BCL3"     "CSF1"     "JUN"     
#>  [55] "F3"       "FOSL2"    "NFKB1"    "CLCF1"    "TRIB1"    "TNIP2"   
#>  [61] "NFE2L2"   "NR4A1"    "NFKBIA"   "EGR2"     "CEBPD"    "F2RL1"   
#>  [67] "LDLR"     "CCRL2"    "TNIP1"    "CD69"     "B4GALT1"  "MXD1"    
#>  [73] "NFAT5"    "ZC3H12A"  "CD44"     "SIK1"     "ZBTB10"   "IFIT2"   
#>  [79] "GADD45A"  "IER2"     "G0S2"     "BTG1"     "PMEPA1"   "STAT5A"  
#>  [85] "PLAU"     "MSC"      "EFNA1"    "KLF6"     "IRF1"     "FJX1"    
#>  [91] "TNC"      "GFPT2"    "LITAF"    "TUBB2A"   "DNAJB4"   "DRAM1"   
#>  [97] "NR4A3"    "RELA"     "PTGER4"   "IER5"     "MAP3K8"   "IL7R"    
#> [103] "IL12B"    "IL1B"     "PLEK"     "CCNL1"    "ATF3"     "SNN"     
#> [109] "SOCS3"    "MCL1"     "PPP1R15A" "IL18"     "GCH1"     "PHLDA1"  
#> [115] "TNFAIP2"  "ID2"      "SLC2A3"   "GPR183"   "MYC"      "HBEGF"   
#> [121] "SGK1"     "CD80"     "PHLDA2"   "TNFAIP8"  "IL6"      "LAMB3"   
#> [127] "IFNGR2"   "SERPINB2" "LIF"      "TNFSF9"   "BIRC3"    "IFIH1"   
#> [133] "CCL4"     "PDE4B"    "BIRC2"    "CXCL2"    "TRAF1"    "CCL2"    
#> [139] "MARCKS"   "SERPINB8" "NFKB2"    "TNFRSF9"  "ICOSLG"   "RNF19B"  
#> [145] "TNF"      "MAFF"     "TNFAIP3"  "TLR2"     "CXCL11"   "ETS2"    
#> [151] "CCL5"     "EHD1"     "PANX1"    "EDN1"     "SAT1"     "KYNU"    
#> [157] "SPSB1"    "TAP1"     "IL1A"     "DUSP2"    "CD83"     "SLC2A6"  
#> [163] "IL23A"    "TRIP10"   "TANK"     "BMP2"     "YRDC"     "B4GALT5" 
#> [169] "PTGS2"    "RELB"     "IL15RA"   "VEGFA"    "ICAM1"    "SPHK1"   
#> [175] "CXCL10"   "CSF2"     "NAMPT"    "MAP2K3"   "PLAUR"    "FUT4"    
#> [181] "CXCL6"    "CFLAR"    "CXCL3"    "RCAN1"    "NFKBIE"   "BCL2A1"  
#> [187] "PNRC1"    "PTX3"     "SOD2"     "CCL20"    "NFIL3"    "RIPK2"   
#> [193] "CXCL1"    "FOSL1"    "CEBPB"    "BTG3"     "ACKR3"    "CCN1"    
#> [199] "PLPP3"    "RIGI"    
#> 
#> $le$HALLMARK_MITOTIC_SPINDLE
#>   [1] "PREX1"    "NUMA1"    "ARFIP2"   "TAOK2"    "KIF3B"    "RASA1"   
#>   [7] "SHROOM2"  "FLNB"     "SHROOM1"  "RAPGEF6"  "ARHGEF12" "APC"     
#>  [13] "PCM1"     "ARHGEF3"  "RHOT2"    "NF1"      "CYTH2"    "ARFGEF1" 
#>  [19] "WASL"     "BCL2L11"  "TSC1"     "ABR"      "RICTOR"   "FGD6"    
#>  [25] "LATS1"    "EZR"      "RFC1"     "ITSN1"    "DYNLL2"   "TRIO"    
#>  [31] "CDC42"    "TUBGCP2"  "SPTAN1"   "ARHGAP29" "CLIP1"    "MYH10"   
#>  [37] "CD2AP"    "TUBGCP6"  "CSNK1D"   "TLK1"     "ARF6"     "PALLD"   
#>  [43] "CDC27"    "ABL1"     "CTTN"     "PXN"      "PPP4R2"   "PAFAH1B1"
#>  [49] "PCGF5"    "TUBGCP5"  "CNTROB"   "ARHGAP5"  "RABGAP1"  "PKD2"    
#>  [55] "PDLIM5"   "ALMS1"    "MID1IP1"  "ROCK1"    "ALS2"     "DST"     
#>  [61] "EPB41"    "NIN"      "CDC42EP2" "KIFAP3"   "ARL8A"    "MARK4"   
#>  [67] "STAU1"    "RAB3GAP1" "SSH2"     "CDK5RAP2" "HDAC6"    "TBCD"    
#>  [73] "SYNPO"    "KPTN"     "AKAP13"   "NEDD9"    "CEP250"   "RANBP9"  
#>  [79] "ARHGDIA"  "CDC42BPA" "KIF22"    "DOCK4"    "NOTCH2"   "RAPGEF5" 
#>  [85] "KIF5B"    "CDC42EP4" "BCAR1"    "CAPZB"    "OPHN1"    "ARHGEF11"
#>  [91] "TIAM1"    "LLGL1"    "KLC1"     "VCL"      "TUBD1"    "NET1"    
#>  [97] "DLG1"     "SMC3"     "DOCK2"    "FARP1"    "ATG4B"    "DYNC1H1" 
#> [103] "SMC1A"    "YWHAE"    "SORBS2"   "CEP192"   "ARHGAP4"  "MAP1S"   
#> [109] "HOOK3"    "ARAP3"    "KNTC1"    "GSN"      "WASF2"    "SPTBN1"  
#> [115] "RALBP1"   "ARHGEF7"  "PCNT"     "FGD4"     "UXT"      "MAP3K11" 
#> [121] "RASA2"    "TOP2A"    "CEP72"    "MAPRE1"   "MYH9"     "CKAP5"   
#> [127] "KIF20B"   "CLASP1"   "ARHGEF2"  "LRPPRC"   "RACGAP1"  "KATNB1"  
#> [133] "BRCA2"    "EPB41L2"  "ESPL1"    "ECT2"     "GEMIN4"   "SAC3D1"  
#> [139] "ARHGAP10" "KIF3C"    "RHOF"     "NEK2"     "TUBGCP3"  "WASF1"   
#> [145] "CLIP2"    "RASAL2"   "MARCKS"   "SASS6"    "FLNA"     "NUSAP1"  
#> [151] "PLEKHG2"  "LMNB1"    "BCR"      "CCDC88A"  "ARHGAP27" "SMC4"    
#> [157] "KIF1B"    "CDK1"     "BIN1"     "MYO9B"    "PRC1"     "ACTN4"   
#> [163] "STK38L"   "CEP57"    "KIF11"    "FBXO5"    "MYO1E"    "SUN2"    
#> [169] "KIF15"    "AURKA"    "CENPE"    "ABI1"     "KIF4A"    "KIF23"   
#> [175] "SOS1"     "CENPF"    "BIRC5"    "TUBA4A"   "INCENP"   "DLGAP5"  
#> [181] "NCK1"     "TPX2"     "PIF1"     "PLK1"     "BUB1"     "CCNB2"   
#> [187] "KATNA1"   "NDC80"    "ANLN"     "KIF2C"    "FSCN1"    "CDC42EP1"
#> [193] "MID1"     "TTK"      "NCK2"     "CEP131"   "CNTRL"    "CPAP"    
#> [199] "SEPTIN9" 
#> 
#> $le$HALLMARK_SPERMATOGENESIS
#>   [1] "PCSK4"   "IFT88"   "RAD17"   "STAM2"   "GSTM3"   "CLGN"    "NPHP1"  
#>   [8] "PHF7"    "TUBA3C"  "SPATA6"  "HSPA1L"  "HSPA2"   "TEKT2"   "JAM3"   
#>  [15] "ZC3H14"  "SIRT1"   "SLC12A2" "IP6K1"   "HOXB1"   "CHRM4"   "STRBP"  
#>  [22] "PHKG2"   "GRM8"    "SHE"     "IDE"     "TNNI3"   "COIL"    "PRKAR2A"
#>  [29] "GAD1"    "CCT6B"   "ACRBP"   "PGK2"    "PAPOLB"  "MAP7"    "GSG1"   
#>  [36] "PEBP1"   "ALOX15"  "GPR182"  "GMCL1"   "ACE"     "NPY5R"   "NEFH"   
#>  [43] "OAZ3"    "SCG5"    "PACRG"   "DPEP3"   "CAMK4"   "DDX25"   "PGS1"   
#>  [50] "BRAF"    "MTOR"    "YBX2"    "DMC1"    "PIAS2"   "CLVS1"   "ADCYAP1"
#>  [57] "NOS1"    "IL13RA2" "VDAC3"   "MAST2"   "MLF1"    "HSPA4L"  "SYCP1"  
#>  [64] "ARL4A"   "SCG3"    "TSN"     "PARP2"   "ELOVL3"  "LDHC"    "MLLT10" 
#>  [71] "TULP2"   "TCP11"   "TOPBP1"  "GAPDHS"  "TKTL1"   "DCC"     "NEK2"   
#>  [78] "CFTR"    "CLPB"    "CRISP2"  "CNIH2"   "GFI1"    "SNAP91"  "TALDO1" 
#>  [85] "CHFR"    "CCNA1"   "NF2"     "MTNR1A"  "AGFG1"   "CDK1"    "PCSK1N" 
#>  [92] "TLE4"    "CDKN3"   "AURKA"   "SLC2A5"  "RPL39L"  "ACRV1"   "POMC"   
#>  [99] "NCAPH"   "CSNK2A2" "EZH2"    "BUB1"    "CCNB2"   "DBF4"    "RFC4"   
#> [106] "KIF2C"   "PSMG1"   "TTK"     "DMRT1"   "IL12RB2" "LPIN1"   "ART3"   
#> [113] "ACTL7B"  "ADAD1"   "ADAM2"   "AKAP4"   "CST8"    "DDX4"    "DNAJB8" 
#> [120] "H1-6"    "HBZ"     "HTR5A"   "MEP1B"   "NAA11"   "ODF1"    "PDHA2"  
#> [127] "PRM2"    "SEPTIN4" "SPMAP2"  "TNP1"    "TNP2"    "TSSK2"   "ZC2HC1C"
#> [134] "ZNRF4"   "ZPBP"   
#> 
#> $le$HALLMARK_GLYCOLYSIS
#>   [1] "TFF3"     "FUT8"     "TPBG"     "STC2"     "COG2"     "CACNA1H" 
#>   [7] "CYB5A"    "ALG1"     "XYLT2"    "IDUA"     "SLC35A3"  "PFKFB1"  
#>  [13] "ANG"      "ARPP19"   "GLCE"     "POLR3K"   "MPI"      "FKBP4"   
#>  [19] "BIK"      "IER3"     "GPC1"     "GPC4"     "ALDOA"    "IL13RA1" 
#>  [25] "CITED2"   "RBCK1"    "ECD"      "GPC3"     "AGL"      "ENO2"    
#>  [31] "PAXIP1"   "PGAM2"    "B4GALT7"  "DCN"      "CHST1"    "GNPDA1"  
#>  [37] "STC1"     "COL5A1"   "GOT2"     "PYGL"     "NOL3"     "GMPPB"   
#>  [43] "VCAN"     "B3GALT6"  "QSOX1"    "ARTN"     "ALDH9A1"  "GALE"    
#>  [49] "IRS2"     "FAM162A"  "GCLC"     "PHKA2"    "PPP2CB"   "CTH"     
#>  [55] "CASP6"    "HK2"      "KIF2A"    "FBP2"     "AGRN"     "GUSB"    
#>  [61] "SDHC"     "CLDN3"    "MED24"    "MXI1"     "SOD1"     "EGLN3"   
#>  [67] "GFPT1"    "PGLS"     "COPB2"    "SRD5A3"   "AK3"      "B4GALT1" 
#>  [73] "SPAG4"    "TPST1"    "CD44"     "SDC2"     "GMPPA"    "ALDOB"   
#>  [79] "P4HA2"    "MERTK"    "GLRX"     "BPNT1"    "AKR1A1"   "P4HA1"   
#>  [85] "LHX9"     "CHPF"     "NDUFV3"   "PMM2"     "PRPS1"    "SLC25A10"
#>  [91] "DPYSL4"   "B3GAT3"   "ALDH7A1"  "ABCB6"    "CHPF2"    "PAM"     
#>  [97] "PGM2"     "NT5E"     "PC"       "LHPP"     "HDLBP"    "PKP2"    
#> [103] "ELF3"     "GNE"      "MIF"      "LDHA"     "GYS1"     "NDST3"   
#> [109] "HSPA5"    "ZNF292"   "SLC37A4"  "IDH1"     "HOMER1"   "HAX1"    
#> [115] "TGFBI"    "GYS2"     "CAPN5"    "GALK2"    "PPFIA4"   "RPE"     
#> [121] "SDC3"     "HS6ST2"   "B4GALT4"  "CHST12"   "GALK1"    "LDHC"    
#> [127] "SLC16A3"  "DLD"      "GOT1"     "EFNA3"    "ANKZF1"   "GAL3ST1" 
#> [133] "GPR87"    "SDC1"     "G6PD"     "MIOX"     "PSMC4"    "GAPDHS"  
#> [139] "KDELR3"   "CLN6"     "TKTL1"    "PPIA"     "TXN"      "PGAM1"   
#> [145] "PYGB"     "PLOD2"    "ISG20"    "IGFBP3"   "HS2ST1"   "PDK3"    
#> [151] "ANGPTL4"  "CXCR4"    "TALDO1"   "HMMR"     "NSDHL"    "ME1"     
#> [157] "PGK1"     "NANP"     "CLDN9"    "EXT2"     "MDH2"     "SAP30"   
#> [163] "CDK1"     "VEGFA"    "DDIT4"    "SOX9"     "CHST6"    "AURKA"   
#> [169] "B3GAT1"   "B4GALT2"  "EXT1"     "KIF20A"   "ME2"      "TPI1"    
#> [175] "ADORA2B"  "B3GNT3"   "MET"      "RRAGD"    "TGFA"     "VLDLR"   
#> [181] "SLC25A13" "DEPDC1"   "MDH1"     "CHST4"    "PFKP"     "CHST2"   
#> [187] "EGFR"     "STMN1"    "NASP"     "CENPA"    "ENO1"     "UGP2"    
#> [193] "LCT"      "PLOD1"    "DSC2"     "AK4"      "ERO1A"    "GFUS"    
#> [199] "PKM"      "RARS1"   
#> 
#> $le$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY
#>  [1] "HHEX"    "EGLN2"   "NQO1"    "GPX4"    "CAT"     "SRXN1"   "ATOX1"  
#>  [8] "GSR"     "PRDX2"   "TXNRD2"  "NDUFB4"  "JUNB"    "PDLIM1"  "GCLC"   
#> [15] "ERCC2"   "HMOX2"   "SOD1"    "GPX3"    "NDUFS2"  "LSP1"    "SBNO2"  
#> [22] "GLRX"    "OXSR1"   "MSRA"    "FES"     "GCLM"    "CDKN2D"  "MBP"    
#> [29] "IPCEF1"  "FTL"     "MGST1"   "G6PD"    "PRDX6"   "PRDX1"   "ABCC1"  
#> [36] "TXN"     "TXNRD1"  "NDUFA6"  "STK25"   "MPO"     "GLRX2"   "PRNP"   
#> [43] "SOD2"    "PFKP"    "PRDX4"   "LAMTOR5" "PTPA"    "SCAF4"   "SELENOS"
#> 
#> $le$HALLMARK_UNFOLDED_PROTEIN_RESPONSE
#>   [1] "XBP1"     "WFS1"     "STC2"     "SLC30A5"  "SLC1A4"   "SPCS1"   
#>   [7] "DNAJA4"   "EXOC2"    "EDEM1"    "CNOT2"    "BAG3"     "DCP2"    
#>  [13] "HSPA9"    "ERN1"     "PARN"     "DCP1A"    "EEF2"     "TSPYL2"  
#>  [19] "ATF6"     "CNOT6"    "DCTN1"    "RPS14"    "GOSR2"    "SPCS3"   
#>  [25] "DNAJB9"   "EIF4A2"   "IMP3"     "EIF4E"    "EIF2AK3"  "IFIT1"   
#>  [31] "SHC1"     "IGFBP1"   "NOP14"    "SERP1"    "NHP2"     "CNOT4"   
#>  [37] "DNAJC3"   "WIPI1"    "SEC31A"   "ARFGAP1"  "LSM1"     "ATP6V0D1"
#>  [43] "NFYB"     "HERPUD1"  "KIF5B"    "EXOSC10"  "EIF4A3"   "EXOSC1"  
#>  [49] "NPM1"     "SEC11A"   "FUS"      "ZBTB17"   "SDAD1"    "CXXC1"   
#>  [55] "EIF2S1"   "YIF1A"    "EXOSC9"   "TUBB2A"   "FKBP14"   "PAIP1"   
#>  [61] "HSPA5"    "EXOSC5"   "ATF3"     "EXOSC2"   "EIF4A1"   "EDC4"    
#>  [67] "RRP9"     "EXOSC4"   "NOLC1"    "ALDH18A1" "BANF1"    "LSM4"    
#>  [73] "TATDN2"   "KDELR3"   "GEMIN4"   "YWHAZ"    "POP4"     "NOP56"   
#>  [79] "SSR1"     "CCL2"     "SRPRB"    "PREB"     "EIF4G1"   "XPOT"    
#>  [85] "HSP90B1"  "DDX10"    "NFYA"     "VEGFA"    "DDIT4"    "HYOU1"   
#>  [91] "KHSRP"    "EIF4EBP1" "PDIA5"    "ATF4"     "MTHFD2"   "CALR"    
#>  [97] "CEBPG"    "CKS1B"    "ASNS"     "CHAC1"    "SLC7A5"   "CEBPB"   
#> [103] "PDIA6"    "DKC1"     "PSAT1"    "ERO1A"    "H2AX"     "IARS1"   
#> [109] "MTREX"    "NABP1"    "SKIC3"    "SRPRA"    "TARS1"   
#> 
#> $le$HALLMARK_BILE_ACID_METABOLISM
#>   [1] "NUDT12"   "AR"       "HSD17B4"  "PEX12"    "EFHC1"    "ABCA3"   
#>   [7] "SLC27A2"  "ABCA2"    "PEX11A"   "PEX19"    "ABCD3"    "SLC23A1" 
#>  [13] "CROT"     "PEX11G"   "PEX7"     "SERPINA6" "SULT2B1"  "DIO1"    
#>  [19] "RXRA"     "DHCR24"   "FDXR"     "NEDD4"    "HSD3B7"   "PECR"    
#>  [25] "RETSAT"   "GNMT"     "LIPE"     "AMACR"    "LONP2"    "EPHX2"   
#>  [31] "PAOX"     "ISOC1"    "BCAR3"    "ABCA8"    "ATXN1"    "HACL1"   
#>  [37] "ABCA9"    "PFKM"     "PEX1"     "APOA1"    "ABCA6"    "PXMP2"   
#>  [43] "SCP2"     "CYP46A1"  "CAT"      "RXRG"     "ABCA5"    "HSD17B6" 
#>  [49] "PEX16"    "ABCA1"    "HSD17B11" "ALDH9A1"  "ALDH1A1"  "TTR"     
#>  [55] "AGXT"     "ALDH8A1"  "SLC27A5"  "MLYCD"    "SLC29A1"  "PNPLA8"  
#>  [61] "HSD3B1"   "SOAT2"    "CH25H"    "ABCA4"    "NR3C2"    "NR0B2"   
#>  [67] "SOD1"     "DIO2"     "AKR1D1"   "GSTK1"    "PEX6"     "GNPAT"   
#>  [73] "PEX26"    "SLC23A2"  "PHYH"     "PIPOX"    "CYP7A1"   "ABCG8"   
#>  [79] "SULT1B1"  "ABCG4"    "IDH1"     "PRDX5"    "FADS1"    "ABCD2"   
#>  [85] "GCLM"     "ACSL1"    "ACSL5"    "BMP6"     "FADS2"    "CYP27A1" 
#>  [91] "CYP8B1"   "IDI1"     "RBP1"     "TFCP2L1"  "NR1I2"    "KLF1"    
#>  [97] "LCK"      "IDH2"     "SLCO1A2"  "PEX13"    "SLC35B2"  "ABCD1"   
#> [103] "AQP9"     "OPTN"     "NPC1"     "CYP7B1"   "CYP39A1"  "BBOX1"   
#> [109] "HAO1"     "GC"       "NR1H4"    "SLC67A1" 
#> 
#> $le$HALLMARK_WNT_BETA_CATENIN_SIGNALING
#>  [1] "FRAT1"  "HDAC11" "AXIN2"  "LEF1"   "PSEN2"  "NOTCH4" "NCOR2"  "MAML1" 
#>  [9] "NUMB"   "JAG1"   "HEY1"   "TP53"   "DLL1"   "KAT2A"  "HDAC5"  "CCND2" 
#> [17] "HEY2"   "FZD8"   "RBPJ"   "DKK4"   "NCSTN"  "AXIN1"  "CTNNB1" "JAG2"  
#> [25] "NKD1"   "FZD1"   "WNT1"   "DVL2"   "MYC"    "WNT5B"  "CSNK1E" "TCF7"  
#> [33] "CUL1"   "PTCH1"  "GNAI1"  "ADAM17" "DKK1"   "NOTCH1" "PPARD"  "HDAC2" 
#> [41] "SKP2"   "WNT6"  
#> 
#> $le$HALLMARK_HYPOXIA
#>   [1] "CA12"     "FBP1"     "BCL2"     "BHLHE40"  "TPBG"     "STC2"    
#>   [7] "SIAH2"    "CCNG2"    "TGFB3"    "PPP1R3C"  "SELENBP1" "MAP3K1"  
#>  [13] "SULT2B1"  "NEDD4L"   "STBD1"    "TIPARP"   "BNIP3L"   "ATP7A"   
#>  [19] "IER3"     "GPC1"     "AKAP12"   "SDC4"     "GPC4"     "DUSP1"   
#>  [25] "ALDOA"    "NDST2"    "CITED2"   "GPC3"     "TPST2"    "HK1"     
#>  [31] "CDKN1B"   "IDS"      "ENO2"     "FOS"      "TPD52"    "PGAM2"   
#>  [37] "NDST1"    "DCN"      "PCK1"     "KIF5A"    "PFKFB3"   "STC1"    
#>  [43] "COL5A1"   "CDKN1A"   "PDGFB"    "SLC6A6"   "GCK"      "B3GALT6" 
#>  [49] "RORA"     "SERPINE1" "ENO3"     "NR3C1"    "CAV1"     "ZFP36"   
#>  [55] "VHL"      "ILVBL"    "IGFBP1"   "IRS2"     "FAM162A"  "PYGM"    
#>  [61] "GAA"      "CASP6"    "HK2"      "AMPD3"    "BCAN"     "PKLR"    
#>  [67] "GRHPR"    "JUN"      "F3"       "BGN"      "FOSL2"    "WSB1"    
#>  [73] "INHA"     "DTNA"     "KLF7"     "SRPX"     "PGF"      "MXI1"    
#>  [79] "GBE1"     "HEXA"     "RBPJ"     "HMOX1"    "LOX"      "SDC2"    
#>  [85] "FOXO3"    "ALDOB"    "P4HA2"    "LXN"      "GLRX"     "B4GALNT2"
#>  [91] "P4HA1"    "SLC25A1"  "CDKN1C"   "MT1E"     "BTG1"     "TGM2"    
#>  [97] "DPYSL4"   "DDIT3"    "HAS1"     "PAM"      "EFNA1"    "KLF6"    
#> [103] "PGM2"     "HDLBP"    "MIF"      "ALDOC"    "LDHA"     "GYS1"    
#> [109] "HSPA5"    "ZNF292"   "NAGK"     "SLC37A4"  "SCARB1"   "KDM3A"   
#> [115] "PRDX5"    "TGFBI"    "ATF3"     "PPFIA4"   "ANXA2"    "PPP1R15A"
#> [121] "HOXB9"    "MT2A"     "SDC3"     "SLC2A3"   "GALK1"    "LDHC"    
#> [127] "EFNA3"    "MYH9"     "EDN2"     "ANKZF1"   "IL6"      "GAPDHS"  
#> [133] "KDELR3"   "TKTL1"    "PLAC8"    "S100A4"   "ETS1"     "PFKL"    
#> [139] "ISG20"    "IGFBP3"   "PDK3"     "ANGPTL4"  "ERRFI1"   "PHKG1"   
#> [145] "CXCR4"    "MAFF"     "TNFAIP3"  "JMJD6"    "SLC2A1"   "PPARGC1A"
#> [151] "HS3ST1"   "PGK1"     "KLHL24"   "SAP30"    "TMEM45A"  "VEGFA"   
#> [157] "CP"       "DDIT4"    "PRKCA"    "ADM"      "SLC2A5"   "GPI"     
#> [163] "PKP1"     "PLAUR"    "EXT1"     "NCAN"     "NDRG1"    "GAPDH"   
#> [169] "TPI1"     "ADORA2B"  "PLIN2"    "PGM1"     "RRAGD"    "PNRC1"   
#> [175] "VLDLR"    "CHST3"    "XPNPEP1"  "PIM1"     "GCNT2"    "TES"     
#> [181] "NFIL3"    "PFKP"     "CHST2"    "EGFR"     "PDK1"     "CSRP2"   
#> [187] "ENO1"     "UGP2"     "ACKR3"    "AK4"      "BRS3"     "CAVIN1"  
#> [193] "CAVIN3"   "CCN1"     "CCN2"     "CCN5"     "ERO1A"    "LALBA"   
#> [199] "LARGE1"   "NOCT"    
#> 
#> $le$HALLMARK_UV_RESPONSE_UP
#>   [1] "NAT1"     "RHOB"     "TMBIM6"   "IL6ST"    "BTG2"     "CYB5R1"  
#>   [7] "RET"      "SULT1A1"  "IGFBP2"   "SPR"      "PRKCD"    "DLG4"    
#>  [13] "OLFM1"    "HSPA2"    "PLCL1"    "FKBP4"    "MAOA"     "MAPK8IP2"
#>  [19] "PPT1"     "BCL2L11"  "FGF18"    "ALDOA"    "HTR7"     "ENO2"    
#>  [25] "FOS"      "EIF5"     "GRINA"    "MSX1"     "FOSB"     "ATP6V1C1"
#>  [31] "SQSTM1"   "TACR3"    "EPHX1"    "PPP1R2"   "SLC25A4"  "ABCB1"   
#>  [37] "POLG2"    "CDO1"     "APOM"     "TGFBRAP1" "ACAA1"    "CA2"     
#>  [43] "CYP1A1"   "JUNB"     "PTPRD"    "DGAT1"    "C4BPB"    "NPTXR"   
#>  [49] "HYAL2"    "TST"      "MRPL23"   "NR4A1"    "COL2A1"   "SLC6A12" 
#>  [55] "PRKACA"   "GRPEL1"   "UROD"     "NFKBIA"   "GPX3"     "STARD3"  
#>  [61] "NXF1"     "RASGRP1"  "HSPA13"   "HNRNPU"   "CLCN2"    "RRAD"    
#>  [67] "NTRK3"    "DNAJB1"   "BSG"      "HMOX1"    "MMP14"    "TCHH"    
#>  [73] "EIF2S3"   "AQP3"     "PDLIM3"   "CDKN1C"   "BTG1"     "FURIN"   
#>  [79] "RAB27A"   "KCNH2"    "SHOX2"    "IRF1"     "CDKN2B"   "NPTX2"   
#>  [85] "RXRB"     "MARK2"    "CLTB"     "AP2S1"    "FMO1"     "PPAT"    
#>  [91] "PARP2"    "ATF3"     "CREG1"    "GCH1"     "CYB5B"    "POLE3"   
#>  [97] "WIZ"      "CDK2"     "MGAT1"    "CASP3"    "DDX21"    "POLR2H"  
#> [103] "PSMC3"    "ATP6V1F"  "ALAS1"    "CCND3"    "IL6"      "E2F5"    
#> [109] "DNAJA1"   "ARRB2"    "TFRC"     "CHKA"     "PPIF"     "STK25"   
#> [115] "STIP1"    "CNP"      "CXCL2"    "LHX2"     "SLC6A8"   "YKT6"    
#> [121] "CDC34"    "HLA-F"    "SIGMAR1"  "TAP1"     "CCK"      "CDC5L"   
#> [127] "BMP2"     "ICAM1"    "PRPF3"    "NKX2-5"   "FEN1"     "RPN1"    
#> [133] "PDAP1"    "GGH"      "TUBA4A"   "TYRO3"    "EPCAM"    "BID"     
#> [139] "BAK1"     "CEBPG"    "GLS"      "SOD2"     "CHRNA5"   "LYN"     
#> [145] "RFC4"     "ASNS"     "KLHDC3"   "AMD1"     "BTG3"     "GAL"     
#> [151] "CCNE1"    "AGO2"     "CTSV"     "H2AX"     "NUP58"    "ONECUT1" 
#> [157] "SELENOW"  "TARS1"   
#> 
#> $le$HALLMARK_CHOLESTEROL_HOMEOSTASIS
#>  [1] "CPEB2"     "TP53INP1"  "SEMA3B"    "ABCA2"     "ALCAM"     "TM7SF2"   
#>  [7] "FASN"      "ATXN2"     "GSTM2"     "HSD17B7"   "HMGCR"     "SCD"      
#> [13] "TRIB3"     "PMVK"      "CYP51A1"   "CLU"       "JAG1"      "PPARG"    
#> [19] "AVPR1A"    "CD9"       "GPX8"      "ANXA13"    "ECH1"      "MAL2"     
#> [25] "GUSB"      "ANTXR2"    "ETHE1"     "FDFT1"     "TMEM97"    "LDLR"     
#> [31] "STX5"      "FBXO6"     "LPL"       "ACSS2"     "CTNNB1"    "ACTG1"    
#> [37] "MVK"       "LGALS3"    "ADH4"      "SREBF2"    "SQLE"      "ALDOC"    
#> [43] "HMGCS1"    "ATF3"      "LSS"       "GLDC"      "FADS2"     "TNFRSF12A"
#> [49] "IDI1"      "ACAT2"     "CHKA"      "PCYT2"     "ATF5"      "ANXA5"    
#> [55] "PDK3"      "ERRFI1"    "LGMN"      "PLSCR1"    "NSDHL"     "STARD4"   
#> [61] "EBP"       "GNAI1"     "DHCR7"     "MVD"       "CXCL16"    "PLAUR"    
#> [67] "FABP5"     "FDPS"      "S100A11"   "CBS"       "PNRC1"     "NFIL3"    
#> [73] "NIBAN1"    "SC5D"     
#> 
#> $le$HALLMARK_KRAS_SIGNALING_UP
#>   [1] "TSPAN1"   "TSPAN13"  "SEMA3B"   "PTCD2"    "KIF5C"    "CROT"    
#>   [7] "FUCA1"    "CAB39L"   "SCN1B"    "DNMBP"    "PLAT"     "CBR4"    
#>  [13] "MAP3K1"   "MTMR10"   "FBXO4"    "CPE"      "SPARCL1"  "GADD45G" 
#>  [19] "CFB"      "ANO1"     "SERPINA3" "ITGBL1"   "IGF2"     "AKAP12"  
#>  [25] "F13A1"    "CDADC1"   "ITGA2"    "NAP1L2"   "EPB41L3"  "MMP11"   
#>  [31] "GNG11"    "MMP10"    "VWA5A"    "TPH1"     "EVI5"     "ABCB1"   
#>  [37] "SPON1"    "JUP"      "CIDEA"    "TSPAN7"   "KLF4"     "ANXA10"  
#>  [43] "USP12"    "APOD"     "ALDH1A2"  "MAP7"     "CA2"      "AKT2"    
#>  [49] "FLT4"     "AVL9"     "PRKG2"    "SNAP25"   "PIGR"     "PRRX1"   
#>  [55] "INHBA"    "CBX8"     "PECAM1"   "PLVAP"    "MAFB"     "CFH"     
#>  [61] "ETV1"     "ACE"      "NIN"      "C3AR1"    "SDCCAG8"  "PTPRR"   
#>  [67] "RBP4"     "RBM4"     "SCG5"     "TOR1AIP2" "IL33"     "TRIB1"   
#>  [73] "DUSP6"    "RGS16"    "NRP1"     "CMKLR1"   "CCND2"    "NR0B2"   
#>  [79] "PEG3"     "ATG10"    "F2RL1"    "LAT2"     "ARG1"     "MMD"     
#>  [85] "ENG"      "TFPI"     "RELN"     "ZNF277"   "TMEM100"  "NGF"     
#>  [91] "RABGAP1L" "GLRX"     "G0S2"     "STRN"     "PLAU"     "PLEK2"   
#>  [97] "PSMB8"    "DOCK2"    "SPP1"     "ADAM8"    "TMEM176A" "EMP1"    
#> [103] "IRF8"     "LAPTM5"   "BPGM"     "CD37"     "GFPT2"    "TMEM176B"
#> [109] "CLEC4A"   "IL10RA"   "EREG"     "LY96"     "TLR8"     "HSD11B1" 
#> [115] "SCG3"     "IL7R"     "IL1B"     "PRDM1"    "FCER1G"   "BTC"     
#> [121] "MMP9"     "PPP1R15A" "CSF2RA"   "ID2"      "RETN"     "AMMECR1" 
#> [127] "GYPC"     "GPNMB"    "PCP4"     "HBEGF"    "LCP1"     "CBL"     
#> [133] "WNT7A"    "ITGB2"    "SATB1"    "PPBP"     "ANKH"     "TNFRSF1B"
#> [139] "IKZF1"    "TRIB2"    "MAP4K1"   "CTSS"     "SPRY2"    "EPHB2"   
#> [145] "HOXD11"   "LIF"      "GABRA3"   "BIRC3"    "ETS1"     "WDR33"   
#> [151] "HKDC1"    "PTBP2"    "IGFBP3"   "DCBLD2"   "IL2RG"    "TRAF1"   
#> [157] "USH1C"    "HDAC9"    "ANGPTL4"  "PDCD1LG2" "ZNF639"   "CXCR4"   
#> [163] "SNAP91"   "TNFAIP3"  "MYCN"     "ST6GAL1"  "ADAMDEC1" "FGF9"    
#> [169] "BMP2"     "YRDC"     "PCSK1N"   "PTGS2"    "IL1RL2"   "SOX9"    
#> [175] "GALNT3"   "CXCL10"   "ALDH1A3"  "CSF2"     "ETV5"     "ETV4"    
#> [181] "ADAM17"   "PLAUR"    "TNNT2"    "BTBD3"    "TMEM158"  "GPRC5B"  
#> [187] "CCL20"    "SLPI"     "MPZL2"    "MALL"     "KCNN4"    "ADGRA2"  
#> [193] "ADGRL4"   "CCSER2"   "CFHR2"    "ERO1A"    "GUCY1A1"  "H2BC3"   
#> [199] "NR1H4"    "PRELID3B"
#> 
#> $le$HALLMARK_NOTCH_SIGNALING
#>  [1] "SKP1"    "CCND1"   "LFNG"    "ARRB1"   "FBXW11"  "PSEN2"   "HEYL"   
#>  [8] "WNT5A"   "JAG1"    "HES1"    "DLL1"    "PSENEN"  "KAT2A"   "NOTCH2" 
#> [15] "NOTCH3"  "FZD1"    "WNT2"    "DTX1"    "APH1A"   "FZD5"    "RBX1"   
#> [22] "DTX2"    "SAP30"   "CUL1"    "FZD7"    "DTX4"    "TCF7L2"  "PRKCA"  
#> [29] "NOTCH1"  "MAML2"   "PPARD"   "ST3GAL6"
#> 
#> $le$HALLMARK_COMPLEMENT
#>   [1] "GATA3"    "F7"       "TMPRSS6"  "CD46"     "DUSP5"    "PLAT"    
#>   [7] "PRKCD"    "AKAP10"   "CTSO"     "ANG"      "USP8"     "CFB"     
#>  [13] "PSEN1"    "ZEB1"     "PCLO"     "CASP7"    "PRDM4"    "SERPINC1"
#>  [19] "GRB2"     "GMFB"     "LIPA"     "ZFPM2"    "CASP9"    "RNF4"    
#>  [25] "LTA4H"    "GPD2"     "LRP1"     "NOTCH4"   "DOCK10"   "IRF7"    
#>  [31] "BRPF3"    "F10"      "CD36"     "PDGFB"    "ATOX1"    "F8"      
#>  [37] "IRF2"     "HSPA1A"   "GNAI2"    "VCPIP1"   "PRCP"     "S100A13" 
#>  [43] "KCNIP2"   "SERPINE1" "DOCK9"    "FN1"      "SERPINA1" "CDH13"   
#>  [49] "CLU"      "STX4"     "MMP13"    "CA2"      "TFPI2"    "USP14"   
#>  [55] "PPP2CB"   "OLR1"     "TIMP2"    "CTSD"     "C4BPB"    "CFH"     
#>  [61] "PRSS36"   "TIMP1"    "SH2B3"    "KIF2A"    "F3"       "HNF4A"   
#>  [67] "CALM3"    "ITGAM"    "DUSP6"    "JAK2"     "ERAP2"    "GNG2"    
#>  [73] "PPP4C"    "CDA"      "USP16"    "RASGRP1"  "CALM1"    "PIK3CG"  
#>  [79] "KCNIP3"   "PIK3R5"   "DOCK4"    "SIRT6"    "GNB2"     "USP15"   
#>  [85] "MMP14"    "APOC1"    "KLKB1"    "DPP4"     "APOBEC3F" "LAP3"    
#>  [91] "RCE1"     "ADAM9"    "CTSH"     "ITIH1"    "CD55"     "CSRP1"   
#>  [97] "LGALS3"   "GZMK"     "RABIF"    "PDP1"     "LTF"      "DGKH"    
#> [103] "GP1BA"    "C2"       "APOBEC3G" "IRF1"     "CPM"      "SERPING1"
#> [109] "DYRK2"    "PHEX"     "C1QC"     "CDK5R1"   "CTSB"     "HSPA5"   
#> [115] "MT3"      "GNGT2"    "HPCAL4"   "LCP2"     "C3"       "ACTN2"   
#> [121] "SCG3"     "PIK3CA"   "C1QA"     "CD40LG"   "MMP8"     "PLEK"    
#> [127] "FCER1G"   "GCA"      "MMP15"    "WAS"      "RAF1"     "CASP10"  
#> [133] "CBLB"     "KLK1"     "GZMA"     "C1S"      "CASP3"    "SRC"     
#> [139] "COL4A2"   "IL6"      "FCN1"     "CASP4"    "CR1"      "CASP1"   
#> [145] "C1R"      "S100A12"  "PSMB9"    "RHOG"     "CTSS"     "SERPINB2"
#> [151] "PCSK9"    "CD59"     "ADRA2B"   "ANXA5"    "F5"       "GNAI3"   
#> [157] "FYN"      "LAMP2"    "LCK"      "MAFF"     "LGMN"     "TNFAIP3" 
#> [163] "CCL5"     "PLSCR1"   "EHD1"     "ME1"      "KYNU"     "CR2"     
#> [169] "PLA2G7"   "GNB4"     "SPOCK2"   "FDX1"     "CP"       "CASP5"   
#> [175] "PLAUR"    "GZMB"     "PFN1"     "DGKG"     "MMP12"    "PREP"    
#> [181] "PRSS3"    "PLA2G4A"  "CTSC"     "XPNPEP1"  "PIM1"     "S100A9"  
#> [187] "LYN"      "L3MBTL4"  "CXCL1"    "CEBPB"    "APOA4"    "C9"      
#> [193] "CPQ"      "CTSL"     "CTSV"     "F2"       "GP9"      "MSRB1"   
#> [199] "PLG"      "RBSN"    
#> 
#> $le$HALLMARK_IL2_STAT5_SIGNALING
#>   [1] "XBP1"      "RHOB"      "BCL2"      "ENPP1"     "IGF1R"     "BHLHE40"  
#>   [7] "LRIG1"     "CISH"      "AHNAK"     "P2RX4"     "SERPINB6"  "FAH"      
#>  [13] "RHOH"      "ALCAM"     "BCL2L1"    "SNX9"      "MUC1"      "IKZF4"    
#>  [19] "SPRED2"    "ECM1"      "F2RL2"     "GATA1"     "AMACR"     "SNX14"    
#>  [25] "BATF"      "SERPINC1"  "PTH1R"     "IL1RL1"    "RNH1"      "TWSG1"    
#>  [31] "SOCS2"     "GPX4"      "MYO1C"     "DENND5A"   "MAP6"      "SYNGR2"   
#>  [37] "SHE"       "BMPR2"     "MAPKAPK2"  "TNFRSF18"  "SH3BGRL2"  "NCOA3"    
#>  [43] "ABCB1"     "PTGER2"    "TNFSF10"   "ITIH5"     "SELP"      "RORA"     
#>  [49] "WLS"       "ITGAV"     "ENO3"      "PRAF2"     "COL6A1"    "CA2"      
#>  [55] "GADD45B"   "PLEC"      "NFKBIZ"    "TLR7"      "POU2F1"    "AHR"      
#>  [61] "HK2"       "IKZF2"     "SCN9A"     "HIPK2"     "FLT3LG"    "HOPX"     
#>  [67] "FGL2"      "CSF1"      "SPRY4"     "PRKCH"     "PENK"      "SYT11"    
#>  [73] "CYFIP1"    "SMPDL3A"   "IRF6"      "RGS16"     "NRP1"      "CCND2"    
#>  [79] "PHTF2"     "IFITM3"    "GPR65"     "IL4R"      "SWAP70"    "IL3RA"    
#>  [85] "CCR4"      "GALM"      "TNFSF11"   "CDC42SE2"  "MXD1"      "PTRH2"    
#>  [91] "CD44"      "RABGAP1L"  "GSTO1"     "P4HA1"     "CDKN1C"    "TIAM1"    
#>  [97] "TGM2"      "FURIN"     "CD81"      "SLC29A2"   "APLP1"     "SPP1"     
#> [103] "KLF6"      "EMP1"      "NT5E"      "IRF8"      "SELL"      "SLC39A8"  
#> [109] "IL10RA"    "AGER"      "GBP4"      "CAPN3"     "IL10"      "MAP3K8"   
#> [115] "ARL4A"     "DHRS3"     "CAPG"      "CD48"      "SLC1A5"    "EOMES"    
#> [121] "TTC39B"    "CD86"      "IGF2R"     "GABARAPL1" "GPR83"     "PHLDA1"   
#> [127] "IL2RB"     "AHCY"      "SLC2A3"    "HUWE1"     "MYC"       "CST7"     
#> [133] "CTSZ"      "CASP3"     "UMPS"      "TNFRSF1B"  "CCND3"     "TNFRSF4"  
#> [139] "CDCP1"     "IL13"      "LTB"       "ITGAE"     "LIF"       "LCLAT1"   
#> [145] "ADAM19"    "ST3GAL4"   "ANXA4"     "S100A1"    "TRAF1"     "CDC6"     
#> [151] "IL18R1"    "TNFRSF8"   "TNFRSF9"   "MAFF"      "CKAP4"     "PLSCR1"   
#> [157] "ICOS"      "PUS1"      "LRRC8C"    "COCH"      "IRF4"      "CTLA4"    
#> [163] "CD83"      "CD79B"     "ITGA6"     "BMP2"      "BATF3"     "PDCD2L"   
#> [169] "PTCH1"     "CXCL10"    "SOCS1"     "CSF2"      "MYO1E"     "IFNGR1"   
#> [175] "ETV4"      "IL2RA"     "PLAGL1"    "PRNP"      "NDRG1"     "PLIN2"    
#> [181] "PNP"       "RRAGD"     "IL1R2"     "PIM1"      "DCPS"      "NOP2"     
#> [187] "NFIL3"     "TNFRSF21"  "GLIPR2"    "NCS1"      "ODC1"      "UCK2"     
#> [193] "CCNE1"     "DRC1"      "EEF1AKMT1" "ETFBKMT"   "GUCY1B1"   "HYCC2"    
#> [199] "PLPP1"    
#> 
#> $le$HALLMARK_INTERFERON_ALPHA_RESPONSE
#>  [1] "IFITM2"   "ELF1"     "PSME1"    "DHX58"    "LPAR6"    "IFI35"   
#>  [7] "TRAFD1"   "TDRD7"    "TXNIP"    "UBA7"     "BST2"     "SP110"   
#> [13] "IRF7"     "RNF31"    "IRF9"     "DDX60"    "IRF2"     "PARP9"   
#> [19] "PROCR"    "IFITM1"   "LGALS3BP" "UBE2L6"   "OGFR"     "SAMD9"   
#> [25] "OAS1"     "CSF1"     "PARP14"   "ADAR"     "CASP8"    "IFITM3"  
#> [31] "TMEM140"  "IL4R"     "TRIM14"   "STAT2"    "TRIM21"   "ISG15"   
#> [37] "IFI27"    "CCRL2"    "PSME2"    "TRIM5"    "MOV10"    "SAMD9L"  
#> [43] "IFIT3"    "SLC25A28" "IFIT2"    "LAP3"     "CD74"     "NUB1"    
#> [49] "TRIM25"   "GMPR"     "PSMB8"    "HERC6"    "MX1"      "B2M"     
#> [55] "IRF1"     "OASL"     "SELL"     "IFI44L"   "GBP2"     "RSAD2"   
#> [61] "HLA-C"    "GBP4"     "PARP12"   "PSMA3"    "CMPK2"    "C1S"     
#> [67] "RTP4"     "EPSTI1"   "USP18"    "IL15"     "BATF2"    "CASP1"   
#> [73] "EIF2AK2"  "PSMB9"    "IL7"      "IFI44"    "IFI30"    "IFIH1"   
#> [79] "LY6E"     "ISG20"    "CNP"      "CD47"     "CXCL11"   "PLSCR1"  
#> [85] "TRIM26"   "TAP1"     "PNPT1"    "CXCL10"   "NMI"      "LAMP3"   
#> [91] "RIPK2"    "NCOA7"    "CMTR1"    "HELZ2"    "MVB12A"   "TENT5A"  
#> [97] "WARS1"   
#> 
#> $le$HALLMARK_PI3K_AKT_MTOR_SIGNALING
#>   [1] "VAV3"     "PRKAG1"   "PLA2G12A" "CAB39L"   "MAPK9"    "PTEN"    
#>   [7] "TBK1"     "GNA14"    "MKNK2"    "TSC2"     "ATF1"     "IRAK4"   
#>  [13] "DUSP3"    "GRB2"     "CDKN1B"   "CLTC"     "YWHAB"    "ACACA"   
#>  [19] "RALB"     "SQSTM1"   "ITPR2"    "CDKN1A"   "PRKAR2A"  "TRIB3"   
#>  [25] "PRKAA2"   "PIN1"     "RPTOR"    "RIPK1"    "THEM4"    "PIKFYVE" 
#>  [31] "EIF4E"    "CAB39"    "UBE2D3"   "AKT1"     "PIK3R3"   "AKT1S1"  
#>  [37] "ECSIT"    "MAPK10"   "PAK4"     "HRAS"     "UBE2N"    "GSK3B"   
#>  [43] "ARPC3"    "CAMK4"    "PLCB1"    "STAT2"    "ARHGDIA"  "MAPK8"   
#>  [49] "PITX2"    "PTPN11"   "GNGT1"    "RIT1"     "NGF"      "TRAF2"   
#>  [55] "PPP1CA"   "SMAD2"    "RAC1"     "MAPK1"    "TIAM1"    "DDIT3"   
#>  [61] "TNFRSF1A" "SLA"      "PRKCB"    "ADCY2"    "PLCG1"    "ARF1"    
#>  [67] "RAF1"     "MKNK1"    "CDK2"     "MAPKAP1"  "MAP2K6"   "SFN"     
#>  [73] "FGF17"    "FASLG"    "RPS6KA1"  "NOD1"     "PPP2R1B"  "MAP3K7"  
#>  [79] "CFL1"     "IL2RG"    "MYD88"    "CXCR4"    "NFKBIB"   "LCK"     
#>  [85] "CDK4"     "HSP90B1"  "SLC2A1"   "DAPP1"    "AP2M1"    "CDK1"    
#>  [91] "CSNK2B"   "E2F1"     "ACTR2"    "RPS6KA3"  "MAP2K3"   "PFN1"    
#>  [97] "CALR"     "NCK1"     "ACTR3"    "EGFR"     "PDK1"     "FGF22"   
#> [103] "FGF6"     "GRK2"     "IL4"     
#> 
#> $le$HALLMARK_UV_RESPONSE_DN
#>   [1] "APBB2"    "INPP4B"   "IGF1R"    "BHLHE40"  "AGGF1"    "IRS1"    
#>   [7] "SPOP"     "CDON"     "SMAD3"    "KCNMA1"   "RXRA"     "PTEN"    
#>  [13] "GJA1"     "AMPH"     "ZMIZ1"    "DBP"      "IGFBP5"   "MAGI2"   
#>  [19] "CAP2"     "ICA1"     "PDGFRB"   "SYNJ2"    "YTHDC1"   "MAP2K5"  
#>  [25] "DUSP1"    "PMP22"    "MMP16"    "MGMT"     "CITED2"   "LPAR1"   
#>  [31] "ATXN1"    "CDKN1B"   "RUNX1"    "ATRX"     "PIAS3"    "COL1A2"  
#>  [37] "RBPMS"    "COL1A1"   "PTPRM"    "DDAH1"    "COL3A1"   "SFMBT1"  
#>  [43] "MIOS"     "ATRN"     "MGLL"     "NIPBL"    "COL5A2"   "SYNE1"   
#>  [49] "PRKAR2B"  "TGFBR3"   "DLC1"     "BDNF"     "PRDM2"    "ATP2B1"  
#>  [55] "RGS4"     "PRKCE"    "SMAD7"    "TJP1"     "PHF3"     "SERPINE1"
#>  [61] "FZD2"     "VAV2"     "NR3C1"    "BCKDHB"   "NR1D2"    "PPARG"   
#>  [67] "CAV1"     "NEK7"     "PDLIM5"   "CDK13"    "ERBB2"    "ITGB3"   
#>  [73] "COL11A1"  "TGFBR2"   "DAB2"     "SRI"      "EFEMP1"   "PIK3R3"  
#>  [79] "FHL2"     "SCN8A"    "INSIG1"   "F3"       "ATP2B4"   "NFKB1"   
#>  [85] "PTPN21"   "NRP1"     "MRPS31"   "PEX14"    "LDLR"     "SNAI2"   
#>  [91] "CDC42BPA" "NOTCH2"   "SIPA1L1"  "TFPI"     "HAS2"     "SDC2"    
#>  [97] "DYRK1A"   "MT1E"     "GRK5"     "BMPR1A"   "GCNT1"    "DLG1"    
#> [103] "ACVR2A"   "LAMC1"    "MAP1B"    "WDR37"    "ID1"      "ADD3"    
#> [109] "AKT3"     "MAPK14"   "RASA2"    "ANXA2"    "CELF2"    "FBLN5"   
#> [115] "MYC"      "MTA1"     "ABCC1"    "KIT"      "KALRN"    "ANXA4"   
#> [121] "PTGFR"    "FYN"      "PIK3CD"   "CACNA1A"  "LTBP1"    "SLC7A1"  
#> [127] "PLCB4"    "PRKCA"    "RND3"     "ARHGEF9"  "ATP2C1"   "ADORA2B" 
#> [133] "MET"      "VLDLR"    "SCHIP1"   "NFIB"     "ADGRL2"   "CCN1"    
#> [139] "DMAC2L"   "PLPP3"    "SCAF8"    "SLC67A1"  "TENT4A"   "TOGARAM1"
#> 
#> $le$HALLMARK_PROTEIN_SECRETION
#>  [1] "KRT18"    "SCAMP1"   "COG2"     "ARFIP1"   "AP2B1"    "MON2"    
#>  [7] "TOM1L1"   "TSG101"   "AP3B1"    "SNX2"     "SNAP23"   "NAPA"    
#> [13] "GALC"     "ICA1"     "RAB14"    "ARFGEF1"  "ARFGEF2"  "GBF1"    
#> [19] "ATP7A"    "PPT1"     "COPB1"    "VAMP4"    "NAPG"     "CLN5"    
#> [25] "TMED10"   "CD63"     "VAMP3"    "TMX1"     "ATP6V1H"  "TPD52"   
#> [31] "CLTC"     "SEC24D"   "CLCN3"    "YIPF6"    "GOSR2"    "VPS4B"   
#> [37] "GLA"      "ATP6V1B1" "ABCA1"    "ADAM10"   "STX16"    "STX12"   
#> [43] "GOLGA4"   "RER1"     "AP1G1"    "RAB22A"   "SGMS1"    "OCRL"    
#> [49] "AP3S1"    "DST"      "RAB5A"    "TMED2"    "ERGIC3"   "SSPN"    
#> [55] "RAB2A"    "SEC31A"   "LMAN1"    "USO1"     "SOD1"     "SCRN1"   
#> [61] "SCAMP3"   "COPB2"    "SEC22B"   "VPS45"    "BET1"     "COPE"    
#> [67] "MAPK1"    "STX7"     "ARCN1"    "M6PR"     "PAM"      "VAMP7"   
#> [73] "BNIP3"    "SH3GL2"   "TSPAN8"   "AP2S1"    "IGF2R"    "ARF1"    
#> [79] "CAV2"     "DNM1L"    "STAM"     "GNAS"     "LAMP2"    "YKT6"    
#> [85] "RAB9A"    "AP2M1"    "ARFGAP3"  "KIF1B"    "CLTA"     "RPS6KA3" 
#> [91] "ATP1A1"   "ZW10"     "CTSC"     "ANP32E"   "EGFR"     "DOP1A"   
#> 
#> $le$HALLMARK_OXIDATIVE_PHOSPHORYLATION
#>   [1] "ACADSB"   "GLUD1"    "ALDH6A1"  "COX6C"    "CYB5A"    "MRPS30"  
#>   [7] "PRDX3"    "SURF1"    "NDUFS4"   "PDHB"     "NDUFA2"   "ETFDH"   
#>  [13] "COX15"    "RHOT2"    "CPT1A"    "ATP6AP1"  "RETSAT"   "UQCRQ"   
#>  [19] "RHOT1"    "ATP6V1G1" "NDUFA5"   "SLC25A20" "ISCU"     "CASP7"   
#>  [25] "PDK4"     "COX7C"    "HSPA9"    "SUCLA2"   "ATP6V1H"  "NDUFA7"  
#>  [31] "COX11"    "COX17"    "BCKDHA"   "GPX4"     "ATP6V1D"  "ATP6V1C1"
#>  [37] "ATP6V0C"  "GOT2"     "ATP6V0E1" "HADHB"    "MAOB"     "UQCR11"  
#>  [43] "SLC25A4"  "NDUFV2"   "ACADVL"   "MRPL35"   "UQCRC2"   "OXA1L"   
#>  [49] "COX6A1"   "ISCA1"    "NDUFC1"   "ACAT1"    "ACAA1"    "MTRR"    
#>  [55] "NNT"      "NDUFB1"   "ECHS1"    "NDUFS7"   "NDUFB4"   "NDUFC2"  
#>  [61] "TIMM9"    "SLC25A12" "CS"       "NDUFA3"   "MFN2"     "TIMM13"  
#>  [67] "ATP6V1E1" "IDH3A"    "TCIRG1"   "ECH1"     "ABCB7"    "DLST"    
#>  [73] "VDAC1"    "PDHX"     "ATP1B1"   "NDUFS1"   "NDUFS8"   "MTRF1"   
#>  [79] "BDH2"     "COX8A"    "ACADM"    "SDHC"     "ETFB"     "GRPEL1"  
#>  [85] "MTX2"     "SLC25A11" "NDUFS2"   "OGDH"     "ATP6V0B"  "NDUFB8"  
#>  [91] "PMPCA"    "SLC25A6"  "NDUFA8"   "NDUFAB1"  "MRPS11"   "IDH3B"   
#>  [97] "MRPL34"   "UQCRC1"   "NDUFB7"   "HADHA"    "NDUFB6"   "COX7A2"  
#> [103] "ACO2"     "PDP1"     "PHYH"     "AIFM1"    "NDUFB5"   "NQO2"    
#> [109] "CYC1"     "COX7A2L"  "NDUFB2"   "TIMM10"   "SLC25A3"  "NDUFS6"  
#> [115] "TIMM17A"  "VDAC3"    "LDHA"     "DECR1"    "OAT"      "VDAC2"   
#> [121] "NDUFA4"   "SDHA"     "IDH1"     "UQCRB"    "NDUFS3"   "UQCR10"  
#> [127] "NDUFA1"   "SDHB"     "MGST3"    "NDUFB3"   "DLD"      "OPA1"    
#> [133] "ETFA"     "ATP6V1F"  "ALAS1"    "LRPPRC"   "COX6B1"   "COX5B"   
#> [139] "IMMT"     "BAX"      "NDUFV1"   "ACAA2"    "HSD17B10" "COX10"   
#> [145] "HTRA2"    "NDUFA6"   "SDHD"     "SUCLG1"   "CYB5R3"   "FXN"     
#> [151] "FH"       "COX7B"    "IDH3G"    "PHB2"     "DLAT"     "MRPS15"  
#> [157] "IDH2"     "MRPS12"   "AFG3L2"   "COX5A"    "COX4I1"   "TIMM8B"  
#> [163] "UQCRFS1"  "HCCS"     "TIMM50"   "MDH2"     "MRPS22"   "FDX1"    
#> [169] "CYCS"     "POR"      "GPI"      "SUPV3L1"  "SLC25A5"  "MRPL11"  
#> [175] "MRPL15"   "TOMM22"   "NDUFA9"   "POLR2F"   "MDH1"     "LDHB"    
#> [181] "PDHA1"    "UQCRH"    "ATP5F1A"  "ATP5F1B"  "ATP5F1C"  "ATP5F1D" 
#> [187] "ATP5F1E"  "ATP5MC1"  "ATP5MC2"  "ATP5MC3"  "ATP5ME"   "ATP5MF"  
#> [193] "ATP5MG"   "ATP5PB"   "ATP5PD"   "ATP5PF"   "ATP5PO"   "ECI1"    
#> [199] "MPC1"     "TOMM70"  
#> 
#> $le$HALLMARK_DNA_REPAIR
#>   [1] "DCTN4"   "BCAM"    "GMPR2"   "ADCY6"   "XPC"     "NME3"    "POLD4"  
#>   [8] "SMAD5"   "TAF9"    "TSG101"  "SURF1"   "AAAS"    "POLL"    "CANT1"  
#>  [15] "VPS37D"  "RRM2B"   "TK2"     "DDB2"    "IMPDH2"  "CCNO"    "POLR2A" 
#>  [22] "GTF2F1"  "ERCC4"   "MPG"     "NUDT9"   "COX17"   "TARBP2"  "SUPT4H1"
#>  [29] "GPX4"    "GTF2H3"  "VPS28"   "NME4"    "ERCC8"   "CETN2"   "NT5C"   
#>  [36] "DAD1"    "POLR2E"  "EIF1B"   "AK1"     "ARL6IP1" "EDF1"    "TP53"   
#>  [43] "POLR3GL" "BRF2"    "POLB"    "GTF2H5"  "RPA3"    "RAD52"   "ERCC2"  
#>  [50] "PRIM1"   "SNAPC5"  "MRPL40"  "ERCC1"   "TAF10"   "USP11"   "RAE1"   
#>  [57] "DUT"     "TMED2"   "NFX1"    "ERCC5"   "NPR2"    "DDB1"    "GTF2H1" 
#>  [64] "SNAPC4"  "ZNF707"  "CDA"     "POLR2G"  "NCBP2"   "AK3"     "ELL"    
#>  [71] "RFC5"    "TAF12"   "POLA1"   "LIG1"    "REV3L"   "NME1"    "POLR1D" 
#>  [78] "POLR2K"  "POLE4"   "RNMT"    "BOLA2"   "ERCC3"   "POLR2I"  "PDE6G"  
#>  [85] "STX3"    "BCAP31"  "HCLS1"   "ITPA"    "POLR2J"  "RALA"    "RPA2"   
#>  [92] "APRT"    "POLH"    "GTF2B"   "CMPK2"   "DGCR8"   "CSTF3"   "GTF2A2" 
#>  [99] "UMPS"    "POLR2H"  "GUK1"    "NUDT21"  "POLR2C"  "TAF13"   "ZWINT"  
#> [106] "POLR3C"  "POM121"  "RBX1"    "ADRM1"   "SAC3D1"  "SDCBP"   "TAF6"   
#> [113] "TAF1C"   "HPRT1"   "PDE4B"   "POLA2"   "POLD3"   "PCNA"    "CLP1"   
#> [120] "SEC61A1" "RFC3"    "RAD51"   "ADA"     "POLD1"   "GTF3C5"  "SUPT5H" 
#> [127] "DGUOK"   "FEN1"    "TYMS"    "PNP"     "POLR1C"  "SSRP1"   "POLR2D" 
#> [134] "VPS37B"  "SF3A3"   "UPF3B"   "POLR2F"  "RFC2"    "RFC4"    "AGO4"   
#> [141] "ALYREF"  "ELOA"    "GSDME"   "MPC2"    "NELFB"   "NELFCD"  "NELFE"  
#> [148] "NT5C3A"  "POLR1H"  "SRSF6"  
#> 
#> $le$HALLMARK_KRAS_SIGNALING_DN
#>   [1] "GAMT"      "SIDT1"     "GP2"       "BTG2"      "RGS11"     "CPEB3"    
#>   [7] "CCDC106"   "CACNG1"    "SLC29A3"   "LFNG"      "CPB1"      "BMPR1B"   
#>  [13] "CAPN9"     "IGFBP2"    "GPRC5C"    "CACNA1F"   "THRB"      "CELSR2"   
#>  [19] "PDK2"      "FGFR3"     "IDUA"      "C5"        "MYO15A"    "PDE6B"    
#>  [25] "BRDT"      "ZBTB16"    "NR4A2"     "EFHD1"     "HNF1A"     "MTHFR"    
#>  [31] "NRIP2"     "ADRA2C"    "SPHK2"     "TAS2R4"    "FGF16"     "P2RY4"    
#>  [37] "COPZ2"     "MAST3"     "EGF"       "EPHA5"     "TNNI3"     "ARPP21"   
#>  [43] "SERPINA10" "ITIH3"     "LYPD3"     "THNSL2"    "TFAP2B"    "SLC25A23" 
#>  [49] "KRT1"      "ASB7"      "ATP6V1B1"  "FGGY"      "SKIL"      "MYOT"     
#>  [55] "HTR1D"     "P2RX6"     "PTPRJ"     "MFSD6"     "GDNF"      "PCDHB1"   
#>  [61] "PNMT"      "MAGIX"     "ABCB11"    "PRODH"     "NTF3"      "NR6A1"    
#>  [67] "SYNPO"     "KRT13"     "OXT"       "HSD11B2"   "COL2A1"    "HTR1B"    
#>  [73] "ENTPD7"    "ARHGDIG"   "RIBC2"     "NPHS1"     "KCNQ2"     "SLC6A3"   
#>  [79] "BARD1"     "SLC5A5"    "ACTC1"     "TG"        "PAX3"      "SPRR3"    
#>  [85] "NGB"       "CLSTN3"    "YBX2"      "SCGB1A1"   "CD207"     "SLC16A7"  
#>  [91] "DLK2"      "AKR1B10"   "ALOX12B"   "KCND1"     "CNTFR"     "GP1BA"    
#>  [97] "IRS4"      "NOS1"      "SHOX2"     "MX1"       "SNCB"      "MSH5"     
#> [103] "UPK3B"     "IFI44L"    "RSAD2"     "ABCG4"     "WNT16"     "YPEL1"    
#> [109] "RYR2"      "CD40LG"    "IL12B"     "PLAG1"     "KCNE2"     "MYH7"     
#> [115] "SNN"       "MEFV"      "TGFB2"     "KLHDC8A"   "SLC12A3"   "CCR8"     
#> [121] "LGALS7"    "NUDT11"    "SGK1"      "CD80"      "EDN2"      "TFF2"     
#> [127] "KCNN1"     "KRT15"     "ITGB1BP2"  "KCNMB1"    "CAMK1D"    "CLDN8"    
#> [133] "DCC"       "SERPINB2"  "TCL1A"     "CLPS"      "TFCP2L1"   "PDCD1"    
#> [139] "PTGFR"     "KRT5"      "CKM"       "CCNA1"     "SPTBN2"    "EDN1"     
#> [145] "IFNG"      "CALML5"    "SLC38A3"   "EDAR"      "RYR1"      "KRT4"     
#> [151] "SLC30A3"   "GTF3C5"    "PKP1"      "SOX10"     "STAG3"     "KLK7"     
#> [157] "CALCB"     "CYP39A1"   "CPA2"      "SMPX"      "GPR19"     "TEX15"    
#> [163] "DTNB"      "KLK8"      "CDKAL1"    "TGM1"      "TLX1"      "GPR3"     
#> [169] "CHST2"     "CLDN16"    "TCF7L1"    "SLC6A14"   "AMBN"      "ATP4A"    
#> [175] "CDH16"     "CHRNG"     "COQ8A"     "CYP11B2"   "FGF22"     "FSHB"     
#> [181] "GRID2"     "IL5"       "INSL5"     "KMT2D"     "MACROH2A2" "NPY4R"    
#> [187] "PAX4"      "PRKN"      "PROP1"     "SCN10A"    "SELENOP"   "SSTR4"    
#> [193] "TENM2"     "TENT5C"    "TSHB"      "UGT2B17"   "VPREB1"    "VPS50"    
#> [199] "ZC2HC1C"   "ZNF112"   
#> 
#> $le$HALLMARK_PEROXISOME
#>   [1] "ABCC8"    "HSD17B4"  "ELOVL5"   "SEMA3C"   "HMGCL"    "SLC27A2" 
#>   [7] "PEX11A"   "ABCD3"    "SERPINA6" "SULT2B1"  "DIO1"     "DLG4"    
#>  [13] "DHCR24"   "CDK7"     "CRAT"     "MVP"      "ABCC5"    "HSD3B7"  
#>  [19] "RETSAT"   "LONP2"    "CRABP2"   "EPHX2"    "ABCB4"    "ISOC1"   
#>  [25] "ATXN1"    "SCP2"     "FIS1"     "PEX2"     "CAT"      "RXRG"    
#>  [31] "IDE"      "STS"      "SLC25A4"  "ABCB1"    "ALB"      "VPS4B"   
#>  [37] "PEX11B"   "ACAA1"    "HSD17B11" "ALDH9A1"  "RDH11"    "ALDH1A1" 
#>  [43] "TTR"      "HAO2"     "MLYCD"    "ABCB9"    "ERCC1"    "ECH1"    
#>  [49] "CLN8"     "ACOT8"    "HRAS"     "SLC25A17" "HSD11B2"  "CTBP1"   
#>  [55] "SOD1"     "PEX14"    "EHHADH"   "GSTK1"    "PEX6"     "ACOX1"   
#>  [61] "GNPAT"    "SCGB1A1"  "PEX5"     "BCL10"    "CADM1"    "SIAH1"   
#>  [67] "SLC23A2"  "CEL"      "ERCC3"    "SMARCC1"  "IDH1"     "PRDX5"   
#>  [73] "DHRS3"    "FADS1"    "ABCD2"    "CNBP"     "TOP2A"    "ACSL1"   
#>  [79] "FABP6"    "ACSL5"    "PRDX1"    "CLN6"     "CACNA1B"  "IDI1"    
#>  [85] "NR1I2"    "NUDT19"   "YWHAH"    "ESR2"     "ITGB1BP1" "ACSL4"   
#>  [91] "IDH2"     "PEX13"    "SLC35B2"  "ABCD1"    "PABPC1"   "SLC25A19"
#>  [97] "FDPS"     "TSPO"     "CRABP1"   "SOD2"     "MSH2"     "CTPS1"   
#> [103] "ECI2"     "UGT2B17" 
#> 
#> $le$HALLMARK_APICAL_JUNCTION
#>   [1] "GAMT"      "EVL"       "WNK4"      "IRS1"      "LIMA1"     "TAOK2"    
#>   [7] "RASA1"     "SHROOM2"   "PCDH1"     "CADM2"     "PTEN"      "AMIGO2"   
#>  [13] "CTNNA1"    "ITGA3"     "NEGR1"     "CRAT"      "NF1"       "JAM3"     
#>  [19] "BAIAP2"    "LAMA3"     "GTF2F1"    "WASL"      "TMEM8B"    "KRT31"    
#>  [25] "TSC1"      "CERCAM"    "CDH11"     "ITGA2"     "FBN1"      "CLDN5"    
#>  [31] "CDH6"      "SLIT2"     "AMIGO1"    "PKD1"      "CRB3"      "VWF"      
#>  [37] "CD34"      "ADAM23"    "CLDN19"    "SORBS3"    "MYH10"     "TRO"      
#>  [43] "GNAI2"     "JUP"       "ARHGEF6"   "NRXN2"     "INPPL1"    "TJP1"     
#>  [49] "LAYN"      "VCAN"      "TSPAN4"    "PPP2R2C"   "THY1"      "CDSN"     
#>  [55] "ALOX15B"   "NLGN3"     "EXOC4"     "VAV2"      "PARVA"     "COL16A1"  
#>  [61] "NLGN2"     "MMP2"      "STX4"      "COL17A1"   "DMP1"      "AKT2"     
#>  [67] "CTNND1"    "CNTN1"     "FLNC"      "CDH8"      "CLDN11"    "SHC1"     
#>  [73] "CLDN7"     "SYMPK"     "AMH"       "PIK3R3"    "ATP1A3"    "PECAM1"   
#>  [79] "LDLRAP1"   "INSIG1"    "CLDN18"    "CDH1"      "MADCAM1"   "HRAS"     
#>  [85] "TIAL1"     "HADH"      "ITGA10"    "TNFRSF11B" "CADM3"     "NRAP"     
#>  [91] "RRAS"      "PIK3CB"    "MYL9"      "TUBG1"     "ADRA1B"    "MAPK11"   
#>  [97] "ACTC1"     "PARD6G"    "B4GALT1"   "ACTA1"     "MAPK13"    "NEXN"     
#> [103] "ACTG1"     "PTK2"      "CD99"      "ADAM9"     "SKAP2"     "CLDN15"   
#> [109] "NFASC"     "ACTN3"     "VCL"       "DLG1"      "KCNH2"     "CDH15"    
#> [115] "DSC1"      "ADAMTS5"   "GRB7"      "ITGA9"     "CD276"     "THBS3"    
#> [121] "RAC2"      "IKBKG"     "PTPRC"     "SGCE"      "CD209"     "AKT3"     
#> [127] "BMP1"      "ACTN2"     "MAPK14"    "MYL12B"    "CD86"      "TGFBI"    
#> [133] "PLCG1"     "MMP9"      "ITGB1"     "CNN2"      "SDC3"      "PBX2"     
#> [139] "CAP1"      "MDK"       "CD274"     "MYH9"      "ACTG2"     "SRC"      
#> [145] "SYK"       "DHX16"     "MAP4K2"    "LAMB3"     "VASP"      "EPB41L2"  
#> [151] "ADAM15"    "VCAM1"     "CLDN8"     "CDH4"      "ICAM2"     "ACTN1"    
#> [157] "RHOF"      "LAMC2"     "CLDN14"    "ITGB4"     "TRAF1"     "YWHAH"    
#> [163] "ICAM5"     "CLDN4"     "NF2"       "ZYX"       "ARPC2"     "DSC3"     
#> [169] "CLDN9"     "ICAM4"     "SPEG"      "ACTN4"     "ICAM1"     "SLC30A3"  
#> [175] "GNAI1"     "MVD"       "ACTB"      "CALB2"     "PFN1"      "COL9A1"   
#> [181] "CDK8"      "SIRPA"     "NRTN"      "MPZL1"     "CX3CL1"    "CLDN6"    
#> [187] "MPZL2"     "MSN"       "EGFR"      "CDH3"      "FSCN1"     "RSU1"     
#> [193] "FYB1"      "MAP3K20"   "NECTIN1"   "NECTIN2"   "NECTIN3"   "NECTIN4"  
#> [199] "NHERF4"    "PALS1"    
#> 
#> $le$HALLMARK_APICAL_SURFACE
#>  [1] "GATA3"    "ATP8B1"   "RTN4RL1"  "GSTM3"    "SCUBE1"   "SHROOM2" 
#>  [7] "HSPB1"    "SLC34A3"  "MDGA1"    "ADIPOR2"  "TMEM8B"   "AFAP1L2" 
#> [13] "SULF2"    "BRCA1"    "LYPD3"    "CROCC"    "THY1"     "ADAM10"  
#> [19] "NTNG1"    "SLC2A4"   "NCOA6"    "SRPX"     "CD160"    "FLOT2"   
#> [25] "B4GALT1"  "PKHD1"    "GAS1"     "AKAP7"    "IL2RB"    "ATP6V0A4"
#> [31] "PCSK9"    "DCBLD2"   "IL2RG"    "GHRL"     "EPHB4"    "MAL"     
#> [37] "APP"      "PLAUR"    "EFNA5"    "CX3CL1"   "LYN"      "RHCG"    
#> [43] "CRYBG1"   "SLC22A12"
#> 
#> $le$HALLMARK_HEME_METABOLISM
#>   [1] "BCAM"     "BTG2"     "BTRC"     "TMEM9B"   "ALAD"     "PIGQ"    
#>   [7] "ADD1"     "VEZF1"    "EZH1"     "ALDH6A1"  "NUDT4"    "LRP10"   
#>  [13] "SELENBP1" "DCAF10"   "RBM5"     "TRIM58"   "BLVRA"    "EPOR"    
#>  [19] "SLC22A4"  "RNF123"   "HEBP1"    "SLC25A38" "HAGH"     "ABCG2"   
#>  [25] "ARHGEF12" "GDE1"     "SIDT2"    "DCAF11"   "GAPVD1"   "KHNYN"   
#>  [31] "MINPP1"   "GATA1"    "ATP6V0A1" "BLVRB"    "TRAK2"    "SLC30A1" 
#>  [37] "BNIP3L"   "NCOA4"    "FN3K"     "TNS1"     "CLCN3"    "DAAM1"   
#>  [43] "FECH"     "YPEL5"    "CAST"     "TSPAN5"   "HBB"      "ADIPOR1" 
#>  [49] "ACP5"     "KLF3"     "TAL1"     "CAT"      "NFE2L1"   "ELL2"    
#>  [55] "PPOX"     "ERMAP"    "SYNJ1"    "PRDX2"    "ISCA1"    "CCDC28A" 
#>  [61] "CDC27"    "FBXO9"    "RAP1GAP"  "NR3C1"    "FBXO34"   "NNT"     
#>  [67] "CA2"      "KAT2B"    "NFE2"     "NEK7"     "EPB42"    "BMP2K"   
#>  [73] "GYPE"     "TOP1"     "LMO2"     "GCLC"     "PSMD9"    "MKRN1"   
#>  [79] "SNCA"     "CTNS"     "LPIN2"    "ALAS2"    "UROS"     "EPB41"   
#>  [85] "HBD"      "UCP2"     "MYL4"     "TCEA1"    "CIR1"     "EIF2AK1" 
#>  [91] "SLC4A1"   "SLC11A2"  "ATG4A"    "ANK1"     "GLRX5"    "MXI1"    
#>  [97] "UROD"     "OSBP2"    "SPTA1"    "ARL2BP"   "PGLS"     "FOXJ2"   
#> [103] "MBOAT2"   "MARK3"    "CDR2"     "RHCE"     "USP15"    "BSG"     
#> [109] "RHD"      "FOXO3"    "P4HA2"    "AQP3"     "RANBP10"  "ALDH1L1" 
#> [115] "MPP1"     "CLIC2"    "ABCB6"    "PPP2R5B"  "PC"       "RNF19A"  
#> [121] "ENDOD1"   "BPGM"     "FBXO7"    "IGSF3"    "KEL"      "CTSB"    
#> [127] "C3"       "MGST3"    "TNRC6B"   "SPTB"     "RAD23A"   "GCLM"    
#> [133] "GYPC"     "XPO7"     "SLC7A11"  "CTSE"     "UBAC1"    "MOCOS"   
#> [139] "BACH1"    "CCND3"    "SEC14L1"  "SDCBP"    "RCL1"     "TFRC"    
#> [145] "DCUN1D1"  "ACSL6"    "HTRA2"    "MOSPD1"   "PICALM"   "KLF1"    
#> [151] "RIOK3"    "SLC6A9"   "SLC6A8"   "LAMP2"    "HTATIP2"  "TFDP2"   
#> [157] "CPOX"     "SLC2A1"   "HBQ1"     "ICAM4"    "PDZK1IP1" "MFHAS1"  
#> [163] "XK"       "HMBS"     "OPTN"     "AGPAT4"   "MAP2K3"   "ADD2"    
#> [169] "FTCD"     "TRIM10"   "E2F2"     "SLC30A10" "SLC10A3"  "NARF"    
#> [175] "SMOX"     "HDGF"     "ASNS"     "SLC25A37" "TMCC2"    "RBM38"   
#> [181] "GMPS"     "ACKR1"    "AHSP"     "CA1"      "CROCCP2"  "DMTN"    
#> [187] "GYPA"     "GYPB"     "H1-0"     "H4C3"     "HBBP1"    "HBZ"     
#> [193] "KDM7A"    "MARCHF2"  "MARCHF8"  "RHAG"     "SLC66A2"  "TENT5C"  
#> [199] "TSPO2"    "TYR"     
#> 
#> $le$HALLMARK_ANGIOGENESIS
#>  [1] "SERPINA5" "THBD"     "POSTN"    "COL3A1"   "MSX1"     "LUM"     
#>  [7] "FGFR1"    "COL5A2"   "LRPAP1"   "STC1"     "VTN"      "VCAN"    
#> [13] "VAV2"     "ITGAV"    "JAG1"     "SLCO2A1"  "KCNJ8"    "OLR1"    
#> [19] "TIMP1"    "FSTL1"    "NRP1"     "CCND2"    "PDGFA"    "LPL"     
#> [25] "PTK2"     "JAG2"     "SPP1"     "S100A4"   "VEGFA"    "APP"     
#> [31] "CXCL6"    "PRG2"     "TNFRSF21" "APOH"     "PF4"      "PGLYRP1" 
#> 
#> $le$HALLMARK_XENOBIOTIC_METABOLISM
#>   [1] "ESR1"      "FBP1"      "MCCC2"     "TMBIM6"    "IGFBP4"    "CSAD"     
#>   [7] "ELOVL5"    "MPP2"      "ACOX2"     "UGDH"      "NINJ1"     "CROT"     
#>  [13] "CYB5A"     "FAH"       "SERPINA6"  "ENTPD5"    "ETFDH"     "SAR1B"    
#>  [19] "ARPP19"    "CFB"       "SLC35D1"   "ACSM1"     "BLVRB"     "CYFIP2"   
#>  [25] "RETSAT"    "GNMT"      "DHPS"      "MAOA"      "CYP26A1"   "ACOX3"    
#>  [31] "PTGR1"     "DCXR"      "NQO1"      "PDK4"      "GSTM4"     "DHRS7"    
#>  [37] "MAN1A1"    "HACL1"     "CYP2E1"    "NFS1"      "PINK1"     "TAT"      
#>  [43] "IGF1"      "ASL"       "SLC46A3"   "AOX1"      "GSS"       "LEAP2"    
#>  [49] "PTGES"     "CYP4F2"    "CAT"       "PTGES3"    "NMT1"      "F10"      
#>  [55] "EPHX1"     "CD36"      "VTN"       "DHRS1"     "GAD1"      "CES1"     
#>  [61] "JUP"       "GSR"       "CDO1"      "GSTT2"     "SLC6A6"    "SERPINE1" 
#>  [67] "RAP1GAP"   "ADH1C"     "CA2"       "ABHD6"     "ENPEP"     "SERTAD1"  
#>  [73] "PDLIM5"    "GSTA3"     "CYP17A1"   "ALDH9A1"   "DDAH2"     "IGFBP1"   
#>  [79] "CYP1A1"    "ITIH4"     "ATOH8"     "GCLC"      "ALDH2"     "ACP2"     
#>  [85] "CASP6"     "CNDP2"     "ABCC3"     "LPIN2"     "FBLN1"     "ECH1"     
#>  [91] "BPHL"      "RBP4"      "HNF4A"     "CYP2S1"    "SPINT2"    "SLC12A4"  
#>  [97] "MBL2"      "IL1R1"     "SLC6A12"   "PEMT"      "TMEM97"    "CDA"      
#> [103] "ARG1"      "CYP2J2"    "HES6"      "PROS1"     "HMOX1"     "ACOX1"    
#> [109] "BCAR1"     "TPST1"     "HGFAC"     "ATP2A2"    "APOE"      "PMM1"     
#> [115] "GSTO1"     "GCKR"      "ITIH1"     "SLC35B1"   "PAPSS2"    "ACO2"     
#> [121] "PTGDS"     "AKR1C3"    "UPB1"      "SSR3"      "TTPA"      "TNFRSF1A" 
#> [127] "MTHFD1"    "PC"        "IRF8"      "ANGPTL3"   "FMO3"      "TMEM176B" 
#> [133] "LONP1"     "PSMB10"    "ALDH3A1"   "IDH1"      "FMO1"      "HSD11B1"  
#> [139] "AKR1C2"    "SLC1A5"    "CRP"       "ABCD2"     "GABARAPL1" "PYCR1"    
#> [145] "TGFB2"     "MT2A"      "GCH1"      "LCAT"      "AHCY"      "ID2"      
#> [151] "ADH5"      "ACP1"      "ALAS1"     "COMT"      "CYP27A1"   "FAS"      
#> [157] "SLC22A1"   "BCAT1"     "HPRT1"     "ABCC2"     "CYP2C18"   "ARG2"     
#> [163] "AP4B1"     "NDRG2"     "DDT"       "CCL25"     "ETS2"      "CBR1"     
#> [169] "KYNU"      "TDO2"      "PGD"       "XDH"       "PTS"       "VNN1"     
#> [175] "AQP9"      "NPC1"      "PGRMC1"    "POR"       "SHMT2"     "GART"     
#> [181] "EPHA2"     "HSD17B2"   "DDC"       "PPARD"     "GCNT2"     "SMOX"     
#> [187] "FETUB"     "UPP1"      "ADH7"      "CYP1A2"    "F11"       "FABP1"    
#> [193] "G6PC1"     "HRG"       "KARS1"     "MARCHF6"   "PLG"       "REG1A"    
#> [199] "TKFC"      "TYR"      
#> 
#> $le$HALLMARK_FATTY_ACID_METABOLISM
#>   [1] "SLC22A5"   "HSD17B4"   "REEP6"     "ELOVL5"    "GLUL"      "ALAD"     
#>   [7] "HMGCL"     "UGDH"      "HSDL2"     "BMPR1B"    "BLVRA"     "SERINC1"  
#>  [13] "SUCLG2"    "HMGCS2"    "PDHB"      "ACOT2"     "ETFDH"     "LTC4S"    
#>  [19] "DHCR24"    "ALDH3A2"   "CRAT"      "PSME1"     "CPT1A"     "RETSAT"   
#>  [25] "CYP4A11"   "MAOA"      "ADIPOR2"   "PTPRG"     "FASN"      "AUH"      
#>  [31] "D2HGDH"    "GSTZ1"     "ALDOA"     "HSD17B7"   "RDH16"     "ACADS"    
#>  [37] "SUCLA2"    "AOC3"      "ENO2"      "INMT"      "RAP1GDS1"  "AQP7"     
#>  [43] "ACSS1"     "GPD2"      "CRYZ"      "MGLL"      "XIST"      "CYP4A22"  
#>  [49] "CPT2"      "EPHX1"     "CD36"      "GPD1"      "HADHB"     "CA4"      
#>  [55] "ACADVL"    "CIDEA"     "MCEE"      "NTHL1"     "ENO3"      "ADH1C"    
#>  [61] "ACAA1"     "HPGD"      "BCKDHB"    "CA2"       "HSD17B11"  "PCBD1"    
#>  [67] "ECHS1"     "ALDH9A1"   "RDH11"     "ALDH1A1"   "GCDH"      "CYP1A1"   
#>  [73] "UBE2L6"    "NBN"       "HAO2"      "MLYCD"     "UROS"      "ECH1"     
#>  [79] "GRHPR"     "BPHL"      "DLST"      "ACOT8"     "TP53INP2"  "ACADM"    
#>  [85] "SDHC"      "HADH"      "HIBCH"     "UROD"      "ERP29"     "EHHADH"   
#>  [91] "APEX1"     "ACOX1"     "IDH3B"     "G0S2"      "LGALS1"    "ACO2"     
#>  [97] "HSP90AA1"  "CD1D"      "HSPH1"     "CEL"       "MIF"       "ALDH3A1"  
#> [103] "LDHA"      "DECR1"     "SDHA"      "IDH1"      "HMGCS1"    "FMO1"     
#> [109] "ACADL"     "ACSM3"     "GABARAPL1" "DLD"       "SMS"       "ACSL1"    
#> [115] "ACSL5"     "PRDX6"     "GAPDHS"    "ACAA2"     "HSD17B10"  "IDI1"     
#> [121] "ACAT2"     "OSTC"      "SDHD"      "SUCLG1"    "METAP1"    "FH"       
#> [127] "IDH3G"     "YWHAH"     "CPOX"      "NCAPH2"    "CBR1"      "NSDHL"    
#> [133] "ME1"       "ACSL4"     "TDO2"      "HCCS"      "PTS"       "MDH2"     
#> [139] "IL4I1"     "VNN1"      "CBR3"      "S100A10"   "AADAT"     "ADSL"     
#> [145] "MDH1"      "CA6"       "PPARA"     "ODC1"      "PDHA1"     "ADH7"     
#> [151] "ECI1"      "ECI2"      "FABP1"     "FABP2"     "GAD2"      "H2AZ1"    
#> [157] "KMT5A"     "MIX23"    
#> 
#> $le$HALLMARK_P53_PATHWAY
#>   [1] "ABAT"     "ANKRA2"   "WWP1"     "GLS2"     "HEXIM1"   "BTG2"    
#>   [7] "SLC19A2"  "CGRRF1"   "PRKAB1"   "ACVR1B"   "XPC"      "NINJ1"   
#>  [13] "KIF13B"   "SERTAD3"  "FUCA1"    "SP1"      "TOB1"     "RB1"     
#>  [19] "MDM2"     "CTSF"     "RAB40C"   "RXRA"     "RPS27L"   "SESN1"   
#>  [25] "FDXR"     "SLC35D1"  "CSRNP2"   "ZBTB16"   "ABCC5"    "CYFIP2"  
#>  [31] "DDB2"     "RETSAT"   "MKNK2"    "BAIAP2"   "ABHD4"    "CCNG1"   
#>  [37] "ZNF365"   "MXD4"     "DCXR"     "IER3"     "CDKN2AIP" "ISCU"    
#>  [43] "F2R"      "TRAFD1"   "STOM"     "RRP8"     "TXNIP"    "NHLH2"   
#>  [49] "FOS"      "PLK2"     "BLCAP"    "TGFB1"    "VWA5A"    "IP6K2"   
#>  [55] "TSPYL2"   "PVT1"     "PHLDA3"   "EPHX1"    "GPX2"     "CDKN1A"  
#>  [61] "PTPRE"    "TRIB3"    "PPM1D"    "TSC22D1"  "MAPKAPK3" "KLF4"    
#>  [67] "PROCR"    "RCHY1"    "ALOX15B"  "AK1"      "TP63"     "CDH13"   
#>  [73] "TP53"     "RAD51C"   "CEBPA"    "NOL8"     "INHBB"    "HDAC3"   
#>  [79] "FAM162A"  "APAF1"    "CTSD"     "DGKA"     "JUN"      "HRAS"    
#>  [85] "ERCC5"    "HINT1"    "ZMAT3"    "RGS16"    "CCND2"    "OSGIN1"  
#>  [91] "RAD9A"    "NUPR1"    "RPL18"    "PDGFA"    "RRAD"     "PLXNB2"  
#>  [97] "PLK3"     "EPS8L2"   "VAMP8"    "ZFP36L1"  "MXD1"     "RALGDS"  
#> [103] "HMOX1"    "GM2A"     "TCHH"     "FOXO3"    "PMM1"     "GADD45A" 
#> [109] "BTG1"     "JAG2"     "DEF6"     "CD81"     "DDIT3"    "RPL36"   
#> [115] "CLCA2"    "VDR"      "CDKN2B"   "CDK5R1"   "TM7SF3"   "DRAM1"   
#> [121] "TRIAP1"   "HSPA4L"   "IER5"     "FGF13"    "TRAF4"    "CCNK"    
#> [127] "POLH"     "ATF3"     "PPP1R15A" "TPD52L1"  "SLC7A11"  "HBEGF"   
#> [133] "PITPNC1"  "SLC3A2"   "TNNI1"    "SFN"      "SDC1"     "CCND3"   
#> [139] "CASP1"    "BAX"      "AEN"      "FAS"      "TPRKB"    "POM121"  
#> [145] "EI24"     "LIF"      "TNFSF9"   "S100A4"   "IFI30"    "RPS12"   
#> [151] "ITGB4"    "FBXW7"    "RNF19B"   "PCNA"     "SAT1"     "TAX1BP3" 
#> [157] "TCN2"     "PTPN14"   "TAP1"     "KRT17"    "IL1A"     "BMP2"    
#> [163] "PRMT2"    "NUDT15"   "SEC61A1"  "DDIT4"    "SPHK1"    "ADA"     
#> [169] "APP"      "DNTTIP2"  "SOCS1"    "TM4SF1"   "SERPINB5" "S100A10" 
#> [175] "CD82"     "NDRG1"    "NOTCH1"   "BAK1"     "EPHA2"    "TGFA"    
#> [181] "PERP"     "RHBDF2"   "STEAP3"   "KLK8"     "CDKN2A"   "ST14"    
#> [187] "IRAK1"    "UPP1"     "RAP2B"    "LDHB"     "CCP110"   "COQ8A"   
#> [193] "ELP1"     "H1-2"     "H2AC25"   "H2AJ"     "IRAG2"    "PIDD1"   
#> [199] "RACK1"    "WRAP73"  
#> 
#> $le$HALLMARK_HEDGEHOG_SIGNALING
#>  [1] "CELSR1" "TLE3"   "RASA1"  "NF1"    "UNC5C"  "RTN1"   "GLI1"   "NRCAM" 
#>  [9] "SCG2"   "LDB1"   "SLIT1"  "THY1"   "HEY1"   "CRMP1"  "NKX6-1" "NRP1"  
#> [17] "HEY2"   "SHH"    "OPHN1"  "CNTFR"  "NRP2"   "CDK5R1" "AMOT"   "MYH9"  
#> [25] "ACHE"   "DPYSL2" "ETS2"   "VEGFA"  "L1CAM"  "PTCH1"  "PML"    "TLE1"  
#> [33] "VLDLR"  "CDK6"   "ADGRG1" "PLG"   
#> 
#> $le$HALLMARK_ANDROGEN_RESPONSE
#>   [1] "SPDEF"    "INPP4B"   "GPD1L"    "ELOVL5"   "CCND1"    "STK39"   
#>   [7] "ABHD2"    "BMPR1B"   "SORD"     "APPBP2"   "KRT8"     "DHCR24"  
#>  [13] "IQGAP2"   "ZMIZ1"    "KRT19"    "AKAP12"   "AZGP1"    "NCOA4"   
#>  [19] "LIFR"     "MAK"      "HMGCR"    "TPD52"    "SEC24D"   "NKX3-1"  
#>  [25] "HERC3"    "STEAP4"   "ELL2"     "SCD"      "GSR"      "TSC22D1" 
#>  [31] "SLC38A2"  "VAPA"     "KLK3"     "CAMKK2"   "ITGAV"    "SPCS3"   
#>  [37] "DNAJB9"   "HPGD"     "MAP7"     "ACSL3"    "HSD17B14" "PDLIM5"  
#>  [43] "TMEM50A"  "AKT1"     "KLK2"     "SRP19"    "INSIG1"   "NGLY1"   
#>  [49] "PTK2B"    "PTPN21"   "HOMER2"   "PIAS1"    "XRCC5"    "TARP"    
#>  [55] "LMAN1"    "MAF"      "PA2G4"    "ARID5B"   "B4GALT1"  "ZBTB10"  
#>  [61] "MERTK"    "PMEPA1"   "UBE2I"    "B2M"      "ELK4"     "FKBP5"   
#>  [67] "RRP12"    "PGM3"     "HMGCS1"   "RAB4A"    "FADS1"    "TMPRSS2" 
#>  [73] "SLC26A2"  "SGK1"     "ADAMTS1"  "ANKH"     "SMS"      "TNFAIP8" 
#>  [79] "CCND3"    "CDC14B"   "ADRM1"    "IDI1"     "ACTN1"    "UAP1"    
#>  [85] "GNAI3"    "SAT1"     "DBI"      "ABCC4"    "MYL12A"   "UBE2J1"  
#>  [91] "XRCC6"    "RPS6KA3"  "ALDH1A3"  "NDRG1"    "SRF"      "CDK6"    
#>  [97] "CENPN"    "GUCY1A1"  "H1-0"     "PLPP1"    "SELENOP" 
#> 
#> $le$HALLMARK_APOPTOSIS
#>   [1] "RHOB"      "RARA"      "ERBB3"     "KRT18"     "BTG2"      "CCND1"    
#>   [7] "ADD1"      "MADD"      "RNASEL"    "PLAT"      "BCL2L1"    "TIMP3"    
#>  [13] "HSPB1"     "EGR3"      "FDXR"      "PSEN1"     "RHOT2"     "CREBBP"   
#>  [19] "HGF"       "RETSAT"    "PDCD4"     "BNIP3L"    "XIAP"      "LEF1"     
#>  [25] "PDGFRB"    "BIK"       "IER3"      "PPT1"      "CASP7"     "BCL2L11"  
#>  [31] "F2R"       "PSEN2"     "PMAIP1"    "MGMT"      "BCL2L2"    "CDKN1B"   
#>  [37] "TXNIP"     "ENO2"      "CASP9"     "IGFBP6"    "DCN"       "SLC20A1"  
#>  [43] "GPX4"      "LUM"       "SQSTM1"    "TGFBR3"    "GUCY2D"    "SPTAN1"   
#>  [49] "BRCA1"     "CDKN1A"    "SMAD7"     "FEZ1"      "GSR"       "AIFM3"    
#>  [55] "TNFSF10"   "CLU"       "MMP2"      "CAV1"      "GSTM1"     "GADD45B"  
#>  [61] "ERBB2"     "AVPR1A"    "ROCK1"     "DPYD"      "ETF1"      "DIABLO"   
#>  [67] "CTH"       "CASP6"     "TIMP2"     "WEE1"      "TIMP1"     "NEFH"     
#>  [73] "JUN"       "BGN"       "GPX1"      "CASP8"     "DNAJC3"    "CCND2"    
#>  [79] "IFITM3"    "NEDD9"     "SOD1"      "GPX3"      "PAK1"      "CD69"     
#>  [85] "HMOX1"     "CD44"      "CYLD"      "CTNNB1"    "PTK2"      "GADD45A"  
#>  [91] "PEA15"     "LGALS3"    "BCL10"     "BMF"       "DDIT3"     "PPP2R5B"  
#>  [97] "EMP1"      "IRF1"      "DAP3"      "PLCB2"     "EREG"      "VDAC2"    
#> [103] "RELA"      "GSN"       "BCAP31"    "IL1B"      "IGF2R"     "ATF3"     
#> [109] "MCL1"      "TGFB2"     "IL18"      "GCH1"      "GNA15"     "CDK2"     
#> [115] "TOP2A"     "CD2"       "SATB1"     "CASP3"     "ANKH"      "FASLG"    
#> [121] "IL6"       "CASP4"     "CASP1"     "BAX"       "FAS"       "TNFRSF12A"
#> [127] "DNM1L"     "IFNB1"     "DNAJA1"    "HMGB2"     "DAP"       "CD14"     
#> [133] "LMNA"      "BIRC3"     "ISG20"     "DFFA"      "TNF"       "CCNA1"    
#> [139] "SAT1"      "PRF1"      "TAP1"      "IL1A"      "EBP"       "BMP2"     
#> [145] "APP"       "IFNGR1"    "CD38"      "TSPO"      "CASP2"     "ANXA1"    
#> [151] "CFLAR"     "PPP3R1"    "CDC25B"    "BID"       "BCL2L10"   "SOD2"     
#> [157] "BTG3"      "F2"        "H1-0"      "PLPPR4"    "SC5D"     
#> 
#> $le$HALLMARK_ADIPOGENESIS
#>   [1] "CMBL"     "REEP6"    "REEP5"    "CYP4B1"   "CCNG2"    "PTGER3"  
#>   [7] "ADCY6"    "BAZ2A"    "COQ5"     "FAH"      "TOB1"     "SULT1A1" 
#>  [13] "PRDX3"    "QDPR"     "OMD"      "SLC27A1"  "HSPB8"    "ELMOD3"  
#>  [19] "LTC4S"    "SPARCL1"  "CRAT"     "RETSAT"   "LIPE"     "PDCD4"   
#>  [25] "UQCRQ"    "ADIPOR2"  "EPHX2"    "GPHN"     "FZD4"     "NDUFA5"  
#>  [31] "ABCB8"    "RREB1"    "MAP4K3"   "ALDOA"    "CD302"    "STOM"    
#>  [37] "DHRS7"    "NKIRAS1"  "LIFR"     "ACADS"    "BCKDHA"   "SCP2"    
#>  [43] "ADIPOQ"   "GPX4"     "GPD2"     "MGLL"     "FABP4"    "ITSN1"   
#>  [49] "ARAF"     "ESYT1"    "PPP1R15B" "DHRS7B"   "CPT2"     "GHITM"   
#>  [55] "CAT"      "GPAM"     "PFKFB3"   "SNCG"     "NMT1"     "CD36"    
#>  [61] "UQCR11"   "CIDEA"    "MRAP"     "ACLY"     "PHLDB1"   "CD151"   
#>  [67] "COX6A1"   "ITIH5"    "CHUK"     "AGPAT3"   "ABCA1"    "DNAJB9"  
#>  [73] "SORBS1"   "LPCAT3"   "PPARG"    "LEP"      "ITGA7"    "ECHS1"   
#>  [79] "LAMA4"    "RAB34"    "VEGFB"    "G3BP2"    "ELOVL6"   "DBT"     
#>  [85] "DGAT1"    "CS"       "BCL6"     "COL15A1"  "UCK1"     "APLP2"   
#>  [91] "ALDH2"    "DRAM2"    "IDH3A"    "UCP2"     "RNF11"    "ECH1"    
#>  [97] "UBC"      "PPM1B"    "TST"      "RTN3"     "SSPN"     "COX8A"   
#> [103] "UBQLN1"   "ACADM"    "SDHC"     "HADH"     "BCL2L13"  "ENPP2"   
#> [109] "HIBCH"    "ETFB"     "PEMT"     "GRPEL1"   "SOD1"     "GBE1"    
#> [115] "GPX3"     "PEX14"    "DNAJC15"  "NDUFAB1"  "LPL"      "ACOX1"   
#> [121] "MYLK"     "APOE"     "UQCRC1"   "GADD45A"  "PIM3"     "SLC25A1" 
#> [127] "NDUFB7"   "STAT5A"   "ACO2"     "SLC25A10" "PHYH"     "AIFM1"   
#> [133] "COL4A1"   "CYC1"     "DECR1"    "C3"       "IDH1"     "SCARB1"  
#> [139] "SAMM50"   "ARL4A"    "ACADL"    "NDUFS3"   "SLC1A5"   "UQCR10"  
#> [145] "SDHB"     "MGST3"    "RETN"     "DLD"      "JAGN1"    "ORM1"    
#> [151] "IMMT"     "YWHAG"    "ACAA2"    "SLC19A1"  "PFKL"     "CHCHD10" 
#> [157] "ANGPT1"   "MTCH2"    "SUCLG1"   "ANGPTL4"  "RIOK3"    "COX7B"   
#> [163] "IDH3G"    "DDT"      "PREB"     "PTCD3"    "TALDO1"   "DLAT"    
#> [169] "CDKN2C"   "ME1"      "COQ9"     "ESRRA"    "MDH2"     "TKT"     
#> [175] "TANK"     "POR"      "DHCR7"    "IFNGR1"   "COQ3"     "AK2"     
#> [181] "MRPL15"   "PLIN2"    "CMPK1"    "PGM1"     "ATP1B3"   "SLC5A6"  
#> [187] "MCCC1"    "ATL2"     "ADIG"     "ATP5PO"   "CAVIN1"   "CAVIN2"  
#> [193] "GPAT4"    "MIGA2"    "MTARC2"   "NABP1"    "RMDN3"    "SLC66A3" 
#> [199] "SOWAHC"   "SQOR"    
#> 
#> $le$HALLMARK_PANCREAS_BETA_CELLS
#>  [1] "ABCC8"   "SPCS1"   "SRP14"   "HNF1A"   "G6PC2"   "STXBP1"  "SYT13"  
#>  [8] "INSM1"   "SRP9"    "GCK"     "PCSK2"   "SCGN"    "ELP4"    "FOXO1"  
#> [15] "LMO2"    "MAFB"    "CHGA"    "NKX6-1"  "PKLR"    "PCSK1"   "NKX2-2" 
#> [22] "DPP4"    "ISL1"    "SEC11A"  "VDR"     "AKT3"    "DCX"     "FOXA2"  
#> [29] "SRPRB"   "PAX6"    "PAK3"    "PDX1"    "GCG"     "IAPP"    "INS"    
#> [36] "NEUROD1" "NEUROG3" "PAX4"    "SLC2A2"  "SST"    
#> 
#> $le$HALLMARK_COAGULATION
#>   [1] "TMPRSS6"  "ACOX2"    "RAPGEF3"  "HPN"      "PLAT"     "PRSS23"  
#>   [7] "HMGCS2"   "TIMP3"    "CTSO"     "CFD"      "MST1"     "ANG"     
#>  [13] "CFB"      "F2RL2"    "HTRA1"    "THBD"     "ISCU"     "CRIP2"   
#>  [19] "SERPINC1" "COMP"     "SPARC"    "ITGA2"    "THBS1"    "FBN1"    
#>  [25] "MASP2"    "CASP9"    "APOA1"    "MMP11"    "MMP10"    "LTA4H"   
#>  [31] "LRP1"     "VWF"      "GNG12"    "PROZ"     "F10"      "CTSK"    
#>  [37] "PROC"     "PDGFB"    "F8"       "S100A13"  "SERPINE1" "FN1"     
#>  [43] "SERPINA1" "CLU"      "MMP2"     "ITGB3"    "CD9"      "TFPI2"   
#>  [49] "PECAM1"   "OLR1"     "CFH"      "PEF1"     "USP11"    "TIMP1"   
#>  [55] "F3"       "HNF4A"    "FGA"      "WDR1"     "DUSP6"    "APOC2"   
#>  [61] "KLF7"     "DUSP14"   "F12"      "MBL2"     "FGG"      "PROS1"   
#>  [67] "GNB2"     "ARF4"     "MSRB2"    "MMP14"    "APOC1"    "KLKB1"   
#>  [73] "DPP4"     "C8B"      "RAC1"     "A2M"      "ADAM9"    "CTSH"    
#>  [79] "ITIH1"    "CSRP1"    "FURIN"    "PLAU"     "RABIF"    "GP1BA"   
#>  [85] "C2"       "MMP3"     "CFI"      "SERPING1" "P2RY1"    "DCT"     
#>  [91] "CTSB"     "MEP1A"    "GSN"      "BMP1"     "C3"       "CPB2"    
#>  [97] "C1QA"     "MMP8"     "CAPN5"    "PLEK"     "MMP9"     "MMP15"   
#> [103] "LEFTY2"   "SIRT2"    "CTSE"     "GDA"      "C1S"      "CAPN2"   
#> [109] "RGN"      "C1R"      "TF"       "SERPINB2" "C8G"      "S100A1"  
#> [115] "FYN"      "LAMP2"    "MAFF"     "LGMN"     "MMP1"     "ANXA1"   
#> [121] "PREP"     "SH2B2"    "MMP7"     "KLK8"     "APOC3"    "C8A"     
#> [127] "C9"       "CPN1"     "CPQ"      "CTSV"     "F11"      "F13B"    
#> [133] "F2"       "F9"       "GP9"      "HRG"      "PF4"      "PLG"     
#> 
#> $le$HALLMARK_MYOGENESIS
#>   [1] "SPDEF"     "ADCY9"     "ERBB3"     "BHLHE40"   "KCNH1"     "OCEL1"    
#>   [7] "STC2"      "REEP1"     "CACNG1"    "PPP1R3C"   "CACNA1H"   "ITGB5"    
#>  [13] "HRC"       "VIPR1"     "RB1"       "MEF2D"     "CASQ1"     "HSPB8"    
#>  [19] "CFD"       "PVALB"     "CAMK2B"    "CRAT"      "TNNT1"     "ATP6AP1"  
#>  [25] "BAG1"      "LAMA2"     "SGCD"      "TSC2"      "NQO1"      "PSEN2"    
#>  [31] "COX7A1"    "CASQ2"     "MYOG"      "SGCG"      "MYH11"     "COL6A3"   
#>  [37] "FKBP1B"    "SPARC"     "AGL"       "PFKM"      "MB"        "IGF1"     
#>  [43] "MAPRE3"    "DMPK"      "BDKRB2"    "PGAM2"     "COL1A1"    "LDB3"     
#>  [49] "COL3A1"    "TGFB1"     "MYL3"      "MYO1C"     "AEBP1"     "SYNGR2"   
#>  [55] "TNNT3"     "SORBS3"    "CD36"      "SPTAN1"    "CDKN1A"    "SCD"      
#>  [61] "FHL1"      "FST"       "ADAM12"    "SGCA"      "MYL7"      "IGFBP7"   
#>  [67] "APOD"      "MYH1"      "MYH8"      "FXYD1"     "AK1"       "ENO3"     
#>  [73] "GABARAPL2" "CDH13"     "CLU"       "MYH4"      "CHRNA1"    "EIF4A2"   
#>  [79] "SORBS1"    "AKT2"      "GADD45B"   "ITGA7"     "FOXO4"     "SH2B1"    
#>  [85] "APLNR"     "TNNC1"     "MYOM2"     "PYGM"      "COL15A1"   "GAA"      
#>  [91] "FABP3"     "GNAO1"     "COX6A2"    "CTF1"      "COL6A2"    "MYH2"     
#>  [97] "MYL4"      "HDAC5"     "ANKRD2"    "MYOM1"     "SH3BGR"    "SSPN"     
#> [103] "AGRN"      "DTNA"      "CKB"       "CHRNB1"    "MEF2C"     "CKMT2"    
#> [109] "MYH3"      "GPX3"      "SLN"       "TAGLN"     "PDLIM7"    "FLII"     
#> [115] "PLXNB2"    "LSP1"      "PKIA"      "ACTC1"     "ACTA1"     "RIT1"     
#> [121] "MYLK"      "MEF2A"     "DAPK2"     "TNNC2"     "MYL6B"     "TCAP"     
#> [127] "FGF2"      "PTGIS"     "ACTN3"     "PDE4DIP"   "GJA5"      "NOS1"     
#> [133] "KCNH2"     "MYOZ1"     "TPM3"      "ATP2A1"    "MYBPH"     "MAPK12"   
#> [139] "PC"        "NAV2"      "DES"       "HSPB2"     "TPM2"      "GSN"      
#> [145] "ACTN2"     "PPFIA4"    "MYH7"      "ITGB1"     "TPD52L1"   "MYBPC3"   
#> [151] "PICK1"     "SIRT2"     "NCAM1"     "HBEGF"     "TNNI1"     "MYH9"     
#> [157] "SVIL"      "ACSL1"     "ABLIM1"    "COL4A2"    "CNN3"      "ACHE"     
#> [163] "PTP4A3"    "IGFBP3"    "SMTN"      "ITGB4"     "SLC6A8"    "CKM"      
#> [169] "SOD3"      "SPEG"      "EFS"       "BIN1"      "TNNI2"     "RYR1"     
#> [175] "SPHK1"     "APP"       "FDPS"      "KIFC3"     "TNNT2"     "DMD"      
#> [181] "PRNP"      "TEAD4"     "WWTR1"     "NOTCH1"    "KLF5"      "EPHB3"    
#> [187] "SCHIP1"    "CRYAB"     "MRAS"      "LPIN1"     "IFRD1"     "CAV3"     
#> [193] "CHRNG"     "CSRP3"     "DENND2B"   "LARGE1"    "MYF6"      "MYL1"     
#> [199] "MYL11"     "MYL2"     
#> 
#> $le$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
#>   [1] "RHOB"      "IGFBP4"    "MATN3"     "FUCA1"     "IGFBP2"    "ITGB5"    
#>   [7] "TIMP3"     "GJA1"      "FSTL3"     "ECM1"      "CAP2"      "HTRA1"    
#>  [13] "NTM"       "LAMA2"     "SGCD"      "LAMA3"     "LOXL1"     "PDGFRB"   
#>  [19] "ECM2"      "GPC1"      "SDC4"      "SLIT3"     "EDIL3"     "PMP22"    
#>  [25] "POSTN"     "MAGEE1"    "COMP"      "SGCG"      "SPOCK1"    "FMOD"     
#>  [31] "AREG"      "COL6A3"    "SPARC"     "CDH11"     "ENO2"      "ITGA2"    
#>  [37] "SFRP4"     "COL12A1"   "THBS1"     "COL1A2"    "FBN1"      "CDH6"     
#>  [43] "DCN"       "COL1A1"    "SLIT2"     "COL3A1"    "TGFB1"     "COL8A2"   
#>  [49] "LRRC15"    "CXCL12"    "MSX1"      "LRP1"      "LUM"       "COL5A2"   
#>  [55] "FAP"       "COL5A1"    "TGFBR3"    "BDNF"      "COL5A3"    "RGS4"     
#>  [61] "SCG2"      "EFEMP2"    "ABI3BP"    "PTHLH"     "MFAP5"     "GEM"      
#>  [67] "ADAM12"    "NID2"      "WNT5A"     "VCAN"      "QSOX1"     "THBS2"    
#>  [73] "SERPINE1"  "THY1"      "OXTR"      "FN1"       "ITGAV"     "COL16A1"  
#>  [79] "MXRA5"     "MMP2"      "COL7A1"    "PCOLCE"    "GADD45B"   "BASP1"    
#>  [85] "ITGB3"     "COL11A1"   "DAB2"      "PRRX1"     "TFPI2"     "ELN"      
#>  [91] "INHBA"     "VEGFC"     "COL6A2"    "DST"       "FBLN1"     "ACTA2"    
#>  [97] "MATN2"     "TIMP1"     "SGCB"      "MEST"      "PDLIM4"    "FSTL1"    
#> [103] "JUN"       "BGN"       "GREM1"     "ANPEP"     "TNFRSF11B" "FERMT2"   
#> [109] "MYL9"      "CTHRC1"    "MGP"       "TAGLN"     "ITGA5"     "SNAI2"    
#> [115] "FZD8"      "NOTCH2"    "DPYSL3"    "LOXL2"     "LOX"       "MMP14"    
#> [121] "CD44"      "WIPF1"     "MYLK"      "FBN2"      "EMP3"      "TPM1"     
#> [127] "GADD45A"   "LGALS1"    "PMEPA1"    "FGF2"      "NNMT"      "TGM2"     
#> [133] "CADM1"     "LAMC1"     "FBLN2"     "CDH2"      "COL4A1"    "APLP1"    
#> [139] "SPP1"      "MMP3"      "NT5E"      "GAS1"      "TNC"       "LAMA1"    
#> [145] "TPM2"      "BMP1"      "CAPG"      "TGFBI"     "COPA"      "ITGB1"    
#> [151] "ID2"       "FBLN5"     "CALD1"     "COL4A2"    "SDC1"      "IL15"     
#> [157] "TPM4"      "GLIPR1"    "IL6"       "CRLF1"     "FAS"       "TNFRSF12A"
#> [163] "VIM"       "VCAM1"     "PLOD2"     "PCOLCE2"   "FOXC2"     "CD59"     
#> [169] "LAMC2"     "IGFBP3"    "SNTB1"     "SLC6A8"    "FLNA"      "TNFAIP3"  
#> [175] "PFN2"      "SAT1"      "PPIB"      "VEGFA"     "MMP1"      "SERPINH1" 
#> [181] "IL32"      "PLAUR"     "SERPINE2"  "GPX7"      "DKK1"      "CXCL6"    
#> [187] "CALU"      "PLOD3"     "PTX3"      "SFRP1"     "MCM7"      "CXCL1"    
#> [193] "PVR"       "PLOD1"     "CCN1"      "CCN2"      "COLGALT1"  "CXCL8"    
#> [199] "P3H1"      "PRSS2"    
#> 
#> $le$HALLMARK_TGF_BETA_SIGNALING
#>  [1] "PPM1A"    "SMAD3"    "APC"      "RAB31"    "XIAP"     "BCAR3"   
#>  [7] "THBS1"    "LTBP2"    "TGFB1"    "SLC20A1"  "RHOA"     "ARID4B"  
#> [13] "NCOR2"    "BMPR2"    "CDK9"     "SMAD7"    "TJP1"     "SMAD1"   
#> [19] "SERPINE1" "KLF10"    "TGIF1"    "SKIL"     "UBE2D3"   "ACVR1"   
#> [25] "JUNB"     "HIPK2"    "CDH1"     "NOG"      "TRIM33"   "ENG"     
#> [31] "SMAD6"    "ID3"      "HDAC1"    "PPP1CA"   "CTNNB1"   "CDKN1C"  
#> [37] "PMEPA1"   "TGFBR1"   "FURIN"    "BMPR1A"   "ID1"      "SPTBN1"  
#> [43] "PPP1R15A" "LEFTY2"   "ID2"      "FNTA"     "IFNGR2"   "MAP3K7"  
#> [49] "SKI"      "SMURF1"   "SMURF2"   "FKBP1A"   "BMP2"     "WWTR1"
```

Once extracted, we can append this information directly into our results
dataframe using
[`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md).

``` r

if (!is.null(gene_lists)) {
  pa_annotated <- addgenesPA(
    pa_data    = pa_single,
    gene_lists = gene_lists
  )

  # Check the new columns added
  print(head(pa_annotated[, c("NAME", "NES", "FDR", "le_genes")], 3))
}
#>                                  NAME       NES FDR
#> 3899             HALLMARK_E2F_TARGETS -3.229399   0
#> 3900 HALLMARK_ESTROGEN_RESPONSE_EARLY  3.026576   0
#> 3901          HALLMARK_G2M_CHECKPOINT -2.979752   0
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             le_genes
#> 3899                                                 RAD50,WDR90,PAN2,UBR7,GSPT1,MLH1,PDS5B,CDKN1B,RFC1,LUC7L3,BRCA1,CDKN1A,PPM1D,RBBP7,PMS2,CTCF,TP53,IPO7,RAD51C,DCTPP1,CBX5,RPA3,NBN,WEE1,DUT,SHMT1,NAP1L1,TUBG1,TIMELESS,STAG1,DCK,BARD1,PA2G4,KIF22,LIG1,RAD1,PRPS1,PNN,NME1,NUP107,SMC3,POLE4,EIF2S1,POLE,SMC1A,PSMC3IP,RAD21,NAA38,HUS1,RPA2,CSE1L,POLD2,ING3,MYC,TOP2A,SNRPB,PRKDC,NOLC1,NUDT21,RACGAP1,PAICS,MXD3,BRCA2,RPA1,ESPL1,SLBP,ATAD2,HMGB2,UNG,NUP153,HELLS,TFRC,SPAG5,NOP56,BRMS1L,PSIP1,CIT,POLA2,POLD3,RNASEH2A,TRA2B,CENPM,PPP1R8,CDK4,ASF1B,XPO1,TIPIN,TK1,E2F8,SPC25,PCNA,SPC24,EXOSC8,TMPO,CDKN2C,HMMR,HNRNPD,LMNB1,DNMT1,PLK4,ILF3,SMC4,CDK1,RAN,POP7,HMGB3,RFC3,UBE2T,XRCC6,GINS1,RRM2,CDKN3,TACC3,BUB1B,POLD1,ZW10,AURKA,MCM4,KPNA2,MCM2,CENPE,GINS4,MAD2L1,TBRG4,CKS2,NUP205,SUV39H1,GINS3,KIF4A,DCLRE1B,EED,UBE2S,TCF19,AK2,ASF1A,KIF18B,MKI67,DSCC1,PTTG1,MTHFD2,BIRC5,CDC25B,SMC6,MCM3,SSRP1,DIAPH3,MYBL2,DLGAP5,USP1,TRIP13,ANP32E,MCM6,RAD51AP1,MSH2,EZH2,CHEK2,PLK1,SYNCRIP,CDCA3,DEPDC1,RFC2,CCNB2,CKS1B,MELK,CDKN2A,DONSON,PHF5A,PRIM2,KIF2C,AURKB,LYAR,PRDX4,CDC25A,MCM7,MCM5,NCAPD2,CHEK1,STMN1,DEK,RANBP1,NASP,CDCA8,LBR,TUBB,CDC20,HMGA1,CCNE1,CCP110,CNOT9,CTPS1,DDX39A,H2AX,H2AZ1,JPT1,MMS22L,MRE11,ORC2,ORC6,SRSF1,SRSF2
#> 3900 CA12,THSD4,MLPH,XBP1,ANXA9,TFF1,GREB1,TFF3,MYB,SLC22A5,ABAT,CELSR1,GFRA1,SLC39A6,MAPT,KCNK15,KDM4B,SYBU,ADCY9,WFS1,PGR,AR,SLC7A2,BCL2,UGCG,RARA,IL6ST,TTC39A,IGF1R,BHLHE40,SCNN1A,MAST4,MED13L,IGFBP4,RAB17,LRIG1,KRT18,MREG,SEMA3B,TPBG,RET,SLC19A2,ARL3,CISH,ELOVL5,CCND1,STC2,REEP1,TJP3,SIAH2,PDZK1,ABCA3,ELOVL2,SLC27A2,ADCY1,NRIP1,AFF1,PEX11A,AMFR,MYOF,ABHD2,DHRS2,ITPK1,TOB1,PRSS23,SULT2B1,SLC1A4,ASB13,CELSR2,CBFA2T3,FLNB,HSPB8,SEC14L2,KRT8,MUC1,ELF1,SNX24,GJA1,OLFM1,CANT1,DYNLT3,TIPARP,EGR3,RAB31,NPY1R,BLVRB,FKBP4,BAG1,KRT19,OVOL2,FASN,P2RY2,PMAIP1,SLC1A1,MSMB,CALCR,AREG,OLFML3,FOS,RAPGEFL1,NADSYN1,UNC119,NBL1,CXCL12,PTGES,MPPED2,SYNGR1,RHOBTB3,NCOR2,DLC1,CHPT1,RBBP8,GLA,KLF4,IL17RB,SOX3,RHOD,KLF10,RPS6KA2,TSKU,HES1,MICB,ESRP2,INHBB,CLDN7,FHL2,SLC24A3,SYT12,GAB2,TMEM164,WWC1,JAK2,KRT13,FDFT1,MED24,RASGRP1,SH3BP5,TMPRSS3,INPP5F,SLC37A1,B4GALT1,FRK,ALDH3B1,CD44,ZNF185,AQP3,PDLIM3,AKAP1,PAPSS2,TIAM1,TGM2,FARP1,NAV2,ENDOD1,ELF3,FKBP5,ADD3,RRP12,SCARB1,DHRS3,MYBL1,ISG20L2,CLIC3,TPD52L1,SLC26A2,MYC,TGIF2,SVIL,TUBB2B,SFN,ABLIM1,KRT15,MYBBP1A,PPIF,CYP26B1,PODXL,TFAP2C,SLC2A1,NXT1,BCL11B,KLK10,DHCR7,CALB2,SLC16A1,OPN3,HR,LAD1,SLC7A5,KCNK5,FOXC1,CCN5,DEPTOR,EEIG1,FCMR,KAZN,MINDY1,NHERF1,PLAAT3,RETREG1,TBC1D30
#> 3901                                                     TLE3,CCND1,NUMA1,ARID4A,SMAD3,CUL3,PURA,RPS6KA5,GSPT1,MNAT1,SLC38A1,BUB3,YTHDC1,SLC12A2,CCNT1,PDS5B,CDKN1B,ATRX,TGFB1,EGF,PRMT5,LIG3,CDC27,ABL1,CTCF,PAFAH1B1,ODF2,TOP1,NUP98,BCL3,RAD23B,G3BP1,SS18,CUL5,HMGN2,STAG1,HOXC10,HNRNPU,BARD1,UPF1,KIF22,NOTCH2,HSPA8,KPNB1,CUL4A,KIF5B,WRN,MEIS1,GINS2,CBX1,FOXN3,POLE,SMC1A,MEIS2,SQLE,SMARCC1,RAD21,MAPK14,HUS1,RBM14,RPA2,TNPO2,CCNF,MT2A,MYC,SMC2,DTYMK,TOP2A,EWSR1,KIF20B,NOLC1,DR1,RACGAP1,TRAIP,BRCA2,ESPL1,RBL1,PBK,CASP8AP2,NEK2,ATF5,POLA2,RASAL2,MARCKS,CDC6,TRA2B,CHMP1A,CDK4,XPO1,NUSAP1,CHAF1A,NUP50,TMPO,CDKN2C,HMMR,HNRNPD,LMNB1,PTTG3P,HIF1A,SFPQ,PLK4,ILF3,SMC4,SAP30,CDK1,POLQ,CUL1,FANCC,SLC7A1,E2F1,MTF2,HMGB3,PRC1,HIRA,TROAP,KIF11,CDKN3,TACC3,NCL,FBXO5,KIF15,AURKA,KPNA2,MCM2,CENPE,MAD2L1,CKS2,SUV39H1,KIF4A,KIF23,EFNA5,UBE2S,CENPF,MKI67,PTTG1,DMD,UBE2C,PML,CDC7,BIRC5,SNRPD1,CDC25B,INCENP,EXO1,CDC45,MCM3,MYBL2,E2F2,MCM6,TFDP1,TPX2,E2F4,EZH2,PLK1,BUB1,CCNA2,SYNCRIP,RAD54L,CCNB2,CKS1B,KATNA1,NDC80,DBF4,STIL,PRIM2,KIF2C,AURKB,CDC25A,MCM5,CHEK1,SLC7A5,STMN1,TTK,ODC1,NASP,LBR,AMD1,UCK2,CENPA,DKC1,CDC20,HMGA1,E2F3,DDX39A,H2AX,H2AZ1,H2AZ2,H2BC12,JPT1,KMT5A,KNL1,MAP3K20,NSD2,ORC5,ORC6,PRP4K,SRSF1,SRSF10,SRSF2,TENT4A
```

## 4. Visualizing Gene Expression inside Pathways (Heatmaps)

OmicsKit’s
[`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md)
function automates the process of drawing expression heatmaps for genes
associated with enriched pathways.

For this vignette, we display the leading-edge genes from
`"HALLMARK_ESTROGEN_RESPONSE_LATE"` in the `"ERpositive_vs_ERnegative"`
comparison. Samples are grouped by estrogen receptor status (`ER_group`)
and ordered as ER Positive and ER Negative.

**Note:** The heatmap is displayed from a pre-rendered static image to
keep vignette compilation lightweight and reproducible.

``` r

# Example code to generate the heatmap
heatmap_PA(
  expression_data = brca_rna_vst_or_logexpr_small,
  metadata        = brca_rna_metadata_tumor_normal,
  pa_data_annot   = pa_annotated,
  ranked_genes    = brca_ranked_genes$ERpositive_vs_ERnegative,
  plot_genes      = "le_genes",
  sample_col      = "sampleID",
  group_col       = "ER_group",
  out_dir         = "figures/heatmaps"
)
```

![Leading-edge expression heatmap for HALLMARK_ESTROGEN_RESPONSE_LATE
grouped by ER
status.](figures/heatmaps/HALLMARK_ESTROGEN_RESPONSE_LATE_heatmap.jpg)

Leading-edge expression heatmap for HALLMARK_ESTROGEN_RESPONSE_LATE
grouped by ER status.

## 5. Reducing Redundancy: Similarity and Clustering

Pathway analysis often returns redundant terms. OmicsKit reduces
redundancy by computing pairwise pathway similarity and clustering
related gene sets.

The code below illustrates how the silhouette plot can be inspected. The
displayed figure is pre-rendered from the ER positive vs ER negative
pathway similarity subset.

``` r

# Example code used during figure generation
brca_pa_clustering_er <- do_clust(brca_pa_similarity_er)
brca_pa_clustering_er$silhouette_plot
```

![Silhouette plot used to evaluate pathway clustering for the ER
positive vs ER negative subset.](figures/pathway_silhouette_plot.png)

Silhouette plot used to evaluate pathway clustering for the ER positive
vs ER negative subset.

## 6. Network Communities and Super-terms

To assign biological meaning to pathway clusters, OmicsKit can detect
communities within a pathway similarity network and extract
representative biological super-terms.

The original similarity object contains multiple comparisons. For this
vignette, the network figures were generated from the
`ERpositive_vs_ERnegative` subset only. This keeps the examples focused
on one biological contrast and avoids mixing nodes from multiple
comparisons.

``` r

# Example code used during figure generation
comparison <- "ERpositive_vs_ERnegative"

brca_pa_similarity_er <- brca_pa_similarity
keep <- grepl(paste0("^", comparison, "::"), rownames(brca_pa_similarity$jaccard_sim))

brca_pa_similarity_er$jaccard_sim <- brca_pa_similarity$jaccard_sim[
  keep,
  keep,
  drop = FALSE
]

# Remove the comparison prefix from node names for cleaner labels
node_names <- sub(
  paste0("^", comparison, "::"),
  "",
  rownames(brca_pa_similarity_er$jaccard_sim)
)

rownames(brca_pa_similarity_er$jaccard_sim) <- make.unique(node_names)
colnames(brca_pa_similarity_er$jaccard_sim) <- make.unique(node_names)
brca_pa_similarity_er$dist_mat <- as.dist(1 - brca_pa_similarity_er$jaccard_sim)

brca_pa_clustering_er <- do_clust(brca_pa_similarity_er)

net_results_er <- get_network_communities(
  x             = brca_pa_similarity_er,
  threshold     = 0.3,
  method        = "louvain",
  superterms    = TRUE,
  n_terms       = 3,
  remove_prefix = TRUE,
  seed          = 174
)

head(net_results_er$superterms$summary)

network_clust_gg(
  x                 = brca_pa_similarity_er,
  clust_result      = brca_pa_clustering_er,
  jaccard_threshold = 0.3,
  min_degree        = 1,
  superterms        = TRUE,
  superterm_data    = net_results_er$superterms,
  type              = "superterms",
  seed              = 174
)
```

### Network without super-term labels

This version emphasizes the topology of the ER positive vs ER negative
pathway similarity network without adding text labels.

![ER positive vs ER negative pathway similarity network without
super-term labels.](figures/pathway_network_er_clean_no_superterms.png)

ER positive vs ER negative pathway similarity network without super-term
labels.

### Network with community super-terms

This version labels pathway communities using automatically derived
super-terms.

![ER positive vs ER negative pathway similarity network with
community-level super-term
labels.](figures/pathway_network_er_superterms.png)

ER positive vs ER negative pathway similarity network with
community-level super-term labels.

### Individual community view

The individual view is useful when detailed community-level inspection
is preferred over a single combined network representation.

![Individual ER positive vs ER negative pathway community
view.](figures/pathway_network_er_individual_with_superterms.png)

Individual ER positive vs ER negative pathway community view.

### Combined network representation

The combined view provides a compact summary of the ER positive vs ER
negative pathway communities.

    #> Figure not found: figures/pathway_network_communities.png

## Session Info

``` r

sessionInfo()
#> R version 4.4.2 (2024-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26200)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> time zone: America/Bogota
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] OmicsKit_1.0.0.0000
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10         generics_0.1.4      tidyr_1.3.2        
#>  [4] lattice_0.22-6      digest_0.6.39       magrittr_2.0.5     
#>  [7] pROC_1.19.0.1       evaluate_1.0.5      grid_4.4.2         
#> [10] RColorBrewer_1.1-3  fastmap_1.2.0       Matrix_1.7-1       
#> [13] jsonlite_2.0.0      backports_1.5.1     survival_3.8-6     
#> [16] purrr_1.2.2         scales_1.4.0        textshaping_1.0.5  
#> [19] jquerylib_0.1.4     cli_3.6.6           rlang_1.2.0        
#> [22] splines_4.4.2       cachem_1.1.0        yaml_2.3.12        
#> [25] otel_0.2.0          tools_4.4.2         dplyr_1.2.1        
#> [28] ggplot2_4.0.3       BiocGenerics_0.52.0 broom_1.0.13       
#> [31] vctrs_0.7.3         R6_2.6.1            stats4_4.4.2       
#> [34] lifecycle_1.0.5     S4Vectors_0.44.0    fs_2.1.0           
#> [37] htmlwidgets_1.6.4   ragg_1.5.2          pkgconfig_2.0.3    
#> [40] desc_1.4.3          pkgdown_2.2.0       pillar_1.11.1      
#> [43] bslib_0.10.0        gtable_0.3.6        glue_1.8.1         
#> [46] Rcpp_1.1.1-1.1      systemfonts_1.3.2   xfun_0.54          
#> [49] tibble_3.3.1        tidyselect_1.2.1    rstudioapi_0.18.0  
#> [52] knitr_1.51          dichromat_2.0-0.1   farver_2.1.2       
#> [55] htmltools_0.5.9     patchwork_1.3.2     rmarkdown_2.31     
#> [58] compiler_4.4.2      S7_0.2.2
```
