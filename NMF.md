## Data Read and Quality Control

    muv41 = read.table("GSM3905417_MUV41.txt")
    muv41 = data_process(muv41)
    muv41[muv41<0]=0.000001
    dim(muv41)

    ## [1] 17140   332

## rank selection

This function is very time consuming and always cause problem when
rendering the md file. However I run it individually and get the
ultimate rank 3. Therefore I will skip this chunk when rendering the md
file

    my_rank = rank_selection(muv41, a = 2, b = 5)
    my_rank

    estimate = nmf(muv41, rank=2:5, method="brunet", nrun=10, seed=888)
    plot(estimate)

![](NMF_files/figure-markdown_strict/unnamed-chunk-3-1.png)

## gene metaprograms

    metaprogram = nmf_process(matrix = muv41, rank = 3, method = "brunet")
    metaprogram

    ## [[1]]
    ##  [1] "SERPINF1"   "CXCR4"      "MGP"        "PGF"        "MFAP4"     
    ##  [6] "ALDH1A3"    "HNRNPA1P10" "S100A11"    "HMGN2"      "BOC"       
    ## [11] "B2M"        "DNAJB1"     "FOS"        "FAP"        "NTRK3"     
    ## [16] "UQCRH"      "VIM"        "CCNG1"      "DDIT4"      "C6orf48"   
    ## [21] "TMEM123"    "SNHG1"      "SNRPB"      "TMEM98"     "FBL"       
    ## [26] "CCNB1IP1"   "ZFP36L1"    "SLCO2A1"    "SNRPD2"     "GNG5"      
    ## 
    ## [[2]]
    ##  [1] "STMN2"     "KCNIP4"    "CNTN2"     "ASIC1"     "ID2"       "NHLH1"    
    ##  [7] "PDZRN3"    "BSCL2"     "GNG3"      "CAMK4"     "CKB"       "TAGLN3"   
    ## [13] "CRIP2"     "NDUFB8"    "SCHIP1"    "PBX3"      "RIMS3"     "PTCHD2"   
    ## [19] "NARF"      "BTBD17"    "ELAVL2"    "DCTN1"     "SCG5"      "MIAT"     
    ## [25] "AFAP1"     "GPR56"     "TUBB2A"    "LINC00599" "PCBP4"     "NEFM"     
    ## 
    ## [[3]]
    ##  [1] "STMN2"     "CHGB"      "TSPAN7"    "NNAT"      "FEZ1"      "HMP19"    
    ##  [7] "SH3GL2"    "STMN4"     "KLHL13"    "ID2"       "TUBB2A"    "TMEM178A" 
    ## [13] "ARL6IP5"   "NEUROD6"   "MEG3"      "BLCAP"     "ENO2"      "GPM6B"    
    ## [19] "GABARAPL2" "GNG3"      "BEX2"      "GALNT18"   "IMPDH1"    "ATP6V1G2" 
    ## [25] "KIFAP3"    "GSTA4"     "RCAN2"     "APLP2"     "MAP1LC3B"  "ASNS"

## GO Enrichment analysis

    res1 = GOenrich(allgene = rownames(muv41), diffgene = metaprogram[[1]])

    ## 
    ## Building most specific GOs .....

    ## Loading required package: org.Hs.eg.db

    ## 

    ##  ( 11258 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 14896 GO terms and 33851 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 12931 genes annotated to the GO terms. )

    ## 
    ##           -- Classic Algorithm -- 
    ## 
    ##       the algorithm is scoring 1462 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    res2 = GOenrich(allgene = rownames(muv41), diffgene = metaprogram[[2]])

    ## 
    ## Building most specific GOs .....

    ##  ( 11258 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 14896 GO terms and 33851 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 12931 genes annotated to the GO terms. )

    ## 
    ##           -- Classic Algorithm -- 
    ## 
    ##       the algorithm is scoring 832 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    res3 = GOenrich(allgene = rownames(muv41), diffgene = metaprogram[[3]])

    ## 
    ## Building most specific GOs .....

    ##  ( 11258 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 14896 GO terms and 33851 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 12931 genes annotated to the GO terms. )

    ## 
    ##           -- Classic Algorithm -- 
    ## 
    ##       the algorithm is scoring 935 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    res <- list(res1,res2,res3)
    names(res) <- c(1,2,3)

    head(res1)

    ##        GO.ID                                     Term Annotated Significant
    ## 1 GO:0051384               response to glucocorticoid       103           5
    ## 2 GO:0031960               response to corticosteroid       113           5
    ## 3 GO:0007417       central nervous system development       829          10
    ## 4 GO:0007399               nervous system development      1944          14
    ## 5 GO:0010712 regulation of collagen metabolic process        24           3
    ## 6 GO:0009605            response to external stimulus      1783          13
    ##   Expected classicFisher         FDR        FC
    ## 1     0.22       2.5e-06 0.009858333 22.418516
    ## 2     0.24       3.9e-06 0.009858333 20.434576
    ## 3     1.80       5.0e-06 0.009858333  5.570825
    ## 4     4.21       1.4e-05 0.020702500  3.325874
    ## 5     0.05       1.8e-05 0.021294000 57.727679
    ## 6     3.86       3.1e-05 0.030560833  3.367178

    head(res2)

    ##         GO.ID                                                   Term Annotated
    ## 2  GO:0048710                regulation of astrocyte differentiation        19
    ## 75 GO:0010968                   regulation of microtubule nucleation        10
    ## 76 GO:0030397                                   membrane disassembly        10
    ## 77 GO:0048557                embryonic digestive tract morphogenesis        10
    ## 78 GO:0051081                           nuclear membrane disassembly        10
    ## 80 GO:0001973 G protein-coupled adenosine receptor signaling pathway        11
    ##    Significant Expected classicFisher FDR       FC
    ## 2            2     0.04       0.00065   1 52.35223
    ## 75           1     0.02       0.01993   1 49.73462
    ## 76           1     0.02       0.01993   1 49.73462
    ## 77           1     0.02       0.01993   1 49.73462
    ## 78           1     0.02       0.01993   1 49.73462
    ## 80           1     0.02       0.02190   1 45.21329

    head(res3)

    ##        GO.ID                                     Term Annotated Significant
    ## 2 GO:0006995 cellular response to nitrogen starvation        10           2
    ## 3 GO:0043562     cellular response to nitrogen levels        10           2
    ## 1 GO:0022411           cellular component disassembly       377           6
    ## 4 GO:0048699                    generation of neurons      1160           9
    ## 6 GO:0000045                   autophagosome assembly        97           3
    ## 7 GO:0051261                 protein depolymerization        97           3
    ##   Expected classicFisher       FDR        FC
    ## 2     0.02       0.00020 0.3943333 92.364286
    ## 3     0.02       0.00020 0.3943333 92.364286
    ## 1     0.82       0.00013 0.3943333  7.349943
    ## 4     2.51       0.00052 0.7689500  3.583097
    ## 6     0.21       0.00117 0.8609611 14.283137
    ## 7     0.21       0.00117 0.8609611 14.283137

## GO heatmap

    GOheatmap = plotGOEnrich (res,n = 3,fdr.cutoff = 2,fc.cutoff = 2)

    ## Using FC as value column: use value.var to override.

    GOheatmap

![](NMF_files/figure-markdown_strict/unnamed-chunk-7-1.png)
