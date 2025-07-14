
cols_samples <- c("79 yo"="#C6C6C6",
                  "116 yo"="#F90F59"
                  )

levels_samples <- c("79 yo", "116 yo")

cols_sample_type <- c("PB" = "#C20454", "BM" = "#6DA0F4")
cols_doublets <- c("Singlet" = "#C5C5C5", "Doublet" = "#F34444")


cols_celltype.l1 <- c("#8294D4", "#CFDF56", "#D4B7CE", "#74E483", "#DD6EBF", "#D3D4A3", "#A94BE0", "#D97C60", "#8BD7D3") %>% setNames(c("B","Mono","NK","CD8 T","CD4 T","DC","other T","other","HSPC"))


cols_celltype.l2 <- c("#E9EC9F", "#ABEACD", "#E3CFBA", "#64E4A0", "#E294B5", "#AA86B2", "#EBDA3C", "#91556A", "#67A552", "#E0B7EB", "#D95126", "#EC3BA9", "#F18B93",
  "#AEEEAC", "#9D883A", "#E1AC3B", "#EAC3D4", "#AA3EA5", "#4F9996", "#E86DBB", "#65EBEA", "#C0BAD2", "#609CDE", "#D32CE9", "#D7E374", "#9D9D89",
  "#6C5EDB", "#5EEA49", "#E2A991", "#5880E2", "#E83F74", "#E9E1F0", "#E3EBC0", "#A8B8F0", "#A2CCDA", "#B65F9B", "#7F9CB7", "#8DE66D", "#C85550",
  "#E08858", "#6F36E4", "#A78BDE", "#E89CEA", "#DF7CE9", "#9ABA82", "#65D8B9", "#B58B85", "#A9E9EE", "#4AEE8F", "#65C0E6", "#545E8B", "#B7EE41",
  "#DDEFE5", "#E4C07E", "#A448D7", "#E655E1") %>% setNames(
c("HSC","Memory B" ,"CD4 Memory","CD4 Naive" ,"Naive B" ,"NK" ,"transitional B","EMP","CD14 Mono", "pre B", "Plasma" ,"CD8 Effector_2" ,"MAIT", "CD8 Memory","ILC", "CD16 Mono" , "cDC2", "CD8 Naive",      
"Early Eryth","CD8 Effector_1","pDC" ,"Late Eryth" ,"Macrophage","pro B","pre-mDC" ,"GMP", "CD4 Effector","BaEoMa" ,  "LMPP", "Prog Mk" ,"T Proliferating" ,"NK CD56+","pre-pDC" ,"ASDC", "NK Proliferating","CLP" ,             
"cDC1" , "Platelet" ,"CD8 Effector_3","Stromal", "B intermediate" ,"CD8 TCM" ,"CD4 CTL", "CD4 TCM","CD8 TEM","Treg", "Plasmablast","gdT"   ,"B naive", "B memory" ,"CD4 TEM" ,"NK_CD56bright","dnT","CD8 Proliferating","HSPC" , "CD4 Proliferating"))



cols_celltype.l3 <- c("#AAECA3", "#59E944", "#E887E8", "#E83F74", "#65BB70", "#66EBEA", "#CE9D90", "#E86830", "#E6BF83", "#C6E771", "#B532E1", "#C14842", "#50EB9A",
  "#E9E7ED", "#EAD737", "#DDB7E9", "#A955DE", "#E0EAC3", "#ED3CAA", "#DA9EE6", "#EFDA6D", "#7CEB6F", "#E9ECA0", "#5C56CF", "#A0A9EE", "#A47699",
  "#E165B3", "#AC1EF7", "#B5ED40", "#ADEFCA", "#9E7FE3", "#66E8C1", "#EDC3D5", "#663BE4", "#B98B32", "#9AB7DC", "#447769", "#E82BE4", "#E156DC",
  "#C8BFD7", "#E4CDB7", "#BEE9E1", "#E18F69", "#5DAEC1", "#DF8FB3", "#99A95E", "#676D99", "#9FD7E6", "#B14EA3", "#6EC0A8", "#5FA7DE", "#EC8994",
  "#A1A093", "#4F81DB") %>% setNames(
c("HSPC",                  "B intermediate kappa" , "CD4 TCM_1",             "B naive kappa" ,        "NK_2" ,                
  "B intermediate lambda", "Treg Memory",           "CD4 TCM_3" ,            "CD14 Mono" ,            "B naive lambda" ,      
  "CD8 TEM_4" ,            "CD8 TEM_5" ,            "NK Proliferating" ,     "MAIT" ,                 "NK_1",                 
  "CD8 TEM_1" ,            "NK_CD56bright",         "CD16 Mono" ,            "CD8 TEM_2" ,            "Plasma",               
  "Plasmablast" ,          "CD8 Naive" ,            "CD8 TCM_2",             "CD4 CTL",               "CD4 TEM_2",            
  "CD8 TEM_6",             "CD4 TEM_1",             "pDC" ,                  "cDC2_1" ,               "NK_4" ,                
  "ILC",                   "CD4 TCM_2" ,            "Platelet" ,             "CD8 TCM_3",             "CD4 Naive",            
  "CD4 Proliferating",     "B memory lambda",       "CD8 Proliferating",     "NK_3",                  "cDC1" ,                
  "dnT_2" ,                "cDC2_2",                "ASDC_mDC",              "B memory kappa" ,       "Treg Naive",           
  "CD8 TCM_1" ,            "CD4 TEM_3",             "CD8 TEM_3",             "Eryth" ,                "CD8 Naive_2",          
  "gdT_4",                 "gdT_3",                 "gdT_1" ,                "gdT_2" ))


cols_diagnosis <- c("#8BCD98", "#5FB6DF", "#C28BCD") %>% setNames(c("AWM", "MGUS", "SWM"))


cols_majority_voting_Immune_All_High <- c("HSC or MPP" = "#EF072A",
                                          "B cells" = "#40D2E1",
                                          "T cells" = "#52A768",
                                          "ILC" = "#216332",
                                          "Early MK" = "#CB8AF9",
                                          "Monocytes" = "#DA8621",
                                          "B-cell lineage" = "#0C93BE",
                                          "Plasma cells" = "#135DA2",
                                          "Erythroid" = "#560334",
                                          "pDC" = "#F0EA76",
                                          "DC" = "#DBD323",
                                          "Megakaryocytes or platelets" = "#6812A5")


cols_majority_voting_Immune_All_Low <- c(
  "HSC or MPP" = "#FF000F", # Brighter red for clear visibility
  "CMP" = "#B21821", # Slightly more purple-toned red
  "MEMP" = "#F75660", # Pinkish red for differentiation
  "Megakaryocyte precursor" = "#BE93F9", # Light purple, distinguishable from mature megakaryocytes
  "Megakaryocytes or platelets" = "#8A2BE2", # Deeper purple for clear contrast
  "Early erythroid" = "#F693BF", # Darker burgundy for early stage
  "Mid erythroid" = "#BC4A7D", # Brighter red as transition
  "Late erythroid" = "#7A0D3E", # Rosy red for late stage, improved visibility
  "Classical monocytes" = "#FFB74B", # Bright orange, slightly lighter than non-classical
  "Non-classical monocytes" = "#E88D05", # Vivid orange for distinction
  "DC2" = "#FFF77F", # Olive green, contrast with pDC
  "pDC" = "#EEDF0A", # Bright yellow, for clear visibility
  "MAIT cells" = "#89BCAA", # Medium sea green, lighter than other T cells
  "CRTAM+ gamma-delta T cells" = "#84F0CB", # Sea green for visibility
  "CD16+ NK cells" = "#A25D13", # Forest green, slightly darker for ILC
  "CD16- NK cells" = "#754007", # Dark green, ensuring visibility
  "Regulatory T cells" = "#B8FEB8", # Dark green for contrast
  "Tcm or Naive helper T cells" = "#79E579", # Different shade of green from MAIT
  "Tem or Effector helper T cells" = "#11B911", # Brighter green, for distinction
  "Tcm or Naive cytotoxic T cells" = "#057E05", # Office green, distinguishable among T cells
  "Tem or Temra cytotoxic T cells" = "#075F07", # Dark green, for cytotoxic lineage visibility
  "Tem or Trm cytotoxic T cells" = "#105010", # Very dark green, maximum contrast in T cells
  "ILC3" = "#194D33", # Dark slate green, standing out within the ILC group
  "Pro-B cells" = "#7FEADD", # Lighter blue for early B-cell stage
  "Small pre-B cells" = "#1CDCC4", # Medium blue, maintaining the gradient
  "Large pre-B cells" = "#0ABEA8", # As provided, solid blue
  "Age-associated B cells" = "firebrick1", # Darker blue for aged B cells
  "Naive B cells" = "#58B2F5", # Deep blue for naive B cells
  "Memory B cells" = "#008EF7", # Turquoise, as provided
  "Plasma cells" = "#0E2E65", # As provided, dark blue
  "Plasmablasts" = "#414F68" # Dark blue, close to plasma cells but distinguishable
  
)

