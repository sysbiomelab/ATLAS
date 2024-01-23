

source("~/rscripts/2019_funcType/funcType.analysis.functions.r")
source("~/rscripts/2019_funcType/funcType.data.files.r")




funcGemJaccNetMode = F
if (funcGemJaccNetMode) {
  #funcGemModNetNodeTableFile = "~/sunjae.mAtlas.func.annotation/funcGemModeNet.nodeTable.20190828.txt"
  #funcGemModNetEdgeTableFile = "~/sunjae.mAtlas.func.annotation/funcGemModeNet.edgeTable.20190828.txt"
  
  funcGemModNetNodeTableFile = "~/sunjae.mAtlas.func.annotation/funcGemModeNet.nodeTable.20191227.txt"
  funcGemModNetEdgeTableFile = "~/sunjae.mAtlas.func.annotation/funcGemModeNet.edgeTable.20191227.txt"
  
  load(funcGemJaccMatRData)
  funcGemJaccNet = makeJaccNet(funcGemJaccMat)
  save(funcGemJaccNet, file = funcGemJaccNetRData)
  
  funcGemModules = makeModuleList(funcGemJaccNet)
  save(funcGemModules, file = funcGemModulesRData)
  
  funcGemModuleAnnots = annotateModulesByCC(funcGemJaccNet, funcGemModules, 3)
  writeModuleNetworkAnnot(funcGemModuleAnnots, funcGemModNetNodeTableFile, funcGemModNetEdgeTableFile)
  save(funcGemModuleAnnots, file=funcGemModuleAnnotsRData)
  
}

rm(list=ls())

source("~/rscripts/2019_funcType/funcType.analysis.functions.r")
source("~/rscripts/2019_funcType/funcType.data.files.r")

funcJaccNetMode = T
if (funcJaccNetMode) {
  #funcModNetNodeTableFile = "~/sunjae.mAtlas.func.annotation/funcModeNet.nodeTable.20190828.txt"
  funcModNetNodeTableFile = "~/sunjae.mAtlas.func.annotation/funcModeNet.nodeTable.20191227.txt"
  #funcModNetEdgeTableFile = "~/sunjae.mAtlas.func.annotation/funcModeNet.edgeTable.20190828.txt"
  funcModNetEdgeTableFile = "~/sunjae.mAtlas.func.annotation/funcModeNet.edgeTable.20191227.txt"
  
  load(funcJaccMatRData)
  
  if (T) {
    funcJaccNet = makeJaccNet(funcJaccMat)
    save(funcJaccNet, file = funcJaccNetRData)
    
    funcModules = makeModuleList(funcJaccNet)
    save(funcModules, file = funcModulesRData)
  }
  if (F) {
    load(funcJaccNetRData)
    load(funcModulesRData)
  }
  
  funcModuleAnnots = annotateModulesByCC(funcJaccNet, funcModules, 3)
  writeModuleNetworkAnnot(funcModuleAnnots, funcModNetNodeTableFile, funcModNetEdgeTableFile)
  save(funcModuleAnnots, file=funcModuleAnnotsRData)
  
}

rm(list=ls())
