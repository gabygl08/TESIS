library(dplyr)
library(mixOmics)

gaby_function_PLSDA <- function(X, Y, max_, step_) {
  ## -------------------------------------------------------------------------------------------------------------------
  srbct.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later
  
  ## ---- fig.cap = "FIGURE 4: Tuning the number of components in PLS-DA on the SRBCT gene expression data. For each component, repeated cross-validation (10 × 3−fold CV) is used to evaluate the PLS-DA classification performance (OER and BER), for each type of prediction distance; `max.dist`, `centroids.dist` and `mahalanobis.dist`."----
  # undergo performance evaluation in order to tune the number of components to use
  perf.splsda.srbct <- perf(
    srbct.splsda,
    validation = "Mfold",
    folds = 5,
    nrepeat = 10,
    # use repeated cross-validation
    progressBar = FALSE,
    auc = TRUE
  ) # include AUC values
  
  # plot the outcome of performance evaluation across all ten components
  plot(perf.splsda.srbct)
  
  ## -------------------------------------------------------------------------------------------------------------------
  perf.splsda.srbct$choice.ncomp # what is the optimal value of components according to perf()
  
  ## ---- fig.cap = "FIGURE 5:  Tuning keepX for the sPLS-DA performed on the SRBCT gene expression data. Each coloured line represents the balanced error rate (y-axis) per component across all tested keepX values (x-axis) with the standard deviation based on the repeated cross-validation folds. As sPLS-DA is an iterative algorithm, values represented for a given component (e.g. comp 1 to 2) include the optimal keepX value chosen for the previous component (comp 1)."----
  # grid of possible keepX values that will be tested for each component
  
  ### SI EN PERFOREMCNE TEN SALTAN MAS DE 2 COPONENTES TENDRAS QUE ONER UNA LISTA poniendo tantos vectores como ponentes haya
  
  
  # undergo the tuning process to determine the optimal number of variables
  tune.splsda.srbct <- tune.splsda(
    X,
    Y,
    ncomp = 4,
    # calculate for first 4 components
    validation = 'Mfold',
    folds = 5,
    nrepeat = 10,
    # use repeated cross-validation
    dist = 'max.dist',
    # use max.dist measure
    measure = "BER",
    # use balanced error rate of dist measure
    test.keepX = seq(1, max_, step_),
    cpus = 2
  ) # allow for paralleliation to decrease runtime
  
  ## -------------------------------------------------------------------------------------------------------------------
  optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
  optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]
  
  ## -------------------------------------------------------------------------------------------------------------------
  # form final model with optimised values for component and variable count
  final.splsda <- splsda(X, Y, ncomp = optimal.ncomp, keepX = optimal.keepX)
  return(final.splsda)
}

gaby_function_PLS <- function(X, Y, min_, max_, step_) {
  ## -------------------------------------------------------------------------------------------------------------------
  spls.liver <- spls(
    X = X,
    Y = Y,
    ncomp = 10,
    mode = 'regression'
  )
  
  
  
  tune.spls.liver <- tune.spls(
    X,
    Y,
    ncomp = 2 ,
    test.keepX = seq(min_, max_, step_),
    test.keepY = 1:2,
    nrepeat = 1,
    folds = 10,
    # use 10 folds
    mode = 'regression',
    measure = 'cor'
  )
  
  optimal.keepX <- tune.spls.liver$choice.keepX
  
  # extract optimal number of variables for Y datafram
  optimal.keepY <- tune.spls.liver$choice.keepY
  
  optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components
  
  ## -------------------------------------------------------------------------------------------------------------------
  # form final model with optimised values for component and variable count
  final.spls.liver <- spls(
    X,
    Y,
    ncomp = 2,
    keepX = optimal.keepX,
    keepY = optimal.keepY,
    mode = "regression"
  ) # explanitory approach being used,
  # hence use regression mode
  
  # evaluate the prediction accuracy for the first two components
  
  return(
    list(
      modelo = final.spls.liver,
      comp = optimal.ncomp,
      optimal.keepX = optimal.keepX,
      optimal.keepY = optimal.keepY
    )
  )
  
  
}
#_--------------------------------------------------------
preprocess_data <- function(A) {
  set.seed(123)
  idx <- sample(1:nrow(A), nrow(A) * 0.7)
  train <- A[idx, ]
  test <- A[-idx, ]
  X <- read.csv("./volumenes_eppendorf.csv", sep = ";") ## leemos y eliminamos la primera columna q la anade el R
  muestras <- apply(X[, 1:2], 1, function(x)
    paste0(paste0(x, collapse = "."), collapse = ".CSV"))
  muestras <- paste0(muestras, ".CSV")
  head(X)
  X$Codigo <- muestras
  Y <- X[, -2]
  A$Codigo  <- rownames(A)
  dats <- inner_join(Y, A, by = "Codigo")
  rownames(dats) <- dats$Codigo
  dats_ <- dats[, -1]
  list(X_ = dats_[, 3:ncol(dats_)], Y_ = dats_[, 1:2])
}


## Primero esto con limpio
#
# limpio2<-bind_rows(etanol_df2,acetico_df2)
# Xpda<-as.data.frame(limpio2)
# Ypda <- condicion
# idxpda <- sample(1:nrow(Xpda), nrow(Xpda) * 0.7)
# length(idxpda)
# train.Xpda <- Xpda[idxpda, ]
# train.Ypda <- Ypda[idxpda]
# test.Xplsda <- Xpda[-idxpda, ]
# test.Yplsda<- Ypda[-idxpda]
# condicion
# length(condicion)
#
#
# length(colnames(limpio2))
#
# modelo_plsda <- gaby_function_PLSDA(train.Xpda,train.Ypda,max_ =7469,step_ = 1000 ) # para limpio
#
# modelo_train_plsda <- splsda(X = train.Xpda, Y = train.Ypda, ncomp = modelo_plsda$ncomp, keepX = modelo_plsda$optimal.keepX)
#
# predict.splsda.srbct <- predict(modelo_train_plsda, test.Xplsda, dist = "mahalanobis.dist")
# predict.comp1 <- predict.splsda.srbct$class$mahalanobis.dist[,1]
# predict.comp1.fact<-as.factor(predict.comp1)
# table(factor(predict.comp1.fact, levels = levels(Ypda)), test.Yplsda)
#
# auc.splsda = auroc(modelo_train_plsda, roc.comp = 1, print = FALSE) # AUROC for the first component
# #ROC 100% - justificar la respuesta . diseño experimental analitos por separado


#CONTAMINADA
vRAa <- readRDS("./variables.rds")
Y <- vRAa[[2]]
muestracruda <- vRAa[[1]]
set.seed(123)
# muestracruda.sin<-as.matrix(muestracruda)
# idx2.sin <- sample(1:nrow(muestracruda.sin),nrow(muestracruda.sin)*0.6)
# train2.sin <- muestracruda.sin[idx2.sin,]
# test2.sin <- muestracruda.sin[-idx2.sin,]
# muestracruda.sin_df<-as.data.frame(muestracruda.sin)
# dim(muestracruda.sin_df)
# colnames(muestracruda.sin_df) <- colnames(muestracruda)
# rownames(muestracruda.sin_df)<- rownames(muestracruda)
muestracruda$Codigo  <- rownames(muestracruda)

dats2.sin <- inner_join(Y, muestracruda, by = "Codigo")
rownames(dats2.sin) <- dats2.sin$Codigo
dats_2.sin <- dats2.sin[, -1]
idx2.sin <- sample(1:nrow(dats_2.sin), nrow(dats_2.sin) * 0.6)

X_2.sin <- dats_2.sin[, 3:ncol(dats_2.sin)]
Y_2.sin <- dats_2.sin[, 1:2]
train.X2.sin <- X_2.sin[idx2.sin, ]
train.Y2.sin <- Y_2.sin[idx2.sin, ]
test.X2.sin <- X_2.sin[-idx2.sin, ]
test.Y2.sin <- Y_2.sin[-idx2.sin, ]
n1.sin <- nrow(test.X2.sin)

modelo_pls <- gaby_function_PLS(
  train.X2.sin,
  train.Y2.sin,
  min_ = 2000,
  max_ = 3000,
  step_ = 100
)

modelo_train_pls <- spls(
  train.X2.sin,
  train.Y2.sin ,
  ncomp = modelo_pls$comp,
  keepX = modelo_pls$optimal.keepX,
  keepY = modelo_pls$optimal.keepY
)
prediccion.final <-  predict(modelo_pls$modelo, newdata = test.X2.sin)

funrmse <- function(pred, actual) {
  n <- length(pred)
  rmse <- sqrt(sum((pred - actual)^2) / n)
  
  
  return(rmse)
}

Comp1 <- prediccion.final$predict[, 1:2, 1]
Comp2 <- prediccion.final$predict[, 1:2, 2]

et <- Comp1[, 2]
ac <- Comp1[, 1]
et2 <- Comp2[, 2]
ac2 <- Comp2[, 1]

SME_comp1_et <- funrmse(et, test.Y2.sin[, 2])
SME_comp1_et

SME_comp1_ac <- funrmse(ac, test.Y2.sin[, 1])
SME_comp1_ac

SME_comp2_et <- funrmse(et2, test.Y2.sin[, 2])
SME_comp2_et

SME_comp2_ac <- funrmse(ac2, test.Y2.sin[, 1])
SME_comp2_ac


minimo <- seq(1000, 2000, 100)
maximo <- seq(3000, 4000, 100)

lista <- list()
lista.nombres <- list()
i <- 1
for (m in seq_along(minimo)) {
  for (ma in seq_along(maximo)) {
    lista.nombres[[i]] <- paste(minimo[m], maximo[ma], sep = "_")
    lista[[i]] <- gaby_function_PLS(
      train.X2.sin,
      train.Y2.sin,
      min_ = minimo[m],
      max_ = maximo[ma],
      step_ = 50
    )
    i <- i * 1
    
  }
  
}

names(lista) <- unlist(lista.nombres)


lista.modelos <- lapply(lista, function(x) x$modelo)
names(lista.modelos) <- unlist(lista.modelos)

fun_extract_bst_model <- function(X,test_data,test_ti){
  
  prediccion.final <-  predict(X, newdata = test_data)
  
  Comp1 <- prediccion.final$predict[, 1:2, 1]
  Comp2 <- prediccion.final$predict[, 1:2, 2]
  
  et <- Comp1[, 2]
  ac <- Comp1[, 1]
  et2 <- Comp2[, 2]
  ac2 <- Comp2[, 1]
  
  SME_comp1_et <- funrmse(et, test_ti[, 2])

  SME_comp1_ac <- funrmse(ac, test_ti[, 1])
  
  
  SME_comp2_et <- funrmse(et2, test_ti[, 2])
  
  
  SME_comp2_ac <- funrmse(ac2, test_ti[, 1])
  
  
  
  vector_respuesta <- c(et_1 = SME_comp1_et,
                        et_2 = SME_comp2_et,
                        ac_1 = SME_comp1_ac,
                        ac_2 = SME_comp2_ac)
  
  return(vector_respuesta)
  
  
  
}




respuesta_final_tesis <- lapply(lista.modelos,function(x) fun_extract_bst_model(x,test_data =test.X2.sin,test_ti = test.Y2.sin))
respuesta_final_tesis <- as.data.frame(Reduce("rbind",respuesta_final_tesis))
colnames(respuesta_final_tesis) <- unlist(lista.nombres)
rownames(respuesta_final_tesis) <- c("et1","et2","ac1","ac2")# colocar el nombre de las columnas 
  
#pasar la mejor comb. pls
modelo_pls_final <- lista.modelos$`2000_4000`

#best_model_para_Etanol <- rownames(respuesta_final_tesis)[which.min(respuesta_final_tesis[,2])]
#print(best_model_para_Etanol)

#modelo_final <- modelos.lista[[best_model_para_Etanol]]

#COMPROBAR
  fun_to_spectra <- function(X) {
    spectra <- as.data.frame((ir_flatten(X)))
    wn <- spectra[, 1]
    spectra <- spectra[, -1]
    rownames(spectra) <- wn
    colnames(spectra) <- X$sample_id
    spectra.t <- t(spectra)
    return(spectra.t)
  }
 
  
#MUESTRAS:
  archivos <- list.files("./MUESTRAS VINAGRE/1/",pattern = "\\.CSV$", full.names = T)
  archivos 
  milista <-  lapply(archivos, function(x)
    read.csv(x, sep = ";", header = F))
  archivos_names <- list.files("./MUESTRAS VINAGRE/1/")
  archivos_names
  for (l in seq_along(milista)) {
    write.csv(milista[[l]],
              file = file.path("./MUESTRAS VINAGRE/", paste0(archivos_names[l], ".csv")),
              row.names = F)
  }
  
  PRIMERADERIVADA<-function(X){X %>%  
      ir_bc(method = "rubberband") %>%
      ir_normalise()  %>%
      ir_smooth(method = "sg", p = 3, n = 91, m = 0)%>%   
      ir_smooth(method = "sg", p = 3, n = 9, m = 1) }
  
  SEGUNDADERIVADA_ir<-function(X,p){
    r<-fun_to_spectra(X)
    w<-as.data.frame(r)
    m<-as.matrix(w)
    baseline <- baseline(m, method = "als")
    baseline <- baseline@corrected
    baseline[baseline<0]<-0
    cont <- ir::ir_import_csv(p)
    contaminado_ir<-fun_to_ir(baseline,cont)
    suavizado_ir<-contaminado_ir%>%ir_smooth(method = "sg", p = 3, n = 91, m = 0)
    suavizado_ir.matriz<-fun_to_spectra(suavizado_ir)
    normalized_spectra <- t(apply(suavizado_ir.matriz, 1, function(row) row / sum(row)))
    datos_bc.filtered.normalized <- fun_to_ir(s = normalized_spectra,X = cont)
    datos_bc.spectra.filtered <-  as.data.frame(t(apply(normalized_spectra,1,function(x) sgolayfilt(x, p = 2, n = 7, m = 2))))
    colnames(datos_bc.spectra.filtered) <- colnames(suavizado_ir.matriz)
    datos_ir<-fun_to_ir(datos_bc.spectra.filtered,cont)
    return(datos_bc.spectra.filtered)
  }

madre<-ir_import_csv("./MUESTRAS VINAGRE/2/madre.CSV.csv",sep=",")
comercial<-ir_import_csv("./MUESTRAS VINAGRE/2/comercial.CSV.csv",sep=",")  
VINAGRE <- list.files("./MUESTRAS VINAGRE/2/",full.names = T)
MUESTRASVINAGRE<-bind_rows(madre,comercial)

MUESTRASVINAGRE1<-PRIMERADERIVADA(MUESTRASVINAGRE)
MUESTRASVINAGRE2<-SEGUNDADERIVADA_ir(MUESTRASVINAGRE1,VINAGRE)



prediccion.comp <-  predict(modelo_pls_final, newdata = MUESTRASVINAGRE2)



Comp1 <- prediccion.comp$predict[, 1:2, 1]
Comp2 <- prediccion.comp$predict[, 1:2, 2]

et <- Comp1[, 2]
et
ac <- Comp1[, 1]
ac
et2 <- Comp2[, 2]
et2
ac2 <- Comp2[, 1]
ac2

tabla_resultados <- data.frame(
  Analito = c("Etanol", "Ácido Acético", "Etanol", "Ácido Acético"),
  Componente = c("Componente 1", "Componente 1", "Componente 2", "Componente 2"),
  Madre = c(Comp1[1, 2], Comp1[1, 1], Comp2[1, 2], Comp2[1, 1]),
  Comercial = c(Comp1[2, 2], Comp1[2, 1], Comp2[2, 2], Comp2[2, 1])
)
tabla_resultados

library(knitr)
kable(tabla_resultados)

#1er elemento de ac2 - va a ser madre y el  segundo el comercial



