M14 <- DATAMUESTRAS1 %>% filter(grepl("^M14\\.", sample_id))
plot(M14)+geom_path(aes(color = sample_id ))+
  labs(x = "Longitud de Onda (cm⁻¹)", y = "Absorbancia")

plot(DATAMUESTRAS1)+ geom_path(aes(color = sample_id ))+
  labs(x = "Longitud de Onda (cm⁻¹)", y = "Absorbancia") +  # Agrega títulos a los ejes
  theme(legend.position = "none")

#GRAFICOS
#CORRECCIÓN DE LA LINEA BASE:
contaminado_mx<-as.matrix(to_spectra(sample))
Correcionbbase2<-baseline(contaminado_mx, method = "als")
Correcionbbase2 <- Correcionbbase2@corrected
Correcionbbase2[Correcionbbase2<0]<-0
cont <- ir::ir_import_csv(MUESTRAS)
contaminado_ir<-fun_to_ir(Correcionbbase2,cont)
plot(contaminado_ir) + 
  geom_path(aes(color = sample_id)) +
  labs(x = "Longitud de Onda (cm⁻¹)", y = "Absorbancia") +  # Agrega títulos a los ejes
  theme(legend.position = "none")

#NORMALIZADO


normalized_spectra <- t(apply(Correcionbbase2, 1, function(row) row / sum(row)))
normalizado2 <- fun_to_ir(s = normalized_spectra,X = cont)
plot(normalizado2)+geom_path(aes(color = sample_id ))+
  labs(x = "Longitud de Onda (cm⁻¹)", y = "Absorbancia") +  # Agrega títulos a los ejes
  theme(legend.position = "none")

#SUAVIZADO
suavizado_ir<-normalizado2%>%ir_smooth(method = "sg", p = 3, n = 91, m = 0)
plot(suavizado_ir)+geom_path(aes(color = sample_id ))+
  labs(x = "Longitud de Onda (cm⁻¹)", y = "Absorbancia") +  # Agrega títulos a los ejes
  theme(legend.position = "none")
suavizado_ir.matriz<-fun_to_spectra(suavizado_ir)


#SEGUNDA
SEGUNDA <-  as.data.frame(t(apply(suavizado_ir.matriz,1,function(x) sgolayfilt(x, p = 2, n = 7, m = 2))))
colnames(SEGUNDA) <- colnames(suavizado_ir.matriz)
segderivada<-fun_to_ir(SEGUNDA,cont)
plot(segderivada)+geom_path(aes(color = sample_id ))+
  labs(x = "Longitud de Onda (cm⁻¹)", y = "Absorbancia") +  # Agrega títulos a los ejes
  theme(legend.position = "none")

#ESCALADA
SEGUNDA.MX<-as.matrix(SEGUNDA)
SEGUNDA_ESCALADA <- as.data.frame(scale(SEGUNDA.MX))
segderivada_escalada <- fun_to_ir(SEGUNDA_ESCALADA, cont)
plot(segderivada_escalada) + 
  geom_path(aes(color = sample_id)) +
  labs(x = "Longitud de Onda (cm⁻¹)", y = "Absorbancia Normalizada") +
  theme(legend.position = "none")


#MUESTRAS DE VINAGRE

MUESTRASVINAGRE2<-SEGUNDADERIVADA_ir(MUESTRASVINAGRE1,VINAGRE)
MUESTRASVINAGRE2_long <- MUESTRASVINAGRE2 %>%
  mutate(sample = c("Madre", "Comercial")) %>%  # Etiquetar las dos muestras
  pivot_longer(cols = -sample, names_to = "Wavelength", values_to = "Intensity")
MUESTRASVINAGRE2_long$Wavelength <- as.numeric(MUESTRASVINAGRE2_long$Wavelength)
ggplot(MUESTRASVINAGRE2_long, aes(x = Wavelength, y = Intensity, color = sample)) +
  geom_line(size = 1) +
  labs(title = "Espectros de Segunda Derivada",
       x = "Longitud de Onda (cm⁻¹)",
       y = "Segunda Derivada de la Absorbancia",
       color = "Muestra") +
  theme_minimal()



#----
  M14_VINA<- rbind(M14,MUESTRASVINAGRE)
plot(M14_VINA)+geom_path(aes(color = sample_id ))+
  labs(x = "Longitud de Onda (cm⁻¹)", y = "Absorbancia")

#----
#PLS
summary(modelo_pls_final)
modelo_pls_final
plotIndiv(modelo_pls_final)
plotArrow(modelo_pls_final)
plotLoadings(modelo_pls_final)

#-------
prediccion.comp$predict



