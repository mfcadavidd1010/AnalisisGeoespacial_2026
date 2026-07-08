library(INLA)
> library(sf)
> library(spdep)
> # Y re-ejecutar la lectura + vecindad
  > malla <- st_read("reticula_bym.gpkg")
> vecinos <- poly2nb(malla, queen = TRUE)
> nb2INLA("malla.adj", vecinos)
> g <- inla.read.graph("malla.adj")
> cat("Vecinos promedio:", mean(card(vecinos)), "\n")        # ~7.66
> cat("Celdas sin vecinos:", sum(card(vecinos) == 0), "\n")  # 0
> malla$lito_grupo <- relevel(factor(malla$lito_grupo), ref = "Volcanica")
> malla$idx <- 1:nrow(malla)
> print(levels(malla$lito_grupo))
> table(malla$lito_grupo)
> # Fórmula del BYM2: litología (efectos fijos) + efecto espacial bym2
> formula_bym <- n_au ~ 1 + lito_grupo +
+     f(idx, model = "bym2", graph = g, scale.model = TRUE,
+       hyper = list(
+           phi = list(prior = "pc", param = c(0.5, 0.5)),      # prior sobre la mezcla espacial/iid
+           prec = list(prior = "pc.prec", param = c(1, 0.01))  # prior sobre la precisión
+       ))
> 
> # Ajustar el modelo
> modelo_bym <- inla(
+     formula_bym,
+     family = "poisson",
+     data = malla,
+     control.compute = list(dic = TRUE, waic = TRUE),   # métricas de ajuste
+     control.predictor = list(compute = TRUE)            # para obtener las predicciones
+ )
> 
> cat("¡Modelo ajustado! Tiempo:", modelo_bym$cpu.used["Total"], "segundos\n")
> summary(modelo_bym)
Time used:
    Pre = 21.3, Running = 3.75, Post = 0.391, Total = 25.4 
Fixed effects:
                               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
(Intercept)                  -2.156 0.156     -2.464   -2.155     -1.853 -2.155   0
lito_grupoIntrusiva           0.746 0.185      0.383    0.746      1.109  0.746   0
lito_grupoMafico_Ultramafica  0.154 0.201     -0.240    0.154      0.549  0.154   0
lito_grupoMetamorfica         0.813 0.201      0.419    0.813      1.207  0.813   0
lito_grupoSed_Cobertura       0.307 0.243     -0.169    0.307      0.783  0.307   0

Random effects:
  Name	  Model
    idx BYM2 model

Model hyperparameters:
                   mean    sd 0.025quant 0.5quant 0.975quant  mode
Precision for idx 0.285 0.021      0.246    0.285      0.329 0.283
Phi for idx       0.991 0.009      0.967    0.994      0.999 0.998

Deviance Information Criterion (DIC) ...............: 4351.77
Deviance Information Criterion (DIC, saturated) ....: 2477.20
Effective number of parameters .....................: 767.44

Watanabe-Akaike information criterion (WAIC) ...: 4352.12
Effective number of parameters .................: 565.27

Marginal log-Likelihood:  -1431.34 
 is computed 
Posterior summaries for the linear predictor and the fitted values are computed
(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

> # 1. El efecto espacial estructurado (las primeras n filas del campo bym2)
> n_celdas <- nrow(malla)
> malla$efecto_espacial <- modelo_bym$summary.random$idx$mean[1:n_celdas]
> 
> # 2. El valor ajustado por celda (intensidad esperada de Au) = mapa de prospectividad final
> malla$prospectividad <- modelo_bym$summary.fitted.values$mean[1:n_celdas]
> 
> # 3. Guardar para mapear (en Python o donde prefieras)
> st_write(malla, "reticula_bym_resultados.gpkg", delete_dsn = TRUE)
Deleting source `reticula_bym_resultados.gpkg' failed
Writing layer `reticula_bym_resultados' to data source 
  `reticula_bym_resultados.gpkg' using driver `GPKG'
Writing 2548 features with 6 fields and geometry type Polygon.
> cat("Resultados guardados. Resumen de prospectividad:\n")
Resultados guardados. Resumen de prospectividad:
> summary(malla$prospectividad)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
 0.00185  0.13747  0.29118  0.78322  0.69255 26.02294 
> # 4. Mapa rápido en R para ver el resultado de una vez
> library(ggplot2)
> ggplot(malla) +
+     geom_sf(aes(fill = prospectividad), color = NA) +
+     scale_fill_viridis_c(option = "inferno") +
+     labs(title = "Mapa de prospectividad de Au (BYM2)",
+          fill = "Intensidad\nesperada") +
+     theme_minimal()
> library(ggplot2)
> 
> ggplot(malla) +
+     geom_sf(aes(fill = prospectividad), color = NA) +
+     scale_fill_viridis_c(option = "inferno", trans = "log1p",   # <- escala log
+                          breaks = c(0, 0.5, 1, 2, 5, 10, 25)) +
+     labs(title = "Mapa de prospectividad de Au (BYM2)",
+          subtitle = "Intensidad esperada de ocurrencias por celda (escala log)",
+          fill = "Intensidad\nesperada") +
+     theme_minimal()
> ggplot(malla) +
+     geom_sf(aes(fill = efecto_espacial), color = NA) +
+     scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
+     labs(title = "Efecto espacial latente (BYM2)",
+          subtitle = "Prospectividad no explicada por la litología",
+          fill = "Efecto\nespacial") +
+     theme_minimal()
> # Residuos del BYM2 = observado - ajustado
> malla$resid_bym <- malla$n_au - modelo_bym$summary.fitted.values$mean[1:nrow(malla)]
> 
> # Moran sobre los residuos, con la MISMA vecindad (vecinos)
> library(spdep)
> pesos <- nb2listw(vecinos, style = "W")          # W estandarizada por filas
> moran_resultado <- moran.test(malla$resid_bym, pesos)
> print(moran_resultado)

	Moran I test under randomisation

data:  malla$resid_bym  
weights: pesos    

Moran I statistic standard deviate = -5.531, p-value = 1
alternative hypothesis: greater
sample estimates:
Moran I statistic       Expectation          Variance 
    -0.0567750999     -0.0003926188      0.0001039141 

> # 1. Guarda el modelo ajustado en disco (si no lo hiciste antes)
> saveRDS(modelo_bym, "modelo_bym.rds")
> 
> # 2. Asegúrate de tener tu código en un SCRIPT guardado (.R), no solo en la consola