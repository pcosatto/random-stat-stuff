#Librerias varias que vamos a usar

librerias <- c("formatR", "latex2exp", "scales", "emmeans",
               "splines", "mgcv", "modelr", "vcdExtra",
               "pROC", "cfcdae", "latex2exp", "nlme", "lme4",
               "fda")

# Comprobar si las librerías están instaladas y si no, instalarlas
for (libreria in librerias) {
  if (!requireNamespace(libreria, quietly = TRUE)) {
    install.packages(libreria, dependencies = TRUE)
  }
  
  # Cargar la librería
  library(libreria, character.only = TRUE)
}

#Parametros graficos (base)
graph_par <- function(mfrow = NULL){
  if(is.null(mfrow)){mfrow <- c(1,1)}
  par(family = 'sans', cex.axis=0.9, cex.lab=1, 
      cex.main = 1.1, 
      mar=c(4,4,4,3) - 1.5,
      mgp=c(1.1,0.25,0), tcl=0, mfrow = mfrow)
}

#Scatter de pares
mv.pairs <- function(data, col=NULL){
  panel.hist <- function(x, ...)
  {
    usr <- par("usr")
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "lightgreen", ...)
  }
  
  panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
  {
    r <- cor(x, y)
    txt <- round(r,2)
    text(mean(range(x)),mean(range(y)), txt, cex = 2,col="gray40", font=2)
  }
  
  if(is.null(col)){
    col <- scales::alpha('black',0.3)
  } else {
    col <- as.factor(col)
    palette <- paletteer_c("grDevices::rainbow", nlevels(col))
    col <- palette[col]
  }
  
  panel.points <- function(x, y, ...) {
    points(x, y, pch = 20, 
           col = col, cex = 1)
  }
  
  pairs(data, lower.panel = panel.points,
        upper.panel= panel.cor, pch=20,
        gap = 0.2, diag.panel = panel.hist)
}


#Plot tridimensional (con Plotly)
mv.plot3d <- function(X,col=NULL,id=NULL, size=6, cube=FALSE){
  
  library(plotly)
  library(RColorBrewer)
  
  n <- nrow(X)
  p <- ncol(X)
  
  data <- data.frame(scale(X, scale=FALSE))
  names(data) <- c('x1','x2','x3')
  
  if(is.null(col)==TRUE){
    data$col <- rep('black',n)
  } else {
    data$col <- col}
  
  if(is.null(id)==TRUE){
    data$id<-1:n
  } else {data$id <- id}
  
  fig <- plot_ly(data,
                 x = ~data[,1], y = ~data[,2], z = ~data[,3],
                 colors = brewer.pal(p,'Set1'), text=~id,
                 marker = list(size=size))
  fig <- fig %>% add_markers(color = ~col)
  fig <- fig %>% layout(scene = list(xaxis = list(title = colnames(X)[1],
                                                  range = c(min(data$x1),max(data$x1))),
                                     yaxis = list(title = colnames(X)[2],
                                                  range = c(min(data$x2),max(data$x2))),
                                     zaxis = list(title = colnames(X)[3],
                                                  range = c(min(data$x3),max(data$x3))),
                                     aspectmode = ifelse(cube==TRUE,'cube','auto')))
  fig
  
}



#Correlograma a color
mv.correlogram <- function(data){
  library(ellipse)
  par(mfrow=c(1,1))
  
  cor <- cor(data)
  
  # Build a Pannel of 100 colors with Rcolor Brewer
  my_colors <- RColorBrewer::brewer.pal(5, "Spectral")
  my_colors <- colorRampPalette(my_colors)(100)
  
  # Order the correlation matrix
  ord <- order(cor[1, ])
  data_ord <- cor[ord, ord]
  plotcorr(data_ord , 
           col=my_colors[data_ord*50+50] , 
           mar=c(1,1,1,1), type='upper')
}

#Distribuciones marginales (con opción a curvas de nivel)
mv.scatter <- function(x, y, contour = FALSE) {
  # library
  library(ggplot2)
  library(ggExtra)
  
  df <- data.frame(x, y)
  
  # Classic plot :
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(color = "gray70") +
    theme_minimal() +
    theme(panel.background = element_blank(),  # Remove gray background
          panel.grid.major = element_blank(),   # Remove grid lines
          panel.grid.minor = element_blank(),   # Remove grid lines
          legend.position = "none") +
    xlab(deparse1(substitute(x))) + 
    ylab(deparse1(substitute(y)))
  
  # Add density contour lines if contour is TRUE
  if (contour) {
    p <- p + geom_density_2d(color = "forestgreen") 
  }
  
  # with marginal boxplot
  p1 <- ggMarginal(p, type = 'boxplot', margins = 'both', 
                   fill = "lightgreen") 
  
  p1
}

#Scatter con gradiente de colores equidistantes
mv.equiplot <- function(x, y, distance, p = NULL, robust = FALSE, ...) {
  
  #Graphical settings
  graph_par <- function(){
    par(family = "Verdana", cex.axis=0.7, cex.lab=0.7, mar=c(4,4,2,3) - 1.5,
        mgp=c(1.1,0.25,0), tcl=0)
  }
  
  #Function to extract color gradient
  equidistance_gradient <- function(x, y, distance = c('minkowski', 'mahalanobis'), p = NULL, robust = FALSE) {
    
    X <- cbind(x, y)
    
    # Definition of palette
    palette <- colorRampPalette(c('palegreen1', 'palegreen3', 'khaki2', 'orange', 'firebrick2', 'maroon3', 'purple'))
    
    # Mean and covariance calculation (classical or robust)
    if (robust == TRUE) {
      covrob <- MASS::cov.rob(X)
      mean <- covrob$center
      covar <- covrob$cov
    } else { 
      mean <- colMeans(X) 
      covar <- cov(X)
    }
    
    # Calculation of distance
    switch(distance, minkowski = {
      d <- apply(X, 1, function(x) {
        sum(abs(x - mean)^p)^(1/p)
      })
    },
    mahalanobis = {
      d <- sqrt(mahalanobis(X, center = mean, cov = covar))
    })
    
    # Return color scale
    index <- round((d - min(d)) / (max(d) - min(d)) * 20 + 1)
    return(list('colors'=palette(21)[index],
                'index'=index))
  }
  
  gradient <- equidistance_gradient(x, y, distance, p, robust)
  
  #Main plot
  graph_par()
  plot(x, y, pch = 20, col = gradient$colors, 
       xlab = deparse1(substitute(x)), 
       ylab = deparse1(substitute(y)),...)
  
  outliers <-  seq_along(x)[gradient$index > quantile(gradient$index, 0.95)]
  text(x[outliers], y[outliers], outliers, pos = 2,
       offset = 0.2, cex = 0.6, col = 'grey')
  
  
}

#Box-cox multivariada aplicada a una matriz de datos
mv.boxcox <- function(X){
  power <- powerTransform(X)
  return(bcPower(X,power$lambda)) 
}

#Biplot
mv.biplot <- function(puntos, direcciones,
                      color_puntos = NULL,
                      id_puntos = NULL,
                      color_ejes = NULL){
  
  n <- nrow(puntos)
  p <- nrow(direcciones)
  
  graph_par <- function(){
    par(family = "Verdana", cex.axis=0.7, cex.lab=0.7, mar=c(4,4,2,3) - 1.5,
        mgp=c(1.1,0.25,0), tcl=0)
  }
  graph_par()
  
  if(is.null(color_puntos)==TRUE){ color_puntos <- rep('grey20',n) } else {
    
    palette <- paletteer::paletteer_c("ggthemes::Temperature Diverging", 
                                      length(unique(color)))
    color <- palette[as.numeric(color)]
    
  }
  
  if(is.null(color_ejes)==TRUE) { color_ejes <- 'firebrick'} 
  
  if(is.null(id_puntos) == TRUE){
    
    plot(scale(puntos, scale = FALSE), pch = 16, col = color_puntos)
    arrows(rep(0,p), rep(0,p), direcciones[,1], direcciones[,2],
           length = 0.1, col = color_ejes, lwd = 2)
    text(direcciones, pos = 1, cex = 0.6, col = color_ejes)
    
  } else {
    
    plot(scale(puntos, scale = FALSE), type = 'n')
    text(scale(puntos,scale=FALSE), 
         labels = id_puntos,
         col = color_puntos, cex = 0.7)
    arrows(rep(0,p), rep(0,p), direcciones[,1], direcciones[,2],
           length = 0.1, col = color_ejes, lwd = 2)
    text(direcciones, 
         labels = rownames(direcciones),
         pos = 1, cex = 0.6, col = color_ejes)
  }
  
}

## Funcion para hacer scatterplots en single factor ANOVA
mv.factorscatter <- function(formula, data, col = 'maroon',
                             width = 0.25){
  boxplot(formula, data, col = 'grey90', cex.axis = 0.8)
  beeswarm::beeswarm(formula, add = TRUE,
                     pch = 20, col = scales::alpha(col, 0.8),
                     corral = 'gutter', corralWidth = width)
}







