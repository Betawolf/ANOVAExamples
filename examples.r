
# Creates 3D matrix of dimensions a*b*r.
#
# Can take in a supplied vector `data`, which it will transform into the right 
# shape by placing the first `r` values into the first cell of the first row,
# and the second `r` values into the second cell of the first row, wrapping
# down row-by-row.
# 
# If the `data=` is left unset, the function will generate data of the required
# dimensions by using the function `rfunc`, which by default is `rnorm`. For
# example, call `3ddata(2,3,4, rfunc=runif)` to generate a 2x3 matrix with each
# cell holding 4 uniformly selected values.
data3d <- function(a,b,r, data=c(), rfunc=rnorm){
  if (length(data) < 1){
    data <- rfunc(a*b*r) 
  }
  outdata <- c()
  ir <- 1
  for (i in 1:a) {
    for (j in 1:b) {
      outdata <- c(outdata, list(data[ir:(ir+(r-1))]))
      ir <- ir + r
    }
  }
  return(matrix(outdata, nrow=a,ncol=b, byrow=T))
}

# Displays a 3D matrix as a latex table.
display3D <- function(data_matrix, alabel="Factor A", blabel="Factor B",
    alevels=c(), blevels=c()){   
  dims <- dim(data_matrix)
  a <- dims[1]
  b <- dims[2]
  
  if (length(alevels) < a){
    alevels <- 1:a
  }
  
  if (length(blevels) < b){
    blevels <- 1:b
  }

  #preamble
  cat("\\begin{table}[ht]\n")
  cat("\\centering\n")
  cat("\\begin{tabular}{|")
  cat('r', rep('r|', b+1))
  cat("}\n")
  cat("\\hline\n")
  cat(rep('&', b-1))
  cat("\\multicolumn{2}{r}{", blabel, "} &")
  cat("\\\\ \\hline \n & & ")
  cat(paste(1:b, collapse=" & "))   
  cat("\\\\ \\hline \n")
  cat("\\multirow{",a,"}{*}{",alabel,"}")
  for (i in 1:a){
    cat(" & ", i, " & ", paste(lapply(data_matrix[i,], function(x) {
    paste(round(x,2), collapse=",")}),
    collapse=" & "))
    if (i < a){
      cat("\\\\ \n \\cline{3 -",b+2,"}  \n")
    }
    else{
      cat("\\\\ \n \\hline  \n")
    } 
  }
  cat("\\end{tabular}\n\\end{table}\n")
}

workedanova <- function(data, alabel="Factor A", blabel="Factor B",
    alevels=c(), blevels=c()){
  library(xtable)

  cat("\n\\section{Example}\n")
  cat("The data are:\n\n")
  display3D(data, alabel, blabel, alevels, blevels)

  dimensions <- c(dim(data), length(unlist(data[1,1])))
  a <- dimensions[1]
  b <- dimensions[2]
  r <- dimensions[3]

  if (length(alevels) < a){
    alevels <- sapply(1:a, as.character)
  }
  
  if (length(blevels) < b){
    blevels <- sapply(1:b, as.character)
  }

  sample_mean <- mean(unlist(data))
  row_means <- apply(data, 1, function(i) {mean(unlist(i))})
  col_means <- apply(data, 2, function(j) {mean(unlist(j))})
  cel_means <- matrix(apply(data, c(1,2), function(k) { mean(unlist(k)) } ), nrow=a, ncol=b)

  means_table <- cbind(rbind(cel_means, col_means), c(row_means, sample_mean))
  means_table <- data.frame(means_table, row.names=c(alevels,
  "col"))
  names(means_table) <- c(blevels, "row")
  
  cat("\n\nThe row means ($\\bar{y}_{i..}$), column means ($\\bar{y}_{.j.}$),
  cell means ($\\bar{y}_{ij.}$) and overall mean ($\\bar{y}_{...}$)\n\n")

  print(xtable(means_table))

  SSA <- b * r * sum((row_means - sample_mean) ^ 2)
  SSB <- a * r * sum((col_means - sample_mean) ^ 2)

  MSA <- SSA/(a-1)
  MSB <- SSB/(b-1)

  SSAB <- 0
  SSerr <- 0
  SStot <- 0
  
  for (i in 1:a) {
    for (j in 1:b){
      SSAB <- SSAB + (cel_means[i,j] - row_means[i] - col_means[j] + sample_mean) ^ 2
      SSerr <- SSerr + sum(sapply(data[i,j], function(ijk) { (ijk - cel_means[i, j]) ^ 2 }))
      SStot <- SStot + sum(sapply(data[i,j], function(ijk) { (ijk - sample_mean) ^ 2 }))
    }
  }

  SSAB <- SSAB * r
  MSAB <- SSAB / ((a - 1)*(b - 1))

  MSerr <- SSerr/(length(unlist(data)) - (a*b))

  ssts <- sapply(unlist(data), function(ijk) { ijk - sample_mean })
  cat("\n\n The table of $y_{ijk} - \\bar{y}_{...}$: \n")
  display3D(data3d(a,b,r, data=ssts))

  ssts2 <- sapply(unlist(data), function(ijk) { (ijk - sample_mean) ^ 2 })
  cat("\n\n The table of $(y_{ijk} - \\bar{y}_{...}) ^ 2$: \n")
  display3D(data3d(a,b,r, data=ssts2))

  cat("\n \\[ \\mbox{SSTotal} = \\sum_{i=1}^a\\sum_{j=1}^b\\sum_{k=1}^r(y_{ijk} -
  \\bar{y}_{...})^2 = ", SStot, "\\] \n")

  cat("\n \\[ \\mbox{SSA} = br\\sum_{i=1}^a(\\bar{y}_{i..} - \\bar{y}_{...})^2
  = ", b, "\\times", r, "\\] \n")
      
#  print("Parameters: ")
#  print(paste("mu hat =",sample_mean))
#  print("alpha hats: ")
#  print(row_means - sample_mean)
#  print("beta hats: ")
#  print(col_means - sample_mean)

#  print("alpha/beta hats: ")
#  for (i in 1:a){
#    for (j in 1:b){
#      print(cel_means[i,j] - row_means[i] - col_means[j] + sample_mean)
#    }
#  }
}
