library("jpeg")
image = readJPEG("forest_medium_sharp.jpg")
head(image)
par(mfcol=c(1,1), mar=c(0,0,2,0))
plot.new()
rasterImage(image, 0, 0, 1, 1, interpolate = FALSE)

# step 0 - initial colorspace clustering
# use simple k-means clustering
#transform image into rgb vectors
data = data.frame(list(r = c(image[,,1]), g = c(image[,,2]), b = c(image[,,3])))

head(data)
#lets try different k values
for(k in 2:16) {
   #choose initial centers
  mu <- runif(k * 3, 0, 1)
  dim(mu) <- c(k, 3)
  
  # lets do clustering
  result = kmeans(data, k)
  # restore image
  rimage = array(NA, dim=c(dim(image)[1], dim(image)[2], 3))
  dim(rimage)
  rimage[,,1] = result$centers[result$cluster, 1]
  rimage[,,2] = result$centers[result$cluster, 2]
  rimage[,,3] = result$centers[result$cluster, 3]
  plot.new()
  rasterImage(rimage, 0, 0, 1, 1, interpolate = FALSE)  
  mtext(paste("Number of clusters:", k))
}

#choose k=6 and recalculate cluster centers
k = 6
result = kmeans(data, k)
  
rimage = array(NA, dim=c(dim(image)[1], dim(image)[2], 3))
dim(rimage)
rimage[,,1] = result$centers[result$cluster, 1]
rimage[,,2] = result$centers[result$cluster, 2]
rimage[,,3] = result$centers[result$cluster, 3]
plot.new()
rasterImage(rimage, 0, 0, 1, 1, interpolate = FALSE)  
mtext(paste("Number of clusters:", k))



counts <- rep(1, k)
barplot(counts, main="Cluster Centers", 
        col=rgb(result$centers[,1],result$centers[,2],result$centers[,3]))

## 1 split image into 16x16 blocks

blocks.size = 16
blocks.dim = dim(image)[1:2] / blocks.size

blocks = matrix(NA, ncol=3*blocks.size**2, nrow=blocks.dim[1]*blocks.dim[2])
for (i in 1:(blocks.dim[1])) {
  for (j in 1:(blocks.dim[2])) {
    RGB = image[((i-1)*blocks.size+1):(i*blocks.size), ((j-1)*blocks.size+1):(j*blocks.size), ]
    blocks[(i-1)*blocks.dim[2]+j, ] = RGB
  }  
}

# make sure we can restore image
blocks.plot <- function (blocks, title) {
  rimage <- array(0, dim=c(blocks.dim[1]*blocks.size, blocks.dim[2]*blocks.size,3))
  dim(rimage)
  for(i in 1:blocks.dim[1]){
    for(j in 1:blocks.dim[2]){
      rimage[(blocks.size*(i-1) + 1:blocks.size),(blocks.size*(j-1) + 1:blocks.size), ] <- 
        blocks[(i-1)*blocks.dim[2] + j, ]
    }
  }
  plot.new()
  rasterImage(rimage, 0, 0, 1, 1, interpolate = FALSE)  
  mtext(paste(title))
}

blocks.plot(blocks, "restored image")

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# assign blocks into initial clusters that we computed using k-means
z <- rep(-1, blocks.dim[1] * blocks.dim[2])

# cluster indexes for each pixel
zimage = result$cluster
dim(zimage) <- c(dim(image)[1], nrow= dim(image)[2])

# find dominant pixel for each block
for(i in 1:blocks.dim[1]){
  for(j in 1:blocks.dim[2]){
    pixels <- zimage[(blocks.size*(i-1) + 1:blocks.size),(blocks.size*(j-1) + 1:blocks.size)]
    mode = pixels[which.max(tabulate(match(pixels,1:4)))] # mode function
    z[(i-1) * blocks.dim[2] + j] <- mode
  }
}



errors <- c()

for (m in 1:20) {
  
  #plot inital clustering of blocks
  bimage <- array(0, dim=c(blocks.dim[1]*blocks.size, blocks.dim[2]*blocks.size,3))
  dim(rimage)
  for(i in 1:blocks.dim[1]){
    for(j in 1:blocks.dim[2]){
      bimage[(blocks.size*(i-1) + 1:blocks.size),(blocks.size*(j-1) + 1:blocks.size),1] <- 
        rep(result$centers[z[(i-1) * blocks.dim[2] + j],1], blocks.size**2)
      bimage[(blocks.size*(i-1) + 1:blocks.size),(blocks.size*(j-1) + 1:blocks.size),2] <- 
        rep(result$centers[z[(i-1) * blocks.dim[2] + j],2], blocks.size**2)
      bimage[(blocks.size*(i-1) + 1:blocks.size),(blocks.size*(j-1) + 1:blocks.size),2] <- 
        rep(result$centers[z[(i-1) * blocks.dim[2] + j],3], blocks.size**2)
    }
  }
  plot.new()
  rasterImage(bimage, 0, 0, 1, 1)  
  mtext(paste("Initial clustering"))

  
  pca.vectors = array(NaN, dim=c(k, ncol(blocks), 16))

  old.par <- par(mfrow=c(4,4), mar = c(0, 0, 0, 0), oma = c(0, 0, 2, 0))
  # 
  for(i in seq_len(k))
  {
    # covariance matrix
    Sigma <- t(blocks[z == i, ]) %*% blocks[z == i, ] # cov(blocks[z == i,]) 
    # Use eigenvector decomposition to find the directions
    W <- eigen(Sigma)
    pca.vectors[i,,]  = W$vectors[,1:16]

    for(j in 1:16){
      block = W$values[j] %*% t(W$vectors[,j])
      block = range01(block)
      dim(block) <- c(16,16,3)
      plot.new()
      rasterImage(block, 0, 0, 1, 1)  
    }  
    mtext(paste("cluster ",i," components, iteration =", m), outer=TRUE)
  }
  par(old.par)
  
  # restore each block by all cluster components
  blocks.masks = array(NaN, dim=c(k, nrow(blocks), ncol(blocks)))
  
  for(i in seq_len(k))
  {
    for (a in 1:nrow(blocks)) {
      weights = t(pca.vectors[i,,]) %*% blocks[a,]
      mask = pca.vectors[i,,] %*% weights
      blocks.masks[i,a, ] = mask
    }
    blocks.plot(range01(blocks.masks[i,, ]), paste("Restored image from cluster ",i,", iteration =", m))
    blocks.plot(range01(abs(blocks - blocks.masks[i,, ])), paste("Errors from cluster ",i,", iteration =", m))
  }

  
  # calculate errors
  blocks.errors <- matrix(NA, nrow(blocks), k)
  for(i in seq_len(k))
  {
    for (a in 1:nrow(blocks)) {
      blocks.errors[a, i] <- mean((blocks[a,] - blocks.masks[i,a, ])**2)
    }
  }
  
  # select cluster with smaller error
  z.new <- apply(blocks.errors, 1, which.min)

  # restore image from best cluster
  composite = matrix(NA, ncol=3*blocks.size**2, nrow=blocks.dim[1]*blocks.dim[2])
  dim(composite)
  for (a in 1:nrow(blocks)) {
    composite[a,] = blocks.masks[z.new[a],a, ]
  }
  blocks.plot(range01(composite), paste("Restored image, iteration =", m))   
  blocks.plot(range01(abs(blocks - composite)), paste("Restored image error, iteration =", m))

  #overall error
  errors <- c(errors, mean((blocks - composite)**2))
  plot(errors, type = "l", lwd = 1, xlab = "", ylab = "", main = "MSE")
  points(errors, pch = 16)
  

  bimage <- array(0, dim=c(blocks.dim[1]*blocks.size, blocks.dim[2]*blocks.size,3))
  dim(rimage)
  for(i in 1:blocks.dim[1]){
    for(j in 1:blocks.dim[2]){
      bimage[(blocks.size*(i-1) + 1:blocks.size),(blocks.size*(j-1) + 1:blocks.size),1] <- 
        rep(result$centers[z[(i-1) * blocks.dim[2] + j],1], blocks.size**2)
      bimage[(blocks.size*(i-1) + 1:blocks.size),(blocks.size*(j-1) + 1:blocks.size),2] <- 
        rep(result$centers[z[(i-1) * blocks.dim[2] + j],2], blocks.size**2)
      bimage[(blocks.size*(i-1) + 1:blocks.size),(blocks.size*(j-1) + 1:blocks.size),2] <- 
        rep(result$centers[z[(i-1) * blocks.dim[2] + j],3], blocks.size**2)
    }
  }
  plot.new()
  rasterImage(bimage, 0, 0, 1, 1)  
  mtext(paste("Clustered blocks, iteration =", m))
  
  if (all(z==z.new)){
    break
  }
  # repeat
  z = z.new
  
}