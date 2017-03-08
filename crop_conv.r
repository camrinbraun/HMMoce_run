load('C:/Users/benjamin.galuardi/Downloads/mako_example.RData')

library(raster)
library(OpenImageR)

inv2 <- r
#inv2 = inv

#ex2 = extent(-1, -90, 45, 55)
# using ex from other script
ex2 <- ex
ridx = cellsFromExtent(inv2, ex2)

c1 = as.matrix(crop(inv2, ex2))

k1 = gausskern(5,2)

conv1 = convolution(c1, k1)

inv2 = (inv2*0)+1e-15
inv2[ridx] = t(conv1)

pdf('try convolve.pdf',width=12,height=8)
par(mfrow=c(1,2))
plot(r)
plot(ex2, add=T, col = 2)
plot(inv2)
dev.off()
