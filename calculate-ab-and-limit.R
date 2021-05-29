

# mutation count at different local GC% quantiles  (local = 2kb)
T2a <- read.table("table2a-2kb-window-intergenic-control.input", header=T)

# 2kb gc% calculated from the human genome
T2b <- read.table("table2b-2kb-intergenic.input", header=T)

######### subroutines -------------------------------------------
################# function: calculate alpha's and beta's
ab <- function(count6){
 nws <- count6[,1]+count6[,2]
 nsw <- count6[,3]+count6[,5]
 a1 <- nsw/(nws+nsw)
 a2 <- nws/(nws+nsw)
 b1 <- count6[,5]*count6[,3] + count6[,4]*count6[,5]+count6[,6]*count6[,3]
 b2 <- count6[,1]*count6[,6] + count6[,2]*count6[,6]+count6[,5]*count6[,1]
 b3 <- count6[,2]*count6[,4] + count6[,1]*count6[,4]+count6[,3]*count6[,2]
 bsum <- b1+b2 +b3
 b1 <- b1/bsum
 b2 <- b2/bsum
 b3 <- b3/bsum
 tmp <- cbind(a1, a2, b1,b2, b3)
 tmp
}

################# function: return limiting x 
xlimit2 <- function(A, x){
 tmp <- 1/(A[,1]/A[,2]* ((1-x)/x) +1)
 tmp
}

################# function: return limiting x, limiting y
xylimit3 <- function(B, x, y){
 w <- B[,1]*(1-x)
 s <- B[,2]*x*(1-y)
 s2 <- B[,3]*x*y
 sum3 <- w+s +s2
 w <- w/sum3
 s <- s/sum3
 s2 <- s2/sum3
 cg <- s+s2
 cpg <- s2/(s+s2)
 tmp <- cbind(cg, cpg)
 tmp
} 
######### end of subroutines -------------------------------------------


ab.intergenic.control.2kb <- ab(T2a[, 4:9])

x.intergenic.2kb <- (T2b[,3]+T2b[,4])/2
y.intergenic.2kb <- (T2b[,5]+T2b[,6]+T2b[,7])/3

xlimit2.intergenic.control.2kb <- xlimit2(ab.intergenic.control.2kb[, 1:2], x.intergenic.2kb)
xlimit3.intergenic.control.2kb <- xylimit3(ab.intergenic.control.2kb[, 3:5], x.intergenic.2kb,y.intergenic.2kb)


## show the result
## last 5 columns of Table 2(i)

cbind(T2a[, 1:2],  ab.intergenic.control.2kb )


#	   gc0  gc1        a1        a2        b1        b2        b3
#	1 0.00 0.35 0.5206176 0.4793824 0.3548048 0.3282431 0.3169521
#	2 0.35 0.38 0.5478261 0.4521739 0.3871115 0.3287313 0.2841571
#	3 0.38 0.40 0.5759504 0.4240496 0.4047506 0.2982819 0.2969675
#	4 0.40 0.43 0.6120916 0.3879084 0.4512356 0.2982248 0.2505396
#	5 0.43 0.47 0.6410438 0.3589562 0.4747930 0.2693332 0.2558739
#	6 0.47 1.00 0.6415902 0.3584098 0.4729347 0.2648238 0.2622415
#	7 0.47 0.50 0.6384892 0.3615108 0.4696151 0.2668568 0.2635281
#	8 0.50 0.54 0.6295249 0.3704751 0.4595317 0.2705762 0.2698921
#	9 0.54 1.00 0.6566265 0.3433735 0.4900213 0.2571874 0.2527913
	

	
## last 3 columns of Table 2(ii)

cbind(T2a[, 1:2], xlimit2.intergenic.control.2kb , xlimit3.intergenic.control.2kb)

#	   gc0  gc1 xlimit2.intergenic.control.2kb        cg        cpg
#	1 0.00 0.35                      0.3100764 0.3108963 0.02454773
#	2 0.35 0.38                      0.3222476 0.3275790 0.02627048
#	3 0.38 0.40                      0.3191483 0.3193197 0.03395483
#	4 0.40 0.43                      0.3092631 0.3169395 0.03308305
#	5 0.43 0.47                      0.3120239 0.3143185 0.04377013
#	6 0.47 1.00                      0.3620974 0.3625020 0.06387820
#	7 0.47 0.50                      0.3455085 0.3461685 0.05497753
#	8 0.50 0.54                      0.3860006 0.3860861 0.06464676
#	9 0.54 1.00                      0.4113758 0.4118747 0.09165449


