###### Forming clusters using data frame  ###
library(factoextra)
df <- USArrests
df <- na.omit(df)
df <- scale(df)
d <- dist(df, method = "euclidean")
hc <- hclust(d)
sub_grp <- cutree(hc, k = 4)         ### forming 4 subgroups
fviz_cluster(list(data = df, cluster = sub_grp))


#### solving & plotting a 3d eqnN #####
library(rgl)
dat <- replicate(2, 1:3)        # Create some dummy data
plot3d(dat, type = 'n', xlim = c(-1, 8), ylim = c(-1, 8), zlim = c(-10, 20), xlab = '', ylab = '', zlab = '')  # Initialize the scene, no data plotted
planes3d(2, 3, -1, 0, col = 'red', alpha = 0.6) # Solving 2*x + 3*y - 1*z = 0
points3d(x=0, y=0, z=0) # Define the origin


#### finding row space, column space and nullspace of a matrix
library(pracma)
library(matlib)
A <- matrix(c(1, -1, 4, 2, 0, -1, -1, -1, 5), nrow=3, ncol=3, byrow=TRUE)
rref(A)            # rowspace of a matrix
orth(t(A))         # rowspace of a matrix
t(rref(t(A)))      # columnspace of a matrix
orth(A)            # columnspace of a matrix
Null(t(A))         # Nullspace of a matrix
nullspace(A)       # Nullspace of a matrix


### Principal Component Analysis(PCA)
library(mvtnorm)
set.seed(123)
sigma <- matrix(c(4,2,2,3), ncol=2)
x <- rmvnorm(n=500, mean=c(0,0), sigma=sigma)
pca <- princomp(x)
pca$loadings
Tr <- pca$loadings
D <- inv(Tr) %*% t(x)
y <- t(D)
plot(y)


### Using Linear Programming
library(lpSolve)
A <- matrix(c(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 1), nrow = 11, ncol = 18, byrow = TRUE)

f.con <- rbind(A, diag(18))
f.rhs <- c(1298, 1948, 465, 605, 451, 338, 260, 183, 282, 127, 535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

f.obj <- c(39.99, 126.27, 102.70, 81.68, 38.81, 71.99, 31.21, 22.28, 321.04, 145.36, 33.82, 154.05, 64.19, 87.90, 107.98, 65.45, 39.08, 167.38)
f.dir <- c("=", "=", "=", "=", "=", "=", "=", "=", "=", "=", "=", 
           ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", 
           ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=")
Sol <- lp ("min", f.obj, f.con, f.dir, f.rhs)
Sol$objval               ## optimal solution of the condition
Sol$solution             ## values for which the soln is optimal


#### Network Analysis
library(pracma)
library(igraph)
library(igraphdata)

G <- sample_gnp(10, 0.5)       #### undirected graph
plot(G)
as_adjacency_matrix(G)       ### always a symmetric matrix
N <- barabasi.game(10)        ### directed graph
plot(N)
as_adjacency_matrix(N)

D = diag(c(3,3,3,2,1))
L <- matrix(c(3, -1, -1, -1, 0, -1, 3, -1, 0, -1, -1, -1, 3, -1, 0, -1, 0, -1, 2, 0, 0, -1, 0, 0, 1), 5, 5)
dim(nullspace(L))

# normalized laplacian
LL <- L         
LL[1,] <- LL[1,]/LL[1,1]
LL[2,] <- LL[2,]/LL[2,2]
LL[3,] <- LL[3,]/LL[3,3]
LL[4,] <- LL[4,]/LL[4,4]
LL[5,] <- LL[5,]/LL[5,5]     ### LL is normalized laplacian of L
E <- eigen(LL)
E$values
sum(E$values)

data("macaque")

LM <- graph.laplacian(macaque)
LM <- as.numeric(LM)
LLM <- matrix(LM, 45, 45)
dim(nullspace(LLM))


#### bivariate data
PID<-read.table("http://stat4ds.rwth-aachen.de/data/PartyID.dat", header=TRUE)
table(PID$race,PID$id) #forms contingency table (not shown here;see Table 1.1)
options(digits=2)
prop.table(table(PID$race,PID$id),margin=1) #Formargin=1, proportions
mosaicplot(table(PID$race,PID$id)) #graphical portrayalof cell sizes


#### binomial distribution
y=seq(0, 12, 1) # y values between 0 and 12 with increment of 1
plot(y,dbinom(y,12,0.2),type="h") # plots binomial probabilities when n=12,pi=0.2
mu<-function(n,pi){n*pi} # function: binomial mean
sigma<-function(n,pi){sqrt(n*pi*1-pi)} # function: binomialstandard deviation
Psd<-function(n,pi,k){pbinom(mu(n,pi)+k*sigma(n,pi),n, pi) - pbinom(mu(n,pi)-k*sigma(n,pi),n,pi)} # function:prob.within k std.dev.ofmean


n=1500;pi=0.60
Psd(n,pi,2) # probability within k=2 standard deviations of mean


probNormal <- function(a,b,mu,sigma){
  prob <- pnorm(b,mu,sigma)-pnorm(a,mu,sigma)
  low <- min(mu-4*sigma,a); up <- max(mu+4*sigma,b)
  curve(dnorm(x,mu,sigma), xlim=c(low,up), main=" ", xlab="y",
        ylab=expression(phi(y)), col="blue", lwd=2)
  x=seq(a,b,length=200)
  y=dnorm(x,mu,sigma)
  polygon(c(a,x,b),c(0,y,0),col="coral2")
  result <- paste("N(",mu,",",sigma^2,"): ",
                  "P(",a,"< Y <",b,") =",signif(prob, digits=3))
  mtext(result,3, col="coral2")
  # text(low,0.15,paste(prob),cex=1,col="grey60")
  curve(dnorm(x,mu,sigma),col="blue", lwd=2,add=T)
  return(prob)}
#-----------------------------------------------------------------------
probNormal(2,5,3,3) # example

#####################################


#### multinomial distribution
dmultinom(c(0,1,11),prob=c(0.20, 0.30, 0.50))



#### Joint probability
GSS <- read.table("http://stat4ds.rwth-aachen.de/data/GSS2018.dat", header=T)
gender <- factor(GSS$SEX, levels=c(1,2), labels = c("Male","Female"))
smallgap <- factor(GSS$SMALLGAP, levels=c(1:5), labels = c("strongly agree", "agree","neutral","disagree","strongly disagree"))
fairsociety <- table(gender,smallgap) # frequency table
fairsociety
joint.prob <- prop.table(fairsociety) # derives proportion table, # from a frequency table
cond.prob1 <- prop.table(fairsociety, 1) # cond. prop. within rows
cond.prob2 <- prop.table(fairsociety, 2) # cond. prop. within columns

barplot(cond.prob2, density=40, main="In a fair society, differences
in people's standard of living should be small",xlab="", ylab="conditional proportions",ylim=c(0,1))
abline(h=0.5, col="blue")


### WEIGHTED CORRELATION'
library(wCorr)
probabilities <- c(0.2,0.1,0.0,0.1,0.2,0.1,0.0,0.1,0.2)
x <- c(1,1,1,2,2,2,3,3,3)
y <- c(1,2,3,1,2,3,1,2,3) #scores
weightedCorr(x, y, weights=probabilities, method="polyserial")


### BINOMIAL CLT
CLT_binom <- function(B,n,p) {
  # B: number of iterations used to approximate the distribution of Xmean
  # n: sample size
  # p: success probability pi
  Y <- rbinom(B,n,p)
  Ymean <- Y/n # vector (length B) with the p-estimates: Algorithm 1 (2)
  var.mean <-p*(1-p)/n # variance of the estimator of p
  p.MC <- mean(Ymean) # Monte Carlo estimate of p
  varp.MC <- var(Ymean) # MC variance estimate of var.mean
  h <- hist(Ymean, col = "gray", probability=TRUE, main=paste("n=",n))
  xfit<-seq(0, 1,length=5000)
  yfit<-dnorm(xfit,mean=p,sd=sqrt(p*(1-p)/n))
  gr <- lines(xfit, yfit, col="blue",lwd=2)
  list(var.mean=var.mean, p.MC=p.MC, varp.MC=varp.MC) }

par(mfrow=c(2,2)) # multiple graphs layout in a 2x2 table format
CLT_binom(100000, 10, 0.3)


### POISSON CLT
pois_CLT <- function(n, mu, B) {
  # n: vector of 2 sample sizes [e.g. n <- c(10, 100)]
  # mu: mean parameter of Poisson distribution
  # B: number of simulated random samples from the Poisson
  par(mfrow = c(2, 2))
  for (i in 1:2){
    Y <- numeric(length=n[i]*B)
    Y <- matrix(rpois(n[i]*B, mu), ncol=n[i])
    Ymean <- apply(Y, 1, mean) # or, can do this with rowMeans(Y)
    barplot(table(Y[1,]), main=paste("n=", n[i]), xlab="y",
            col="lightsteelblue") # sample data dist. for first sample
    hist(Ymean, main=paste("n=",n[i]), xlab=expression(bar(y)),
         col="lightsteelblue") # histogram of B sample mean values
  } }

    # implement:with 100000 random sample sizes of 10 and 100, mean = 0.7
n <- c(10, 100)
pois_CLT(n, 0.7, 100000)
 

### cOMPARING MEDIANS
Y<-matrix(rpois(10*100000,0.7), ncol=10) # simulate Poisson rv's with mean 0.7
Ymed<-apply(Y,1, median) # find median for each sample of size 10
hist(Ymed,freq=FALSE) # histogram of the 100,000 sample medians

sdmed <- function(B, n, mu, sigma){
  medians <- rep(0, B)
  for(i in 1:B){
    y <- rnorm(n, mu, sigma)
    medians[i] <- median(y) }
  sd(medians)
}

sdmed(1000000, 100, 100, 16) # B=1000000, n=100, mu=100, sigma=16



