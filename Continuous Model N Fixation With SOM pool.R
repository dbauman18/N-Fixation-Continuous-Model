# We create the grid and set the locations of our plants
library(ReacTran)


#_____________________________________________________________________________#
# This sets up the grid, the arrangements of trees, and computes the 
# rooting kernels/litterfall kernels
L <- 20; # length of square grid in meters
areaM2 <- L^2 # area of grid in m^2
N <- 20; # number of grid cells per side
cellArea <- areaM2 / N^2 # area of single cell in m^2

x.grid <- setup.grid.1D(0, L, N = N)
y.grid <- setup.grid.1D(0, L, N = N)

Grid <- setup.grid.2D(x.grid = x.grid, y.grid = y.grid)


# set tree locations randomly

j <- 5 #number NF
k <- 5 #number F
r <- 3

xNF <- runif(n = j, min = 0, max = L)
yNF <- runif(n = j, min = 0, max = L)
xF <- runif(n=k, min = 0, max = L)
yF <- runif(n=k, min = 0, max = L)
rMaxNF <- runif(n = j, min = r/2, max = 3*r/2)
rMaxF <- runif(n = k, min = r/2, max = 3*r/2)

x <- Grid$x.mid
y <- Grid$y.mid



# store the rooting kernel and litter arrays for each plant
lambda <- 3# rate of litter transfer. Lower lambda means more dispersal

rootingKernelNF <- list()
litterKernelNF <- list()
for (i in 1:j){
  root <- outer(x,y, function(x,y){
    v <- rMaxNF[i]^2 - ((x-xNF[i] + 3*L/2) %% L - L/2)^2 - ((y-yNF[i]+3*L/2) %% L - L/2)^2
    return(v * (v>= 0))
  })
  rootingKernelNF[[i]] <- root/sum(root)*rMaxNF[i]^2 # rMax is the benefit multiplier of larger roots
  lit <- outer(x,y, function(x,y){
    return(lambda^2*exp(-lambda*sqrt(((x-xNF[i] + 3*L/2) %% L - L/2)^2 + ((y-yNF[i]+3*L/2) %% L - L/2)^2)/(2*pi)))
  })
  litterKernelNF[[i]] <- lit/sum(lit)
}

rootingKernelF <- list()
litterKernelF <- list()
for (i in 1:k){
  root <- outer(x,y, function(x,y){
    v <- rMaxF[i]^2 - ((x-xF[i] + 3*L/2) %% L - L/2)^2 - ((y-yF[i]+3*L/2) %% L - L/2)^2
    return(v * (v>= 0))
  })
  rootingKernelF[[i]] <- root/sum(root)*rMaxF[i]^2
  lit <- outer(x,y, function(x,y){
    return(lambda^2*exp(-lambda*sqrt(((x-xF[i] + 3*L/2) %% L - L/2)^2 + ((y-yF[i]+3*L/2) %% L - L/2)^2)/(2*pi)))
  })
  litterKernelF[[i]] <- lit/sum(lit)
}




#______________________________________________________________________________#
# next, we will define our functions for growth and compute the PDE
# growth of non-fixer as function of available N and P (AN/AP)
gB1 <- function(AN,AP){
  g <- min(wN1*vN1*AN, wP1*vP1*AP)
  return(g)
}


# growth of fixer as a function of available N and P (AN/AP)
gB2 <- function(AN, AP){
  f <- max(0, min(fMax, wP2*vP2*AP/wN2 - vN2*AN))
  return(min(wN2*(vN2*AN + f), wP2*vP2*AP))
}

# growth of an obligate fixer at rate O with biomass B 
gB2Obligate <- function(AN,AP,B) {
  return(min(wN2*(vN2*AN + O*B), wP2*vP2*AP))
}



pdeSpatialModel <- function(t, y, parms) {
  AN = y[1:N^2]
  AP = y[(N^2 + 1): (2*N^2)]
  LN = y[(2*N^2 + 1): (3*N^2)]
  LP = y[(3*N^2 + 1): (4*N^2)]
  ON = y[(4*N^2 + 1): (5*N^2)]
  OP = y[(5*N^2 + 1): (6*N^2)]
  NonFixers = y[(6*N^2 + 1): (6*N^2 + j)]
  Fixers = y[(6*N^2 + j+1) : (6*N^2 + j+k)]
  
  # now we compute advection/diffusion for available nitrogen
  availableNMatrix <- matrix(nr = N, nc = N, AN)
  top <- availableNMatrix[1,]
  left <- availableNMatrix[,1]
  bottom <- availableNMatrix[N,]
  right <- availableNMatrix[,N]
  tranAN <- tran.2D(C = availableNMatrix, A.x = 1, A.y = 1, D.x = DA, D.y = DA, 
                    grid = Grid, v.x = uxA, v.y = uyA, C.x.down = right, 
                    C.x.up = left, C.y.down = top, C.y.up = bottom)$dC
  
  # do same for available phosphorus
  availablePMatrix <- matrix(nr = N, nc = N, AP)
  top <- availablePMatrix[1,]
  left <- availablePMatrix[,1]
  tranAP <- tran.2D(C = availablePMatrix, A.x = 1, A.y = 1, D.x = DA, D.y = DA, 
                    grid = Grid, v.x = uxA, v.y = uyA, C.x.down = left, 
                    C.x.up = left, C.y.down = top, C.y.up = top)$dC
  
  #Do same for litter nitrogen and phosphorus
  litterNMatrix <- matrix(nr = N, nc = N, LN)
  top <- litterNMatrix[1,]
  left <- litterNMatrix[,1]
  tranLN <- tran.2D(C = litterNMatrix, A.x = 1, A.y = 1, D.x = DL, D.y = DL, 
                    grid = Grid, v.x = uxL, v.y = uyL, C.x.down = left, 
                    C.x.up = left, C.y.down = top, C.y.up = top)$dC
  
  litterPMatrix <- matrix(nr = N, nc = N, LP)
  top <- litterPMatrix[1,]
  left <- litterPMatrix[,1]
  tranLP <- tran.2D(C = litterPMatrix, A.x = 1, A.y = 1, D.x = DL, D.y = DL, 
                    grid = Grid, v.x = uxL, v.y = uyL, C.x.down = left, 
                    C.x.up = left, C.y.down = top, C.y.up = top)$dC
  
  # Next, we can compute the impact of each of the tree at the same time as we
  # compute the growth function for each
  x <- Grid$x.mid
  y <- Grid$y.mid
  
  if (seasonal_litterfall){
    u1 <- 0.05 + 1.5*(t%%1 > 0.6)*(t%%1 < 0.7)
    u2 <- 0.05 + 1.5*(t%%1 > 0.6)*(t%%1 < 0.7)
  }
  
  # Last, we need to compute the amount of litter nutrients added, taking
  # account of the distribution the place falls
  # ADD THIS LATER!!!!! Right now we are assuming litter does not get added back
  for (m in 1:j)
  {
    drop <- u1*NonFixers[m]*litterKernelNF[[m]]
    # add litter to pool
    tranLN <- tranLN + drop/wN1
    tranLP <- tranLP + drop/wP1
  }
  
  # Do it for the fixers
  for (n in 1:k)
  {
    drop <- u2*Fixers[n]*litterKernelF[[n]]
    # add litter to pool
    tranLN <- tranLN + drop/wN2
    tranLP <- tranLP + drop/wP2
  }
  
  # We do rooting now
  # go through non-fixers
  for (m in 1:j)
  {
    roots <- rootingKernelNF[[m]]
    nUptake <- roots*availableNMatrix
    pUptake <- roots*availablePMatrix
    
    ANm <- sum(nUptake)
    APm <- sum(pUptake)
    
    # growth of plant
    g <- gB1(ANm,APm)
    # loss of available nutrients
    tranAN <- tranAN - NonFixers[m]*g*nUptake/sum(nUptake)/wN1
    tranAP <- tranAP - NonFixers[m]*g*pUptake/sum(pUptake)/wP1
    NonFixers[m] <- NonFixers[m]*(g - u1)
  }
  
  # go through fixers
  for (n in 1:k)
  {
    roots <- rootingKernelF[[n]]
    
    nUptake <- roots*availableNMatrix
    pUptake <- roots*availablePMatrix
    ANm <- sum(nUptake)
    APm <- sum(pUptake)
    # growth of plant
    g <- gB2(ANm, APm)
    # loss of available nutrients
    tranAN <- tranAN - Fixers[n]*g*nUptake/sum(nUptake)/wN2
    tranAP <- tranAP - Fixers[n]*g*pUptake/sum(nUptake)/wP2
    Fixers[n] <- Fixers[n]*(g - u2)
  }
  
  dAN <-  c(tranAN) + IN + gammaN*kN*ON -mN*AN
  dAP <- c(tranAP) + IP + gammaP*kP*OP - mP*AP
  dLN <- c(tranLN) - dN*LN
  dLP <- c(tranLP) - dP*LP
  dON <- eN*dN*LN - kN*ON
  dOP <- eP*dP*LP - kP*OP
  
  
  # finally, we add all of the derivatives together and turn into a vector
  return(list(c(dAN, dAP, dLN, dLP, dON, dOP, NonFixers, Fixers)))
}



#_____________________________________________________________________________#
# Now, we are ready to set all the parameters and run the model
wN1 <- 50 # g C/g N
wN2 <- 38.5 # g C/g N
vN1 <- 0.0002 # m^2/g C/y
vN2 <- 0.0002 # m^2/g C/y



fMax <- 0.05 # g N/g C/y
wP1 <- 600 # g C/g P
wP2 <- 600 # g C/g P
vP1 <- 0.0002 # m^2/g C/y
vP2 <- 0.0002 # m^2/g C/y
u1 <- 0.1 # /y
u2 <- 0.1 # /y

DA <- 4 # diffusion constant of available
DL <- 4 # diffusion constant of litter
uxA <- -14 # x advection of available
uyA <- 0 # y advection of available
uxL <- -10 # x advection of litter
uyL <- 0 # y advection of litter


IN <- 0.6 # g N/m^2/y input of N 
IP <- 0.06 # g P/m^2/y abiotic input of P

# rate at which SOM becomes available, and amount which actually becomes available (rest mineralized)
kN <- 0.1
kP <- 0.1
gammaN <- 0.98
gammaP <- 0.98


dN <- 4 # rate of litter decomposition
dP <- 4 # rate of litter decomposition
eN <- 0.9 # proportion of litter N decomp that goes to organic
eP <- 0.9 # proportion of litter P decomp that goes to organic


mN <- 0.04 # rate which available N lost to environment
mP <- 0.03 # rate which available P lost to environment


#_____________________________________________________________________________#
# we need to set our initial conditions and run the code
# soil organic matter pool
SOM_N <- 163 # g N/m^2
SOM_P <- 15.8 # g P/m^2
proportionAvailable <- 0.03 # 3% of organic pool is available for initial conditions
iaN <- SOM_N*proportionAvailable * cellArea
iaP <- SOM_P*proportionAvailable * cellArea
ilN <- 25*cellArea
ilP <- 2*cellArea
iBNF <- 5000
iBF <-5000

seasonal_litterfall <- FALSE #If true, then u1, u2 become periodic (most litter drops at a certain time)
ilN <- 0
ilP <- 0


yInitial = c(rep(iaN, N^2), rep(iaP, N^2), rep(ilN, N^2), rep(ilP, N^2), 
             rep(SOM_N, N^2), rep(SOM_P, N^2), rMaxNF*iBNF/r, rMaxF*iBF/r)
require(deSolve)

t <- 1000
times <- seq (from = 0, to = t, by = 0.2) #t year simulation


print(system.time(out <- lsodes(y = yInitial, times = times, func = pdeSpatialModel, parms = NULL)))

sum(out < 0) # if this is positive, the simulation failed (negative values)

