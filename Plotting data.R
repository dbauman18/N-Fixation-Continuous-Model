# we make plots of what occurs over time

data <- as.data.frame(out)


AN <- data[,2:(N^2+1)]
AP <- data[,(N^2 + 2):(2*N^2 + 1)]
LN <- data[,(2*N^2 +2):(3*N^2 + 1)]
LP <- data[,(3*N^2 + 2):(4*N^2 + 1)]
ON <- data[,(4*N^2 + 2):(5*N^2 + 1)]
OP <- data[,(5*N^2 + 2):(6*N^2 + 1)]
NonFixers <- data[,(6*N^2 + 2):(6*N^2 + 1 + j)]
Fixers <- data[,(6*N^2 + 2 + j):(6*N^2 + 1 + j + k)]

data$totalAN <- rowSums(AN)
data$totalAP <- rowSums(AP)
data$totalLN <- rowSums(LN)
data$totalLP <- rowSums(LP)
data$totalON <- rowSums(ON)
data$totalOP <- rowSums(OP)
data$nonFixerN <- rowSums(NonFixers)/wN1
data$nonFixerP <- rowSums(NonFixers)/wP1
data$fixerN <- rowSums(Fixers)/wN2
data$fixerP <- rowSums(Fixers)/wP2

Nitrogen <- c("totalAN", "totalLN","nonFixerN", "fixerN", "totalON")
data$totalN <- rowSums(data[Nitrogen])
Phosphorus <- c("totalAP", "totalLP", "nonFixerP", "fixerP", "totalOP")
data$totalP <- rowSums(data[Phosphorus])



# plot all of the pools over time (totalled over the entire region)
layout_matrix_pools <- matrix(1:8, nr = 4, ncol = 2)
layout(layout_matrix_pools)
plot(times, data$totalN, 'l', main = "Total Nitrogen",xlab = "Time [yr]", ylab = "g N")
plot(times, data$totalAN, 'l', main = "Available Nitrogen",xlab = "Time [yr]", ylab = "g N")
plot(times, data$totalON, 'l', main = "Organic Nitrogen",xlab = "Time [yr]", ylab = "g N")
plot(times, data$totalLN, 'l', main = "Litter Nitrogen",xlab = "Time [yr]", ylab = "g N")
plot(times, data$totalP, 'l', main = "Total Phosphorus",xlab = "Time [yr]", ylab = "g P")
plot(times, data$totalAN, 'l', main = "Available Phosphorus",xlab = "Time [yr]", ylab = "g P")
plot(times, data$totalOP, 'l', main = "Organic Phosphorus",xlab = "Time [yr]", ylab = "g P")
plot(times, data$totalLP, 'l', main = "Litter Phosphorus",xlab = "Time [yr]", ylab = "g P")



layout_matrix_biomass <- matrix(1:2, nr = 1, ncol = 2)
layout(layout_matrix_biomass)
plot(times, data$nonFixerN*wN1, 'l', main = "Non-Fixer Biomass",xlab = "Time [yr]", ylab = "g C")
plot(times, data$fixerN*wP1, 'l', main = "Fixer Biomass",xlab = "Time [yr]",ylab = "g C")

