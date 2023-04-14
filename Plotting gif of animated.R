library(tidyverse)
library(manipulate)
library(gganimate)
library(gifski)

# # j = 1 -> AN, j = 2 -> AP, j = 3 -> LN, j = 4 -> LP
# make_plot <- function(t, j) {
#   
#   aN <- matrix(nr = N, nc = N, out[t, 2 + N^2*(j-1):(N^2*j+1)])
#   aN <- cbind(x,as.data.frame(aN))
#   colnames(aN) <- c("x", y)
#   
#   return(pivot_longer(aN, cols=-x, names_to="y", values_to="value") %>%
#            mutate(y=as.numeric(y)) %>%
#            ggplot() +
#            geom_tile(aes(x=x,y=y,fill=value)) +
#            scale_fill_gradientn(limits=c(min(out[, 2:(N^2+1)]),max(out[, 2:(N^2+1)])),
#                                 colors=c("blue","lightblue")))
#   
# }
# 
# # require(manipulate)
# # manipulate(plot(1:x), x = slider(5,10))
# # manipulate(make_plot(t), t = slider(1,200,step=1))

x <- Grid$x.mid
y <- Grid$y.mid

out[,1:(N^2+1)] %>%
  as.data.frame %>% 
  pivot_longer(cols=-time, names_to="row", values_to="value") %>%
  mutate(
    row = as.numeric(row),
    x = x[(row-1) %% N + 1],
    y = y[(row-1) %/% N +1]) %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=value)) +
  scale_fill_gradientn(limits=c(min(out[, 2:(N^2+1)]),max(out[, 2:(N^2+1)])),
                       colors=c("blue","green")) +
  transition_time(time) + 
  labs(title = "Available nitrogen pool over 20 years")



