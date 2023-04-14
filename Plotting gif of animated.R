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

times <- out[,1]
AN <- out[,2:(N^2+1)]
AP <- out[,(N^2 + 2):(2*N^2 + 1)]
LN <- out[,(2*N^2 +2):(3*N^2 + 1)]
LP <- out[,(3*N^2 + 2):(4*N^2 + 1)]
ON <- out[,(4*N^2 + 2):(5*N^2 + 1)]
OP <- out[,(5*N^2 + 2):(6*N^2 + 1)]
NonFixers <- out[,(6*N^2 + 2):(6*N^2 + 1 + j)]
Fixers <- out[,(6*N^2 + 2 + j):(6*N^2 + 1 + j + k)]


makeGif <- function(time, data, min_time, name){
  new_data <-  
    cbind(time=time, data) %>%
    as.data.frame() %>% 
    pivot_longer(cols=-time, names_to="row", values_to="value") %>%
    mutate(
      row = as.numeric(row)-min(as.numeric(row))+1,
      x = x[(row-1) %% N + 1],
      y = y[(row-1) %/% N +1]) %>%
    filter(time>min_time)
  return(ggplot(new_data) +
    geom_tile(aes(x=x,y=y,fill=value)) +
    scale_fill_gradientn(limits=c(min(new_data$value),max(new_data$value)),
                         colors=rev(rainbow(7)[1:6])) +
    transition_time(time) + 
    labs(title = paste0(name),
         subtitle="Year: {as.integer(frame_time)}") +
    scale_x_continuous(expand=c(0,0), breaks=seq(0,N,5)) +
    scale_y_continuous(expand=c(0,0), breaks=seq(0,N,5)))
}

# data <- cbind(var=list(rep(c("AN", "AP", "LN", "LP"), time=rep(times,4),each=length(times))),rbind(AN, AP, LN, LP))
# data %>%
#   as.data.frame() %>% 
#   pivot_longer(cols=-time, names_to="row", values_to="value") %>% 
#   mutate(row = as.numeric(row),
#     x = x[(row-1) %% N + 1],
#     y = y[(row-1) %/% N +1]) %>%
#   filter(time > min_time)

makeGif(times, LN, 15, "litter N")
# makeGif(times, LP, 15, "litter P")
makeGif(times, AN, 15, "available N")
# makeGif(times, AP, 15, "available P")

# out[,1:(N^2+1)] %>%
#   as.data.frame %>% 
#   pivot_longer(cols=-time, names_to="row", values_to="value") %>%
#   mutate(
#     row = as.numeric(row),
#     x = x[(row-1) %% N + 1],
#     y = y[(row-1) %/% N +1]) %>%
#   ggplot() +
#   geom_tile(aes(x=x,y=y,fill=value)) +
#   scale_fill_gradientn(limits=c(min(out[, 2:(N^2+1)]),max(out[, 2:(N^2+1)])),
#                        colors=c("blue","green")) +
#   transition_time(time) + 
#   labs(title = "Available nitrogen pool over 20 years")
# 
# 

