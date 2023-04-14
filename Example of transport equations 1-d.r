library("ReacTran")

# we make our grid
L <- 20;
N <- 100;
Grid <- setup.grid.1D(-L/2, L/2, N = N)

# we define a property of the grid cells
# such as the rooting density

tree <- 2;
rMax <- 1;


# we are gonna make a pde which solves the equation with diffusion and advection
# and has the plant suck up nutrients according to rooting
pdeAvailable <- function(t, y, parms, A = 1) {
  end <- y[N] #this enforces the 
  tran <- tran.1D(C = y, A = A, D = D, v = ux, C.down = end, C.up = end, dx = Grid)$dC
  g <- Grid$x.mid
  reac <- -1*(0 + (((rMax^2-(Grid$x.mid -tree)^2))> 0))*1/rMax^2 * (rMax^2-(Grid$x.mid -tree)^2)*v*y + I
  return(list(tran + reac))
}


#initial condition for available nitrogen
ini <- rep(1,N)
D <-0.3
v <- 0.4
I <- 0.01
ux <--0.1
#now we solve the PDE equation

require(deSolve)
times <- seq (from = 0, to = 100, by = 1)
system.time(
  out <- ode.1D(y = ini, times = times, func = pdeAvailable, parms = NULL, nspec = 1)
)
image(out, grid = Grid$x.mid, xlab = "time, days",
      ylab = "Distance, cm", main = "PDE", add.contour = TRUE)