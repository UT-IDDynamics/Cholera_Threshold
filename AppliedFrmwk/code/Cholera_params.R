
h <- 100 #number of stochastic runs #number of points
set.seed(3) #set a single seed so runif() is the same each time
seed <- runif(h)
#create additional empty storage objects needed for over prop.m.abx
# k<-seq(0,1.0,by=0.1)
j<-1
labels<-c()
S_0 <- 500000
lhs <- maximinLHS(n=h,k=22) #n=num of points or simulations. #k= number of variables or columns
# Draws a Latin Hypercube Sample from a set of uniform distributions 
#   for use in creating a Latin Hypercube Design. 
#   This function attempts to optimize the sample by maximizing the 
#   minium distance between design points (maximin criteria).
# R0.min <- 2.0 #reproductive number. use to calc beta (transmission rate)
# R0.max <- 2.5 #1.3-2.0; 2.3-2.8
sigma.min <- 0.62 #rate at which symptoms develop; 1/latent or incubation #X
sigma.max <- 0.77 #X
v.a.min <- 0.15 #modifier; reduction in transmission for asymptomatic vs severe symptomatic #X
v.a.max <- 0.35 #X
v.m.min <- 0.45 #X
v.m.max <- 0.75 #X
v.sh.min <- 0.3 #X
v.sh.max <- 0.6 #X
v.abx.min <- 0.3 #X
v.abx.max <- 0.7 #X

epsilon.a.min <- 0.65 #X
epsilon.a.max <- 0.85 #X
epsilon.s.min <- 0.2 #X
epsilon.s.max <- 0.4 #X
# epsilon.m.T.min <- 1 #X
# epsilon.m.T.max <- 1 #X
epsilon.s.T.min <- 0.5 #X
epsilon.s.T.max <- 0.9 #X

alpha.m.min <- 0.2 #X
alpha.m.max <- 0.4 #X
alpha.s.min <- 0.1 #X
alpha.s.max <- 0.3 #X
tau.min <- 1.33 #X
tau.max <- 4.0 #X
theta.m.min <- 0.05 #X
theta.m.max <- 0.3 #X
theta.s.min <- 0.4 #X
theta.s.max <- 0.8 #X

gamma.a.min <- 0.1 #X
gamma.a.max <- 0.3 #X
gamma.m.min <- 0.1 #X
gamma.m.max <- 0.25 #X
gamma.s.min <- 0.05 #LINDSAY CHANGED
gamma.s.max <- 0.2 #X
gamma.m.abx.min <- 1.0 #LINDSAY CHANGED
gamma.m.abx.max <- 1.0 #LINDSAY CHANGED
gamma.s.abx.min <- 1.0 #X
gamma.s.abx.max <- 1.0#X

mu.m.min <- 0.5#X
mu.m.max <- 1.0#X
mu.s.min <- 1.0#X
mu.s.max <- 2.0#X

