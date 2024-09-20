library(dplyr)
library(ggplot2)
library(lhs)
library(dplyr) #need for bind_rows for incidences_sum
library(reshape2) #need for melt for total.inc_long

# runSEIR takes the SEIR model parameters and initial condition 
# and returns a time series for each group
SEIIRfunct <- function(beta, sigma, v.a, v.m, v.sh, v.abx,
                       epsilon.a, epsilon.s, epsilon.m.T, epsilon.s.T,
                       alpha.m, alpha.s,
                       tau,q,delta,
                       theta.m, theta.s,
                       gamma.a, gamma.m, gamma.s,
                       gamma.m.abx, gamma.s.abx,
                       mu.m, mu.s,
                       initial.state,
                       step.size = 1,
                       freq.dependent=TRUE,
                       final.only=FALSE){
  
  #If the model is frequency dependent we modify beta based on the total population size
  beta.divisor <- ifelse(freq.dependent,
                         initial.state["S"]+
                           initial.state["E"]+
                           initial.state["Ia_sh"]+
                           initial.state["Im_syU"]+initial.state["Im_sh"]+
                           initial.state["Im_syT"]+initial.state["Im_abx"]+
                           initial.state["Is_syU"]+initial.state["Is_sh"]+
                           initial.state["Is_syT"]+initial.state["Is_abx"]+
                           initial.state["R"]+initial.state["R_abx"]+initial.state["D"],
                         1)
  
  #create the parameter vector.
  param <- c(beta.scaled=beta/beta.divisor, sigma=sigma, 
             v.a=v.a, v.m=v.m, v.sh=v.sh, v.abx=v.abx,
             epsilon.a=epsilon.a, epsilon.s=epsilon.s, epsilon.m.T=epsilon.m.T, epsilon.s.T=epsilon.s.T,
             alpha.m=alpha.m, alpha.s=alpha.s,
             tau=tau,q=q,delta=delta,
             theta.m=theta.m, theta.s=theta.s,
             gamma.a=gamma.a, gamma.m=gamma.m, gamma.s=gamma.s,
             gamma.m.abx=gamma.m.abx, gamma.s.abx=gamma.s.abx,
             mu.m=mu.m, mu.s=mu.s)
  
  #Since we are not using a fancy solver we will need to run this on our own.
  #note that the simulation ends once there are 0 people in groups E or I
  t <- 0
  y <- initial.state
  
  SEIIR.output <- matrix(ncol=length(initial.state)+1, nrow=1)
  colnames(SEIIR.output) <- c("time", "S", "E",
                              "Ia_sh",
                              "Im_syU", "Im_sh", "Im_syT", "Im_abx",
                              "Is_syU", "Is_sh", "Is_syT", "Is_abx",
                              "R", "D", "R_abx")
  SEIIR.output[1,] <- c(t,y)
  SEIIR.output
  incidences <- matrix(ncol=18+1, nrow=1)
  colnames(incidences) <- c("time",
                            "incident.exposed",
                            "incident.cases.a",
                            "incident.cases.symp.mU", "incident.cases.symp.mT",
                            "incident.cases.symp.sU", "incident.cases.symp.sT",
                            "shedding.cases.mU", "shedding.cases.mT", "abx.cases.m",
                            "shedding.cases.s", "abx.cases.s",
                            "nat.recovered.a",
                            "nat.recovered.m", "died.m", "abx.recovered.m",
                            "nat.recovered.s", "died.s", "abx.recovered.s")
  incidences[1,] <- rep(0,18+1)
  
  while(
    y["E"]>0 | 
    y["Ia_sh"]>0 |
    y["Im_syU"]>0 | y["Im_sh"]>0 | y["Im_syT"]>0 | y["Im_abx"]>0 |
    y["Is_syU"]>0 | y["Is_sh"]>0 | y["Is_syT"]>0 | y["Is_abx"]>0) {
    
    t <- t+step.size 
    
    lambda <- (param["beta.scaled"]*(y["Is_syU"] + y["Is_syT"])) + 
      (v.sh*param["beta.scaled"]*y["Is_sh"]) + 
      (v.a*param["beta.scaled"]*y["Ia_sh"]) + 
      ((v.m*param["beta.scaled"])*(y["Im_syU"]+y["Im_syT"])) + 
      (v.sh*v.m*param["beta.scaled"]*y["Im_sh"]) + 
      (v.m*v.abx*param["beta.scaled"]*y["Im_abx"])
    lambda
    #at this point, lambda is basically beta, looks complicated since is so many different infected types
    #is a number, should look like a number. print to look at
    
    #calculate the probability of infection and recovery in this time step
    #probability that thing happens in a given time step, updates every time step. range 0-1
    #need one of these per arrow, moderate/severe specific etc
    Pr.expose <- 1-exp(-step.size*lambda)
    
    Pr.infect.all <- 1-exp(-step.size*param["sigma"])
    
    Pr.mU.shed <- 1-exp(-step.size*param["alpha.m"]) #the rate at which moderates who don't seek care move to shed
    Pr.mT.shed <- 1-exp(-step.size*param["alpha.m"]) #this is only rate of move to shed
    Pr.mT.abx <- 1-exp(-step.size*param["delta"]*param["tau"]) #this is only rate of moderates who seek care move to abx
    # Pr.mT.shed <- 1-exp(-step.size*(1-param["q"])*param["alpha.m"])
    # Pr.mT.abx <- 1-exp(-step.size*param["q"]*param["delta"]*param["tau"])
    Pr.sU.shed <- 1-exp(-step.size*param["alpha.s"]) #the rate at which severes who don't seek care move to shed
    Pr.sT.abx <- 1-exp(-step.size*param["tau"]) #the rate at which severes who do seek care move to abx
    
    Pr.a.recover <- 1-exp(-step.size*param["gamma.a"]) #rate at which asymptomatics recover
    Pr.m.shed.recover <- 1-exp(-step.size*param["gamma.m"]) #this only incorporates the rate at which shedding moderates recover
    Pr.m.shed.die <- 1-exp(-step.size*param["mu.m"]) #this only incorporates the rate at which shedding moderates die
    # Pr.m.shed.recover <- 1-exp(-step.size*(1-param["theta.m"])*param["gamma.m"])
    # Pr.m.shed.die <- 1-exp(-step.size*param["theta.m"]*param["mu.m"])
    Pr.m.abx.recover <- 1-exp(-step.size*param["gamma.m.abx"]) #rate at which moderates who receive abx recover
    Pr.s.shed.recover <- 1-exp(-step.size*param["gamma.s"]) #this only incorporates the rate of recover among shedding severes
    Pr.s.shed.die <- 1-exp(-step.size*param["mu.s"]) #this only incorporates the rate of die among shedding severes
    # Pr.s.shed.recover <- 1-exp(-step.size*(1-param["theta.s"])*param["gamma.s"])
    # Pr.s.shed.die <- 1-exp(-step.size*param["theta.s"]*param["mu.s"])
    Pr.s.abx.recover <- 1-exp(-step.size*param["gamma.s.abx"]) #rate of recover among severes who receive abx
    
    #draw random variable from binomial distribution for new number of
    #using prob from above, is the actual number based on prob and denom (num ppl in that box)
    #should be whole numbers. if aren't something's wrong
    #rbinom(num observations, num trials, prob success each trial)
    incident.exposed <- rbinom(1, y["S"], Pr.expose) #the number who go S to E
    
    incident.cases.all <- rbinom(1, y["E"], Pr.infect.all) #this is the number that leave E
    incident.cases.a <- round(incident.cases.all*param["epsilon.a"],0) #this the number that enter I_a
    #now only have (number that leave)-(number that go to A) left to go to symptomatics
    incident.cases.symp <- (incident.cases.all-incident.cases.a) #this is the number that go to symptomatic I's
    incident.cases.symp.s <- round(incident.cases.symp*param["epsilon.s"],0) #this is the number that go to severe symptomatic I's
    incident.cases.symp.m <- (incident.cases.symp-incident.cases.symp.s) #this is the number that goes to moderate symptomatic I's
    incident.cases.symp.mT <- round(incident.cases.symp.m*param["epsilon.m.T"],0) #this is the number that go to moderate seek care
    incident.cases.symp.mU <- (incident.cases.symp.m-incident.cases.symp.mT) #this is the number that go to moderate don't seek care
    incident.cases.symp.sT <- round(incident.cases.symp.s*param["epsilon.s.T"],0) #this is the number that go to severe seek care
    incident.cases.symp.sU <- (incident.cases.symp.s-incident.cases.symp.sT) #this is the number that go to severe don't seek care
    # incident.cases.all
    # sum(incident.cases.a,incident.cases.symp.mU,incident.cases.symp.mT,incident.cases.symp.sU,incident.cases.symp.sT)
    
    shedding.cases.mU <- rbinom(1, y["Im_syU"], Pr.mU.shed) #the number who move to moderate shedding
    m.who.will.abx <- round(y["Im_syT"]*param["q"],0) #the number of moderates who seek care who will eventually receive abx#are treating q like a proportion, though this not actually correct. the math works though and conceptually is ok via prop.m.abx equation
    shedding.cases.mT <- rbinom(1, y["Im_syT"]-m.who.will.abx, Pr.mT.shed) #the number of moderates who seek care who don't receive abx and move to shed
    abx.cases.m <- rbinom(1, m.who.will.abx, Pr.mT.abx) #the number of moderates who receive abx at this step
    # shedding.cases.mT <- rbinom(1, y["Im_syT"], Pr.mT.shed)
    # abx.cases.m <- rbinom(1, y["Im_syT"], Pr.mT.abx)
    shedding.cases.s <- rbinom(1, y["Is_syU"], Pr.sU.shed) #the number of severes who don't seek care and move to shed
    abx.cases.s <- rbinom(1, y["Is_syT"], Pr.sT.abx) #the number of severes who seek care and receive abx
    
    nat.recovered.a <- rbinom(1, y["Ia_sh"], Pr.a.recover) #the number of asymptomatics who recover
    m.who.will.die <- round(y["Im_sh"]*param["theta.m"],0) #this is the whole number of moderates who will eventually die
    nat.recovered.m <- rbinom(1, y["Im_sh"]-m.who.will.die, Pr.m.shed.recover) #this is the actual number of I_sh who naturally recover. this incorporates the stochastic prob that will actually move in this step #rbinom requires whole numbers
    died.m <- rbinom(1, m.who.will.die, Pr.m.shed.die) #the number of moderates who move from I_sh to dead
    # nat.recovered.m <- rbinom(1, y["Im_sh"], Pr.m.shed.recover)
    # died.m <- rbinom(1, y["Im_sh"], Pr.m.shed.die)
    abx.recovered.m <- rbinom(1, y["Im_abx"], Pr.m.abx.recover) #number of moderates who move from abx to recover in this step
    s.who.will.die <- round(y["Is_sh"]*param["theta.s"],0) #this is the whole number of severes who will die
    nat.recovered.s <- rbinom(1, y["Is_sh"]-s.who.will.die, Pr.s.shed.recover) #number of severes who move from Is_sh to recover w/o abx
    died.s <- rbinom(1, s.who.will.die, Pr.s.shed.die) #number of Is_sh who die in this time step
    abx.recovered.s <- rbinom(1, y["Is_abx"], Pr.s.abx.recover) #number of severes who recover after abx in this time step
    
    #Find the deltas for each compartment
    dS <- -incident.exposed
    dE <- incident.exposed - incident.cases.a - 
      incident.cases.symp.mU - incident.cases.symp.mT - 
      incident.cases.symp.sU - incident.cases.symp.sT
    dIa.sh <- incident.cases.a - nat.recovered.a
    
    dIm.syU <- incident.cases.symp.mU - shedding.cases.mU
    dIm.sh <- shedding.cases.mU + shedding.cases.mT - nat.recovered.m - died.m
    dIm.syT <- incident.cases.symp.mT - shedding.cases.mT - abx.cases.m
    dIm.abx <- abx.cases.m - abx.recovered.m
    
    dIs.syU <- incident.cases.symp.sU - shedding.cases.s
    dIs.sh <- shedding.cases.s - nat.recovered.s - died.s
    dIs.syT <- incident.cases.symp.sT - abx.cases.s
    dIs.abx <- abx.cases.s - abx.recovered.s
    
    dR <- nat.recovered.a + nat.recovered.m + nat.recovered.s
    dD <- died.m + died.s
    dR.abx <- abx.recovered.m + abx.recovered.s
    dE
    sum(dIa.sh,dIm.syU,dIm.syT,dIs.syU,dIs.syT)
    
    deltas <- c(dS,dE,dIa.sh,dIm.syU,dIm.sh,dIm.syT,dIm.abx,dIs.syU,dIs.sh,dIs.syT,dIs.abx,dR,dD,dR.abx) # calculate step sizes. this is net change for each box
    z <- c(incident.exposed,
           incident.cases.a,
           incident.cases.symp.mU, incident.cases.symp.mT,
           incident.cases.symp.sU, incident.cases.symp.sT,
           shedding.cases.mU, shedding.cases.mT, abx.cases.m,
           shedding.cases.s, abx.cases.s,
           nat.recovered.a,
           nat.recovered.m, died.m, abx.recovered.m,
           nat.recovered.s, died.s, abx.recovered.s)
    #this is the incident (new) number of ppl along each arrow at each step
    
    y <- y+deltas 
    # print(c(t,y))
    
    #trick to speed up the code
    if (!final.only){
      SEIIR.output<-rbind(SEIIR.output, c(t,y))
      incidences<-rbind(incidences, c(t,z))
    }
  }
  
  if(final.only){
    SEIRR.output[1,]<-c(t,y)
    incidences[1,]<-c(t,z)
  }
  
  combo<-list(SEIIR.output=as.data.frame(SEIIR.output),incidences=as.data.frame(incidences))
  return(combo)
}
