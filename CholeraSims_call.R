###########################################################################
# Quantifying the impact of expanded antibiotic treatment on cholera 
# outbreaks
# Code authors: Sharia M. Ahmed, Lindsay T. Keegan
# September 17, 2024
###########################################################################


# Load Libraries ----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(lhs)
library(dplyr) 
library(reshape2)
library(stringr)


# Define Functions --------------------------------------------------------

## create.variable.params: function that takes the parameter values we vary and
## puts them into a dataframe that can be used to loop over when calling the model
create.variable.params <- function(R0.min.values, R0.max.values, epsilon.values) {
  # Create a data frame with all combinations of the parameters
  variable.params <- expand.grid(
    R0.min.values = R0.min.values,
    epsilon.values = epsilon.values
  )
  
  # Assign R0.max.values with the correct R0.min.values
  variable.params <- variable.params %>%
    mutate(R0.max.values = R0.max.values[match(R0.min.values, R0.min.values)])
      
  # Use mutate and case_when to create the categorical columns
  variable.params <- variable.params %>%
    mutate(
      R0.categorical = case_when(
        R0.min.values <= 1.4 ~ "low",
        R0.min.values <= 2.0 ~ "intermediate",
        R0.min.values > 2.0 ~ "high"
      ),
      treatment.categorical = case_when(
        epsilon.values <= 0.1 ~ "verylow",
        epsilon.values > 0.1 & epsilon.values <= 0.3 ~ "low",
        epsilon.values > 0.3 & epsilon.values <= 0.6 ~ "intermediate",
        epsilon.values > 0.6 & epsilon.values <= 0.99 ~ "high",
        epsilon.values > 0.99 ~ "perfect"
      )
    )
  
  return(variable.params)
}

## run.model: function that calls the SEIR model and runs LHS for given parameter ranges
## Cholera_params.R is the script with all the unchanging parameters
## run.model takes in a range for R0, the values of proportion of moderates 
## treated (k), and the range for the proportion of moderates who seek treatment 
## (epsilon)
run.model <- function(R0.min, R0.max, k, epsilon.m.T.min, epsilon.m.T.max) {

  source("Cholera_Function.R")
  source("Cholera_params.R")

  lhs.params <- cbind(
    R0 = lhs[,1]*(R0.max-R0.min)+R0.min,
    sigma = lhs[,2]*(sigma.max-sigma.min)+sigma.min,
    v.a = lhs[,3]*(v.a.max-v.a.min)+v.a.min,
    v.m = lhs[,4]*(v.m.max-v.m.min)+v.m.min,
    v.sh = lhs[,5]*(v.sh.max-v.sh.min)+v.sh.min,
    v.abx = lhs[,6]*(v.abx.max-v.abx.min)+v.abx.min,
    
    epsilon.a = lhs[,7]*(epsilon.a.max-epsilon.a.min)+epsilon.a.min,
    epsilon.s = lhs[,8]*(epsilon.s.max-epsilon.s.min)+epsilon.s.min,
    epsilon.m.T = lhs[,9]*(epsilon.m.T.max-epsilon.m.T.min)+epsilon.m.T.min,
    epsilon.s.T = lhs[,10]*(epsilon.s.T.max-epsilon.s.T.min)+epsilon.s.T.min,
    
    alpha.m = lhs[,11]*(alpha.m.max-alpha.m.min)+alpha.m.min,
    alpha.s = lhs[,12]*(alpha.s.max-alpha.s.min)+alpha.s.min,
    tau = lhs[,13]*(tau.max-tau.min)+tau.min,
    theta.m = lhs[,14]*(theta.m.max-theta.m.min)+theta.m.min,
    theta.s = lhs[,15]*(theta.s.max-theta.s.min)+theta.s.min,
    
    gamma.a = lhs[,16]*(gamma.a.max-gamma.a.min)+gamma.a.min,
    gamma.m = lhs[,17]*(gamma.m.max-gamma.m.min)+gamma.m.min,
    gamma.s = lhs[,18]*(gamma.s.max-gamma.s.min)+gamma.s.min,
    gamma.m.abx = lhs[,19]*(gamma.m.abx.max-gamma.m.abx.min)+gamma.m.abx.min,
    gamma.s.abx = lhs[,20]*(gamma.s.abx.max-gamma.s.abx.min)+gamma.s.abx.min,
    
    mu.m = lhs[,21]*(mu.m.max-mu.m.min)+mu.m.min,
    mu.s = lhs[,22]*(mu.s.max-mu.s.min)+mu.s.min
  )
  
  # Initialize lists to store outputs for each value of k
  SEIR.results <- list()
  INCID.results <- list()
  
  # Loop over k
  for (l in 1:length(k)) {
    prop.abx <- k[l]
    
    # Initialize dataframes to accumulate results for the current k
    SEIR.output.final <- NULL
    INCID.output.final <- NULL
    
    for (i in 1:h) {
      # Define initial state and parameters (same as your original code)
      n.seed.events <- 50
      initial.state <- c(S=S_0-(n.seed.events), 
                         E=n.seed.events,
                         Ia_sh=0, Im_syU=0, Im_sh=0, Im_syT=0, Im_abx=0, 
                         Is_syU=0, Is_sh=0, Is_syT=0, Is_abx=0, 
                         R=0, D=0, R_abx=0)
      set.seed <- seed[[i]]
      
      # Define SEIR parameters
      R0 <- as.numeric(lhs.params[i,"R0"])
      sigma <- as.numeric(lhs.params[i,"sigma"])
      v.a <- as.numeric(lhs.params[i,"v.a"])
      v.m <- as.numeric(lhs.params[i,"v.m"])
      v.sh <- as.numeric(lhs.params[i,"v.sh"])
      v.abx <- as.numeric(lhs.params[i,"v.abx"])
      epsilon.a <- as.numeric(lhs.params[i,"epsilon.a"])
      epsilon.s <- as.numeric(lhs.params[i,"epsilon.s"])
      epsilon.m.T <- as.numeric(lhs.params[i,"epsilon.m.T"])
      epsilon.s.T <- as.numeric(lhs.params[i,"epsilon.s.T"])
      alpha.m <- as.numeric(lhs.params[i,"alpha.m"])
      alpha.s <- as.numeric(lhs.params[i,"alpha.s"])
      tau <- as.numeric(lhs.params[i,"tau"])
      theta.m <- as.numeric(lhs.params[i,"theta.m"])
      theta.s <- as.numeric(lhs.params[i,"theta.s"])
      gamma.a <- as.numeric(lhs.params[i,"gamma.a"])
      gamma.m <- as.numeric(lhs.params[i,"gamma.m"])
      gamma.s <- as.numeric(lhs.params[i,"gamma.s"])
      gamma.m.abx <- as.numeric(lhs.params[i,"gamma.m.abx"])
      gamma.s.abx <- as.numeric(lhs.params[i,"gamma.s.abx"])
      mu.m <- as.numeric(lhs.params[i,"mu.m"])
      mu.s <- as.numeric(lhs.params[i,"mu.s"])
      
      #calculate beta given R0 etc
      beta <- R0/(
        ((epsilon.a*v.a)/(gamma.a))+
          (((((1-epsilon.a)*(1-epsilon.s)*(1-epsilon.m.T))*v.m)/(alpha.m))+
             ((((1-epsilon.a)*(1-epsilon.s)*(1-epsilon.m.T))*v.sh*v.m)/(((1-theta.m)*gamma.m)+(theta.m*mu.m)))+
             ((((1-epsilon.a)*(1-epsilon.s)*epsilon.m.T)*v.m)/(alpha.m))+
             ((((1-epsilon.a)*(1-epsilon.s)*epsilon.m.T)*v.sh*v.m)/(((1-theta.m)*gamma.m)+(theta.m*mu.m))))+
          ((((1-epsilon.a)*epsilon.s*(1-epsilon.s.T)))/(alpha.s))+
          (((1-epsilon.a)*epsilon.s*epsilon.s.T)/(tau))+
          ((((1-epsilon.a)*epsilon.s*(1-epsilon.s.T))*v.sh)/(((1-theta.s)*gamma.s)+(theta.s*mu.s)))+
          ((((1-epsilon.a)*epsilon.s*epsilon.s.T)*v.abx)/(gamma.s.abx))
      ) 
      
      delta <- 0.5
      q <- (alpha.m * prop.abx) / ((alpha.m * prop.abx) - (delta * tau * prop.abx) + (delta * tau))
      
      # Run the SEIIR model
      temp.outputs <- SEIIRfunct(beta, sigma, v.a, v.m, v.sh, v.abx,
                                 epsilon.a, epsilon.s, epsilon.m.T, epsilon.s.T,
                                 alpha.m, alpha.s,
                                 tau, q, delta,
                                 theta.m, theta.s,
                                 gamma.a, gamma.m, gamma.s,
                                 gamma.m.abx, gamma.s.abx,
                                 mu.m, mu.s,
                                 initial.state,
                                 step.size = 1,
                                 freq.dependent = TRUE,
                                 final.only = FALSE)
      
      # Extract outputs and add prop.abx column
      temp.seir.out <- as.data.frame(temp.outputs$SEIIR.output)
      temp.seir.out$run <- i
      temp.seir.out$prop.m.abx <- prop.abx
      
      temp.incidences.out <- as.data.frame(temp.outputs$incidences)
      temp.incidences.out$run <- i
      temp.incidences.out$prop.m.abx <- prop.abx
      
      # Accumulate the results for this value of k
      SEIR.output.final <- rbind(SEIR.output.final, temp.seir.out)
      INCID.output.final <- rbind(INCID.output.final, temp.incidences.out)
    }
    
    # Store the results for this value of k in the lists
    SEIR.results[[l]] <- SEIR.output.final
    INCID.results[[l]] <- INCID.output.final
    
  }
  
  # Combine all the results into a single dataframe for further analysis
  SEIR.final <- do.call(rbind, SEIR.results)
  INCID.final <- do.call(rbind, INCID.results)
  return(list(SEIR.final = SEIR.final, INCID.final = INCID.final))
  
}


## summarize.incidences: function that takes the output of run_model and summarizes the
## output (mean, quantiles)
summarize.incidences <- function(INCID.final) {
  
  # Sum values in each column except prop.m.abx, grouped by run
  INCID.summed <- INCID.final %>%
    group_by(run, prop.m.abx) %>%
    summarize(across(-time, sum))
  
  # Add additional calculated columns
  INCID.summed <- INCID.summed %>%
    mutate(untreated.cases.m = shedding.cases.mU + shedding.cases.mT) %>%
    mutate(incident.cases.symp.m = incident.cases.symp.mU + incident.cases.symp.mT) %>%
    mutate(incident.cases.symp.s = incident.cases.symp.sU + incident.cases.symp.sT) %>%
    mutate(incident.cases.all = incident.exposed)
  
  # Calculate mean and quantiles for each column, grouped by prop.m.abx
  INCID.summary <- INCID.summed %>%
    group_by(prop.m.abx) %>%
    summarize(across(-run, 
                     list(mean = mean, q25 = ~ quantile(.x, 0.25), 
                          q50 = ~ quantile(.x, 0.5), q75 = ~ quantile(.x, 0.75)), 
                     .names = "{.col}_{.fn}"))
  
  # Pivot the summary data into long format
  INCID.summary.long <- INCID.summary %>%
    pivot_longer(
      cols = -prop.m.abx, 
      names_to = c("compartment", "statistic"), 
      names_sep = "_", 
      values_to = "value"
    )
  
  # Filter for specific compartment values
  filtered.INCID <- INCID.summary.long %>%
    filter(compartment %in% c("incident.cases.a", 
                              "incident.cases.symp.m",
                              "untreated.cases.m", 
                              "abx.cases.m", 
                              "incident.cases.symp.sU", 
                              "incident.cases.symp.sT", 
                              "incident.cases.symp.s",
                              "incident.cases.all"))
  
  # Pivot summary to wide format for easier ggplot usage with geom_ribbon
  INCID.summary.wide <- filtered.INCID %>%
    pivot_wider(names_from = statistic, values_from = value)
  
  return(INCID.summary.wide)
}

## retrieve.datasets: function that loops over the names of datasets and extracts them
retrieve.datasets <- function(names) {
  # Create an empty list to store the datasets
  datasets <- list()
  
  # Loop through the summary.names and use get() to retrieve the objects
  for (i in seq_along(names)) {
    datasets[[names[i]]] <- get(names[i])
  }
  
  # Return the list of datasets
  return(datasets)
}

## combine.incident.cases: function that takes in multiple run.model simulations and 
## combines them, extracting only the total number of cases to be used in making 
## figures. Takes in all of the datasets that need to be combined as well as the 
## corresponding R0 and treatment seeking values
combine.incident.cases <- function(datasets, R0.values, seek.treatment.values) {
  # Check if the lengths of inputs match
  if (length(datasets) != length(R0.values) || length(datasets) != length(seek.treatment.values)) {
    stop("The number of datasets, R0 values, and seek treatment values must match.")
  }
  
  # Create an empty list to store the modified datasets
  modified_datasets <- list()
  
  # Loop through each dataset and corresponding R0 and seek_treatment values
  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    R0.value <- R0.values[[i]]
    seek.treatment.value <- seek.treatment.values[[i]]
    
    # Mutate to add R0 and seek treatment columns
    modified.dataset <- dataset %>%
      mutate(R0 = R0.value,
             seektreatment = seek.treatment.value) %>%
      # Filter for rows where compartment == "incident.cases.all"
      filter(compartment == "incident.cases.all")
    
    # Add the modified dataset to the list
    modified_datasets[[i]] <- modified.dataset
  }
  
  # Combine all modified datasets into one
  combined.dataset <- bind_rows(modified_datasets)
  
  # Return the combined dataset
  return(combined.dataset)
}

## combine.incident.cases.all: function that takes in multiple run.model simulat0ions and 
## combines them and makes the formatting wide (for plotting) of cases to be used in  
## making figures. Takes in all of the datasets that need to be combined as well 
## as the corresponding R0 and treatment seeking values
combine.incident.cases.all <- function(datasets, R0.values, seek.treatment.values) {
  # Check if the lengths of inputs match
  if (length(datasets) != length(R0.values) || length(datasets) != length(seek.treatment.values)) {
    stop("The number of datasets, R0 values, and seek treatment values must match.")
  }
  
  # Create an empty list to store the modified datasets
  modified.datasets <- list()
  
  # Loop through each dataset and corresponding R0 and seek treatment values
  for (i in seq_along(datasets)) {
    dataset <- datasets[[i]]
    R0.value <- R0.values[[i]]
    seek.treatment.value <- seek.treatment.values[[i]]
    
    # Mutate to add R0 and seek treatment columns
    modified.dataset <- dataset %>%
      mutate(R0 = R0.value,
             seektreatment = seek.treatment.value)
    
    # Add the modified dataset to the list
    modified.datasets[[i]] <- modified.dataset
  }
  
  # Combine all modified datasets into one
  combined.dataset <- bind_rows(modified.datasets)
  

  # Return the combined dataset
  return(combined.dataset)
}

## averted.incidence.processing: function to compare the expanded guidelines to the same
## parameter run under current guidlines. Calculates averted incidence and additional doses
averted.incidence.processing <- function(dataset) {
  # Sum all columns for each run and prop.m.abx combination, except for 'time'
  INCID.summed <- dataset %>%
    group_by(run, prop.m.abx) %>%
    summarize(across(-time, sum, na.rm = TRUE))  # Sum all columns except 'time'
  
  # Extract values where prop.m.abx == 0 for each run
  INCID.base <- INCID.summed %>%
    filter(prop.m.abx == 0) %>%
    select(-prop.m.abx)  # Drop the prop.m.abx column as it's not needed in the base
  
  # Rename columns to indicate they are base values (for clarity)
  colnames(INCID.base)[-1] <- paste0(colnames(INCID.base)[-1], ".base")
  
  # Join the base values with the original dataset
  INCID.adjusted <- INCID.summed %>%
    left_join(INCID.base, by = "run") %>%
    mutate(across(ends_with(".base"), ~ get(sub("\\.base$", "", cur_column())) - .x, .names = "adjusted.{.col}")) %>%
    select(run, prop.m.abx, starts_with("adjusted."))
  
  INCID.adjusted <- INCID.adjusted %>%
    mutate(total.abx.cases = adjusted.abx.cases.m.base + adjusted.abx.cases.s.base) %>%
    mutate(adjusted.incident.exposed.base = -adjusted.incident.exposed.base) %>%
    mutate(inf.averted = adjusted.incident.exposed.base / S_0) %>%
    mutate(addnl.abx = total.abx.cases / S_0) 
  
  # Create the color group column 
  INCID.adjusted <- INCID.adjusted %>%
    mutate(color.group = case_when(
      inf.averted < addnl.abx ~ "between.zero.and.line",  # Points between y = 0 and y = x
      inf.averted >= addnl.abx & addnl.abx > 0 ~ "between.line.and.zero",  # Points between y = x and x = 0
      addnl.abx < 0 ~ "x.below.zero"  # Points where x < 0
    ))
  
  return(INCID.adjusted)
}

# Cholera Final Size  -----------------------------------------------------
## Run the model -----------------------------------------------------------

R0.min.values <- c(1.2, 1.5, 2.3)  # R0 min values
R0.max.values <- c(1.4, 2.0, 2.8)  # R0 max values
epsilon.values <- c(0.05, 0.25,0.5, 0.75, 1)  # Epsilon values

variable.params <- create.variable.params(R0.min.values, R0.max.values, epsilon.values)

# Initialize output.names
output.names <- NULL

# Loop through all rows of variable.params (which now contains all combinations of R0 and epsilon)
for (i in 1:nrow(variable.params)) {
  
  # Dynamically create the output name using the values from the data frame
  output.name <- paste("outputR.", variable.params$R0.categorical[i], "R0.", 
                       variable.params$treatment.categorical[i], "treat", sep = "")
  
  # Use assign to run the model and store it under the dynamically generated name
  assign(output.name, 
         run.model(R0.min = variable.params$R0.min.values[i], 
                   R0.max = variable.params$R0.max.values[i], 
                   k = seq(0, 1.0, by = 0.1), 
                   epsilon.m.T.min = 0.5, #variable.params$epsilon.values[i], 
                   epsilon.m.T.max = variable.params$epsilon.values[i])
  )
  
  # Add the generated name to the output.names vector
  output.names <- c(output.names, output.name)
}


## Summarize the data ------------------------------------------------------
summary.names <- NULL

for(i in 1:length(output.names)){
  # Create the output summary name using the values from the output.names
  summary.name <- paste("INCID.summary", str_extract(output.names[i], "(?<=outputR\\.)[^R0]+"), "R0.", 
                        str_extract(output.names[i], "(?<=R0\\.)[a-zA-Z]+(?=treat)"), "treat", sep = "")
  
  assign(summary.name,
         summarize.incidences(INCID.final = get(output.names[i])$INCID.final))
  
  # Add the generated name to the output.names vector
  summary.names <- c(summary.names, summary.name)
}

datasets <- retrieve.datasets(summary.names)


R0.values<-variable.params$R0.categorical
seek.treatment.values <- variable.params$epsilon.values

## Combine the datasets 
combined.dataset <- combine.incident.cases(datasets, R0.values, seek.treatment.values)



## Final size plot ---------------------------------------------------------

ggplot(combined.dataset, 
       aes(x = prop.m.abx, y = mean/S_0, group = R0, color = R0)) +
  geom_ribbon(aes(ymin = q25/S_0, ymax = q75/S_0, fill = R0), alpha = 0.4, colour = NA) +
  geom_line(aes(linetype = "expanded"), size = 1) +
  geom_hline(data = combined.dataset %>% filter(prop.m.abx ==0),
             aes(yintercept = mean/S_0, group = R0, color = R0, linetype = "current"), size =1) +
  facet_wrap(~ factor(seektreatment, levels=c("0.05", "0.25", "0.5", "0.75", "1")), 
             ncol = 4, labeller=labeller(seektreatment=labels)) +
  labs(x = "Proportion of moderates treated", 
       y = "Proportion of total population infected") +
  theme_bw()  +
  labs(color = expression(R[0]), 
       fill = expression(R[0])) +
  scale_color_manual(values = c("low" = "#FAC484FF", 
                                "intermediate" = "#CE6693FF",
                                "high" = "#5C53A5FF")) +
  scale_fill_manual(values = c("low" = "#FAC484FF", 
                               "intermediate" = "#CE6693FF",
                               "high" = "#5C53A5FF")) +
  scale_linetype_manual(values = c("expanded" = "solid", "current" = "dashed"), 
                        labels = c("current" = "Current", "expanded" = "Expanded"), 
                        name = "Treatment guidelines") + 
  theme(strip.background=element_rect(colour="black",
                                      fill="white"), 
        strip.placement = "outside") +
  theme(panel.grid.major = element_blank()) 




# Relationship between treatment seeking and proportion treated -----------

### Using the same model runs as above
## Summarize the data ------------------------------------------------------

INCID.summary.all <-combine.incident.cases.all(datasets, R0.values, seek.treatment.values)

# Pivot wider 

INCID.summary.wider.all <- INCID.summary.all %>%
  pivot_wider(names_from = "compartment", values_from = c("mean", "q25", "q50", "q75")) %>%
  mutate(
    # Proportion of population treated
    prop_population_treated = (mean_abx.cases.m + mean_incident.cases.symp.sT) / S_0,
    
    # Proportion of moderates treated
    prop_moderates_treated = mean_abx.cases.m / S_0
  )

## Order the factor levels
INCID.summary.wider.all <- INCID.summary.wider.all %>%
  mutate(seektreatment_ordered = factor(seektreatment, levels = c("0.05", "0.25", "0.5","0.75", "1")))


# Proportion of the population treated plot -------------------------------

ggplot(INCID.summary.wider.all, aes(x = prop.m.abx,  
                                    y = prop_population_treated, 
                                    group = interaction(R0, seektreatment_ordered), 
                                    color = seektreatment_ordered)) +
  geom_line(size = 1) +
  theme_bw() +
  xlab("Proportion of moderates treated") + 
  facet_wrap(~ factor(R0, levels=c("low", "intermediate", "high")), 
             ncol = 4, labeller=labeller(R0=labels)) +
  scale_color_manual(values = c("0.05" = "#D8D97AFF", 
                                "0.25" = "#95C36EFF",
                                "0.5" ="#74C8C3FF", 
                                "0.75" = "#5A97C1FF", 
                                "1" = "#0A2E57FF"),
                     name = "Seeking\ntreatment") +
  ylab("Proportion of the population treated") + 
  theme(strip.background=element_rect(colour="black",
                                      fill="white"), 
        text = element_text(size=16),
        strip.placement = "outside") +
  theme(panel.grid.major = element_blank()) 



# Additional Doses Used ---------------------------------------------------

## Run the model for new values of R0
R0.min.values <- c(1.2, 1.5, 2.3)  # R0 min values
R0.max.values <- c(1.4, 2.0, 2.8)  # R0 max values
epsilon.values <- c(1)  # Epsilon values

variable.params <- create.variable.params(R0.min.values, R0.max.values, epsilon.values)

# Initialize output.names
output.quarter.names <- NULL

# Loop through all rows of variable.params (which now contains all combinations of R0 and epsilon)
for (i in 1:nrow(variable.params)) {
  
  # Dynamically create the output name using the values from the data frame
  output.quarter.name <- paste("outputR.quarter.", variable.params$R0.categorical[i], "R0.", 
                       variable.params$treatment.categorical[i], "treat", sep = "")
  
  # Use assign to run the model and store it under the dynamically generated name
  assign(output.quarter.name, 
         run.model(R0.min = variable.params$R0.min.values[i], 
                   R0.max = variable.params$R0.max.values[i], 
                   k = c(0, 0.05, 0.25,0.5, 0.75, 1), 
                   epsilon.m.T.min = variable.params$epsilon.values[i], 
                   epsilon.m.T.max = variable.params$epsilon.values[i])
  )
  
  # Add the generated name to the output.names vector
  output.quarter.names <- c(output.quarter.names, output.quarter.name)
}


quarter.names <- NULL

for(i in 1:length(output.quarter.names)){
  # Dynamically create the output summary name using the values from the output.quarter.names
  quarter.name <- paste("INCID.quarter", str_extract(output.quarter.names[i], "(?<=outputR\\.)[^R0]+"), "R0.", 
                        str_extract(output.quarter.names[i], "(?<=R0\\.)[a-zA-Z]+(?=treat)"), "treat", sep = "")
  
  assign(quarter.name,
         averted.incidence.processing(dataset = get(output.quarter.names[i])$INCID.final))
  
  # Add the generated name to the output.quarter.names vector
  quarter.names <- c(quarter.names, quarter.name)
}


quarter.datasets <- retrieve.datasets(quarter.names)

R0.values<-variable.params$R0.categorical
seek.treatment.values <- variable.params$epsilon.values


combined.adjusted.all <-combine.incident.cases.all(quarter.datasets, R0.values, seek.treatment.values)



# Plot of additional doses used --------------------------------------------


ggplot(combined.adjusted.all %>%
         filter(prop.m.abx != 0), aes(x = addnl.abx, y = inf.averted, color = color.group)) +
  geom_abline(slope=1, intercept=0,linetype="dashed", size = 0.9, color = "black") +
  geom_vline(xintercept = 0,linetype="dashed", size = 0.9, color = "black") +
  geom_hline(yintercept = 0,linetype="dashed", size = 0.9, color = "black") +
  geom_point(size = 3) +
  xlab("Additional antibiotic doses used") +
  ylab("Infections averted") +
  xlim(-0.05, 0.3) + 
  scale_color_manual(values = c("below.zero" = "#FDE725FF", 
                                "between.zero.and.line" = "#FDE725FF", 
                                "between.line.and.zero" = "#20A387FF", 
                                "x.below.zero" = "#481567FF"), 
                     labels = c("below.zero" = "Only individual benefits",
                                "between.zero.and.line" = "\nOnly individual benefits\n",
                                "between.line.and.zero" = "\nEach dose prevents\n>1 cases", 
                                "x.below.zero" = "Fewer doses used"),
                     name = "Expanded eligibility\nimpact") + 
  theme_bw() +
  theme(text=element_text(size=20))+
  facet_grid(factor(R0, 
                    levels=c("low", "intermediate", "high")) ~ factor(prop.m.abx, 
                                                                      levels = c("0.05", "0.25", "0.5", "0.75", "1")), 
             labeller=labeller(prop.m.abx=labels)) + 
  theme_bw(base_size=13) + 
  theme(strip.background=element_rect(colour="black",
                                      fill="white"), 
        strip.placement = "outside") +
  theme(panel.grid.major = element_blank(), 
        legend.position = "bottom") 

