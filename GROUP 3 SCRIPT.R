## Load the package

library(deSolve)
library(tidyverse)
library(ggplot2)
library(DiagrammeR)
library(plotly)
library(dplyr)
library(ggplotify)
library(gridExtra)
# flow diagram-------------------

graph <- grViz("
digraph model {
  graph [layout = dot, rankdir = LR]

  # Human compartments
  Sh [label = 'Sh\nSusceptible Humans', shape = box, style=filled, fillcolor=lightblue]
  Ih [label = 'Ih\nInfected Humans', shape = box, style=filled, fillcolor=lightpink]

  # Mosquito compartments
  Sm [label = 'Sm\nSusceptible Mosquitoes', shape = box, style=filled, fillcolor=lightgreen]
  Em [label = 'Em\nExposed Mosquitoes', shape = box, style=filled, fillcolor=khaki]
  Im [label = 'Im\nInfected Mosquitoes', shape = box, style=filled, fillcolor=orange]

  # Birth/Death node
  BD [label = 'Birth / Death', shape = ellipse, style=filled, fillcolor=grey90]

  # Human transitions
  Sh -> Ih [label = 'λh * Sh']
  Ih -> Sh [label = 'γ * Ih']

  # Mosquito transitions
  Sm -> Em [label = 'λm * Sm']
  Em -> Im [label = 'σ * Em']
  Im -> BD [label = 'μm * Im']
  Em -> BD [label = 'μm * Em']
  Sm -> BD [label = 'μm * Sm']
  BD -> Sm [label = 'birth']

}
")

graph

#====================================================================
#Baseline prevalence------------------------------------------------
#====================================================================
# baseline without intervention 

Nh<- 15000 #total number of people 
Ih<- 25    # number of infected people 
Nm<-Nh*8   # 8 mosquito per peerson 
Im<-0.12*Nm # 12% of the mosquito are infectious

# mosquito cxcs

bitting<- 0.3
seeking<- 0.1
transmission<-0.45

alpha= bitting*seeking*transmission

# parameters 
parameters<-list(alpha= alpha,
                 beta= 0.5,
                 gamma= 1/25,
                 sigma= 1/7,
                 mu_m= 1/14,
                 Nh= 15000,
                 Nm= 15000*8
)

# initial state
initial_state_values <- c(
    Sh = (parameters$Nh - 25),
    Ih = 25,
    Sm = (parameters$Nh*8 - 0.12*parameters$Nh),
    Em = 0,
    Im = (0.12*parameters$Nh)
)


## Time and solve

times <- seq(0, 365, by =1)  # days

#model function 
base_model<- function(time, state, parameters){
    with(as.list(c(state, parameters)),{
        
        Nh <- Sh + Ih   # human population
        Nm <- Sm +Em + Im  # mosquito population
        
        # Forces of infection (Ross-Macdonald style)
        
        lambda_h <- alpha * (Im / Nh)       # per-susceptible-human per day
        lambda_m <- beta * (Ih / Nh)       # per-susceptible-mosquito per day
        
        dSh <- - lambda_h * Sh + gamma * Ih
        dIh <-   lambda_h * Sh - gamma * Ih
        
        dSm <- mu_m *Nm-lambda_m * Sm - mu_m * Sm
        dEm <- lambda_m * Sm - sigma * Em - mu_m * Em
        dIm <- sigma * Em - mu_m * Im
        
        list(c(dSh, dIh, dSm, dEm, dIm))
    })
}

# change it to data frame 
out_put_baseline<- as.data.frame(ode(
    y     = initial_state_values,
    times = times,
    func  = base_model,
    parms = parameters
))


# ---- Tidy + Plot ----
out_put_long_baseline <- out_put_baseline%>%
    pivot_longer(-time, 
                 names_to = "Compartment", 
                 values_to = "Count")


# --- Human plot ---

# Questions 
# a. Calculate the prevalence of the disease in the population at the end of
#the simulation period?
# Prevalence at day 365: 0.5918367 ( 59.18 % )
# -----------------------------
# Prevalence calculation
# -----------------------------

out_put_baseline$Nh <- out_put_baseline$Sh + out_put_baseline$Ih
out_put_baseline$prevalence <- out_put_baseline$Ih / out_put_baseline$Nh

# period prevalence 
final_prev <- tail(out_put_baseline$prevalence, 1)
cat("Prevalence at day 365:", final_prev, 
    "(", round(final_prev*100, 2), "% )\n")

# point prevalence 

final_prev_point <- (out_put_baseline$prevalence)

final_prev_point

# PLOT-----------------------------------------------
# Final point prevalence (day 365)
final_prev <- tail(out_put_baseline$prevalence, 1)

# Plot prevalence with dotted line

ggplot(out_put_baseline, aes(x = time, y = prevalence)) +
    geom_line(color = "blue", linewidth = 1.2) +
    geom_hline(yintercept = final_prev, linetype = "dotted", color = "red", size = 1) +
    labs(
        x = "Time (days)",
        y = "Prevalence (proportion infected)",
        title = "Prevalence of Infection Over Time"
    ) +
    annotate("text", x = max(out_put$time)*0.8, y = final_prev + 0.01, 
             label = paste0("Final point prevalence = ", round(final_prev*100, 2), "%"),
             color = "red", hjust = 0) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )


#b. What is your interpretation of the computed prevalence?
# At the end of the follow-up period 59.2% of the people has been infected 

#Baseline plot=============================================================
#Humans --------------------------------------------------------------

human_plot <- out_put_long_baseline %>%
    filter(Compartment %in% c("Sh", "Ih")) %>%
    plot_ly(x = ~time, y = ~Count, color = ~Compartment, 
            colors = c("green", "orange"),
            type = 'scatter', mode = 'lines') %>%
    layout(title = "Human Compartments",
           xaxis = list(title = "Days"),
           yaxis = list(title = "Number of human"))


# --- Mosquitoes ---
vector_plot <- out_put_long_baseline %>%
    filter(Compartment %in% c("Sm", "Im")) %>%
    plot_ly(x = ~time, y = ~Count, color = ~Compartment, 
            colors = c("green", "orange"),
            type = 'scatter', mode = 'lines') %>%
    layout(title = "Mosquito Compartments",
           xaxis = list(title = "Days"),
           yaxis = list(title = "Number of Mosquitoes"))


# Show them
human_plot
vector_plot

##################################################################################
# Vector control================================================================
##################################################################################

Nh<- 15000
Ih<- 25
Nm<-Nh*8
Im<-0.12*Nm

bitting<- 0.3
seeking<- 0.1
transmission<-0.45

alpha= bitting*seeking*transmission


parameters<-list(alpha= alpha,
                 beta= 0.5,
                 gamma= 1/25,
                 sigma= 1/7,
                 mu_m= 1/14*1.6, # the vector control has killed 60% of the mosquito
                 Nh= 15000,
                 Nm= 15000*8
)

initial_state_values <- c(
    Sh = (parameters$Nh - 25),
    Ih = 25,
    Sm = (parameters$Nh*8 - 0.12*parameters$Nh),
    Em = 0,
    Im = (0.12*parameters$Nh)
)


## Time and solve
times <- seq(0, 365, by =1)  # days

vector_model<- function(time, state, parameters){
    with(as.list(c(state, parameters)),{
        
        Nh <- Sh + Ih   # human population
        Nm <- Sm +Em + Im  # mosquito population
        
        # Forces of infection (Ross-Macdonald style)
        
        lambda_h <- alpha * (Im / Nh)       # per-susceptible-human per day
        lambda_m <- beta * (Ih / Nh)       # per-susceptible-mosquito per day
        
        dSh <- - lambda_h * Sh + gamma * Ih
        dIh <-   lambda_h * Sh - gamma * Ih
        
        dSm <- mu_m *Nm-lambda_m * Sm - mu_m * Sm
        dEm <- lambda_m * Sm - sigma * Em - mu_m * Em
        dIm <- sigma * Em - mu_m * Im
        
        list(c(dSh, dIh, dSm, dEm, dIm))
    })
}


out_put_vector <- as.data.frame(ode(
    y     = initial_state_values,
    times = times,
    func  = vector_model,
    parms = parameters
))


# ---- Tidy + Plot ----
out_put_long_vector <- out_put%>%
    pivot_longer(-time, 
                 names_to = "Compartment", 
                 values_to = "Count")



#Vector control plot=============================================================
#Humans --------------------------------------------------------------
Prevalence at day 365: 0.5085714 ( 50.86 % )
# Questions 
# a. Calculate the prevalence of the disease in the population at the end of
#the simulation period?

# -----------------------------
# Prevalence calculation
# -----------------------------
out_put_vector$Nh <- out_put_vector$Sh + out_put_vector$Ih
out_put_vector$prevalence <- out_put_vector$Ih / out_put_vector$Nh

final_prev <- tail(out_put_vector$prevalence, 1)
cat("Prevalence at day 365:", final_prev, 
    "(", round(final_prev*100, 2), "% )\n")

#b. What is your interpretation of the computed prevalence?

#the vector control has reduced the prevalence by  8.32% compared to without intervention.

# --- Humans ---
human_plot <- out_put_long_vector%>%
    filter(Compartment %in% c("Sh", "Ih")) %>%
    plot_ly(x = ~time, y = ~Count, color = ~Compartment, 
            colors = c("green", "orange"),
            type = 'scatter', mode = 'lines') %>%
    layout(title = "Human Compartments",
           xaxis = list(title = "Days"),
           yaxis = list(title = "Number of human"))


# --- Mosquitoes ---
vector_plot <- out_put_long_vector%>%
    filter(Compartment %in% c("Sm", "Im")) %>%
    plot_ly(x = ~time, y = ~Count, color = ~Compartment, 
            colors = c("green", "orange"),
            type = 'scatter', mode = 'lines') %>%
    layout(title = "Mosquito Compartments",
           xaxis = list(title = "Days"),
           yaxis = list(title = "Number of Mosquitoes"))


# Show them
human_plot
vector_plot


#plot for basline vs vector control----------------------- ---------------------

# Add a label column to distinguish scenarios
out_put_baseline <- out_put_baseline %>%
    mutate(scenario = "Baseline")

out_put_vector <- out_put_vector %>%
    mutate(scenario = "Vector Control")

# Combine both outputs
out_put_combined <- bind_rows(out_put_baseline, out_put_vector)

# Calculate final prevalence for each scenario
final_prev_baseline <- tail(out_put_baseline$prevalence, 1)
final_prev_vector <- tail(out_put_vector$prevalence, 1)

# Plot
ggplot(out_put_combined, aes(x = time, y = prevalence, color = scenario)) +
    geom_line(linewidth = 1.2) +
    
    # Add dotted horizontal lines for each scenario
    geom_hline(yintercept = final_prev_baseline, linetype = "dotted", color = "blue") +
    geom_hline(yintercept = final_prev_vector, linetype = "dotted", color = "red") +
    
    # Labels
    labs(
        x = "Time (days)",
        y = "Prevalence (proportion infected)",
        title = "Comparison of Infection Prevalence: Baseline vs Vector Control"
    ) +
    
    # Annotate final prevalence values
    annotate("text", x = max(out_put_combined$time)*0.8, y = final_prev_baseline + 0.01, 
             label = paste0("Baseline final = ", round(final_prev_baseline*100, 2), "%"),
             color = "blue", hjust = 0) +
    annotate("text", x = max(out_put_combined$time)*0.8, y = final_prev_vector + 0.03, 
             label = paste0("Vector control final = ", round(final_prev_vector*100, 2), "%"),
             color = "red", hjust = 0) +
    
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

#==============================================================================
# Vaccination------------------------------------------------------------------
#==============================================================================

# Flow diagram________________________


graph <- grViz("
digraph model {
  graph [layout = dot, rankdir = LR]

  # Human compartments
  Sh [label = 'Sh\\nSusceptible Humans', shape = box, style=filled, fillcolor=lightblue]
  Ih [label = 'Ih\\nInfected Humans', shape = box, style=filled, fillcolor=lightpink]
  Vh [label = 'Vh\\nVaccinated Humans', shape = box, style=filled, fillcolor=lightgoldenrod]

  # Mosquito compartments
  Sm [label = 'Sm\\nSusceptible Mosquitoes', shape = box, style=filled, fillcolor=lightgreen]
  Em [label = 'Em\\nExposed Mosquitoes', shape = box, style=filled, fillcolor=khaki]
  Im [label = 'Im\\nInfected Mosquitoes', shape = box, style=filled, fillcolor=orange]

  # Birth/Death node
  BD [label = 'Birth / Death', shape = ellipse, style=filled, fillcolor=grey90]

  # Human transitions
  Sh -> Ih [label = 'λh * Sh']
  Ih -> Sh [label = 'γ * Ih']
  Vh -> Sh [label = 'waning * Vh']
  Vh -> Ih [label = 'λh * (1 - efficacy) * Vh']

  # Mosquito transitions
  Sm -> Em [label = 'λm * Sm']
  Em -> Im [label = 'σ * Em']
  Sm -> BD [label = 'μm * Sm']
  Em -> BD [label = 'μm * Em']
  Im -> BD [label = 'μm * Im']
  BD -> Sm [label = 'birth']
}
")

graph

# -----------------------------------------------------------------------------
# Parameters--------------------------------------------------------------------
# -----------------------------------------------------------------------------
Nh <- 15000
Ih <- 25
Nm <- Nh * 8


Im <- 0.12 * Nm       # 12% of mosquitoes infected initially
Sm <- Nm - Im
Em <- 0

efficacy <- 0.65
waning   <- 1/30      # waning immunity rate (per day)
coverage <- 0.30
vacc_cost <- 20       # USD per vaccinated person
bitting<- 0.3
seeking<- 0.1
transmission<-0.45

Sh0 <- Nh - Ih
Vh  <- coverage * Sh0        # 30% vaccination at t=0
Sh  <- Sh0 - Vh

alpha <- bitting * seeking * transmission

# -----------------------------
# Initial states
# -----------------------------
initial_state_values <- c(
    Sh = Sh,
    Ih = Ih,
    Vh = Vh,
    Sm = Sm,
    Em = Em,
    Im = Im
)

# -----------------------------
# Parameters
# -----------------------------
parameters <- list(
    alpha    = alpha,
    beta     = 0.5,      # human → mosquito transmission
    gamma    = 1/25,     # human recovery
    sigma    = 1/7,      # mosquito incubation
    mu_m     = 1/14,     # mosquito death
    nu       = 0.3,      # vaccination rate
    efficacy = efficacy,
    waning   = waning
)

# Vaccination cost
total_cost <- Vh * vacc_cost
cat("Total vaccination cost (USD):", total_cost, "\n")

# -----------------------------
# Model function
# -----------------------------
vaccination_model <- function(time, state, parameters){
    with(as.list(c(state, parameters)), {
        
        Nh <- Sh + Ih + Vh
        Nm <- Sm + Em + Im
        
        # Forces of infection
        lambda_h <- alpha * (Im / Nh)   # mosquito → human
        lambda_m <- beta  * (Ih / Nh)   # human → mosquito
        
        # Humans
        
        dSh <- -lambda_h * Sh + gamma * Ih + waning * Vh
        dIh <-  lambda_h * Sh + (1 - efficacy) * lambda_h * Vh - gamma * Ih
        dVh <- - (1 - efficacy) * lambda_h * Vh - waning * Vh
        
        # Mosquitoes
        dSm <- -lambda_m * Sm - mu_m * Sm + mu_m * Nm
        dEm <- lambda_m * Sm - sigma * Em - mu_m * Em
        dIm <- sigma * Em - mu_m * Im
        
        list(c(dSh, dIh, dVh, dSm, dEm, dIm))
    })
}

# -----------------------------
# Solve ODE
# -----------------------------
times <- seq(0, 365, by = 1)

out_put_vaccine<- as.data.frame(ode(
    y = initial_state_values,
    times = times,
    func = vaccination_model,
    parms = parameters
))


# ---- Tidy + Plot ----
out_put_long_vaccine<- out_put_vaccine%>%
    pivot_longer(-time, 
                 names_to = "Compartment", 
                 values_to = "Count")
# -----------------------------
# Prevalence calculation
# -----------------------------
out_put_vaccine$Nh <- out_put_vaccine$Sh + out_put_vaccine$Ih + 
    out_put_vaccine$Vh

out_put_vaccine$prevalence <- out_put_vaccine$Ih / out_put_vaccine$Nh

final_prev <- tail(out_put_vaccine$prevalence, 1)
cat("Prevalence at day 365:", final_prev,
    "(", round(final_prev*100, 2), "% )\n")

#Prevalence at day 365: 0.5918367 ( 59.18 % )

# ------------------------------------------------------------------------------
# Plot results
# ------------------------------------------------------------------------------

vaccine_plot <- out_put_long_vaccine %>%
    filter(Compartment %in% c("Sm", "Im", "Vh")) %>%
    ggplot(aes(x = time, y = Count, color = Compartment)) +
    geom_line(size = 1) +
    scale_color_manual(
        values = c("Sm" = "green", "Im" = "orange", "Vh" = "red"),
        labels = c("Sm" = "Susceptible",
                   "Im" = "Infectious", 
                   "Vh" = "Vaccinated")
    ) +
    labs(
        x = "Days",
        y = "Population",
        color = "Compartments",
        title = "Human compartments with vaccination"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right"
    )

vaccine_plot


#vaccine vs vector vs baseline----------------------------------

# Make sure each dataset has prevalence already computed:
out_put_baseline$prevalence <- out_put_baseline$Ih / out_put_baseline$Nh
out_put_vector$prevalence   <- out_put_vector$Ih / out_put_vector$Nh
out_put_vaccine$prevalence  <- out_put_vaccine$Ih / out_put_vaccine$Nh

# Add a column to identify the scenario
out_put_baseline$scenario <- "Baseline"
out_put_vector$scenario   <- "Vector Control"
out_put_vaccine$scenario  <- "Vaccine"

# Combine into one dataframe
all_out <- bind_rows(out_put_baseline, out_put_vector, out_put_vaccine)

# Calculate final prevalence per scenario
final_prevs <- all_out %>%
    group_by(scenario) %>%
    summarise(final_prev = tail(prevalence, 1), .groups = "drop")

# Plot
ggplot(all_out, aes(x = time, y = prevalence, color = scenario)) +
    geom_line(linewidth = 1.2) +
    # Add final prevalence dotted lines for each scenario
    geom_hline(data = final_prevs,
               aes(yintercept = final_prev, color = scenario),
               linetype = "dotted", size = 1, show.legend = FALSE) +
    # Annotate final prevalence
    geom_text(data = final_prevs,
              aes(x = max(all_out$time)*0.9, 
                  y = final_prev + 0.01,
                  label = paste0("Final = ", round(final_prev*100, 2), "%"),
                  color = scenario),
              hjust = 0) +
    labs(
        x = "Time (days)",
        y = "Prevalence (proportion infected)",
        title = "Prevalence of Infection Over Time (Baseline vs Vector Control vs Vaccine)"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
#===================================================================
##############################################################
#Combined vector control + vaccination______________________
# -----------------------------


Nh <- 15000
Ih <- 25
Nm <- Nh * 8


Im <- 0.12 * Nm       # 12% of mosquitoes infected initially
Sm <- Nm - Im
Em <- 0

efficacy <- 0.65
waning   <- 1/30      # waning immunity rate (per day)
coverage <- 0.30
vacc_cost <- 20       # USD per vaccinated person
bitting<- 0.3
seeking<- 0.1
transmission<-0.45

Sh0 <- Nh - Ih
Vh  <- coverage * Sh0        # 30% vaccination at t=0
Sh  <- Sh0 - Vh

alpha <- bitting * seeking * transmission

# -----------------------------
# Initial states
# -----------------------------
initial_state_values <- c(
    Sh = Sh,
    Ih = Ih,
    Vh = Vh,
    Sm = Sm,
    Em = Em,
    Im = Im
)

# -----------------------------
# Parameters
# -----------------------------
parameters <- list(
    alpha    = alpha,
    beta     = 0.5,      # human → mosquito transmission
    gamma    = 1/25,     # human recovery
    sigma    = 1/7,      # mosquito incubation
    mu_m     = 1/14*1.6,     # mosquito death
    nu       = 0.3,      # vaccination rate
    efficacy = efficacy,
    waning   = waning
)

# Vaccination cost
total_cost <- Vh * vacc_cost
cat("Total vaccination cost (USD):", total_cost, "\n")

# -----------------------------
# Model function
# -----------------------------
vaccination_model <- function(time, state, parameters){
    with(as.list(c(state, parameters)), {
        
        Nh <- Sh + Ih + Vh
        Nm <- Sm + Em + Im
        
        # Forces of infection
        lambda_h <- alpha * (Im / Nh)   # mosquito → human
        lambda_m <- beta  * (Ih / Nh)   # human → mosquito
        
        # Humans
        
        dSh <- -lambda_h * Sh + gamma * Ih + waning * Vh
        dIh <-  lambda_h * Sh + (1 - efficacy) * lambda_h * Vh - gamma * Ih
        dVh <- - (1 - efficacy) * lambda_h * Vh - waning * Vh
        
        # Mosquitoes
        dSm <- -lambda_m * Sm - mu_m * Sm + mu_m * Nm
        dEm <- lambda_m * Sm - sigma * Em - mu_m * Em
        dIm <- sigma * Em - mu_m * Im
        
        list(c(dSh, dIh, dVh, dSm, dEm, dIm))
    })
}

# Run simulation
out_put_combined  <- ode(y = initial_state_values, 
                         times = times, 
                         func = vaccination_model, 
                         parms = parameters)

out_put_combined <- as.data.frame(out_put_combined)

# ---- Tidy + Plot ----
out_put_long_combined <- out_put_combined%>%
    pivot_longer(-time, 
                 names_to = "Compartment", 
                 values_to = "Count")

# Prevalence at end of simulation
out_put_combined$prevalence <- out_put_combined$Ih / (out_put_combined$Sh + out_put_combined$Ih + out_put_combined$Vh)

final_prev <- tail(out_put_combined$prevalence, 1)
cat("Prevalence at day 365:", round(final_prev*100,2), "%\n")



# Combined plot------------------------------------------------------------------
# Make sure each dataset has prevalence already computed:

out_put_baseline$prevalence <- out_put_baseline$Ih / out_put_baseline$Nh
out_put_vector$prevalence   <- out_put_vector$Ih / out_put_vector$Nh
out_put_vaccine$prevalence  <- out_put_vaccine$Ih / out_put_vaccine$Nh
out_put_combined$prevalence <- out_put_combined$Ih / (out_put_combined$Sh +
                                                          out_put_combined$Ih + 
                                                          out_put_combined$Vh)

# Add a column to identify the scenario
out_put_baseline$scenario <- "Baseline"
out_put_vector$scenario   <- "Vector Control"
out_put_vaccine$scenario  <- "Vaccine"
out_put_combined$scenario <- "Combined"

# Combine into one dataframe
all_out <- bind_rows(out_put_baseline, out_put_vector, out_put_vaccine, out_put_combined)

# Calculate final prevalence per scenario
final_prevs <- all_out %>%
    group_by(scenario) %>%
    summarise(final_prev = tail(prevalence, 1), .groups = "drop")

# Plot
ggplot(all_out, aes(x = time, y = prevalence, color = scenario)) +
    geom_line(linewidth = 1.2) +
    # Add final prevalence horizontal dotted lines
    geom_hline(data = final_prevs,
               aes(yintercept = final_prev, color = scenario),
               linetype = "dotted", size = 1, show.legend = FALSE) +
    # Annotate final prevalence
    geom_text(data = final_prevs,
              aes(x = max(all_out$time)*0.9, 
                  y = final_prev + 0.01,
                  label = paste0("Final = ", round(final_prev*100, 2), "%"),
                  color = scenario),
              hjust = 0) +
    labs(
        x = "Time (days)",
        y = "Prevalence (proportion infected)",
        title = "Prevalence of Infection Over Time (Baseline vs Vector Control vs Vaccine vs Combined)"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )


#===================================================================
#Model Assumptions
#===================================================================
##1.	Closed human and mosquito populations; no births/deaths except mosquito natural death.
## 2.	Homogeneous mixing: each mosquito can bite any human equally.
## 3.	Vaccine is leaky: reduces infection probability by 65%.
## 4.	Vaccine immunity wanes after ~30 days (return to susceptibility).
## 5.	Human recovery is immediate after infectious period (~25 days), returning to susceptible.
## 6.	Mosquitoes have a latent period of 7 days; after that they remain infectious until death (~14 days).
## 7.	No seasonal variation in transmission other than the initial rainy season assumption.
## 8.	Population sizes: humans = 15,000, mosquitoes = 8 × 15,000 = 120,000 initially.

#==============================================================================

# residual spraying 
# Interpretation
# final_prev gives the proportion of infectious humans at the end of 1 year.
# For a vaccination with 30% coverage and 65% efficacy, the prevalence should be lower than baseline (without vaccination), but not zero due to:
## partial coverage
## leaky vaccine efficacy
## waning immunity
## This illustrates the impact of vaccination on reducing transmission, while showing that combining with vector control could further lower prevalence.

##__________________________________________________________________________________
# Narrative description_______________________________________________________

# Populations (t = 0)
##	Humans: N_h(0)=150000
##	Mosquitoes per human: 8*Nm
##	Infectious humans: I_h(0)=25
##	Infectious mosquitoes:
##	Vaccination at baseline: 30% of susceptibles 


#Epidemiology and life history
##	Human infectious period: 25 days 
##	Mosquito lifespan: 14 days 
##	Mosquito latent period (E→I): 7 days 
##	Biting rate: a = 0.3 bites/mosquito/day
##	Mosquito→human transmission prob per bite: alpha
##	Human→mosquito transmission prob per bite: beta

#Vaccination
## Coverage at baseline: 30% of susceptibles (encoded in initial state)
## Vaccine efficacy (leaky): varepsilon = 0.65
## Waning: average 30 days :omega = 1/30


#Vector control
##	Increases mosquito mortality by 60% 
##	(No change to biting or recruitment unless you want to add them)

#Costs
##	Vaccination cost: $20 per person vaccinated at t=0
##	Vector control annual budget: $80,000

#===================================================================
# Infections averted with interevetion_____________________________
#===================================================================

# -----------------------------------------------------------------------------
# Initial state values
# ------------------------------------------------------------------------------
Nh <- 15000

initial_state_values <- c(
    Sh = Nh - 25,
    Ih = 25,
    Vh = Nh * 0.30,
    Sm = Nh * 8 - 0.12 * Nh,
    Em = 0,
    Im = 0.12 * Nh,
    CumInf = 0   # cumulative infections
)

# ------------------------------------------------------------------------------
# Malaria transmission model----------------------------------------------------
# ------------------------------------------------------------------------------

malaria_model <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
        Nh <- Sh + Ih + Vh
        Nm <- Sm + Em + Im
        
        mu_m_eff <- (1 + vc_mort_increase) * mu_m
        alpha <- bites * seeking * transmission
        
        lambda_h <- alpha * (Im / Nh)
        lambda_m <- beta * (Ih / Nh)
        
        # New infections in humans (incidence flow)
        new_inf <- lambda_h * Sh + (1 - epsilon) * lambda_h * Vh
        
        # Humans
        dSh <- -lambda_h * Sh + gamma * Ih + omega * Vh
        dIh <- new_inf - gamma * Ih
        dVh <- -(1 - epsilon) * lambda_h * Vh - omega * Vh
        
        # Mosquitoes
        dSm <- -lambda_m * Sm - mu_m_eff * Sm
        dEm <-  lambda_m * Sm - sigma * Em - mu_m_eff * Em
        dIm <-  sigma * Em - mu_m_eff * Im
        
        # Cumulative infections
        dCumInf <- new_inf
        
        list(c(dSh, dIh, dVh, dSm, dEm, dIm, dCumInf))
    })
}

# ------------------------------------------------------------------------------
# Run one scenario
# ------------------------------------------------------------------------------
run_scenario <- function(params, times) {
    ode(y = initial_state_values, times = times,
        func = malaria_model, parms = params) %>%
        as.data.frame()
}

# ------------------------------------------------------------------------------
# Base parameters
# ------------------------------------------------------------------------------

base_params <- c(
    bites  = 0.3,
    seeking   = 0.1,
    transmission = 0.45,
    beta  = 0.5,
    gamma  = 1/25,
    mu_m   = 1/14,
    sigma  = 1/7,
    epsilon= 0.0,       # vaccine efficacy
    omega  = 0.0,       # waning
    nu     = 0.0,       # not used
    vc_mort_increase = 0.0 # vector control mortality increase
)

times <- seq(0, 365, by = 1)

#------------------------------------------------------------------------------
# Define scenarios
#------------------------------------------------------------------------------
scenarios <- list(
    baseline = base_params,
    vaccine  = {p <- base_params; p["epsilon"] <- 0.65; p},
    vector   = {p <- base_params; p["vc_mort_increase"] <- 0.60; p},
    combined = {p <- base_params; p["epsilon"] <- 0.65; p["vc_mort_increase"] <- 0.60; p}
)

#------------------------------------------------------------------------------
# Process results
#------------------------------------------------------------------------------

process_result <- function(df) {
    df$Nh <- df$Sh + df$Ih + df$Vh
    
    # daily incidence (new cases each day)
    df$daily_incidence <- c(NA, diff(df$CumInf))
    
    # total incident cases over the simulation period
    total_incident <- tail(df$CumInf, 1) - head(df$CumInf, 1)
    
    # mean population (for denominator)
    mean_pop <- mean(df$Nh, na.rm = TRUE)
    
    # incidence per 1000 persons over the year
    incidence_per_1000 <- (total_incident / mean_pop) * 1000
    
    list(
        df = df,
        total_incident = total_incident,
        incidence_per_1000 = incidence_per_1000,
        mean_pop = mean_pop
    )
}

#------------------------------------------------------------------------------
# Run all scenarios
#------------------------------------------------------------------------------
results <- lapply(scenarios, run_scenario, times = times)
processed <- lapply(results, process_result)

# Extract totals
totals <- sapply(processed, function(x) x$total_incident)
inc_per_1000 <- sapply(processed, function(x) x$incidence_per_1000)
mean_pops <- sapply(processed, function(x) x$mean_pop)

# Cases averted relative to baseline
averted <- totals["baseline"] - totals
percent_averted <- (averted / totals["baseline"]) * 100

#------------------------------------------------------------------------------
# Output summary table
#------------------------------------------------------------------------------

output_table <- data.frame(
    Scenario = names(totals),
    Total_Incident = round(totals, 0),
    Incidence_per_1000 = round(inc_per_1000, 1),
    Mean_Pop = round(mean_pops, 1),
    Cases_Averted = round(averted, 0),
    Percent_Averted = round(percent_averted, 1),
    row.names = NULL
)

print(output_table)

#------------------------------------------------------------------------------
# Plot incidence curves
#------------------------------------------------------------------------------

plot_df <- bind_rows(
    lapply(names(processed), function(name) {
        df <- processed[[name]]$df
        df$Scenario <- name
        df
    })
)

# Daily incidence plot
ggplot(plot_df, aes(x = time, y = daily_incidence, color = Scenario)) +
    geom_line(size = 1) +
    labs(title = "Daily Incidence of Malaria",
         x = "Time (days)", y = "New cases per day") +
    theme_minimal()

# Cumulative incidence plot
ggplot(plot_df, aes(x = time, y = CumInf, color = Scenario)) +
    geom_line(size = 1) +
    labs(title = "Cumulative Incidence of Malaria",
         x = "Time (days)", y = "Cumulative cases") +
    theme_minimal()

#======================================================
#cost---------------------------------------------------------------------------
#======================================================

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------

Nh <- 15000         # human population
Ih0 <- 25           # initial infectious humans
Sh0 <- Nh - Ih0
Vh0 <- 0            # vaccinated initially

Nm <- Nh * 8        # mosquito population
Im0 <- round(0.12 * Nm)   # infectious mosquitoes
Em0 <- 0            # latent mosquitoes
Sm0 <- Nm - Im0

# Epidemiological parameters
bitting<- 0.3
seeking<- 0.1
transmission<-0.45
gamma <- 1/25       # human recovery (1/25 days)
mu_m <- 1/14        # mosquito mortality (1/14 days)
b <- 0.3            # mosquito biting rate
p_h <- bitting*seeking*transmission # prob mosquito→human transmission
p_m <- 0.5          # prob human→mosquito transmission
sigma_m <- 1/7      # mosquito latent rate


# Vaccination
vacc_cov <- 0.30    # coverage
vacc_eff <- 0.65    # efficacy
waning <- 1/30      # waning immunity (30 days)
vacc_cost <- 20     # cost per vaccinated person

# Vector control
vc_eff <- 0.60      # 60% increase in mortality
vc_budget <- 80000

#------------------------------------------------------------------------------
# Model equations
#------------------------------------------------------------------------------
model <- function(t, state, parms) {
    with(as.list(c(state, parms)), {
        
        Nh <- Sh + Ih + Vh
        Nm <- Sm + Em + Im
        
        # Forces of infection
        lambda_h <- b * p_h * (Im / Nh)
        lambda_m <- b * p_m * (Ih / Nh)
        
        # Humans
        dSh <- -lambda_h * Sh + gamma * Ih + waning * Vh
        dIh <- lambda_h * Sh + lambda_h * (1 - vacc_eff) * Vh - gamma * Ih
        dVh <- -lambda_h * (1 - vacc_eff) * Vh - waning * Vh
        
        # Mosquitoes
        dSm <- -lambda_m * Sm - mu_m * Sm
        dEm <- lambda_m * Sm - sigma_m * Em - mu_m * Em
        dIm <- sigma_m * Em - mu_m * Im
        
        # Track incidence
        incidence <- lambda_h * Sh + lambda_h * (1 - vacc_eff) * Vh
        
        return(list(c(dSh, dIh, dVh, dSm, dEm, dIm),
                    incidence = incidence))
    })
}

#------------------------------------------------------------------------------
# Simulation wrapper
#------------------------------------------------------------------------------
simulate_scenario <- function(scenario) {
    
    # Adjust parameters
    parms <- list()
    parms$mu_m <- mu_m
    parms$vacc_eff <- ifelse(scenario %in% c("vacc", "both"), vacc_eff, 0)
    parms$waning <- ifelse(scenario %in% c("vacc", "both"), waning, 0)
    parms$gamma <- gamma
    parms$b <- b; parms$p_h <- p_h; parms$p_m <- p_m; parms$sigma_m <- sigma_m
    
    # Vector control
    parms$mu_m <- ifelse(scenario %in% c("vc", "both"), mu_m * (1 + vc_eff), mu_m)
    
    # Initial states
    Sh <- Sh0
    Ih <- Ih0
    Vh <- 0
    if (scenario %in% c("vacc", "both")) {
        Sh <- round((1 - vacc_cov) * Sh0)
        Vh <- round(vacc_cov * Sh0)
    }
    
    state <- c(Sh = Sh, Ih = Ih, Vh = Vh, Sm = Sm0, Em = Em0, Im = Im0)
    
    times <- seq(0, 365, by = 1)
    
    out <- ode(y = state, times = times, func = model, parms = parms)
    out <- as.data.frame(out)
    
    out$cum_incidence <- cumsum(out$incidence)
    return(out)
}

#------------------------------------------------------------------------------
# Run scenarios
#------------------------------------------------------------------------------
baseline <- simulate_scenario("base")
vc <- simulate_scenario("vc")
vacc <- simulate_scenario("vacc")
both <- simulate_scenario("both")

#------------------------------------------------------------------------------
# Cost-effectiveness
#------------------------------------------------------------------------------
cases_base <- tail(baseline$cum_incidence, 1)
cases_vc <- tail(vc$cum_incidence, 1)
cases_vacc <- tail(vacc$cum_incidence, 1)
cases_both <- tail(both$cum_incidence, 1)

# Cases averted
ca_vc <- cases_base - cases_vc
ca_vacc <- cases_base - cases_vacc
ca_both <- cases_base - cases_both

# Costs
cost_vc <- vc_budget
cost_vacc <- vacc_cov * Sh0 * vacc_cost
cost_both <- cost_vc + cost_vacc

# Cost per case averted
cea <- tibble(
    Scenario = c("Vector Control", "Vaccination", "Combined"),
    Cases_Averted = c(ca_vc, ca_vacc, ca_both),
    Cost_USD = c(cost_vc, cost_vacc, cost_both),
    Cost_per_Case_Averted = Cost_USD / Cases_Averted
)

print(cea)

#------------------------------------------------------------------------------
# Plot incidence curves
#------------------------------------------------------------------------------

incidence_plot <- baseline %>%
    select(time, base = incidence) %>%
    left_join(vc %>% select(time, vc = incidence), by = "time") %>%
    left_join(vacc %>% select(time, vacc = incidence), by = "time") %>%
    left_join(both %>% select(time, both = incidence), by = "time") %>%
    pivot_longer(-time, names_to = "Scenario", values_to = "Incidence")

ggplot(incidence_plot, aes(x = time, y = Incidence, color = Scenario)) +
    geom_line(size = 1) +
    labs(title = "Incidence over 365 days", x = "Time (days)", y = "New infections per day")

# 

#------------------------------------------------------------------------------
#Cost assumptions
#------------------------------------------------------------------------------
cost_vaccine <- 20 * (0.30 * Sh)   # 20 USD per vaccinated person
cost_vector  <- 80000              # assumed cost for vector control
cost_combined <- cost_vaccine + cost_vector

costs <- c(
    baseline = 0,
    vaccine = cost_vaccine,
    vector = cost_vector,
    combined = cost_combined
)

#------------------------------------------------------------------------------
# Cost-effectiveness results
#------------------------------------------------------------------------------
cea <- data.frame(
    Scenario = names(totals),
    Total_Incident = round(totals, 0),
    Cases_Averted = round(averted, 0),
    Cost_USD = round(costs[names(totals)], 0),
    Cost_per_Case_Averted = round(costs[names(totals)] / pmax(averted, 1), 2) # avoid div by zero
)

print(cea)

#------------------------------------------------------------------------------
# Case averted vs cost-effectiveness plot
#------------------------------------------------------------------------------

cea_plot <- cea %>%
    filter(Scenario != "baseline")  # exclude baseline (0 cost, 0 averted)

ggplot(cea_plot, aes(x = Cases_Averted, y = Cost_per_Case_Averted, color = Scenario)) +
    geom_point(size = 6, alpha = 0.8) +
    geom_text(aes(label = paste0("US$", Cost_per_Case_Averted)),
              vjust = -1, size = 4.2, show.legend = FALSE) +
    geom_line(aes(group = 1), color = "grey40", linetype = "dashed") +
    labs(
        title = "Cost-Effectiveness of Malaria Interventions",
        subtitle = "Costs vs. cases averted (365 days)",
        x = "Cases Averted",
        y = "Cost per case averted"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")

