###### Malaria Network Transmission Model ####
##Network Analysis and Malaria Epidemic Model##

####Project Description####
##This project models malaria transmission through social networks incorporating household clustering, 
##geographic proximity, seasonal patterns, and vector-mediated transmission. 
##Various intervention strategies (ITNs, IRS, case management, mass screening) are evaluated.

#Packages for this work
library(igraph) # For graphs
library(ggplot2) # For plots
library(dplyr) # For manipulation
library(RColorBrewer) # For colours

##Create enhanced dataset with household and age information
data <- data.frame(
  Participant = c("KAM", "SIN", "RAH", "DAW", "LUC", "TUH", "BAS", "MBA", "MUS", "DAV", "CHA", "EMM", "LYD", "VAL", "GOD", "AGN", "ROS", "FRA", "JAM", "DGB", "ELI", "WAN", "MUT", "GEO", "GEF", "SAH", "PAT", "TAB", "BRI", "STA", "OUM", "STE", "THU", "EMI", "CAM", "ZEN", "BIL", "JOH", "ALI", "PUR", "WAI", "MIL", "BRL", "TRI"),  Country = c("Cameroon", "Ethiopia", "Tanzania", "Ethiopia", "Zambia", "Uganda", "Nigeria", "Cameroon", "Zimbabwe", "Zambia", "Kenya", "Kenya", "Kenya", "Kenya", "Malawi", "Malawi", "Tanzania", "Tanzania", "Kenya", "Ethiopia", "Tanzania", "Kenya", "Kenya", "Kenya", "Kenya", "Tunisia", "Kenya", "Kenya", "Kenya", "Kenya", "Tunisia", "Kenya", "Kenya", "Switzerland", "Switzerland", "Ghana", "Switzerland", "Kenya", "Kenya", "Kenya", "Kenya", "Kenya", "Kenya", "Kenya"),
  # Add household information (grouped by accommodation/country)
  Household = c(1,2,3,2,4,5,6,1,7,4,8,8,9,9,10,10,3,3,11,2,3,12,13,14,15,16,17,18,19,20,
                16,21,22,23,23,24,23,25,26,27,28,29,30,31),
  # Add age groups: Child (5-17), Adult (18-64), Elder (65+)
  Age_Group = c("Adult","Adult","Adult","Child","Adult","Adult","Adult","Child","Adult","Adult","Adult","Child","Adult","Adult","Adult","Child","Adult","Adult","Adult","Adult", "Child","Adult","Adult","Adult","Adult","Adult","Adult","Child","Adult","Adult","Adult","Adult","Child","Adult","Adult","Adult","Adult","Adult","Child","Adult", "Adult","Adult","Adult","Adult"),
  # Add village clustering (based on country but with sub-villages)
  Village = c(1,2,3,2,4,5,6,1,7,4,8,8,8,8,9,9,3,3,8,2,3,8,8,8,8,10,8,8,8,8,
              10,8,8,11,11,12,11,8,8,8,8,8,8,8)
)

#####PART 1. MALARIA NETWORK CONSTRUCTION ####

##Create malaria-specific network with household clustering and geographic proximity
create_malaria_network <- function(data, household_prob = 0.7, village_prob = 0.3, age_mixing = TRUE, seasonal_factor = 1.0) {
  n <- nrow(data)
  adj_matrix <- matrix(0, n, n)
  rownames(adj_matrix) <- colnames(adj_matrix) <- data$Participant
  
  set.seed(123)
  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j) {
        base_prob <- 0.05  # Base connection probability
        
        # Household clustering - much higher within households
        if(data$Household[i] == data$Household[j]) {
          base_prob <- household_prob
        }
        # Village proximity - moderate increase for same village
        else if(data$Village[i] == data$Village[j]) {
          base_prob <- village_prob
        }
        
        # Age-structured mixing (children connect more with children, adults with adults)
        if(age_mixing) {
          if(data$Age_Group[i] == data$Age_Group[j]) {
            base_prob <- base_prob * 1.5
          }
          # Children connect more overall (school, play)
          if(data$Age_Group[i] == "Child" | data$Age_Group[j] == "Child") {
            base_prob <- base_prob * 1.2
          }
        }
        
        # Seasonal mobility factor (farming, markets increase connections)
        base_prob <- base_prob * seasonal_factor
        
        if(runif(1) < base_prob) {
          adj_matrix[i,j] <- 1
        }
      }
    }
  }
  return(adj_matrix)
}

##Create the malaria network with mosquito breeding sites as hubs
adj_matrix <- create_malaria_network(data)
g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

##Add attributes
V(g)$country <- data$Country
V(g)$name <- data$Participant
V(g)$household <- data$Household
V(g)$age_group <- data$Age_Group
V(g)$village <- data$Village

##Enhanced network visualization
plot_malaria_network <- function(g, data) {
  # Color nodes by age group
  age_colors <- c("Child" = "red", "Adult" = "blue", "Elder" = "green")
  V(g)$color <- age_colors[V(g)$age_group]
  
  # Size nodes by household size (larger households = larger nodes)
  household_sizes <- table(data$Household)
  V(g)$size <- household_sizes[as.character(V(g)$household)] * 2
  
  set.seed(123)
  layout <- layout_with_fr(g)
  
  plot(g,
       layout = layout,
       vertex.label = V(g)$name,
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       vertex.frame.color = "white",
       edge.color = "gray70",
       edge.width = 0.8,
       main = "Malaria Transmission Network\n(Node size = Household size, Color = Age group)")
  
  legend("topright", legend = names(age_colors), fill = age_colors, 
         cex = 0.8, title = "Age Group")
}

plot_malaria_network(g, data)

####Part 2 - Malaria-Specific SEIR Model ####
malaria_seir_model <- function(adj_matrix, beta, sigma, gamma, immunity_waning, 
                               initial_infected, days, interventions = list(), 
                               environmental_factor = 1.0, asymptomatic_rate = 0.3) {
  n <- nrow(adj_matrix)
  participants <- rownames(adj_matrix)
  
  # States: 0=S, 1=E, 2=I_symptomatic, 3=I_asymptomatic, 4=R, 5=Immune_waning
  states <- array(0, dim = c(days + 1, n, 6))
  dimnames(states) <- list(NULL, participants, c("S", "E", "I_symp", "I_asymp", "R", "Immune"))
  
  # Initialize - everyone starts susceptible
  states[1, , "S"] <- 1
  infected_idx <- which(participants %in% initial_infected)
  states[1, infected_idx, "S"] <- 0
  states[1, infected_idx, "I_symp"] <- 1
  
  # Track interventions per individual
  itn_protection <- rep(1, n)  # 1 = no protection, <1 = protected
  case_management <- rep(FALSE, n)
  
  for(day in 1:days) {
    # Environmental effects (rainfall, temperature - simplified as seasonal factor)
    current_beta <- beta * environmental_factor * (1 + 0.3 * sin(2 * pi * day / 365))
    
    # Current states
    S_current <- states[day, , "S"]
    E_current <- states[day, , "E"]
    I_symp_current <- states[day, , "I_symp"]
    I_asymp_current <- states[day, , "I_asymp"]
    R_current <- states[day, , "R"]
    Immune_current <- states[day, , "Immune"]
    
    # Next day states
    S_next <- S_current
    E_next <- E_current
    I_symp_next <- I_symp_current
    I_asymp_next <- I_asymp_current
    R_next <- R_current
    Immune_next <- Immune_current
    
    # Apply interventions
    if("ITN" %in% names(interventions) && day >= interventions$ITN$start_day) {
      itn_coverage <- interventions$ITN$coverage
      itn_efficacy <- interventions$ITN$efficacy
      # Randomly assign ITNs based on coverage
      if(day == interventions$ITN$start_day) {
        itn_users <- sample(1:n, size = round(n * itn_coverage))
        itn_protection[itn_users] <- 1 - itn_efficacy
      }
    }
    
    if("IRS" %in% names(interventions) && day >= interventions$IRS$start_day) {
      # IRS affects everyone in treated villages
      current_beta <- current_beta * (1 - interventions$IRS$efficacy * interventions$IRS$coverage)
    }
    
    # SEIR transitions
    for(i in 1:n) {
      # S -> E (Susceptible to Exposed)
      if(S_current[i] == 1) {
        # Force of infection from both symptomatic and asymptomatic
        infectious_neighbors <- which((I_symp_current == 1 | I_asymp_current == 1) & adj_matrix[i, ] == 1)
        if(length(infectious_neighbors) > 0) {
          total_force <- length(infectious_neighbors) * current_beta * itn_protection[i]
          prob_infection <- 1 - exp(-total_force)
          
          if(runif(1) < prob_infection) {
            S_next[i] <- 0
            E_next[i] <- 1
          }
        }
      }
      
      # E -> I (Exposed to Infected)
      else if(E_current[i] == 1) {
        if(runif(1) < sigma) {
          E_next[i] <- 0
          # Decide if symptomatic or asymptomatic
          if(runif(1) < asymptomatic_rate) {
            I_asymp_next[i] <- 1
          } else {
            I_symp_next[i] <- 1
          }
        }
      }
      
      # I -> R (Infected to Recovered)
      else if(I_symp_current[i] == 1 || I_asymp_current[i] == 1) {
        recovery_rate <- gamma
        
        # Case management increases recovery rate for symptomatic cases
        if("case_management" %in% names(interventions) && I_symp_current[i] == 1) {
          if(day >= interventions$case_management$start_day) {
            if(runif(1) < interventions$case_management$coverage) {
              recovery_rate <- recovery_rate * 3  # Faster recovery with treatment
            }
          }
        }
        
        if(runif(1) < recovery_rate) {
          I_symp_next[i] <- 0
          I_asymp_next[i] <- 0
          R_next[i] <- 1
        }
      }
      
      # R -> S (Immunity waning)
      else if(R_current[i] == 1) {
        if(runif(1) < immunity_waning) {
          R_next[i] <- 0
          S_next[i] <- 1
        }
      }
    }
    
    # Mass screening intervention
    if("mass_screening" %in% names(interventions) && day >= interventions$mass_screening$start_day) {
      if(day %% interventions$mass_screening$frequency == 0) {
        # Screen high-degree nodes first
        degrees <- rowSums(adj_matrix)
        high_risk <- which(degrees >= quantile(degrees, 0.8))
        
        for(idx in high_risk) {
          if(I_asymp_next[idx] == 1) {
            # Detect and treat asymptomatic cases
            I_asymp_next[idx] <- 0
            R_next[idx] <- 1
          }
        }
      }
    }
    
    # Seasonal Malaria Chemoprevention (SMC) for children
    if("SMC" %in% names(interventions)) {
      current_month <- ((day - 1) %/% 30) + 1
      if(current_month %in% interventions$SMC$months) {
        child_indices <- which(data$Age_Group == "Child")
        for(idx in child_indices) {
          if(runif(1) < interventions$SMC$coverage) {
            # SMC prevents infection
            if(E_next[idx] == 1) {
              E_next[idx] <- 0
              S_next[idx] <- 1
            }
          }
        }
      }
    }
    
    # Update states
    states[day + 1, , "S"] <- S_next
    states[day + 1, , "E"] <- E_next
    states[day + 1, , "I_symp"] <- I_symp_next
    states[day + 1, , "I_asymp"] <- I_asymp_next
    states[day + 1, , "R"] <- R_next
    states[day + 1, , "Immune"] <- Immune_next
  }
  
  return(list(states = states))
}

####Part 3 - MALARIA SIMULATIONS ####

# Malaria-specific parameters
beta <- 0.4           # Transmission rate (higher for malaria)
sigma <- 0.14         # Incubation rate (7-day incubation period)
gamma <- 0.07         # Recovery rate (14-day illness)
immunity_waning <- 0.003  # Immunity wanes over ~1 year
days <- 180           # 6 months simulation

initial_infected <- c("ZEN", "AGN", "BAS")  # Start with 3 infected

##Scenario 1: No interventions
set.seed(42)
result_baseline <- malaria_seir_model(adj_matrix, beta, sigma, gamma, immunity_waning, 
                                      initial_infected, days)

##Scenario 2: ITNs only
interventions_itn <- list(
  ITN = list(coverage = 0.6, efficacy = 0.5, start_day = 30)
)
set.seed(42)
result_itn <- malaria_seir_model(adj_matrix, beta, sigma, gamma, immunity_waning, 
                                 initial_infected, days, interventions_itn)

##Scenario 3: ITNs + IRS
interventions_combined <- list(
  ITN = list(coverage = 0.6, efficacy = 0.5, start_day = 30),
  IRS = list(coverage = 0.8, efficacy = 0.6, start_day = 60)
)
set.seed(42)
result_combined <- malaria_seir_model(adj_matrix, beta, sigma, gamma, immunity_waning, initial_infected, days, interventions_combined)

##Scenario 4: Full intervention package
interventions_full <- list(
  ITN = list(coverage = 0.7, efficacy = 0.5, start_day = 30),
  IRS = list(coverage = 0.8, efficacy = 0.6, start_day = 60),
  case_management = list(coverage = 0.8, start_day = 1),
  mass_screening = list(start_day = 90, frequency = 30),
  SMC = list(coverage = 0.9, months = c(6, 7, 8, 9))
)
set.seed(42)
result_full <- malaria_seir_model(adj_matrix, beta, sigma, gamma, immunity_waning, initial_infected, days, interventions_full)

####Part 4 - ANALYSIS AND VISUALIZATION ####

##Helper function for total infections (symptomatic + asymptomatic)
get_total_infected <- function(states, day) {
  sum(states[day + 1, , "I_symp"] + states[day + 1, , "I_asymp"])
}

##Analysis
cat("MALARIA TRANSMISSION ANALYSIS\n")
cat("=============================\n\n")

# Peak infections
baseline_peak <- max(rowSums(result_baseline$states[, , "I_symp"] + result_baseline$states[, , "I_asymp"]))
itn_peak <- max(rowSums(result_itn$states[, , "I_symp"] + result_itn$states[, , "I_asymp"]))
combined_peak <- max(rowSums(result_combined$states[, , "I_symp"] + result_combined$states[, , "I_asymp"]))
full_peak <- max(rowSums(result_full$states[, , "I_symp"] + result_full$states[, , "I_asymp"]))

cat("Peak Infections:\n")
cat("Baseline (no intervention):", baseline_peak, "\n")
cat("ITNs only:", itn_peak, "(", round((baseline_peak-itn_peak)/baseline_peak*100,1), "% reduction)\n")
cat("ITNs + IRS:", combined_peak, "(", round((baseline_peak-combined_peak)/baseline_peak*100,1), "% reduction)\n")
cat("Full package:", full_peak, "(", round((baseline_peak-full_peak)/baseline_peak*100,1), "% reduction)\n\n")

# Asymptomatic burden
asymp_baseline <- sum(result_baseline$states[, , "I_asymp"])
asymp_full <- sum(result_full$states[, , "I_asymp"])
cat("Asymptomatic burden:\n")
cat("Baseline:", asymp_baseline, "person-days\n")
cat("With full interventions:", asymp_full, "person-days\n\n")

##Visualization
plot_malaria_epidemic <- function(states_list, scenario_names) {
  days <- dim(states_list[[1]])[1] - 1
  
  # Create data frame for all scenarios
  df_list <- list()
  for(i in seq_along(states_list)) {
    total_infected <- rowSums(states_list[[i]][, , "I_symp"] + states_list[[i]][, , "I_asymp"])
    symptomatic <- rowSums(states_list[[i]][, , "I_symp"])
    asymptomatic <- rowSums(states_list[[i]][, , "I_asymp"])
    
    df_list[[i]] <- data.frame(
      Day = 0:days,
      Total_Infected = total_infected,
      Symptomatic = symptomatic,
      Asymptomatic = asymptomatic,
      Scenario = scenario_names[i]
    )
  }
  
  combined_df <- do.call(rbind, df_list)
  
  # Plot total infections
  p1 <- ggplot(combined_df, aes(x = Day, y = Total_Infected, color = Scenario)) +
    geom_line(size = 1.2) +
    labs(title = "Malaria Transmission: Impact of Interventions",
         subtitle = "Total infected individuals (symptomatic + asymptomatic)",
         x = "Day", y = "Number of Infected") +
    theme_minimal() +
    scale_color_manual(values = c("red", "blue", "green", "purple"))
  
  print(p1)
  
  # Plot symptomatic vs asymptomatic for baseline
  baseline_df <- df_list[[1]]
  baseline_long <- reshape2::melt(baseline_df[,c("Day", "Symptomatic", "Asymptomatic")], 
                                  id.vars = "Day")
  
  p2 <- ggplot(baseline_long, aes(x = Day, y = value, fill = variable)) +
    geom_area(alpha = 0.7) +
    labs(title = "Malaria Cases: Symptomatic vs Asymptomatic (Baseline)",
         x = "Day", y = "Number of Cases", fill = "Case Type") +
    theme_minimal() +
    scale_fill_manual(values = c("red", "orange"))
  
  print(p2)
}

# Plot results
states_list <- list(result_baseline$states, result_itn$states, 
                    result_combined$states, result_full$states)
scenario_names <- c("Baseline", "ITNs Only", "ITNs + IRS", "Full Package")

plot_malaria_epidemic(states_list, scenario_names)

##Network metrics analysis
cat("NETWORK IMPACT ANALYSIS\n")
cat("=======================\n")

degrees <- degree(g)
household_degrees <- tapply(degrees, data$Household, mean)
age_degrees <- tapply(degrees, data$Age_Group, mean)

cat("Average connections by age group:\n")
for(age in names(age_degrees)) {
  cat(age, ":", round(age_degrees[age], 2), "connections\n")
}

cat("\nHouseholds with highest connectivity (top 5):\n")
top_households <- sort(household_degrees, decreasing = TRUE)[1:5]
print(top_households)

##Cost-effectiveness analysis (simplified)
cat("\nCOST-EFFECTIVENESS ANALYSIS\n")
cat("===========================\n")

# Simplified costs (relative units)
costs <- list(
  ITN = 100,
  IRS = 200, 
  case_management = 50,
  mass_screening = 150,
  SMC = 80
)

scenarios_cost <- c(0, costs$ITN, costs$ITN + costs$IRS, 
                    sum(unlist(costs)))
infections_averted <- baseline_peak - c(baseline_peak, itn_peak, combined_peak, full_peak)

cost_effectiveness <- data.frame(
  Scenario = scenario_names,
  Cost = scenarios_cost,
  Infections_Averted = infections_averted,
  Cost_per_Infection_Averted = ifelse(infections_averted > 0, 
                                      scenarios_cost / infections_averted, NA)
)

print(cost_effectiveness)

cat("\n=== EXERCISES ===\n")
cat("1. Modify seasonal_factor to 1.5 and rerun - what happens?\n")
cat("2. Increase asymptomatic_rate to 0.6 - how does this affect control?\n") 
cat("3. What if ITN coverage was 0.9 instead of 0.6?\n")
cat("4. Identify the 3 highest-degree households - would targeting them be effective?\n")
cat("5. How would you modify the model for urban vs rural transmission patterns?\n")

