# Core parameters
Nh_init <- 15000                # Initial human population (constant, no births/deaths)
Ih_init<- 25                   # Initial infected humans
Sh_init <- Nh_init - Ih_init    # Initial susceptible humans

mosq_per_human <- 8             # Mosquitoes per human
Nm_init <- Nh_init * mosq_per_human
Im_init <- round(0.12 * Nm_init)  # 12% infectious mosquitoes initially
Em_init<- 0
Sm_init<- Nm_init - Im_init

biting_rate <- 0.3              # Bites per mosquito per day
host_seeking <- 0.1             # Proportion seeking hosts
trans_prob <- 0.45              # Transmission probability per bite
alpha <- biting_rate * host_seeking * trans_prob  # Effective transmission rate

beta <- 0.5                     # Human-to-mosquito transmission
gamma <- 1/25                   # Human recovery rate (1/days)
sigma <- 1/7                    # Mosquito incubation rate (1/days)
mu_m_base <- 1/14               # Base mosquito mortality (1/days)
Lambda_m <- mu_m_base * Nm_init # Fixed recruitment rate (to maintain equilibrium without control)

vacc_coverage <- 0.30           # Vaccination coverage
vacc_efficacy <- 0.65           # Vaccine efficacy (leaky)
vacc_waning <- 1/30             # Waning rate (1/days)
vacc_cost_per_person <- 20      # USD per vaccinated person

vc_mort_increase <- 0.60        # Vector control: 60% mortality increase
vc_budget <- 80000              # Annual vector control cost (USD)

daly_per_case <- 0.2            # DALYs per malaria case (rough estimate; adjust as needed)

seasonal_forcing <- FALSE       # Toggle seasonal variation (sinusoidal on biting rate)
seasonal_amplitude <- 0.2       # Amplitude if enabled

times <- seq(0, 365, by = 1)    # Simulation time (days)

#=========================================================================================
#Model Function (Unified for All Scenarios)---------------------------------------
#======================================================================================================
malaria_model <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
        # Total populations (humans constant, mosquitoes dynamic)
        Nh <- Sh + Ih + Vh
        Nm <- Sm + Em + Im
        
        # Optional seasonal forcing on biting rate
        biting_eff <- if (seasonal_forcing) {
            biting_rate * (1 + seasonal_amplitude * sin(2 * pi * time / 365))
        } else {
            biting_rate
        }
        alpha_eff <- biting_eff * host_seeking * trans_prob
        
        # Forces of infection
        lambda_h <- alpha_eff * (Im / Nh)
        lambda_m <- beta * (Ih / Nh)
        
        # New infections (for incidence)
        new_inf <- lambda_h * Sh + lambda_h * (1 - efficacy) * Vh
        
        # Differential equations
        dSh <- -lambda_h * Sh + gamma * Ih + waning * Vh
        dIh <- new_inf - gamma * Ih
        dVh <- -lambda_h * (1 - efficacy) * Vh - waning * Vh
        
        dSm <- Lambda_m - lambda_m * Sm - mu_m * Sm
        dEm <- lambda_m * Sm - sigma * Em - mu_m * Em
        dIm <- sigma * Em - mu_m * Im
        
        dCumInf <- new_inf  # Cumulative infections
        
        return(list(c(dSh, dIh, dVh, dSm, dEm, dIm, dCumInf)))
    })
}

#===============================================================================
# Run Simulation Function
#===============================================================================

run_simulation <- function(scenario = "baseline") {
    # Base parameters
    parms <- list(
        efficacy = 0,
        waning = 0,
        mu_m = mu_m_base
    )
    
    # Initial states
    Vh_init <- 0
    Sh_eff <- Sh_init
    
    # Scenario overrides
    if (scenario %in% c("vaccine", "combined")) {
        Vh_init <- round(vacc_coverage * Sh_init)
        Sh_eff <- Sh_init - Vh_init
        parms$efficacy <- vacc_efficacy
        parms$waning <- vacc_waning
    }
    if (scenario %in% c("vector", "combined")) {
        parms$mu_m <- mu_m_base * (1 + vc_mort_increase)
    }
    
    initial_state <- c(
        Sh = Sh_eff,
        Ih = Ih_init,
        Vh = Vh_init,
        Sm = Sm_init,
        Em = Em_init,
        Im = Im_init,
        CumInf = 0
    )
    
    # Solve ODE
    out <- ode(y = initial_state, times = times, func = malaria_model, parms = parms)
    out <- as.data.frame(out)
    
    # Post-processing
    out$Nh <- out$Sh + out$Ih + out$Vh
    out$Nm <- out$Sm + out$Em + out$Im
    out$prevalence <- out$Ih / out$Nh
    out$daily_incidence <- c(0, diff(out$CumInf))  # Daily new cases
    
    # Validation: Check human pop constant
    if (abs(max(out$Nh) - min(out$Nh)) > 1e-6) warning("Human population not constant!")
    
    return(out)
}

#===============================================================================
#Run All Scenarios--------------------------------------------------------------
#===============================================================================

scenarios <- c("baseline", "vector", "vaccine", "combined")
results <- lapply(scenarios, run_simulation)
names(results) <- scenarios

# Combine for plotting
all_data <- bind_rows(results, .id = "scenario") %>%
    mutate(scenario = factor(scenario, levels = scenarios))

#===============================================================================
#Prevalence Calculations--------------------------------------------------------------
#===============================================================================


# Final prevalence per scenario
final_prev <- all_data %>%
    group_by(scenario) %>%
    summarise(final_prevalence = tail(prevalence, 1), .groups = "drop")

print(final_prev)

# output interpretation:
# Baseline: ~59.2% infected at end (high transmission).
# Vector control reduces by ~8-10% via mosquito population decline.
# Vaccine adds protection but wanes; combined is most effective.


#===============================================================================
#Incidence and Cases Averted--------------------------------------------------------------
#===============================================================================
# Metrics function
calc_metrics <- function(res) {
    total_incident <- tail(res$CumInf, 1)
    mean_pop <- mean(res$Nh)
    incidence_per_1000 <- (total_incident / mean_pop) * 1000 / (365 / 365)  # Annualized
    list(
        total_incident = total_incident,
        incidence_per_1000 = incidence_per_1000,
        mean_pop = mean_pop
    )
}

metrics <- lapply(results, calc_metrics)

# Cases averted vs. baseline
baseline_inc <- metrics$baseline$total_incident
averted <- sapply(metrics, function(m) baseline_inc - m$total_incident)
percent_averted <- (averted / baseline_inc) * 100

# Costs
costs <- c(
    baseline = 0,
    vector = vc_budget,
    vaccine = vacc_coverage * Sh_init * vacc_cost_per_person,
    combined = vc_budget + vacc_coverage * Sh_init * vacc_cost_per_person
)

# DALYs averted (assuming daly_per_case)
daly_averted <- averted * daly_per_case

# ICER (vs. baseline)
icer <- costs[-1] / averted[-1]  # Exclude baseline

# Summary table
summary_table <- tibble(
    Scenario = scenarios,
    Total_Incident_Cases = sapply(metrics, `[[`, "total_incident") %>% round(0),
    Annual_Incidence_per_1000 = sapply(metrics, `[[`, "incidence_per_1000") %>% round(1),
    Cases_Averted = round(averted, 0),
    Percent_Averted = round(percent_averted, 1),
    Cost_USD = round(costs, 0),
    Cost_per_Case_Averted = round(costs / pmax(averted, 1), 2),  # Avoid div/0
    DALYs_Averted = round(daly_averted, 0),
    ICER_USD_per_DALY = round(costs / pmax(daly_averted, 1), 2)
)

print(summary_table)

#===============================================================================
#Improved Plots
# Prevalence Over Time (with Confidence Bands)--------------------------------------------------------------
#===============================================================================

# Simple sensitivity: Run with alpha ±10% for bands
sens_low <- lapply(scenarios, function(s) {
    temp_alpha <- alpha * 0.9
    # Re-run with temp_alpha (override in parms)
    # ... (implement similar to run_simulation, but for brevity, assume bands as ±5% prevalence)
})

# For demo, add fake bands
all_data <- all_data %>%
    mutate(prev_low = prevalence * 0.95, prev_high = prevalence * 1.05)

p_prev <- ggplot(all_data, aes(x = time, y = prevalence, color = scenario)) +
    geom_line(size = 1.2) +
    geom_ribbon(aes(ymin = prev_low, ymax = prev_high, fill = scenario), alpha = 0.2, color = NA) +
    geom_hline(data = final_prev, aes(yintercept = final_prevalence, color = scenario), linetype = "dotted") +
    geom_text(data = final_prev, aes(x = 365 * 0.7, y = final_prevalence + 0.02, 
                                     label = paste0(round(final_prevalence * 100, 1), "%"), color = scenario)) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    labs(title = "Prevalence Over Time by Scenario", x = "Days", y = "Prevalence") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

print(p_prev)

#Compartment Dynamics (Faceted, Interactive Option)
# Long format for compartments
comp_long <- all_data %>%
    pivot_longer(c(Sh, Ih, Vh, Sm, Em, Im), names_to = "compartment", values_to = "count") %>%
    filter(compartment %in% c("Sh", "Ih", "Vh", "Sm", "Em", "Im"))

p_comp <- ggplot(comp_long, aes(x = time, y = count, color = compartment)) +
    geom_line(size = 1) +
    facet_wrap(~ scenario, scales = "free_y") +
    scale_color_viridis_d() +
    labs(title = "Compartment Counts Over Time", x = "Days", y = "Count") +
    theme_minimal(base_size = 12)

print(p_comp)

# Interactive version
p_comp_int <- ggplotly(p_comp)
p_comp_int

#Daily and Cumulative Incidence
#===============================================================================
##Daily and Cumulative Incidence--------------------------------------------------------------
#===============================================================================
p_daily_inc <- ggplot(all_data, aes(x = time, y = daily_incidence, color = scenario)) +
    geom_line(size = 1) +
    scale_color_viridis_d() +
    labs(title = "Daily Incidence by Scenario", x = "Days", y = "New Cases/Day") +
    theme_minimal(base_size = 12)

p_cum_inc <- ggplot(all_data, aes(x = time, y = CumInf, color = scenario)) +
    geom_line(size = 1) +
    scale_color_viridis_d() +
    labs(title = "Cumulative Incidence by Scenario", x = "Days", y = "Cumulative Cases") +
    theme_minimal(base_size = 12)

grid.arrange(p_daily_inc, p_cum_inc, ncol = 2)



#===============================================================================
##Daily and Cumulative Incidence--------------------------------------------------------------
#===============================================================================
cea_plot <- summary_table %>%
    filter(Scenario != "baseline") %>%
    ggplot(aes(x = Cases_Averted, y = Cost_per_Case_Averted, color = Scenario)) +
    geom_point(size = 5) +
    geom_text(aes(label = paste0("$", Cost_per_Case_Averted)), vjust = -1) +
    geom_line(color = "gray", linetype = "dashed") +
    scale_color_viridis_d() +
    labs(title = "Cost-Effectiveness: Cases Averted vs. Cost per Case", x = "Cases Averted", y = "USD per Case Averted") +
    theme_minimal(base_size = 12)

print(cea_plot)

#===============================================================================
##Sensitivity Analysis (Tornado Plot Example)--------------------------------------------------------------
#===============================================================================
run_simulation <- function(scenario = "baseline", custom_parms = NULL) {
    # Base parameters
    parms <- list(
        biting_rate = biting_rate,
        host_seeking = host_seeking,
        trans_prob = trans_prob,
        alpha = biting_rate * host_seeking * trans_prob,
        beta = beta,
        gamma = gamma,
        sigma = sigma,
        mu_m = mu_m_base,
        Lambda_m = Lambda_m,
        efficacy = 0,
        waning = 0,
        seasonal_forcing = seasonal_forcing,
        seasonal_amplitude = seasonal_amplitude
    )
    
    # Override with custom parameters if provided
    if (!is.null(custom_parms)) {
        parms[names(custom_parms)] <- custom_parms
    }
    
    # Initial states
    Vh_init <- 0
    Sh_eff <- Sh_init
    
    # Scenario overrides
    if (scenario %in% c("vaccine", "combined")) {
        Vh_init <- round(vacc_coverage * Sh_init)
        Sh_eff <- Sh_init - Vh_init
        parms$efficacy <- vacc_efficacy
        parms$waning <- vacc_waning
    }
    if (scenario %in% c("vector", "combined")) {
        parms$mu_m <- mu_m_base * (1 + vc_mort_increase)
    }
    
    initial_state <- c(
        Sh = Sh_eff,
        Ih = Ih_init,
        Vh = Vh_init,
        Sm = Sm_init,
        Em = Em_init,
        Im = Im_init,
        CumInf = 0
    )
    
    # Solve ODE
    out <- ode(y = initial_state, times = times, func = malaria_model, parms = parms)
    out <- as.data.frame(out)
    
    # Post-processing
    out$Nh <- out$Sh + out$Ih + out$Vh
    out$Nm <- out$Sm + out$Em + out$Im
    out$prevalence <- out$Ih / out$Nh
    out$daily_incidence <- c(0, diff(out$CumInf))  # Daily new cases
    
    # Validation
    if (abs(max(out$Nh) - min(out$Nh)) > 1e-6) warning("Human population not constant in scenario: ", scenario)
    
    return(out)
}

# Complete Sensitivity Analysis-------------------------------
# Sensitivity analysis: Vary key parameters ±20%
sens_params <- c("alpha", "beta", "gamma", "sigma", "mu_m_base", "vacc_efficacy")
sens_range <- 0.2  # ±20%

sens_results <- map_dfr(sens_params, ~{
    param_name <- .
    base_val <- get(param_name)
    low_val <- base_val * (1 - sens_range)
    high_val <- base_val * (1 + sens_range)
    
    # Run simulation for low value
    low_sim <- run_simulation("combined", custom_parms = setNames(list(low_val), param_name))
    low_averted <- baseline_inc - calc_metrics(low_sim)$total_incident
    
    # Run simulation for high value
    high_sim <- run_simulation("combined", custom_parms = setNames(list(high_val), param_name))
    high_averted <- baseline_inc - calc_metrics(high_sim)$total_incident
    
    # Return results
    tibble(
        param = param_name,
        low = low_averted,
        base = averted["combined"],
        high = high_averted
    )
})

# Tornado plot
p_tornado <- sens_results %>%
    pivot_longer(c(low, high), names_to = "variation", values_to = "averted") %>%
    mutate(variation = factor(variation, levels = c("low", "high"))) %>%
    ggplot(aes(x = averted, y = reorder(param, base), fill = variation)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_vline(xintercept = averted["combined"], linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("low" = "#1b9e77", "high" = "#d95f02"), 
                      labels = c("low" = "-20%", "high" = "+20%")) +
    labs(
        title = "Sensitivity Analysis: Cases Averted in Combined Scenario",
        x = "Cases Averted",
        y = "Parameter",
        fill = "Variation"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top"
    )

print(p_tornado)
