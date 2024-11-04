library(tidyr)
library(dplyr)
library(ggplot2)
library(glue)
library(stringr)

######---------------------------------------------------------------------#####
######                  Visualization of Simulation Result                 #####
######---------------------------------------------------------------------#####

Scenarios = c("KN", "UN", "UC")

Labels = c("Known PS & No Censoring", 
           "Unknown PS & No Censoring", 
           "Unknown PS & Censoring") 

# Initialize empty data frames for both plots
policy_value_plot_data = data.frame()  # For policy value plot
beta_bias_plot_data = data.frame()  # For beta bias comparison plot

for(i in seq_along(Scenarios)){
  
  Scenario = Scenarios[i]
  Label = Labels[i]
  
  # List all CSV files in the directory
  files <- list.files(glue("data"), pattern = "^Output_.*_m\\d+_.*\\.csv$", full.names = FALSE)
  
  # Extract unique m values using regular expressions
  m_values <- sort(unique(as.numeric(sub("Output_.*_m(\\d+)_.*\\.csv", "\\1", files))))
  
  # Loop over different sample sizes
  for(m in m_values){  # Include multiple m values as needed
    
    # Get a list of all CSV files in the folder for the given sample size m
    csv_files <- list.files(path = glue("data"), pattern = glue("Output_{Scenario}_m{m}_.*.csv"), full.names = TRUE)
    
    # Read all CSV files into a list of data frames
    data_list <- lapply(csv_files, read.csv)
    
    # Combine all data frames into one
    combined_data <- bind_rows(data_list, .id = "simul")
    colnames(combined_data)[2] = "Method"
    combined_data$Method = factor(combined_data$Method, 
                                  levels = c("NoIntIPW", "ClusIPW", "addIPW", "Oracle"), 
                                  ordered = TRUE)
    
    combined_data$m = m  # Add the sample size column
    combined_data$Scenario = Label
    
    # Add to policy value plot data
    policy_value_plot_data = rbind(policy_value_plot_data, combined_data)
    
    
    ###----------- Beta Bias Comparison Plot ---------------###
    
    # Step 1: Identify beta columns dynamically
    beta_cols <- grep("^beta", names(combined_data), value = TRUE)
    
    # Step 2: Calculate bias by grouping by `simul`, `m`, and `Scenario`
    final_data <- combined_data %>%
      group_by(simul, m, Scenario) %>%
      # Calculate bias by subtracting Oracle values within each group
      mutate(across(all_of(beta_cols), 
                    ~ . - .[Method == "Oracle"], 
                    .names = "bias_{col}")) %>%
      ungroup()
    
    # Step 3: Reshape to long format for easy analysis
    final_data <- final_data %>%
      pivot_longer(cols = all_of(beta_cols), names_to = "beta", values_to = "beta_value") %>%
      pivot_longer(cols = starts_with("bias_beta"), names_to = "bias", values_to = "bias_value") %>%
      filter(str_remove(bias, "bias_") == beta) %>%
      filter(Method != "Oracle") %>%
      select(Method, beta, bias = bias_value, m, Scenario)
    
    # Add to beta bias plot data
    beta_bias_plot_data = rbind(beta_bias_plot_data, final_data)
    
  }
}

policy_value_plot_data$Scenario = factor(policy_value_plot_data$Scenario, levels = Labels)
beta_bias_plot_data$Scenario = factor(beta_bias_plot_data$Scenario, levels = Labels)


###----------- Policy Value Plot ---------------###

ggplot(policy_value_plot_data, 
       aes(x = as.factor(m), y = Value, fill = Method)) + 
  geom_boxplot(outlier.shape = 1, outlier.alpha = 0.7) + 
  theme_bw() + 
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("NoIntIPW" = "darkgreen", 
                               "ClusIPW" = "blue", 
                               "addIPW" = "red", 
                               "Oracle" = "purple")) +
  labs(title = "logT Value Maximization",
       x = "Sample Size", y = "Policy Value") +
  facet_wrap(~Scenario)

ggsave("logT_value_Maximization_Scenarios.jpeg", width = 12, height = 6)


###----------- Beta Bias Plot ---------------###

ggplot(beta_bias_plot_data, aes(x = beta, y = bias, fill = Method)) +
  geom_boxplot(outlier.shape = 1, outlier.alpha = 0.7) +
  scale_fill_manual(values = c("NoIntIPW" = "darkgreen", 
                               "ClusIPW" = "blue", 
                               "addIPW" = "red")) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
  labs(title = "Bias of Beta Estimates by Method and Sample Size",
       x = "Beta",
       y = "Bias") +
  facet_grid(Scenario ~ m) +  # Facet by `Scenario`
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity

ggsave("Beta_Bias_Plot_logT_Maximization.jpeg", width = 12, height = 6)
