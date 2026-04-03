# Multiple-Comparisons-Method



# My Comparison Function
# The adjustment takes the number of pairwise comparisons (m) and divides it by a slightly larger number based off m
# Alpha is divided by the number of comparisons to the power of the adjustment
# The more pairwise comparisons, the smaller the alpha adjusted will be

my_comparison <- function(means, J, MSE, alpha = 0.05) { 
  I <- length(means) # Gets amount of groups
  df_error <- I*(J-1)
  m <- I*(I-1)/2  
  adjustment <- m/(m+log10(m)) 
  alpha_adjusted <- alpha / m^(adjustment) 
  t_crit <- qt(1 - alpha_adjusted / 2, df_error) 
  critical_diff <- t_crit * sqrt(2 * MSE / J) 
  diff_matrix <- outer(means, means, "-") 
  sig_matrix <- (diff_matrix > critical_diff) | (diff_matrix < -critical_diff)
  return (as.data.frame(sig_matrix)) 
}

# Simulate Family Wise Error
  # First calculate summary statistics then apply my method

simulate_fwer <- function(I, J, sigma = 1, n_sim = sim, alpha = 0.05) { 
  false_rejections <- 0 # Will be used to count times a significant difference was found
  for (sim in 1:n_sim) { # Repeat for as many simulations as instructed

    data <- matrix(rnorm(I * J, mean = 0, sd = sigma), nrow = I, ncol = J)  # Generate data under H0 (all means = 0)
    means <- rowMeans(data) # Find means of populations
    MSE <- mean(apply(data, 1, var)) # Use pooled variance to find mean standard error
    
    results <- my_comparison(means, J, MSE, alpha) # Calls my function and saves true/false matrix
    # Check if any rejection occurred
    if (any(results == TRUE)) { # Loop will check if there were any TRUEs
      false_rejections <- false_rejections + 1 # If there is a significant difference recorded, add it to counter
    }
  }
  return(false_rejections / n_sim) # Returns the proportion of times a simulation (incorrectly) reported a significant difference
}

# Create FWER data frame
fwerDF <- data.frame() 
fwerCalc <- function(i){
  fwerI <- c() 
  for(j in 5:20){
    fwer <- simulate_fwer (2, j, sigma = 1, n_sim = 1000, alpha = 0.05) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
    fwerI <- c(fwerI, fwer) 
  }
  return <- fwerI 
}

# Fill in FWER data frame with each I and J FWER value
I2 <- fwerCalc(2)
fwerDF <- rbind(I2)
I3 <- fwerCalc(3)
fwerDF <- rbind(fwerDF, I3)
I4<- fwerCalc(4)
fwerDF <- rbind(fwerDF,I4)
I5 <- fwerCalc(5)
fwerDF <- rbind(fwerDF,I5)
I6 <- fwerCalc(6)
fwerDF <- rbind(fwerDF,I6)
I7 <- fwerCalc(7)
fwerDF <- rbind(fwerDF,I7)
I8 <- fwerCalc(8)
fwerDF <- rbind(fwerDF,I8)
I9 <- fwerCalc(9)
fwerDF <- rbind(fwerDF,I9)
I10 <- fwerCalc(10)
fwerDF <- rbind(fwerDF,I10)
colnames(fwerDF) <- c("J5", "J6","J7","J8","J9","J10","J11","J12","J13","J14","J15","J16","J17","J18","J19","J20") 

# Simulate Power
simulate_pow <- function(I, J, sigma = 1, n_sim, alpha = 0.05, delta = 2) { 
  true_rejections <- 0   
  for (sim in 1:n_sim) {   
    row_means <- c(0, 0, 0, delta, delta) 
    mean_vector <- rep(row_means, each = J) 
    data <- matrix(rnorm(I * J, mean = mean_vector, sd = sigma), nrow = I,  ncol = J, byrow = TRUE)
    MSE <- mean(apply(data, 1, var)) 
    
    results <- my_comparison(row_means, J, MSE, alpha) 
    if (sum(results == TRUE) >0 ) {          
      true_rejections <- true_rejections + 1 
    }
  }
  return(true_rejections / n_sim) 
}

# Create power data frame
powDF <- data.frame() 
powCalc <- function(i){
  powI <- c() 
  for(j in 5:20){
    pow <- simulate_pow (2, j, sigma = 1, n_sim = 1000, alpha = 0.05)
    powI <- c(powI, pow) 
  }
  return <- powI 
}

# Fill in power data frame with each I and J power value
I2 <- powCalc(2)
powDF <- rbind(I2)
I3 <- powCalc(3)
powDF <- rbind(powDF, I3)
I4<- powCalc(4)
powDF <- rbind(powDF,I4)
I5 <- powCalc(5)
powDF <- rbind(powDF,I5)
I6 <- powCalc(6)
powDF <- rbind(powDF,I6)
I7 <- powCalc(7)
powDF <- rbind(powDF,I7)
I8 <- powCalc(8)
powDF <- rbind(powDF,I8)
I9 <- powCalc(9)
powDF <- rbind(powDF,I9)
I10 <- powCalc(10)
powDF <- rbind(powDF,I10)
colnames(powDF) <- c("J5", "J6","J7","J8","J9","J10","J11","J12","J13","J14","J15","J16","J17","J18","J19","J20")
