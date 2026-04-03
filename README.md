# Multiple-Comparisons-Method
This code provides an alternative to methods such as Tukey HSD by adjusting alpha based off the number of pairwise comparisons necessary. This code also shows the family wise error and power of my method, so it is clear how it functions.


my_comparison <- function(means, J, MSE, alpha = 0.05) { 
  #gets amount of groups
  I <- length(means) 
  df_error <- I*(J-1)
  #m is number of pairwise comparisons
  m <- I*(I-1)/2
  #the adjustment takes the number of pairwise comparisons and divides it by a slightly larger number based off m
  adjustment <- m/(m+log10(m)) 
  #alpha is divided by the number of comparisons to the power of the adjustment
  #the more pairwise comparisons, the smaller the alpha adjusted will be
  alpha_adjusted <- alpha / m^(adjustment)  
  #find critical t value using adjusted alpha
  t_crit <- qt(1 - alpha_adjusted / 2, df_error) 
  #find the critical difference for significance
  critical_diff <- t_crit * sqrt(2 * MSE / J)
  #creates a matrix of all the differences of means from the samples
  diff_matrix <- outer(means, means, "-") 
  #creates true false matrix based on if the differences of means are outside critical region
  sig_matrix <- (diff_matrix > critical_diff) | (diff_matrix < -critical_diff)
  #returns the true/false matrix
  return (as.data.frame(sig_matrix)) 
}

simulate_fwer <- function(I, J, sigma = 1, n_sim = sim, alpha = 0.05) { #simulates Family Wise Error
  false_rejections <- 0 #will be used to count times a significant difference was found
  for (sim in 1:n_sim) { #repeat for as many simulations as instructed

    data <- matrix(rnorm(I * J, mean = 0, sd = sigma), nrow = I, ncol = J)  # Generate data under H0 (all means = 0)
    # Calculate summary statistics
    means <- rowMeans(data) #Find means of populations
    MSE <- mean(apply(data, 1, var)) #Use pooled variance to find mean standard error
    
    # Apply my method
    results <- my_comparison(means, J, MSE, alpha) #calls my function and saves true/false matrix
    # Check if any rejection occurred
    if (any(results == TRUE)) { #loop will check if there were any TRUEs
      false_rejections <- false_rejections + 1 #if there is a significant difference recorded, add it to counter
    }
  }
  return(false_rejections / n_sim) #returns the proportion of times a simulation (incorrectly) reported a significant difference
}

fwerDF <- data.frame() #create an empty data frame
fwerCalc <- function(i){
  fwerI <- c() #create empty vector
  for(j in 5:20){
    fwer <- simulate_fwer (2, j, sigma = 1, n_sim = 1000, alpha = 0.05) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
    fwerI <- c(fwerI, fwer) #save the fwer for a given I and J
  }
  return <- fwerI #returns the FWER for a given I and J's 5-20
}

#fill in FWER data frame with each I and J FWER value
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
colnames(fwerDF) <- c("J5", "J6","J7","J8","J9","J10","J11","J12","J13","J14","J15","J16","J17","J18","J19","J20") #name columns for clarity

#POWER
simulate_pow <- function(I, J, sigma = 1, n_sim, alpha = 0.05, delta = 2) { #function to simulate power
  true_rejections <- 0    #set base count for rejections
  for (sim in 1:n_sim) {      #runs 1,000 simulations
    row_means <- c(0, 0, 0, delta, delta)   #creates the population means, with two that are different from the other three
    mean_vector <- rep(row_means, each = J) #creates vector of the means
    data <- matrix(rnorm(I * J, mean = mean_vector, sd = sigma), nrow = I,  ncol = J, byrow = TRUE)    #creates an I by J matrix of generated data using means and sigma = 1,
    MSE <- mean(apply(data, 1, var)) # Simplified; use pooled variance
    
    # Applying my method
    results <- my_comparison(row_means, J, MSE, alpha)      #runs my comparison function with given means and MSE, alpha =.05
    if (sum(results == TRUE) >0 ) {           #sees if my comparison correctly identified a difference
      true_rejections <- true_rejections + 1      #counts how many times there's a TRUE value in results matrix aka every time the null was correctly rejected
    }
  }
  return(true_rejections / n_sim) #returns the proportion of simulations that
}
powDF <- data.frame() #creates empty power data frame
powCalc <- function(i){
  powI <- c() #creates empty vector
  for(j in 5:20){
    pow <- simulate_pow (2, j, sigma = 1, n_sim = 1000, alpha = 0.05) #I = i, J =j sigma = 1, 100 simulations, alpha is .05
    powI <- c(powI, pow) #adds power for I J to the power I vector
  }
  return <- powI #returns vector with the power for I and for j 5-20
}

#fill in power data frame with each I and J power value
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
