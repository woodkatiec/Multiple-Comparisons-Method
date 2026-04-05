set.seed (000)
zewari_comparison <- function(means, J, MSE, alpha = 0.05) {

  I <- length(means) #gets amount of groups
  df_error <- I*(J-1)
  m <- I*(I-1)/2  # number of pairwise comparisons
  adjustment <- m/(m+log10(m)) #adjustment will be a number slgihtly smaller than m itself
  alpha_adjusted <- alpha / m^adjustment     #divide alpha by m^adjustment to shrink alpha
  t_crit <- qt(1 - alpha_adjusted / 2, df_error) #calculate critical t-value
  critical_diff <- t_crit * sqrt(2 * MSE / J) #calculate critical difference
  diff_matrix <- outer(means, means, "-")  #creates a matrix of all the differences of means from the samples
  sig_matrix <- (diff_matrix > critical_diff) | (diff_matrix < -critical_diff) #creates true false matrix based on if the differences of means are outside critical region
  return (as.data.frame(sig_matrix)) #returns the true/false matrix
}

#FWER
simulate_fwer <- function(I, J, sigma = 1, n_sim = 10000, alpha = 0.05) {
    false_rejections <- 0 #will be used to count times a significant difference was found
    for (sim in 1:n_sim) { #repeat for as many simulations as instructed
      # Generate data under H0 (all means = 0)
      data <- matrix(rnorm(I * J, mean = 0, sd = sigma), nrow = I, ncol = J)
      # Calculate summary statistics
      means <- rowMeans(data) # Find means of populations
      MSE <- mean(apply(data, 1, var)) # Simplified; use pooled variance to find mean standard error
      # Apply my method
      results <- zewari_comparison(means, J, MSE, alpha) #calls my function and saves true/false matrix
      # Check if any rejection occurred
      if (any(results == TRUE)) { # Loop will add to counter if there were any TRUEs
        false_rejections <- false_rejections + 1 # If there is a significant difference recorded, add it to counter
      }
    }
    return(false_rejections / n_sim) # Returns the proportion of times a simulation (incorrectly) reported a significant difference
  }

fwerDF <- data.frame()
fwerCalc <- function(i){ #calcuates fwer for J 6-20 for a given I
  fwerI <- c()
  for(j in 5:20){
    fwer <- simulate_fwer (i, j, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
    fwerI <- c(fwerI, fwer)
    }
    return <- fwerI # Returns vector that contains all FWER's for I and J 6-20
  }

#Build a matrix of FWER for I 2-9 and J 6-20 using my comparison method
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
  colnames(fwerDF) <- c("J5", "J6","J7","J8","J9","J10","J11","J12","J13","J14","J15","J16","J17","J18","J19","J20") #Columns are renamed for clarity

# Checks if FWER is between .3 and .7
  fwerCheck <- data.frame() #creates empty data frame to fill
  l_37 <- 0
  h_37<- 0
  r_37 <- 0
  l_46 <- 0
  h_46 <- 0
  r_46 <- 0
  for (i in 1:9){
    for(j in 1:16){
      if (fwerDF[i,j] > .07){
        fwerCheck[i,j] <- "H"
        h_37 <- h_37+1
      }
      if (fwerDF[i,j] < .03){
        fwerCheck[i,j] <- "L"
        l_37 <- l_37+1
      }
      if(fwerDF[i,j] >=.03 && fwerDF[i,j] <= .07){
        fwerCheck[i,j] <- "R"
        r_37 <- r_37+1
      }
    }
  }
  colnames(fwerCheck) <- c("J5", "J6","J7","J8","J9","J10","J11","J12","J13","J14","J15","J16","J17","J18","J19","J20")
  rownames(fwerCheck) <- c("I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9", "I10")

# Checks if FWER is between .4 and .6
  fwerCheck <- data.frame()
  for (i in 1:9){
    for(j in 1:16){
      if (fwerDF[i,j] > .06){
        fwerCheck[i,j] <- "H"
        h_46 <- h_46+1
      }
      if (fwerDF[i,j] < .04){
        fwerCheck[i,j] <- "L"
        l_46 <- l_46+1
      }
      if(fwerDF[i,j] >=.04 && fwerDF[i,j] <= .06){
        fwerCheck[i,j] <- "R"
        r_46 <- r_46+1
      }
    }
  }

#Calculates average FWER for scenarios I 2-9, J 5-20
  avgFWER <- sum(fwerDF)/(15*9)
#Creates a vector that has how many FWERS were lower, between, higher than each range, and the average FWER
  results <- c(l_37, h_37, r_37, l_46, h_46, r_46, avgFWER)

##############################################################################################################
pow <- 0

simulate_pow <- function(I, J, sigma = 1, n_sim, alpha = 0.05, delta = 1) { #function to simulate power
    true_rejections <- 0    #set base count for rejections
    for (sim in 1:n_sim) {      #runs 1,000 simulations
      row_means <- c()
      for (i in 1:(I-1))
        {
          row_means <- c(row_means, 0)
        }
    row_means <- c(row_means, delta)
    mean_vector <- rep(row_means, each = J) # Creates vector of the means
    data <- matrix(rnorm(I * J, mean = mean_vector, sd = sigma), nrow = I,  ncol = J, byrow = TRUE) #Creates data with one populatoin whose mean is different by delta      MSE <- mean(apply(data, 1, var)) # Simplified; use pooled variance

    # Applying Zewari method
      results <- my_comparison(row_means, J, MSE, alpha)      #runs my comparison function with given means and MSE, alpha =.05
      if (sum(results == TRUE) >0 ) {           #used sum(results == TRUE) > 0 because the original function was returning false even when I printed and the results matrix had trues in it
        true_rejections <- true_rejections + 1      #counts how many times there's a true value in results matrix aka every time the null (all pops equal) was rejected
        }
      }
      return(true_rejections / n_sim) #returns the proportion of simulations that
    }

powDF <- data.frame()
powCalc <- function(i){
  powI <- c()
  ival <- i
  for(j in 5:20){
    #this is where delta can be adjusted
      pow <- simulate_pow (ival, j, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1.5) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
      powI <- c(powI, pow)
      }
      return (powI)
    }

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


  avgPow <- sum(powDF)/(16*9)
  avgPow

#####################################################################################################################################
#running specific scenarios for FWER (can also display fwerDF and take values from that)
   fwerA <- simulate_fwer (3, 10, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  fwerB <- simulate_fwer (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  fwerC <- simulate_fwer (7, 10, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  fwerD <- simulate_fwer (5, 5, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  fwerE <- simulate_fwer (5, 20, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
#show output of scenarios
  fwerA
  fwerB
  fwerC
  fwerD
  fwerE

#simulate power for when three populations have same mean, and two have different mean
  simulate_pow <- function(I, J, sigma = 1, n_sim, alpha = 0.05, delta) { #function to simulate power
    true_rejections <- 0    #set base count for rejections
    for (sim in 1:n_sim) {      #runs 1,000 simulations
      row_means <- c(0, 0, 0, delta, delta)   #I manually inputed data/changed means for each test
      mean_vector <- rep(row_means, each = J) #creates vector of the means
      data <- matrix(rnorm(I * J, mean = mean_vector, sd = sigma), nrow = I,  ncol = J, byrow = TRUE)    #creates an I by J (5 by 10) matrix of generated data using means and sd = 1,
      MSE <- mean(apply(data, 1, var)) # Simplified; use pooled variance

      # Applying my method
      results <- zewari_comparison(row_means, J, MSE, alpha)      #runs my comparison function with given means and MSE, alpha =.05
      if (sum(results == TRUE) >0 ) {           #used sum(results == TRUE) > 0 because the original function was returning false even when I printed and the results matrix had trues in it
        true_rejections <- true_rejections + 1      #counts how many times there's a true value in results matrix aka every time the null (all pops equal) was rejected
      }
    }
    return(true_rejections / n_sim) #returns the proportion of simulations that
  }

#simulate power for deltas 0.5-2.0
  pow05 <- simulate_pow (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = .5) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  pow75 <- simulate_pow (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = .75) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  pow1 <- simulate_pow (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 1) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  pow125 <- simulate_pow (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 1.25) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  pow15 <- simulate_pow (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 1.5) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  pow175 <- simulate_pow (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 1.75) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  pow2 <- simulate_pow (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 2) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05

# Show output of specific power simulations
  pow05
  pow75
  pow1
  pow125
  pow15
  pow175
  pow2

# Create vector for Zewari method's power
  myPow <- c(pow05,pow75,pow1,pow125, pow15,  pow175, pow2)

#Tukey's HSD power simulation
simulate_powT <- function(I, J, sigma = 1, n_sim, alpha = 0.05, delta) {
  true_rejections <- 0 #set base count for rejections
  for (sim in 1:n_sim) {      #runs 1,000 simulations
    row_means <- c(0, 0, 0, delta, delta)   #input data for each test
    mean_vector <- rep(row_means, each = J) #get means
    data <- matrix(rnorm(I * J, mean = mean_vector, sd = sigma), nrow = I,  ncol = J, byrow = TRUE)    #creates an I by J (5 by 10) matrix of generated data that has the correct means and sd = 1,
    dataAnova <- data.frame(
      Y = c(data[1,1],data[1,2],data[1,3],data[1,4],data[1,5],data[1,6],data[1,7],data[1,8],data[1,9],data[1,10],
            data[2,1],data[2,2],data[2,3],data[2,4],data[2,5],data[2,6],data[2,7],data[2,8],data[2,9],data[2,10],
            data[3,1],data[3,2],data[3,3],data[3,4],data[3,5],data[3,6],data[3,7],data[3,8],data[3,9],data[3,10],
            data[4,1],data[4,2],data[4,3],data[4,4],data[4,5],data[4,6],data[4,7],data[4,8],data[4,9],data[4,10],
            data[5,1],data[5,2],data[5,3],data[5,4],data[5,5],data[5,6],data[5,7],data[5,8],data[5,9],data[5,10]),
      Population = factor(rep(c("P1","P2","P3","P4", "P5"), each = 10))
      )
    aov_model <- aov(Y ~ Population, data = dataAnova) #runs anova test
    results <-TukeyHSD(aov_model, conf.level = 0.95)   #runs Tukey test
    p_vals <- results$Population[,4] #gets just the p values from the Tukey test
    tf <- p_vals < .05
    if (sum(tf == TRUE) > 0) {        #adds to the counter if it was significant
      true_rejections <- true_rejections + 1      #counts how many times there's a true value in results matrix
      }
    }
    return(true_rejections / n_sim) #returns proportion of times difference was successfully found
  }


#Tukey's power simulated for specific scenarios and deltas
  powT05 <- simulate_powT (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = .5) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  powT75 <- simulate_powT (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = .75) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  powT1 <- simulate_powT (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 1) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  powT125 <- simulate_powT (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 1.25) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  powT15 <- simulate_powT (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 1.5) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  powT175 <- simulate_powT (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 1.75) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05
  powT2 <- simulate_powT (5, 10, sigma = 1, n_sim = 1000, alpha = 0.05, delta = 2) #I = i, J =j sigma = 1, 1000 simulations, alpha is .05

#Show output for Tukey's simulated power
  powT05
  powT75
  powT1
  powT125
  powT15
  powT175
  powT2

  tukPow <- c(powT05,powT75,powT1,powT125, powT15,  powT175, powT2) #vector of tukeys power



#Plotting Zewari vs Tukey's power
  delta <- c(.05, .75, 1, 1.25, 1.5, 1.75, 2)
    plot(x = delta, y = myPow,
       type = "l",
       col = "pink",
       ylim = c(-1,2),
       main = "Tukey and Zewari Power",
       xlab = "Delta",
       ylab = "Power")

  lines(x = delta, y = tukPow,
        col = "blue")

  legend("topleft",
         legend = c("Zewari Method", "Tukey"),
         col = c("pink", "blue"),
         lty = 1,
         cex = 0.8)
