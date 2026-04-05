# ZEWARI COMPARISON
# KATIE WOOD
# Zewari comparison provides an alternative to methods such as Tukey's HSD
# Comparison functions by dividing the alpha level by m^(m/(m+log10(m)) where m is the number of pariwise comparisons

set.seed (0)
zewari_comparison <- function(means, J, MSE, alpha = 0.05) {

  I <- length(means) #gets amount of groups
  df_error <- I*(J-1)
  m <- I*(I-1)/2  # number of pairwise comparisons
  adjustment <- m/(m+log10(m)) #adjustment will be a number slightly smaller than m itself
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
  return (fwerI) # Returns vector that contains all FWER's for I and J 6-20
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
fwerCheck37 <- data.frame() #creates empty data frame to fill
l_37 <- 0
h_37<- 0
r_37 <- 0
l_46 <- 0
h_46 <- 0
r_46 <- 0
# Fills fwerCheck37 with "H" for a value higher than .07, "L" for a value lower than .03, and "R" for a value that is between .03 and.07
for (i in 1:nrow(fwerDF)){
  for(j in 1:ncol(fwerDF)){
    if (fwerDF[i,j] > .07){
      fwerCheck37[i,j] <- "H"
      h_37 <- h_37+1
    }
    if (fwerDF[i,j] < .03){
      fwerCheck37[i,j] <- "L"
      l_37 <- l_37+1
    }
    if(fwerDF[i,j] >=.03 && fwerDF[i,j] <= .07){
      fwerCheck37[i,j] <- "R"
      r_37 <- r_37+1
    }
  }
}
colnames(fwerCheck37) <- c("J5", "J6","J7","J8","J9","J10","J11","J12","J13","J14","J15","J16","J17","J18","J19","J20")
rownames(fwerCheck37) <- c("I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9", "I10")

# Checks if FWER is between .3 and .7
fwerCheck46 <- data.frame()

# Fills fwerCheck46 with "H" for a value higher than .06, "L" for a value lower than .06, and "R" for a value that is between .04 and.06
for (i in 1:nrow(fwerDF)){
  for(j in 1:ncol(fwerDF)){
    if (fwerDF[i,j] > .06){
      fwerCheck46[i,j] <- "H"
      h_46 <- h_46+1
    }
    if (fwerDF[i,j] < .04){
      fwerCheck46[i,j] <- "L"
      l_46 <- l_46+1
    }
    if(fwerDF[i,j] >=.04 && fwerDF[i,j] <= .06){
      fwerCheck46[i,j] <- "R"
      r_46 <- r_46+1
    }
  }
}

#Calculates average FWER for scenarios I 2-9, J 5-20
avgFWER <- sum(fwerDF)/(nrow(fwerDF)*ncol(fwerDF))
#Creates a vector that has how many FWERS were lower, between, higher than each range, and the average FWER
results <- c(l_37, h_37, r_37, l_46, h_46, r_46, avgFWER)

##############################################################################################################
# Switch from simulating FWER to simulating Power
pow <- 0
simulate_pow <- function(I, J, sigma = 1, n_sim, alpha = 0.05, delta = 1) { #function to simulate power
  true_rejections <- 0    #set base count for rejections
  for (sim in 1:n_sim) {      #runs n_sim simulations
    row_means <- c()
    for (i in 1:(I-1))
    {
      row_means <- c(row_means, 0)
    }
    row_means <- c(row_means, delta)
    mean_vector <- rep(row_means, each = J) # Creates vector of the means
    data <- matrix(rnorm(I * J, mean = mean_vector, sd = sigma), nrow = I,  ncol = J, byrow = TRUE) #Creates data with one population whose mean is different by delta
    MSE <- mean(apply(data, 1, var)) # Simplified; using pooled variance

    # Applying my method
    means <- rowMeans(data)
    results <- zewari_comparison(means, J, MSE, alpha)      #runs my comparison function with given means and MSE, alpha =.05
    if (sum(results == TRUE) >0 ) {           #checks if matrix has TRUE(s) in it and thus found a difference
      true_rejections <- true_rejections + 1      #counts how many times there's a true value in results matrix aka every time the null (all populations equal) was correctly rejected
    }
  }
  return(true_rejections / n_sim) #returns the proportion of simulations that correctly found a significant difference
}

#Runs power simulations for I 2:9 and J 5:20
powSimAll <- function(){
  powDF <- data.frame() #creates empty power data frame
  powCalc <- function(i){
    powI <- c()
    for(j in 5:20){
      # Here is where you can adjust delta
      pow <- simulate_pow (i, j, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1.5) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
      powI <- c(powI, pow)
    }
    return (powI)
  }

  # Fills power data frame with the average power for each value of I and J
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

  return(powDF)
}

powdataFrame <- powSimAll()
avgPow <- sum(powdataFrame)/(nrow(powdataFrame)*ncol(powdataFrame)) #finds average power for all scenarios I: 2-9 and J: 5-20

#####################################################################################################################################
#Runs specific scenarios for FWER (one could also display fwerDF and pull values from that)
fwerA <- simulate_fwer (3, 10, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
fwerB <- simulate_fwer (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
fwerC <- simulate_fwer (7, 10, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
fwerD <- simulate_fwer (5, 5, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
fwerE <- simulate_fwer (5, 20, sigma = 1, n_sim = 10000, alpha = 0.05) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
#show output of scenarios
fwerA
fwerB
fwerC
fwerD
fwerE


#simulate power for deltas 0.5-2.0, I=5, J=10
pow25 <- simulate_pow (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = .25) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
pow05 <- simulate_pow (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = .5) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
pow75 <- simulate_pow (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = .75) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
pow1 <- simulate_pow (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
pow125 <- simulate_pow (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1.25) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
pow15 <- simulate_pow (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1.5) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
pow175 <- simulate_pow (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1.75) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
pow2 <- simulate_pow (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 2) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05

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
  for (sim in 1:n_sim) {      #runs n_sim simulations
    row_means <- c(rep(0, I-1), delta)   #input data for each test
    mean_vector <- rep(row_means, each = J) #get means
    data <- matrix(rnorm(I * J, mean = mean_vector, sd = sigma), nrow = I,  ncol = J, byrow = TRUE)    #creates an I by J matrix of generated data that has the correct means and sd = 1,
    #convert data to data frame
    Y <- as.vector(t(data))
    Population <- factor(rep(paste0("P",1:I), each = J))
    dataAnova <- data.frame(Y, Population)
    #run ANOVA and Tukey
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
powT25 <- simulate_powT (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = .25) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
powT05 <- simulate_powT (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = .5) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
powT75 <- simulate_powT (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = .75) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
powT1 <- simulate_powT (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
powT125 <- simulate_powT (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1.25) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
powT15 <- simulate_powT (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1.5) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
powT175 <- simulate_powT (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 1.75) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05
powT2 <- simulate_powT (5, 10, sigma = 1, n_sim = 10000, alpha = 0.05, delta = 2) #I = i, J =j sigma = 1, 10000 simulations, alpha is .05

#Show output for Tukey's simulated power
powT05
powT75
powT1
powT125
powT15
powT175
powT2

#Create vector of Tukey's power
tukPow <- c(powT05,powT75,powT1,powT125, powT15,  powT175, powT2)



#Plotting Zewari's vs Tukey's power
delta <- c(0.5, .75, 1, 1.25, 1.5, 1.75, 2)
plot(x = delta, y = myPow,
     type = "l",
     col = "red",
     ylim = c(-1,2),
     main = "Tukey and Zewari Power",
     xlab = "Delta",
     ylab = "Power")

lines(x = delta, y = tukPow,
      col = "blue")

legend("topleft",
       legend = c("Zewari Method", "Tukey"),
       col = c("red", "blue"),
       lty = 1,
       cex = 0.8)
