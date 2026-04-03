# Multiple-Comparisons-Method
# My Comparison Function: alpha/m^(m/(m+log10(m)))

my_comparison <- function(means, J, MSE, alpha = 0.05) {

  I <- length(means) #gets amount of groups
  df_error <- I*(J-1)
  m <- I*(I-1)/2  # number of pairwise comparisons
  adjustment <- m/(m+log10(m))
  alpha_adjusted <- alpha / m^(adjustment)  # divide alpha by number of comparisons
  t_crit <- qt(1 - alpha_adjusted / 2, df_error)
  critical_diff <- t_crit * sqrt(2 * MSE / J)
  diff_matrix <- outer(means, means, "-") # creates a matrix of all the differences of means from the samples
  sig_matrix <- (diff_matrix > critical_diff) | (diff_matrix < -critical_diff) #creates true false matrix based on if the differences of means are outside critical region
  return (as.data.frame(sig_matrix)) #returns the true/false matrix
}
