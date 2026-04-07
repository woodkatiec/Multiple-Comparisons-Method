# Multiple-Comparisons-Method
# My comparison function uses standard t value and critical difference calculation, but adjusts alpha to using the formula: 
      number of pairwise comparisons = m= I*(I-1)/2)
      adjusted alpha = alpha/m^(m/(m+log10(m)))
      # Method:

# To use my function call the function "zewari_comparison" with the following parameters
  1. a vector with the means of each population
  2. the number of observations per population
  3. the mean standard error of the populations
  4. the desired alpha/type one error (will automatically be set to .05 if not given)

# My function will return a true or false matrix comparing the populations
  if the results matrix at [i,j] is FALSE it means no difference was detected between populations I and J
  if the results matrix at [i,j] is TRUE it means a difference was detcted between populations I and J

# Example results for comparing petal lengths of different species of irises
> irisResults
     V1    V2    V3
1 FALSE  TRUE  TRUE
2  TRUE FALSE  TRUE
3  TRUE  TRUE FALSE

# Zewari method is best for when false positives want to be avoided and when there are 8 or more observations per population
   The Zewari method tends to be slightly conservative, check full simulation code for data frame with FWER for I:2-9, J:5-20
   The Zewari method tends to have almost identical powr to Tukey for I:5, J:10, check power table in full simulation code to see if       power for a given I and J is desired for a given scenario 
