# Multiple-Comparisons-Method
# My comparison function adjusts alpha to using the formula: 
      adjusted alpha = alpha/m^(m/(m+log10(m)))

# To use my function call the function "wood_comparison" with the following parameters
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



