# Multiple-Comparisons-Method
# Zewari function uses standard t value and critical difference calculation, but adjusts alpha using the following formula: 
      number of pairwise comparisons = m= I*(I-1)/2)
      adjusted alpha = alpha/m^(m/(m+log10(m)))

# To use the function call "zewari_comparison" with the following parameters:
  1. a vector with the means of each population
  2. the number of observations per population
  3. the mean standard error of the populations
  4. the desired alpha/type one error (will automatically be set to .05 if not given)

# Zewari function will return a true or false matrix comparing the populations
  If the results matrix at [i,j] is FALSE it means no difference was detected between populations I and J
  If the results matrix at [i,j] is TRUE it means a difference was detected between populations I and J

# Example results for comparing petal lengths of different species of irises
> irisResults
     V1    V2    V3
1 FALSE  TRUE  TRUE
2  TRUE FALSE  TRUE
3  TRUE  TRUE FALSE

# Zewari Method Usage and Limitations
The Zewari method tends to be slightly conservative, with lower FWER than alpha. It has almost identical, but slightly less power than Tukey's method. However, difference in power between Tukey and Zewari method is much smaller than the difference in FWER. Thus, the Zewari method should be used when false positives need to be avoided. The ideal scenario for the Zewari method (the scenario with has highest power and lowest FWER combination) is for I=10, J=19. In this scenario, Zewari's power is only .0002 less than Tukey, whereas the error is lowered by .01002. 
   
# Checking Power and FWER
Code for checking Zewari's power and FWER for specific I, J, standard deviation, and alpha is available in repository. There is also a set of simulations already run in the full_simulations file that provides a data frame of power and FWER for I: 2-9, J: 5-20.
