# Multiple-Comparisons-Method
# Zewari function uses standard t value and critical difference calculation, but adjusts alpha using the following formula: 
      number of pairwise comparisons = m= I*(I-1)/2
      adjusted alpha = alpha/m^(m/(m+log10(m)))

# To use the function call "zewari_comparison" with the following parameters:
  1. a vector with the means of each population
  2. the number of observations per population
  3. the pooled MSE of the populations
  4. the desired alpha/Type I error (will automatically be set to .05 if not given)

# Zewari function will return a TRUE/FALSE matrix comparing the populations
  If the results matrix at [i,j] is FALSE it means no difference was detected between populations i and j
  If the results matrix at [i,j] is TRUE it means a difference was detected between populations i and j

# Example results for comparing petal lengths of different species of irises
The Zewari method identified a difference between each of the different populations, which is accurate because iris petal length does depends on the species. The full example code is found in zewari_comparison_example

# Zewari Method Usage and Limitations
The Zewari method tends to be slightly conservative, and generally has lower FWER than alpha. It has almost identical, but slightly less power than Tukey's method. However, difference in power between Tukey and Zewari method is much smaller than the difference in FWER. Thus, the Zewari method is generally preferable and should be used when false positives need to be avoided. However, Tukey should be used when it is crucial for power to be absolutely maximized.
Considering I:2-9 and J:5-20, the scenario in which the Zewari method produces the best combination of power-FWER at alpha = .05, delta = 1.5, is when I=10 and J=19. In this scenario, Zewari's power is only .0002 less than Tukey, whereas the error is lowered by .01002. 
The scenario in which the Zewari method produces the highest power/FWER ratio at alpha = .05, delta = 1.5 is when I=9 and J=17. In this scenario, the Zewari method has a 25.8:1 Power:FWER ratio, whereas Tukey has a 19.9:1 Power:FWER ratio.
   
# Checking Power and FWER
Code for checking Zewari's power and FWER for specific I, J, standard deviation, and alpha is available in the file "Test FWER and Power for specific data." There is also a set of simulations already run in the file "Full_Simulations_Zewari.R" that provides data frames of power and FWER for I: 2-9, J: 5-20.
