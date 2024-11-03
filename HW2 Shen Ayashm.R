##################################################
#    Principles of Statistical Data Analysis     #
#    Home Work 2                                 #
#                                                #
#    Authors:                                    #
#     - Jiacheng Shen                            #
#     - Sajjad Ayashm                            #
#                                                #
##################################################

# Question 1

# function of median.test(x, y , n)

median.test <- function(x, y, n = 1000){
  # Check is the inputs are vector or not.
  stopifnot("the input of function should be a vector" = is.vector(x))
  stopifnot("the input of function should be a vector" = is.vector(y))
  
  # Calculate observed statistic
  observed_stat <- median(x) - median(y)
  
  # Combine data for permutation
  combine <- c(x, y)
  
  extreme_count = 0
  
  # Permutation test loop
  for (i in 1:n){
    shuffle <- sample(combine)
    perm_stat <- median(head(shuffle,length(x))) - median(tail(shuffle,length(y)))
    
    if (abs(perm_stat) >= abs(observed_stat)){
      extreme_count <- extreme_count + 1
    }
  }
  
  # Calculate p-value
  p_value <-  extreme_count / n
  return(p_value)
}



# Question 2

# Load the 'coin' package for permutation tests
library('coin')

# Initialize vectors to store p-values for each test
p.t <- p.wmw <- p.median <- c()

# Set sample size and effect size (delta) for simulations
n <- 20
delta <- sqrt(3)/2

# Run 1,000 simulations to compare test performance
for (i in 1:1000) {
  # Generate two samples from the t-distribution, one shifted by delta
  Y1 <- rt(n, 3)
  Y2 <- rt(n, 3) + delta
  
  # Create a group indicator and combine samples
  X <- factor(c(rep("A", n), rep("B", n)))
  Y <- c(Y1, Y2)
  
  # Compute p-values for each test
  p.t[i] <- pvalue(oneway_test(Y ~ X, distribution = approximate(B = 9999)))
  p.wmw[i] <- wilcox.test(Y1, Y2, exact = TRUE)$p.value
  p.median[i] <- median.test(Y1, Y2, n)
}



mean(p.median < 0.05)
mean(p.t < 0.05)
mean(p.wmw < 0.05)













