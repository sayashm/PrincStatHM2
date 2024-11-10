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



# Question 3
library(extraDistr)
library(coin)

n <- 20
delta <- sqrt(3)/2

distribution_list <- list(
  exp = function(n) rexp(n, rate = 1),
  t3 = function(n) rt(n, df = 3),
  laplace = function(n) rlaplace(n, mu = 0, sigma = 1),
  t5 = function(n) rt(n, df = 5),
  logistic = function(n) rlogis(n, location = 0, scale = 1),
  normal = function(n) rnorm(n, mean = 0, sd = 1),
  uniform = function(n) runif(n, min = 0, max = 1)
)


# Initialize matrices to store p-values for each test and distribution
p_values_t <- matrix(0, nrow = 1000, ncol = length(distribution_list))
p_values_wmw <- matrix(0, nrow = 1000, ncol = length(distribution_list))
p_values_median <- matrix(0, nrow = 1000, ncol = length(distribution_list))
colnames(p_values_t) <- colnames(p_values_wmw) <- colnames(p_values_median) <- names(distribution_list)

# Run the simulations
for (i in 1:1000) {
  for (dist_name in names(distribution_list)) {
    # Generate samples
    Y1 <- distribution_list[[dist_name]](n)
    Y2 <- distribution_list[[dist_name]](n) + delta
    
    # Combine samples and create group indicator
    X <- factor(c(rep("A", n), rep("B", n)))
    Y <- c(Y1, Y2)
    
    # Perform each test and store p-values
    p_values_t[i, dist_name] <- pvalue(oneway_test(Y ~ X, distribution = approximate(B = 9999)))
    p_values_wmw[i, dist_name] <- wilcox.test(Y1, Y2, exact = TRUE)$p.value
    p_values_median[i, dist_name] <- median.test(Y1, Y2, n)
  }
}

# Calculate power (proportion of p-values below alpha) for each test and distribution
power_t <- colMeans(p_values_t < 0.05) * 100
power_wmw <- colMeans(p_values_wmw < 0.05) * 100
power_median <- colMeans(p_values_median < 0.05) * 100

# Combine power results into a data frame for plotting
power_data <- data.frame(
  Distribution = rep(names(distribution_list), 3),
  Power = c(power_t, power_wmw, power_median),
  Test = factor(rep(c("Permutation t-test", "WMW test", "Median test"), each = length(distribution_list)))
)

# Plot the results with horizontal bars
library(ggplot2)
ggplot(power_data, aes(x = Distribution, y = Power, fill = Test)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Power of Different Tests across Distributions",
       x = "Distribution",
       y = "Power (%)") +
  scale_fill_manual(values = c("gray4", "gray37", "gray83")) +
  theme_minimal() +
  coord_flip() # This flips the axes, making the bars horizontal

# Question 4

library(extraDistr)
library(coin)

n <- 20
delta_new <- 0

distribution_list_new <- list(
  exp = function(n) rexp(n, rate = 1),
  t3 = function(n) rt(n, df = 3),
  laplace = function(n) rlaplace(n, mu = 0, sigma = 1),
  t5 = function(n) rt(n, df = 5),
  logistic = function(n) rlogis(n, location = 0, scale = 1),
  normal = function(n) rnorm(n, mean = 0, sd = 1),
  uniform = function(n) runif(n, min = 0, max = 1)
)


# Initialize new matrices to store p-values for each test and distribution
p_values_t_new <- matrix(0, nrow = 1000, ncol = length(distribution_list))
p_values_wmw_new <- matrix(0, nrow = 1000, ncol = length(distribution_list))
p_values_median_new <- matrix(0, nrow = 1000, ncol = length(distribution_list))
colnames(p_values_t_new) <- colnames(p_values_wmw_new) <- colnames(p_values_median_new) <- names(distribution_list_new)

# Run the simulations
for (i in 1:1000) {
  for (dist_name in names(distribution_list_new)) {
    # Generate samples
    Y1 <- distribution_list_new[[dist_name]](n)
    Y2 <- distribution_list_new[[dist_name]](n) + delta_new
    
    # Combine samples and create group indicator
    X <- factor(c(rep("A", n), rep("B", n)))
    Y <- c(Y1, Y2)
    
    # Perform each test and store p-values
    p_values_t_new[i, dist_name] <- pvalue(oneway_test(Y ~ X, distribution = approximate(B = 9999)))
    p_values_wmw_new[i, dist_name] <- wilcox.test(Y1, Y2, exact = TRUE)$p.value
    p_values_median_new[i, dist_name] <- median.test(Y1, Y2, n)
  }
}

# Calculate power (proportion of p-values below alpha) for each test and distribution
power_t_new <- colMeans(p_values_t_new < 0.05) * 100
power_wmw_new <- colMeans(p_values_wmw_new < 0.05) * 100
power_median_new <- colMeans(p_values_median_new < 0.05) * 100

# Combine power results into a data frame for plotting
power_data_new <- data.frame(
  Distribution_new = rep(names(distribution_list_new), 3),
  Power_new = c(power_t_new, power_wmw_new, power_median_new),
  Test_new = factor(rep(c("Permutation t-test", "WMW test", "Median test"), each = length(distribution_list_new)))
)

# Plot the results with horizontal bars
library(ggplot2)
ggplot(power_data_new, aes(x = Distribution_new, y = Power_new, fill = Test_new)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Type I Errors across Distributions",
       x = "Distribution",
       y = "Type I Error (%)") +
  scale_fill_manual(values = c("gray4", "gray37", "gray83")) +
  theme_minimal() +
  coord_flip() # This flips the axes, making the bars horizontal
