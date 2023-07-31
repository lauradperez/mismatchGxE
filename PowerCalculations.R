#!/usr/bin/env Rscript

#===================================================================================
#--Power calculations for gene-by-environment analysis
#===================================================================================

# Set the parameters
sample_size <- 
allele_frequency <- # population specific
effect_size <- 
effect_size_se <- 
disease_prevalence <- # population specific
bazar_food <-
food_secure <- 
sex <- # population specific

# Generate the genotype variable
genotype <- rbinom(sample_size, 2, allele_frequency)

# Generate the effect size for each genotype, incorporating standard error
effect_size <- rep(effect_size, sample_size)
effect_size_se <- rep(effect_size_se, sample_size)
effect_size[genotype == 0] <- 0

# Generate the environmental variables
bazarfood <- rbinom(sample_size, 1, bazar_food)
foodsecure <- rbinom(sample_size, 1, food_secure)
sex <- rbinom(sample_size, 1, sex)
age <- abs(rnorm(sample_size, mean = , sd = )) # population specific 
bmi <- abs(rnorm(sample_size, mean = , sd = )) # population specific

# Calculate the log-odds based on the genotype and adjusted effect size
log_odds <- log(disease_prevalence / (1 - disease_prevalence)) + (effect_size * genotype)

# Calculate the standard error of the log-odds
log_odds_se <- sqrt((effect_size_se^2) * genotype)

# Generate the log-odds with the added standard error
log_odds_with_se <- log_odds + abs(rnorm(sample_size, mean = 0, sd = log_odds_se))

# Calculate the multiplier based on the interaction between genotype and environmental variables
multiplier <- ifelse(genotype == 1 & bazarfood == 1 & foodsecure == 1, XXX,
                     ifelse(genotype == 2 & bazarfood == 1 & foodsecure == 1, XXX,
                            ifelse(genotype == 1 & bazarfood == 0 & foodsecure == 1, XXX,
                                   ifelse(genotype == 1 & bazarfood == 1 & foodsecure == 0, XXX,
                                          ifelse(genotype == 2 & bazarfood == 1 & foodsecure == 0, XXX,
                                                 ifelse(genotype == 2 & bazarfood == 0 & foodsecure == 1, XXX, 1))))))

# Apply the multiplier to the log-odds
log_odds_with_multiplier <- log_odds_with_se * multiplier

# Convert log-odds to probabilities using the logistic function
disease_outcome <- rbinom(sample_size, 1, plogis(log_odds_with_multiplier))

# Create the simulated dataset
simulated_data <- data.frame(
  genotype = genotype,
  effect_size = effect_size,
  effect_size_se = effect_size_se,
  disease_outcome = disease_outcome,
  bazarfood = bazarfood,
  foodsecure = foodsecure,
  sex = sex,
  age = age,
  bmi = bmi
)

# Power analysis
n.it <- 1000
sub.size <- 800
big.power <- matrix(0, nrow = sub.size - 49, ncol = 2)

for (j in 50:sub.size) {
  p.out <- numeric(n.it)
  subs <- sample(nrow(simulated_data), size = n.it * j, replace = TRUE)
  
  for (i in 1:n.it) {
    sub_data <- simulated_data[subs[((i - 1) * j + 1):(i * j)], ]
    glm.null <- glm(disease_outcome ~ genotype + age + sex + bazarfood + foodsecure, 
                    family = "binomial", data = sub_data)  # Null model with intercept only
    glm.full <- glm(disease_outcome ~ genotype*bazarfood*foodsecure + age + sex, 
                    family = "binomial", data = sub_data) # Full model with predictor variable
    glm.p <- anova(glm.null, glm.full, test = "LRT")[["Pr(>Chi)"]][2]
    p.out[i] <- glm.p
  }
  
  sub.pow <- length(which(p.out < 0.05))/n.it
  big.power[j - 49, ] <- c(j, sub.pow)
  print(j)
  
}

colnames(big.power) <- c("Sample_Size", "Power")

big.power <- as.data.frame(big.power)



