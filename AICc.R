# Akaike information criterion for comparing phylogenetic models

AICc <- function(LogL_constrained, LogL_unconstrained, K, n) {
  
  # AIC
  
  AIC_constrained <- -2*LogL_constrained + 2*K
  
  # Correct AIC for small sample size
  
  AICc_constrained <- AIC_constrained + (2*K*(K+1))/(n-K-1)
  
  # AIC
  
  AIC_unconstrained <- -2*LogL_unconstrained + 2*K
  
  # Correct AIC
  
  AICc_unconstrained <- AIC_unconstrained + (2*K*(K+1))/(n-K-1)
  
  # Calculate minimum AIC between models
  
  minAICc <- min(AICc_constrained,AICc_unconstrained)
  
  # Delta AIC for unconstrained
  
  dAICc_unconstrained <- AICc_unconstrained - minAICc
  
  # Delta AIC for constrained
  
  dAICc_constrained <- AICc_constrained - minAICc
  
  # Weights
  
  sum <- exp(-dAICc_unconstrained/2) + exp(-dAICc_constrained/2)
  
  # Weight unconstrained
  
  w_unconstrained <- exp(-dAICc_unconstrained/2)/sum
  
  # Weight constrained
  
  w_constrained <- exp(-dAICc_constrained/2)/sum
  
  Constrained <- c(LogL_constrained,AIC_constrained,AICc_constrained,dAICc_constrained,w_constrained)
  Unconstrained <- c(LogL_unconstrained,AIC_unconstrained,AICc_unconstrained,dAICc_unconstrained,w_unconstrained)
  
  result <- t(data.frame(Constrained,Unconstrained,row.names=c("Log-likelihood","AIC","AICc","dAICc","Weight")))
  
  return(result)
  
}

### Data: Concatenated alignments

# Number of free parameters (df)

K <- 190

# Number of characters in alignment

n <- 484

# Log-likelihood of best constrained tree (Adeno/NCLDVs sister groups)

LogL_constrained <- -41618.186771

# Log-likelihood of best unconstrained tree (Adeno/NCLDVs NOT sister groups)

LogL_unconstrained <- -41612.481141

# Compare models based on the concatenated alignment

AICc(LogL_constrained, LogL_unconstrained, K, n)

### Data: Major capsid protein

# Number of free parameters (df)

K <- 125

# Number of characters in alignment

n <- 150

# Log-likelihood of best constrained tree (Adeno/NCLDVs sister groups)

LogL_constrained <- -16305.521480

# Log-likelihood of best unconstrained tree (Adeno/NCLDVs NOT sister groups)

LogL_unconstrained <- -16304.569687

# Compare models based on the major capsid alignment

AICc(LogL_constrained, LogL_unconstrained, K, n)

# Data: Minor capsid protein

# Number of free parameters (df)

K <- 125

# Number of characters in alignment

n <- 115

# Log-likelihood of best constrained tree (Adeno/NCLDVs sister groups)

LogL_constrained <- -10234.557166

# Log-likelihood of best unconstrained tree (Adeno/NCLDVs NOT sister groups)

LogL_unconstrained <- -10231.749148

# Compare models based on the minor capsid alignment

AICc(LogL_constrained, LogL_unconstrained, K, n)

# Data: Protease

# Number of free parameters (df)

K <- 126

# Number of characters in alignment

n <- 111

# Log-likelihood of best constrained tree (Adeno/NCLDVs sister groups)

LogL_constrained <- -6323.233122

# Log-likelihood of best unconstrained tree (Adeno/NCLDVs NOT sister groups)

LogL_unconstrained <- -6317.549935

# Compare models based on the protease capsid alignment

AICc(LogL_constrained, LogL_unconstrained, K, n)
