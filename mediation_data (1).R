library(MASS)
library(extraDistr)
library(R2jags)

setwd("C:/Users/au199986/OneDrive - Aarhus Universitet/Courses/F24Thesis/Julie")

#---- Generate mediated constructs using standard Sobel test ------------------#

# number of simulated subjects
nsubs <- 100

#---- Define parameters of Sobel test ----------------

# Intercepts for models
gamma_M1 <- 1
gamma_M2 <- 1
gamma_M3 <- .4
gamma_total <- -2

# Coefficient - X -> M path
alpha_M1 <- 2 
alpha_M2 <- 3
alpha_M3 <- 4

# Coefficient - M -> Y path
beta_M1 <- .3
beta_M2 <- .1
beta_M3 <- .7

# Coefficient - Direct path
tau_prime <- .01

#Sobel test - magnitude of indirect effect for mediators
psy_M1 <- alpha_M1 * beta_M1 
psy_M2 <- alpha_M2 * beta_M2 
psy_M3 <- alpha_M3 * beta_M3 

#SD of mediators - (this will be doing the work of the error term "epsilon" in frequent frameworks)
sigma_M1 <- .5 
sigma_M2 <- .5
sigma_M3 <- .5
sigma_Y <- .1

#------------Arrays for simulation-------------------
X <- rnorm(nsubs,.4,.1) #here we sample 1000 subjects from a normal dist. mean = 3, sd = 1...
#...that looks like this:
par(mfrow=c(1,2))
plot(density(X)) #.... or this
hist(X)

mu_M1 <- array(NA,nsubs) # we need the "mu" arrays, for the Bayesian way of doing regression (error)
M1 <- array(NA,nsubs)

mu_M2 <- array(NA,nsubs)
M2 <- array(NA,nsubs)

mu_M3 <- array(NA,nsubs)
M3 <- array(NA,nsubs)

mu_Y <- array(NA,nsubs)
Y <- array(NA,nsubs)

#----------- Mediation models -----------------------
for (n in 1:nsubs) {
  
  # ------ Models for mediators  --------------------
  mu_M1[n] <- gamma_M1 + (alpha_M1*X[n]) # model 2 in wiki for Sobel test
  M1[n] <- rnorm(1,mu_M1[n],sigma_M1) # this adds some noise to the MEASURE - measurement noise - sigma from above!!
  
  mu_M2[n] <- gamma_M2 + (alpha_M2*X[n])
  M2[n] <- rnorm(1,mu_M2[n],sigma_M2) # this adds some noise to the MEASURE - measurement noise - sigma from above!!

  mu_M3[n] <- gamma_M3 + (alpha_M3*X[n])
  M3[n] <- rnorm(1,mu_M3[n],sigma_M3) # this adds some noise to the MEASURE - measurement noise - sigma from above!!
    
  # ---------- Total model -----------
  # Model 3 in wiki Sobel test
  mu_Y[n] <- gamma_total + 
    (tau_prime*X[n]) +
    (beta_M1*M1[n]) + 
    (beta_M2*M2[n]) +
    (beta_M3*M3[n])    
  Y[n] <- rnorm(1,mu_Y[n],sigma_Y) # this adds some noise to the MEASURE - measurement noise - sigma from above!!
  
}

# how Y looks
par(mfrow=c(1,2))
plot(density(Y)) #.... or this
hist(Y)

par(mfrow=c(3,1))
plot(density(M1))
plot(density(M2))
plot(density(M3))

#-------------------------------------------------------------
#---------- Test Construct Level Simulation ------------------
#---------- Visualise correlations in paths-------------------

par(mfrow=c(3,3))
plot(X,M1)
plot(X,M2)
plot(X,M3)
plot(M1,Y)
plot(M2,Y)
plot(M3,Y)
plot(X,Y)

#----------- test if LM recovers paths 

# X -> M paths - to check alpha
print(gamma_M1)
print(alpha_M1)
summary(lm(M1~X))

print(gamma_M2)
print(alpha_M2)
summary(lm(M2~X))

print(gamma_M3)
print(alpha_M3)
summary(lm(M3~X))

#rest of the model - to check tau prime and beta
print(gamma_total)
print(tau_prime)
print(beta_M1)
print(beta_M2)
print(beta_M3)
summary(lm(Y~X+M1+M2+M3))

#-------------------------------------------------------------------------------
#------------ Generate Item Level responses from simulated constructs ----------
#------------ Model is Kruschke chapter on ordinal probit regression -----------
#-------------------------------------------------------------------------------

# Survey X questions and answers
nq_X <- 1
nr_X <- 11 #WE ARE COLLAPSING THE 0 < 10% RESPONSES HERE, so there is one category less  

# Survey M1 questions and answers
nq_M1 <- 1
nr_M1 <- 5

# Survey M2 questions and answers
nq_M2 <- 1
nr_M2 <- 5

# Survey M2 questions and answers
nq_M3 <- 1
nr_M3 <- 5

# Survey Y questions and answers
nq_Y <- 1 
nr_Y <- 11 # number of proportion options

#---------- simulate survey X --------------------------------------------------
#---------- CONSTRUCT IN % BUT DATA IN ORDINAL CATEGORIES 1-12 -----------------
sigma_theta_X <- .01 #noise in threshold distances - must be less

#threshold array - one fewer threshold than response options
theta_X <- array(NA,c(nsubs,nq_X,nr_X-1)) 
#arrays for probabilities for response options.
p_X <- array(NA,c(nsubs,nq_X,nr_X))
surveyX <- array(NA,c(nsubs,nq_X))

#-Special for response options-
# response can be anything from 0 to 100, in approximate intervals of %10
# here we will represent the data as proportions, and then have thresholds
# that match proportion intervals, rather than a Likert scale
mu_theta_X <- c(.05,.15,.25,.35,.45,.55,.65,.75,.85,.95)

for (s in 1:nsubs) {
  
  for (q in 1:nq_X) { 
    
    # simulate thresholds/cutoffs for responses
    for (t in 1:(nr_X-1)) { # define one fewer threshold than n questions  
      # find thresholds on response scale
      theta_X[s,q,t] <- rnorm(1,mu_theta_X[t],sigma_theta_X) #note we have the mu from above 
    }
    
    # generate probabilities for response categories from thresholded model
    # See Kruschke Chapter 21 for model
    # mean is underlying construct from mvnorm simulation above - i.e X
    p_X[s,q,1] <- pnorm((theta_X[s,q,1]-X[s])/sigma_theta_X)
    p_X[s,q,2] <- pnorm((theta_X[s,q,2]-X[s])/sigma_theta_X)-pnorm((theta_X[s,q,1]-X[s])/sigma_theta_X)
    p_X[s,q,3] <- pnorm((theta_X[s,q,3]-X[s])/sigma_theta_X)-pnorm((theta_X[s,q,2]-X[s])/sigma_theta_X)
    p_X[s,q,4] <- pnorm((theta_X[s,q,4]-X[s])/sigma_theta_X)-pnorm((theta_X[s,q,3]-X[s])/sigma_theta_X)
    p_X[s,q,5] <- pnorm((theta_X[s,q,5]-X[s])/sigma_theta_X)-pnorm((theta_X[s,q,4]-X[s])/sigma_theta_X)
    p_X[s,q,6] <- pnorm((theta_X[s,q,6]-X[s])/sigma_theta_X)-pnorm((theta_X[s,q,5]-X[s])/sigma_theta_X)
    p_X[s,q,7] <- pnorm((theta_X[s,q,7]-X[s])/sigma_theta_X)-pnorm((theta_X[s,q,6]-X[s])/sigma_theta_X)
    p_X[s,q,8] <- pnorm((theta_X[s,q,8]-X[s])/sigma_theta_X)-pnorm((theta_X[s,q,7]-X[s])/sigma_theta_X)
    p_X[s,q,9] <- pnorm((theta_X[s,q,9]-X[s])/sigma_theta_X)-pnorm((theta_X[s,q,8]-X[s])/sigma_theta_X)
    p_X[s,q,10] <- pnorm((theta_X[s,q,10]-X[s])/sigma_theta_X)-pnorm((theta_X[s,q,9]-X[s])/sigma_theta_X)
    p_X[s,q,11] <- 1-pnorm((theta_X[s,q,10]-X[s])/sigma_theta_X)
    
    # simulate response on item using probabilities
    surveyX[s,q] <- rcat(1,p_X[s,q,])
  }
}

#---------- simulate survey M1 --------------------------------------------------
sigma_theta_M1 <- .1 #noise in threshold distances

#threshold array - one fewer threshold than response options
theta_M1 <- array(NA,c(nsubs,nq_M1,nr_M1-1)) 
#arrays for probabilities for response options.
p_M1 <- array(NA,c(nsubs,nq_M1,nr_M1))
surveyM1 <- array(NA,c(nsubs,nq_M1))

for (s in 1:nsubs) {
  
  for (q in 1:nq_M1) { 
    
    # simulate thresholds/cutoffs for responses
    for (t in 1:(nr_M1-1)) { # define one fewer threshold than n responses  
      # find thresholds on response scale
      theta_M1[s,q,t] <- rnorm(1,(t+.5),sigma_theta_M1) 
    }
    
    # generate probabilities for response categories from thresholded model
    # See Kruschke Chapter 21 for model
    # mean is underlying construct from mvnorm simulation above - i.e M1
    p_M1[s,q,1] <- pnorm((theta_M1[s,q,1]-M1[s])/sigma_theta_M1)
    p_M1[s,q,2] <- pnorm((theta_M1[s,q,2]-M1[s])/sigma_theta_M1)-pnorm((theta_M1[s,q,1]-M1[s])/sigma_theta_M1)
    p_M1[s,q,3] <- pnorm((theta_M1[s,q,3]-M1[s])/sigma_theta_M1)-pnorm((theta_M1[s,q,2]-M1[s])/sigma_theta_M1)
    p_M1[s,q,4] <- pnorm((theta_M1[s,q,4]-M1[s])/sigma_theta_M1)-pnorm((theta_M1[s,q,3]-M1[s])/sigma_theta_M1)
    p_M1[s,q,5] <- 1-pnorm((theta_M1[s,q,4]-M1[s])/sigma_theta_M1)
    
    # simulate response on item using probabilities
    surveyM1[s,q] <- rcat(1,p_M1[s,q,])
  }
}

#---------- simulate survey M2 --------------------------------------------------
sigma_theta_M2 <- .1 #noise in threshold distances

#threshold array - one fewer threshold than response options
theta_M2 <- array(NA,c(nsubs,nq_M2,nr_M2-1)) 
#arrays for probabilities for response options.
p_M2 <- array(NA,c(nsubs,nq_M2,nr_M2))
surveyM2 <- array(NA,c(nsubs,nq_M2))

for (s in 1:nsubs) {
  
  for (q in 1:nq_M2) { 
    
    # simulate thresholds/cutoffs for responses
    for (t in 1:(nr_M2-1)) { # define one fewer threshold than n questions  
      # find thresholds on response scale
      theta_M2[s,q,t] <- rnorm(1,(t+.5),sigma_theta_M2) 
    }
    
    # generate probabilities for response categories from thresholded model
    # See Kruschke Chapter 21 for model
    # mean is underlying construct from mvnorm simulation above - i.e M2
    p_M2[s,q,1] <- pnorm((theta_M2[s,q,1]-M2[s])/sigma_theta_M2)
    p_M2[s,q,2] <- pnorm((theta_M2[s,q,2]-M2[s])/sigma_theta_M2)-pnorm((theta_M2[s,q,1]-M2[s])/sigma_theta_M2)
    p_M2[s,q,3] <- pnorm((theta_M2[s,q,3]-M2[s])/sigma_theta_M2)-pnorm((theta_M2[s,q,2]-M2[s])/sigma_theta_M2)
    p_M2[s,q,4] <- pnorm((theta_M2[s,q,4]-M2[s])/sigma_theta_M2)-pnorm((theta_M2[s,q,3]-M2[s])/sigma_theta_M2)
    p_M2[s,q,5] <- 1-pnorm((theta_M2[s,q,4]-M2[s])/sigma_theta_M2)
    
    # simulate response on item using probabilities
    surveyM2[s,q] <- rcat(1,p_M2[s,q,])
  }
}

#---------- simulate survey M3 --------------------------------------------------
sigma_theta_M3 <- .1 #noise in threshold distances

#threshold array - one fewer threshold than response options
theta_M3 <- array(NA,c(nsubs,nq_M3,nr_M3-1)) 
#arrays for probabilities for response options.
p_M3 <- array(NA,c(nsubs,nq_M3,nr_M3))
surveyM3 <- array(NA,c(nsubs,nq_M3))

for (s in 1:nsubs) {
  
  for (q in 1:nq_M3) { 
    
    # simulate thresholds/cutoffs for responses
    for (t in 1:(nr_M3-1)) { # define one fewer threshold than n questions  
      # find thresholds on response scale
      theta_M3[s,q,t] <- rnorm(1,(t+.5),sigma_theta_M3) 
    }
    
    # generate probabilities for response categories from thresholded model
    # See Kruschke Chapter 21 for model
    # mean is underlying construct from mvnorm simulation above - i.e M3
    p_M3[s,q,1] <- pnorm((theta_M3[s,q,1]-M3[s])/sigma_theta_M3)
    p_M3[s,q,2] <- pnorm((theta_M3[s,q,2]-M3[s])/sigma_theta_M3)-pnorm((theta_M3[s,q,1]-M3[s])/sigma_theta_M3)
    p_M3[s,q,3] <- pnorm((theta_M3[s,q,3]-M3[s])/sigma_theta_M3)-pnorm((theta_M3[s,q,2]-M3[s])/sigma_theta_M3)
    p_M3[s,q,4] <- pnorm((theta_M3[s,q,4]-M3[s])/sigma_theta_M3)-pnorm((theta_M3[s,q,3]-M3[s])/sigma_theta_M3)
    p_M3[s,q,5] <- 1-pnorm((theta_M3[s,q,4]-M3[s])/sigma_theta_M3)
    
    # simulate response on item using probabilities
    surveyM3[s,q] <- rcat(1,p_M3[s,q,])
  }
}


#---------- simulate survey Y --------------------------------------------------
#---------- CONSTRUCT IN % BUT DATA IN ORDINAL CATEGORIES 1-12 -----------------

sigma_theta_Y <- .01 #noise in threshold distances - must be less

#threshold array - one fewer threshold than response options
theta_Y <- array(NA,c(nsubs,nq_Y,nr_Y-1)) 
#arrays for probabilities for response options.
p_Y <- array(NA,c(nsubs,nq_Y,nr_Y))
surveyY <- array(NA,c(nsubs,nq_Y))

#-Special for response options-
# response can be anything from < 50% to >50%, in intervals of %10
# here we will represent the data as proportions, and then have thresholds
# that match proportion intervals, rather than a Likert scale
mu_theta_Y <- seq(-.45,.45,.1)

for (s in 1:nsubs) {
  
  for (q in 1:nq_Y) { 
    
    # simulate thresholds/cutoffs for responses
    for (t in 1:(nr_Y-1)) { # define one fewer threshold than n questions  
      # find thresholds on response scale
      theta_Y[s,q,t] <- rnorm(1,mu_theta_Y[t],sigma_theta_Y) #note we have the mu from above 
    }
    
    # generate probabilities for response categories from thresholded model
    # See Kruschke Chapter 21 for model
    # mean is underlying construct from mvnorm simulation above - i.e Y
    p_Y[s,q,1] <- pnorm((theta_Y[s,q,1]-Y[s])/sigma_theta_Y)
    p_Y[s,q,2] <- pnorm((theta_Y[s,q,2]-Y[s])/sigma_theta_Y)-pnorm((theta_Y[s,q,1]-Y[s])/sigma_theta_Y)
    p_Y[s,q,3] <- pnorm((theta_Y[s,q,3]-Y[s])/sigma_theta_Y)-pnorm((theta_Y[s,q,2]-Y[s])/sigma_theta_Y)
    p_Y[s,q,4] <- pnorm((theta_Y[s,q,4]-Y[s])/sigma_theta_Y)-pnorm((theta_Y[s,q,3]-Y[s])/sigma_theta_Y)
    p_Y[s,q,5] <- pnorm((theta_Y[s,q,5]-Y[s])/sigma_theta_Y)-pnorm((theta_Y[s,q,4]-Y[s])/sigma_theta_Y)
    p_Y[s,q,6] <- pnorm((theta_Y[s,q,6]-Y[s])/sigma_theta_Y)-pnorm((theta_Y[s,q,5]-Y[s])/sigma_theta_Y)
    p_Y[s,q,7] <- pnorm((theta_Y[s,q,7]-Y[s])/sigma_theta_Y)-pnorm((theta_Y[s,q,6]-Y[s])/sigma_theta_Y)
    p_Y[s,q,8] <- pnorm((theta_Y[s,q,8]-Y[s])/sigma_theta_Y)-pnorm((theta_Y[s,q,7]-Y[s])/sigma_theta_Y)
    p_Y[s,q,9] <- pnorm((theta_Y[s,q,9]-Y[s])/sigma_theta_Y)-pnorm((theta_Y[s,q,8]-Y[s])/sigma_theta_Y)
    p_Y[s,q,10] <- pnorm((theta_Y[s,q,10]-Y[s])/sigma_theta_Y)-pnorm((theta_Y[s,q,9]-Y[s])/sigma_theta_Y)
    p_Y[s,q,11] <- 1-pnorm((theta_Y[s,q,10]-Y[s])/sigma_theta_Y)
    
    # simulate response on item using probabilities
    surveyY[s,q] <- rcat(1,p_Y[s,q,])
  }
}

#----------------------------------- TEST --------------------------------
#---- recover relation between modeled construct and simulated survey responses
#---- lm result should be intercept = 0, slope = 1
par(mfrow=c(3,2))

#---- survey X ------------
measure <- mu_theta_X[surveyX]
construct <- X

plot(construct,measure,#xlim=c(0,6),ylim=c(0,6),
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
     pch=16,cex=2,
     abline(lm(measure ~ construct)))

print(summary(lm(measure ~ construct))) #intercept should be 0, slope should be 1
print(cor(measure,construct)) #should get a good correlation

#---- survey M1 ------------
measure <- surveyM1
construct <- M1

plot(construct,rowMeans(measure),xlim=c(0,6),ylim=c(0,6),
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
     pch=16,cex=2,
     abline(lm(rowMeans(measure) ~ construct)))

print(summary(lm(rowMeans(measure) ~ construct)))
print(cor(rowMeans(measure),construct))

#---- survey M2 ------------
measure <- surveyM2
construct <- M2

plot(construct,rowMeans(measure),xlim=c(0,6),ylim=c(0,6),
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
     pch=16,cex=2,
     abline(lm(rowMeans(measure) ~ construct)))

print(summary(lm(rowMeans(measure) ~ construct)))
print(cor(rowMeans(measure),construct))

#---- survey M3 ------------
measure <- surveyM3
construct <- M3

plot(construct,rowMeans(measure),#xlim=c(0,6),ylim=c(0,6),
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
     pch=16,cex=2,
     abline(lm(rowMeans(measure) ~ construct)))

print(summary(lm(rowMeans(measure) ~ construct)))
print(cor(rowMeans(measure),construct))

#---- survey Y ------------
measure <- mu_theta_Y[surveyY]
construct <- Y

plot(construct,measure,#xlim=c(0,6),ylim=c(0,6),
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
     pch=16,cex=2,
     abline(lm(measure ~ construct)))

print(summary(lm(measure ~ construct)))
print(cor(measure,construct))

hist(Y)

#----------------------------------- TEST --------------------------------
#----------- test if LM recovers paths from generated survey responses ---

# Now we are doing a Sobel test, as if the "survey responses" were our only "data"

# X -> M paths - to check alpha
print(gamma_M1) # IS 1, because we defined it right at the start
print(alpha_M1) # IS .7, because we defined it right at the start
summary(lm(rowMeans(surveyM1)~rowMeans(surveyX)))

print(gamma_M2) # IS 1, because we defined it right at the start
print(alpha_M2) # IS .2, because we defined it right at the start
summary(lm(rowMeans(surveyM2)~rowMeans(surveyX)))

print(gamma_M3) # IS 1, because we defined it right at the start
print(alpha_M3) # IS .2, because we defined it right at the start
summary(lm(rowMeans(surveyM3)~rowMeans(surveyX)))

#rest of the model - to check tau prime and beta
print(gamma_total) # IS 1, because we defined it right at the start
print(tau_prime) # IS .3, because we defined it right at the start 
print(beta_M1) # IS .7, because we defined it right at the start
print(beta_M2) # IS .2, because we defined it right at the start
print(beta_M3) # IS .2, because we defined it right at the start
summary(lm(Y~X+M1+M2+M3))

#---------------- Apply jags model ---------------------------
#----- effects of self-efficacy on anxiety ----------------------
data <- list("nsubs","nq_X","nr_X","nq_M1","nr_M1","nq_M2","nr_M2",
             "nq_M3","nr_M3","nq_Y","nr_Y",
             "X","M1","M2","M3","Y")
params <- c("gamma_M1","gamma_M2","gamma_M3","gamma_total",
            "alpha_M1","alpha_M2","alpha_M3",
            "beta_M1","beta_M2","beta_M3",
            "tau_prime",
            "psy_M1","psy_M2","psy_M3")


# need to initialise survey thresholds
initList <- function(){
  list(
    theta_X=array(c(rep(.05,nsubs*nq_X),
                      rep(.15,nsubs*nq_X),
                      rep(.25,nsubs*nq_X),
                      rep(.35,nsubs*nq_X),
                      rep(.45,nsubs*nq_X),
                      rep(.55,nsubs*nq_X),
                      rep(.65,nsubs*nq_X),
                      rep(.75,nsubs*nq_X),
                      rep(.85,nsubs*nq_X),
                      rep(.95,nsubs*nq_X)),
                    c(nsubs,nq_X,nr_X-1)),
    
    theta_M1=array(c(rep(1.5,nsubs*nq_M1),
                      rep(2.5,nsubs*nq_M1),
                      rep(3.5,nsubs*nq_M1),
                      rep(4.5,nsubs*nq_M1)),
                    c(nsubs,nq_M1,nr_M1-1)),
    
    theta_M2=array(c(rep(1.5,nsubs*nq_M2),
                     rep(2.5,nsubs*nq_M2),
                     rep(3.5,nsubs*nq_M2),
                     rep(4.5,nsubs*nq_M2)),
                   c(nsubs,nq_M2,nr_M2-1)),

    theta_M3=array(c(rep(1.5,nsubs*nq_M3),
                     rep(2.5,nsubs*nq_M3),
                     rep(3.5,nsubs*nq_M3),
                     rep(4.5,nsubs*nq_M3)),
                   c(nsubs,nq_M3,nr_M3-1)),    
        
    theta_Y=array(c(rep(-.45,nsubs*nq_Y),
                      rep(-.35,nsubs*nq_Y),
                      rep(-.25,nsubs*nq_Y),
                      rep(-.15,nsubs*nq_Y),
                      rep(-.05,nsubs*nq_Y),
                      rep(.05,nsubs*nq_Y),
                      rep(.15,nsubs*nq_Y),
                      rep(.25,nsubs*nq_Y),
                      rep(.35,nsubs*nq_Y),
                      rep(4.5,nsubs*nq_Y)),
                    c(nsubs,nq_Y,nr_Y-1))
  )
}
initList()


samples <- jags(data, inits=initList, params,
                model.file ="jags_model_solution.txt",
                n.chains=3, n.iter=10000, n.burnin=2000, n.thin=1)
