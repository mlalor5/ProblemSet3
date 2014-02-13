#Lalor PS03
#PS Applied Statistical Programming Spring 2014
#Due Feb 13, 2013

## NO FOR LOOPS ##

### PART I: Sampling Distributions and p-values ###

# Load relevant libraries
library(plyr)
library(doMC)
library(abind) #abind combines multi-dim arrays


## 1) Make a three dimensional array with dim=c(20,5, 1000) and fill it with random data. ##

array1 <- array(data=(rnorm(100*1000, mean=0)) , dim=c(20,5,1000))
#dim(array1) 20 x 5 x 1000


## 2) Make function to create a Y values that are a linear combination of X plus normal error ##

#Here is the vector of covariates
Beta <- matrix(c(1,2,0,4,0), ncol=1)

#Function to create Y vals
yfun <- function(x, Beta, Parallel=F) { #Add in Parallel option for Q #7
  .singley <- function(x, Beta) { #Function to make a single y value
    x%*%Beta + rnorm(1) # 1x5 * 5x1 
  }
  yval <- t(aaply(.data=x, .margins=3, .fun=.singley, Beta, .parallel=Parallel)) #for entire array
  return(yval) #return transpose so requested dimensions (below)
}

#The output should be a 20 by 1000 array.
yval <- yfun(array1,Beta)
#dim(yval) #output is 20 x 1000


## 3) Run 1,000 regressions across all of this simulated data. ##

#combine into one array:
simdata <- abind(array1, yval, along=2)
dim(simdata)
data <- simdata
regfun <- function(data, Parallel=F) { #input of 3-d array where 2nd dim x1,...,xn,y
  
  .extract <- function(data) { #Function to get beta estimates
    dims<- dim(data) # 2nd dim numeric will be different variables
    y<- data[,dims[2]] # breaking up to y and x for readability
    x<- data[,1:(dims[2]-1)]  
    datalm <- lm(y~x)
    return(as.numeric(datalm$coefficients)) #6 num, otherwise would also get names like intercept, etc
  }
  output <- aaply(.data=data, .margins=3, .fun=.extract, .parallel=Parallel)
  return(as.matrix(output))
}

#output a 1000 by 6 matrix of estimated regression coefficients.
  test <- regfun(simdata) 
  #dim(test) # 1000 x 6
  #class(test) #matrix


## 4) Create a density plot for each of the 6 coefficients   ##

plot(density(test[,1]), main="Intercept") #Main describes particular coefficient
plot(density(test[,2]), main="Beta 1")
plot(density(test[,3]), main="Beta 2")
plot(density(test[,4]), main="Beta 3")
plot(density(test[,5]), main="Beta 4")
plot(density(test[,6]), main="Beta 5")
 
#What does this density represent? 
# These plots are the distributions of estimated betas from all regressions,
# aka the proportion of regressions that are found a given value e.g. 1 for a given coefficient e.g. intercept


## 5) Alter your code so that you now collect t-statistics for all 1,000 regressions for all six coefficients ##

tstatsfun <- function(data, , Parallel=F) { #renamed for tstats
  .extract <- function(data) { #Function to get t-estimates for betas
    dims<- dim(data) # 2nd dim numeric will be different variables
    y<- data[,dims[2]] # breaking up to y and x for readability
    x<- data[,1:(dims[2]-1)]  
    datalm <- lm(y~x) 
    sum <- summary(datalm) #summary function to save coefficient t-vals
    return(as.numeric(coef(sum)[, "t value"])) #6 num, otherwise would also get names like intercept, etc
  }
  toutput <- aaply(.data=data, .margins=3, .fun=.extract, .parallel=Parallel)
  return(as.matrix(toutput))
}

#run altered function
test <- tstatsfun(simdata) 
dim(test) # 1000 x 6
class(test) #matrix
#summary(test) #note 4 and 6 have mean t-vals clearly < crit, makes sense since for those beta=0


## 6) For the 1,000 regressions, calculate how many t-statistics are statistically “significant” (p≤.05) for each variable ##

# calculate critical t-value for alpha=0.05 (2.1448 here)
criticalval <- qt(0.975, 14)#0.975 for 2-sided t-test, 20 observations (n) with 5 parameters (k)  so 20-5 -1 = 14 df

#Function to count number of crits
critstats <- function(tstats, critval, Parallel=F) { #input stats and critical t-value (assuming 2-sided)
  #test crit for a given regression
  .crit <- function(x, critval){ #checks a given row
    critcoeff <- ifelse(critval >= abs(as.numeric(x)),0,1) # 0 if less than crit, 1 otherwise
    return(critcoeff)
  }
  #Find for all regressions
  critdata <- aaply(.data=tstats, .margins=1, .fun=.crit, critval=criticalval, .parallel=Parallel) #run by row
  summary(critdata) #1s all down when summary(test) has a large enough min aka all but 4 and 6
  #aaply wouldn't recognize sum as a function, so using standard apply
  numbercrit <- apply (critdata, MARGIN=2, FUN=sum)#return how many t-stats were significiant for each coefficient
  return(numbercrit) #return the number of crit for each coeff
}

#Run function
(crits <- critstats(test, criticalval)) #Number of crits for each coefficient, notice beta=0 were ones with fewer crits (aka # of regressions with Type I Error)


## 7) Re-run that code (facebook answer- any) in parallel. Using the system.time command, estimate how much time is saved
system.time(regfun(simdata)) #one run user: 1.428, system: 0.000, elapsed: 1.423


#Added in parallel option so can rerun in parallel
registerDoMC(cores=4)
system.time(regfun(simdata, Parallel=T)) #one run user: 2.820, system: 0.148, elapsed: 1.021

#From system.time estimates it appears running in parallel saved  0.402 (1.423-1.021)

### Part II: Calculating Fit Statistics ###

## 1) Using the Out of step dataset, randomly subset the data into two partitions. 
incumb2 <- read.table(file="~/Documents/Spring_2014/Montgomery/PS03/incumbents2.txt", header=T, sep="\t")
#summary(incumb2$voteshare) #125 NAs in dependent variable - remove these
#dim(incumb2) #6687 x 20
incumbentdata <- na.omit(incumb2) # dim(incumbentdata) now 3193 x 20
index <- 1:dim(incumbentdata)[1] #1 is number of rows


#set seed so results are reproduceable
set.seed(25)
testind <- sample(index, size = (length(index)/2))
testdata <- incumbentdata[testind,] #3281 randomly subsetted
trainingdata <- incumbentdata[-testind,]

# Use one partition (“training set”) to build three statistical models 
#with incumbent vote share as dependent variable

#Model 1 voteshare~ Party of President + Midterm Election + Incumbent Spending + Quality Challenger + Number Unemployed (logged)
model1 <- lm( voteshare~incparty + midterm + incspend + chalquality + unemployed, data=trainingdata)
summary(model1)

#Model 2 - Model 1 - midterm - unemployment + Incumbent Spending Squared + Vote Share of Presidential candidate (same party, last 2 elections)
model2 <- lm( voteshare~incparty + incspend + incspend2 + presvote + chalquality, data=trainingdata)
summary(model2)

#Model 3 - Kitchen Sink (all variables that are meaningful/not redundant)
model3 <- lm(voteshare~.-x-difflog-congress, data=trainingdata)
summary(model3)

# Use these models to make “predictions” for the partition of the data you did
#NOT use to fit the model. This is your “test” set.
m1pred <- predict(model1, newdata=testdata)
m2pred <- predict(model2, newdata=testdata)
m3pred <- predict(model3, newdata=testdata)


## 2) Write a function that takes as arguments 
#(1) a vector of “true” observed outcomes (y), 
#(2) a matrix of predictions (P), and (3) a vector of naive forecasts (r). 
#The matrix should be organized so that each column represents a single forecasting model and 
#the rows correspond with each observation being predicted.
#From Facebook: "The naive model should be whatever you want. A regression with just a constant, or just including one variable."

#The function should output a matrix where each column corresponds with 
# one of the above (PS03) fit statistics, and each row corresponds to a model.

# THIS VERSION WAS IN PREVIOUS COMMIT, Final version incorporates Q 3


## 3) Alter function so user can specify stats and r optional ##

predstats <- function(y, p, r="", stats="") { #r optional, specify stats to run

  #Pre-stat computing math
  #absolute error e_i = |p_i - y_i|
  abserror <- abs(p-cbind(y,y,y)) # technically don't need cbind, but reminder of dim
  #absolute error percentage a_i =e_i / |y_i| * 100
  abserrorpercent <- abserror/abs(y)*100
  #IF R is specified prevent warning by only looking at first element
  if (r[1] != "") { 
    #baseline b_i = |r_i - y_i|
    baseline <- abs(r-y)
  }

  #Compute Requested Stats
  #RMSE = sqrt( sum e^2/n)  
  if ("RMSE" %in% stats){
    RMSE <- sqrt(colSums(e^2)/length(y))    
  } else {
    RMSE =NULL
  }
  
  #MAD = median(e)
  if ("MAD" %in% stats){
    MAD = aaply(e, .margins=2, .fun=median)    
  } else {
    MAD =NULL
  }
  
  
  #RMSLE = sqrt ( sum (ln(p_i+1) - ln(y_i+1))^2 /n )
  if ("RMSLE" %in% stats){
    numerator <- (log(p+1)-log(y+1))^2
    numerator <- colSums(numerator)
    RMSLE <- sqrt(numerator/length(y))    
  } else {
    RMSLE =NULL
  }
  

  #MAPE = sum a / n
  if ("MAPE" %in% stats){
    MAPE <- colSums(abserrorpercent)/length(y)    
  } else {
    MAPE =NULL
  }

  
  #MEAPE = median(a)
  if ("MEAPE" %in% stats){
    MEAPE = aaply(abserrorpercent, .margins=2, .fun=median)  
  } else {
    MEAPE =NULL
  }

  
  #MRAE = median ( e_i / b_i)
  if ("MRAE" %in% stats & r[1] != ""){
    MRAE = aaply((abserror/baseline), .margins=2, .fun=median) 
  } else {
    MRAE =NULL
  }


  #return matrix where  column corresponds with a fit stat, and each row a model
  returnmatrix <- cbind(RMSE, MAD, RMSLE, MAPE, MEAPE,MRAE)
  return(returnmatrix)
}

#Naive model, incumbent voteshare is 50%, or 0.5
r <- rep(0.5, length(m1pred))
# Y from definitions
y <- testdata$voteshare
# Matrix of predictions
p <- as.matrix(cbind(m1pred,m2pred,m3pred))

#Run function - note MRAE requested but not returned since no r specified
(run <- predstats(y,p, stats=c("MAD","MAPE", "MRAE"))) #prints requested stats with no r
# class(run) #returns matrix

# test function
testfunction <- function(y,p,r) { #Will run with and without r, and various stats
  ok = TRUE
  
  #Run with specified r
  evaluate <- predstats(y,p,r, stats=c("RMSE","MAD","RMSLE", "MAPE","MEAPE","MRAE"))
  if (dim(evaluate)[1] != 3 | dim(evaluate)[2] != 6){
    ok = FALSE
    print("PROBLEM WITH FULL SET OF FIT STATISTICS")
  }
  
  #Run without r
  evaluate <- predstats(y,p, stats=c("RMSE","MAD","RMSLE", "MAPE","MEAPE","MRAE"))
  if (dim(evaluate)[1] != 3 | dim(evaluate)[2] != 5){
    ok = FALSE
    print("PROBLEM WITH FULL SET OF FIT STATISTICS, NO R")
  }
  
  #Run with fewer stats
  evaluate <- predstats(y,p,r, stats=c("MRAE"))
  if (dim(evaluate)[1] != 3 | dim(evaluate)[2] != 1){
    ok = FALSE
    print("PROBLEM WITH MRAE, WITH R")
  }
  evaluate <- predstats(y,p, stats=c("MRAE"))
  if (length(evaluate) != 0){ #SHOULD RETURN EMPTY SINCE NO R
    ok = FALSE
    print("PROBLEM WITH NO R, REQUESTING MRAE")
  }
  return(ok) #RETURNS TRUE IF ALL TESTS PASSED, else false
} #END TEST FUNCTION

#RETURNS TRUE SINCE ALL TESTS PASSED
(teststatfun <- testfunction(y,p,r))


## 4) Evaluate the accuracy of models
(evaluate <- predstats(y,p,r, stats=c("RMSE","MAD","RMSLE", "MAPE","MEAPE","MRAE")))

#From output we see that the more variables included, the smaller the fit statistics.
#Thus, model 3 has the smallest values for ALL fit statistics of the three.  Comparing
#models with an equal number of variables aka models 1 and 2, removing the dummy variable 
#for midterms and the number of unemployed, and replacing with the vote share of the presidential
#candidate (of same party) as well as making incumbent spending quadratic, appears to improve
#all fit statistics across the board (see evaluate output).
