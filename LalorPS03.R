#Lalor PS03
#PS Applied Statistical Programming Spring 2014
#Due Feb 13, 2013

## NO FOR LOOPS ##
### PART I: Sampling Distributions and p-values ###

# Load relevant libraries
library(plyr)

## 1) Make a three dimensional array with dim=c(20,5, 1000) and fill it with random data.

array1 <- array(data=(rnorm(100*1000, mean=2)) , dim=c(20,5,1000))
#dim(array1) 20 x 5 x 1000
## 2) Make function to create a Y values that are a linear combination of X pplus normal error

#Here is the vector of covariates
Beta <- matrix(c(1,2,0,4,0), ncol=1)

#Function to create Y vals
yfun <- function(x, Beta) {
  .singley <- function(x, Beta) { #Function to make a single y value
    x%*%Beta + rnorm(1) # 1x5 * 5x1 
  }
  yval <- t(aaply(.data=x, .margins=3, .fun=.singley, Beta)) #for entire array
  return(yval) #return transpose so requested dimensions (below)
}

#The output should be a 20 by 1000 array.
yval <- yfun(array1,Beta)
#dim(yval) #output is 20 x 1000

## 3) Run 1,000 regressions across all of this simulated data. 
#output a 1000 by 6 matrix of estimated regression coefficients.

#combine into one array:
library(abind) #abind combines multi-dim arrays
simdata <- abind(array1, yval, along=2)
dim(simdata)
data <- simdata
regfun <- function(data) { #input of 3-d array where 2nd dim x1,...,xn,y
  
  .extract <- function(data) { #Function to get beta estimates
    dims<- dim(data) # 2nd dim numeric will be different variables
    y<- data[,dims[2]] # breaking up to y and x for readability
    x<- data[,1:(dims[2]-1)]  
    datalm <- lm(y~x)
    return(as.numeric(datalm$coefficients)) #6 num, otherwise would also get names like intercept, etc
  }
  output <- aaply(.data=data, .margins=3, .fun=.extract)
  return(as.matrix(output))
}

#try out function to run regressions
  test <- regfun(simdata) 
  dim(test) # 1000 x 6
  class(test) #matrix

## 4) Create a density plot for each of the 6 coefficients 
 #What does this density represent? 

plot(density(test[,1]), main="Intercept") #Main describes particular coefficient
plot(density(test[,2]), main="Beta 1")
plot(density(test[,3]), main="Beta 2")
plot(density(test[,4]), main="Beta 3")
plot(density(test[,5]), main="Beta 4")
plot(density(test[,6]), main="Beta 5")
 
#Note: These plots are the distributions of estimated betas from all regressions,
# aka the proportion of regressions that found a given value e.g. 1 for a given coefficient e.g. intercept



 
