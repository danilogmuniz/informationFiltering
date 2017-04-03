#================================================================
library(truncnorm)
#===============================================================
#Auxiliary functions used during the simulations
#===============================================================
#Function that calculates euclidean distances. It receives two 
#matrices of x and y coordinates as parameters.
euclid = function(m1, m2=m1)
{
  sqrt(outer(m1[,1], m2[,1],"-")^2 + outer(m1[,2], m2[,2],"-")^2)
}
#===============================================================
#Function that takes samples from a vector of probabilities
#Given a vector of length n, it samples values from 1 to n, the 
#probability of value i (i <=n )being sampled is n[i].
sampleProb = function(prob, n=1, replace=FALSE)
{
  prob[is.na(prob)] = 0
  positiveProbs = sum(prob!=0)
  if (positiveProbs == 1)
    return(which(prob!=0))
  else if(positiveProbs > 1)
    return(sample(1:length(prob), size = n, replace=replace, prob=prob))
  else return(NA)
}
#===============================================================
#Function that runs a simulation of a mating system with information
#filtering during mate choice.
#It simulates one mating season of a population with 1:1 sex ratio 
#It returns to the user a list containing three data.frames, if a 
#filename is given, all data is also registered in csv files, one 
#for males, one for females and one for copulations (the groups file).

simScramble = function(N=100, w = 1, zmean=4, zsd=1, radius=0.10, 
                       B=2, filename=NA, append=TRUE, seed=NA,
                       choice=1, zmin=1, zmax=7)
{  

  #if a random seed is given, it must be used
  if(!is.na(seed))
    set.seed(seed)  
  
  #a ttick to facilitate indexing
  xy = c("x","y")
  
  #creating the data.frames for males and females
  females = data.frame(wid=seed, id=1:N, x=runif(N, min=0, max=w),
                       y=runif(N, min=0, max=w))
  males = data.frame(wid=seed, id=1:N, x=runif(N, min=0, max=w), y=runif(N, min=0, max=w), ms=0)
  
  zmatrix = NA
  difs = NA
  
  #creating data.frame for males
  zmean = as.numeric(zmean);zsd = as.numeric(zsd);
  zmin = as.numeric(zmin);zmax = as.numeric(zmax);
  males$trait=rtruncnorm(n = N, a = zmin, b = zmax, mean=zmean, sd=zsd)  
  
  if(choice==1)#directional choice
  {
    #male trait matrix
    zbase = matrix(males$trait, nrow=N, ncol=N, byrow=TRUE)
    zmatrix = zbase^B
    
  }
  else #assortative choice
  {
    #then females also have a trait value(from the same distribution) 
    females$trait=rtruncnorm(n = N, a = zmin, b = zmax, mean=zmean, sd=zsd)
    
    #calculates all pairwise differences
    difs = outer(females$trait, males$trait, "-")
    
    #divides by female trait
    fzmatrix = matrix(females$trait, nrow=N, ncol=N, byrow=TRUE)
    difs = difs/fzmatrix
    difs = abs(difs)
    
    zbase = difs
    zmatrix = exp(-B*(difs))
  }  
  
  #In all the matrices, females are lines, males are columns
  #distance matrix
  smatrix = euclid(females[,xy], males[,xy])
  
  #adjacency matrix of neighborhood (males available)
  #females are connected to males within r distance
  nmatrix = ifelse(smatrix<=radius, 1, 0)
  
  #also, by default females are connected to the closest male
  nmatrix[cbind(1:N, apply(smatrix, 1, which.min))] = 1 
  
  #choice matrix
  cmatrix = nmatrix * zmatrix
  
  #choosing the males
  chosenMales = apply(cmatrix, 1, sampleProb, replace=FALSE)
  
  groups = data.frame(wid=seed, female = females$id, male = chosenMales)
  
  groups$z = zbase[cbind(groups$female,groups$male)]
  
  #adding mating success to the males
  males$ms=0  
  tms = tabulate(chosenMales, nbins = N)#table of mating success
  
  
  #If there is a filename, files must be written
  #Columns are separated by ";" and "," is used as decimal separator
  if(!(is.na(filename)))
  {
    write.table(x = females, file = paste(filename, "Females.csv", sep=""),
                append = append, sep = ";", dec=",", row.names = FALSE, col.names = !append)
    write.table(x = males, file = paste(filename, "Males.csv", sep=""),
                append = append, sep = ";", dec=",", row.names = FALSE, col.names = !append)
    write.table(x = groups, file = paste(filename, "Groups.csv", sep=""),
                append = append, sep = ";", dec=",", row.names = FALSE, col.names = !append)
  }
  
  invisible(list(females=females, males=males,groups=groups))
  
}
#====================================================================