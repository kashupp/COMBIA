#' @import gdata
#' @import hash
#' @import oro.nifti
#' @import latticeExtra
#' @import lattice 
#' @importFrom 'grDevices'  'colorRampPalette'  'dev.off' 'png'
#' @importFrom 'stats' 'optimize' 'quantile' 'sd' 'update'
#' @importFrom 'utils' 'combn' 'read.csv' 'read.table' 'write.csv'
{}


####Load libraries
library(gdata)
library(hash)
library(oro.nifti)
library(latticeExtra)
library(lattice)


# Function to calculate predicted survival values
calculateSiPredicted <- function( conc, para, h){
    predictedSIs <- (1 / (1 + ((conc/ exp(para)  )^h)))
    return(predictedSIs) 
  }

# Function to calculate sum of squared error
J <- function(sis, conc, para, h){
    siPredicted <- apply( matrix(conc, nrow=1 , ncol=length(conc)), 2, calculateSiPredicted, para, h )
    j <- (1/(2*length(conc))) * sum((sis-(siPredicted))^2) 
    return(j)
  }

# Function to calculate gradiant of change of IC50 values
gradJ_a <- function(sis, conc, para, h){
    siPredicted <- apply( matrix(conc, nrow=1, ncol=length(conc)), 2, calculateSiPredicted, para, h )
    errors <- sis - siPredicted
    gradj_a <- -1 * (1/length(conc)) *   sum((errors * ((siPredicted)^2) ) * h * ((conc/exp(para))^h)  ) 
    return(gradj_a) 
  }

# Function to calculate gradiant of change of h values
gradJ_h <- function(sis, conc, para, h){
    siPredicted <- apply( matrix(conc, nrow=1, ncol=length(conc)), 2, calculateSiPredicted, para, h )
    errors <- sis - siPredicted
    gradj_h <- -1 * ((1/length(conc)) *  sum( (errors- ( ((siPredicted)^2) ) *h *((conc/exp(para))^(h-1) ))  ) )  
    #gradj_h <- (1/length(conc)) *  sum( (errors* ((siPredicted)^2) ) * ((conc/a)^h) * log(conc/a) )  
    return(gradj_h)
  }

# Function to extract maximum error
Jmax <- function(sis, conc, para, h){
    siPredicted <- apply( matrix(conc, nrow=1, ncol=length(conc)), 2, calculateSiPredicted, para, h )
    return( max(abs(sis-siPredicted)) )
  }


nlsComIntAct <- function(sis, conc, a, h) {
  if (a<=0 ||h<=0 ){
    print(sis)
    print(conc)
    print(a)
    print(h)
    stop("Reasonable prediction of IC50 or Hill Coefficient can not be made on this data ")
  }
  
  
  hstartX <- 0    # Variable to store h values
  c50startX <-  0 # Variable to store IC50 values
  EX <- 0         # Variable to store mean of Error sequare
  EXmax <- 0      # Variable to store maximum of mean of Error sequare
  gXa <- 0        # Variable to store gradiant  of IC50
  gXh <- 0        # Variable to store gradiant of h
  ErrorValues <- 0 # For diagnostics
  ErrorValuesMax <- 0 # For diagnostics
  alpha_aXStore <- 0 # For diagnostics
  alpha_hStore <- 0 # For diagnostics
  alpha_c50 <-0      # Variable to store the alpha values of the parameter IC50
  alpha_c50Store <- 0 # For diagnostics
  
  # Starting guess of h and ICs
  i <- 1
  hstartX[i] <- h
  c50startX[i] <- a
  alpha_c50[i] <- log(a)
  
  # gradient of alpha_c50
  gXa[i] <- gradJ_a(sis, conc, alpha_c50[i], hstartX[i])
  # gradient of h
  gXh[i] <- gradJ_h(sis, conc, alpha_c50[i], hstartX[i])
  # root mean sequare error
  EX[i] <- J(sis, conc, alpha_c50[i], hstartX[i])
  # maximum error
  EXmax[i] <- Jmax(sis, conc, alpha_c50[i], hstartX[i])
  ErrorValues[i] <- EX[i] #Diagnostics
  ErrorValuesMax[i] <- EXmax[i] #Diagnostics
  
  # StepSIze
  alpha_aX <- alpha_c50[i]/ (1000 * abs(gXa[i]) )
  alpha_hX <-   hstartX[i]/(1000 * abs(gXh[i]))
  alpha_aXStore[i] <- alpha_aX # For diagnostics
  alpha_hStore[i] <- alpha_hX # For diagnostics
  
  flagContinue <- TRUE
  cStore <- 0
  hStore <- 0
  c50startX[i] <- exp(alpha_c50[i])
  cStore[i] <- c50startX[i]
  hStore[i] <- hstartX[i]
  alpha_c50Store[i] <- alpha_c50[i]
    
  maxItr <- 1000
  
  while(flagContinue==TRUE){
    i <- i+1
    alpha_c50[i] <- alpha_c50[i-1] - alpha_aX * gXa[i-1]
    hstartX[i] <- hstartX[i-1] - alpha_hX * gXh[i-1]
    
    EX[i] <- J(sis, conc, alpha_c50[i], hstartX[i])
    EXmax[i] <- Jmax(sis, conc, alpha_c50[i], hstartX[i])
    ErrorValues[i] <- EX[i] # Diagnostics
    ErrorValuesMax[i] <- EXmax[i] # Diagnostics
    
    # Stop Criteria
    if (EXmax[i] < 0.05){
      flagContinue <- FALSE
    }
  
    if ( (i-1) >= maxItr ){
      flagContinue <- FALSE
    }
  
  # Adaptive selection of step size
  cStore[i] <- exp(alpha_c50[i])
  hStore[i] <- hstartX[i]
  alpha_c50Store[i] <- alpha_c50[i] 
    
    if (EX[i] < EX[i-1]){
      
      alpha_aX <- 1.05 * alpha_aX
      alpha_hX <- 1.05 * alpha_hX 
      alpha_aXStore[i] <- alpha_aX  # For diagnostics
      alpha_hStore[i] <- alpha_hX  # For diagnostics
      
      gXa[i] <- gradJ_a(sis, conc, alpha_c50[i], hstartX[i])
      gXh[i] <- gradJ_h(sis, conc, alpha_c50[i], hstartX[i])
      
    } else {
    
      c50startX[i] <- c50startX[i-1] 
      hstartX[i] <- hstartX[i-1] 
      alpha_c50[i] <- alpha_c50[i-1]
      EX[i] <- EX[i-1]
      alpha_aX <- 0.95 * alpha_aX
      alpha_hX <- 0.95 * alpha_hX 
      alpha_aXStore[i] <- alpha_aX  # For diagnostics
      alpha_hStore[i] <- alpha_hX  # For diagnostics
      gXa[i] <- gXa[i-1]
      gXh[i] <- gXh[i-1]
    }
  
  } # End of while
return(list(c( exp(alpha_c50[i]) , hstartX[i], EXmax[i])) )
}

# Function to select values of SI and H for estimating starting parameters of IC50 and h 
startingGuessIC50nH <- function(drugObs_Mean, ConcentrationCleaned )
  {
    if (any(drugObs_Mean < 0.5) )
    {
      drugObs_Mean2Index <-  which(max( drugObs_Mean[which(drugObs_Mean < 0.50)]) == drugObs_Mean)      
      drugObs_Mean2 <- drugObs_Mean[drugObs_Mean2Index ] 
      drugObs_MeanIC2 <- ConcentrationCleaned[drugObs_Mean2Index ]
      if (drugObs_Mean2 <= 0 ) {drugObs_Mean2 <- 0.01}
      } else{
        drugObs_Mean2 <- min(drugObs_Mean) 
        if (drugObs_Mean2 <= 0 ){drugObs_Mean2 <- 0.01}
        drugObs_MeanIC2 <- ConcentrationCleaned[ which(min(drugObs_Mean) == drugObs_Mean) ]
     }
    
  
  if ((drugObs_Mean2 < 0.5)  )
  {
    searchAbleValues <- drugObs_Mean[1: (drugObs_Mean2Index-1) ]
    ConcentrationCleanedSearchable <- ConcentrationCleaned[1: (drugObs_Mean2Index-1)]
    minInd80  <- order(abs(searchAbleValues- 0.8))[1:2]
    drugObs_Mean1 <-   mean(searchAbleValues[ minInd80[!is.na(minInd80)==TRUE ] ])
    if (drugObs_Mean1 >= 1 ){drugObs_Mean1 <- 0.99}
    drugObs_MeanIC1 <- mean(ConcentrationCleanedSearchable[  minInd80[!is.na(minInd80)==TRUE ]  ] )
    } else {
      drugObs_Mean1 <- mean(drugObs_Mean[3:4] )
      if (drugObs_Mean1 >= 1 ){drugObs_Mean1 <- 0.99}
      drugObs_MeanIC1 <- mean(ConcentrationCleaned[3:4])
    }
  
  
  r <- 1
  while (drugObs_MeanIC2 <= drugObs_MeanIC1)
  { #second best 
    secondbestind <- which(rank(drugObs_Mean)==r+1)
    drugObs_Mean2 <- drugObs_Mean[secondbestind] 
    if (drugObs_Mean2 <= 0 ){drugObs_Mean2 <- 0.01}
    drugObs_MeanIC2 <- ConcentrationCleaned[secondbestind]
    r <- r+1
    if (r > 3)
    {
      break;
    }
  }
  
  
  if ((drugObs_Mean1 < drugObs_Mean2 ) & (drugObs_Mean2 > 0.5))
  {
    drugObs_Mean1 <- mean(drugObs_Mean[1:2] )
    if (drugObs_Mean1 >= 1 ){drugObs_Mean1 <- 0.99}
    drugObs_MeanIC1 <- mean(ConcentrationCleaned[1:2])
  }
  
    return(c(drugObs_Mean2,drugObs_MeanIC2, drugObs_Mean1,drugObs_MeanIC1  ))
  }


#' This function applies Loewe Model
#' @param xConcentration X drug concentrations
#' @param yConcentration Y drug concentrations
#' @param drugYObs_Mean Mean of y drug observations
#' @param drugXObs_Mean Mean of x drug observations
#' @return Loewe Model values 
#' @examples 
#' xConcentration <- c(0.00,0.20, 0.39,  0.78,  1.56, 3.12, 6.25, 12.50, 25.00, 50.0) 
#' yConcentration <- c(128,  64,  32,  16,   8,   4,   2,   0)
#' drugXObs_Mean <- c(0.9747255, 0.9197924, 0.9520692, 0.9517162, 0.9032701, 0.7892114,
#'                      0.6768190, 0.6524227, 0.4561164)
#' drugYObs_Mean <- rev( c( 0.93, 0.89, 0.73, 0.42, 0.24, 0.21, 0.11) )
#' rslt <- loeweModel( xConcentration, yConcentration, drugYObs_Mean, drugXObs_Mean)
#' @author Muhammad kashif 
#' @export
loeweModel<- function ( xConcentration, yConcentration,
                        drugYObs_Mean , drugXObs_Mean)
  {
    noOfRows <- length(yConcentration)
    noOfCols <- length(xConcentration)
  
    # Variables to store calculated parameters
    cx50 <- 0 
    hx   <- 0
    cy50 <- 0
    hy   <- 0
  
    # Indices of the concentrations to be search for nls implementation
    xConcentrationCleaned <- xConcentration[2: length(xConcentration)] 
    yConcentrationRevCleaned <- rev(yConcentration[1: (length(yConcentration) - 1) ] )
    drugYObs_Mean_Rev <- rev(drugYObs_Mean)

    drugXObs_StartingGuess <- startingGuessIC50nH(drugXObs_Mean,xConcentrationCleaned )

   drugXObs_Mean2 <- drugXObs_StartingGuess[1] # si that is just below IC50 or minimum
   drugXObs_MeanIC2 <- drugXObs_StartingGuess[2] #  Concentration for above mentioned SI
   drugXObs_Mean1 <- drugXObs_StartingGuess[3]  # si that is close to 80
   drugXObs_MeanIC1 <- drugXObs_StartingGuess[4] # Cocnetration for the above mentioned SI
   
   drugYObs_StartingGuess <- startingGuessIC50nH(drugYObs_Mean_Rev, yConcentrationRevCleaned )
   
   drugYObs_Mean_Rev2 <- drugYObs_StartingGuess[1]
   drugYObs_Mean_RevIC2 <- drugYObs_StartingGuess[2]
   drugYObs_Mean_Rev1 <- drugYObs_StartingGuess[3]
   drugYObs_Mean_RevIC1 <- drugYObs_StartingGuess[4]
  
  
  # Starting guess of h and ICs
  hstX <- log10( (drugXObs_Mean1/drugXObs_Mean2) * (( 1- drugXObs_Mean2)/(1-drugXObs_Mean1 ) )  ) / 
          log10(drugXObs_MeanIC2/drugXObs_MeanIC1)
  # Incase the slope of dose response is almost constant and proper initialization is not possible then
  # it is changed to 0.1
  if( (hstX <  0) | is.nan(hstX) ){
    #print("hstx<0 was true")
    hstX<-0.1}
  
  hstY <- suppressWarnings( log10((drugYObs_Mean_Rev1/drugYObs_Mean_Rev2) * (( 1- drugYObs_Mean_Rev2)/(1-drugYObs_Mean_Rev1 ) )  ) / 
          log10(drugYObs_Mean_RevIC2/drugYObs_Mean_RevIC1) )
  # Incase the slope of dose response is almost constant and proper initialization is not possible then
  # it is changed to 0.1
  if( (hstY <  0) | is.nan(hstY) ){
    #print("hsty<0 was true")
    hstY <- 0.1}
  
  denoX <-  (( (1-drugXObs_Mean1)/drugXObs_Mean1 ) ^ (1/hstX) )
  if (denoX==0){denoX <- 0.0001} 
  c50stX <- drugXObs_MeanIC1/denoX 
  
  denoY <- (( (1-drugYObs_Mean_Rev1)/drugYObs_Mean_Rev1 ) ^ (1/hstY) )
  if (denoY==0){denoY <- 0.0001} 
  c50stY <- drugYObs_Mean_RevIC1/ denoY

  paramValuesX <- 0
  paramValuesY <- 0
  
  paramValuesX <- nlsComIntAct(drugXObs_Mean, xConcentrationCleaned, c50stX, hstX)
  paramValuesY <- nlsComIntAct(drugYObs_Mean_Rev, yConcentrationRevCleaned, c50stY, hstY) 
  cx50 <- paramValuesX[[1]][1]
  hx   <- paramValuesX[[1]][2]
  cy50 <- paramValuesY[[1]][1]
  hy   <- paramValuesY[[1]][2]
  

  concentration.product <- expand.grid(xConcentration[2:length(xConcentration)], yConcentration[1: (length(yConcentration)-1)] )  # X and Y cons with PBS well  removed 
  
  # Inverse survival function
  Inv_Survival <- function(s, h, c50,alpha) {
      ( ( (((1- ( s ^  (1/ alpha))    )/  (s ^ (1/ alpha))     )) ^ (1/h  ) )* c50 ) 
  }
  
  # Function SurvivalAB implements the equation of loewe response surface 
  SurvivalAB <- function(s, da, db, cah, ca50, cbh, cb50,alpha){ 
    abs(((da/ Inv_Survival(s, cah, ca50,alpha)) + 
           (db/ Inv_Survival(s, cbh, cb50, alpha)))- 1)
  }
  
  survivalindex.normal <- 0
  for (i in 1: nrow(concentration.product)){#i=1  
    survivalindex.normal.intermediate <- optimize(SurvivalAB, c(0.0001,1), concentration.product[i,1], concentration.product[i,2],
                                                  hx, cx50, hy, cy50, alpha=1, tol= 1e-16)
    
    survivalindex.normal[i] <- as.numeric(survivalindex.normal.intermediate[1])
  }
  
  # concentration.product and  combXYObs_Mean have same dimension
  # and by default matrix will create it bycol and therefore it will be 
  # in wrong order  therefore
  loeweSynObs_Model <- matrix(survivalindex.normal, ncol=noOfCols-1, nrow=noOfRows-1, byrow=TRUE)
  
  return( list(loeweSynObs_Model) )
  
} # End of loeweModel function



#' This function calculates Loewe synergy/antagonism and associated BIs
#' @param rawDataPreProcessed Raw preprocessed experimental data 
#' @param xConcentration X drug concentrations
#' @param yConcentration Y drug concentrations
#' @param nBoot Number of times to bootstrap
#' @return Three lists, first list consisting of Loewe Synergy/Antagonism, lower bound of BI and  
#' upper bound of BI. 2nd list consists of global BI for maximum synergy and 3rd list 
#' consists of global BI of maximum antagonistic combination.
#' @examples
#'\dontrun{
#' dataFile <- system.file("extdata", "rawDataPreProcessed.csv", package="COMBIA")
#' dataSample <- read.csv(dataFile, header=FALSE )
#' xConc  <- c(0.00,  0.20,  0.39,  0.78,  1.56,3.12,  6.25, 12.50, 25.00, 50) 
#' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
#' noOFBoot <- 500 # a large number is recomended
#' rslt <- applyLoewe(as.matrix(dataSample), xConc, yConc, noOFBoot)
#' }
#' @author Muhammad kashif
#' @export
applyLoewe <- function( rawDataPreProcessed, xConcentration, yConcentration, nBoot)
  {
  
    # calculates total number of rows and colmuns
    noOfRows <- length(yConcentration); # calculate directly from  data;
    noOfCols <- length(xConcentration); # calculate directly from  data;
 
    # calculates number of replicates 
    rawDataPreProcessed_NA <- rawDataPreProcessed;
    rawDataPreProcessed_NA[which(rawDataPreProcessed == 0, arr.ind=TRUE)] <- NA
    replicateCount_individual <-   apply(rawDataPreProcessed_NA, 2, function(x) length (which(!is.na(x)) ))
    
    # Prepare data for Loewe analysis
    
    # Mean of data
    totalNumberofReplicates <- nrow(rawDataPreProcessed)
    rawDataPreProcessedMean <- apply(rawDataPreProcessed_NA,2,mean, na.rm=TRUE) 
    
    rawDataPreProcessed_mat_temp <- matrix(rawDataPreProcessedMean, noOfRows, noOfCols)
    
    drugYObsMean <- as.vector(rawDataPreProcessed_mat_temp[1:noOfRows-1 ,1])  # all rows of first column, descending order
    drugXObsMean <- as.vector(rawDataPreProcessed_mat_temp[noOfRows,2:noOfCols]) # all columns of last row ascending oredr
    combXYObsMean <- rawDataPreProcessed_mat_temp[1:noOfRows-1, 2:noOfCols ]# x ascending and y descending
  
    # Start of applying Loewe  
    loeweSynObsModelAll <- loeweModel( xConcentration, yConcentration,
                                    drugYObsMean,drugXObsMean )
    loeweSynObsModel <- as.matrix( loeweSynObsModelAll[[1]] )
     # Note: Data in loeweSynObs_Model is formated simialr to combXYObs_Mean 
    
    # Calculate Loewe Synergy
    # data is in matrix loeweSynergy with same order as loeweSynObs_Model 
    
    loeweSynergy <-  loeweSynObsModel - combXYObsMean # positive is synergy
    loeweSynergy_vector <- as.vector(loeweSynergy)
    # End of loewe Application
    # calculate BIs 

    # Non Heteroscedastic residual selection based on sliding window of 10% nearest neighbours.
    residuesWindowsList <- residualSelection(totRepl=totalNumberofReplicates, nr=noOfRows , 
                                             nc=noOfCols, rawDataPPNA=rawDataPreProcessed_NA )
    
    # CREATE NEW RESIDUES THAT IS RESIDUES DUE TO VARIABILITY IN OBSERVED VALUES
    reps <- 1000
    residualVar <- matrix(nrow = ( (noOfRows-1) * (noOfCols-1) ) , ncol =  reps+1 ) ; # variable to hold the residuals due to variability
    
    residualVar[,1]     <- as.vector(combXYObsMean)
    
  # TO EXRACT INDEX OF RESIDUES LISTS OF THE X DRUG, Y DRUG AND COMBINATIONS
    matForInd <- matrix(1:(noOfRows* noOfCols), noOfRows, noOfCols)
    drugXInds <- matForInd[noOfRows, 2:noOfCols] # all columns of last row ascending order of concentration
    drugYInds <- matForInd[1:noOfRows-1, 1] # all rows of first column in the descending order of concentration
    drugXYInds <- matForInd[1:noOfRows-1, 2:noOfCols ]# x ascending and y descending
    
    for( i in (1:reps)) # 
    { 
      drugYObsMeanBar <- 0
      drugXObsMeanBar <- 0
      
      for(yl in 1:length(drugYObsMean))# yl<- 1
      { 
        drugYObsMeanBar[yl] <- drugYObsMean[yl] + sample( unlist( residuesWindowsList[  drugYInds[yl] ]) , 1 , replace=TRUE)
      }
      
      for(xl in 1:length( drugXObsMean))
      { 
        drugXObsMeanBar[xl] <- drugXObsMean[xl] + sample(unlist(residuesWindowsList[ drugXInds[xl] ]) , 1 , replace=TRUE)
      }
        if ((i %% 10)==0)
        {
        print(paste(" Bootstrap no=", i))
        }
        #options(warn=2)
        loeweSynObs_Model_residualVar_All <- loeweModel(xConcentration, yConcentration,
                                                        drugYObsMeanBar, drugXObsMeanBar)
        
        
        loeweSynObs_Model_residualVar<- as.matrix(loeweSynObs_Model_residualVar_All[[1]] )
        
        loeweSynergy_Boot <-  loeweSynObs_Model_residualVar - loeweSynObsModel # residuals
    
        residualVar[,i+1] <- as.vector( loeweSynergy_Boot )
        
    }
    
    
    # calculate the eroor distrbution of the every combination by selecting errors from two different 
    
    # index of combination residues in the list
    drugXYIndsV <- as.vector(drugXYInds)
    
    error.loewe.bootstrap <- matrix(rep(NA, ((noOfCols-1) * (noOfRows-1))  * nBoot ),
                                      nrow=(noOfCols-1) * (noOfRows-1) , ncol=nBoot) # varaible to store the loewe error distribution 
  
    quantilesCI_ConcComb <- matrix(rep(0, 3*((noOfCols-1) * (noOfRows-1)) ), 3, ((noOfCols-1) * (noOfRows-1)) ) # variable to store the bis
    
    for(i in  (1: ((noOfCols-1) * (noOfRows-1)) ) ) #xDrug Mean  i=1
    {    
          # sample El and Eb
          # Note: eb and el should correspond to a combinations and thats why indices are used in residualVar
          error.loewe.bootstrap[ i, (1:nBoot) ] <- sample( unlist(residuesWindowsList[  drugXYIndsV[i] ] ) , nBoot , replace=TRUE)
                                                  + sample(residualVar[ i , 2:(reps+1) ] , nBoot , replace=TRUE)                    
    }
    
    # 0.025= Lower Bound and 0.975 is Upper bound
    # the first row is the loewe synergy, second in the lower bound and 3rd is the upper boound of synergy.
    quantilesCI_ConcComb <- matrix(rep(0, 3*((noOfCols-1) * (noOfRows-1)) ), 3, ((noOfCols-1) * (noOfRows-1)) )
    
    for (pv in (1:((noOfCols-1) * (noOfRows-1))) )
    {
      # 0.025= Lower Bound and 0.975 is Upper bound
     quantilesCI_ConcComb[2:3,pv] <-
     apply(matrix(error.loewe.bootstrap[pv,], 1, nBoot) , 1, quantile, c(.025, 0.975  ))
    }  
    
      quantilesCI_ConcComb[1,] <- loeweSynergy_vector
      quantilesCI_Max <-    apply(matrix(apply(error.loewe.bootstrap,2, max), 1, nBoot) , 1, quantile, c(.025, 0.975  ))
      quantilesCI_Min <-    apply(matrix(apply(error.loewe.bootstrap,2, min), 1, nBoot) , 1, quantile, c(.025, 0.975  ))
      
    return(list(quantilesCI_ConcComb, quantilesCI_Max, quantilesCI_Min))
    
      
 }# end of applyLoewe function



#' This function extract numerical indices of a given range e-g B2  
#' @param range Range e-g B2
#' @param excelFormate TRUE if range is in spreadsheet formate 
#' @return Number of starting row, ending row, starting column and ending column
#' @examples
#' rng <- c("B2")
#' exclF <- TRUE
#' rslt <-  extractValuesFromRange(rng, exclF)
#' @author Muhammad kashif
#' @export
extractValuesFromRange <- function(range, excelFormate)
{
  # Given a list a range("B2") This function will extract 
  # the numeric starting row, ending row, starting column and ending column.
  if (excelFormate==FALSE){
      columnstartindex <- as.integer(substr(range[1], 2, 4)) #extract strating column
      rowstartindex <- which( letters[1 : 26] == tolower(substr(range[1], 0, 1))) #extract strating row
      columnendindex <- as.integer(substr(range[2], 2, 4))  #extract ending column
      rowendindex <- which( letters[1 : 26] ==tolower( substr(range[2], 0, 1))) #extract ending row
  
      # Update variables for reformating data
      noOfRows <- (rowendindex - rowstartindex) + 1 # Calcualte No of rows in a single data
      noOfCols <- (columnendindex - columnstartindex) + 1 # Calculate No of cols in a single data
      return(c(rowstartindex, rowendindex, columnstartindex, columnendindex ))  
  
  } else { #  if input of range is similar to excel formate
      rowstartindex <- as.integer(substr(range[1], 2, 4)) #extract strating column
      columnstartindex <- which( letters[1 : 26] == tolower( substr(range[1], 0, 1)) )#extract strating row
      rowendindex <- as.integer(substr(range[2], 2, 4))  #extract ending column
      columnendindex <- which( letters[1 : 26] == tolower( substr(range[2], 0, 1)) )#extract ending row
      
      # Update variables for reformating data
      noOfRows <- (rowendindex- rowstartindex) + 1 # Calcualte No of rows in a single data
      noOfCols <- (columnendindex- columnstartindex) + 1 # Calculate No of cols in a single data
      return(c(rowstartindex, rowendindex, columnstartindex, columnendindex ))    
  }
  
}


#' This function plots the synergy analysis 2D and 3D graphs
#' @param processedData A matrix to plot
#' @param xConcentration X drug concentrations
#' @param yConcentration Y drug concentrations
#' @param xDrug X drug name
#' @param yDrug Y drug name
#' @param cellLine Cell line name
#' @return Plot the values 
#' @examples
#' dataFile <- system.file("extdata", "processedData.csv", package="COMBIA")
#' procData <- read.csv( dataFile, header=FALSE)
#' xConc <- c(0.00,  0.20, 0.39, 0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50) 
#' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
#' xD <- "X_Drug"
#' yD <- "Y_Drug"
#' clN <- "myCell"
#' rslt <- synAntPlot(as.matrix(procData),xConc,yConc, xD, yD, clN)  
#' @author Muhammad kashif
#' @export
synAntPlot<- function(processedData, xConcentration, yConcentration, xDrug, yDrug
                      , cellLine) 
  {
    cellLine <- cellLine;
    noOfRows <- length(yConcentration); # calcualte directly from  data;
    noOfCols <- length(xConcentration); # calcualte directly from  data;
  
    # Color Key
    ColorFun <- colorRampPalette(tim.colors(255))
    #Levelplot plots the columns and first column at botto so change arrangement of data as:
    #Plot only those are significantat = do.breaks( c(-60 , 60), 255)
    xConcentration_label <- xConcentration
    yConcentration_label <- yConcentration
  
    objectSynAntPlot <- levelplot( 100* (t(processedData[nrow(processedData):1, ])) , col.regions= ColorFun(255), 
                                    at = do.breaks( c(-115 , 115), 255),
                                    scales=list( x= list( at=c(1: (noOfCols-1) ), labels=xConcentration_label[2:length(xConcentration_label)] ,
                                                cex=1.3 ),  y= list( at=c(1:(noOfRows-1) ), 
                                                labels=rev(yConcentration_label[1:(length(yConcentration_label)-1)]), cex=1.3)
                                                ),
                                    xlab= as.vector(xDrug) , ylab= as.vector(yDrug), main=""
                                    ,colorkey=list(labels= list(cex=1.3,at=c(-100,-50,0,50,100), 
                                                          labels= as.character(c(-100,-50,0,50,100)))
                                             ) )
   
   
    # Update the font size                           
    objectSynAntPlotUpdated <- update(objectSynAntPlot,par.settings = list(fontsize = list(text = 14, points = 14),
                                                                            par.ylab.text = list(cex = 1.5) ,
                                                                            par.xlab.text = list(cex = 1.5) 
                                                                           ))
    # Print Formated plot
    print(objectSynAntPlotUpdated)
    
    #Data saving code is removed due to change in the Cran policy 2018-10-07  
    myPathRoot<- system.file( package="COMBIA")
    #Add prefiX Combined_AnalyzedData
    myPath <- paste(myPathRoot,   "AnalyzedData_", sep="")
    myPath <- paste(myPath,   xDrug, sep="")
    myPath <- paste(myPath,   "_", sep="")
    myPath <- paste(myPath,   yDrug, sep="")
    myPath <- paste(myPath,   "_", sep="")
    myPath <- paste(myPath,   cellLine, sep="")
    myPath2 <- paste(myPath,   "_3D.png", sep="")
    myPath <- paste(myPath,   ".png", sep="")
    #png(filename=myPath, width=800, height=600) #Data saving code is removed 2018-10-07
    #print(objectSynAntPlotUpdated) #Data saving code is removed 2018-10-07
    #dev.off() #Data saving code is removed 2018-10-07

  # For 3D graphs
      processedDataMat<-as.matrix(processedData)
      colnames(processedDataMat) <- c(1:ncol(processedDataMat))
      rownames(processedDataMat) <- c(1:nrow(processedDataMat))
      h <- cloud( processedDataMat , 
              col.facet= level.colors(  as.matrix(processedDataMat) , at = do.breaks( c(-1.15 , 1.15), 255), colors =TRUE,
                        col.regions=ColorFun(255)), 
              zlim=c(-1.15, 1.15),
              colorkey=list(col=ColorFun(255), at = do.breaks( c(-1.15 , 1.15), 255),
                         labels= list(cex=1.3, at=c(-1,-0.5,0,0.5,1), 
                                      labels= as.character(c(-100,-50,0,50,100)))
              ),
              xlab= paste(as.vector(yDrug),"  ") , ylab= paste("        ", as.vector(xDrug)), main="", zlab="",
              scales=list( y= list(arrows=FALSE, at=c(1: (noOfCols-1) ), lab= xConcentration_label[2:length(xConcentration_label)], cex=1.1),  
                           x= list(arrows=FALSE, at=c(1:(noOfRows-1) ), lab=  
                                     c(yConcentration_label[1:(length(yConcentration_label)-2)], paste(yConcentration_label[(length(yConcentration_label)-1)], "   ", sep="") ), 
                                   cex=1.1),
                           z= list(arrows=FALSE, at=c(-1,-0.5,0,0.5,1),lab=c(-100,-50,0,50,100), cex=1.2)
              ),
              panel.3d.cloud=panel.3dbars, type="h",reference=TRUE,
              screen = list(z = -40, x = -25)
            )
    
  hUpdated <- update(h, par.settings = list(fontsize = list(text = 14, points = 14),
                                         par.ylab.text = list(cex = 1.5) ,
                                         par.xlab.text = list(cex = 1.5) 
                    ))
    
  print(hUpdated)
  
  #png(myPath2, width=850, height=600) #Data saving code is removed 2018-10-07
  #print(hUpdated) #Data saving code is removed 2018-10-07
  #dev.off() #Data saving code is removed 2018-10-07
  
  print(as.matrix(processedData))
  
  
  # End of data saving
  }

#' Function calculates significant synergy/antagonism 
#' @param synergyCalculationLists List of synergy antagonism calculations
#' @param noOfRows Number of rows
#' @param noOfCols Number of columns
#' @param xDrug Name of drug on x-axis
#' @param yDrug Name of drug on y-axis
#' @param cellLine Cell Line
#' @return Processed data 
#' @examples
#' dataFile <- system.file("extdata", "rawDataPreProcessed.csv", package="COMBIA")
#' dataSample <- read.csv(dataFile, header=FALSE)
#' nR <- 8
#' nC <- 10
#' rslt <- applyBliss(nR, nC,  as.matrix(dataSample ), 100) 
#' synergySignificant(rslt, nR, nC,"A", "B", "Cell")
#' @author Muhammad kashif
#' @export
synergySignificant <- function(synergyCalculationLists, noOfRows, noOfCols, xDrug, yDrug, cellLine )
{
  #Separate the lists synergyCalculationLists<-synergyBlissCalculationLists
  
  synergy_Calculation_Matrix <- matrix(unlist(synergyCalculationLists[[1]]), nrow=3, ncol=(noOfRows-1)* (noOfCols-1))
  CIForBestSynergy <- unlist(synergyCalculationLists[[2]])
  CIForBestSynergy[3] <- max(synergy_Calculation_Matrix[1,])
  CIForBestAntagonism <- unlist(synergyCalculationLists[[3]]) 
  CIForBestAntagonism[3] <- min(synergy_Calculation_Matrix[1,])
  
  # Formating analyzed data same way as of the input data
  # Note: While calclating Synergy  first row counterpart in synergy_Calculation_Matrix
  # should be postive and should be nagative for antagonism.
  
  #plot graph of those values that are synergistic or antagonistic and are 
  # at 95% significant and bonferroni corrected
  levelOut <-  rep(0, ( (noOfRows-1)* (noOfCols-1) ) )   # variable to hold the data to print out
  indecesOfSynergy  <- which (synergy_Calculation_Matrix[1,] > 0)  #IndecesOfSynergy
  indecesOfSynergy_Significant95  <- which( synergy_Calculation_Matrix[1,indecesOfSynergy] > synergy_Calculation_Matrix[3,indecesOfSynergy] ) #only those indeces which are significant 95% and bonferroni corrected
  
  indecesOfAntagonism <- which (synergy_Calculation_Matrix[1,] < 0)
  indecesOfAntagonism_Significant95  <- which( synergy_Calculation_Matrix[1,indecesOfAntagonism] < synergy_Calculation_Matrix[2,indecesOfAntagonism] ) #only those indeces which are significant and bonferroni corrected 95%
  #Update outputvariable
  levelOut[indecesOfSynergy[indecesOfSynergy_Significant95] ] <- 
    synergy_Calculation_Matrix[1,indecesOfSynergy[indecesOfSynergy_Significant95] ]
  
  levelOut[indecesOfAntagonism[indecesOfAntagonism_Significant95] ] <- 
    synergy_Calculation_Matrix[1,indecesOfAntagonism[indecesOfAntagonism_Significant95] ]
  
  processedData <- matrix(levelOut, nrow=noOfRows-1, ncol=noOfCols-1 )
  
  
  
  myPathRoot <- system.file(  package="COMBIA")
  #Data saving code is removed due to change in crean policy 2018-10-07
  myPath <- paste(myPathRoot,   "AnalyzedData_", sep="")
  myPath <- paste(myPath,   xDrug, sep="")
  myPath <- paste(myPath,   "_", sep="")
  myPath <- paste(myPath,   yDrug, sep="")
  myPath <- paste(myPath,   "_", sep="")
  myPath <- paste(myPath,   cellLine, sep="")
  myPath <- paste(myPath,   ".csv", sep="")
  #write.csv(100*processedData, file=myPath ) #Data saving code is removed 2018-10-07
  
  
  #Save Global CIS for Synergy and Antagonism
  myPath <- 0
  myPath <- paste(myPathRoot,   "AnalyzedData_", sep="")
  myPath <- paste(myPath,   xDrug, sep="")
  myPath <- paste(myPath,   "_", sep="")
  myPath <- paste(myPath,   yDrug, sep="")
  myPath <- paste(myPath,   "_", sep="")
  myPath <- paste(myPath,   cellLine, sep="")
  myPath <- paste(myPath,   "_", sep="")
  myPath <- paste(myPath,"CIs", sep="")
  myPath <- paste(myPath,   ".csv", sep="")
  #Data saving code is removed 2018-10-07
  # write.csv(100*c(CIForBestSynergy, CIForBestAntagonism), file=myPath ) #Data saving code is removed 2018-10-07
  
  
  #Save local well statistical data for Synergy and Antagonism
  myPath <- 0
  myPath <- paste(myPathRoot,   "AnalyzedData_", sep="")
  myPath <- paste(myPath,   xDrug, sep="")
  myPath <- paste(myPath,   "_", sep="")
  myPath <- paste(myPath,   yDrug, sep="")
  myPath <- paste(myPath,   "_", sep="")
  myPath <- paste(myPath,   cellLine, sep="")
  myPath <- paste(myPath,   "_", sep="")
  myPath <- paste(myPath,"LocalData", sep="")
  myPath <- paste(myPath,   ".csv", sep="")
  #write.csv(100*synergy_Calculation_Matrix, file=myPath ) #Data saving code is removed 2018-10-07
  
  return(processedData)
}# End of synergySignificant


#' This function will takes a list of ranges removes case wells and extract replicate values separately
#' @param rawDataUnProcessed A data matrix
#' @param wellRanges Ranges of wells
#' @param wellplace Place of treated (case) well range
#' @param simple TRUE if survival values are already calculated otherwise it is FALSE
#' @param excelFormate True if ranges are in excel formate
#' @return Replicate values
#' @examples 
#' dataFile <- system.file("extdata", "testData.csv", package="COMBIA")
#' rData <- read.csv( dataFile, skip=0, sep=",", nrows=41, 
#'                     fill=TRUE, header=FALSE,
#'                     blank.lines.skip = FALSE)[,1:13]
#' wellR= c( "l3:l10","m3:m10","b3:k10",  "l13:l20","m13:m20","b13:k20", 
#'             "l23:l30","m23:m30","b23:k30",  "l33:l40","m33:m40","b33:k40")
#' rslt <-  extractReplicateValues(rData, wellR, excelFormate=TRUE )
#' @author Muhammad kashif
#' @export
extractReplicateValues <- function(rawDataUnProcessed, wellRanges, wellplace=3, simple=FALSE, excelFormate=FALSE)
  { 
    # This function will takes a list of ranges and then 
    # removes case wells and extract replicate values separately
    if (simple==FALSE){
        cnt <- 0
        replicateCounter <- 0
        # List to hold the replicated data
        replicatedData <- vector("list", length(wellRanges)/3)
          while (cnt < (length(wellRanges)) ){
            #cnt=0, wellplace=1
            range <- unlist(strsplit( wellRanges[ cnt + wellplace ] ,":")) #split the range i-e from A2:D6 in to A2 and D6
            numericValues <- extractValuesFromRange(range, excelFormate)
            replicateCounter <- replicateCounter + 1;
            replicatedData[[replicateCounter]] <- rawDataUnProcessed[c(numericValues[1]:numericValues[2] ),
                                                      c(numericValues[3]:numericValues[4] )]
          cnt <- cnt + 3;
          }
  return(replicatedData)
  }
  else
  {
    cnt <- 0
    replicateCounter <- 0
    # List to hold the replicated data
    replicatedData<- vector("list", length(wellRanges)/3)
    while (cnt < (length(wellRanges)) ){
      #cnt=1
      range <- unlist(strsplit( wellRanges[ cnt + wellplace ] ,":")) #split the range i-e from A2:D6 in to A2 and D6
      numericValues <- extractValuesFromRange(range, excelFormate)
      replicateCounter <- replicateCounter + 1;
      replicatedData[[replicateCounter]] <- rawDataUnProcessed[c(numericValues[1]:numericValues[2] ),
                                                        c(numericValues[3]:numericValues[4] )]
      cnt <- cnt + 3;
    }
    
    return(replicatedData)
  }
  
  }



#' This function calculates CV
#' @param vals Values
#' @return cv of input values
#' @examples 
#' mData <- matrix(1:10, 2,5)
#' rslt <-  cVCal(mData)
#' @author Muhammad kashif
#' @export
cVCal <- function(vals){
        # Function to calculate CV[sd/mean ] of a vector
         cv <- sd(vals,na.rm=TRUE)/mean(vals,na.rm=TRUE);
    
         return(cv)
        }



#' Function to make unique perturbations of the replicates these will be used incase if CV is greater than threshold. 
#' @param totalNumberofReplicates Total replicate number
#' @return unique possible perturbations
#' @examples
#' rslt <- createUniquePertbs(5) 
#' @author Muhammad kashif
#' @export
createUniquePertbs <- function(totalNumberofReplicates){
    # Code to make unique perturbations of replicates that will be used incase if CV 
    # will be greater than 30%. 
    #totalNumberofReplicates<-8
  
    n <- totalNumberofReplicates; # Number of replicates will be decreasing
    perturblist <- vector("list",n) # 1 and 2
    for(i in (1:n) ){ #n=1  i=8
        perturblist[[i]] <- combn(n,i,unique )
    }
 #Always remove first and  last
  return(perturblist)
  }  



#' This function  Remove Outliers
#' @param arrangeReplicates A data matrix
#' @param minThersholdForCVCal Threshold for value removal in CV
#' @param minThersholdForCV  Values to be excluded
#' @return Replicate values
#' @examples
#' dataFile <- system.file("extdata", "rawDataPreProcessed.csv", package="COMBIA")
#' dataSample <- read.csv(dataFile, header=FALSE )
#' minThersholdForCV <- 0.3
#' minThersholdForCVCal <- 0.1
#' removeOutliers( as.matrix(dataSample ), minThersholdForCV,
#'        minThersholdForCVCal) 
#' @author Muhammad kashif
#' @export
removeOutliers <- function(arrangeReplicates, minThersholdForCVCal, minThersholdForCV)
  {
    # Add NA to data that was empty
    arrangeReplicates_NA <- arrangeReplicates;
    arrangeReplicates_NA[which(arrangeReplicates == 0)] <- NA
    # Count replicates for each data value
    replicateCount_individual <-   apply(arrangeReplicates_NA, 2, function(x) length (which(!is.na(x)) ))
    #Count which replicates should be used, Each list represents each row
    if (class(apply(arrangeReplicates_NA, 2, function(x) which(!is.na(x)) ))=="matrix")
        {
          replicateCount_indices <-   as.list(data.frame(apply(arrangeReplicates_NA, 2, function(x) which(!is.na(x)) )))
        } else {
          replicateCount_indices <-   apply(arrangeReplicates_NA, 2, function(x) which(!is.na(x)) )  
        }
    totalCombinationsnSingle <- (ncol(arrangeReplicates_NA))  
    for(i in  1:totalCombinationsnSingle  ){ 
        #calculate Coefficient of variance
        # If values are less than 10 SIs then dont do cv
        if (length(which( arrangeReplicates_NA[,i] > minThersholdForCVCal) ) > 0)
            {
              currCV <- cVCal(c(arrangeReplicates_NA[,i]) )
              # if CV is greater than 30 then apply function of removal 
              if ( (currCV > minThersholdForCV)==TRUE ){
                  cvfound <- 0
                  currentNoOfPerts <- 0
                  currentPertList <- 0
        
                  if (replicateCount_individual[i] > 2){
                      perturblist <- createUniquePertbs(replicateCount_individual[i])
                      rangeOfPerturbation <- 2:(length(perturblist)-1)
                      cntPertb <- length(rangeOfPerturbation)
                      indecesOfReplicatesIn_arrangeReplicates <- replicateCount_indices[[i]]
                      while(cntPertb > 0){
                            currentPertList <- perturblist[[rangeOfPerturbation[cntPertb]]]
                            pertSIval <- matrix(arrangeReplicates_NA[  indecesOfReplicatesIn_arrangeReplicates[currentPertList] ,i ], 
                                                      nrow=nrow(currentPertList), ncol=ncol(currentPertList) )
                            if(cntPertb > 1)
                              {
                                if(min(apply(pertSIval, 2, cVCal)) < minThersholdForCV)
                                  {
                                    indPert <- which(min(apply(pertSIval, 2, cVCal))==apply(pertSIval, 2, cVCal))
                                    indicesToBeZerod <- setdiff(c(replicateCount_indices[[i]]),
                                                            c(indecesOfReplicatesIn_arrangeReplicates[currentPertList[,indPert]]))
                                    arrangeReplicates[indicesToBeZerod,i] <- 0
                                    cvfound <- 1
                                break 
                              }
            
                              } else { #if(cntPertb==1)
                                    indPert <- which(min(apply(pertSIval, 2, cVCal))==apply(pertSIval, 2, cVCal))
                                    indicesToBeZerod <- setdiff(c(replicateCount_indices[[i]]),
                                                                c(indecesOfReplicatesIn_arrangeReplicates[currentPertList[,indPert]]))
                                    arrangeReplicates[indicesToBeZerod,i] <- 0
                              }
                      cntPertb <- cntPertb-1 
                      } # End of while
        
                    }# End of if replciates >2
                  } # End of if (length(which(currSIvals > 0.1) ) > 0)
                } # End of perturb loop
            } # End of for(i in 80)
   
    return(arrangeReplicates)
  }# End OF VARIABLITY CHECKING FUNCTION


#' Function calculates Bliss Synergy, associated BIs and global BIs
#' @param noOfRows Number of rows in the experiment
#' @param noOfCols Number of columns in the experiment
#' @param rawDataPreProcessed Data matrix 
#' @param nBoot Number of bootstrap
#' @return Three lists, first list consists of Bliss Synergy/Antagonism, lower bound of BI and 
#' upper bound of BI. 2nd list consists of global BI of Maximun synergistic combiantion and 3rd list 
#' consists of global BI of maximum antagonistic combination.
#' @examples 
#' dataFile <- system.file( "extdata", "rawDataPreProcessed.csv", package="COMBIA" )
#' dataSample <- read.csv(dataFile, header=FALSE )
#' nR <- 8
#' nC <- 10
#' rslt <- applyBliss(nR, nC,  as.matrix(dataSample ), 500) 
#' @author Muhammad kashif
#' @export
applyBliss <- function(noOfRows, noOfCols, rawDataPreProcessed , nBoot)
{
  rawDataPreProcessed_NA <- rawDataPreProcessed;
  rawDataPreProcessed_NA[which(rawDataPreProcessed == 0)] <- NA
  replicateCount_individual <-   apply(rawDataPreProcessed_NA, 2, function(x) length (which(!is.na(x)) ))
  
  totalNumberofReplicates <- nrow(rawDataPreProcessed)
  rawDataPreProcessedMean <- apply(rawDataPreProcessed_NA, 2, mean, na.rm=TRUE) 
  rawDataPreProcessed_mat_temp <- matrix(rawDataPreProcessedMean, noOfRows, noOfCols)
  drugYObsMean <- as.vector(rawDataPreProcessed_mat_temp[1:noOfRows-1 ,1])  # all rows of first column, descending order
  drugXObsMean <- as.vector(rawDataPreProcessed_mat_temp[noOfRows, 2:noOfCols]) # all columns of last row ascending oredr
  combXYObsMean <- rawDataPreProcessed_mat_temp[1:noOfRows-1, 2:noOfCols ]# x ascending and y descending
  combXYObsMeanVector <- as.vector(combXYObsMean)
  
  # Apply Bliss on experimental data
  # Application of Bliss model only for mean values
  comb_xy_model_Bliss_inter <- expand.grid( drugXObsMean, drugYObsMean )  # use all values , same formate as of  input data
  comb_xy_model_Bliss <-  (comb_xy_model_Bliss_inter[,1] * comb_xy_model_Bliss_inter[,2]) # generate model data for normal cells
  # comb_xy_model_Bliss_Mat is of formate as in the raw data fetch accordoring to wellRange parameter
  comb_xy_model_Bliss_Mat <- matrix(comb_xy_model_Bliss, noOfRows-1, noOfCols-1, byrow=TRUE) # formate as in inputdata
  
  # Calculation of Bliss Synergy
  synergy_Bliss_obs <-   comb_xy_model_Bliss_Mat - combXYObsMean  # if positive then synergy
  synergy_Bliss_obsVector <- as.vector(synergy_Bliss_obs)
  
# Sliding window based residuals
  
  # Non Heteroscedastic residual selection based on sliding window of 10% nearest neighbours.
  residuesWindowsList <- residualSelection(totRepl=totalNumberofReplicates, nr=noOfRows,
                                           nc=noOfCols, rawDataPPNA=rawDataPreProcessed_NA )

  replicateCount_individual_temp <- matrix(replicateCount_individual, noOfRows, noOfCols)
  drugYObsReplicates <- replicateCount_individual_temp[1:noOfRows-1, 1]  # all rows of first column, descending order
  drugXObsReplicates <- replicateCount_individual_temp[noOfRows, 2:noOfCols] # all columns of last row ascending oredr
  combXYObsReplicates <- replicateCount_individual_temp[1:noOfRows-1, 2:noOfCols ]# x ascending and y descending
  combXYObsReplicatesVector <- as.vector(combXYObsReplicates)

  #Index matrix
  matForInd <- matrix(1:(noOfRows*noOfCols), noOfRows, noOfCols)
  drugXInds <- matForInd[noOfRows, 2:noOfCols] # all columns of last row ascending oredr
  
  # Bootstrap
  synergy_Best_Bar <- matrix(rep(NA, ((noOfCols-1) * (noOfRows-1))  * nBoot ),
                             nrow=(noOfCols-1) * (noOfRows-1), ncol=nBoot)
  for (bootCntr in (1:nBoot) )
  {
    counter=1;
    for(i in  (1: length(drugXObsMean)) ) #xDrug Mean  i=1
    {
      for (j in (1:length(drugYObsMean))) #yDrug mean j=1 
      {
        xBar <- 0
        yBar <- 0
        xyBar <- 0
        
        
        currYind <- j; 
        currXind <- drugXInds[i];
        currXYind <- matForInd[j, i+1] 
          
        # sample Ea, Eb and Eab
        residualsEa <- sample(residuesWindowsList[[currXind]] , 1, replace=TRUE)
        residualsEb <- sample(residuesWindowsList[[currYind]] , 1, replace=TRUE)
        residualsEab <- sample(residuesWindowsList[[currXYind]] , 1, replace=TRUE)
        
        xBar <- drugXObsMean[i] * residualsEb
        yBar <- drugYObsMean[j] * residualsEa
        
        # according to equation of bliss
        synergy_Best_Bar[counter, bootCntr] <-   (xBar+yBar+(residualsEa * residualsEb) ) - residualsEab

        counter <- counter + 1;
      }
    }
  }
  
  
  quantilesCI_ConcComb <- matrix(rep(0, 3*length(combXYObsMeanVector)), 3, length(combXYObsMeanVector) )
  for (pv in (1:length(combXYObsMeanVector)) )
  {
    #0.025= Lower Bound and 0.975 is Upper bound
    quantilesCI_ConcComb[2:3, pv] <- apply(matrix(synergy_Best_Bar[pv, ], 1, nBoot) , 1, quantile, c(.025, 0.975  ))
  }  
  
  quantilesCI_ConcComb[1,] <- synergy_Bliss_obsVector
  #Maximum max(synergy_Bliss_obsVector)
  quantilesCI_Max <-    apply(matrix(apply(synergy_Best_Bar, 2, max), 1, nBoot), 1, quantile, c(.025, 0.975  ))
  #Minimum min(synergy_Bliss_obsVector)
  quantilesCI_Min <-    apply(matrix(apply(synergy_Best_Bar, 2, min), 1, nBoot) , 1, quantile, c(.025, 0.975  ))
  
  return( list(quantilesCI_ConcComb, quantilesCI_Max, quantilesCI_Min))
}  


 # Non Heteroscedastic residual selection based on sliding window of 10% nearest neighbours.
 residualSelection <- function(totRepl, nr, nc, rawDataPPNA )
 {
   #Start of Sliding code
   # to make a function follwing arguments are required
   # 1-totalNumberofReplicates*noOfRows* noOfCols
   # 2-rawDataPreProcessed_NA
   # and output will be the residuesWindowsList as per order  of rawDataPreProcessed_NA.
   
   # Extract all residues and store them in myResidualMat, where last row of it 
   # contains the respective means.
   myResidualMat <-  matrix(rep(NA, (totRepl + 1)*nr* nc),
                            nrow = totRepl + 1, ncol = nr* nc)
   rawDataPPNAMean <- apply(rawDataPPNA, 2, mean, na.rm=TRUE) 
   
   for(k in (1:(nr * nc) )){ #k=1
     myResidualMat[ (1: (nrow(myResidualMat)-1) ) ,k] <- rawDataPPNA[,k]- 
       rawDataPPNAMean[k]
     myResidualMat[ nrow(myResidualMat), k ] <- rawDataPPNAMean[k]
   }
   
   
   # Make a sliding window
   # 10% of nearest neighbours
   nNeig= ceiling(0.1 * ncol(myResidualMat))
   # store sliding residual windows j wise
   residuesWindowsL <- list( );
   
   for (j in 1:ncol(myResidualMat))
   {
     currentMean <-  myResidualMat[ nrow(myResidualMat),  j] 
     diffCurrent <- (  abs(myResidualMat[ nrow(myResidualMat), ] - currentMean)    )
     sortedAbsDiffs <- sort(diffCurrent)
     selectedWindowDiffs <- sortedAbsDiffs[2:(k + 1) ]
     indicesofWindows <- 0
     for(cntr in 1:nNeig )
     {
       indicesofWindows[cntr] <- which( selectedWindowDiffs[cntr] == diffCurrent  )
     }
     # Add the index of current combination also
     indicesofWindows[cntr+1] <- j
     # Extract sliding residues
     residuesWindows <- 0
     # residual window indices
     residuesWindows <- myResidualMat[ 1: (nrow(myResidualMat)-1), indicesofWindows  ] 
     # Clean for NAs
     residuesWindows_cleaned <- residuesWindows[!is.na(residuesWindows)]
     # make residue window symmetrical
     residuesWindowsL[[j]] <-  c( c(residuesWindows_cleaned), -1 * c(residuesWindows_cleaned)  )
     
   }
   # End of Sliding window
   return(residuesWindowsL)
 }

#' Combine data from multiple files
#' @param yConcentration Y drug concentrations 
#' @param xConcentration X drug concentrations
#' @param replNo Number of Replicates in all files
#' @param file  File name 
#' @param totalNumberofReplicates Total number of replicates per files
#' @param siReplicates data
#' @return Combined data of replicate survival indices from multiple experiments 
#' @examples 
#' xConc <- c(0.00,  0.20,  0.39,  0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50) 
#' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
#' rN <- 4
#' fN <- 1
#' trN <- 4
#' dataFile <- system.file("extdata", "rawDataPreProcessed.csv", package="COMBIA")
#' dataSample <- read.csv(dataFile, header=FALSE )
#' replList <- list(vector, 4)
#' for( i in 1:4)
#' { replList[[i]] <- dataSample[i,] }
#' rslt <- combineDataFromMultipleFiles(list(yConc), 
#' list(xConc), rN,fN,trN, replList )
#' @author Muhammad kashif
#' @export
combineDataFromMultipleFiles <- function(yConcentration, xConcentration, replNo, 
                                        file, totalNumberofReplicates, siReplicates )
   {
    # Make a big matrix that contain all data
  
    allY  <-   unique(as.vector(sapply(yConcentration, function(x){as.numeric(x)})))
    allY  <-  allY[order(allY, decreasing=TRUE)] # keep y cons in proper order, large to small
    
    allX  <-  unique(as.vector(sapply(xConcentration, function(x){as.numeric(x)})))
    allX  <-  allX[order(allX, decreasing=FALSE)] # keep x cons in proper order, small to large
    bigDataList <- list("vector", replNo-1)
    bigMatrix <- matrix(rep(0,(length(allX)*length(allY))), nrow= length(allY), ncol=length(allX) )
    # Add data of all replicates of all files at proper location in bigMatrix
    replNoNew <- 0
    for (noOfFiles in (1:length(file) )){
        curr_xConcentration <- as.vector(xConcentration[[noOfFiles]])
        curr_yConcentration <- as.vector(yConcentration[[noOfFiles]])
        for (repl in (1:totalNumberofReplicates[noOfFiles] )){ #repl=1
          replNoNew <- replNoNew + 1;
          # Convert SiData intomatrix column wise
          currentSIData <- matrix(siReplicates[[replNoNew]], nrow=length( unlist(yConcentration[[noOfFiles]] ) ) , ncol=length(xConcentration[[noOfFiles]]) )
         
          for (i in (1:length(allX)) ){
            for (j in (1:length(allY)) ){
                xIndex <- which(curr_xConcentration==allX[i])
                yIndex <- which(curr_yConcentration==allY[j]) 
                if( (length(xIndex)==0) ||(length(yIndex)==0)   ) 
                    {
                      bigMatrix[j, i] <- 0
                    }else{
                      bigMatrix[j, i] <- as.numeric(currentSIData[yIndex, xIndex])
                    }# end of descision structure
            }# end of allY
          }# end of allX
        bigDataList[[replNoNew]] <- bigMatrix
        bigMatrix <- matrix(rep(0,(length(allX)*length(allY))), nrow = length(allY), ncol = length(allX) )
      }# end of replciates in a file
    }# end of files
  
    arrangeReplicates <- matrix(rep(0,sum(totalNumberofReplicates) * length(allX) *length(allY) ),
                            nrow=sum(totalNumberofReplicates), ncol=length(allY)*length(allX))
    for ( i in 1:sum(totalNumberofReplicates)){#i=1
        arrangeReplicates[i,] <- as.numeric(bigDataList[[i]])
        }
    
  return(arrangeReplicates)
  }# End of combineDataFromMultipleFiles



#' Read data from macsynergyII formate and clean for outliers
#' @param file Name of fiele to be read
#' @param sheet Sheet Number
#' @param nrow Number of rows in the sheet
#' @param wellRangesExcel TRUE if wells in excel formate
#' @param minThersholdForCVCal Thresolld for data outliears in CV
#' @param minThersholdForCV Thresold of values in CV not to remove
#' @param survivalFunc <-  function (x,y,z) {(x-z)/(y-z)} # It can be any function
#' @return Matrix of replicated values
#' @examples
#' fl <- system.file("extdata", "testData.csv", package="COMBIA")
#' sh <- 1
#' wellR <- list(c( "l3:l10","m3:m10","b3:k10",  "l13:l20","m13:m20","b13:k20", 
#'            "l23:l30","m23:m30","b23:k30",  "l33:l40","m33:m40","b33:k40"))
#' minThersholdForCV <- 0.3
#' minThersholdForCVCal <- 0.1
#' survivalFunc <-  function (x,y,z) {(x-z)/(y-z)}
#' rslt <- readMacSynergyValues(fl, sh, nrow=41, wellR,  
#' minThersholdForCVCal, minThersholdForCV, survivalFunc)
#' @author Muhammad kashif
#' @export
readMacSynergyValues <- function(file, sheet, nrow=41, wellRangesExcel,
                                minThersholdForCVCal, minThersholdForCV, survivalFunc){
  # It will be vector that conatianing total number of replciates per file.
    totalNumberofReplicates <- as.numeric(lapply(wellRangesExcel, length))/3 
    siReplicates <- vector("list",sum(totalNumberofReplicates))
    replNo <- 1
    yConcentration <- vector("list", length(totalNumberofReplicates))
    xConcentration <- vector("list", length(totalNumberofReplicates))
    for(cn in (1:length(file)) )
      {
      
      # extract extension of file from file path
      filenamechunks <- unlist( strsplit( file[cn], ".", fixed = TRUE))
      if ( (filenamechunks[length(filenamechunks)] == "xls") | (filenamechunks[length(filenamechunks)] == "xlsx") ){
        #print("XLSSSS")
        #library(gdata)
        plateData <- read.xls(  file[cn], sheet=sheet, skip=0, sep=",", nrows=41, fill=TRUE, header=FALSE, blank.lines.skip = FALSE)[,1:13]
      } else if ( filenamechunks[length(filenamechunks)] == "csv"){
        #print("CSV")
        plateData <- read.csv(  file[cn], skip=0, sep=",", nrows=41, fill=TRUE, header=FALSE, blank.lines.skip = FALSE)[,1:13]
      } else {
        #print("TAB")
        plateData <- read.table(  file[cn], skip=0, sep=",", nrows=41, fill=TRUE, header=FALSE, blank.lines.skip = FALSE)[,1:13]
      } 
      
        # Extracct replicate data for every control, case and empty well
        # create variables to hold data
        controlValues <- vector("list", totalNumberofReplicates[[cn]])
        emptyValues <-  vector("list", totalNumberofReplicates[[cn]])
        caseValues <-    vector("list", totalNumberofReplicates[[cn]]) 
        yConcentration[[cn]] <- round(as.numeric(as.matrix(plateData[3:10, 1])),2)
        xConcentration[[cn]] <- round(as.numeric(as.matrix(plateData[11, 2:11])),2)  
        controlValues <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=1, excelFormate=TRUE) # extract control wells
        emptyValues   <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=2, excelFormate=TRUE) # extract empty wells
        caseValues    <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=3, excelFormate=TRUE) # extract case wells
        sur <- survivalFunc
        for (macI in (1:totalNumberofReplicates[cn]) )
            { print(paste("Ratio between empty and control",  ( mean( as.numeric(as.vector(controlValues[[macI]])), na.rm=TRUE  ) /mean( as.numeric(as.vector(emptyValues[[macI]])), na.rm=TRUE)  ) ))
              print( paste(   paste(paste("CV for control:", macI  ), ":")  ,  100 * (sd( as.numeric(as.vector(  controlValues[[macI]]) ), na.rm=TRUE )/mean( as.numeric(as.vector(controlValues[[macI]] )), na.rm=TRUE) ) ) )
              # Calculate Sis
              # A matrix is decomposed by column wise
              siReplicates[[replNo]] <- lapply(as.numeric(as.matrix(caseValues[[macI]])),
                                      sur, y=mean( as.numeric(as.vector(controlValues[[macI]]) ), na.rm=TRUE),
                                      z=mean( as.numeric(as.vector(emptyValues[[macI]])), na.rm=TRUE)  )
              replNo <- replNo + 1;
            }  
      }# loop no of files
  
    
    
    arrangeReplicates <- combineDataFromMultipleFiles(yConcentration, xConcentration, replNo, 
                                        file, totalNumberofReplicates, siReplicates )
    arrangeReplicatesSave <- arrangeReplicates
    arrangeReplicatesSave[which(arrangeReplicatesSave==0)] <- NA 
    arrangeReplicatesSaver <- rbind(arrangeReplicatesSave, apply(arrangeReplicatesSave, 2, cVCal))
     
   # Print No of datapoints with CV > minThresholdForCV
    print(paste("Number of the variable datapoints before data removal=", 
                length(which(arrangeReplicatesSaver[nrow(arrangeReplicatesSaver), ] > minThersholdForCV)  ) ))
    # Function to remove Outliers:
    rawDataPreProcessed <- removeOutliers(arrangeReplicates, minThersholdForCVCal, minThersholdForCV )
    # Add cv in the last column
    rawDataPreProcessedSave <- rawDataPreProcessed
    rawDataPreProcessedSave[which(rawDataPreProcessedSave==0)] <- NA 
    rawDataPreProcessedSaver <- rbind(rawDataPreProcessedSave,apply(rawDataPreProcessedSave, 2, cVCal))
    
    print(rawDataPreProcessed)
    # Print No of datapoints with CV > minThresholdForCV after correction
    print(paste("Number of the variable datapoints after data removal=", 
                length(which(rawDataPreProcessedSaver[nrow(rawDataPreProcessedSaver), ] > minThersholdForCV)  ) ))
    
    
  return(rawDataPreProcessed);
  } # End of MacSynergy function




#' Read data from raw FMCA format and clean for outliers
#' @param file Name of file to be read
#' @param platetype 384 etc 
#' @param keyposition Bar code position
#' @param selectionkey 65000
#' @param platekey Barcode
#' @param wells Wells ranges
#' @param minThersholdForCVCal Thresolld for data outliears in CV
#' @param minThersholdForCV Thresold of values in CV not to remove
#' @param yConcentration Concentrations of y drug
#' @param xConcentration Concentrations of x drug
#' @return Matrix of replicated survival values 
#' @examples 
#' fl <- system.file("extdata","FluoOptima_384_2014-03-28test.txt", package="COMBIA")
#' wls <- list(c("A11:H11", "A12:H12","A1:H10",   "I11:P11", "I12:P12","I1:P10", 
#'         "A23:H23", "A24:H24","A13:H22",   "I23:P23", "I24:P24","I13:P22")
#'                         )
#' pltype <- "384"
#' keypos <- 2     
#' seleckey <- "65000"
#' barCode <- 7049
#' minThersholdForCVCal <- 0.1 
#' minThersholdForCV <- 0.3
#' xConc <- c(0.00,  0.20,  0.39,0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50.00) 
#' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
#' readFMCAValues(fl, pltype, keypos, seleckey, barCode,
#'               wls, minThersholdForCVCal, minThersholdForCV, xConc, yConc   )
#' @author Muhammad kashif
#' @export
readFMCAValues <- function(file, platetype, keyposition,      
                          selectionkey, platekey, wells,
                          minThersholdForCVCal,
                          minThersholdForCV, 
                          yConcentration,
                          xConcentration
                          ){
  
  # Extract no of replicates
  # It will be vector that conatianing total number of replciates per file.
  totalNumberofReplicates <- as.numeric(lapply(wells, length))/3 
  siReplicates <- vector("list", sum(totalNumberofReplicates))
  replNo <- 1
  for(cn in (1:length(file)) )
    {
    # Call to function that read raw fmca data and convert  them into survival values
     rawDataUnProcessed <- readFluostarPlates(filename = unlist(file[cn]), platetype = platetype[cn], keyposition = keyposition,      
                                                  selectionkey = selectionkey, platekey = platekey[cn], wells = wells[[cn]]                  
                                            )
     replicateValues <- extractReplicateValues(rawDataUnProcessed, wells[[cn]], wellplace=3, excelFormate=FALSE) # extract control wells
     for(cm in 1:length(replicateValues) )
        {
          siReplicates[[(replNo) ]] <- replicateValues[[cm]]
          replNo <- replNo + 1;
        }
    }# loop no of files

  arrangeReplicates <- combineDataFromMultipleFiles(list(yConcentration), list(xConcentration),replNo, 
                                                    file,totalNumberofReplicates, siReplicates
                                                    )
  # Function to remove Outliers:
  rawDataPreProcessed <- removeOutliers(arrangeReplicates,  minThersholdForCVCal, minThersholdForCV
                                        )
  
  return(rawDataPreProcessed);
  } # End of READ fmca FUNCTION



#' Read data from raw format and clean for outliers
#' @param file Name of fiele to be read
#' @param sheet Sheet 
#' @param rskip Number of rows to skip before reading data, default rskip=0 
#' @param cStart Number of column to start reading data, default cStart=1
#' @param wellRangesExcel well ranges in excel formate
#' @param platetype 384 or 96 
#' @param minThersholdForCVCal Thresolld for data outliears in CV
#' @param minThersholdForCV Thresold of values in CV not to remove
#' @param survivalFunc A function to calculate survival values
#' @param xConcentration Concentrations of drug at x-axis
#' @param yConcentration Concentrations of drugs at y-axis
#' @return Matrix of survival values of experimental replicates
#' @examples 
#' fl <- system.file("extdata", "FluoOptima_384_2014-03-28test.csv", package="COMBIA")
#' wls <- list(  c(  "K1:K8", "L1:L8","A1:J8",     "K9:K16", "L9:L16","A9:J16", 
#'                   "W1:W8", "X1:X8","M1:V8",     "W9:W16", "X9:X16","M9:V16")
#'                   )
#' sh <- 1
#' rskip <- 0 
#' cStart <- 1
#' pltype <- "384"
#' minThersholdForCVCal <- 0.1
#' minThersholdForCV<- 0.3
#' survivalFunc <- function (x,y,z) {(x-z)/(y-z)}
#' xConc <- c(0.00,  0.20,  0.39,  0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50.00) 
#' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
#' rslt <- readOtherValues(fl, sh, rskip, cStart, wls, pltype, minThersholdForCVCal, 
#'                 minThersholdForCV, survivalFunc, xConc, yConc )
#' @author Muhammad kashif
#' @export
readOtherValues <- function(file, sheet, rskip=0, cStart=1, wellRangesExcel, platetype,
                            minThersholdForCVCal, minThersholdForCV,  survivalFunc,
                           xConcentration, yConcentration)
                          {
      if (platetype == "384"){
          rowsPerPlate <- 16; colsPerPlate <- 24  
          } else{
          rowsPerPlate <- 8; colsPerPlate <- 12 
        }
  
  
    # Extract no of replicates
    # It will be vector that conatianing total number of replicates per file.
    totalNumberofReplicates <- as.numeric(lapply(wellRangesExcel, length))/3 
    siReplicates <- vector("list", sum(totalNumberofReplicates))
    replNo <- 1
    for(cn in (1:length(file)) )
        {
          # In raw format
          
      filenamechunks <- unlist( strsplit( unlist(file[cn]), ".", fixed = TRUE))
      if ( (filenamechunks[length(filenamechunks)] == "xls") | (filenamechunks[length(filenamechunks)] == "xlsx") ){ 
        #library(gdata)
          plateData <- read.xls( unlist(file[cn]), sheet=sheet, skip=rskip, sep=",", nrows=rowsPerPlate, fill=TRUE, header=FALSE, blank.lines.skip = FALSE)[ ,cStart:( (cStart-1) + colsPerPlate)]
      } else if ( filenamechunks[length(filenamechunks)] == "csv"){  
        plateData <- read.csv( unlist(file[cn]), skip=rskip, sep=",", nrows=rowsPerPlate, fill=TRUE, header=FALSE, blank.lines.skip = FALSE)[ ,cStart:( (cStart-1) + colsPerPlate)]
        
        } else {  
          plateData <- read.table( unlist(file[cn]), skip=rskip, sep=",", nrows=rowsPerPlate, fill=TRUE, header=FALSE, blank.lines.skip = FALSE)[ ,cStart:( (cStart-1) + colsPerPlate)]
      }
            # Extract replicate data for every control, case and empty well
          # Create variables to hold data
          controlValues <- vector("list", totalNumberofReplicates[[cn]])
          emptyValues <-   vector("list", totalNumberofReplicates[[cn]])
          caseValues <-    vector("list", totalNumberofReplicates[[cn]]) 
    
          controlValues <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=1, excelFormate=TRUE) # extract control wells
          emptyValues   <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=2, excelFormate=TRUE) # extract empty wells
          caseValues    <-  extractReplicateValues(plateData, wellRangesExcel[[cn]], wellplace=3, excelFormate=TRUE) # extract case wells
    
          for (macI in (1:totalNumberofReplicates[cn]) )
              { 
                
                print(paste("Ratio between empty and control", ( mean( as.integer(as.vector(controlValues[[macI]])), na.rm=TRUE  ) / mean( as.integer(as.vector(emptyValues[[macI]])), na.rm=TRUE)   ) ))
                print( paste(   paste(paste("CV for control:", macI  ), ":"),  100 * (sd( as.integer(as.vector(  controlValues[[macI]]) ), na.rm=TRUE )/mean( as.integer(as.vector(controlValues[[macI]] )), na.rm=TRUE) ) ) )
      
                # Calculate survial valures
                # A matrix is decomposed by column wise
                siReplicates[[replNo]] <- lapply(as.numeric(as.matrix(caseValues[[macI]])),
                                      survivalFunc, y=mean( as.integer(as.vector(controlValues[[macI]]) )), z=mean( as.integer(as.vector(emptyValues[[macI]])), na.rm=TRUE)  )
                
                replNo <- replNo + 1;
                
              }  
    
        }# loop no of files
  
  
      arrangeReplicates <- combineDataFromMultipleFiles(list(yConcentration), list(xConcentration), replNo, 
                                                   file, totalNumberofReplicates, siReplicates )
     
    # Function to remove Outliers:
      rawDataPreProcessed <- removeOutliers(arrangeReplicates,
                                             minThersholdForCVCal, minThersholdForCV
                                          )
      return(rawDataPreProcessed);
    } # End of other function



#' This function calculates significant synergy/antagonism according to Bliss or Loewe model 
#' and creates scientific publication ready graphs.    
#' @param filename Name of file containing experimental data.
#' For MS Excel files, working version of Perl must be present in the executable search path.
#' @param sheet Optional, sheet number if excel file is used for input.
#' @param model  bliss or loewe.
#' @param inputFormates Any of these three formates "fmca", "macsynergy" and "others" are supported.
#' Example is provided with macsynergy format and test data for this example can be found in 
#' installation directory ("extdata")  of COMBIA.  See files FluoOptima_384_2014-03-28test_M and testDataM in 
#' directory "extdata" for format details of "fmca" and "macsynergy". "others" can be any other format.
#' @param platetype Optional default is 384. Only 384 and 96 well plates are supported. 
#' @param keyposition Optional default is 2. Usefull for automated barcoded data.
#' @param selectionkey Optional default is 65000.
#' @param platekey Optional barcode.
#' @param minThersholdForCVCal Optional default is 0.15. 
#' @param minThersholdForCV Optional default is 0.3.
#' @param wells  wells argument should  be in triplet form that is 
#' 1-Untreated control wells range, 2-empty wells range and 3-case wells range.
#' Thus in example below (see well argument) experiment has four replicates.
#' "l3:l10","m3:m10","b3:k10" is first replicate. Where "l3:l10" is the location
#' of untreated control values in the testData.csv, "m3:m10" is the background/
#' empty well well values and "b3:k10" are values after treatment. 
#' @param yConcentration Y drug Concentrations. 
#' @param xConcentration X drug Concentrations.
#' @param xDrug X drug name.
#' @param yDrug Y drug name.
#' @param cellLine Cell/Experiment name.
#' @param survivalFunc Optional default is function (x,y,z) {(x-z)/(y-z)} 
#' i.e (treated - background)/ (untreated - background).
#' @param nBoot Optional Number of time to bootstrap default is 5000  
#' @return Stores and show graph/data of synergy/antagonism analyses
#' @examples
#' fl <- system.file("extdata", "testData.csv", package="COMBIA")
#' wellR <- list(c("l3:l10","m3:m10","b3:k10", "l13:l20","m13:m20","b13:k20", 
#'            "l23:l30","m23:m30","b23:k30", "l33:l40","m33:m40","b33:k40") )
#' mdl <- "bliss"
#' xConc <- c(0.00,  0.20,  0.39,  0.78,  1.56,  3.12,  6.25, 12.50, 25.00, 50) 
#' yConc <- c(128,  64,  32,  16,   8,   4,   2,   0)
#' xDrug <- "A"
#' yDrug <- "B"
#' cellLine <-"Cell"
#' analyzeCOMBO(filename = c(fl), model = "bliss", inputFormates = "macsynergy", 
#'                  wells = wellR, yConcentration = yConc, xConcentration = xConc,
#'                  xDrug = xDrug, yDrug=yDrug, cellLine = cellLine, nBoot=500)
#'
#' @author Muhammad kashif
#' @export
analyzeCOMBO <- function( filename, sheet=1, model, inputFormates, platetype="384", keyposition = 2,
                              selectionkey="65000", platekey=7051, minThersholdForCVCal=0.15, minThersholdForCV=0.3,
                              wells, yConcentration,xConcentration,
                              xDrug, yDrug,cellLine,     survivalFunc= function (x,y,z) {(x-z)/(y-z)}, nBoot=5000)
  {
    # Local variable declaration
    replicateValues <- 0;
    rawDataPreProcessed <- 0;  
    noOfRows <- 0   #  Variable store No of Rows to be used to formate data from matrix to linear and vice versa 
    noOfCols <- 0   #  Variable store No of Cols to be used to formate data from matrix to linear and vice versa
    inputFileFormates <- inputFormates
    # Extract number of rows and columns of an experiment
    noOfRows <- length(yConcentration)  # Calcualte No of rows in a single data
    noOfCols <- length(xConcentration)  # Calculate No of cols in a single data
  
    # Apply function and calculate values
    if (inputFileFormates == "fmca")
      {
        rawDataPreProcessed <- readFMCAValues( file=filename, platetype=platetype, keyposition = keyposition,      
                                                selectionkey = selectionkey, platekey= platekey, wells= wells,   
                                                minThersholdForCVCal = minThersholdForCVCal,
                                                minThersholdForCV = minThersholdForCV,  yConcentration, xConcentration
                                            )
      } else if(inputFileFormates == "macsynergy"){
          # wells in excel format
          wellRangesExcel <- wells
          rawDataPreProcessed <- readMacSynergyValues(filename, sheet, nrow=41, wellRangesExcel,
                                                      minThersholdForCVCal, minThersholdForCV, survivalFunc)
       
          
      } else if(inputFileFormates == "others"){
          # wells in excel format
          wellRangesExcel <- wells
          rawDataPreProcessed <- readOtherValues( file=filename, sheet,  wellRangesExcel, platetype,
                                            minThersholdForCVCal, minThersholdForCV,  
                                            survivalFunc, yConcentration=yConcentration, xConcentration=xConcentration )
       
      }

    
   if (model=="bliss")
      { 
      
        # synergyBlissCalculation consists of three lists are resturned one is
        # first row of raw Bliss Synergy/Antagonism
        # 2nd row lower bound of BootStrap Interval
        # 3rd row upper bound of BootStrap Interval
        # 2nd list is CI of the maximum significant concentration combination 
        # 3rd list is CI of the minimum significant concentration combination 
        print("Applying Bliss model")
        synergyBlissCalculationLists <-  applyBliss(noOfRows, noOfCols, rawDataPreProcessed, nBoot)
    
        # Process lists and save them in user readale formate 
        processedData <- synergySignificant(synergyBlissCalculationLists, noOfRows, noOfCols, xDrug, yDrug, cellLine )
        # Synergy/Analysis plot 
        synAntPlot(processedData, xConcentration, yConcentration, xDrug, yDrug, cellLine)
    
      } else if(model=="loewe"){   
          # synergyBlissCalculation consists of three lists are resturned one is
          # first row of raw Bliss Synergy/Antagonism
          # 2nd row lower bound of BootStrap Interval
          # 3rd row upper bound of BootStrap Interval
          # 2nd list is CI of the maximum significant concentration combination 
          # 3rd list is CI of the maximum significant concentration combination 
          synergyLoeweCalculationLists <- applyLoewe(rawDataPreProcessed, xConcentration, yConcentration, nBoot)
          # Process lists and save them in user readale formate 
          processedData_Loewe <- synergySignificant(synergyLoeweCalculationLists, noOfRows, noOfCols, xDrug, yDrug, cellLine )
          # Synergy/Analysis plot 
          synAntPlot(processedData_Loewe, xConcentration, yConcentration, xDrug, yDrug, cellLine )
      }
  } # End of analyzeCOMBO Main function



################################Floustar data reading##############################################
# This function can read data from thefile that may or may not be generated by automated 
# wet lab experimental plateforms, e-g fluostar, automated  robotics etc. 
# It calculates survival indices of the specified wells, you only need to call "readFluostarPlates" 
# function with proper arguments to execute everything (Example are provided with sample values).
# It is also important to note that wells argument of readFluostarPlates should always be in triplet form that is 
# 1-control wells range, 2-empty wells range and 3-case wells range.
# It can also read the data from the file where a plate is read only one time, still it cope with variations if an experiment is
# repeated twice or many time in adjacent rows in the file. Another flexibilty of is its ability to calculate S.I of those
# experiments in which single plate and no repeated row is used. 
# 
#
####################################################################################################

#' Reads experimental data from a file.
#' This function reads the data from specified (excel,log, txt etc) file and store it in a data frame. 
#' @param filename Filename.ext.
#' @param separator Any character(, ; ' etc) that is used as a separator in specified file.
#' @param sheet Need to use only when reading excel files. It is the number of the excel sheet to be read in a worksheet.
#' @param noofrows_skip Number of the rows in the file that should be skipped before starting the data reading.
#' @param readplates Number of the plates that you want to read from a set of plates in a file.This parameter can only 
#' be used with excel files. Otherwise it will be ignored.
#' @param numberofrowsperplate It is calculated on the basis of type of plates i-e number of rows per plates are 17 for 
#' 384 well plates(16 lines from plates + 1 header lines) and 9 for 96 well plates (8 lines from plates + 1 header lines).
#' @param platetype type of plate used i-e 384 or 96 well plate. 
#' @return Data frame of file data.
#' @examples 
#' f <- system.file("extdata", "optima.log", package="COMBIA")
#' fileDF <- readFile(filename = f, separator ="\t", sheet=1, noofrows_skip=0,   
#' readplates=1, numberofrowsperplate=17, platetype="384")  
#' @author Muhammad Kashif
#' @export
readFile <- function(filename, separator, sheet, noofrows_skip, readplates, numberofrowsperplate, platetype)
  {
    # Change for Version3 Date 2011-05-27 
    # an argument platetype is added to specify the colulmns
    columnrange <- 0 
    if (platetype == "384"){
        columnrange <- 1:24
      } else if (platetype =="96"){
          columnrange <- 1:12
      }  
      
    # Change for Version2 Date 2011-02-17
    filenamechunks <- unlist( strsplit( filename, ".", fixed = TRUE))
    if ( (filenamechunks[length(filenamechunks)] == "xls") | (filenamechunks[length(filenamechunks)] == "xlsx") ){ 
        # If plate types are not specified
        if (platetype == ""){ 
          #library(gdata)
            rawdata <- read.xls( filename, sheet = sheet, skip = noofrows_skip, nrows = readplates * numberofrowsperplate, header = FALSE, fill = TRUE)
        } else { 
          #library(gdata)
            rawdata <- ( read.xls( filename, sheet = sheet, skip = noofrows_skip, nrows = readplates * numberofrowsperplate, header = FALSE, fill = TRUE)[,columnrange])
        }
        return(rawdata)   
    }else{
    
    # This part of code can read .log file and store in a data structure####
    
    # 1.classical parser was not written because we were interested in barcode only that is always at second position of 
    # header and it can be extracted with methods provided by R. I guess we dont have any well designed and fairly
    # complex grammer for parser. Therefore regular expressions and regular expression like structures are the logical choice.
    
    # 2. Second choice of storing plate data was, one vector and one matrix (As per concept of database tables), 
    # one(vector) for storing headers and 
    # other one for storing the data. As we already know each plate has a line of header immediately followed by 16 lines of 
    # data(following the geometry of plates), therefore, it is not needed to have a link between these vector and matric through 
    # pointer like structure because 
    # 1st header * 1 and header > 1 * 17 give the same thing.
    
    
    # 3. Still another way of reading the data from .txt or .log is to use the built in functions of R that are related to file  
    # reading and string operations: they are
    # srcfile() for file name
    # getsrclines() read source file line by line
    # unlist(strsplit()) to split the line 
    
    # But it is clear from the above functions that they will still need datastructure to store their results 
    # after some  # manipulations.
    
    # 4. Dataframe was the choice made because of there ability to store the numeric amd non numeric data. Here only one built in 
    # function (read.table) is sufficient to read and store data.
    
    # If plate types are not specified
    if (platetype == ""){
        rawdata <- data.frame(read.table( filename, sep = separator, fill = TRUE, skip = noofrows_skip))
      } else {  
        rawdata <- data.frame(read.table( filename, sep = separator, fill = TRUE, skip = noofrows_skip)[ , c(columnrange) ])
      }
      return(rawdata)
    }
  }

 
#' Extracts the keyvalues (Barcode) from a dataset, every plate needs barcode.
#' Keyvalues are extracted from the header of the plates at the position specified by keyposition argument.
#' @param keyposition Position of keyvalue in the header of plate.
#' @param rawdata An object(dataframe) of rawdata.
#' @param numberofrowsperplate  This argument is not needed when you call function "readFluostarPlates". The number of rows depend upon the
#' geometry of the plates. These are 16 in case of 384well paltes.
#' @param doubleplateexperiment This parameter can have TRUE & FALSE values only. It is set to TRUE when an experiment is performed
#' twice and we only want to choose only one of them. 
#' @return A complete set of keyvalues.
#' @examples
#' f <- system.file("extdata", "optima.log", package="COMBIA")
#' fileDF <- readFile(filename = f,  separator = "\t", noofrows_skip=0,
#'                    platetype="384")  
#' Generatedbarcode <- extractKey(keyposition = 2, rawdata = fileDF, 
#'                               numberofrowsperplate = 17, 
#'                               doubleplateexperiment = TRUE) 
#' @author Muhammad Kashif
#' @export
extractKey <- function(keyposition, rawdata, numberofrowsperplate, doubleplateexperiment){
              # Generation of the missing barcode is not possible mainly there are limitless
              # ways plates can be arranged on FMCA day.               
  numberofrows  <- nrow(rawdata)  # total number of rows
  rawbarcode    <- seq(length = (numberofrows / numberofrowsperplate), from = 0, to = 0)   # initializing vector of raww barcodes
  rawbarcode[1] <- as.vector(rawdata [1, keyposition] )  # raw barcode of first plate      
  
  cnt <- 1     			
  while (cnt <= ((numberofrows / numberofrowsperplate) - 1) )
    {  
      # loop to extract raw barcodes
      rawbarcode[cnt + 1] <- as.vector(rawdata[ (cnt * numberofrowsperplate) + 1, keyposition ])
      cnt <- cnt + 1 					
    }
  
    # Perform empty barcode check 
    if( any(rawbarcode == "NOREAD") )  
      {
        stop("Barcode is missing")
      }
  return(rawbarcode)
  }

#' Select one of the two read plates and built a hashtable.
#' One plate from each pair of the read plate is selected in case of double plate experinment on the basis of presence 
#' of minimum selection key and if none have maxed out values then one with highest mean value is picked.
#' @param rawdata An object(dataframe) of rawdata.
#' @param processedbarcode A vector of regenerated missing keyvalues. In this case it is the output of function "extractKey".
#' @param numberofrowsperplate  This argument is not needed when you call function "readFluostarPlates". The number of rows depends upon the
#' geometry of the plates. These are 16 in case of 384well paltes.
#' @param selectionkey keyvalue on basis of which a plate is slected from a pair of plates read in double plate experiment. 
#' @param doubleplateexperiment This parameter can have TRUE & FALSE values only. It is set to TRUE when an experiment is read twice.
#' @return A hashtable of picked plates.
#' @examples
#' f <- system.file("extdata", "optima.log", package="COMBIA")
#' fileDF <- readFile(filename = f, separator = "\t", noofrows_skip=0,
#'                    platetype = "384") 
#' Generatedbarcode <- Generatedbarcode <- extractKey(keyposition = 2,
#'     rawdata = fileDF, numberofrowsperplate = 17, doubleplateexperiment = TRUE) 
#' hashedplates <-  selectPlate(rawdata = fileDF,
#'     processedbarcode = Generatedbarcode, numberofrowsperplate=17,
#'     selectionkey="65000", doubleplateexperiment = TRUE  )
#' @author Muhammad Kashif
##' @export
selectPlate <- function(rawdata, processedbarcode, numberofrowsperplate, selectionkey, doubleplateexperiment)
  {
    startindex <- 0 # variable starting from zero 
    dataplate1 <- 0 # variable storing data of plate1
    dataplate2 <- 0 # variable storing data of plate2
  
    # Variables for hashing function
    barcodekey <- 0
    hashedplates <- hash(keys = unique(processedbarcode), values = seq(length = length(unique(processedbarcode)), 0, 0))
    while(startindex < nrow(rawdata) ){
    # Selecting one of the two consecutive plates ......
      dataplate1 <- as.integer(as.matrix(rawdata[(startindex + 2) : (startindex + numberofrowsperplate), ]))
      startindex <- startindex + numberofrowsperplate
      # Change for Version2 Date 2011-02-17 it is introduced here to control double and single plate experiments
      if (doubleplateexperiment == TRUE){
          dataplate2 <- as.integer(as.matrix(rawdata[(startindex + 2) : (startindex + numberofrowsperplate), ]))
          startindex <- startindex + numberofrowsperplate
            # Selection on the basis of number of 65000 values
            if (length(dataplate1[dataplate1 == selectionkey] ) == length(dataplate2[dataplate2 == selectionkey])){
                if(mean(dataplate1) > mean(dataplate2)){
                  # Plate1 selected
                  print("Plate 1 selected as per mean criteria")
                  barcodekey <- processedbarcode[(startindex / numberofrowsperplate)]
                  hashedplates[[as.character( barcodekey )]] <- dataplate1
                  } else{
                    # Plate 2 selected
                    print("Plate 2 selected as per mean criteria")
                    barcodekey <- processedbarcode[(startindex / numberofrowsperplate)]
                    hashedplates[[as.character( barcodekey )]] <- dataplate2
                  }  
            } else {     
                    if( length(dataplate1[dataplate1 == selectionkey]) < length(dataplate2[dataplate2 == selectionkey])){
                        print("Plate 1 selected as per least 65000 values")
                        barcodekey<- processedbarcode[(startindex / numberofrowsperplate)]
                        hashedplates[[as.character(barcodekey)]] <- dataplate1
                      } else{
                          # Plate2 selected 
                          print("Plate 2 selected as per least 65000 values")
                          barcodekey <- processedbarcode [(startindex / numberofrowsperplate)]
                          hashedplates[[as.character(barcodekey)]] <- dataplate2
                      }  
            }
        } else{ 
            # Change for Version2 Date 2011-02-17 #if single plate data is read
            # print("Plate 1 selected in single plate dataread")
            barcodekey <- processedbarcode[(startindex / numberofrowsperplate)]
            hashedplates[[as.character(barcodekey)]] <- dataplate1
        }
      }
    return(hashedplates)
  }



rangemean <- function(platebarcode, range, platetype, hashedplates, printcv )
  {
    # Start of function to extract range meanings, it works on the basis of plate labels and don't follow the excel style.
    plateforSI <- hashedplates[[as.character(platebarcode)]]
    if (platetype == "384"){
        dim(plateforSI) <- c(16,24)
      } else if (platetype == "96"){
      dim(plateforSI) <- c(8,12)
     }
  
    rangelist <- unlist(strsplit(range, ":"))
    rowstart <- substr(rangelist[1], 0, 1)
    columnstartindex <- as.integer(substr(rangelist[1], 2, 4))
    rowstartindex <- which( letters[1 : 26] == rowstart)
  
    rowend <- substr(rangelist[2 ], 0, 1)
    columnendindex <- as.integer(substr(rangelist[ 2 ], 2, 4 ))
    rowendindex <- which( letters[1:26] == rowend)
  
    welldata <- plateforSI[rowstartindex : rowendindex, columnstartindex : columnendindex]
    # FMCA quality checks
    # Relation between empty and control   >5
    # Control well CV <30%
    # Calculate Cv of wells and print it
  
    if (printcv==TRUE)
      {
        print(paste("CV% for control well is :", 100 * cVCal(welldata) ) )
      }
    return (mean( welldata ))
  }





#' Calculates survival indices (S.Is) for a range of wells (casewells).
#' S.Is for a range of wells are calculated, that range is specified at the third place of wells argument list. This function call the rangemean function 
#' to calculate the mean of the  range of the specified range. S.I is calculated by (Case well- meanofemptyrange/mean    
#' of controlwell- meanofemptyrange). In the wells argument one should provide arguments in the triplet form that is first one is control 
#' data range, second one is the empty data range  while third one is the control range.  
#' @param hashedplates A hash table of picked plates. It is the output of function "selectPlate".
#' @param platekey It is the key of the plate whose S.I is needed to be calculated.
#' @param platetype It is the type of plate (386 and 96). 
#' @param rowsperexperiment It is the argument that specifies if the same experiment is reptead and how many times in a plate. If an experiment is 
#' repeated twice in adjacent rows then average of its values will be used in the SI calculation. 
#' @param wells This argument can take a list of arguments in the triplet form. Where first argument of triplet is the range of control wells,      
#' second argument is the range of empty wells while third one is the range of case wells. It is made so that in labs plates layouts can differ  
#' greatly. By using this  triplet scheme one can handel a number of palte layouts.
#' @return A matrix with S.I showing values where they are actually exist on the plate.
#' @examples  
#' f <- system.file("extdata", "optima.log", package="COMBIA")
#' fileDF <- readFile(filename = f, separator = "\t", noofrows_skip=0,
#'                     platetype="384") 
#' Generatedbarcode <- extractKey(keyposition = 2,
#'                             rawdata = fileDF, numberofrowsperplate = 17, 
#'                             doubleplateexperiment=TRUE) 
#' hashedplates <-  selectPlate(rawdata = fileDF,
#'                             processedbarcode = Generatedbarcode, 
#'                             numberofrowsperplate = 17,
#'                             selectionkey = "65000", 
#'                             doubleplateexperiment = TRUE  )
#' survivalindeces <- calculateSi(hashedplates = hashedplates, 
#'                                 platekey = "7051", platetype = "384",rowsperexperiment=1,
#'                                 wells = c( "c8:h8","c1:n1","c3:c7",    "c8:h8","c1:n1","c9:c11", 
#'                                 "c8:h8","c1:n1","e3:e7",     "c8:h8","c1:n1","e9:e11",
#'                                 "c8:h8","c1:n1","g3:g7",     "c8:h8","c1:n1","g9:g11") 
#'                               )    
#' @author Muhammad Kashif
#' @export
calculateSi <- function(hashedplates, platekey, platetype, rowsperexperiment, wells )
  {
    # Change for Version2 Date 2011-02-17
    # three parameter cntrlrange, emptrange and caserange were removed and new "wells" is introduced to handel the 
    # variable number of these arguments.
    # rowsperexperiment argument stored the no of rows repeated per experiment 
  
    # Change for Version3 Date 2011-05-27 :: tolower function is added
    wellranges <- tolower(wells)
    if ( (length(wellranges)==0) | ( (length(wellranges)%%3)!=0 )){
        stop("Incorrect control, empty and case well ranges")
     }
  
    plateforSI <- 0 ## variable will store values of the plate underprocessing
    SI <- 0
    if (platetype == "384"){ ### this if statement will be helpful when we will be using same code for other plates may be 96 well
        plateforSI <- hashedplates[[as.character(platekey)]]  # extracting values of the plate based on barcode
        dim(plateforSI) <- c(16, 24)
        SI <- seq(length = 384, 0, 0)
        dim(SI) <- c(16, 24)
    } else if (platetype == "96" ){
          plateforSI <- hashedplates[[as.character(platekey)]]  # extracting values of the plate based on barcode
          dim(plateforSI) <- c(8, 12)
          SI <- seq(length = 96, 0, 0)
        dim(SI) <- c(8, 12)
    }  
    
      # Change for Version2 Date 2011-02-17
    cnt <- 0;
    cntrlMeanStore <- 0;
    emptyMeanStore <- 0;
    meanStoreCounter <- 0
    while (cnt < (length(wellranges)) ){
        # mean of the specified rangesControl and empty are calculated here
        cntrlmean <- rangemean(platekey, wellranges[cnt + 1 ], platetype, hashedplates, printcv=TRUE)
        print(paste("Controlmean:", cntrlmean))
        # Change 280514, print  CV of CONTROLS
        # FMCA Quality checks
        # Ration between empty and control   >5
        # Control well CV <30%
        meanStoreCounter <- meanStoreCounter + 1
        cntrlMeanStore[meanStoreCounter] <- cntrlmean
        emptmean  <- rangemean(platekey, wellranges[ cnt +2 ], platetype, hashedplates, printcv=FALSE)
        emptyMeanStore <- emptmean
        # These lines can be a part of the function but for clearity they are written seprate.
        range <- unlist(strsplit( wellranges[ cnt +3 ] ,":")) # split the range i-e from A2:D6 in to A2 and D6
        columnstartindex <- as.integer(substr(range[1], 2, 4)) # extract strating column
        rowstartindex <- which( letters[1 : 26] == substr(range[1], 0, 1)) #extract strating row
        columnendindex <- as.integer(substr(range[2], 2, 4))  # extract ending column
        rowendindex<- which( letters[1 : 26] == substr(range[2], 0, 1)) #extract ending row
          while(rowstartindex <= rowendindex)
            {
              if (rowsperexperiment > 1){
                SI[rowstartindex, columnstartindex : columnendindex] <- 
                    (( colMeans( matrix(plateforSI[rowstartindex: (rowstartindex + rowsperexperiment - 1),  columnstartindex : columnendindex],
                              nrow = length(rowstartindex: (rowstartindex + rowsperexperiment - 1)), ncol = length(columnstartindex : columnendindex) )) 
                              - emptmean) / (cntrlmean - emptmean))
              }else{
                SI[rowstartindex, columnstartindex : columnendindex] <- (  plateforSI[rowstartindex,  
                                                                              columnstartindex : columnendindex]  -                      
                                                                     emptmean) / (cntrlmean - emptmean)
              }
          rowstartindex <- rowstartindex + rowsperexperiment
          } # end of the while
      cnt <- cnt + 3 
    } # end of while loop of Change for Version2 Date 2011-02-17
  
    # FMCA Quality checks
    # Ratios between empty and control   >5
    # Control well CV <30%
    cntrlMeanStore;
    emptyMeanStore;
    print(paste("Ratio between empty and control=", 100* (emptyMeanStore/ mean(cntrlMeanStore))  ) )
    return(SI)
  }


#' Read a file and process it to calculate the Survival indeces(S.I).
#' This function calls other functions to complete its task. It reads a file to separate and regenerate the missing platekeys. 
#' Checks are performed to keep regenerated missing keyvalues in sync with data. It calculates survival indeces of the provided 
#' control wells, where wells should always be in triplet form that is  control well range, empty well range and case well range.
#' It can also handle the double plate experiments in which one plate is read twice and only one of them is selected in S.I calculations. 
#' Secondly it can also read the data from the file where a plate is read only one time, still it cope with variations if an experiment is
#' repeated twice or many time in adjacent rows in the file.
#' @param filename value of this argument should be path and filename.ext e=g "e:/optima.txt".
#' @param separator is the sepration character within the file assigned to filename.  
#' @param noofrows_skip Number of the rows in the file that should be skipped before starting the data reading.
#' @param sheet Need to use only when reading excel files. It is the number of the excel sheet to be read in a worksheet.
#' @param readplates Number of the plates to read from a set of plates from an excel file, This feature is only workable with xls files.
#' @param platetype Two types of plate formates are supported 384 and 96 wells. 
#' @param doubleplateexperiment This parameter can have TRUE & FALSE values only. It is set to TRUE when an experiment is read twice.
#' @param keyposition It is the position of key in the header. Currently it is located at the second position
#' but it can be at any position in the header.
#' @param selectionkey value, that will be used during the selection of plate. Current value is 65000.
#' @param platekey barcode of the plate whose wells you want to measure for Survival index    
#' @param rowsperexperiment It is the argument that specifies if the same experiment is repeated and how many times in a plate. If an experiment is 
#' repeated twice in adjacent rows then average of its values will be used in the SI calculation. 
#' @param wells This argument can take a list of arguments in the triplet form. Where first argument of triplet is the range of control wells,      
#' second argument is the range of empty wells while third one is the range of case wells. It is made so that in labs plates layouts can differ  
#' greatly. By using this  triplet scheme one can handel a number of palte layouts. Values should be given in the according to plate range e-g a4:d5   
#' means start from the  a(1) row and first column and continue to d(4) row 5th column. 
#' @return Matrix of S.I.
#' @examples
#' f <- system.file("extdata", "optima.log", package = "COMBIA")
#' platematrix <- readFluostarPlates(filename = f, platetype = "384", 
#'                                    keyposition=2, separator= "\t",     
#'                                    selectionkey = "65000", platekey = 7051,
#' wells = c( "c8:h8","c1:n1","c3:c7",    "c8:h8","c1:n1","c9:c11", 
#'            "c8:h8","c1:n1","e3:e7",     "c8:h8","c1:n1","e9:e11",
#'            "c8:h8","c1:n1","g3:g7",     "c8:h8","c1:n1","g9:g11" )
#'                                  )
#' @author Muhammad Kashif
#' @export  
readFluostarPlates <- function(filename, separator=",", noofrows_skip=0, sheet="1", readplates=1, platetype, doubleplateexperiment=TRUE, keyposition,    
                               selectionkey, platekey, rowsperexperiment=1, wells )
    {
      numberofrowsperplate <- 0  
      if (platetype == "384"){
          numberofrowsperplate <- 17
      } else if (platetype == "96") {
      numberofrowsperplate <- 9
    }
  
    rawdata <- readFile(filename, separator, sheet, noofrows_skip, readplates, numberofrowsperplate, platetype)
    Generatedbarcode <- extractKey(keyposition, rawdata, numberofrowsperplate, doubleplateexperiment)
    hashedplates <-  selectPlate(rawdata, processedbarcode = Generatedbarcode, numberofrowsperplate, selectionkey, doubleplateexperiment)
    survivalindeces <- calculateSi(hashedplates, platekey, platetype, rowsperexperiment, wells )    
  
  }





