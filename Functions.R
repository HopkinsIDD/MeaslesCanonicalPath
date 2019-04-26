

scl <- function(x){ (x - min(x,na.rm = T))/(max(x,na.rm = T) - min(x,na.rm = T)) }


####################################

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}



#' Main function to generate the canonical path data
#'
#' @param window.length - how many years is are the local cvs calculated over
#' @param regions - which WHO regions we consider
#' @param gaussian.st.dev - standarad deviation of the gaussian weighting function
#' @param cutoff - the number of years in the past at which we stop including data. Generally 
#' this is set at 50, meaning that all past data is used.
#' @param interp.resolution - how many segments we split each year into for animation purposes
#' @param year.shift.inc - how many years in the past are the gaussian weights centered?
#'

generate.data <- function(window.length, regions, 
                          gaussian.st.dev, cutoff = 50, 
                          interp.resolution = 20,
                          year.shift.inc = 3){
  
  list[subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data] = get.data.for.animation(regions)
  
  x = seq(1980, 2017) 
  
  ##' interpolate the datasets to have entries for all points in time once the interpolation is done.
  list[interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop] = interp.datasets(subset.data, 
                                                                                                      subset.vaccination, 
                                                                                                      subset.birth.rates, 
                                                                                                      subset.pop.by.year,
                                                                                                      x,
                                                                                                      x)
  
  ##' output matrices the correct size for our animation
  list[mean.cases, coeff.var.cases, incidence.per.1000, mean.br, mean.vac] = 
    prepare.matrices.for.animation(interp.subset.data, subset.data)
  
  ##' number of unique years that we will have data for. The longer the window, the less unique years of data.
  num.windows = length(x) - window.length + 1
  
  ##' first year of data
  year = 1980
  
  ##' setting up the datasets
  coeff.var = matrix(0, length(subset.data[ , 1]), num.windows)
  incidence.per.1000 = matrix(0, length(subset.data[ , 1]), num.windows)
  mean.br = matrix(0, length(subset.data[ , 1]), num.windows)
  mean.vac = matrix(0, length(subset.data[ , 1]), num.windows)
  
  ##' do calculations that calculate the coefficient of variation, incidence per 100, mean birth rate and
  ##' mean vaccination rate over periods of length given by the window length.
  for ( j in 1 : num.windows){
    for ( i in 1 : length(subset.data[ , 1])){
      coeff.var[i, j]  =  sd(interp.subset.data[i, paste(seq(year, year + window.length - 1))], na.rm = TRUE) /  
        mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      if(is.na( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)) == FALSE ){
        if( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) == 0) {
          coeff.var[i, j]  =  0
        } 
      }
      incidence.per.1000[i, j]  =  sum(as.numeric(interp.subset.data[i,  paste(seq(year, year + window.length - 1))]) / 
                                         as.numeric(interp.subset.pop[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) * 1000
      if(length(which(is.na(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))])))) < window.length){
        mean.br[i, j]  =  mean(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      }
      if(length(which(is.na(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))])))) < window.length){
        mean.vac[i, j]  =  mean(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      }
    }
    year = year + 1
  }
  
  incidence.per.1000.each.year = matrix(0, nrow(subset.data), length(seq(1980, 2017)))
  for(i in 1 : nrow(subset.data)){
    incidence.per.1000.each.year[i, ] = 1000 * as.numeric(interp.subset.data[i,  paste(seq(1980, 2017))]) / 
      as.numeric(interp.subset.pop[i, paste(seq(1980, 2017))])
    
  }
  ##' for calculating the weighted coefficient of variation, we need to take the weighted
  ##' mean of locally calculated coefficient of variations. The following function will
  ##' generate the matrix of weights used to calculate this weighted average.
  ##' the weights are gaussian, centred on a specific year, with
  ##' specified number of years for the standard deviation
  
  w1 = generate.cv.weights(coeff.var,
                           gaussian.st.dev,
                           cutoff)
  
  ##' for incidence, birth rate and vaccination rate, we can weight slightly differently
  ##' as there is no issue with the weighting of yearly calculated values.
  ##' the following function produces the weights used for weighting these variables
  
  w2 = generate.other.weights(d = incidence.per.1000.each.year,
                              window.length,
                              gaussian.st.dev,
                              cutoff, year.shift.inc)
  
  
  ##' make a set of matrices that are the same size as the matrices containing the data.
  coeff.2 = coeff.var
  incidence.2 = incidence.per.1000
  mbr2 = mean.br
  mvacc2 = mean.vac
  for(i in 1 : length(coeff.var[1, ])){
    for(j in 1 : length(coeff.var[, 1])){
      ##' make the entries of these newly created matrices to be the weighted averages of the originally calculated datasets
      coeff.2[j, i] = sum(coeff.var[j, ] * w1[i, ], na.rm = T)
      incidence.2[j, i] = sum(incidence.per.1000.each.year[j, ] * w2[i, ], na.rm = T)
      mbr2[j, i] =  sum(mean.br[j, i] * w2[i, ], na.rm = T)
      mvacc2[j, i] = sum(mean.vac[j, i]* w2[i, ], na.rm = T)
      
      ##' Should we do weighted average of birth rate and vaccination rate?
      ##' If so uncomment the next two lines
      
      #mbr2[j, i] = sum(mean.br[j, ] * w1[i, ], na.rm = T)
      #mvacc2[j, i] = sum(mean.vac[j, i] * w1[i, ], na.rm = T)
    }
  }
  
  ##' set the original data to be equal to the weighted data
  coeff.var.cases = coeff.2
  incidence.per.1000 = incidence.2
  mean.br = mbr2
  mean.vac = mvacc2
  
  
  ##' set up the timeline on which we do the interpolation. 
  ##' The number of sections that the yearly data is split up to is given by interp.resolution  
  x = seq(1980 + (window.length - 1), 2017)
  xout = seq(1980 + (window.length - 1), 2017, 1/interp.resolution)
  
  ##' interpolate the data and add columns that contain the corresponding country and WHO region of each line
  
  coeff.var.cases = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(coeff.var.cases, x, xout))
  incidence.per.1000 = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(incidence.per.1000, x, xout))
  mean.br = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.br, x, xout))
  mean.vac = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.vac, x, xout))
  
  ##' round the data to 2 decimal places for ease of reading.
  mean.vac[, -(1:2)] = as.numeric(mean.vac[, -(1:2)])
  mean.br[, -(1:2)] = as.numeric(mean.br[, -(1:2)])
  incidence.per.1000[, -(1:2)] = as.numeric(incidence.per.1000[, -(1:2)])
  coeff.var.cases[, -(1:2)] = as.numeric(coeff.var.cases[, -(1:2)])
  
  
  ##' set up the output to be the appropriate size and add column labels.
  ##' Additionally add enough room to include additional data for each year that will be used 
  ##' to calibrate the data for each year, so that the minimum and maximum of vaccination rate is 0 and 100 each time.
  ##' This ensures that the colour scale is constant
  
  output.data = matrix(0, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + (2 * length(regions) * length(xout)), 7)
  output.data  =  data.frame(output.data)
  colnames(output.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
  
  ##' input the appropriate data to the outputs 
  output.data$Country[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 1], length(coeff.var.cases[1, -(1:2)]))
  output.data$WHO_REGION[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 2], length(coeff.var.cases[1, -(1:2)]))
  count = 1
  for(i in 3 : length(coeff.var.cases[1, ])){
    output.data$Coefficient.of.Variation[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = coeff.var.cases[, i]
    output.data$Incidence[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = incidence.per.1000[, i]
    output.data$Mean.vaccination[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.vac[, i]
    output.data$Mean.birth.rate[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.br[, i]
    output.data$Year[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = xout[i - 2]  
    count = count + 1
  }
  
  ##' add in the dummy data for each year to keep the scales constant.
  year.mins = matrix(0, length(xout), 2)
  
  for(i in 1 : length(xout)){
    t  =  subset(output.data, output.data$Year ==  unique(xout)[i])
    year.mins[i, 1] = xout[i]
    year.mins[i, 2] = as.numeric(min(t$Mean.birth.rate,na.rm = T)  )
  }
  
  l =  expand.grid("", -1, 0, c(0,100), 0, xout, regions)
  
  for(i in 1 : (2 * length(regions) * length(xout))){
    y = l[i, 6]
    j = which(year.mins[, 1] == y)
    l[i, 5]  =  year.mins[j, 2]
  }
  
  output.data[((length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + 1 ): length(output.data[, 1]), ]  =  l
  
  for( i in 1 : (2 * length(regions) * length(xout))){
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] = regions[as.numeric(output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] )]
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)])  + i, 1] = ""
  }
  
  ##' make sure that each column that should be numeric is numeric.
  output.data$Coefficient.of.Variation = as.numeric(output.data$Coefficient.of.Variation)
  output.data$Incidence  =  as.numeric(output.data$Incidence)
  output.data$Mean.vaccination  =  as.numeric(output.data$Mean.vaccination)
  output.data$Mean.birth.rate  =  as.numeric(output.data$Mean.birth.rate)
  output.data$Year   =  as.numeric(output.data$Year)
  output.data$Coefficient.of.Variation[which(output.data$Coefficient.of.Variation == "Inf")] = 0
  return(output.data)
}

########################
### interpolate datasets for state space model shiny animation
########################


interp.datasets.state.space <- function(subset.data, 
                                        subset.vaccination, 
                                        subset.birth.rates, 
                                        subset.pop.by.year, 
                                        x,
                                        xout){
  
  interp.subset.data = matrix(0, length(subset.data[, 1]), length(xout) + 2)
  interp.subset.data[, 1] = as.character(subset.data$Country)
  interp.subset.data[, 2] = as.character(subset.data$WHO_REGION)
  
  interp.subset.vacc = matrix(0, length(subset.vaccination[, 1]), length(xout) + 2)
  interp.subset.vacc[, 1] = as.character(subset.data$Country)
  interp.subset.vacc[, 2] = as.character(subset.data$WHO_REGION)
  
  interp.subset.br = matrix(0, length(subset.birth.rates[, 1]), length(xout) + 2)
  interp.subset.br[, 1] = as.character(subset.data$Country)
  interp.subset.br[, 2] = as.character(subset.data$WHO_REGION)
  
  interp.subset.pop = matrix(0, length(subset.pop.by.year[, 1]), length(xout) + 2)
  interp.subset.pop[, 1] = as.character(subset.data$Country)
  interp.subset.pop[, 2] = as.character(subset.data$WHO_REGION)
  
  for ( i in 1 : length(subset.pop.by.year[, 1])){
    y = as.numeric(as.character(subset.data[i, paste(x)]))
    list[qq,ww] =  approx (as.numeric(x), as.numeric(y),  method = "linear", xout )
    interp.subset.data[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
    colnames(interp.subset.data) = c("Country", "WHO_REGION", paste(x))
    
    y1 = as.numeric(as.character(subset.vaccination[i, paste("X", x, sep = "")]))
    list[qq,ww] =  approx (as.numeric(x), as.numeric(y1),  method = "linear", xout )
    interp.subset.vacc[i, 3: length(interp.subset.vacc[1, ])]  =  round(ww, 2)
    colnames(interp.subset.vacc) = c("Country", "WHO_REGION", x)
    
    y2 =as.numeric( c(subset.birth.rates[i, paste("X", x, sep = "")]))
    if(length(y2) < length(y1)){
      y2 =as.numeric( c(subset.birth.rates[i, paste("X", x, sep = "")], subset.birth.rates[i, paste("X", 2012, sep = "")]))
    }
    if(length(which(is.na(y2) == FALSE)) < 2) 
    {interp.subset.br[i, 3: length(interp.subset.data[1, ])]  = 0} else{
      list[qq,ww] =  approx (as.numeric(x), as.numeric(y2),  method = "linear", xout )
      interp.subset.br[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
    }
    colnames(interp.subset.br) = c("Country", "WHO_REGION", x)
    
    
    y3 = subset.pop.by.year[i, paste("X", x, sep = "")]
    list[qq,ww] =  approx (as.numeric(x), as.numeric(y3),  method = "linear", xout )
    interp.subset.pop[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
    colnames(interp.subset.pop) = c("Country", "WHO_REGION",x)
  }
  
  return(list(interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop))
  
}




############################

##' Import data needed to generate animations

############################
get.data.for.animation <- function(regions){
  
  ##' import the case, birth rate, population and vaccination data for all countries
  cases.by.country.by.year = read.csv("data/Measles_cases_by_year.csv", stringsAsFactors = FALSE)
  Birth.rates = read.csv("data/Birth_rates.csv", stringsAsFactors = FALSE)
  pop.by.year = read.csv("data/All_populations.csv", stringsAsFactors = FALSE)
  vacc.rates = read.csv("data/Measles_vac_all.csv", stringsAsFactors = FALSE)
  
  
  ##' only include data for the specified WHO regions (generally use all regions)
  subset.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION %in% regions)
  subset.vaccination = subset(vacc.rates, vacc.rates$WHO_REGION %in% regions)
  subset.birth.rates = subset(Birth.rates, Birth.rates$WHO_REGION %in% regions)
  subset.data = subset(cases.by.country.by.year, cases.by.country.by.year$WHO_REGION %in% regions)
  
  ##' rename the Cname variable to Country
  colnames(subset.data)[3] = "Country"
  
  ##' edit the data, so that only countries which are in all data sets are included 
  missing1 = setdiff(subset.vaccination$Country, subset.pop.by.year$Country.Name)
  if(length(missing1) > 0){
    j = which(subset.vaccination$Country %in% missing1)
    subset.vaccination = subset.vaccination[-j, ]
  }
  missing2 = setdiff(subset.birth.rates$Country, subset.vaccination$Country)
  if(length(missing2) > 0){
    j = which(subset.birth.rates$Country %in% missing2)
    subset.birth.rates = subset.birth.rates[-j, ]
  }
  missing3 = setdiff(subset.data$Country, subset.vaccination$Country)
  if(length(missing3) > 0){
    j = which(subset.data$Country %in% missing3)
    subset.data = subset.data[-j, ]
  }
  
  ##' rename data sets for editing
  p1  =  subset.pop.by.year
  p2  =  subset.vaccination
  p3  =  subset.birth.rates
  p4  =  subset.data
  
  ##' edit data so that the entry for each country, is on the same line in each data set.
  for ( i in 1 : length(subset.vaccination[, 1])){
    C  =  subset.vaccination$Country[i]
    p2[i, ]  =  subset.vaccination[i, ]
    j = which(subset.pop.by.year$Country.Name == C)
    p1[i, ]  =  subset.pop.by.year[j, ]
    j = which(subset.birth.rates$Country == C)
    p3[i, ]  =  subset.birth.rates[j, ]
    j = which(subset.data$Country == C)
    p4[i, ]  =  subset.data[j, ]
  }
  
  ##" rename back to original names
  subset.pop.by.year = p1
  subset.vaccination = p2
  subset.birth.rates = p3
  subset.data = p4
  
  return(list(subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data))
}




############################

##' Get datasets into the same format

############################

interp.datasets <- function(subset.data, 
                            subset.vaccination, 
                            subset.birth.rates, 
                            subset.pop.by.year, 
                            x,
                            xout){
  
  interp.subset.data = matrix(0, length(subset.data[, 1]), length(xout) + 2)
  interp.subset.data[, 1] = subset.data$Country
  interp.subset.data[, 2] = subset.data$WHO_REGION
  
  interp.subset.vacc = matrix(0, length(subset.vaccination[, 1]), length(xout) + 2)
  interp.subset.vacc[, 1] = subset.data$Country
  interp.subset.vacc[, 2] = subset.data$WHO_REGION
  
  interp.subset.br = matrix(0, length(subset.birth.rates[, 1]), length(xout) + 2)
  interp.subset.br[, 1] = subset.data$Country
  interp.subset.br[, 2] = subset.data$WHO_REGION
  
  interp.subset.pop = matrix(0, length(subset.pop.by.year[, 1]), length(xout) + 2)
  interp.subset.pop[, 1] = subset.data$Country
  interp.subset.pop[, 2] = subset.data$WHO_REGION
  
  for ( i in 1 : length(subset.pop.by.year[, 1])){
    y = subset.data[i, paste("X", x, sep = "")]
    list[qq,ww] =  approx (as.numeric(x), as.numeric(y),  method = "linear", xout )
    interp.subset.data[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
    colnames(interp.subset.data) = c("Country", "WHO_REGION", x)
    
    y1 = subset.vaccination[i, paste("X", x, sep = "")]
    list[qq,ww] =  approx (as.numeric(x), as.numeric(y1),  method = "linear", xout )
    interp.subset.vacc[i, 3: length(interp.subset.vacc[1, ])]  =  round(ww, 2)
    colnames(interp.subset.vacc) = c("Country", "WHO_REGION", x)
    
    y2 = subset.birth.rates[i, paste("X", x, sep = "")]
    if(length(which(is.na(y2) == FALSE)) < 2) 
    {interp.subset.br[i, 3: length(interp.subset.data[1, ])]  = 0} else{
      list[qq,ww] =  approx (as.numeric(x), as.numeric(y2),  method = "linear", xout )
      interp.subset.br[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
    }
    colnames(interp.subset.br) = c("Country", "WHO_REGION", x)
    
    
    y3 = subset.pop.by.year[i, paste("X", x, sep = "")]
    list[qq,ww] =  approx (as.numeric(x), as.numeric(y3),  method = "linear", xout )
    interp.subset.pop[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
    colnames(interp.subset.pop) = c("Country", "WHO_REGION", x)
  }
  
  return(list(interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop))
  
}




############################

##' Continue preparing matrices

############################

prepare.matrices.for.animation <- function(interp.subset.data, subset.data){
  
  mean.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.cases[, 1] = subset.data$Country
  mean.cases[, 2] = subset.data$WHO_REGION
  
  coeff.var.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  coeff.var.cases[, 1] = subset.data$Country
  coeff.var.cases[, 2] = subset.data$WHO_REGION
  
  incidence.per.1000 = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  incidence.per.1000[, 1] = subset.data$Country
  incidence.per.1000[, 2] = subset.data$WHO_REGION
  
  mean.br = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.br[, 1] = subset.data$Country
  mean.br[, 2] = subset.data$WHO_REGION
  
  mean.vac = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
  mean.vac[, 1] = subset.data$Country
  mean.vac[, 2] = subset.data$WHO_REGION
  
  return(list(mean.cases, coeff.var.cases, incidence.per.1000, mean.br, mean.vac))
}





############################

##' output the gaussian weight for the averaging process

############################
output.weights.gaussian.with.cutoff <- function(x, st.dev, cutoff, neg.only = F){
  
  weights = dnorm(x, mean = 0, sd = st.dev)
  if(neg.only == T){
    weights[which(x > cutoff)] = 0
  } else{
    weights[which(abs(x) > cutoff)] = 0
  }
  return(weights)
}



############################

##' Interpolate datasets to the 

############################

#' Interpolate a given data set over a specified time frame.
#'
#' @param data - data to be interpolated
#' @param x - original time frame of the data
#' @param xout - time frame of the iinterpolated data
#'
#' @return - interpolated data set

interpolate.give.dataset <- function(data,  
                                     x,
                                     xout){
  ##' make a matrix to hold the interpolated data, which has number of rows equal to number
  ##' of countries in the data set, and number of colunms equal to the length of the chosen time frame
  interp.data = matrix(0, length(data[, 1]), length(xout))
  for ( i in 1 : length(data[, 1])){
    y = data[i, ]
    ##' if there are less than 2 non-NA entries in the data, then report 0, as there
    ##' will be no way to interpolate the data
    if(length(which(!is.na(y)) == F) < 2){
      interp.data[i, ]  =  0} else{
        ##' otherwise interpolate the data y, over the time frame xout, and capture the x
        ##' and interpolated y in qq and ww respectively
        list[qq,ww] =  approx (as.numeric(x), as.numeric(y),  method = "linear", xout )
        ##' put the interpolated data into the interpolated data matrix
        interp.data[i, ]  =  ww
      }
  }
  
  return(interp.data)
  
}




##' Estimate the number of infectious individuals by age for a given country in a given year.
##' Code and method via Takahashi et al. 2015
getLexisVaccPopStructSpecifyYear<- function(country="Sierra Leone",overlap=1,
                                            proportion.pop=FALSE, barchart.susc=FALSE,
                                            add.disrupted=FALSE, resurgence=FALSE,
                                            max.cover=0.95,
                                            cover.from.DHS=NULL,
                                            cover.from.DHS.disrupted=NULL, 
                                            year,
                                            max.x = 200000, output.plot = T){
  
  ages <- 1:60
  
  # 1. Routine Coverage - bring WHO data and find right country
  df <- read.csv("data/UNICEF.WHO.Measles.AdminCoverage.downloadedMay2014.csv")
  df <- read.csv("data/Measles.vac.WHO.estimates.csv")
  mtch <- match(df$Cname, country)
  
  coverage <- as.numeric(df[which(mtch==1 & !is.na(mtch),arr.ind=TRUE),paste("X",1980:min(year, 2013),sep="")]/100)
  if (country=="Liberia") coverage[2] <- 0.40 #Probable mistage Take the average for the decade
  #print(coverage)
  if (country=="Liberia") coverage[2] <- 0.40
  # -- Data ends in 2013 - extend from 2013 to 2015 - by assuming same level as 2013
  if(year > 2013){
    coverage <-c( coverage[1:length(coverage)],rep(coverage[length(coverage)],year - 2013))
  }
  
  if (is.na(coverage[1])) coverage[1] <- 0  #if 1980 value is NA, set it to zero
  #missing data fix
  for (j in 2:length(coverage)) {if (is.na(coverage[j])) coverage[j] <- coverage[j-1]} #set subsequent values to previous years value if NA
  #express as coverage in every cohort by reversing
  prop.vacc <- c(coverage[length(coverage):1],rep(0,60-length(coverage))) #get coverage in every cohort
  prop.vacc <- pmin(prop.vacc,max.cover)#constrain max to max.cover
  
  #adjust for 90% vaccine efficacy for all routine vaccination [note: Saki does not have this in]
  prop.vacc <- prop.vacc*0.95
  
  #for the 1-5 year olds, bring in data based on Saki's and multiply by vaccine efficacy
  #  prop.vacc[1:5] <- cover.from.DHS*0.95
  prop.vacc = c(prop.vacc,rep(tail(prop.vacc,1),60-length(prop.vacc)))
  #create disrupted - and bring in Saki's numbers - again multiply by vaccine efficacy
  prop.vacc.disrupted <- prop.vacc
  #prop.vacc.disrupted[1:5] <- cover.from.DHS.disrupted[1:5]*0.95
  
  # 2. SIAs - bring in WHO data, find right country, age range, year
  country1 = country
  #if( country == "Tanzania"){country1 = "United Republic of Tanzania"}
  # if(country == "Congo, Republic of the"){country1 = "Congo (the)"}
  #if(country == "Gambia, The"){country1 = "Gambia (the)"}
  #if(country == "Congo, Democratic Republic of the"){country1 = "Democratic Republic of the Congo (the)"}
  if(country == "Niger"){
    coverage.file <- "data/All_SIA_Routine_May2014_2.csv"
    year.now <- year
    df <- read.csv(coverage.file,stringsAsFactors =FALSE)
    mtch <- which(df$country=="Niger",arr.ind=TRUE)
    if (sum(is.na(mtch))==length(mtch)) print("could not match country name to data-file")  else df <- df[mtch,]
    years <- as.numeric(substring(df$date,5,nchar(df$date)))
    df <- df[years<year.now,]
  }else{
    coverage.file <- "data/All_SIA_Routine_May2014_2.csv"
    year.now <- year
    df <- read.csv(coverage.file,stringsAsFactors =FALSE)
    mtch <- which(df$country==country1,arr.ind=TRUE)
    if (sum(is.na(mtch))==length(mtch)) print("could not match country name to data-file")  else df <- df[mtch,]
    years <- as.numeric(substring(df$date,5,nchar(df$date)))
    df <- df[years<year.now,]
  }
  
  if(grepl("(the)", country, fixed = TRUE)) country =  gsub(" (the)","", country, fixed = T)
  if(grepl("(the)", country1, fixed = TRUE)){ country1 =  gsub(" (the)","", country1, fixed = T)}
  # 3. set up the SIAs - here also, we just need yearly timing
  # also assume 97% efficacy
  years.sia <- as.numeric(substring(df$date,5,nchar(df$date))[df$is.SIA==1])
  if (length(years.sia)>0) {
    coverage.sia <- pmin(df$percent.cov[df$is.SIA==1]*0.97,max.cover)    #adjust for 97% vaccine efficacy for all SIA delivery
    coverage.sia[is.na(coverage.sia)] <- 0 #if no coverage reported, assume 0%
    age.range <- cbind(df$age.low[df$is.SIA==1],df$age.high[df$is.SIA==1])
    
    for (j in 1:length(years.sia)) {
      if(!is.na(age.range[j,1]) & !is.na(age.range[j,2])){
        cov.ages.SIA <- year.now-years.sia[j]+1+c(age.range[j,1],age.range[j,2])/12
        find.bins.SIA <- findInterval(cov.ages.SIA,ages,all.inside=TRUE)
        
        # complete dependence
        p.overlap <- pmax(pmin(pmax(prop.vacc[find.bins.SIA[1]:find.bins.SIA[2]],coverage.sia[j]),1),0)
        # complete independence
        p.indep <- pmax(pmin((prop.vacc[find.bins.SIA[1]:find.bins.SIA[2]] +
                                (1-prop.vacc[find.bins.SIA[1]:find.bins.SIA[2]])*(coverage.sia[j])),1),0)
        # new cover
        prop.vacc[find.bins.SIA[1]:find.bins.SIA[2]] <- overlap * p.overlap + (1-overlap) * p.indep
        
        # complete dependence for disrupted
        p.overlap <- pmax(pmin(pmax(prop.vacc.disrupted[find.bins.SIA[1]:find.bins.SIA[2]],coverage.sia[j]),1),0)
        # complete independence  for disrupted
        p.indep <- pmax(pmin((prop.vacc.disrupted[find.bins.SIA[1]:find.bins.SIA[2]] +
                                (1-prop.vacc[find.bins.SIA[1]:find.bins.SIA[2]])*(coverage.sia[j])),1),0)
        # new cover  for disrupted
        prop.vacc.disrupted[find.bins.SIA[1]:find.bins.SIA[2]] <- overlap * p.overlap + (1-overlap) * p.indep
        
      }
      
    }
  }else{
    age.range = NA
    ages.vacc = NA
    coverage.sia = NA
  }
  
  # 4. Population structure
  
  
  ##' we look at any years between 2000 and 2015 to construct Lexis diagram. 
  ##' We can calculate the population at any multiple of 5 from 2000 to 2015, so to 
  ##' calculate population for in between years, we consider the population at either end of
  ##' these 5 year intervals and take the point at which we lie e.g. 2002 lies 2/5th's of 
  ##' the way between 2000 and 2005.
  
  if(year < 2015){
    possible.years = seq(1950, 2015, 5)
    y = which(possible.years > year)[1]
    z = possible.years[y]
    z1 = possible.years[y-1]
  } else{
    z = year
    z1 = 2010
  }
  
  
  p <- read.csv("data/Population_through_time.csv")
  country2 = country
  pop.adjust = 1
  if(country == "Andorra") {
    country2 = "Spain"
    pop.adjust = 77281/46570000}
  if(country == "Monaco"){
    country2 = "Japan"
    pop.adjust = 38695/126800000
  }
  if(country == "Tuvalu"){
    country2 = "Syria"
    pop.adjust = 11192/18270000
  }
  if(country == "Dominica"){
    country2 = "Albania"
    pop.adjust = 73925/2873000
  }
  if(country == "Saint Kitts and Nevis") {
    country2 = "Albania"
    pop.adjust = 55345/2873000
  }
  if(country == "Palau") {
    country2 = "Albania"
    pop.adjust = 21729/2873000
  }
  if(country == "Marshall Islands") {
    country2 = "Syria"
    pop.adjust = 53127/18270000
  }
  if(country == "San Marino") {
    country2 = "Denmark"
    pop.adjust = 33400/5770000
  }
  # if(country == "Korea, South") {country2 = "Republic of Korea"}
  # if(country == "Laos"){country2 = "Lao People's Democratic Republic"}
  # if(country == "Vietnam"){country2 = "Viet Nam"}
  # if(country == "Brunei"){country2 = "Brunei Darussalam"}
  # if(country == "Macedonia") {country2 = "TFYR Macedonia"}
  # if(country == "Venezuela"){country2 = "Venezuela (Bolivarian Republic of)"}
  # if(country == "United States") {country2 = "United States of America"}
  # if(country == "Micronesia, Federated States of"){country2 = "Micronesia (Fed. States of)"}
  # if(country == "Burma"){country2 = "Myanmar"}
  pop.struct.past = p[which(p$Region == country2 & p$Data == z1), 4:ncol(p)] * 1000 * pop.adjust
  if(z>2015){
    z=2015
  }
  pop.struct.future = p[which(p$Region == country2 & p$Data == z), 4:ncol(p)] * 1000 * pop.adjust
  
  fraction = year - z1
  if(fraction > 5){
    fraction = 5
  }
  pop.struct = (5 - fraction)/5 * pop.struct.past + fraction / 5 * pop.struct.future
  pop.struct[14] = sum(pop.struct[14:length(pop.struct)])
  pop.struct = pop.struct[1:14]
  age <- seq(5,70,by=5)
  sp <- smooth.spline(age,pop.struct)
  pred <- predict(sp,ages)$y
  #!make sure pop size stays the same
  pred <- sum(pop.struct)*pred/sum(pred)
  
  
  # 5. Natural Immunity 
  
  # This is the state space model adjusted disease burden.
  df <- read.csv("data/State_space_cases_new.csv",stringsAsFactors=FALSE,row.names=1)
  country.codes <- read.csv("data/country_codes.csv")
  iso <- country.codes[which(country.codes[,1]==country1),"ISO3_code"]
  if(country == "Cabo Verde"){iso = "CPV"}
    est.burden <- as.numeric(df[which(rownames(df)==iso), which(colnames(df) %in% paste("X", 1981:year,sep = ""))])
  base.burden <- mean(est.burden[5:15],na.rm=T)
  rel.burden <- est.burden / base.burden
  j = which(is.na(rel.burden))
  if(length(j) > 0){
    for(k in 1 : length(j)){
      if((j[k] > 1) & (j[k] < length(rel.burden))){
        p = which(!is.na(rel.burden[(j[k]+1): length(rel.burden)] ))[1]
        if(is.na(p) ){
          rel.burden[j[k]] = 0
        }else{
          rel.burden[j[k]] = mean(rel.burden[j[k]-1], rel.burden[j[k]+p])  
        }
      }
      if(j[k] == length(rel.burden)){
        rel.burden[j[k]] = rel.burden[j[k]-1]
      }
    }
  }
  # this needs to be considered -- what is the magnitude of measles resurgence ????
  # needs to be some kind of best/worst case scenario
  # currently set at the "baseline" which is the average of the 10 years centered around 1990
  rel.burden <- c(rep(1,60 - length(rel.burden[10:length(rel.burden)]) ),
                  rel.burden[10:length(rel.burden)])
  rel.burden[is.na(rel.burden)] = 1
  
  #     ifelse(resurgence==FALSE,rel.burden <- c(rep(1,60 - length(rel.burden[10:length(rel.burden)]) ),
  #                                              rel.burden[10:length(rel.burden)]),
  #            rel.burden <- c(rep(1,35),rel.burden[10:length(rel.burden)],rep(1,
  #                                                                            60 - 35 - length(rel.burden[10:length(rel.burden)]))))
  #     
  
  # ifelse(resurgence==FALSE,rel.burden <- c(rep(1,35),rel.burden[10:32],rep(rel.burden[32],2)),
  #         rel.burden <- c(rep(1,35),rel.burden[10:32],rep(1,2)))
  
  
  
  foi.base <- log(.05)/-20  				# baseline foi -- 95% immune by 20y -- should be tunable
  rel.foi <- rev(rel.burden)*foi.base
  prop.nat.imm <- 1-exp(-cumsum(rel.foi))		# cumulative probability of natural immunity
  #browser()
  
  
  # 6. Maternal immunity - remove those protected by maternal immunity from the susceptibles [approx 1/2 kids aged < 1]
  mat.protect <-c(pred[1]*0.5,rep(0,length(pred)-1))
  
  # 7. Figure
  if(output.plot ){
    
    layout(matrix(c(1,1,1,1,2,2),2,3))
    par(mar=c(4,4,3,2))
    
    plot(1980:year,type="n",ylim=c(0,40),xlim=c(1980,year - 1), xlab="", ylab="Age (years)",cex.lab = 1.6, cex.axis = 1.6)
    points(1980:year,0.5+0:(year - 1980),pch=15,col="skyblue", cex=coverage[1])  
    for (j in 0:(year - 1980)) { points((1980:year)+j,0.5+(0:(year - 1980)),pch=15,col="skyblue", cex=coverage[j+1]*1.5)}
    
    if (length(years.sia)>0) {
      for (j in 1:length(years.sia)) {
        ages.vacc <- seq(round(age.range[j,1]/12),round(age.range[j,2]/12),by=1)
        points(rep(years.sia[j],length(ages.vacc))+0.1,ages.vacc+
                 0.6,cex=coverage.sia[j]*1.5,col=4, pch=15)
      }}
    
    legend("topleft",legend=country,bty="n",cex=2.5)
    legend("topright",legend=c("routine","camapign"),col=c("skyblue","blue"),pch=15,bty="n", cex = 2)
    par(mar=c(4,0,3,2))
    
  }
  
  
  # plot out the lexis
  
  
  tot.imm <- mat.protect+(pred-mat.protect)*prop.vacc + (pred-mat.protect)*(1-prop.vacc)*(prop.nat.imm)
  tot.imm.disrupted <- mat.protect+ (pred-mat.protect)*prop.vacc.disrupted  +
    (pred-mat.protect)*(1-prop.vacc.disrupted)*(prop.nat.imm)
  #print(tot.imm.disrupted)
  
  if (!barchart.susc)  {
    if (proportion.pop) denom <- pred else denom <- 1
    if(output.plot == T){
      plot(pred/denom,ages,type="l",xlab="Population size", ylim=c(0,40),
           xlim=range(c(0,pred,pred*(1-prop.vacc)),na.rm=TRUE), cex.axis = 1.6, cex.lab= 1.6)
      points(pred*prop.vacc/denom,ages, type="l")
      points(tot.imm/denom,ages, type="l")
      
      polygon(c(pred/denom,rep(0,length(pred))),c(ages,ages[length(ages):1]),col="red",border="NA")
      polygon(c(tot.imm/denom,rep(0,length(pred))),c(ages,ages[length(ages):1]),col="thistle",border="NA")
      polygon(c(pred*prop.vacc/denom,rep(0,length(pred))),c(ages,ages[length(ages):1]),col="skyblue",border="NA")
      if (add.disrupted) polygon(c(pred*prop.vacc.disrupted/denom,rep(0,length(pred))),c(ages,ages[length(ages):1]),
                                 col="grey",border="NA")
      
      legend("topright",legend=c("vaccinated","naturally immune","unvaccinated"),col=c("skyblue","thistle","red"),
             pch=15,bty="n", cex = 2)
    }
  } else {
    if (proportion.pop) denom <- pred else denom <- rep(1,length(pred))
    tot.susc <- pred-tot.imm
    tot.susc.disrupted <- pred-tot.imm.disrupted
    
    xlims <- range(c(0,tot.susc/denom),na.rm=TRUE)
    if (add.disrupted) xlims <- range(c(0,tot.susc/denom,tot.susc.disrupted/denom),na.rm=TRUE)
    
    #quick and dirty for display
    xlims <- c(0,max.x)
    if(output.plot == T){
      if(proportion.pop == T){
        plot(tot.susc/denom,ages,type="n", xlab=" ", ylim=c(0,40))
      }
      else{plot(tot.susc/denom,ages,type="n", xlab=" ", ylim=c(0,40), xlim=xlims)}
      for (j in 1:length(tot.susc)) {
        if (add.disrupted) points(c(0,tot.susc.disrupted[j]/denom[j]),
                                  rep(ages[j],2)-0.5,lwd=6,type="l", col="coral", pch=15)
        points(c(0,tot.susc[j]/denom[j]),rep(ages[j],2)-0.5,lwd=6,type="l", col="seagreen4", pch=15)
        
      }
      if (add.disrupted & !proportion.pop)
        legend("topright",legend=c("Number of susceptibles","with 18 months disruption"),
               bty="n", pch=15,col=c("red","coral"))
      if (add.disrupted & proportion.pop)
        legend("topright",legend=c("Proportion of susceptibles","with 18 months disruption"),
               bty="n", pch=15,col=c("red","coral"))
    }
    
  }
  
  if (barchart.susc) return(list(tot.susc=tot.susc,tot.susc.disrupted=tot.susc.disrupted,
                                 popsize=pred,ages=ages, rel.foi=rel.foi,
                                 prop.vacc=prop.vacc,prop.vacc.disrupted=prop.vacc.disrupted,
                                 prop.nat.imm=prop.nat.imm, coverage, years.sia,
                                 age.range, ages.vacc, coverage.sia, mat.protect))
  
}



#########################
### Do multiple ggplots in a grid
#########################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}








#########################
### Prepare anim.data for analysis by adding columns for MCV2 use and Mean age of susceptibles.
### Also output African only data as we use this for most of the analysis.
#########################


prepare.anim.data.for.analysis <- function(input.data, mean.age.sus){
  require(dplyr)
  input.data = add.mcv2.data.to.anim.data (input.data, non.mcv2.alpha = 0.4, mcv2.alpha = 0.4)
  years = colnames(mean.age.sus)
  names = colnames(input.data)
  anim.data = cbind(input.data, matrix(0, nrow(input.data), 1))
  colnames(anim.data) = c(names, "Mean.Age.Sus")
  anim.data = anim.data %>% dplyr::filter(Year %in% years)
  Af.d = dplyr::filter(anim.data, WHO_REGION == "AFR", Country != "")
  Af.countries = unique(Af.d$Country)
  for( i in 1 : nrow(anim.data)){
    y = anim.data$Year[i]
    c = anim.data$Country[i]
    k = which(rownames(mean.age.sus)==c)
    if(length(k) > 0){
      anim.data$Mean.Age.Sus[i] = mean.age.sus[k, paste(y)]
    }
  }
  
  
  years = seq(min(anim.data$Year) + 1, max(anim.data$Year))
  require(dplyr)
  anim.data$MCV2 = as.numeric(as.character(anim.data$MCV2))
  anim.data$shape = as.numeric(as.character(anim.data$shape))
  
  Af.data = dplyr::filter(anim.data, WHO_REGION == "AFR", Year %in% years, !is.na(Mean.Age.Sus),
                   Country != "")
  Af.data = dplyr::filter(anim.data, WHO_REGION == "AFR", Year %in% years, !is.na(Mean.Age.Sus),
                   Country != "South Sudan")
  Af.data = Af.data[order(Af.data$Country), ]
  Af.data$Country = as.factor(Af.data$Country)
  Af.data = dplyr::filter(Af.data, Country != "")
  Af.data$br.vacc = (100-Af.data$Mean.vaccination) * Af.data$Mean.birth.rate / 100
  levels(Af.data$Country)[levels(Af.data$Country)=="Congo, Democratic Republic of the"] <- "DRC"
  levels(Af.data$Country)[levels(Af.data$Country)=="Congo, Republic of the"] <- "Republic of the Congo"
  levels(Af.data$Country)[levels(Af.data$Country)=="Gambia, The"] <- "The Gambia"
  Af.data = Af.data[order(Af.data$Country), ]
  
  
  Amr.data = dplyr::filter(anim.data, WHO_REGION == "AMR", Year %in% years, !is.na(Mean.Age.Sus),
                    Country != "" , Mean.birth.rate > 0 , Mean.vaccination > 0)
  Amr.data = Amr.data[order(Amr.data$Country), ]
  Amr.data$Country = as.factor(Amr.data$Country)
  Amr.data = dplyr::filter(Amr.data, Country != "")
  Amr.data$br.vacc = (100-Amr.data$Mean.vaccination) * Amr.data$Mean.birth.rate / 100
  Amr.data = Amr.data[order(Amr.data$Country), ]
  
  Rest.data = dplyr::filter(anim.data, WHO_REGION %in% c("WPR","EUR", "SEAR", "EMR"), Year %in% years, !is.na(Mean.Age.Sus),
                     Country != "" , Mean.birth.rate > 0 , Mean.vaccination > 0)
  Rest.data = Rest.data[order(Rest.data$Country), ]
  Rest.data$Country = as.factor(Rest.data$Country)
  Rest.data = dplyr::filter(Rest.data, Country != "")
  Rest.data$br.vacc = (100-Rest.data$Mean.vaccination) * Rest.data$Mean.birth.rate / 100
  Rest.data = Rest.data[order(Rest.data$Country), ]
  
  anim.data$br.vacc = (100-anim.data$Mean.vaccination) * anim.data$Mean.birth.rate / 100
  return(list(anim.data, Af.data, Amr.data, Rest.data))
}





first.sia.year.by.country.correct <- function (){
  ##' read in data
  SIA = read.csv("data/output-sia_age_corrected.csv", stringsAsFactors = F)
  SIA$Country = ""
  SIA$WHO_REGION = ""
  country.codes <- read.csv("data/country_codes.csv")
  k = which(country.codes$Report_country_name == "Serbia future")
  if(length(k) > 0){
    country.codes = country.codes[-(k),]
  }
  k = which(country.codes$Report_country_name == "Uruguay")
  country.codes$Region_Code[k] = "AMRO"
  for(i in 1 : nrow(SIA)){
    a = rownames(SIA)[i]
    j = which(country.codes$ISO3_code == a)
    if(length(j) > 0 ){
      SIA$Country[i] = as.character(country.codes$Report_country_name[j])
      pp = as.character(country.codes$Region_Code[j])
      SIA$WHO_REGION[i] = substring(pp,1, nchar(pp) -1)
    }
  }
  cs = (SIA$Country)
  reg = SIA$WHO_REGION
  ##' set up matrix
  first.sia = matrix(NA, length(cs), 3)
  first.sia [, 1:2] = cbind(cs, reg)
  ##' load packages for filter and str_sub
  require(dplyr)
  require(stringr)
  
  for(i in 1 : length(cs)){
    c = cs[i]
    a = SIA %>% filter(., Country == c)
    a = a[paste("X", 1980:2012, sep = "")]
    if(length(a) > 0){
      first.sia[i, 3] = which(a>0)[1] + 1979
    }
  }
  first.sia = data.frame(first.sia)
  colnames(first.sia) = c("Country", "Region", "Year")
  return(first.sia)
}


##########################################
##' work out when mean position in taxonomy space when sia introduced
##########################################
mean.sia.pos <- function(first.sia.by.country, countries = Af.countries,
                         data = Af.data, region){
  d = first.sia.by.country %>% filter(., Region == region)
  if(region == "AFR"){
    d$Country = as.character(d$Country)
    d$Country[which(d$Country == "C\xf4te d\x92Ivoire")] =
      "Cote d'Ivoire"
    d$Country[which(d$Country == "Democratic Republic of the Congo")] =
      "Congo, Democratic Republic of the"
    d$Country[which(d$Country == "Cape Verde")] =
      "Cabo Verde"
    d$Country[which(d$Country == "Gambia")] = "Gambia, The"
    d$Country[which(d$Country == "United Republic of Tanzania")] = "Tanzania"
    d$Country[which(d$Country == "Congo")] = "Congo, Republic of the"
  }
  
  sia = filter(d, Country %in% countries)
  
  first.sia.inc.cv = matrix(NA, nrow(sia), 4)
  for(i in 1 : nrow(sia)){
    c = as.character(sia$Country[i])
    y = as.numeric(as.character(sia$Year[i]))
    if(!is.na(y)){
      p = filter(data, Country == c, Year == y)
      first.sia.inc.cv[i, 1] = p$Coefficient.of.Variation
      first.sia.inc.cv[i, 2] = p$Incidence
      first.sia.inc.cv[i, 3] = p$Mean.birth.rate
      first.sia.inc.cv[i, 4] = p$Mean.vaccination
    }
    
  }
  
  mean.sia.intro = data.frame(cbind(mean(first.sia.inc.cv[, 1],na.rm = T), mean(first.sia.inc.cv[, 2], na.rm = T),
                                    mean(first.sia.inc.cv[, 3], na.rm = T),mean(first.sia.inc.cv[, 4], na.rm = T)))
  
  colnames(mean.sia.intro) = c("CV", "inc", "br", "vacc")
  return(mean.sia.intro)
}





##########################################
##' work out when median position in taxonomy space when sia introduced
##########################################
median.sia.pos <- function(first.sia.by.country, countries = Af.countries,
                           data = Af.data){
  sia = filter(first.sia.by.country, Country %in% countries)
  
  first.sia.inc.cv = matrix(NA, nrow(sia), 4)
  for(i in 1 : nrow(sia)){
    c = as.character(sia$Country[i])
    y = as.numeric(as.character(sia$Year[i]))
    if(!is.na(y)){
      p = filter(data, Country == c, Year == y)
      first.sia.inc.cv[i, 1] = p$Coefficient.of.Variation
      first.sia.inc.cv[i, 2] = p$Incidence
      first.sia.inc.cv[i, 3] = p$Mean.birth.rate
      first.sia.inc.cv[i, 4] = p$Mean.vaccination
    }
    
  }
  
  mean.sia.intro = data.frame(cbind(median(first.sia.inc.cv[, 1],na.rm = T), median(first.sia.inc.cv[, 2], na.rm = T),
                                    mean(first.sia.inc.cv[, 3], na.rm = T),mean(first.sia.inc.cv[, 4], na.rm = T)))
  
  colnames(mean.sia.intro) = c("CV", "inc", "br", "vacc")
  return(mean.sia.intro)
}



##########################################
##' work out when position in taxonomy space and br, vacc 
##########################################

#' Title
#'
#' @param MCV2.countries 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
calculate.mcv2.intro.pos <- function(MCV2.countries, data = Af.data){
  mcv2.intro.pos = matrix(0, nrow(MCV2.countries), 4)
  for(i in 1 : nrow(MCV2.countries)){
    if(as.numeric(MCV2.countries[i, 2]) < 1990){
      AAA = data[which(data$Country == as.character(MCV2.countries[i, 1]) &
                         data$Year == 1990), ]
      
    }else{
      AAA = data[which(data$Country == as.character(MCV2.countries[i, 1]) &
                         data$Year == as.numeric(MCV2.countries[i, 2])), ]
      
    }
    
    mcv2.intro.pos[i, ] = c(AAA$Coefficient.of.Variation, AAA$Incidence, 
                            AAA$Mean.vaccination, AAA$Mean.birth.rate)
  }
  
  mean.pos.br.vacc.mcv2 = cbind(mean(mcv2.intro.pos[,1 ]), mean(mcv2.intro.pos[, 2]),
                                mean(mcv2.intro.pos[, 3]), mean(mcv2.intro.pos[, 4]))
  mean.pos.br.vacc.mcv2 = data.frame(mean.pos.br.vacc.mcv2)
  colnames(mean.pos.br.vacc.mcv2) = c('CV', 'inc', 'vacc', 'br')
  return(mean.pos.br.vacc.mcv2)
}





##########################################
##' work out when position in taxonomy space and br, vacc - median
##########################################

#' Title
#'
#' @param MCV2.countries 
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
calculate.mcv2.intro.pos.median <- function(MCV2.countries, data = Af.data){
  mcv2.intro.pos = matrix(0, nrow(MCV2.countries), 4)
  for(i in 1 : nrow(MCV2.countries)){
    if(as.numeric(MCV2.countries[i, 2]) < 1990){
      AAA = data[which(data$Country == as.character(MCV2.countries[i, 1]) &
                         data$Year == 1990), ]
      
    }else{
      AAA = data[which(data$Country == as.character(MCV2.countries[i, 1]) &
                         data$Year == as.numeric(MCV2.countries[i, 2])), ]
      
    }
    
    mcv2.intro.pos[i, ] = c(AAA$Coefficient.of.Variation, AAA$Incidence, 
                            AAA$Mean.vaccination, AAA$Mean.birth.rate)
  }
  
  mean.pos.br.vacc.mcv2 = cbind(median(mcv2.intro.pos[,1 ]), median(mcv2.intro.pos[, 2]),
                                mean(mcv2.intro.pos[, 3]), mean(mcv2.intro.pos[, 4]))
  mean.pos.br.vacc.mcv2 = data.frame(mean.pos.br.vacc.mcv2)
  colnames(mean.pos.br.vacc.mcv2) = c('CV', 'inc', 'vacc', 'br')
  return(mean.pos.br.vacc.mcv2)
}




##########################################
##' calculate the mean position each year 
##########################################
arrow.position.by.year <- function (years, data){
  data = filter(data, Country != "")
  av.pos.by.year = matrix(0, length(years), 5)
  for(i in 1 : length(years)){
    d = data[which(data$Year == years[i]), ]
    av.pos.by.year[i, ] = cbind(years[i], mean(d$Coefficient.of.Variation), mean(d$Incidence),
                                mean(d$Mean.birth.rate), mean(d$Mean.vaccination))
  }
  
  av.pos.by.year = data.frame(av.pos.by.year)
  colnames(av.pos.by.year) = c('year', 'x', 'y','br','vacc')
  return(av.pos.by.year)   
}



##########################################
##' calculate the mean position each year 
##########################################
arrow.position.by.year.median <- function (years, data){
  data = filter(data, Country != "")
  av.pos.by.year = matrix(0, length(years), 5)
  for(i in 1 : length(years)){
    d = data[which(data$Year == years[i]), ]
    av.pos.by.year[i, ] = cbind(years[i], median(d$Coefficient.of.Variation), median(d$Incidence),
                                mean(d$Mean.birth.rate), mean(d$Mean.vaccination))
  }
  
  av.pos.by.year = data.frame(av.pos.by.year)
  colnames(av.pos.by.year) = c('year', 'x', 'y','br','vacc')
  return(av.pos.by.year)   
}



#############################################
##' Produce tufte style arrow figure
#############################################

#' Title
#'
#' @param arrow.data 
#' @param output.plot.points 
#' @param pos.mcv2 
#' @param pos.sia 
#' @param p.title 
#' @param normal.t 
#' @param mcv2.t 
#' @param sia.t 
#' @param p.font 
#' @param sia.col 
#' @param normal.col 
#' @param do.arrows 
#'
#' @return
#' @export
#'
#' @examples
plot.arrow.with.labels <- function(arrow.data, output.plot.points, pos.mcv2, pos.sia,
                                   p.title = "Mean path for countries in Africa, along with mean vaccination and birth rates
                                   over time", 
                                   normal.t="Mean position in taxonomy space when mean vaccination and birth rates first reach given values" ,
                                   mcv2.t = "Mean position upon introduction of MCV2",
                                   sia.t = "Mean position upon introduction of SIA",
                                   p.font = "sans",
                                   mcv.col = "rosybrown1",
                                   sia.col = "plum",
                                   normal.col = "skyblue1",
                                   do.arrows = F, x.axis.ticks = seq(0.75, 2, 0.25),
                                   y.axis.ticks = seq(0, 2.5, 0.5),
                                   arrow.size = 10,
                                   title = "Taxonomic regime change: Africa 1990-2014"){
  
  
  a <- ggplot() + 
    geom_path(data = arrow.data, aes(x, y, label = NULL ), size = arrow.size,
              arrow = arrow(type="closed"), color = "wheat2", alpha = 0.7) + 
    geom_point(data=output.years.to.plot, aes(x = x, y = y, size = 8), pch=1, col = normal.col, size = 8) +
    geom_text(data = output.years.to.plot, color = "gray30",family=p.font,
              mapping=aes(x=x + 0.045, y=y, label=paste(round(vacc, 1), ", ", round(br), sep = "")), size=4) +
    geom_point(aes( x = CV, y = inc), pch=1,data = pos.mcv2, colour = mcv.col, size = 8) +
    geom_text(data = pos.mcv2, color = "gray30",family=p.font,
              mapping=aes(x=CV + 0.045, y=inc, label=paste(round(vacc, 1), ", ", round(br), sep = "")), size=4) +
    geom_point(aes(x = CV, y = inc), pch = 1, data = pos.sia, colour = sia.col, size = 8) +
    geom_text(data = pos.sia, color = "gray30",family=p.font,
              mapping=aes(x=CV + 0.045, y=inc, label=paste(round(vacc, 1), ", ", round(br), sep = "")), size=4) 
  
  desc <- p.title %>%
    strwrap(width = 0.9 * getOption("width")) %>%
    paste0(collapse = "\n")
  
  normal.text <- normal.t %>%
    strwrap(width = 0.25 * getOption("width")) %>%
    paste0(collapse = "\n")
  
  
  mcv2.text <- mcv2.t %>%
    strwrap(width = 0.27 * getOption("width")) %>%
    paste0(collapse = "\n")
  
  sia.text <- sia.t %>%
    strwrap(width = 0.27 * getOption("width")) %>%
    paste0(collapse = "\n")
  
  
  if(do.arrows == T){
    qqq1<- a + 
      
      theme(axis.ticks = element_blank(),
            text=element_text(family=p.font),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.text = element_text(color = "gray30", size = 11,family=p.font),
            plot.title = element_text(face = "bold", hjust = 0.012,
                                      vjust = 0.8, color = "#3C3C3C", size = 20)) +
      
      labs(x = NULL, y = NULL, title = title) +
      
      # annotate("text", x = 0.73, y = 2.53, size = 4, p.fontface = "bold", color = "gray30",
      #          hjust = 0.38, vjust = -1, label = "Incidence") +
      
      geom_vline(xintercept = x.axis.ticks,
                 color = "wheat4", linetype = "dotted", size = 0.5)+
      
      geom_hline(yintercept = y.axis.ticks,
                 color = "wheat4", linetype = "dotted", size = 0.5)  +
      
      
      #annotate("text", x = 0.73, y = 2.53, size = 5, color = "gray30",
      #         hjust = 0.045, vjust = -0.51, label = desc,family=p.font) +
      
      geom_segment(aes(x = pos.mcv2$CV+.013, xend = pos.mcv2$CV + .12,
                       y = pos.mcv2$inc + .03, yend =  pos.mcv2$inc + 0.32),
                   color = mcv.col, 
                   alpha = 0.7, size = 0.7, show.legend = FALSE) +
      
      annotate("text", x =  pos.mcv2$CV + .13, y =  pos.mcv2$inc + .32, size = 4, color = "gray30",
               label = mcv2.text, hjust = 0, family=p.font) +
      
      geom_segment(aes(x = pos.sia$CV - .013, xend = pos.sia$CV - .12,
                       y = pos.sia$inc - .03, yend = pos.sia$inc - .32),color = sia.col, 
                   alpha = 0.7, size = 0.7, show.legend = FALSE) +
      
      annotate("text", x = pos.sia$CV - .18, y = pos.sia$inc - .40, size = 4, color = "gray30",
               label = sia.text, hjust = 0, family=p.font) +
      
      geom_segment(aes(x = arrow.data$x[1]+.013, xend = arrow.data$x[1] + .12,
                       y = arrow.data$y[1]- .03, yend = arrow.data$y[1] - .32),color = normal.col, 
                   alpha = 0.7, size = 0.7, show.legend = FALSE) +
      
      annotate("text", x = arrow.data$x[1] + .13, y =  arrow.data$y[1] - .32, size = 4, color = "gray30",
               label = normal.text, hjust = 0, family=p.font) 
    
  } else {
    qqq1<- a + 
      
      theme(axis.ticks = element_blank(),
            text=element_text(family=p.font),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.text = element_text(color = "gray30", size = 11,family=p.font),
            plot.title = element_text(face = "bold", hjust = 0.012,
                                      vjust = 0.8, color = "#3C3C3C", size = 20)) +
      
      labs(x = NULL, y = NULL, title = title) +
      
      # annotate("text", x = 0.73, y = 2.53, size = 4, fontface = "bold", color = "gray30",
      #          hjust = 0.38, vjust = -1, label = "Incidence") +
      
      geom_vline(xintercept = x.axis.ticks,
                 color = "wheat4", linetype = "dotted", size = 0.5)+
      
      geom_hline(yintercept = y.axis.ticks,
                 color = "wheat4", linetype = "dotted", size = 0.5)  
    
    
    
  }
  
  
  
  return(qqq1)
  
  
  
}








#' Make figure 1a from paper
#'
#' @param anim.data 
#' @param years 
#' @param regions 
#' @param shapes 
#' @param countries.of.interest 
#' @param colors 
#' @param text.size 
#' @param arrow.size 
#' @param xint 
#' @param yint 
#' @param line.color 
#' @param breaks 
#'
#' @return
#' @export
#'
#' @examples
incidence.space.fig <- function(anim.data = anim.data, years = c(1990, 2017), 
                                state.space = 0,
                                regions, shapes = c(1,2,16,17),
                                countries.of.interest = c("Malawi", "United States", "Brazil",
                                                          "Congo, Democratic Republic of the",
                                                          "Zambia", "Tanzania",
                                                          "Colombia"),
                                colors = c("deeppink", "royalblue1"), 
                                text.size = 4,
                                arrow.size = 3, xint = 1.7, yint = 7,
                                line.color = 'grey2',
                                breaks = c(1,5,10,20,40,80,120,160)){
  
  if(state.space == 1){
    years = c(min(anim.data$Year), max(anim.data$Year))
  }
  
  D = anim.data[which(anim.data$WHO_REGION %in% regions & 
                        anim.data$Year %in% years & 
                        anim.data$Country != ""), ]
  D$colour = matrix(0, nrow(D), 1)
  count = 1
  
  
  for(i in 1 : length(years)){
    for(j in 1 : length(regions)){
      k = which(D$Year == years[i] & D$WHO_REGION == regions[j])
      D$shape[k] = count
      D$colour[k] = count
      count = count + 1
    }
  }
  max.incidence = max(D$Incidence)
  D$alpha[which(D$Country %in% countries.of.interest)] = 1
  D = rbind(D, D[1, ])
  D$alpha[nrow(D)] = 0
  D$size = matrix(0, nrow(D), 1)
  D$size[which(D$Country %in% countries.of.interest)] = text.size
  if("Congo, Democratic Republic of the" %in% countries.of.interest){
    D$Country[which(D$Country == "Congo, Democratic Republic of the")] = "DRC"
  }
  if("United States" %in% countries.of.interest){
    D$Country[which(D$Country == "United States")] = "USA"
  }
  a <- ggplot(D, aes(x = Coefficient.of.Variation, y = 100*Incidence, label = Country)) 
  
  
  a <-  a + geom_point(aes( alpha = alpha, shape = factor(shape), color = factor(colour)), size = 10) +
    scale_shape_manual(values=shapes, name  ="",
                       breaks=c(1, 2, 3, 4),
                       labels=c("Africa 1990", "Americas 1990", "Africa 2017","Americas 2017"))+
    scale_colour_manual(values = rep(colors, 2),name  ="",
                        breaks=c(1, 2, 3, 4),
                        labels=c("Africa 1990", "Americas 1990", "Africa 2017","Americas 2017")) + 
    scale_y_continuous(limits = c(0, 100*round(max.incidence + 2)), trans= 'sqrt' ,
                       breaks = 100*breaks) + 
    theme_classic() + scale_alpha(guide = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=20), #, color = 'darkgrey'), 
          axis.title=element_text(size=20), #, color = 'darkgrey'),
          legend.position = "bottom", legend.text=element_text(size=16)) + geom_text(size = D$size) +
    labs(x = "coefficient of variation",
         y = paste("mean incidence per 100,000"), 
         color = "% not vaccinated" )
  
  
  y = seq(min(years), max(years))
  Mean_X = matrix(0, length(y), length(regions))
  Mean_Y = matrix(0, length(y), length(regions))
  for ( i in 1:length(y)){
    for(j in 1 : length(regions)){
      Mean_X[i, j] = mean(anim.data$Coefficient.of.Variation[which(anim.data$Year == y[i] &
                                                                     anim.data$WHO_REGION == regions[j] & 
                                                                     anim.data$Country != "")], na.rm = T)
      Mean_Y[i, j] = mean(anim.data$Incidence[which(anim.data$Year == y[i] &
                                                      anim.data$WHO_REGION == regions[j] & 
                                                      anim.data$Country != "")], na.rm = T)
    }
    
  }
  
  
  df <- data.frame(x = Mean_X, y  = 100*Mean_Y)
  qqq1<- a + geom_path(data = df, aes(x.1, y.1, label = NULL ), size = arrow.size,
                       arrow = arrow(), color = colors[1]) +
    geom_path(data = df, aes(x.2, y.2, label = NULL ), size = arrow.size,
              arrow = arrow(), color = colors[2])  +
    #geom_vline(xintercept = xint,  colour = line.color, linetype = 'dashed') +
    #geom_hline(yintercept = 100*yint,  colour = line.color, linetype = 'dashed') +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  return(qqq1)
}


incidence.space.fig.dashed.arrow <- function(anim.data = anim.data, years = c(1990, 2017), 
                                             state.space = 0,
                                             regions, shapes = c(1,2,16,17),
                                             countries.of.interest = c("Malawi", "United States", "Brazil",
                                                                       "Congo, Democratic Republic of the",
                                                                       "Zambia", "Tanzania",
                                                                       "Colombia"),
                                             colors = c("deeppink", "royalblue1"), 
                                             text.size = 4,
                                             arrow.size = 3, xint = 1.7, yint = 7,
                                             line.color = 'grey2',
                                             breaks = c(1,5,10,20,40,80,120,160)){
  
  if(state.space == 1){
    years = c(min(anim.data$Year), max(anim.data$Year))
  }
  
  years2 = years
  years2[2] = years[2] - 4
  D = anim.data[which(anim.data$WHO_REGION %in% regions & 
                        anim.data$Year %in% years2 & 
                        anim.data$Country != ""), ]
  D$colour = matrix(0, nrow(D), 1)
  count = 1
  
  
  for(i in 1 : length(years)){
    for(j in 1 : length(regions)){
      k = which(D$Year == years2[i] & D$WHO_REGION == regions[j])
      D$shape[k] = count
      D$colour[k] = count
      count = count + 1
    }
  }
  max.incidence = max(D$Incidence)
  D$alpha[which(D$Country %in% countries.of.interest)] = 1
  D = rbind(D, D[1, ])
  D$alpha[nrow(D)] = 0
  D$size = matrix(0, nrow(D), 1)
  D$size[which(D$Country %in% countries.of.interest)] = text.size
  if("Congo, Democratic Republic of the" %in% countries.of.interest){
    D$Country[which(D$Country == "Congo, Democratic Republic of the")] = "DRC"
  }
  if("United States" %in% countries.of.interest){
    D$Country[which(D$Country == "United States")] = "USA"
  }
  a <- ggplot(D, aes(x = Coefficient.of.Variation, y = 100*Incidence, label = Country)) 
  
  
  a <-  a + geom_point(aes( alpha = alpha, shape = factor(shape), color = factor(colour)), size = 10) +
    scale_shape_manual(values=shapes, name  ="",
                       breaks=c(1, 2, 3, 4),
                       labels=c("Africa 1990", "Americas 1990", "Africa 2014","Americas 2014"))+
    scale_colour_manual(values = rep(colors, 2),name  ="",
                        breaks=c(1, 2, 3, 4),
                        labels=c("Africa 1990", "Americas 1990", "Africa 2014","Americas 2014")) + 
    scale_y_continuous(limits = c(0, 100*round(max.incidence + 2)), trans= 'sqrt' ,
                       breaks = 100*breaks) + 
    theme_classic() + scale_alpha(guide = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=20), #, color = 'darkgrey'), 
          axis.title=element_text(size=20), #, color = 'darkgrey'),
          legend.position = "bottom", legend.text=element_text(size=16)) + geom_text(size = D$size) +
    labs(x = "coefficient of variation",
         y = paste("mean incidence per 100,000"), 
         color = "% not vaccinated" )
  
  
  y = seq(min(years), max(years))
  Mean_X = matrix(0, length(y), length(regions))
  Mean_Y = matrix(0, length(y), length(regions))
  for ( i in 1:length(y)){
    for(j in 1 : length(regions)){
      Mean_X[i, j] = mean(anim.data$Coefficient.of.Variation[which(anim.data$Year == y[i] &
                                                                     anim.data$WHO_REGION == regions[j] & 
                                                                     anim.data$Country != "")], na.rm = T)
      Mean_Y[i, j] = mean(anim.data$Incidence[which(anim.data$Year == y[i] &
                                                      anim.data$WHO_REGION == regions[j] & 
                                                      anim.data$Country != "")], na.rm = T)
    }
    
  }
  
  
  df <- data.frame(x = Mean_X, y  = 100*Mean_Y)
  qqq1<- a + geom_path(data = df, aes(x.1, y.1, label = NULL ), size = arrow.size,
                       arrow = arrow(), color = colors[1]) +
    geom_path(data = df[1:(nrow(df)-4),], aes(x.2, y.2, label = NULL ), size = arrow.size,
              arrow = arrow(), color = colors[2])  +
    geom_path(data = df[(nrow(df)-4):nrow(df),], aes(x.2, y.2, label = NULL ), size = arrow.size,
              arrow = arrow(), color = colors[2], linetype = 2)+
    #geom_vline(xintercept = xint,  colour = line.color, linetype = 'dashed') +
    #geom_hline(yintercept = 100*yint,  colour = line.color, linetype = 'dashed') +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  return(qqq1)
}




incidence.space.fig.median <- function(anim.data = anim.data, years = c(1990, 2017), 
                                       regions, shapes = c(1,2,16,17),
                                       countries.of.interest = c("Malawi", "United States", "Brazil",
                                                                 "Congo, Democratic Republic of the",
                                                                 "Zambia", "Tanzania",
                                                                 "Colombia"),
                                       colors = c("deeppink", "royalblue1"), 
                                       text.size = 4,
                                       arrow.size = 3, xint = 1.7, yint = 7,
                                       line.color = 'grey2',
                                       breaks = c(1,5,10,20,40,80,120,160)){
  
  D = anim.data[which(anim.data$WHO_REGION %in% regions & 
                        anim.data$Year %in% years & 
                        anim.data$Country != ""), ]
  D$colour = matrix(0, nrow(D), 1)
  count = 1
  for(i in 1 : length(years)){
    for(j in 1 : length(regions)){
      k = which(D$Year == years[i] & D$WHO_REGION == regions[j])
      D$shape[k] = count
      D$colour[k] = count
      count = count + 1
    }
  }
  max.incidence = max(D$Incidence)
  D$alpha[which(D$Country %in% countries.of.interest)] = 1
  D = rbind(D, D[1, ])
  D$alpha[nrow(D)] = 0
  D$size = matrix(0, nrow(D), 1)
  D$size[which(D$Country %in% countries.of.interest)] = text.size
  if("Congo, Democratic Republic of the" %in% countries.of.interest){
    D$Country[which(D$Country == "Congo, Democratic Republic of the")] = "DRC"
  }
  if("United States" %in% countries.of.interest){
    D$Country[which(D$Country == "United States")] = "USA"
  }
  a <- ggplot(D, aes(x = Coefficient.of.Variation, y = 100*Incidence, label = Country)) 
  
  
  a <-  a + geom_point(aes( alpha = alpha, shape = factor(shape), color = factor(colour)), size = 10) +
    scale_shape_manual(values=shapes, name  ="",
                       breaks=c(1, 2, 3, 4),
                       labels=c("Africa 1990", "Americas 1990", "Africa 2017","Americas 2017"))+
    scale_colour_manual(values = rep(colors, 2),name  ="",
                        breaks=c(1, 2, 3, 4),
                        labels=c("Africa 1990", "Americas 1990", "Africa 2017","Americas 2017")) + 
    scale_y_continuous(limits = c(0, 100*round(max.incidence + 2)), trans= 'sqrt' ,
                       breaks = 100*breaks) + 
    theme_classic() + scale_alpha(guide = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=20), #, color = 'darkgrey'), 
          axis.title=element_text(size=20), #, color = 'darkgrey'),
          legend.position = "bottom", legend.text=element_text(size=16)) + geom_text(size = D$size) +
    labs(x = "coefficient of variation",
         y = paste("median incidence per 100,000"), 
         color = "% not vaccinated" )
  
  
  y = seq(min(years), max(years))
  Median_X = matrix(0, length(y), length(regions))
  Median_Y = matrix(0, length(y), length(regions))
  for ( i in 1:length(y)){
    for(j in 1 : length(regions)){
      Median_X[i, j] = median(anim.data$Coefficient.of.Variation[which(anim.data$Year == y[i] &
                                                                         anim.data$WHO_REGION == regions[j] & 
                                                                         anim.data$Country != "")], na.rm = T)
      Median_Y[i, j] = median(anim.data$Incidence[which(anim.data$Year == y[i] &
                                                          anim.data$WHO_REGION == regions[j] & 
                                                          anim.data$Country != "")], na.rm = T)
    }
    
  }
  
  
  df <- data.frame(x = Median_X, y  = 100*Median_Y)
  qqq1<- a + geom_path(data = df, aes(x.1, y.1, label = NULL ), size = arrow.size,
                       arrow = arrow(), color = colors[1]) +
    geom_path(data = df, aes(x.2, y.2, label = NULL ), size = arrow.size,
              arrow = arrow(), color = colors[2])  +
    #geom_vline(xintercept = xint,  colour = line.color, linetype = 'dashed') +
    #geom_hline(yintercept = 100*yint,  colour = line.color, linetype = 'dashed') +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  return(qqq1)
}






incidence.space.fig.median.dashed.arrow <- function(anim.data = anim.data, years = c(1990, 2017), 
                                       regions, shapes = c(1,2,16,17),
                                       countries.of.interest = c("Malawi", "United States", "Brazil",
                                                                 "Congo, Democratic Republic of the",
                                                                 "Zambia", "Tanzania",
                                                                 "Colombia"),
                                       colors = c("deeppink", "royalblue1"), 
                                       text.size = 4,
                                       arrow.size = 3, xint = 1.7, yint = 7,
                                       line.color = 'grey2',
                                       breaks = c(1,5,10,20,40,80,120,160)){
  
  years2 = years
  years2[2] = years[2] - 4
  D = anim.data[which(anim.data$WHO_REGION %in% regions & 
                        anim.data$Year %in% years2 & 
                        anim.data$Country != ""), ]
  D$colour = matrix(0, nrow(D), 1)
  count = 1
  
  
  for(i in 1 : length(years)){
    for(j in 1 : length(regions)){
      k = which(D$Year == years2[i] & D$WHO_REGION == regions[j])
      D$shape[k] = count
      D$colour[k] = count
      count = count + 1
    }
  }
  max.incidence = max(D$Incidence)
  D$alpha[which(D$Country %in% countries.of.interest)] = 1
  D = rbind(D, D[1, ])
  D$alpha[nrow(D)] = 0
  D$size = matrix(0, nrow(D), 1)
  D$size[which(D$Country %in% countries.of.interest)] = text.size
  if("Congo, Democratic Republic of the" %in% countries.of.interest){
    D$Country[which(D$Country == "Congo, Democratic Republic of the")] = "DRC"
  }
  if("United States" %in% countries.of.interest){
    D$Country[which(D$Country == "United States")] = "USA"
  }
  a <- ggplot(D, aes(x = Coefficient.of.Variation, y = 100*Incidence, label = Country)) 
  
  
  a <-  a + geom_point(aes( alpha = alpha, shape = factor(shape), color = factor(colour)), size = 10) +
    scale_shape_manual(values=shapes, name  ="",
                       breaks=c(1, 2, 3, 4),
                       labels=c("Africa 1990", "Americas 1990", "Africa 2014","Americas 2014"))+
    scale_colour_manual(values = rep(colors, 2),name  ="",
                        breaks=c(1, 2, 3, 4),
                        labels=c("Africa 1990", "Americas 1990", "Africa 2014","Americas 2014")) + 
    scale_y_continuous(limits = c(0, 100*round(max.incidence + 2)), trans= 'sqrt' ,
                       breaks = 100*breaks) + 
    theme_classic() + scale_alpha(guide = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=20), #, color = 'darkgrey'), 
          axis.title=element_text(size=20), #, color = 'darkgrey'),
          legend.position = "bottom", legend.text=element_text(size=16)) + geom_text(size = D$size) +
    labs(x = "coefficient of variation",
         y = paste("median incidence per 100,000"), 
         color = "% not vaccinated" )
  
  
  y = seq(min(years), max(years))
  Median_X = matrix(0, length(y), length(regions))
  Median_Y = matrix(0, length(y), length(regions))
  for ( i in 1:length(y)){
    for(j in 1 : length(regions)){
      Median_X[i, j] = median(anim.data$Coefficient.of.Variation[which(anim.data$Year == y[i] &
                                                                         anim.data$WHO_REGION == regions[j] & 
                                                                         anim.data$Country != "")], na.rm = T)
      Median_Y[i, j] = median(anim.data$Incidence[which(anim.data$Year == y[i] &
                                                          anim.data$WHO_REGION == regions[j] & 
                                                          anim.data$Country != "")], na.rm = T)
    }
    
  }
  
  
  df <- data.frame(x = Median_X, y  = 100*Median_Y)
  qqq1<- a + geom_path(data = df, aes(x.1, y.1, label = NULL ), size = arrow.size,
                       arrow = arrow(), color = colors[1]) +
    geom_path(data = df[1:(nrow(df)-4),], aes(x.2, y.2, label = NULL ), size = arrow.size,
              arrow = arrow(), color = colors[2])  +
    geom_path(data = df[(nrow(df)-4):nrow(df),], aes(x.2, y.2, label = NULL ), size = arrow.size,
              arrow = arrow(), color = colors[2], linetype = 2)+
    #geom_vline(xintercept = xint,  colour = line.color, linetype = 'dashed') +
    #geom_hline(yintercept = 100*yint,  colour = line.color, linetype = 'dashed') +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  return(qqq1)
}









#' @examples
awesome.mega.figure.scaled <- function(anim.data = anim.data, years = c(1990, 2017), 
                                       regions, shapes = c(1,2,16,17),
                                       countries.of.interest = c("Malawi", "United States", "Brazil",
                                                                 "Congo, Democratic Republic of the",
                                                                 "Zambia",
                                                                 "Colombia"),
                                       colors = c("deeppink", "royalblue1"), 
                                       text.size = 4,
                                       arrow.size = 3, xint = 1.7, yint = 7,
                                       line.color = 'grey2',
                                       breaks = c(1,5,10,20,40,80,120,160)){
  
  D = anim.data[which(anim.data$WHO_REGION %in% regions & 
                        anim.data$Year %in% years & 
                        anim.data$Country != ""), ]
  D$colour = matrix(0, nrow(D), 1)
  count = 1
  for(i in 1 : length(years)){
    for(j in 1 : length(regions)){
      k = which(D$Year == years[i] & D$WHO_REGION == regions[j])
      D$shape[k] = count
      D$colour[k] = count
      count = count + 1
    }
  }
  max.incidence = max(D$Incidence)
  D$alpha[which(D$Country %in% countries.of.interest)] = 1
  D = rbind(D, D[1, ])
  D$alpha[nrow(D)] = 0
  D$size = matrix(0, nrow(D), 1)
  D$size[which(D$Country %in% countries.of.interest)] = text.size
  if("Congo, Democratic Republic of the" %in% countries.of.interest){
    D$Country[which(D$Country == "Congo, Democratic Republic of the")] = "DRC"
  }
  if("United States" %in% countries.of.interest){
    D$Country[which(D$Country == "United States")] = "USA"
  }
  a <- ggplot(D, aes(x = Coefficient.of.Variation, y = Incidence, label = Country)) 
  
  
  a <-  a + geom_point(aes( alpha = alpha, shape = factor(shape), color = factor(colour)), size = 10) +
    scale_shape_manual(values=shapes, name  ="",
                       breaks=c(1, 2, 3, 4),
                       labels=c("Africa 1990", "Americas 1990", "Africa 2017","Americas 2017"))+
    scale_colour_manual(values = rep(colors, 2),name  ="",
                        breaks=c(1, 2, 3, 4),
                        labels=c("Africa 1990", "Americas 1990", "Africa 2017","Americas 2017")) + 
    #scale_y_continuous(limits = c(0, round(max.incidence)), trans= 'sqrt' ,
    #                  breaks = breaks) + 
    theme_classic() + scale_alpha(guide = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=20), #, color = 'darkgrey'), 
          axis.title=element_text(size=20), #, color = 'darkgrey'),
          legend.position = "bottom", legend.text=element_text(size=16)) + geom_text(size = D$size) +
    labs(x = "scaled coefficient of variation",
         y = paste("scaled incidence"), 
         color = "% not vaccinated" )
  
  
  y = seq(min(years), max(years))
  Mean_X = matrix(0, length(y), length(regions))
  Mean_Y = matrix(0, length(y), length(regions))
  for ( i in 1:length(y)){
    for(j in 1 : length(regions)){
      Mean_X[i, j] = mean(anim.data$Coefficient.of.Variation[which(anim.data$Year == y[i] &
                                                                     anim.data$WHO_REGION == regions[j] & 
                                                                     anim.data$Country != "")], na.rm = T)
      Mean_Y[i, j] = mean(anim.data$Incidence[which(anim.data$Year == y[i] &
                                                      anim.data$WHO_REGION == regions[j] & 
                                                      anim.data$Country != "")], na.rm = T)
    }
    
  }
  
  
  df <- data.frame(x = Mean_X, y  = Mean_Y)
  qqq1<- a + geom_path(data = df, aes(x.1, y.1, label = NULL ), size = 2,
                       arrow = arrow(), color = colors[1]) +
    geom_path(data = df, aes(x.2, y.2, label = NULL ), size = 2,
              arrow = arrow(), color = colors[2])  +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  return(qqq1)
}



awesome.mega.figure.scaled.dashed.arrow <- function(anim.data = anim.data, years = c(1990, 2017), 
                                       regions, shapes = c(1,2,16,17),
                                       countries.of.interest = c("Malawi", "United States", "Brazil",
                                                                 "Congo, Democratic Republic of the",
                                                                 "Zambia",
                                                                 "Colombia"),
                                       colors = c("deeppink", "royalblue1"), 
                                       text.size = 4,
                                       arrow.size = 3, xint = 1.7, yint = 7,
                                       line.color = 'grey2',
                                       breaks = c(1,5,10,20,40,80,120,160)){
  
  years2 = years
  years2[2] = years[2] - 4
  D = anim.data[which(anim.data$WHO_REGION %in% regions & 
                        anim.data$Year %in% years2 & 
                        anim.data$Country != ""), ]
  D$colour = matrix(0, nrow(D), 1)
  count = 1
  
  
  for(i in 1 : length(years)){
    for(j in 1 : length(regions)){
      k = which(D$Year == years2[i] & D$WHO_REGION == regions[j])
      D$shape[k] = count
      D$colour[k] = count
      count = count + 1
    }
  }
  max.incidence = max(D$Incidence)
  D$alpha[which(D$Country %in% countries.of.interest)] = 1
  D = rbind(D, D[1, ])
  D$alpha[nrow(D)] = 0
  D$size = matrix(0, nrow(D), 1)
  D$size[which(D$Country %in% countries.of.interest)] = text.size
  if("Congo, Democratic Republic of the" %in% countries.of.interest){
    D$Country[which(D$Country == "Congo, Democratic Republic of the")] = "DRC"
  }
  if("United States" %in% countries.of.interest){
    D$Country[which(D$Country == "United States")] = "USA"
  }
  a <- ggplot(D, aes(x = Coefficient.of.Variation, y = Incidence, label = Country)) 
  
  
  a <-  a + geom_point(aes( alpha = alpha, shape = factor(shape), color = factor(colour)), size = 10) +
    scale_shape_manual(values=shapes, name  ="",
                       breaks=c(1, 2, 3, 4),
                       labels=c("Africa 1990", "Americas 1990", "Africa 2014","Americas 2014"))+
    scale_colour_manual(values = rep(colors, 2),name  ="",
                        breaks=c(1, 2, 3, 4),
                        labels=c("Africa 1990", "Americas 1990", "Africa 2014","Americas 2014")) + 
    #scale_y_continuous(limits = c(0, round(max.incidence)), trans= 'sqrt' ,
    #                  breaks = breaks) + 
    theme_classic() + scale_alpha(guide = "none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=20), #, color = 'darkgrey'), 
          axis.title=element_text(size=20), #, color = 'darkgrey'),
          legend.position = "bottom", legend.text=element_text(size=16)) + geom_text(size = D$size) +
    labs(x = "scaled coefficient of variation",
         y = paste("scaled incidence"), 
         color = "% not vaccinated" )
  
  
  y = seq(min(years), max(years))
  Mean_X = matrix(0, length(y), length(regions))
  Mean_Y = matrix(0, length(y), length(regions))
  for ( i in 1:length(y)){
    for(j in 1 : length(regions)){
      Mean_X[i, j] = mean(anim.data$Coefficient.of.Variation[which(anim.data$Year == y[i] &
                                                                     anim.data$WHO_REGION == regions[j] & 
                                                                     anim.data$Country != "")], na.rm = T)
      Mean_Y[i, j] = mean(anim.data$Incidence[which(anim.data$Year == y[i] &
                                                      anim.data$WHO_REGION == regions[j] & 
                                                      anim.data$Country != "")], na.rm = T)
    }
    
  }
  
  
  df <- data.frame(x = Mean_X, y  = Mean_Y)
  qqq1<- a + geom_path(data = df, aes(x.1, y.1, label = NULL ), size = 2,
                       arrow = arrow(), color = colors[1]) +
    geom_path(data = df[1:(nrow(df)-4),], aes(x.2, y.2, label = NULL ), size = arrow.size,
              arrow = arrow(), color = colors[2])  +
    geom_path(data = df[(nrow(df)-4):nrow(df),], aes(x.2, y.2, label = NULL ), size = arrow.size,
              arrow = arrow(), color = colors[2], linetype = 2)+
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  return(qqq1)
}




#' Get data ready for the arrow plot
#'
#' @param anim.data - data from which to get data on countries
#' @param regions - which regions do we want to calculate for
#'
#' @return - MCV2.countries: countries which have introduced MCV2 and the year which they did it,
#' MCV2.year: br and vacc of each country when they introduce MCV2, countries: which countries are in
#' this region, data - output from anim.data for region in country
#' @export
#'
#' @examples
sort.arrow.data <- function(d,
                            regions){
  
  years = seq(min(d$Year), max(d$Year))
  data = d[which(d$Year %in% years & d$WHO_REGION %in% regions &
                   d$Country != ""), ]
  countries= unique(data$Country)
  MCV2.year = data.frame(matrix(NA, length(countries), 4))
  colnames(MCV2.year) = c("Country", "Intro.year", "BR", "Vacc")
  MCV2.year[, 1] = countries
  for(i in 1 : length(countries)){
    g = data[which(data$Country == countries[i]), ]
    if(length(which(g$MCV2 == 1)) > 0){
      MCV2.year[i, 2] = g$Year[which(g$MCV2 == 1)][1]
      MCV2.year[i, 3:4] = c(g$Mean.birth.rate[which(g$MCV2 == 1)][1], 
                            g$Mean.vaccination[which(g$MCV2 == 1)][1])
    }
  }
  
  
  MCV2.countries = MCV2.year[which(!is.na(MCV2.year[,2])), ]
  MCV2.countries = data.frame(MCV2.countries)
  MCV2.countries[, 2] = as.numeric(as.character(MCV2.countries[, 2]))
  MCV2.countries[, 3] = MCV2.countries[,2] - min(data$Year)
  MCV2.countries[, 4] = max(data$Year) - MCV2.countries[,2]
  
  colnames(MCV2.countries) = c("Country", "Intro_Year", "Years_before_Intro", "Years_after_Intro")
  return(list(MCV2.countries, MCV2.year, countries, data))
}







#' Output animation data when we use the state space corrected cases to generate
#' position over time.
#'
#' @param window.length 
#' @param regions 
#' @param gaussian.st.dev 
#' @param cutoff 
#' @param interp.resolution 
#' @param year.shift.inc 
#'
#' @return
#' @export
#'
#' @examples
back.only.method.state.space <- function(window.length, regions, 
                                         gaussian.st.dev, cutoff = 50, 
                                         interp.resolution = 20,
                                         year.shift.inc = 3){
  
  require(stats)
  list[subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data] = 
    get.data.for.animation.state.space(regions)
  
  
  x = seq(1981, 2012)
  colnames(subset.data) = c(seq(1981, 2012),"Country", "WHO_REGION")
  ##' interpolate the datasets to have entries for all points in time once the interpolation is done.
  list[interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop] = 
    interp.datasets.state.space(subset.data, 
                                subset.vaccination, 
                                subset.birth.rates, 
                                subset.pop.by.year,
                                x,
                                x)
  
  # list[mean.cases, coeff.var.cases, incidence.per.1000, mean.br, mean.vac] = 
  #     prepare.matrices.for.animation(interp.subset.data, subset.data)
  
  ##' number of unique years that we will have data for. The longer the window, the less unique years of data.
  num.windows = length(x) - window.length + 1
  
  ##' first year of data
  year = x[1]
  
  ##' setting up the datasets
  coeff.var = matrix(0, length(subset.data[ , 1]), num.windows)
  incidence.per.1000 = matrix(0, length(subset.data[ , 1]), num.windows)
  mean.br = matrix(0, length(subset.data[ , 1]), num.windows)
  mean.vac = matrix(0, length(subset.data[ , 1]), num.windows)
  
  ##' do calculations that calculate the coefficient of variation, incidence per 100, mean birth rate and
  ##' mean vaccination rate over periods of length given by the window length.
  for ( j in 1 : num.windows){
    for ( i in 1 : length(subset.data[ , 1])){
      coeff.var[i, j]  =  sd(interp.subset.data[i, paste(seq(year, year + window.length - 1))], na.rm = TRUE) /  
        mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      if(is.na( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)) == FALSE ){
        if( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) == 0) {
          coeff.var[i, j]  =  0
        } 
      }
      incidence.per.1000[i, j]  =  sum(as.numeric(interp.subset.data[i,  paste(seq(year, year + window.length - 1))]) / 
                                         as.numeric(interp.subset.pop[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) * 1000
      if(length(which(is.na(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))])))) < window.length){
        mean.br[i, j]  =  mean(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      }
      if(length(which(is.na(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))])))) < window.length){
        mean.vac[i, j]  =  mean(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      }
    }
    year = year + 1
  }
  
  incidence.per.1000.each.year = matrix(0, nrow(subset.data), length(x))
  for(i in 1 : nrow(subset.data)){
    incidence.per.1000.each.year[i, ] = 1000 * as.numeric(interp.subset.data[i,  paste(x)]) / 
      as.numeric(interp.subset.pop[i, paste(x)])
    
  }
  
  ##' we take the weighted average of the values that we have calculated for each year, where the weights are gaussian,
  ##' with a specified number of years for the standard deviation
  x1 = seq(1, length(coeff.var[1, ]) + 1)
  w1 = matrix(0, length(x1), length(x1))
  for (i in 1 : length(x1)){
    ##' set up the gaussian weights for averaging
    
    w.input = x1 - x1[i]
    w1[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
    
    ##' make sure that the weights add up to 1 for each of the specific weightings
    j = which(w.input == 0)
    w1[i,j : length(w1[i,])] = 0
    w1[i, ] = w1[i, ] / sum(w1[i, ])
  }
  
  w1 = w1[-(1), ]
  w1 = w1[, -(ncol(w1)) ]
  
  x2 = seq(1, length(incidence.per.1000.each.year[1, ]))
  w2 = matrix(0, length(x2), length(x2))
  for (i in 1 : length(x2)){
    ##' set up the gaussian weights for averaging
    
    w.input = x2 - x2[i] + year.shift.inc
    w2[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
    
    ##' make sure that the weights add up to 1 for each of the specific weightings
    j = which(w.input == (year.shift.inc + 1))
    if(length(j) > 0){
      w2[i,j : length(w2[i,])] = 0
      w2[i, ] = w2[i, ] / sum(w2[i, ])
    }
    
  }
  #w2 = w2[-(1), ]
  #w2 = w2[, -(ncol(w2)) ]
  w2 = w2[-(1:9), ]
  
  ##' make a set of matrices that are the same size as the matrices containing the data.
  coeff.2 = coeff.var
  incidence.2 = incidence.per.1000
  mbr2 = mean.br
  mvacc2 = mean.vac
  for(i in 1 : length(coeff.var[1, ])){
    for(j in 1 : length(coeff.var[, 1])){
      ##' make the entries of these newly created matrices to be the weighted averages of the originally calculated datasets
      coeff.2[j, i] = sum(coeff.var[j, ] * w1[i, ], na.rm = T)
      incidence.2[j, i] = sum(incidence.per.1000.each.year[j, ] * w2[i, ], na.rm = T)
      mbr2[j, i] =  sum(mean.br[j, i] * w2[i, ], na.rm = T)
      mvacc2[j, i] = sum(mean.vac[j, i]* w2[i, ], na.rm = T)
      
      ##' Should we do weighted average of birth rate and vaccination rate?
      ##' If so uncomment the next two lines
      
      # mbr2[j, i] = sum(mean.br[j, ] * w1[i, ], na.rm = T)
      # mvacc2[j, i] = sum(mean.vac[j, i] * w1[i, ], na.rm = T)
    }
  }
  
  ##' set the original data to be equal to the weighted data
  coeff.var.cases = coeff.2
  incidence.per.1000 = incidence.2
  mean.br = mbr2
  mean.vac = mvacc2
  
  
  ##' set up the timeline on which we do the interpolation. 
  ##' The number of sections that the yearly data is split up to is given by interp.resolution  
  x = seq(1981 + (window.length - 1), 2012)
  xout = seq(1981+ (window.length - 1), 2012, 1/interp.resolution)
  
  ##' interpolate the data and add columns that contain the correspondin country and WHO region of each line
  
  coeff.var.cases = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(coeff.var.cases, x, xout))
  incidence.per.1000 = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(incidence.per.1000, x, xout))
  mean.br = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.br, x, xout))
  mean.vac = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.vac, x, xout))
  
  ##' round the data to 2 decimal places for ease of reading.
  mean.vac[, -(1:2)] = as.numeric(mean.vac[, -(1:2)])
  mean.br[, -(1:2)] = as.numeric(mean.br[, -(1:2)])
  incidence.per.1000[, -(1:2)] = as.numeric(incidence.per.1000[, -(1:2)])
  coeff.var.cases[, -(1:2)] = as.numeric(coeff.var.cases[, -(1:2)])
  
  
  ##' set up the output to be the appropriate size and add column labels.
  ##' Additionally add enough room to include additional data for each year that will be used 
  ##' to calibrate the data for each year, so that the minimum and maximum of vaccination rate is 0 and 100 each time.
  ##' This ensures that the colour scale is constant
  
  output.data = matrix(0, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + (2 * length(regions) * length(xout)), 7)
  output.data  =  data.frame(output.data)
  colnames(output.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
  
  ##' input the appropriate data to the outputs 
  output.data$Country[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 1], length(coeff.var.cases[1, -(1:2)]))
  output.data$WHO_REGION[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 2], length(coeff.var.cases[1, -(1:2)]))
  count = 1
  for(i in 3 : length(coeff.var.cases[1, ])){
    output.data$Coefficient.of.Variation[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = coeff.var.cases[, i]
    output.data$Incidence[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = incidence.per.1000[, i]
    output.data$Mean.vaccination[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.vac[, i]
    output.data$Mean.birth.rate[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.br[, i]
    output.data$Year[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = xout[i - 2]  
    count = count + 1
  }
  
  ##' add in the dummy data for each year to keep the scales constant.
  year.mins = matrix(0, length(xout), 2)
  
  for(i in 1 : length(xout)){
    t  =  subset(output.data, output.data$Year ==  unique(xout)[i])
    year.mins[i, 1] = xout[i]
    year.mins[i, 2] = as.numeric(min(t$Mean.birth.rate,na.rm = T)  )
  }
  
  l =  expand.grid("", -1, 0, c(0,100), 0, xout, regions)
  
  for(i in 1 : (2 * length(regions) * length(xout))){
    y = l[i, 6]
    j = which(year.mins[, 1] == y)
    l[i, 5]  =  year.mins[j, 2]
  }
  
  output.data[((length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + 1 ): length(output.data[, 1]), ]  =  l
  
  for( i in 1 : (2 * length(regions) * length(xout))){
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] = regions[as.numeric(output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] )]
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)])  + i, 1] = ""
  }
  
  ##' make sure that each column that should be numeric is numeric.
  output.data$Coefficient.of.Variation = as.numeric(output.data$Coefficient.of.Variation)
  output.data$Incidence  =  as.numeric(output.data$Incidence)
  output.data$Mean.vaccination  =  as.numeric(output.data$Mean.vaccination)
  output.data$Mean.birth.rate  =  as.numeric(output.data$Mean.birth.rate)
  output.data$Year   =  as.numeric(output.data$Year)
  output.data$Coefficient.of.Variation[which(output.data$Coefficient.of.Variation == "Inf")] = 0
  return(output.data)
}




#' Function to calculate the mean path taken by the regions through taxonomy space
#'
#' @param d - data set to calculate from
#' @param regions - regions to perform analysis for
#'
#' @return mean.x, mean.y - the x and y values of the mean position over time.
#' @export
#'
#' @examples

find.mean.path <- function(d, regions){
  
  years = unique(d$Year)
  mean.x = data.frame(matrix(0, length(years), length(regions) + 2))
  mean.y = data.frame(matrix(0, length(years), length(regions) + 2))
  mean.x[, 1] = years
  mean.y[, 1] = years
  
  for(j in 1 : length(regions)){
    mean.x[, j + 1] = d %>%
      filter(WHO_REGION == regions[j]) %>%
      group_by(Year) %>%
      dplyr::summarise(pos = mean(Coefficient.of.Variation)) %>%
      as.data.frame %>%
      dplyr::select(pos) 
    
    
    mean.y[, j + 1] = d %>%
      filter(WHO_REGION == regions[j]) %>%
      group_by(Year) %>%
      dplyr::summarise(pos = mean(Incidence)) %>%
      as.data.frame %>%
      dplyr::select(pos)
  }
  mean.x[, j + 2] = as.matrix(d %>%
                                group_by(Year) %>%
                                dplyr::summarise(pos = mean(Coefficient.of.Variation)) %>%
                                dplyr::select(pos))
  
  mean.y[, j + 2] = as.matrix(d %>%
                                group_by(Year) %>%
                                dplyr::summarise(pos = mean(Incidence)) %>%
                                dplyr::select(pos))
  colnames(mean.x) = c("year", regions, 'all')
  colnames(mean.y) = c("year", regions, 'all')
  return(list(mean.x, mean.y))
}







#' Title
#'
#' @param d 
#' @param regions 
#' @param years 
#' @param use.rep.cases 
#' @param make.inc.cv.scale.same 
#' @param connect.canonical.path.to.zero 
#' @param log.incidence 
#'
#' @return
#' @export
#'
#' @examples

closest.path.point <- function(d, regions, years, use.rep.cases = 1,
                               make.inc.cv.scale.same = TRUE, sqrt.inc = F, sqrt.cv=F,
                               connect.canonical.path.to.zero = 0, log.incidence = F,
                               make.figure.plot = F){
  
  ##' filter the data to only include the regions, and years required, along with only
  ##' entries which have incidence greater than 0.
  d = d %>% filter(., Year %in% years, WHO_REGION %in% regions, Incidence > 0,
                   Country != "")
  
  ##' remove years where the coefficient of variation is 0, and the year is before 1995,
  ##' as before this time, no countries had reached elimination, so the only way to 
  ##' have a 0 here, is if there is no data before this time (unless a country reported
  ##' the exact same non-zero number of cases for every year, which didn't happen)
  j = which(d$Year < 1995 & d$Coefficient.of.Variation == 0)
  if(length(j) > 0){
    d = d[-(j), ]
  }
  
  ##' if make.figure.plot = T then we will produce the plot that is seen in the paper, which is 
  ##' make by putting the incidence data onto the same range as the coefficient of variation data.
  if(make.figure.plot == T){
    d$Incidence = d$Incidence* max(d$Coefficient.of.Variation)/ max(d$Incidence)
  }
  
  ##' if we want to log the incidence before assigning countries to their location
  ##' on the path, then we do that here. Along with adding a small non-zero value to each of the
  ##' incidences, so that there are no -Inf values created.
  if(log.incidence == T){
    d$Incidence = log(d$Incidence + 0.000001)
  }
  
  ##' if sqrt.inc=T then take sqrt of incidence
  if(sqrt.inc == T){
    d$Incidence <-  sqrt(d$Incidence)
  }
  
  ##' if sqrt.cv=T then take sqrt of cv
  if(sqrt.cv == T){
    d$Coefficient.of.Variation = sqrt(d$Coefficient.of.Variation)
  }
  
  ##' if wanted scale the incidence and coefficient of variation so that they are both in the 0-1 range
  if(make.inc.cv.scale.same == T){
    d$Incidence = scl(d$Incidence)
    d$Coefficient.of.Variation = scl(d$Coefficient.of.Variation)
  }
  
  ##' find the mean path for each region by averaging positions in each year
  ##' 
  list[mean.x, mean.y] = find.mean.path(d = filter(d, Year <= 2013), regions)
  
  
  if(use.rep.cases ==  1){
    ##' If we are considering the reported cases to locate the countries, then
    ##' Africa and the Americas mean paths converge at roughly 2007 in Africa and 
    ##' 1993 for the Americas. Join the paths at this point, and use it as the canonical
    ##' path 
    k = which(years == 2008)
    k1 = which(years == 1995)
    #zzz = mean.x$AFR[k] + (c(1,2,3)* (mean.x$AMR[k1] - mean.x$AFR[k]))/4
    #zzz1 = mean.y$AFR[k] + (c(1,2,3)* (mean.y$AMR[k1] - mean.y$AFR[k]))/4
    canonical.path = rbind(cbind(mean.x$AFR[1:k], mean.y$AFR[1:k]),# cbind(zzz,zzz1),
                           cbind(mean.x$AMR[k1:nrow(mean.x)], mean.y$AMR[k1:nrow(mean.x)])) %>%
      data.frame()
    
    
  }else{
    ##' If we are considering the state space cases to locate the countries, then
    ##' Africa and the Americas mean paths converge at roughly 2011 in Africa and 
    ##' 1994 for the Americas. We add in 3 intermediate points here and join the paths at this point, 
    ##' and use it as the canonical path
    k = which(years == 2008)
    k1 = which(years == 1995)
    zzz = mean.x$AFR[k] + (c(1,2,3)* (mean.x$AMR[k1] - mean.x$AFR[k]))/4
    zzz1 = mean.y$AFR[k] + (c(1,2,3)* (mean.y$AMR[k1] - mean.y$AFR[k]))/4
    canonical.path = rbind(cbind(mean.x$AFR[1:k], mean.y$AFR[1:k]), cbind(zzz,zzz1),
                           cbind(mean.x$AMR[k1:nrow(mean.x)], mean.y$AMR[k1:nrow(mean.x)])) %>%
      data.frame()
    canonical.path = canonical.path[-c(8,29,31),]
    
  }
  colnames(canonical.path) = c('x', 'y')
  last.canonical = c(tail(canonical.path$x, 1), tail(canonical.path$y, 1))
  
  end.canonical =c(0,0)
  
  if(connect.canonical.path.to.zero == 1){
    l = cbind(seq(last.canonical[1], end.canonical[1], -last.canonical[1] / 3),
              seq(last.canonical[2], end.canonical[2], -last.canonical[2] / 3))
    
    colnames(l) = colnames(canonical.path)
    canonical.path = rbind(canonical.path, l)
  }
  
  # canonical.path$x = scl(canonical.path$x)
  # canonical.path$y = scl(canonical.path$y)
  ##' add a column to hold the closest position on the mean trajectory
  d$closest = 0
  
  ##' we only need to calculate the closest position for each country once a year, 
  ##' therefore we restrict to only yearly data
  y1 = unique(round(d$Year))
  d = filter(d, Year %in% y1)
  
  ##' for each entry, then calculate which point is the closest on the trajectory
  for(i in 1 : nrow(d)){
    q = c(d$Coefficient.of.Variation[i], d$Incidence[i]) %>%
      rbind( cbind(canonical.path$x, canonical.path$y)) %>%
      dist() %>% as.matrix()
    d$closest[i] = which(q[1, -(1)] == min(q[1, -(1)])) %>% as.matrix() %>% tail(1)
    
  }
  d$canonical.x = canonical.path$x[d$closest]
  d$canonical.y = canonical.path$y[d$closest]
  return(list(d, canonical.path))
}











#' Generate the weights for the averaging of the Coefficient of Variation
#'
#' @param coeff.var - data set which is used for indicating the number of years in the analysis
#' @param gaussian.st.dev - the standard deviation of the gaussian weigthing function
#' @param cutoff - the number of years in the past at which we stop including data. Generally 
#' this is set at 50, meaning that all past data is used.
#'
#' @return - w1, a set of weights
generate.cv.weights <- function(coeff.var, gaussian.st.dev, cutoff){
  x1 = seq(1, length(coeff.var[1, ]) + 1)
  w1 = matrix(0, length(x1), length(x1))
  for (i in 1 : length(x1)){
    ##' set up the gaussian weights for averaging
    
    w.input = x1 - x1[i]
    w1[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
    
    ##' find the position in the time input at which we are at, and make all weights 0
    ##' for any time after that, so that we are only using past values for the weigthing
    j = which(w.input == 0)
    w1[i,j : length(w1[i,])] = 0
    ##' make sure that the weights add up to 1 for each of the specific weightings
    w1[i, ] = w1[i, ] / sum(w1[i, ])
  }
  ##' remove the first row, as this is all NA
  w1 = w1[-(1), ]
  ##' remove the final column, as this is always all 0, due to only including previoulsy observed values.
  w1 = w1[, -(ncol(w1)) ]    
  return(w1)
}



#' Generate the weights for the averaging of the incidence, birth rate and vaccination rate
#'
#' @param d - data set which is used for indicating the number of years in the analysis
#' @param gaussian.st.dev - the standard deviation of the gaussian weighting function. Set at 3 for analysis in the paper.
#' @param cutoff - the number of years in the past at which we stop including data. For analysis in the paper 
#' this is set at 50, meaning that all past data is used.
#' @param year.shift.inc - the number of years in the past that we centre the Gaussian
#' distribution at, rather than including the year in question as the centre of the Gaussian.
#' Was set at 2 for all analysis in the paper

#' @return - w2, a set of weights
generate.other.weights <- function(d, 
                                   window.length,
                                   gaussian.st.dev, 
                                   cutoff,
                                   year.shift.inc){
  x2 = seq(1, length(d[1, ]))
  w2 = matrix(0, length(x2), length(x2))
  for (i in 1 : length(x2)){
    ##' set up the gaussian weights for averaging
    
    w.input = x2 - x2[i] + year.shift.inc
    w2[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
    
    ##' make sure that the weights add up to 1 for each of the specific weightings
    j = which(w.input == (year.shift.inc + 1))
    if(length(j) > 0){
      w2[i,j : length(w2[i,])] = 0
      w2[i, ] = w2[i, ] / sum(w2[i, ])
    }
    
  }
  ##' we only begin weighting the data sets at the point at which we have at least the number of years
  ##' of data as we include each windowed data set. 
  ##' therefore, we delete the entries before this number of years 
  w2 = w2[-(1:(window.length - 1)), ]
  
  return(w2)
}


#' 
#' 
#' #' Create weighted data sets, from original non-weighted data, and the matrices of weights
#' #'
#' #' @param coeff.var - data set with cv data
#' #' @param incidence.per.1000.each.year - data set with incidence data
#' #' @param mean.br - data set with birth data
#' #' @param mean.vacc - data set with vaccination data
#' #' @param w1 - matrix of weights for weighting cv data
#' #' @param w2 - matrix of weights for other data sets
#' #'
#' #' @return the weighted data sets are returned as a list
#' #' list(coeff.var.cases, incidence.per.1000, mean.br, mean.vac)
#' make.weighted.data <- function(coeff.var, 
#'                                incidence.per.1000.each.year,
#'                                mean.br, 
#'                                mean.vac,
#'                                w1,
#'                                w2){
#'     
#'     coeff.2 = coeff.var
#'     incidence.2 = incidence.per.1000.each.year
#'     mbr2 = mean.br
#'     mvacc2 = mean.vac
#'     for(i in 1 : length(coeff.var[1, ])){
#'         for(j in 1 : length(coeff.var[, 1])){
#'             ##' make the entries of these newly created matrices to be the 
#'             ##' weighted averages of the originally calculated datasets
#'             coeff.2[j, i] = sum(coeff.var[j, ] * w1[i, ], na.rm = T)
#'             incidence.2[j, i] = sum(incidence.per.1000.each.year[j, ] * w2[i, ], na.rm = T)
#'             # mbr2[j, i] =  sum(mean.br[j, ] * w2[i, ], na.rm = T)
#'             # mvacc2[j, i] = sum(mean.vac[j, ]* w2[i, ], na.rm = T) 
#'             mbr2[j, i] = mean.br[j, i] 
#'             mvacc2[j, i] = mean.vac[j, i]
#'         }
#'     }
#'     
#'     ##' set the original data to be equal to the weighted data
#'     coeff.var.cases = coeff.2
#'     incidence.per.1000 = incidence.2[1:j,1:i]
#'     # mean.br = mbr2[1:j,10:ncol(mean.br)]
#'     # mean.vac = mvacc2[1:j,10:ncol(mean.vac)]
#'     # 
#'     mean.br = mbr2[1:j,1:28]
#'     mean.vac = mvacc2[1:j,1:28]
#'     
#'     
#'     return(list(coeff.var.cases, incidence.per.1000, mean.br, mean.vac))
#' }

make.weighted.data <- function(coeff.var, 
                               incidence.per.1000.each.year,
                               mean.br, 
                               mean.vac,
                               w1,
                               w2){
  
  coeff.2 = coeff.var
  incidence.2 = incidence.per.1000.each.year
  mbr2 = mean.br
  mvacc2 = mean.vac
  for(i in 1 : length(coeff.var[1, ])){
    for(j in 1 : length(coeff.var[, 1])){
      ##' make the entries of these newly created matrices to be the weighted averages of the originally calculated datasets
      coeff.2[j, i] = sum(coeff.var[j, ] * w1[i, ], na.rm = T)
      incidence.2[j, i] = sum(incidence.per.1000.each.year[j, ] * w2[i, ], na.rm = T)
      mbr2[j, i] =  sum(mean.br[j, i] * w2[i, ], na.rm = T)
      mvacc2[j, i] = sum(mean.vac[j, i]* w2[i, ], na.rm = T) 
      
      
      
      
    }
  }
  
  ##' set the original data to be equal to the weighted data
  coeff.var.cases = coeff.2
  incidence.per.1000 = incidence.2
  mean.br = mbr2
  mean.vac = mvacc2
  
  return(list(coeff.2, incidence.2, mbr2, mvacc2))
}





#' Calculate the coefficient of variation and the mean number of cases, along with the 
#' mean birth rate and vaccination proportion over each period of years specified by the 
#' window length entry of the function
#'
#' @param year - year to begin analysis from. set to 1980 for analysis in the paper
#' @param window.length - the number of years which is considered at once. set to 10 for papaer
#' @param num.windows - the number of windows used in the analysis. this is calculated as the number of windows needed to go from the 
#' starting year of analysis, to the last year in the data
#' @param interp.subset.data - data set which has the cases data in it
#' @param interp.subset.pop - data set which has the population data in it
#' @param interp.subset.br - data set which has the birth rate data in it
#' @param interp.subset.vacc - data set which has the vaccination data in it
#'
#' @return - list(coeff.var, incidence.per.1000, mean.br, mean.vac) - data sets which are averaged over the window length



prepare.data.for.weighting <- function(year, window.length, num.windows,
                                       interp.subset.data,
                                       interp.subset.pop,
                                       interp.subset.br,
                                       interp.subset.vacc,
                                       subset.data){
  
  ##' set up data sets that will contain the desired data
  coeff.var = matrix(0, length(interp.subset.data[ , 1]), num.windows)
  incidence.per.1000 = matrix(0, length(interp.subset.data[ , 1]), num.windows)
  mean.br = matrix(0, length(interp.subset.data[ , 1]), num.windows)
  mean.vac = matrix(0, length(interp.subset.data[ , 1]), num.windows)
  
  
  ##' do calculations that calculate the coefficient of variation, incidence per 100, mean birth rate and
  ##' mean vaccination rate over periods of length given by the window length.
  for ( j in 1 : num.windows){
    for ( i in 1 : length(subset.data[ , 1])){
      coeff.var[i, j]  =  sd(interp.subset.data[i, paste(seq(year, year + window.length - 1))], na.rm = TRUE) /  
        mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      if(is.na( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)) == FALSE ){
        if( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) == 0) {
          coeff.var[i, j]  =  0
        } 
      }
      incidence.per.1000[i, j]  =  sum(as.numeric(interp.subset.data[i,  paste(seq(year, year + window.length - 1))]) / 
                                         as.numeric(interp.subset.pop[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) * 1000
      if(length(which(is.na(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))])))) < window.length){
        mean.br[i, j]  =  mean(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      }
      if(length(which(is.na(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))])))) < window.length){
        mean.vac[i, j]  =  mean(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      }
    }
    year = year + 1
  }
  
  incidence.per.1000.each.year = matrix(0, nrow(subset.data), length(seq(1980, 2017)))
  for(i in 1 : nrow(subset.data)){
    incidence.per.1000.each.year[i, ] = 1000 * as.numeric(interp.subset.data[i,  paste(seq(1980, 2017))]) / 
      as.numeric(interp.subset.pop[i, paste(seq(1980, 2017))])
    
  }
  
  
  return(list(coeff.var, incidence.per.1000.each.year, mean.br, mean.vac))
}


#' Calculate the coefficient of variation and the mean number of cases, along with the 
#' mean birth rate and vaccination proportion over each period of years specified by the 
#' window length entry of the function
#'
#' @param year - year to begin analysis from. set to 1980 for analysis in the paper
#' @param window.length - the number of years which is considered at once. set to 10 for papaer
#' @param num.windows - the number of windows used in the analysis. this is calculated as the number of windows needed to go from the 
#' starting year of analysis, to the last year in the data
#' @param interp.subset.data - data set which has the cases data in it
#' @param interp.subset.pop - data set which has the population data in it
#' @param interp.subset.br - data set which has the birth rate data in it
#' @param interp.subset.vacc - data set which has the vaccination data in it
#'
#' @return - list(coeff.var, incidence.per.1000, mean.br, mean.vac) - data sets which are averaged over the window length

calc.variables.by.window.length <- function(year, window.length, num.windows,
                                            interp.subset.data,
                                            interp.subset.pop,
                                            interp.subset.br,
                                            interp.subset.vacc){
  
  ##' set up data sets that will contain the desired data
  coeff.var = matrix(0, length(interp.subset.data[ , 1]), num.windows)
  incidence.per.1000 = matrix(0, length(interp.subset.data[ , 1]), num.windows)
  mean.br = matrix(0, length(interp.subset.data[ , 1]), num.windows)
  mean.vac = matrix(0, length(interp.subset.data[ , 1]), num.windows)
  
  
  ##' entry paste(seq(year, year + window.length - 1)) specifies which years of data we are considering
  ##' at each position
  for ( j in 1 : num.windows){
    for ( i in 1 : nrow(interp.subset.data)){
      ##' coefficient of variation is the standard deviation of the cases divided by the mean number of cases
      coeff.var[i, j]  =  sd(interp.subset.data[i, paste(seq(year, year + window.length - 1))], na.rm = TRUE) /  
        mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      ##' if cases are missing over a period, then report 0 as the coefficient of variation
      if(is.na( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)) == FALSE ){
        if( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) == 0) {
          coeff.var[i, j]  =  0
        } 
      }
      ##' calculate the 
      # incidence.per.1000[i, j]  =  sum(as.numeric(interp.subset.data[i,  paste(seq(year, year + window.length - 1))]) / 
      #                                      as.numeric(interp.subset.pop[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) * 1000
      # if(length(which(is.na(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))])))) < window.length){
      #     mean.br[i, j]  =  mean(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      # }
      # if(length(which(is.na(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))])))) < window.length){
      #     mean.vac[i, j]  =  mean(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      # }
    }
    ##' increment the year so that we are stepping our analysis forward one year at a time
    year = year + 1
  }
  
  mean.vac = data.frame(interp.subset.vacc[, -(1:2)])
  mean.br = data.frame(interp.subset.br[, -(1:2)]) 
  for(i in 1 : ncol(mean.vac)){
    mean.vac[, i] = mean.vac[, i] %>% as.character %>% as.numeric
    mean.br[, i] = mean.br[, i] %>% as.character %>% as.numeric
  }
  
  
  return(list(coeff.var, incidence.per.1000, mean.br, mean.vac))
}





#' Title
#'
#' @param d1 - incidence and cv data
#' @param d2 - estimated number of susceptibles by country and by year
#' @param d3 - estimated number of individuals by age and country and year
#' @param rep.cases - are we using the reported or estimated cases
#' @param regions - the WHO regions we are considering
#' @param make.inc.cv.scale.same - if 1, then we scale both the incidence and cv data to be on 0-1 range
#' Generally redundant now due to the log.incidence argument
#' @param connect.canonical.path.to.zero - do we attach the canonical path to the (0,0) point
#' @param log.incidence - should we log the incidence
#' 
#' @return
#' @export
#'
#' @examples
plots.for.sus.dist.by.canonical.path <- function(d1, d2, d3, rep.cases,
                                                 regions,
                                                 make.inc.cv.scale.same,
                                                 sqrt.inc,
                                                 connect.canonical.path.to.zero,
                                                 log.incidence,
                                                 make.figure.plot = F){
  
  ##' use the incidence and cv data to create the canonical path
  list[d1, canonical.path] = closest.path.point(d = d1, 
                                                use.rep.cases = rep.cases,
                                                regions = regions, 
                                                years = seq(1990, 2017),
                                                make.inc.cv.scale.same = make.inc.cv.scale.same,
                                                sqrt.inc = sqrt.inc,
                                                connect.canonical.path.to.zero,
                                                log.incidence = log.incidence,
                                                make.figure.plot = make.figure.plot)
  
  d3 <- d3[3:ncol(d3)]
  colnames(d3) <- paste("pop", 0:59, sep = "")
  d4 <- cbind(d2, d3) #data.frame for age profile of susceptiblity
  
  list[subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data] = get.data.for.animation(regions)
  d2$closest = NA
  d2$population = NA
  for(i in 1 : nrow(d2)){
    c = d2$Country[i]
    y = d2$Year[i]
    k = which(d1$Country == c & d1$Year == y )
    k2 = which(subset.pop.by.year$Country.Name == c)
    if(length(k) > 0){
      d2$closest[i] = d1$closest[k]
      if(y == 2014){
        d2$population[i] = d2$population[i-1]
      }else{
        d2$population[i] = subset.pop.by.year[k2, paste("X", y, sep = "")]
      }
      
    }
    
  }
  
  sample.size = 100000
  for(i in 0:59){
    d2[, paste("X", i, sep = "")] = as.numeric(as.character(d2[, paste("X", i, sep = "")]))
  }
  
  d2$population = d2$population %>% as.character %>% as.numeric
  for(i in 1 : max(d2$closest, na.rm = T)){
    qq = d2 %>%
      filter(closest == i) 
    if(nrow(qq)>0){
      qq = qq[, paste("X", 0:59, sep = "")]
      qq = qq %>% colSums(na.rm = T)
      qq = qq %>% divide_by(max(1,min(qq))) %>%
        multiply_by(20) %>%
        round
      
      pp = matrix(0, sum(qq), 1)
      count = 1
      for(j in 1 : length(qq)){
        pp[count:(count + qq[j] - 1)] = j-1
        count = count + qq[j]
      }
      if(i==1){
        set.seed(100)
        dist.by.age = cbind(1, sample(pp, sample.size, replace = T))
      }
      if(i > 1){
        set.seed(100)
        dist.by.age = rbind(dist.by.age, cbind(i, sample(pp, sample.size, replace = T)))
      }
      head(dist.by.age)
    }
  }
  prop.sus <- data.frame(matrix(NA, max(d2$closest, na.rm = T), 2))
  prop.sus[, 1] = seq(1, max(d2$closest, na.rm = T))
  
  for(i in 1 : max(prop.sus[, 1])){
    qq = d2 %>%
      filter(closest == i) %>%
      na.omit
    if(nrow(qq) >0){
      num.sus = qq[ paste("X", 0:59, sep = "")] %>% sum
      total.pop = qq$population %>% sum
      prop.sus[i, 2] = num.sus / total.pop
    }
  }
  #prop.sus %<>% na.omit
  colnames(prop.sus) = c("x", "y")
  
  
  
  dist.by.age = data.frame(dist.by.age)
  colnames(dist.by.age) = c("point", "age")
  num.points = nrow(canonical.path)
  missing.points = setdiff(1:num.points, unique(dist.by.age$point))
  if(length(missing.points) > 0){
    missing.boxes = data.frame(cbind(missing.points, -5))   
    colnames(missing.boxes) = c("point", "age")
    dist.by.age = rbind(dist.by.age, missing.boxes)
    
  }
  
  #age profile of susceptiblity
  d4$closest = NA
  for(i in 1 : nrow(d4)){
    c = d4$Country[i]
    y = d4$Year[i]
    k = which(d1$Country == c & d1$Year == y )
    if(length(k) > 0){
      d4$closest[i] = d1$closest[k]
    }
  } #this adds the point of the canonial path for each country and year in d2 (i.e., row)
  
  for(i in 0:59){
    d4[, paste("X", i, sep = "")] = as.numeric(as.character(d4[, paste("X", i, sep = "")]))
    d4[, paste("pop", i, sep = "")] = as.numeric(as.character(d4[, paste("pop", i, sep = "")]))
  } #just changing numbers to numeric
  
  prop.sus.by.age <- data.frame(matrix(NA, max(d4$closest, na.rm = T), 61))
  prop.sus.by.age[, 1] = seq(1, max(d4$closest, na.rm = T)) #matrix of nrow= 38, all positions, and 2 columns
  for(i in 1 : max(prop.sus.by.age[, 1])){ #for each position
    qq = d4 %>%
      filter(closest == i) %>%
      na.omit  #rows with this canonical point
    if(nrow(qq) >0){
      num.sus = qq[ paste("X", 0:59, sep = "")] %>% colSums
      total.pop = qq[ paste("pop", 0:59, sep = "")] %>% colSums
      prop.sus.by.age[i, 2:61] = num.sus / total.pop
    }
  }
  colnames(prop.sus.by.age) = c("X", paste("prop.sus", 0:59, sep = ""))
  
  
  return(list( dist.by.age, canonical.path, prop.sus, d1, prop.sus.by.age))
}



#' Title
#'
#' @param data - the data set used for the GAM 
#' @param ticks - points on the y-axis that we want our ticks to appear at
#' @param tick.labels - the lables we want at these tick positions
#'
#' @return
#' @export
#'
#' @examples
produce.gams.plot <- function(data, 
                              ticks = c(0,sqrt(10), sqrt(50), 10, sqrt(200)), 
                              tick.labels = c(0,10,50,100,200)){
  
  ##' fit a GAM to the incidence and the coefficient of variation with the independent variable being birth rate
  A <- gam(sqrt(Incidence) ~ s(Mean.birth.rate), data = data )
  A2 <- gam(Coefficient.of.Variation ~s(Mean.birth.rate) , data = data )
  
  ##' find all the (integer) values of birth rate that are in the data set
  brs = unique(c(floor(data$Mean.birth.rate), ceiling(data$Mean.birth.rate)))
  
  ##' order all of these birth rates to be used in prediction
  brs = brs[order(brs)]
  
  ##' predict the incidence, given the birth rate
  g1 = predict.gam(A, newdata = data.frame(Mean.birth.rate = brs))
  
  ##' predict the coefficient of variation given the birth rate
  g2 = predict.gam(A2, newdata = data.frame(Mean.birth.rate = brs))
  
  ##' collect these two predictions, along with the independent variable in a data frame
  h = data.frame(cbind(as.matrix(as.numeric(g2)),as.matrix(as.numeric( g1)), brs))
  colnames(h) = c("x", "y", "num")
  
  ##' fit a GAM to the incidence and the coefficient of variation with the independent variable being vaccination 
  C <- gam(sqrt(Incidence) ~s(Mean.vaccination), data = data )
  C2 <- gam(Coefficient.of.Variation ~s(Mean.vaccination), data = data )
  
  ##' find the minimum and maximum values of the vaccination
  vac.min = min(data$Mean.vaccination, na.rm = T)
  vac.max = max(data$Mean.vaccination, na.rm = T)
  
  ##' predict the incidence, given the vaccination proportion
  s1 = predict.gam(C, newdata = data.frame(Mean.vaccination = vac.min:vac.max))
  
  ##' predict the coefficient of variation given the vaccination proportion
  s2 = predict.gam(C2, newdata = data.frame(Mean.vaccination = vac.min:vac.max))
  
  ##' collect these two predictions, along with the independent variable in a data frame
  z = data.frame(cbind(as.matrix(as.numeric(s2)),as.matrix(as.numeric( s1)), vac.min:vac.max))
  colnames(z) = c("x", "y", "num")
  
  ##' fit a GAM to the incidence and the coefficient of variation with the independent variable being (birth rate)*(1-vacc proportion)
  D <- gam(sqrt(Incidence) ~s(br.vacc), data = data )
  D2 <- gam(Coefficient.of.Variation ~s(br.vacc), data = data )
  
  ##' find the minimum and maximum values of this independent variable
  br.vac.min = min(data$br.vacc, na.rm = T)
  br.vac.max = max(data$br.vacc, na.rm = T)
  
  ##' predict the incidence given the independent variable
  f1 = predict.gam(D, newdata = data.frame(br.vacc = br.vac.min:br.vac.max))
  
  ##' predict the coefficient of variation given the independent variable
  f2 = predict.gam(D2, newdata = data.frame(br.vacc = br.vac.min:br.vac.max))
  
  ##' collect these two predictions, along with the independent variable in a data frame
  j = data.frame(cbind(as.matrix(as.numeric(f2)),as.matrix(as.numeric( f1)), br.vac.min:br.vac.max))
  colnames(j) = c("x", "y", "num")
  
  ##' add a column to each of the prediction data sets we have created to label what the independent variable was
  h = cbind(h, 'br')
  colnames(h) = c("x", "y", "num", 'type')
  z = cbind(z, 'vacc')
  colnames(z) = c("x", "y", "num", "type")
  j = cbind(j, 'br.vacc')
  colnames(j) = c("x", "y", "num", "type")
  
  ##' combine these 3 data sets into 1 for plotting
  z1 = rbind(z, h, j)
  
  ##' generate colors
  colfunc.br <- colorRampPalette(c("mediumaquamarine", "darkgreen"))
  colfunc.vacc = colorRampPalette(c("lightblue", "navyblue"))
  colfunc.br.vacc = colorRampPalette(c("pink", "firebrick"))
  
  
  
  return(list(h, z, j, z1, colfunc.br, colfunc.vacc, colfunc.br.vacc))
  
}




#' @param input.data - data set to add mcv2 data to
#' @param non.mcv2.alpha - tranparency score for countries which have not introduced mcv2
#' @param mcv2.alpha - tranparency score for countries which have introduced mcv2

add.mcv2.data.to.anim.data <- function(input.data, non.mcv2.alpha = 0.4, mcv2.alpha = 0.6){
  ##' read in the MCV2 data. this data set gives the year of introduction of second measles routine vaccination
  mcv2 = read.csv("data/MCV2_introduction.csv", stringsAsFactors = F)
  ##' change entries for the year introduction of second measles routine vaccination which are n/a to NA
  j = which(mcv2$MCV2_Intro_Year == "n/a")
  mcv2$MCV2_Intro_Year[j] = NA
  mcv2$MCV2_Intro_Year = as.numeric(as.character(mcv2$MCV2_Intro_Year))
  ##' change entries for the part of the country mcv2 was introduced in from n/a to NA
  j = which(mcv2$MCV2_Intro_part_of_country == "n/a" |
              mcv2$MCV2_Intro_part_of_country == "n/d")
  mcv2$MCV2_Intro_part_of_country[j] = NA
  mcv2$MCV2_Intro_part_of_country = as.numeric(as.character(mcv2$MCV2_Intro_part_of_country))
  
  ##' set the data set which will be output from this function to be the same as the input data set
  output.data = input.data
  ##' add a column which will hold the year of introduction of the MCV2 for each country
  output.data$MCV2 = as.numeric(matrix(0, nrow(output.data), 1))
  
  ##' which countries are in the MCV2 data set
  aa = unique(mcv2$Country)
  ##' loop over these countries adding the year of introduction of MCV2 by country to the output data 
  for(i in 1 : length(aa)){
    year = mcv2$MCV2_Intro_Year[i]
    if(!is.na(year)){
      j = which(output.data$Country == aa[i])
      k = which(output.data$Year[j] >= year)
      output.data$MCV2[j[k]] = 1
    }
  }
  
  
  
  output.data$alpha = as.numeric(matrix(non.mcv2.alpha, nrow(output.data), 1))
  output.data$shape = as.numeric(matrix(1, nrow(output.data), 1))
  j = which(output.data$MCV2 == 1)
  output.data$alpha[j] = mcv2.alpha
  output.data$shape[j] = 2
  output.data$alpha = as.numeric(as.character(output.data$alpha))
  return(output.data)
}





plot.world.map <- function(d, 
                           use.rep.cases,
                           regions,
                           year,
                           make.inc.cv.scale.same,
                           sqrt.inc,
                           missing.countries,
                           replace.countries,
                           connect.canonical.path.to.zero,
                           log.incidence,
                           with.text =1){
  
  list[d1, canonical.path] = closest.path.point(d = d, 
                                                use.rep.cases = use.rep.cases,
                                                regions = regions, 
                                                years = seq(1990, 2017),
                                                make.inc.cv.scale.same = make.inc.cv.scale.same,
                                                sqrt.inc = sqrt.inc,
                                                connect.canonical.path.to.zero,
                                                log.incidence = log.incidence)
  
  
  A = d1%>%filter(Year == year)    
  data(World)
  World$pos = NA
  for(i in 1 : nrow(A)){
    if(A$Country[i] %in% missing.countries){
      k = which(missing.countries == A$Country[i])
      A$Country[i] = replace.countries[k]
    }
    j = which(World$name == A$Country[i])
    if(length(j) > 0){
      World$pos[j] = A$closest[i]
    }
  }
  #Replace Somaliland with Somalia's position
  World$pos[(which(World$iso_a3=="SOL"))] <- World$pos[(which(World$iso_a3=="SOM"))] 
  World$iso_a3[(which(World$iso_a3=="SOL"))] <-  NA
  #Replace Taiwan with China's position
  World$pos[(which(World$iso_a3=="TWN"))] <- World$pos[(which(World$iso_a3=="CHN"))] 
  World$iso_a3[(which(World$iso_a3=="TWN"))] <-  NA
  #Replace Peurto Rico with USA's position
  World$pos[(which(World$iso_a3=="PRI"))] <- World$pos[(which(World$iso_a3=="USA"))] 
  World$iso_a3[(which(World$iso_a3=="PRI"))] <-  NA
  
  
  cols.1 <- colorRampPalette(c("red", "white", "blue"))(nrow(canonical.path))
  
  kk = which(World$iso_a3 == "ATA" | World$iso_a3 == "GRL")
  World$cols = cols.1[World$pos]
  
  #Fill in color NAs with grey.  These are countries for which we do not have data
  World$cols[is.na(World$cols)] <- "#969696"
  
  if(with.text==1){
    p<-tm_shape(World[-(kk),]) +
      tm_polygons("cols", textNA="No Data", 
                  title="Well-Being Index") +
      tm_text("iso_a3", size="area", root=5) 
  }else{p<-tm_shape(World[-(kk),]) +
    tm_polygons("cols", textNA="No Data", 
                title="Well-Being Index")  }
  
  #print(p)
  return(list(p, A, canonical.path))
}




get.data.for.animation.state.space <- function(regions){
  
  
  cases.by.country.by.year = read.csv("data/Measles_cases_by_year2.csv", stringsAsFactors = FALSE)
  subset.data = read.csv("data/State_space_cases_new.csv")
  Birth.rates = read.csv("data/Birth_rates.csv", stringsAsFactors = FALSE)
  pop.by.year = read.csv("data/All_populations.csv", stringsAsFactors = FALSE)
  vacc.rates = read.csv("data/Measles_vac_all.csv", stringsAsFactors = FALSE)
  
  subset.data$Country = NA
  subset.data$WHO_REGION = NA
  for(i in 1 : nrow(subset.data)){
    iso = subset.data$iso[i]
    k = which(cases.by.country.by.year$ISO_code == iso)
    if(length(k) > 0){
      subset.data$Country[i] = cases.by.country.by.year$Cname[k]
      subset.data$WHO_REGION[i] = cases.by.country.by.year$WHO_REGION[k]
    }
  }
  
  subset.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION %in% regions)
  subset.vaccination = subset(vacc.rates, vacc.rates$WHO_REGION %in% regions)
  subset.birth.rates = subset(Birth.rates, Birth.rates$WHO_REGION %in% regions)
  # 
  # for(i in 1 : nrow(subset.data)){
  #     subset.data$Country[i] = cases.by.country.by.year$Cname[which(cases.by.country.by.year$Cname == subset.data$Country[i])]
  #     subset.data$WHO_REGION[i] = cases.by.country.by.year$WHO_REGION[which(cases.by.country.by.year$Cname == subset.data$Country[i])]
  # }
  # 
  subset.data = subset(subset.data, subset.data$WHO_REGION %in% regions)
  
  missing3 = setdiff(subset.vaccination$Country, subset.data$Country)
  if(length(missing3) > 0){
    j = which(subset.vaccination$Country %in% missing3)
    subset.vaccination = subset.vaccination[-j, ]
  }
  
  missing1 = setdiff(subset.vaccination$Country, subset.pop.by.year$Country.Name)
  if(length(missing1) > 0){
    j = which(subset.vaccination$Country %in% missing1)
    subset.vaccination = subset.vaccination[-j, ]
  }
  missing2 = setdiff(subset.birth.rates$Country, subset.vaccination$Country)
  if(length(missing2) > 0){
    j = which(subset.birth.rates$Country %in% missing2)
    subset.birth.rates = subset.birth.rates[-j, ]
  }
  missing3 = setdiff(subset.data$Country, subset.vaccination$Country)
  if(length(missing3) > 0){
    j = which(subset.data$Country %in% missing3)
    subset.data = subset.data[-j, ]
  }
  
  missing3 = setdiff(subset.pop.by.year$Country.Name, subset.data$Country)
  if(length(missing3) > 0){
    j = which(subset.pop.by.year$Country.Name %in% missing3)
    subset.pop.by.year = subset.pop.by.year[-j, ]
  }
  p1  =  subset.pop.by.year
  p2  =  subset.vaccination
  p3  =  subset.birth.rates
  p4  =  subset.data
  
  
  for ( i in 1 : length(subset.vaccination[, 1])){
    C  =  subset.vaccination$Country[i]
    p2[i, ]  =  subset.vaccination[i, ]
    j = which(subset.pop.by.year$Country.Name == C)
    p1[i, ]  =  subset.pop.by.year[j, ]
    j = which(subset.birth.rates$Country == C)
    p3[i, ]  =  subset.birth.rates[j, ]
    j = which(subset.data$Country == C)
    p4[i, ]  =  subset.data[j, ]
  }
  subset.pop.by.year = p1
  subset.vaccination = p2
  subset.birth.rates = p3
  subset.data = p4
  subset.data = subset.data[, -(1)]
  x = colnames(subset.data)[which(grepl("X", colnames(subset.data)))]
  x = min(gsub("X", "", x)):max(gsub("X", "", x))
  colnames(subset.data) = c(x,"Country", "WHO_REGION")
  
  return(list(subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data, x))
}



#############################################################################################################

generate.state.space.data <- function(window.length, regions, 
                                      gaussian.st.dev, cutoff = 50, 
                                      interp.resolution = 20,
                                      year.shift.inc = 3){
  
  require(stats)
  list[subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data, x] =
    get.data.for.animation.state.space(regions)
  
  
  ##' interpolate the datasets to have entries for all points in time once the interpolation is done.
  list[interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop] = 
    interp.datasets.state.space(subset.data, 
                                subset.vaccination, 
                                subset.birth.rates, 
                                subset.pop.by.year,
                                x,
                                xout = x)
  
  # list[mean.cases, coeff.var.cases, incidence.per.1000, mean.br, mean.vac] = 
  #     prepare.matrices.for.animation(interp.subset.data, subset.data)
  
  ##' number of unique years that we will have data for. The longer the window, the less unique years of data.
  num.windows = length(x) - window.length + 1
  
  ##' first year of data
  year = x[1]
  
  ##' setting up the datasets
  coeff.var = matrix(0, length(subset.data[ , 1]), num.windows)
  incidence.per.1000 = matrix(0, length(subset.data[ , 1]), num.windows)
  mean.br = matrix(0, length(subset.data[ , 1]), num.windows)
  mean.vac = matrix(0, length(subset.data[ , 1]), num.windows)
  
  ##' do calculations that calculate the coefficient of variation, incidence per 100, mean birth rate and
  ##' mean vaccination rate over periods of length given by the window length.
  for ( j in 1 : num.windows){
    for ( i in 1 : length(subset.data[ , 1])){
      coeff.var[i, j]  =  sd(interp.subset.data[i, paste(seq(year, year + window.length - 1))], na.rm = TRUE) /  
        mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      if(is.na( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)) == FALSE ){
        if( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) == 0) {
          coeff.var[i, j]  =  0
        } 
      }
      incidence.per.1000[i, j]  =  sum(as.numeric(interp.subset.data[i,  paste(seq(year, year + window.length - 1))]) / 
                                         as.numeric(interp.subset.pop[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) * 1000
      if(length(which(is.na(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))])))) < window.length){
        mean.br[i, j]  =  mean(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      }
      if(length(which(is.na(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))])))) < window.length){
        mean.vac[i, j]  =  mean(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
      }
    }
    year = year + 1
  }
  
  incidence.per.1000.each.year = matrix(0, nrow(subset.data), length(x))
  for(i in 1 : nrow(subset.data)){
    incidence.per.1000.each.year[i, ] = 1000 * as.numeric(interp.subset.data[i,  paste(x)]) / 
      as.numeric(interp.subset.pop[i, paste(x)])
    
  }
  
  ##' we take the weighted average of the values that we have calculated for each year, where the weights are gaussian,
  ##' with a specified number of years for the standard deviation
  x1 = seq(1, length(coeff.var[1, ]) + 1)
  w1 = matrix(0, length(x1), length(x1))
  for (i in 1 : length(x1)){
    ##' set up the gaussian weights for averaging
    
    w.input = x1 - x1[i]
    w1[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
    
    ##' make sure that the weights add up to 1 for each of the specific weightings
    j = which(w.input == 0)
    w1[i,j : length(w1[i,])] = 0
    w1[i, ] = w1[i, ] / sum(w1[i, ])
  }
  
  w1 = w1[-(1), ]
  w1 = w1[, -(ncol(w1)) ]
  
  x2 = seq(1, length(incidence.per.1000.each.year[1, ]))
  w2 = matrix(0, length(x2), length(x2))
  for (i in 1 : length(x2)){
    ##' set up the gaussian weights for averaging
    
    w.input = x2 - x2[i] + year.shift.inc
    w2[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
    
    ##' make sure that the weights add up to 1 for each of the specific weightings
    j = which(w.input == (year.shift.inc + 1))
    if(length(j) > 0){
      w2[i,j : length(w2[i,])] = 0
      w2[i, ] = w2[i, ] / sum(w2[i, ])
    }
    
  }
  #w2 = w2[-(1), ]
  #w2 = w2[, -(ncol(w2)) ]
  w2 = w2[-(1:9), ]
  
  ##' make a set of matrices that are the same size as the matrices containing the data.
  coeff.2 = coeff.var
  incidence.2 = incidence.per.1000
  mbr2 = mean.br
  mvacc2 = mean.vac
  for(i in 1 : length(coeff.var[1, ])){
    for(j in 1 : length(coeff.var[, 1])){
      ##' make the entries of these newly created matrices to be the weighted averages of the originally calculated datasets
      coeff.2[j, i] = sum(coeff.var[j, ] * w1[i, ], na.rm = T)
      incidence.2[j, i] = sum(incidence.per.1000.each.year[j, ] * w2[i, ], na.rm = T)
      mbr2[j, i] =  sum(mean.br[j, i] * w2[i, ], na.rm = T)
      mvacc2[j, i] = sum(mean.vac[j, i]* w2[i, ], na.rm = T)
      
      ##' Should we do weighted average of birth rate and vaccination rate?
      ##' If so uncomment the next two lines
      
      # mbr2[j, i] = sum(mean.br[j, ] * w1[i, ], na.rm = T)
      # mvacc2[j, i] = sum(mean.vac[j, i] * w1[i, ], na.rm = T)
    }
  }
  
  ##' set the original data to be equal to the weighted data
  coeff.var.cases = coeff.2
  incidence.per.1000 = incidence.2
  mean.br = mbr2
  mean.vac = mvacc2
  
  
  ##' set up the timeline on which we do the interpolation. 
  ##' The number of sections that the yearly data is split up to is given by interp.resolution  
  x = seq(1981 + (window.length - 1), max(x))
  xout = seq(1981+ (window.length - 1), max(x), 1/interp.resolution)
  
  ##' interpolate the data and add columns that contain the correspondin country and WHO region of each line
  
  coeff.var.cases = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(coeff.var.cases, x, xout))
  incidence.per.1000 = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(incidence.per.1000, x, xout))
  mean.br = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.br, x, xout))
  mean.vac = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.vac, x, xout))
  
  ##' round the data to 2 decimal places for ease of reading.
  mean.vac[, -(1:2)] = as.numeric(mean.vac[, -(1:2)])
  mean.br[, -(1:2)] = as.numeric(mean.br[, -(1:2)])
  incidence.per.1000[, -(1:2)] = as.numeric(incidence.per.1000[, -(1:2)])
  coeff.var.cases[, -(1:2)] = as.numeric(coeff.var.cases[, -(1:2)])
  
  
  ##' set up the output to be the appropriate size and add column labels.
  ##' Additionally add enough room to include additional data for each year that will be used 
  ##' to calibrate the data for each year, so that the minimum and maximum of vaccination rate is 0 and 100 each time.
  ##' This ensures that the colour scale is constant
  
  output.data = matrix(0, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + (2 * length(regions) * length(xout)), 7)
  output.data  =  data.frame(output.data)
  colnames(output.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
  
  ##' input the appropriate data to the outputs 
  output.data$Country[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 1], length(coeff.var.cases[1, -(1:2)]))
  output.data$WHO_REGION[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 2], length(coeff.var.cases[1, -(1:2)]))
  count = 1
  for(i in 3 : length(coeff.var.cases[1, ])){
    output.data$Coefficient.of.Variation[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = coeff.var.cases[, i]
    output.data$Incidence[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = incidence.per.1000[, i]
    output.data$Mean.vaccination[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.vac[, i]
    output.data$Mean.birth.rate[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.br[, i]
    output.data$Year[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = xout[i - 2]  
    count = count + 1
  }
  
  ##' add in the dummy data for each year to keep the scales constant.
  year.mins = matrix(0, length(xout), 2)
  
  for(i in 1 : length(xout)){
    t  =  subset(output.data, output.data$Year ==  unique(xout)[i])
    year.mins[i, 1] = xout[i]
    year.mins[i, 2] = as.numeric(min(t$Mean.birth.rate,na.rm = T)  )
  }
  
  l =  expand.grid("", -1, 0, c(0,100), 0, xout, regions)
  
  for(i in 1 : (2 * length(regions) * length(xout))){
    y = l[i, 6]
    j = which(year.mins[, 1] == y)
    l[i, 5]  =  year.mins[j, 2]
  }
  
  output.data[((length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + 1 ): length(output.data[, 1]), ]  =  l
  
  for( i in 1 : (2 * length(regions) * length(xout))){
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] = regions[as.numeric(output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] )]
    output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)])  + i, 1] = ""
  }
  
  ##' make sure that each column that should be numeric is numeric.
  output.data$Coefficient.of.Variation = as.numeric(output.data$Coefficient.of.Variation)
  output.data$Incidence  =  as.numeric(output.data$Incidence)
  output.data$Mean.vaccination  =  as.numeric(output.data$Mean.vaccination)
  output.data$Mean.birth.rate  =  as.numeric(output.data$Mean.birth.rate)
  output.data$Year   =  as.numeric(output.data$Year)
  output.data$Coefficient.of.Variation[which(output.data$Coefficient.of.Variation == "Inf")] = 0
  return(output.data)
}







#' Return Gaussian weighted mean and cv from state-space model estimates
#'
#' @param df - data frame incidene per 100000 from 1981 to 2017 in columns and countries in rows
#'
#' @return list(incs.all, coeff.all)
return.cv.inc.from.ss <- function(df){
  gaussian.st.dev = 3
  cutoff = 50
  
  output.weights.gaussian.with.cutoff <- function(x, st.dev, cutoff, neg.only = F){
    
    weights = dnorm(x, mean = 0, sd = st.dev)
    if(neg.only == T){
      weights[which(x > cutoff)] = 0
    } else{
      weights[which(abs(x) > cutoff)] = 0
    }
    return(weights)
  }
  
  coeff.all <- incs.all <- matrix(NA, nrow(df), 28)
  for (c in 1:nrow(df)){
    
    infs.by.year <- as.numeric(df[c,-1])
    
    cv.inc = matrix(NA, length(10:length(infs.by.year)), 1)
    count = 1
    for(i in 10:length(infs.by.year)){
      cv.inc[count, 1] = sd(infs.by.year[(i-9):i])/mean(infs.by.year[(i-9):i])
      count = count + 1
    }
    
    x1 = seq(1, length(cv.inc[, 1]) + 1)
    w1 = matrix(0, length(x1), length(x1))
    for (i in 1 : length(x1)){
      ##' set up the gaussian weights for averaging
      
      w.input = x1 - x1[i]
      w1[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
      
      ##' make sure that the weights add up to 1 for each of the specific weightings
      j = which(w.input == 0)
      w1[i,j : length(w1[i,])] = 0
      w1[i, ] = w1[i, ] / sum(w1[i, ])
    }
    
    w1 = w1[-(1), ]
    w1 = w1[, -(ncol(w1)) ]
    
    
    x2 = seq(1, length(infs.by.year) + 1)
    w2 = matrix(0, length(x2), length(x2))
    for (i in 1 : length(x2)){
      ##' set up the gaussian weights for averaging
      
      w.input = x2 - x2[i]
      w2[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
      
      ##' make sure that the weights add up to 1 for each of the specific weightings
      j = which(w.input == 1)
      if(length(j) > 0){
        w2[i,j : length(w2[i,])] = 0
      }
      w2[i, ] = w2[i, ] / sum(w2[i, ])
    }
    w2 = w2[-(1), ]
    #w2 = w2[, -(ncol(w2)) ]
    w2 = w2[-(1:9), ]
    
    #colfunc = colorRampPalette(c("grey","red"))
    
    coeff.2 = cv.inc
    incidence.2 = cv.inc
    for(i in 1 : length(coeff.2)){
      ##' make the entries of these newly created matrices to be the weighted averages
      ##'  of the originally calculated datasets
      w.input = seq(-9-(i-1), 0, 1)
      k1 = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)/
        sum(output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff))
      
      w.input = seq(-(i-1), 0, 1)
      k2 = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)/
        sum(output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff))
      
      coeff.2[i] = sum(cv.inc[1:length(k2)] * k2, na.rm = T)
      incidence.2[i] = sum(infs.by.year[1:length(k1)] * k1, na.rm = T)
      
    }
    
    coeff.all[c,] <- coeff.2
    incs.all[c,] <- incidence.2
  }
  
  return(list(incs.all, coeff.all))
}


##' Got online to allow different colors for line and point
##' https://stackoverflow.com/questions/19053440/r-legend-with-points-and-lines-being-different-colors-for-the-same-legend-item
LEGEND <- function (x, y = NULL, legend, fill = NULL, 
                    col = par("col"), pt.col=col, line.col=col,
                    border = "black", lty, lwd, pch, angle = 45, density = NULL,
                    bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"),
                    box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd,
                    xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0,
                                                                                0.5), text.width = NULL, text.col = par("col"), text.font = NULL,
                    merge = do.lines && has.pch, trace = FALSE, plot = TRUE,
                    ncol = 1, horiz = FALSE, title = NULL, inset = 0, xpd, title.col = text.col,
                    title.adj = 0.5, seg.len = 2)
{
  if (missing(legend) && !missing(y) && (is.character(y) ||
                                         is.expression(y))) {
    legend <- y
    y <- NULL
  }
  mfill <- !missing(fill) || !missing(density)
  if (!missing(xpd)) {
    op <- par("xpd")
    on.exit(par(xpd = op))
    par(xpd = xpd)
  }
  title <- as.graphicsAnnot(title)
  if (length(title) > 1)
    stop("invalid 'title'")
  legend <- as.graphicsAnnot(legend)
  n.leg <- if (is.call(legend))
    1
  else length(legend)
  if (n.leg == 0)
    stop("'legend' is of length 0")
  auto <- if (is.character(x))
    match.arg(x, c("bottomright", "bottom", "bottomleft",
                   "left", "topleft", "top", "topright", "right", "center"))
  else NA
  if (is.na(auto)) {
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2)
      stop("invalid coordinate lengths")
  }
  else nx <- 0
  xlog <- par("xlog")
  ylog <- par("ylog")
  rect2 <- function(left, top, dx, dy, density = NULL, angle,
                    ...) {
    r <- left + dx
    if (xlog) {
      left <- 10^left
      r <- 10^r
    }
    b <- top - dy
    if (ylog) {
      top <- 10^top
      b <- 10^b
    }
    rect(left, top, r, b, angle = angle, density = density,
         ...)
  }
  segments2 <- function(x1, y1, dx, dy, ...) {
    x2 <- x1 + dx
    if (xlog) {
      x1 <- 10^x1
      x2 <- 10^x2
    }
    y2 <- y1 + dy
    if (ylog) {
      y1 <- 10^y1
      y2 <- 10^y2
    }
    segments(x1, y1, x2, y2, ...)
  }
  points2 <- function(x, y, ...) {
    if (xlog)
      x <- 10^x
    if (ylog)
      y <- 10^y
    points(x, y, ...)
  }
  text2 <- function(x, y, ...) {
    if (xlog)
      x <- 10^x
    if (ylog)
      y <- 10^y
    text(x, y, ...)
  }
  if (trace)
    catn <- function(...) do.call("cat", c(lapply(list(...),
                                                  formatC), list("\n")))
  cin <- par("cin")
  Cex <- cex * par("cex")
  if (is.null(text.width))
    text.width <- max(abs(strwidth(legend, units = "user",
                                   cex = cex, font = text.font)))
  else if (!is.numeric(text.width) || text.width < 0)
    stop("'text.width' must be numeric, >= 0")
  xc <- Cex * xinch(cin[1L], warn.log = FALSE)
  yc <- Cex * yinch(cin[2L], warn.log = FALSE)
  if (xc < 0)
    text.width <- -text.width
  xchar <- xc
  xextra <- 0
  yextra <- yc * (y.intersp - 1)
  ymax <- yc * max(1, strheight(legend, units = "user", cex = cex)/yc)
  ychar <- yextra + ymax
  if (trace)
    catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra,
                                                   ychar))
  if (mfill) {
    xbox <- xc * 0.8
    ybox <- yc * 0.5
    dx.fill <- xbox
  }
  do.lines <- (!missing(lty) && (is.character(lty) || any(lty >
                                                            0))) || !missing(lwd)
  n.legpercol <- if (horiz) {
    if (ncol != 1)
      warning(gettextf("horizontal specification overrides: Number of columns := %d",
                       n.leg), domain = NA)
    ncol <- n.leg
    1
  }
  else ceiling(n.leg/ncol)
  has.pch <- !missing(pch) && length(pch) > 0
  if (do.lines) {
    x.off <- if (merge)
      -0.7
    else 0
  }
  else if (merge)
    warning("'merge = TRUE' has no effect when no line segments are drawn")
  if (has.pch) {
    if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L],
                                                      type = "c") > 1) {
      if (length(pch) > 1)
        warning("not using pch[2..] since pch[1L] has multiple chars")
      np <- nchar(pch[1L], type = "c")
      pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
    }
    if (!is.character(pch))
      pch <- as.integer(pch)
  }
  if (is.na(auto)) {
    if (xlog)
      x <- log10(x)
    if (ylog)
      y <- log10(y)
  }
  if (nx == 2) {
    x <- sort(x)
    y <- sort(y)
    left <- x[1L]
    top <- y[2L]
    w <- diff(x)
    h <- diff(y)
    w0 <- w/ncol
    x <- mean(x)
    y <- mean(y)
    if (missing(xjust))
      xjust <- 0.5
    if (missing(yjust))
      yjust <- 0.5
  }
  else {
    h <- (n.legpercol + (!is.null(title))) * ychar + yc
    w0 <- text.width + (x.intersp + 1) * xchar
    if (mfill)
      w0 <- w0 + dx.fill
    if (do.lines)
      w0 <- w0 + (seg.len + x.off) * xchar
    w <- ncol * w0 + 0.5 * xchar
    if (!is.null(title) && (abs(tw <- strwidth(title, units = "user",
                                               cex = cex) + 0.5 * xchar)) > abs(w)) {
      xextra <- (tw - w)/2
      w <- tw
    }
    if (is.na(auto)) {
      left <- x - xjust * w
      top <- y + (1 - yjust) * h
    }
    else {
      usr <- par("usr")
      inset <- rep_len(inset, 2)
      insetx <- inset[1L] * (usr[2L] - usr[1L])
      left <- switch(auto, bottomright = , topright = ,
                     right = usr[2L] - w - insetx, bottomleft = ,
                     left = , topleft = usr[1L] + insetx, bottom = ,
                     top = , center = (usr[1L] + usr[2L] - w)/2)
      insety <- inset[2L] * (usr[4L] - usr[3L])
      top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] +
                      h + insety, topleft = , top = , topright = usr[4L] -
                      insety, left = , right = , center = (usr[3L] +
                                                             usr[4L] + h)/2)
    }
  }
  if (plot && bty != "n") {
    if (trace)
      catn("  rect2(", left, ",", top, ", w=", w, ", h=",
           h, ", ...)", sep = "")
    rect2(left, top, dx = w, dy = h, col = bg, density = NULL,
          lwd = box.lwd, lty = box.lty, border = box.col)
  }
  xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1),
                                              rep.int(n.legpercol, ncol)))[1L:n.leg]
  yt <- top - 0.5 * yextra - ymax - (rep.int(1L:n.legpercol,
                                             ncol)[1L:n.leg] - 1 + (!is.null(title))) * ychar
  if (mfill) {
    if (plot) {
      if (!is.null(fill))
        fill <- rep_len(fill, n.leg)
      rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox,
            col = fill, density = density, angle = angle,
            border = border)
    }
    xt <- xt + dx.fill
  }
  if (plot && (has.pch || do.lines)) {
    pt.COL <- rep_len(pt.col, n.leg)
    line.COL <- rep_len(line.col, n.leg)
  }
  if (missing(lwd) || is.null(lwd))
    lwd <- par("lwd")
  if (do.lines) {
    if (missing(lty) || is.null(lty))
      lty <- 1
    lty <- rep_len(lty, n.leg)
    lwd <- rep_len(lwd, n.leg)
    ok.l <- !is.na(lty) & (is.character(lty) | lty > 0) &
      !is.na(lwd)
    if (trace)
      catn("  segments2(", xt[ok.l] + x.off * xchar, ",",
           yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
    if (plot)
      segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len *
                  xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l],
                col = line.COL[ok.l])
    xt <- xt + (seg.len + x.off) * xchar
  }
  if (has.pch) {
    pch <- rep_len(pch, n.leg)
    pt.bg <- rep_len(pt.bg, n.leg)
    pt.cex <- rep_len(pt.cex, n.leg)
    pt.lwd <- rep_len(pt.lwd, n.leg)
    ok <- !is.na(pch)
    if (!is.character(pch)) {
      ok <- ok & (pch >= 0 | pch <= -32)
    }
    else {
      ok <- ok & nzchar(pch)
    }
    x1 <- (if (merge && do.lines)
      xt - (seg.len/2) * xchar
      else xt)[ok]
    y1 <- yt[ok]
    if (trace)
      catn("  points2(", x1, ",", y1, ", pch=", pch[ok],
           ", ...)")
    if (plot)
      points2(x1, y1, pch = pch[ok], col = pt.COL[ok], cex = pt.cex[ok],
              bg = pt.bg[ok], lwd = pt.lwd[ok])
  }
  xt <- xt + x.intersp * xchar
  if (plot) {
    if (!is.null(title))
      text2(left + w * title.adj, top - ymax, labels = title,
            adj = c(title.adj, 0), cex = cex, col = title.col)
    text2(xt, yt, labels = legend, adj = adj, cex = cex,
          col = text.col, font = text.font)
  }
  invisible(list(rect = list(w = w, h = h, left = left, top = top),
                 text = list(x = xt, y = yt)))
}




##' Function to plot boxplots in Figure S9 of paper, of the age of cases in Malawi and Angola.

plot.case.boxplots <- function(Cases, country, col = 'cornflowerblue', alpha = 1){
  
  Cases = filter(Cases, Country == country, !is.na(Cases$AgeInyears))
  
  
  cases.plot <- bwplot(AgeInyears ~ factor(YrOnset), data = Cases, ylim = c(-0.5,42.5))
  #cases.plot <- bwplot(AgeInyears ~ factor(YrOnset), data = Cases)
  
  bw.theme <- trellis.par.get()
  bw.theme$box.dot$pch <- "|"
  bw.theme$box.rectangle$col <- "black"
  bw.theme$box.rectangle$lwd <- 4
  bw.theme$box.rectangle$fill <- "grey90"
  bw.theme$box.rectangle$alpha <- 1
  bw.theme$box.umbrella$lty <- 1
  bw.theme$box.umbrella$col <- "black"
  bw.theme$plot.symbol$col <- "grey40"
  bw.theme$plot.symbol$pch <- "*"
  bw.theme$plot.symbol$cex <- 2
  bw.theme$strip.background$col <- "grey80"
  bw.theme$par.main.text$lineheight = 3
  bw.theme$add.text$cex = 5
  bw.theme$par.ylab.text$cex = 2
  bw.theme$par.xlab.text$cex = 2
  bw.theme$par.main.text$cex = 2
  bw.theme$axis.text$cex=2
  l.bw <- update(cases.plot, par.settings = bw.theme, 
                 xlab = "year", ylab = "age", main = country)
  
  l.bw
  #Cases = filter(Cases, AgeInyears<21)
  cases.plot <- boxplot(AgeInyears~factor(YrOnset), data = Cases, 
                        col = alpha('mediumaquamarine',0.75), ylab = "Age", 
                        boxwex = 0.6, xlab = 'year', main = country, pch = 16, cex = 0.2,
                        frame.plot = FALSE, ylim = c(-0.5,20.5))
  
}




closest.path.point.movement.comparison <- function(d, regions, years, use.rep.cases = 1,
                                                   make.inc.cv.scale.same = TRUE, sqrt.inc=F, sqrt.cv=F,
                                                   connect.canonical.path.to.zero = 0, log.incidence = F,
                                                   number.of.additional.points = 5){
  
  ##' filter the data to only include the regions, and years required, along with only
  ##' entries which have incidence greater than 0.
  d = d %>% filter(., Year %in% years, WHO_REGION %in% regions, Incidence > 0,
                   Country != "")
  
  ##' remove years where the coefficient of variation is 0, and the year is before 1995,
  ##' as before this time, no countries had reached elimination, so the only way to 
  ##' have a 0 here, is if there is no data before this time (unless a country reported
  ##' the exact same non-zero number of cases for every year, which didn't happen)
  j = which(d$Year < 1995 & d$Coefficient.of.Variation == 0)
  if(length(j) > 0){
    d = d[-(j), ]
  }
  
  ##' if we want to log the incidence before assigning countries to their location
  ##' on the path, then we do that here. Along with adding a small non-zero value to each of the
  ##' incidences, so that there are no -Inf values created.
  if(log.incidence == T){
    d$Incidence = log(d$Incidence + 0.000001)
  }
  
  if(sqrt.cv == T){
    d$Coefficient.of.Variation = sqrt(d$Coefficient.of.Variation)
  }
  
  ##' if sqrt.inc=T then take sqrt of incidence
  if(sqrt.inc == T){
    d$Incidence <-  sqrt(d$Incidence)
  }
  
  ##' if wantedscale the incidence and coefficient of variation so that they are both in the 0-1 range
  if(make.inc.cv.scale.same == T){
    d$Incidence = scl(d$Incidence)
    d$Coefficient.of.Variation = scl(d$Coefficient.of.Variation)
  }
  
  
  ##' find the mean path for each region by averaging positions in each year
  ##' 
  list[mean.x, mean.y] = find.mean.path(d = filter(d, Year <= 2013), regions)
  
  
  
  if(use.rep.cases ==  1){
    ##' If we are considering the reported cases to locate the countries, then
    ##' Africa and the Americas mean paths converge at roughly 2007 in Africa and 
    ##' 1993 for the Americas. Join the paths at this point, and use it as the canonical
    ##' path 
    k = which(years == 2008)
    k1 = which(years == 1995)
    #zzz = mean.x$AFR[k] + (c(1,2,3)* (mean.x$AMR[k1] - mean.x$AFR[k]))/4
    #zzz1 = mean.y$AFR[k] + (c(1,2,3)* (mean.y$AMR[k1] - mean.y$AFR[k]))/4
    canonical.path = rbind(cbind(mean.x$AFR[1:k], mean.y$AFR[1:k]),# cbind(zzz,zzz1),
                           cbind(mean.x$AMR[k1:nrow(mean.x)], mean.y$AMR[k1:nrow(mean.x)])) %>%
      data.frame()
    
    
  }else{
    ##' If we are considering the state space cases to locate the countries, then
    ##' Africa and the Americas mean paths converge at roughly 2011 in Africa and 
    ##' 1994 for the Americas. We add in 3 intermediate points here and join the paths at this point, 
    ##' and use it as the canonical path
    k = which(years == 2008)
    k1 = which(years == 1995)
    zzz = mean.x$AFR[k] + (c(1,2,3)* (mean.x$AMR[k1] - mean.x$AFR[k]))/4
    zzz1 = mean.y$AFR[k] + (c(1,2,3)* (mean.y$AMR[k1] - mean.y$AFR[k]))/4
    canonical.path = rbind(cbind(mean.x$AFR[1:k], mean.y$AFR[1:k]), cbind(zzz,zzz1),
                           cbind(mean.x$AMR[k1:nrow(mean.x)], mean.y$AMR[k1:nrow(mean.x)])) %>%
      data.frame()
    canonical.path = canonical.path[-c(8,29,31),]
    
  }
  colnames(canonical.path) = c('x', 'y')
  last.canonical = c(tail(canonical.path$x, 1), tail(canonical.path$y, 1))
  
  end.canonical =c(0,0)
  
  if(connect.canonical.path.to.zero == 1){
    l = cbind(seq(last.canonical[1], end.canonical[1], -last.canonical[1] / 3),
              seq(last.canonical[2], end.canonical[2], -last.canonical[2] / 3))
    
    colnames(l) = colnames(canonical.path)
    canonical.path = rbind(canonical.path, l)
  }
  
  ###########################################
  
  ###########################################
  
  ##' here is where the new higher resolution path is constructed
  
  granular.canonical.path = calculate.granular.path(canonical.path, number.of.additional.points)
  
  ###########################################
  
  ###########################################
  
  # canonical.path$x = scl(canonical.path$x)
  # canonical.path$y = scl(canonical.path$y)
  ##' add a column to hold the closest position on the mean trajectory
  d$closest = 0
  
  ##' we only need to calculate the closest position for each country once a year, 
  ##' therefore we restrict to only yearly data
  y1 = unique(round(d$Year))
  d = filter(d, Year %in% y1)
  
  ##' for each entry, then calculate which point is the closest on the trajectory
  for(i in 1 : nrow(d)){
    q = c(d$Coefficient.of.Variation[i], d$Incidence[i]) %>%
      rbind( cbind(granular.canonical.path$x, granular.canonical.path$y)) %>%
      dist() %>% as.matrix()
    d$closest[i] = which(q[1, -(1)] == min(q[1, -(1)])) %>% as.matrix() %>% tail(1)
    
  }
  d$canonical.x = granular.canonical.path$x[d$closest]
  d$canonical.y = granular.canonical.path$y[d$closest]
  return(list(d, granular.canonical.path))
}









##' create a canonical path where the number.of.additional.points are placed between
##' the 38 canonical path points
calculate.granular.path <- function(canonical.path, number.of.additional.points){
  ##' calculate how many total points will be on the granular path
  number.points = (nrow(canonical.path)-1) * (number.of.additional.points + 1) + 1
  
  ##' set up mdata frame to hold this path
  granular.path = data.frame(matrix(NA, number.points, 2))
  colnames(granular.path) = c("x", "y")
  
  ##' intialise the new path at the first positions
  granular.path[1, ] = c(canonical.path$x[1],canonical.path$y[1])
  
  ##' begin a counter
  count = 2
  
  for(i in 2 : nrow(canonical.path)){
    
    ##' x1 and y1 defines the point on the canonical path that we begin from
    x1 = canonical.path$x[i-1]
    y1 = canonical.path$y[i-1]
    
    ##' x2 and y2 defines the point on the canonical path that we move towards
    x2 = canonical.path$x[i]
    y2 = canonical.path$y[i]
    
    ##' calculate the slope we have to move along
    slope = (y2 - y1)/(x2 - x1)
    
    ##' calculate the distance along the path that we have to move
    path.distance = sqrt((y2-y1)^2 + (x2-x1)^2)
    
    ##' split this distance up so that we get the desired number of points between the two points on the path
    step.distance = path.distance / (number.of.additional.points + 1)
    
    ##' calculate the change in x
    x.change = sqrt(step.distance^2/(1+slope^2))
    
    ##' if we move backwards in the x direction, then amend the direction of this movement
    if(canonical.path$x[i] - canonical.path$x[i-1]  < 0){
      x.change = -x.change
    }
    
    ##' if we move up or down in the y direction, then calculate the change in y appropriately
    if(canonical.path$y[i] - canonical.path$y[i-1]>0){
      y.change = abs(slope * x.change)
    } else{
      y.change = -abs(slope * x.change)
    }
    
    ##' add in the additional points iteratively
    for(j in 1 : number.of.additional.points){
      granular.path$x[count] = granular.path$x[count - 1] + x.change
      granular.path$y[count] = granular.path$y[count -1] + y.change
      count = count + 1
    }
    
    ##' add in the point on the path that we were moving towards
    granular.path$x[count] = x2
    granular.path$y[count] = y2
    count = count + 1
  }
  return(granular.path)
}



GetLineOrthogonalToPath <- function(canonical.path.line, smoothed.scaled.pop){
  
  #starting code from https://stackoverflow.com/questions/53167797/determine-transects-perpendicular-to-a-coastline-in-r
  
  AllTransects <- vector('list', 100000) # DB that should contain all transects
  subset_geometry <- canonical.path.line
  
  dx <- c(0, diff(subset_geometry[,'x'])) # Calculate difference at each cell comapred to next cell
  dy <- c(0, diff(subset_geometry[,'y']))
  
  dseg <- sqrt(dx^2+dy^2)                 # get rid of negatives and transfer to uniform distance per segment (pythagoras)
  dtotal <- cumsum(dseg)                  # cumulative sum total distance of segments
  
  linelength = sum(dseg)                  # total linelength
  sep <- 0.001
  start <- 0
  pos = seq(start,linelength, by=sep)     # Array with postions numbers in meters
  whichseg = unlist(lapply(pos, function(x){sum(dtotal<=x)})) # Segments corresponding to distance
  
  pos=data.frame(pos=pos,                            # keep only 
                 whichseg=whichseg,                  # Position in meters on line
                 x0=subset_geometry[whichseg,1],     # x-coordinate on line
                 y0=subset_geometry[whichseg,2],     # y-coordinate on line
                 dseg = dseg[whichseg+1],            # segment length selected (sum of all dseg in that segment)
                 dtotal = dtotal[whichseg],          # Accumulated length
                 x1=subset_geometry[whichseg+1,1],   # Get X coordinate on line for next point
                 y1=subset_geometry[whichseg+1,2]    # Get Y coordinate on line for next point
  )
  
  pos$further =  pos$pos - pos$dtotal       # which is the next position (in meters)
  pos$f = pos$further/pos$dseg              # fraction next segment of its distance
  pos$x = pos$x0 + pos$f * (pos$x1-pos$x0)  # X Position of point on line which is x meters away from x0
  pos$y = pos$y0 + pos$f * (pos$y1-pos$y0)  # Y Position of point on line which is x meters away from y0
  
  pos$theta = atan2(pos$y0-pos$y1,pos$x0-pos$x1)  # Angle between points on the line in radians
  pos$object = i
  
  ##xxamy - add column to pos with thickness desire per whichseg
  pos$thickness <- smoothed.scaled.pop[whichseg]
  
  ###### Define transects
  pos$thetaT = pos$theta+pi/2         # Get the angle
  dx_poi <- pos$thickness*cos(pos$thetaT) # coordinates of point of interest as defined by position length (sep)
  dy_poi <- pos$thickness*sin(pos$thetaT) 
  
  # transect is defined by x0,y0 and x1,y1 with x,y the coordinate on the line
  output <-     data.frame(pos = pos$pos,
                           x0 = pos$x + dx_poi,       # X coordinate away from line
                           y0 = pos$y + dy_poi,       # Y coordinate away from line
                           x1 = pos$x - dx_poi,       # X coordinate away from line
                           y1 = pos$y - dy_poi,       # X coordinate away from line
                           theta = pos$thetaT,    # angle
                           x = pos$x,             # Line coordinate X
                           y = pos$y,             # Line coordinate Y
                           object = pos$object,
                           nextx = pos$x1,
                           nexty = pos$y1) 
  
  
  pos.opts <- seq(start,linelength, by=sep)   
  #build new data frame
  row.indexes <- rep(NA, nrow(canonical.path.line))
  for (i in 1:nrow(canonical.path.line)){
    row.indexes[i] <- which(abs(pos.opts-dtotal[i])==min(abs(pos.opts-dtotal[i])))
  }
  df.polygon <- output[row.indexes,c("x0","y0","x1","y1","x","y")]  
  
  return(df.polygon)
}



Unscaled.Canonical.Path <- function(canonical.path.data, canonical.path, regions, years){
  
  x <- canonical.path$x
  y <- canonical.path$y
  
  #setting up the data the same as when the canonical path was created
  d <- canonical.path.data
  d = d %>% filter(., Year %in% years, WHO_REGION %in% regions, Incidence > 0,
                   Country != "")
  j = which(d$Year < 1995 & d$Coefficient.of.Variation == 0)
  if(length(j) > 0){
    d = d[-(j), ]
  }
  
  #back transforming incidence and cv
  log.d.inc <- log(d$Incidence + 0.000001)
  new.x <- x*(max(d$Coefficient.of.Variation)-min(d$Coefficient.of.Variation))+(min(d$Coefficient.of.Variation))
  new.y <- exp((y*(max(log.d.inc)-min(log.d.inc)))+(min(log.d.inc)))-0.000001
  
  #plot(d$Coefficient.of.Variation,d$Incidence)
  #points(new.x, new.y, col=2)
  #plot(new.x, new.y, col=2)
  #points(new.x, sqrt(new.y), col=4)
  #axis(4, at=sqrt(c(0.1,0.5,1,1.5)), labels=c(0.1,0.5, 1, 1.5), col="blue", col.ticks="blue", col.lab="blue")
  
  canonical.path2 <- data.frame(x=new.x, y=new.y)
  return(canonical.path2)
  
}



post.mcv2.movements <- function(df.mcv2, d1b){
  for(i in 1 : nrow(df.mcv2)){
    cc = df.mcv2$Country[i]
    post.mcv2.years = max(min(d1b$Year),df.mcv2$MCV2_Intro_Year[i]) : max(d1b$Year)
    post.mcv2.country = d1b %>% filter(Country == cc, Year %in% post.mcv2.years)
    time.since.mcv2 =  0:(length(post.mcv2.years)-1)
    pos.delta = post.mcv2.country$closest - post.mcv2.country$closest[1]
    post.mcv2.country = cbind(post.mcv2.country, time.since.mcv2=time.since.mcv2, pos.delta=pos.delta)
    if(i == 1){
      post.mcv2.positions = post.mcv2.country
    } else{
      post.mcv2.positions = rbind(post.mcv2.positions, post.mcv2.country)
    }
  }
  return(post.mcv2.positions)
}



average.movement.since.mcv2 <- function(d, lower_q = 0.025, upper_q = 0.975){
  average.movement = data.frame(matrix(0, max(d$time.since.mcv2), 5))
  colnames(average.movement) = c("time.since.mcv2", "mean_postion", "median_movement", "lower_quantile","upper_quantile")
  for(i in 1 : max(d$time.since.mcv2)){
    dd = d %>% filter(time.since.mcv2 == i)
    average.movement[i, ] = c(i, mean(dd$closest), median(dd$pos.delta), quantile(dd$pos.delta, lower_q), quantile(dd$pos.delta, upper_q))
  }
  return(average.movement)
}

