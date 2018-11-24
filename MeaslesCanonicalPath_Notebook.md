---
title: "Measles and the Canonical Path to Elimination"
output:
  html_document:
    keep_md: yes
  pdf_document: default
bibliography: bibliography.bib
---



This notebook will construct the figures for the paper by Graham et al., titled "Measles and the Canonical Path to Elimination." First we weighted case data along with population data to construct a dataset which can be used to plot each country through time in 'incidence-space', which has the incidence on the y-axis, and the coefficient of variation of incidence over time on the x-axis.





Following the construction of this data set, we plot Figure 1A from the paper. This figure shows the position of
countries in the Americas and in Africa in two years, along with the mean path taken by these continents from 
the beginning of the data set to the end of it. These mean paths are shown by the green and purple lines in the
figure below, with the green line being for Africa and the purple one for the Americas.




![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

By combining the trajectories, we can create the canonical path towards elimination seen in Fig. 1B in the paper.
This figure was created in Adobe Illustrator, so cannot be reproduced here.

The canonical path seen in Fig. 2A is constructed by putting the incidence and coefficient of variation on the same scale (by multiplying the incidence by max(coefficient of variation) / max(incidence)) and then calculating the mean trajectories of Africa and the Americas over time, and combining these at the point that they intersect, i.e. at the point that the green and purple lines cross in the plot above.




![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

The analysis performed in the paper relies on the fact that when we calculate the position of a country at a given time in the incidence-space, we can calculate which point on the canonical path this country's location is closest to. However, we can see that when we do this there are positions which are close to each other distance wise on this path, but far in terms of progression towards elimination, due to the fact that on the lower end of incidence, there is very little distance between the points. For example if, in a given year, a country lay at the point 1.4 on the x-axis, and has a fairly low incidence, then it could easily be assigned to a point very close to the end of the path to elimination or one which is about half way along the path. This is not a desirable property of the canoncial path. To help distinguish the points more clearly, and help this assignment of nations, we do two things. Firstly, we take the log (natural base) of the incidence, and secondly, we transform both the incidence and the coefficient of variation data to be on the 0-1 scale. This is done in the chunk below.



When we plot the canonical path now, we see that there is a much greater distinction between points at the low incidence part of the path (Fig. S3 in the supplement).  

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

We can now assign countries to the closest point on this path each year that we have data for, using the re-scaled canonical path. This is seen in Fig. 4A in the paper, where we plot all countries position on the path overtime, in addition to population-weighted region averages.  



![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Countries do not necessarily progress smoothly along the path as if, for example, after years of low incidence is interrupted by a year of high incidence, then this will increase the x and y position of the country in incidence-space and hence the position on the canonical path will head backwards. We demonstrate these movements along the path by plotting the % change in path postion between 1990 and 2017 for each country in addition to population-weighted region averages, as see in Fig. 4B in the paper.

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

Given our assignement of all countries to the closest point on the path, for any given year a map can be produced which is colored according to each countries position along the path. This is seen below for 2017 (Fig. 2C of the paper).  We have also made available web application, developed with R package Shiny, to view all countries position on the canonical path between 1980 and 2014 (available at http://iddynamics.jhsph.edu/apps/shiny/measlescanonicalpath/).



![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

We can also capture direction of movement of each country within the incidence space in order to assess "deviations" from the path or trajectories of countries in incidence space that were different than expected given the characterized canonical path.  Fig. 4C in the paper is reproduced here.  



We can plot the expected movement, observed movement, and difference for each region over 10 or 8 year increments.  The heat map captures the frequency of countries movement in a direction grouped by 10 degree angle groups. This plot will look slightly different than the main text because the main text range is from 0 to >20%, where we grouped frequencies between 20 and 25% into one color dark blue color.  Below we show the range 0 to 25%. 



We can plot all regions and years here as shown in Fig. 4C.


```r
all.plots
```

```
## $reg1.year1
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

```
## 
## $reg1.year2
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-2.png)<!-- -->

```
## 
## $reg1.year3
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-3.png)<!-- -->

```
## 
## $reg2.year1
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-4.png)<!-- -->

```
## 
## $reg2.year2
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-5.png)<!-- -->

```
## 
## $reg2.year3
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-6.png)<!-- -->

```
## 
## $reg3.year1
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-7.png)<!-- -->

```
## 
## $reg3.year2
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-8.png)<!-- -->

```
## 
## $reg3.year3
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-9.png)<!-- -->

```
## 
## $reg4.year1
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-10.png)<!-- -->

```
## 
## $reg4.year2
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-11.png)<!-- -->

```
## 
## $reg4.year3
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-12.png)<!-- -->

```
## 
## $reg5.year1
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-13.png)<!-- -->

```
## 
## $reg5.year2
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-14.png)<!-- -->

```
## 
## $reg5.year3
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-15.png)<!-- -->

```
## 
## $reg6.year1
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-16.png)<!-- -->

```
## 
## $reg6.year2
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-17.png)<!-- -->

```
## 
## $reg6.year3
```

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-16-18.png)<!-- -->

We can also plot the expected and observed movements for all countries and years, as shown as the top of Fig. 4C. 

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

Along with this analysis, we used a method for estimating the susceptibility in each year of age which was developed by @tak2015. Each age cohort's susceptiblity is estimated based on its opportunity for immunization (via routine and supplmentary activities) and risk of natural infection. The probabliity of immunization was estimated per WHO reported administrative vaccination coverage estimates (@who). The probability of natural infection by age was estimated by assuming a constant hazard of infection over age that scaled each year relative to the proportional decline in estimated measles incidence corrected for under-reporting. Measles incidence, corrected for under-reporting, was estimated using a state space model per @statespace and @simons2012. The following chunks of code first uses the state-space model to estimate the number of measles cases (corrected for under-reporting) by year and country, then goes on to infer susceptibilty by age for each country and each year.





For ease, we can simply read in the already estimated proportion susceptible by age, year, and country.



We can then use the dataset of the number of susceptibles by country, year, and age, to link each country and year to the canonical path and determine the estimated total proporiton of susceptibles, and age-specific proportion susceptibles at each canonical path point. 



Here we plot the proportion of susceptible individuals by canonical path point (Fig. 2C in the paper). The horizontal dashed lines display the critical level of immunity if the basic reproduction number of measles is 15 or 20 (critical level = 1-(1 / basic reproductive number)).

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

We plot the age distribution of susceptibles by canonical path point 

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

We also have data on the mean age of measles cases in multiple countries from 2000-2016. This data has increasing numbers of data as time goes on.



Again, these data can be linked to canonical path point. Here we make a boxplot of the esimated distribution of ages against the mean age of measles cases, to see how they compare (Fig. 2F). The estimates are colored the same colors as the canonical path points, and the case data is transparent yellow and green. 

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-25-1.png)<!-- -->


The figures below show the estimated proportion of susceptibles who are under 5 (Fig. 2D), and the age at which an SIA would have to go up to in order to cover 90% of all susceptibles (Fig. 2E).




![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-27-1.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-27-2.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-27-3.png)<!-- -->


# Supplementary 

Here we reproduce the figures seen in the supplement of the paper.  The following codes reproduces Fig. S1. 

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-29-1.png)<!-- -->


![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

The following codes reproduces Fig. S2. 

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-31-1.png)<!-- -->




![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-33-1.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-33-2.png)<!-- -->


Figure S4 of the supplement has the location of countries in the WHO Africa and Americas Regions in incidence-space in 1990 and 2014 post scaling of incidence and CV. This is re-created here.




![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-35-1.png)<!-- -->

Figure S5 shows the trajectory of the Americas and Africa when we take the median of these regions paths rather than the mean.





![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-37-1.png)<!-- -->


We used an established discrete time age-structured mathematical model, introduced in @metcalf2012a and @metcalf2012b to simulate measles transmission dynamics for each country in the WHO Americas and Africa Regions.  We used the same gaussian weights as in the empirical analysis to create an incidence-space for each country over time, and compared this to estimates of measles incidence per @statespace and @simons2012.


```
## -------------------------------------------------------------------------
```

```
## You have loaded plyr after dplyr - this is likely to cause problems.
## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
## library(plyr); library(dplyr)
```

```
## -------------------------------------------------------------------------
```

```
## 
## Attaching package: 'plyr'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     arrange, count, desc, failwith, id, mutate, rename, summarise,
##     summarize
```

```
## The following object is masked from 'package:purrr':
## 
##     compact
```

We plotted AFRO countries in incidence-space (Fig. S6-S11).

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-1.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-2.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-3.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-4.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-5.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-6.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-7.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-8.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-9.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-10.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-11.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-12.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-13.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-14.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-15.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-16.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-17.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-18.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-19.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-20.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-21.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-22.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-23.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-24.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-25.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-26.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-27.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-28.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-29.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-30.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-31.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-32.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-33.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-34.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-35.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-36.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-37.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-38.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-39.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-40.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-41.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-42.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-43.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-44.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-45.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-46.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-47.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-39-48.png)<!-- -->

We plotted AMRO countries in incidence-space (Fig. S12-S15).

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-1.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-2.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-3.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-4.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-5.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-6.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-7.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-8.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-9.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-10.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-11.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-12.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-13.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-14.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-15.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-16.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-17.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-18.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-19.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-20.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-21.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-22.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-23.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-24.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-25.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-26.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-27.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-28.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-29.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-30.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-31.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-32.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-40-33.png)<!-- -->

The estimated case data per @statespace and @simons2012 can also be used to plot a figure similar to figure 1 of the paper, demonstrating the path of the Americas and Africa through incidence space. This is seen below and is Fig. S16. The only difference is that Fig. 1A uses reported cases, and Fig. S16 uses cases corrected for under-reporting.






![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-43-1.png)<!-- -->

Figure S17 demonstrates movement along the canonical path for 10 countries of interest from the Americas, Africa, and European.  These 10 countries are highlighted from Fig. 4A-B in the main text to conceptualize an individual country's movement on the path and discuss data biases or true dynamics that may have contributed to these movements.  

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-44-1.png)<!-- -->![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-44-2.png)<!-- -->

Figures S18 and S19, contains the results of modeling incidence and coefficient of variation as dependent variables separately, with birth rate, vaccination proportion and (birth rate) times (1-vacciantion proportion) all used independent variables (one at a time) using generalized additivie models (GAM). The term (birth rate) times (1-vaccination proportion) is an approximation of how quickly individuals who are susceptible to measles are recruited into the population, therefore we term this the rate of susceptible recruitment. To produce the GAM plots, we first need to construct a data set which has the estimated mean age of susceptibles by country and by year.



We can then plot these predicted paths in the incidence-space over values of birth rate, vaccination rate, and the rate of susceptible recruitment with reported cases (Fig. S18).

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-46-1.png)<!-- -->
   
    
![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-47-1.png)<!-- -->

Figure S18C also displays the association between incidence or coefficient of variation with the birth rate and vaccination coverage in Africa.

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-48-1.png)<!-- -->

Similarly, we can use model incidence and coefficient of variation using GAMS, as seen above for the reported case data. This is seen below and Fig. S19 of the paper.





![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-50-1.png)<!-- -->
   
    
![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-51-1.png)<!-- -->

Figure S20 shows the age distribution of cases in Malawi and Angola from 2006-2013. This figure is reproduced below.

![](MeaslesCanonicalPath_Notebook_files/figure-html/unnamed-chunk-52-1.png)<!-- -->


# References










