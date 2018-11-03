##########################################################################################
# Check that the order of the rows in the EVPI table are correct
##########################################################################################

library(shiny)

shinyServer(function(input,output) {
isovalue<-reactive({input$country})

output$Est_Plot <- renderPlot({
iso<-isovalue()
ind<-which(rownames(allcases)==iso)
par(mar=c(2,4,2,4))
par(mfrow=c(2,1))
plot(1980:(1980+dim(allcases)[2]-1),unlist(allcases[ind,]),xlab="year",ylab="observed cases",type="l",lwd=2,col=4)
par(new=T)
plot(1980:(1980+dim(allcases)[2]-1),c(NA,unlist(CASES[ind,])),xlab=NA,ylab=NA,type="l",lwd=2,col=2,axes=F)
mtext("estimated cases",4,line=2.5)
axis(4)
par(new=T)
plot(1978:(1980+dim(allcases)[2]-1),c(1,NA,is.na(unlist(outbreaks[ind,]))),xlab=NA,ylab=NA,type="p",lwd=2,col=1,axes=F,xlim=c(1980,(1980+dim(allcases)[2]-1)))
legend(1980,.8,legend=c("observed cases","estimated cases","outbreaks"),bty="n",lty=c(1,1,NA),col=c(4,2,1),pch=c(NA,NA,1),cex=.75)
par(new=F)

plot(1980:(1980+dim(allcases)[2]-1),unlist(allmcv1[ind,]),xlab="year",ylab="MCV1",type="l",lwd=2,col=4)
par(new=T)
plot(1980:(1980+dim(allcases)[2]-1),c(NA,unlist(measles_deaths[ind,])),xlab=NA,ylab=NA,type="l",lwd=2,col=2,axes=F)
mtext("estimated deaths",4,line=2.5)
axis(4)
par(new=T)
plot(1980:(1980+dim(allcases)[2]-1),unlist(allsia[ind,]),xlab=NA,ylab=NA,type="p",lwd=2,col=1,axes=F)
legend(1980,.8*max(unlist(allsia[ind,])),legend=c("MCV1","estimated deaths","SIA coverage"),lty=c(1,1,NA),col=c(4,2,1),pch=c(NA,NA,1),cex=.75,bty="n")
par(new=F)
})

output$transmission_Plot <- renderPlot({
iso<-isovalue()
ind<-which(rownames(allcases)==iso)
hist(PAR[,1],breaks=seq(min(PAR[,1],na.rm=T),max(PAR[,1],na.rm=T),length=100),xlab="transmission parameter",main="transmission parameter")
abline(v=PAR[ind,1],col=2,lwd=2)
})

output$reporting1_Plot <- renderPlot({
iso<-isovalue()
ind<-which(rownames(allcases)==input$country)
hist(PAR[,3],breaks=seq(min(PAR[,3],na.rm=T),max(PAR[,3],na.rm=T),length=100),xlab="reporting probability",main="reporting probability")
abline(v=PAR[ind,3],col=2,lwd=2)
})

output$reporting2_Plot <- renderPlot({
iso<-isovalue()
ind<-which(rownames(allcases)==input$country)
hist(PAR[,4],breaks=seq(min(PAR[,4],na.rm=T),max(PAR[,4],na.rm=T),length=100),xlab="reporting probability",main="reporting probability")
abline(v=PAR[ind,4],col=2,lwd=2)
})


})