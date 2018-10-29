datainput=function(allbirths.w,allcases.w,alldeaths.w,allpop.w,allsia.w,allmcv1.w,mcv1.12mo.w,allmcv2.w,outbreaks.w){              # number in country list
births<-as.numeric(allbirths.w)   # transform birth rate into numerical form
cases<-as.numeric(allcases.w)     # transform reported rate into numerical form
death<-as.numeric(alldeaths.w) # death data
births<-births-death                # adjusted birth
pop<-as.numeric(allpop.w)         # transform population into numerical form
# sia= supplementary immunization activities
sia<-as.numeric(allsia.w) 		# transform SIA into numerical form
mcv1<-as.numeric(allmcv1.w)       # transform 1st pulse campaign into numerical form
mcv2<-as.numeric(allmcv2.w)       # transform 2nd pulse campaign into numerical form

ifelse(mcv1.12mo.w==0,
X<-as.numeric(births - 0.84*(mcv1)*(1-mcv2)*births - 0.99*(mcv1)*mcv2*births - 0.925*(1-mcv1)*(mcv2)*births),
X<-as.numeric(births - 0.925*(mcv1)*(1-mcv2)*births - 0.99*(mcv1)*mcv2*births - 0.925*(1-mcv1)*(mcv2)*births))
# adjusted birth rate
vecs<-list(X=X,Ir=cases,sia=sia,mcv1=mcv1,mcv2=mcv2,pop=pop,births=births,outbreaks=outbreaks.w)
# sort all necessary data into a list
returns=vecs
# return the list
returns
# function return
}
