# a4a assessment of the Iberian stock to mimic SS3 assessment

library(FLa4a)

load("ib.Rdata")

SS3.fit<-IB.stk # stock file with SS3 output
  
# set a4a to mimic SS3

fmodel <- ~s(replace(age,age %in% c(3:5), 3), k = 4, by = breakpts(year, 1991)) + s(year, k = 20)
qmod <- list(~s(replace(age, age %in% c(2:5), 2),k = 3),~1)      
          
a4a.fit <- sca(IB.stk, IB.idx,fit="assessment",fmodel=fmodel,qmodel=qmod) # final model
            
             




