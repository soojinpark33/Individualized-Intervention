#::::: R code for the application to HSLS:09 in 
#::::: Causal Decomposition Analysis with Time-Varying Mediators:
#::::: Designing Individualized Interventions to Reduce Social Disparities
#::::: This R code is written by Soojin Park and Namhwa Lee.

### Required Package ###
library(geepack); library(DynTxRegime) ; library(boot); library(rpart) 
library(rpart.plot) ; library(parallel); library(dplyr); library(tidyverse)

### Load the data and pre-processing ###
load("~/GSR_project2/data_imputed.RData")
codebook <- read.csv("~/GSR_project2/Codebook.csv")
out.name <- codebook[4:5, 1] # outcome
med.name <- codebook[6:7, 1] # mediator: algebra1 & algebra2
basecov.name <- codebook[9, 1] # female
conf.name1 <- codebook[c(10:29), 1] # S1MUNDERST / X1MTHID / X3THIMATH
conf.name2 <- codebook[3,1] # X1TXMTH: Mathematics score in 9th grade
DATA <- data_impute
DATA$S1M8GRADE <- relevel(DATA$S1M8GRADE, ref = "b") # Final grade in 9th grader's most advanced 8th grade math course
DATA1 <- subset(DATA, DATA$race %in% c("white", "black"))
DATA1$race <- factor(DATA1$race)
DATA1$seqn <- rep(1:nrow(DATA1))
DATA1$GRADE <-ifelse(DATA1$S1M8GRADE=="Below D",1,0)
DATA1$female <- scale(as.numeric(DATA1$female), center=TRUE) # center baseline covariates
DATA1$advanced <- as.factor(as.numeric(as.numeric(DATA1$X3THIMATH)>8))


### Estimate Optimal DTRs through weighting method ###
# Second-Stage Analysis using IPW #
moPropen2 <- buildModelObj(model = ~ race+female+S1M8GRADE+ses
                           +X1PAR1OCC_STEM1+X1PAR2OCC_STEM1+P1EDUEXPECT+P1EDUASPIRE+M1SEX
                           +M1INTEREST+locale+X1SCHOOLCLI+X2PROBLEM+M1UNPREPPCT+requirements
                           +schextracur+friend+X1MTHINT+X1MTHUTI+X1SCHOOLENG+X1SCHOOLBEL+X1TXMTH+algebra1,
                           solver.method = 'glm',
                           solver.args = list('family'='binomial'),
                           predict.method = 'predict.glm',
                           predict.args = list(type='response'))

moClass2 <- buildModelObj(model = ~ X1TXMTH +X1MTHEFF+X1MTHINT,
                          solver.method = 'rpart',
                          solver.args = list(method="class"),
                          predict.args = list(type='class'))

fitSS_IPW <- optimalClass(moPropen = moPropen2,
                          moClass = moClass2,
                          data = DATA1, response = DATA1$X2TXMTH,
                          txName = 'advanced', verbose=FALSE)

### First-Stage Analysis using IPW #
moPropen1 <- buildModelObj(model = ~ race + female+S1M8GRADE+ses+X1PAR1OCC_STEM1
                           +X1PAR2OCC_STEM1+P1EDUEXPECT+P1EDUASPIRE+M1SEX+M1INTEREST
                           +locale+X1SCHOOLCLI+X2PROBLEM+M1UNPREPPCT+requirements+schextracur
                           +friend+X1MTHINT+X1MTHUTI+X1MTHEFF+X1SCHOOLENG+X1SCHOOLBEL,
                           solver.method = 'glm',
                           solver.args = list('family'='binomial'),
                           predict.method = 'predict.glm',
                           predict.args = list(type='response'))

moClass1 <- buildModelObj(model = ~  S1M8GRADE+X1MTHEFF+X1MTHINT,
                          solver.method = 'rpart',
                          solver.args = list(method="class"),
                          predict.args = list(type='class'))

fitFS_IPW <- optimalClass(moPropen = moPropen1,
                          moClass = moClass1,
                          data = DATA1, response = fitSS_IPW,
                          txName = 'algebra1', verbose = FALSE)

# Define optimal DTR #
optFS <- ifelse(is.na(optTx(fitFS_IPW)$optimalTx),1,optTx(fitFS_IPW)$optimalTx)
DATA1$opt.alg1 <- as.factor(optFS)
DATA1$opt.adv <- as.factor(optTx(fitSS_IPW)$optimalTx)

### Implement a Bootstrap method for obtaining proposed estimates ###
boot_ipw <- function(data, indices){
  tryCatch(
    {
      DATA1 <- data[indices,]
      #:::::::::::::::::::::::::::::::::::::::#
      #:::::::::: Initial disparity ::::::::::#
      #:::::::::::::::::::::::::::::::::::::::#
      fit.lm <- lm(X2TXMTH ~ race + female, data = DATA1)
      tau <- coef(fit.lm)[2]
      
      #:::::::::::::::::::::::::::::::::::::::#
      #:::::: Controlled Direct Effects ::::::#
      #:::::::::::::::::::::::::::::::::::::::#
      ### construct weight for med2 ###
      fit.m2 <- glm(advanced ~ race+female+S1M8GRADE+ses+X1PAR1OCC_STEM1+X1PAR2OCC_STEM1+P1EDUEXPECT+P1EDUASPIRE
                    +M1SEX+M1INTEREST+locale+X1SCHOOLCLI+X2PROBLEM+M1UNPREPPCT+requirements+schextracur+friend
                    +X1MTHINT+X1MTHUTI+X1MTHEFF+X1SCHOOLENG+X1SCHOOLBEL+X1TXMTH+algebra1,
                    family = binomial(logit),
                    data = DATA1)
      
      ### Calculate weights ###
      p.med2 <- ifelse(DATA1$advanced == "0",
                       1 - predict(fit.m2, type = "response"),
                       predict(fit.m2, type = "response"))
      w2 <- 1 / p.med2
      
      ### construct weight for med1 ###
      fit.m1 <- glm(algebra1 ~ race + female + S1M8GRADE + ses + X1PAR1OCC_STEM1 + X1PAR2OCC_STEM1
                    + P1EDUEXPECT + P1EDUASPIRE + M1SEX + M1INTEREST + locale + X1SCHOOLCLI
                    + X2PROBLEM + M1UNPREPPCT + requirements + schextracur + friend
                    + X1MTHINT + X1MTHUTI + X1MTHEFF + X1SCHOOLENG + X1SCHOOLBEL,
                    family = binomial(), data = DATA1)
      
      ### Calculate weights ###
      p.med1 <- ifelse(DATA1$algebra1 == "0",
                       1 - predict(fit.m1, type = "response"),
                       predict(fit.m1, type = "response"))
      w1 <- 1 / p.med1
      
      ### Truncate weights ###
      w <- w1*w2
      w.threshold <- quantile(w, probs=c(0.01, 0.99)) # 1.04, 12.41
      w.trunc <- ifelse(w < w.threshold[1], as.numeric(w.threshold[1]), w)
      w.trunc <- ifelse(w.trunc > w.threshold[2], as.numeric(w.threshold[2]), w.trunc)
      DATA1$w.trunc <- w.trunc
      
      ### Calculate CDE ###
      fit.lm <- lm(X2TXMTH ~ race + female + algebra1 + advanced + race*advanced,
                   weights = w.trunc , data=DATA1)
      
      zeta_cde10 <- coef(fit.lm)[2]
      zeta_cde11 <- sum(coef(fit.lm)[c(2,6)])

      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::: Interventional Marginal Effects (VanderWeele):::::#
      #::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
      ### construct a mediator model ###
      fit.m1 <- glm(algebra1 ~ race + female, data = DATA1, family = binomial(logit))
      fit.m2 <- glm(advanced ~ race + female, data = DATA1, family = binomial(logit))
      
      ### Enforcing Race ###
      DATA_R1 <- DATA_R2 <- DATA1
      DATA_R1[,"race"]<-"white" #white
      DATA_R2[,"race"]<-"black" #black
      
      ### P(M1 = 1 | R=0, c) ###
      p_w1 <- predict(fit.m1, newdata=DATA_R1,type = "response")
      ### P(M1 = 1 | R=1, c) ###
      p_b1 <- predict(fit.m1, newdata=DATA_R2,type = "response")
      ### Marginalize ###
      p1 <- mean(p_b1-p_w1)
      
      ### P(M2 = 1 | R=0, c) ###
      p_w2 <- predict(fit.m2, newdata=DATA_R1,type = "response")
      ### P(M2 = 1 | R=1, c) ###
      p_b2 <- predict(fit.m2, newdata=DATA_R2,type = "response")
      ### Marginalize ###
      p2 <- mean(p_b2-p_w2)
      
      ### Calculate IME ###
      delta_ime <- coef(fit.lm)[4]*p1 + sum(coef(fit.lm)[c(5,6)])*p2
      zeta_ime <- tau-delta_ime
      
      #:::::::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::: Individualized Conditional Direct Effects :::::#
      #:::::::::::::::::::::::::::::::::::::::::::::::::::::#
      ##### Regression Method #####
      ### indicator for compliance sample ###
      DATA1$ind2 <- as.numeric(DATA1$advanced==DATA1$opt.adv)
      DATA1$ind1 <- as.numeric(DATA1$algebra1==DATA1$opt.alg1)
      # denominator in W_ICDE is the same with the weight of CDE.
      fit.lm <- lm(X2TXMTH ~ race + female, data=DATA1, weights = (ind1 & ind2)*w.trunc)
      zeta_icde <- coef(fit.lm)[2]
      
      ##### Weighting Method #####
      # Index: female = 0
      DATA1 <- DATA1 %>% 
        mutate(C0.idx = as.numeric(female == unique(female)[2,]))
    
      E_W.icdeY_R1C0 <- DATA1 %>% 
        filter(race == "black", C0.idx==1) %>% 
        summarise(weighted.mean(X2TXMTH, ind1*ind2*w.trunc))
      
      E_W.icdeY_R0C0 <- DATA1 %>% 
        filter(race == "white", C0.idx==1) %>% 
        summarise(weighted.mean(X2TXMTH, ind1*ind2*w.trunc))
      
      zeta_icde.wm <- as.numeric(E_W.icdeY_R1C0 - E_W.icdeY_R0C0)
      
      #:::::::::::::::::::::::::::::::::::::::::::::::::
      #::::: Individualized Interventional Effects :::::
      #:::::::::::::::::::::::::::::::::::::::::::::::::
      DATA_white <- subset(DATA1, race == "white") 
      DATA_black <- subset(DATA1, race == "black")
      
      # E [Y | R=1, C=c] #
      fit.lm1 <- lm(X2TXMTH ~ female, data=DATA_black)
      W_Y1 <- coef(fit.lm1)[1]
      # E [Y | R=0, C=c] #
      fit.lm0 <- lm(X2TXMTH ~ female, data=DATA_white)
      W_Y0 <- coef(fit.lm0)[1]
      
      ##### Weighting Method #####
      ### Step 1 ###
      # Among white (R=0) #
      fit.Im1_0c <- glm(ind1 ~ female, data = DATA_white, family = binomial(logit))
      fit.Im2_0c <- glm(ind2 ~ female, data = DATA_white, family = binomial(logit))
      # P [ I(M1=d1.opt)=1 | R=0, C=c ] # 
      pi_Im11_0c <- plogis(coef(fit.Im1_0c)[1])
      # P [ I(M1=d1.opt)=0 | R=0, C=c ] # 
      pi_Im10_0c <- 1 - pi_Im11_0c
      # P [ I(M2=d2.opt)=1 | R=0, C=c ] # 
      pi_Im21_0c <- plogis(coef(fit.Im2_0c)[1])
      # P [ I(M2=d2.opt)=0 | R=0, C=c ] # 
      pi_Im20_0c <- 1 - pi_Im21_0c
      
      ### Step 2 ###
      ind1 <- DATA1$ind1 ; ind2 <- DATA1$ind2
      W_IIE.11 <- (ind1 & ind2) * w.trunc
      # theta1 = 1, theta2 = 0 : M1=d1.opt & M2!=d2.opt #
      W_IIE.10 <- (ind1 & !ind2) * w.trunc
      # theta1 = 0, theta2 = 1 : M1!=d1.opt & M2=d2.opt #
      W_IIE.01 <- (!ind1 & ind2) * w.trunc
      # theta1 = 0, theta2 = 0 : M1!=d1.opt & M2!=d2.opt #
      W_IIE.00 <- (!ind1 & !ind2) * w.trunc
      
      ### Step 3 ###
      # (theta1, theta1) = (1,1) #
      idx_R1 <- which(DATA1$race == "black")
      fit.lm11 <- lm(X2TXMTH ~ female, data=DATA_black, weights = W_IIE.11[idx_R1])
      term11 <- pi_Im11_0c * pi_Im21_0c * coef(fit.lm11)[1]
      
      # (theta1, theta1) = (1,0) #
      fit.lm10 <- lm(X2TXMTH ~ female, data=DATA_black, weights = W_IIE.10[idx_R1])
      term10 <- pi_Im11_0c * pi_Im20_0c * coef(fit.lm10)[1]
      
      # (theta1, theta1) = (0,1) #
      fit.lm01 <- lm(X2TXMTH ~ female, data=DATA_black, weights = W_IIE.01[idx_R1])
      term01 <- pi_Im10_0c * pi_Im21_0c * coef(fit.lm01)[1]
      
      # (theta1, theta1) = (0,0) #
      fit.lm00 <- lm(X2TXMTH ~ female, data=DATA_black, weights = W_IIE.00[idx_R1])
      term00 <- pi_Im10_0c * pi_Im20_0c * coef(fit.lm00)[1]
      
      W_iie.Y <- term11 + term10 + term01 + term00
      delta_iie1.wm <- W_Y1 - W_iie.Y 
      zeta_iie0.wm <- W_iie.Y - W_Y0
      
      ##### Regression Method #####
      fit.lm <- lm(X2TXMTH ~ race + female + ind1 + ind2, data=DATA1, weights=w.trunc)
      fit.ind1 <- glm(ind1 ~ race + female, data=DATA1, family=binomial(logit))
      fit.ind2 <- glm(ind2 ~ race + female, data=DATA1, family=binomial(logit))
      
      alpha1 <- plogis(sum(coef(fit.ind1)[1:2])) - plogis(coef(fit.ind1)[1])
      alpha2 <- plogis(sum(coef(fit.ind2)[1:2])) - plogis(coef(fit.ind2)[1])
      
      delta_iie1 <- as.numeric(alpha1*coef(fit.lm)["ind1"]+alpha2*coef(fit.lm)["ind2"])
      zeta_iie0 <- as.numeric((W_Y1 - W_Y0) - delta_iie1)
      
      #:::::::::::::::::::::::::::::::::::::::::::::::::#
      #::::::: Individualized Conditional Effects ::::::#
      #:::::::::::::::::::::::::::::::::::::::::::::::::#
      ### Step 1 ###
      ### Numerator ###
      ## Fit a mediator model : advanced ##
      step1.fit2 <- glm(advanced ~ female + X1MTHEFF + X1MTHINT + X1TXMTH,
                        family = binomial(logit), data = DATA_white)
      
      ### Fit a mediator model: algebra1 ###
      step1.fit1 <- glm(algebra1 ~ female + S1M8GRADE + X1MTHEFF + X1MTHINT,
                        family = binomial(), data = DATA_white)
      
      ### algebra1: Predicted Probability after enforcing Black ###
      nu_m1 <- ifelse(DATA_black$algebra1 == "0",
                      1-predict(step1.fit1, newdata=DATA_black, type="response"),
                      predict(step1.fit1, newdata=DATA_black ,type='response'))
      
      ### advanced: Predicted Probability after enforcing Black ###
      nu_m2 <- ifelse(DATA_black$advanced == "0",
                      1-predict(step1.fit2, newdata=DATA_black, type="response"),
                      predict(step1.fit2, newdata=DATA_black ,type='response'))
      
      ### Step 2 ###
      ### Denominator ###
      ## Fit a mediator model: advanced ##
      ## Among Black students (R=1) ##
      fit2 <- glm(advanced ~ female + S1M8GRADE + ses + X1PAR1OCC_STEM1 + X1PAR2OCC_STEM1
                  + P1EDUEXPECT + P1EDUASPIRE + M1SEX + M1INTEREST + locale
                  + X1SCHOOLCLI + X2PROBLEM + M1UNPREPPCT + requirements + schextracur
                  + friend + X1MTHINT + X1MTHUTI + X1MTHEFF + X1SCHOOLENG + X1SCHOOLBEL + X1TXMTH
                  + algebra1, family = binomial(logit), data = DATA_black)
      
      ### Predicted probability of the mediator advanced among R = 1 ###
      pi_m2 <-  ifelse(DATA_black$advanced == "0",
                       1 - predict(fit2, type = "response"),
                       predict(fit2, type = "response"))
      
      ### Fit a mediator model: algebra1 ###
      ## Among Black students ##
      fit1 <- glm(algebra1 ~ female + S1M8GRADE + ses
                  + P1EDUEXPECT + P1EDUASPIRE + M1SEX + M1INTEREST + locale
                  + X1SCHOOLCLI + X2PROBLEM + M1UNPREPPCT + requirements + schextracur
                  + friend + X1MTHINT + X1MTHUTI  + X1SCHOOLENG + X1SCHOOLBEL,
                  family = binomial(logit), data = DATA_black)
      
      ### Predicted probability of the mediator algebra1 among R = 1 ###
      pi_m1 <-  ifelse(DATA_black$algebra1 == "0",
                       1 - predict(fit1, type = "response"),
                       predict(fit1, type = "response"))
      
      ### Truncate ###
      w_ice.denom <- 1/(pi_m1*pi_m2)
      w_ice.denom_thr <- quantile(w_ice.denom, probs=c(0.01, 0.99)) # 1.05, 15.53
      w_ice.denom_tr <- ifelse(w_ice.denom < w_ice.denom_thr[1], as.numeric(w_ice.denom_thr[1]), w_ice.denom)
      w_ice.denom_tr <- ifelse(w_ice.denom_tr > w_ice.denom_thr[2], as.numeric(w_ice.denom_thr[2]), w_ice.denom_tr)
      
      ### Step 3 ###
      ## Calculate weights ##
      W_ice.trunc <- (nu_m1 * nu_m2 * w_ice.denom_tr)
      
      ### Fit a weighted regression of Y on C among black (R=1) ###
      fit.lm.trunc <- lm(X2TXMTH ~ female, data=DATA_black, weights = W_ice.trunc)
      W_ice.Y.trunc <- coef(fit.lm.trunc)[1]
      
      delta_ice1 <- as.numeric(W_Y1 - W_ice.Y.trunc)
      zeta_ice0 <- as.numeric(W_ice.Y.trunc - W_Y0)
      
      ### Return values ###
      method <- c("Initial", "CDE_zeta10", "CDE_zeta11",
                  "IME_delta", "IME_zeta",
                  "ICDE_zeta", 'ICDE_zeta.wm',
                  "IIE_delta", "IIE_zeta", "IIE_delta.w", "IIE_zeta.w", 
                  "ICE_delta", "ICE_zeta")
      
      values <- c(tau, zeta_cde10, zeta_cde11,
                  delta_ime, zeta_ime,
                  zeta_icde, zeta_icde.wm, 
                  delta_iie1, zeta_iie0, delta_iie1.wm, zeta_iie0.wm,
                  delta_ice1, zeta_ice0)
      
      names(values) <- method
      return(values)
    },
    #if an error occurs, tell me the error
    error=function(e) {
      message('An Error Occurred')
      print(e)
      return(NA)
    },
    #if a warning occurs, tell me the warning
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
  )
}

Date_Time <- format(Sys.time(),"%Y-%m-%d-%H-%M")
R <- 1000
boot.result <- boot(data=DATA1, statistic=boot_ipw, R=R, parallel = "multicore",
                    ncpus = getOption("boot.ncpus",detectCores()))  

cc.boot.result <- boot.result$t[complete.cases(boot.result$t),]
result <- data.frame(original=boot.result$t0,
                     boot_mean=colMeans(cc.boot.result),
                     boot_sd=apply(cc.boot.result, 2, sd))

save_dir <- paste("HSLS_Boot.Result(R=", R, "): ", Date_Time, ".RData", sep="")
save(list=c("boot.result", "result"), file = save_dir, version = 2) 