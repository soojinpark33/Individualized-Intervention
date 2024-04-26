#::::: R code for the simulation study in Supplementary Materials for 
#::::: Causal Decomposition Analysis with Time-Varying Mediators:
#::::: Designing Individualized Interventions to Reduce Social Disparities
#::::: This R code is written by Soojin Park and Namhwa Lee.

### Required Package ###
library(DynTxRegime); library(parallel); library(rpart); library(rpart.plot); library(dplyr)

gen.data <- function(n, seed){
  ###################################
  ### Generate Simulation Dataset ###
  ###################################
  if (!is.null(seed)) set.seed(seed = seed)
  # Generate a baseline covariates
  C <- rbinom(n, 1, 0.4)
  # Generate a Race variable
  pR1_C <- exp(1-0.7*C)/(1+exp(1-0.7*C)) 
  R <- rbinom(n, 1, pR1_C) 
  # Generate three covariates at the 1st time point
  for (i in 1:3){
    assign(paste("X1", i, sep=""), rnorm(n))
  }
  # Generate M1 #
  pM1_RX <- plogis(0.5 + 0.5*X11+0.5*X12 - 0.5*R - 0.5*C)
  M1 <- rbinom(n, 1, prob=pM1_RX)
  # Generate M1opt #
  M1opt <- (X11 > -0.54) * (X12 < 0.54)
  # Generate three covariates at the 2nd time points
  X21 <- 0.2*C -1.4*R + 0.8*X11 + 0.6*M1 + rnorm(n)
  X22 <- 0.2*C -0.4*R + 0.8*X12 + 0.6*M1 + rnorm(n)
  X23 <- 0.2*C -1.4*R + 0.8*X13 + 0.6*M1 + rnorm(n)
  # Generate M1 #
  pM2_RX <- plogis(-0.5 + 0.5*X21 + 0.5*X12 - 0.5*R - 0.5*C)
  M2 <- rbinom(n, 1, prob=pM2_RX)
  # Generate M1opt #
  M2opt <- (X21 > 0.3) * (X23 < 0.46)
  # Generate Y #
  y <- 1 + 0.25*X11 + 0.25*X12 - 0.25*X13 - 0.5*(M1 - M1opt)^2 - (M2 - M2opt)^2  + 0.25*C - 0.5*R + rnorm(n)
  y_ice <- 1 + 0.25*X11 + 0.25*X12 - 0.25*X13 - 0.5*M1 + 0.5*M1*X11 - M2 + 0.7*M2*X12   + 0.25*C - 0.5*R + rnorm(n)

    # Generate a dataframe
  data <- data.frame(R=R, C=C, C.centered=scale(C),
                     X11=X11, X12=X12, X13=X13, M1=M1, M1opt=M1opt, 
                     X21=X21, X22=X22, X23=X23, M2=M2, M2opt=M2opt,
                     Y=y, Y_ice=y_ice,
                     true.Ind1=(M1==M1opt), true.Ind2=(M2==M2opt),
                     true.M1prob=pM1_RX, true.M2prob=pM2_RX)
  
  ###########################
  ### Calculate True ICDE ###
  ###########################
  EY11_R1C0 <- (1+0.25*mean((X11 + X12 - X13)[R==1 & C==0])-0.5)*mean((M1opt & M2opt)[R==1 & C==0])
  EY10_R1C0 <- (1+0.25*mean((X11 + X12 - X13)[R==1 & C==0])-0.5)*mean((M1opt & !M2opt)[R==1 & C==0])
  EY01_R1C0 <- (1+0.25*mean((X11 + X12 - X13)[R==1 & C==0])-0.5)*mean((!M1opt & M2opt)[R==1 & C==0])
  EY00_R1C0 <- (1+0.25*mean((X11 + X12 - X13)[R==1 & C==0])-0.5)*mean((!M1opt & !M2opt)[R==1 & C==0])
  
  EY_R1C0 <- EY11_R1C0 + EY10_R1C0 + EY01_R1C0 + EY00_R1C0
  
  EY11_R0C0 <- (1+0.25*mean((X11 + X12 - X13)[R==0 & C==0]))*mean((M1opt & M2opt)[R==0 & C==0])
  EY10_R0C0 <- (1+0.25*mean((X11 + X12 - X13)[R==0 & C==0]))*mean((M1opt & !M2opt)[R==0 & C==0])
  EY01_R0C0 <- (1+0.25*mean((X11 + X12 - X13)[R==0 & C==0]))*mean((!M1opt & M2opt)[R==0 & C==0])
  EY00_R0C0 <- (1+0.25*mean((X11 + X12 - X13)[R==0 & C==0]))*mean((!M1opt & !M2opt)[R==0 & C==0])
  
  EY_R0C0 <- EY11_R0C0 + EY10_R0C0 + EY01_R0C0 + EY00_R0C0
  
  zeta_icde.true <- EY_R1C0 - EY_R0C0
  
  ##########################
  ### Calculate True IIE ###
  ##########################
  EY_R1C <- data %>% filter(R == 1, C == 0) %>% summarise(EY_R1c=mean(Y))
  
  EY_R0C <- data %>% filter(R == 0, C == 0) %>% summarise(EY_R0c=mean(Y)) 
  
  wm.Y11 <- data %>% filter(R == 1, C == 0) %>% summarise(wm.Y11=mean(1 + 0.25*(X11 + X12 - X13) - 0.5))
  
  I11 <- data %>% filter(R == 0 , C == 0) %>% summarise(mean(true.Ind1) * mean(true.Ind2))
  
  wm.Y10 <- data %>% filter(R == 1, C == 0) %>% summarise(wm.Y10=mean(1 + 0.25*(X11 + X12 - X13) - 1 - 0.5))
  
  I10 <- data %>% filter(R == 0 , C == 0) %>% summarise(mean(true.Ind1) * mean(!true.Ind2))
  
  wm.Y01 <- data %>% filter(R == 1, C == 0) %>% summarise(wm.Y11=mean(1 + 0.25*(X11 + X12 - X13) - 0.5 - 0.5))
  
  I01 <- data %>% filter(R == 0 , C == 0) %>% summarise(mean(!true.Ind1) * mean(true.Ind2))
  
  wm.Y00 <- data %>% filter(R == 1, C == 0) %>% summarise(wm.Y11=mean(1 + 0.25*(X11 + X12 - X13) - 0.5 - 1 - 0.5))
  
  I00 <- data %>% filter(R == 0 , C == 0) %>% summarise(mean(!true.Ind1) * mean(!true.Ind2))
  
  EY_K1K2 <- wm.Y11*I11 + wm.Y10*I10 + wm.Y01*I01 + wm.Y00*I00
  delta_iie1.true <- EY_R1C - EY_K1K2
  zeta_iie0.true <- EY_K1K2 - EY_R0C
  
  ##########################
  ### Calculate True ICE ###
  ##########################
  EY.ice_R1C <- data %>% filter(R == 1, C == 0) %>% summarise(mean(Y_ice))
  EY.ice_R0C <- data %>% filter(R == 0, C == 0) %>% summarise(mean(Y_ice))
  
  ### P(X|R=1, C=0) 
  PX11_R1C0 <- data %>% filter(R == 1, C == 0) %>% summarise(mean(X11))
  PX12_R1C0 <- data %>% filter(R == 1, C == 0) %>% summarise(mean(X12))
  PX13_R1C0 <- data %>% filter(R == 1, C == 0) %>% summarise(mean(X13))
  
  ### M1 ###
  fit_m1_rx <- glm(M1 ~  R + C + X11 , data=data, family=binomial(logit))
  m1_rx.coef <- fit_m1_rx$coef
  pm1_R0x <- plogis(as.numeric(m1_rx.coef[1] + m1_rx.coef[3]*C + m1_rx.coef[4]*X11))
  
  ### M2 ###
  fit_m2_rx <- glm(M2 ~  R + C + X12 , data=data, family=binomial(logit))
  m2_rx.coef <- fit_m2_rx$coef
  pm2_R0X <- plogis(as.numeric(m2_rx.coef[1] + m2_rx.coef[3]*C + m2_rx.coef[4]*X12))

  ### Y ###
  term1 <- 1 - 0.5 + 0.25*(PX11_R1C0 + PX12_R1C0 - PX13_R1C0)
  term2 <- - 0.5*mean(pm1_R0x[R==1 & C==0]) + 0.5*mean(X11[R==1 & C==0]*pm1_R0x[R==1 & C==0])
  term3 <- - mean(pm2_R0X[R==1 & C==0]) + 0.7*mean(X12[R==1 & C==0]*pm2_R0X[R==1 & C==0])
  EY_G1G2 <- as.numeric(term1 + term2 + term3)
  delta_ICE.true <- as.numeric(EY.ice_R1C -  EY_G1G2)
  zeta_ICE.true <- as.numeric(EY_G1G2 - EY.ice_R0C)
  
  True <- cbind(ICDE=as.numeric(zeta_icde.true), 
                delta_IIE=as.numeric(delta_iie1.true), zeta_IIE=as.numeric(zeta_iie0.true),
                delta_ICE=as.numeric(delta_ICE.true), zeta_ICE=as.numeric(zeta_ICE.true))
  return(list(pop.data=data, True=True))
}

population <- gen.data(n=1e6, seed=2024)
pop.data <- population$pop.data
true.value <- population$True

N <- 500 # Number of replication
nobs <- 500 # Number of observations: c(500, 1000, 2000)

sample.data <- function(seed, pop.data=pop.data, nobs){
  set.seed(seed)
  n <- nrow(pop.data)
  idx <- sample(1:n, nobs, replace = TRUE)
  samp.data <- pop.data[idx, ]
  return(samp.data)
}

sample_data <- lapply(X=1:N, FUN=sample.data, pop.data=pop.data, nobs=nobs)

parallel.compare <- function(data) {
  ##### Weighting Methods for Estimating Optimal DTR #####
  # IPW for Stage 2 #
  moPropen2 <- buildModelObj(model = ~ R + C + X21 + X12, 
                             solver.method = "glm",
                             solver.args = list("family"="binomial"),
                             predict.method = 'predict.glm',
                             predict.args = list(type="response"))
  
  moClass2 <- buildModelObj(model = ~ X21 + X23,
                            solver.method = "rpart",
                            solver.args = list(method="class"),
                            predict.args = list(type='class'))
  
  fitSS_IPW <- optimalClass(moPropen = moPropen2,
                            moClass = moClass2,
                            data=data, response=data$Y, txName="M2", verbose = F)
  
  optSS <- optTx(fitSS_IPW)$optimalTx
  # IPW for STAGE 1 #
  moPropen1 <- buildModelObj(model = ~ R + C + X11 + X12, 
                             solver.method = "glm",
                             solver.args = list("family"="binomial"),
                             predict.method = 'predict.glm',
                             predict.args = list(type="response"))
  
  moClass1 <- buildModelObj(model = ~X11 + X12, 
                            solver.method = "rpart",
                            solver.args = list(method="class"),
                            predict.args = list(type='class'))
  
  fitFS_IPW <- optimalClass(moPropen = moPropen1, 
                            moClass = moClass1,
                            data=data, response = fitSS_IPW, txName = "M1", verbose = F)
  
  optFS <- ifelse(is.na(optTx(fitFS_IPW)$optimalTx), 0, optTx(fitFS_IPW)$optimalTx)

  ### Assign optimal decision ###
  data$opt_M1 <- optFS
  data$opt_M2 <- optSS
  DATA1 <- data
  
  ##########################
  ########## ICDE ##########
  ##########################
  ### index for compliance sample ###
  DATA1 <- DATA1 %>% mutate(Ind1 = (M1 == opt_M1), Ind2=(M2 == opt_M2))
  ind1 <- (DATA1$M1 == DATA1$opt_M1)
  ind2 <- (DATA1$M2 == DATA1$opt_M2)
  
  ### Construct weight model for med2: Pr(M2 | R, C, X1, M1, X2) ###
  fit.m2 <- glm(M2 ~ R + C + X12 + X21 , family = binomial(logit), data = DATA1)
  p.med2 <- ifelse(DATA1$M2 == 0,
                   1 - predict(fit.m2, type = "response"),
                   predict(fit.m2, type = "response"))
  
  ### construct weight for med1 ###
  fit.m1 <- glm(M1 ~ R + C + X11 + X12 , family = binomial(logit), data = DATA1)
  p.med1 <- ifelse(DATA1$M1 == 0, 
                   1 - predict(fit.m1, type = "response"), 
                   predict(fit.m1, type = "response"))
  
  ### Weights ###
  w2 <- 1/p.med2 ; w1 <- 1/p.med1
  w <- w1 * w2
  
  # Calculate zeta_icde #
  fit <- lm(Y ~ R + C, data=DATA1, weights = (ind1 & ind2)*w)
  zeta_icde <- coef(fit)[2]
  
  #########################################
  ########## ICDE: Weighted Mean ##########
  #########################################
  E_W.icdeY_R1C0 <- DATA1 %>% 
    mutate(w_icde=w) %>%  
    filter(R == 1, C == 0) %>% 
    summarise(weighted.mean(Y, Ind1*Ind2*w_icde))
  
  E_W.icdeY_R0C0 <- DATA1 %>%
    mutate(w_icde=w) %>% 
    filter(R == 0, C == 0) %>% 
    summarise(weighted.mean(Y, Ind1*Ind2*w_icde))
  
  zeta_icde.wm <- as.numeric(E_W.icdeY_R1C0 - E_W.icdeY_R0C0)

  #########################################
  ############# IIE: Method 1 #############
  #########################################
  DATA_R0 <- subset(DATA1, R==0) # white
  DATA_R1 <- subset(DATA1, R==1) # black
  
  ### Step 1 ###
  # Among white (R=0) & Among black (R=1) #
  fit.Im1_0c <- glm(Ind1 ~ C, data = DATA_R0, family = binomial(logit))
  fit.Im2_0c <- glm(Ind2 ~ C, data = DATA_R0, family = binomial(logit))
  
  # P [ I(M1=d1.opt)=1 | R=0, C=c ] # 
  pi_Im11_0c <- plogis(coef(fit.Im1_0c)[1])
  # P [ I(M1=d1.opt)=0 | R=0, C=c ] # 
  pi_Im10_0c <- 1 - pi_Im11_0c
  # P [ I(M2=d2.opt)=1 | R=0, C=c ] # 
  pi_Im21_0c <- plogis(coef(fit.Im2_0c)[1])
  # P [ I(M2=d2.opt)=0 | R=0, C=c ] # 
  pi_Im20_0c <- 1 - pi_Im21_0c
  
  ### Step 2 ###
  W_IIE.11 <- (ind1 & ind2) * w
  # theta1 = 1, theta2 = 0 : M1=d1.opt & M2!=d2.opt #
  W_IIE.10 <- (ind1 & !ind2) * w
  # theta1 = 0, theta2 = 1 : M1!=d1.opt & M2=d2.opt #
  W_IIE.01 <- (!ind1 & ind2) * w
  # theta1 = 0, theta2 = 0 : M1!=d1.opt & M2!=d2.opt #
  W_IIE.00 <- (!ind1 & !ind2) * w
  
  # E [Y | R=1, C=c] #
  fit.lm1 <- lm(Y ~  C, data=DATA_R1)
  W_Y1 <- coef(fit.lm1)[1]
  # E [Y | R=0, C=c] #
  fit.lm0 <- lm(Y ~ C, data=DATA_R0)
  W_Y0 <- coef(fit.lm0)[1]
  
  # (theta1, theta1) = (1,1) #
  idx_R1 <- which(DATA1$R == 1)
  fit.lm11 <- lm(Y ~ C, data=DATA_R1, weights = W_IIE.11[idx_R1])
  term11 <- pi_Im11_0c * pi_Im21_0c * coef(fit.lm11)[1]
  # (theta1, theta1) = (1,0) #
  fit.lm10 <- lm(Y ~ C, data=DATA_R1, weights = W_IIE.10[idx_R1])
  term10 <- pi_Im11_0c * pi_Im20_0c * coef(fit.lm10)[1]
  # (theta1, theta1) = (0,1) #
  fit.lm01 <- lm(Y ~ C, data=DATA_R1, weights = W_IIE.01[idx_R1])
  term01 <- pi_Im10_0c * pi_Im21_0c * coef(fit.lm01)[1]
  # (theta1, theta1) = (0,0) #
  fit.lm00 <- lm(Y ~ C, data=DATA_R1, weights = W_IIE.00[idx_R1])
  term00 <- pi_Im10_0c * pi_Im20_0c * coef(fit.lm00)[1]
  
  W_iie.Y <- term11 + term10 + term01 + term00
  
  delta_iie1 <- W_Y1 - W_iie.Y 
  zeta_iie0 <- W_iie.Y - W_Y0
  
  #########################################
  ############# IIE: Method 2 #############
  #########################################
  fit.lm <- lm(Y ~ R+ C + Ind1 + Ind2, data=DATA1, weights=w)
  fit.Ind1 <- glm(Ind1 ~ R + C, data = DATA1, family = binomial(logit))
  fit.Ind2 <- glm(Ind2 ~ R + C, data = DATA1, family = binomial(logit))
  alpha1 <- plogis(sum(coef(fit.Ind1)[1:2])) - plogis(coef(fit.Ind1)[1])
  alpha2 <- plogis(sum(coef(fit.Ind2)[1:2])) - plogis(coef(fit.Ind2)[1])
  
  delta_iie1_2 <- alpha1*coef(fit.lm)["Ind1TRUE"]+alpha2*coef(fit.lm)["Ind2TRUE"]
  zeta_iie0_2 <- (W_Y1 - W_Y0) - delta_iie1_2
  
  #########################################
  ############# ICE: Method 1 #############
  #########################################
  ### Step 1 ###
  ### Numerator ###
  ## Fit a mediator model : M2 ##
  step1.fit2 <- glm(M2 ~ C + X12, family = binomial(logit), data = DATA_R0)
  ### Fit a mediator model: M1 ###
  step1.fit1 <- glm(M1 ~ C + X11, family = binomial(logit), data = DATA_R0)
  ### algebra1: Predicted Probability after enforcing Black ###
  nu_m1 <- ifelse(DATA_R1$M1 == 0,
                  1-predict(step1.fit1, newdata=DATA_R1, type="response"),
                  predict(step1.fit1, newdata=DATA_R1 ,type='response'))
  
  ### advanced: Predicted Probability after enforcing Black ###
  nu_m2 <- ifelse(DATA_R1$M2 == 0,
                  1-predict(step1.fit2, newdata=DATA_R1, type="response"),
                  predict(step1.fit2, newdata=DATA_R1 ,type='response'))
  
  ### Denominator ###
  ## Fit a mediator model: M2 ##
  ## Among Black students (R=1) ##
  fit2 <- glm(M2 ~ C + X12 + X21, family = binomial(logit), data = DATA_R1)
  ### Predicted probability of the mediator M2 among R = 1 ###
  pi_m2 <-  ifelse(DATA_R1$M2 == 0,
                   1 - predict(fit2, type = "response"),
                   predict(fit2, type = "response"))
  
  ### Fit a mediator model: algebra1 ###
  ## Among Black students ##
  fit1 <- glm(M1 ~ C + X11 + X12, family = binomial(logit), data = DATA_R1)
  ### Predicted probability of the mediator M1 among R = 1 ###
  pi_m1 <-  ifelse(DATA_R1$M1 == 0,
                   1 - predict(fit1, type = "response"),
                   predict(fit1, type = "response"))
  
  w_ice.denom <- 1/(pi_m1*pi_m2)
  
  ### Step 3 ###
  ## Calculate weights ##
  W_ice <- (nu_m1 * nu_m2 * w_ice.denom)
  
  ### Fit a weighted regression of Y on C among black (R=1) ###
  fit.lm <- lm(Y_ice ~ C, data=DATA_R1, weights = W_ice)
  W_ice.Y <- coef(fit.lm)[1]
  
  fit.lm1 <- lm(Y_ice ~ C, data=DATA_R1)
  W_Y1 <- coef(fit.lm1)[1]
  fit.lm0 <- lm(Y_ice ~ C, data=DATA_R0)
  W_Y0 <- coef(fit.lm0)[1] 
  
  ### Step 4 ### 
  delta_ice1 <- as.numeric(W_Y1 - W_ice.Y)
  zeta_ice0 <- as.numeric(W_ice.Y - W_Y0)
  
  return(list(zeta_icde.wm=zeta_icde.wm, zeta_icde=zeta_icde,
              delta_iie1=delta_iie1, zeta_iie0=zeta_iie0,
              delta_iie1_2=delta_iie1_2, zeta_iie0_2=zeta_iie0_2,
              delta_ice1=delta_ice1, zeta_ice0=zeta_ice0))
}

## Date-Time ####
Date_Time <- format(Sys.time(),"%Y-%m-%d-%H-%M") # Used as a unique identifier
ncores <- 1
### Parallel Computing ###
start <- Sys.time()# Used for timing process
parallel.results <- mclapply(sample_data, parallel.compare , mc.cores = ncores) # Setting the number of cores to use
print("Parallel: mclapply")
(comp.time <- Sys.time()-start) # Time it took

zeta_icde.wm = zeta_icde <- matrix(ncol=1, nrow=N)
delta_iie1 = zeta_iie0 = delta_iie1_2 = zeta_iie0_2 <- matrix(ncol=1, nrow=N)
delta_ice1 = zeta_ice0  <- matrix(ncol=1, nrow=N)

for (i in 1:N){
  # ICDE #
  zeta_icde.wm[i,] <- parallel.results[[i]]$zeta_icde.wm # Weighted mean 
  zeta_icde[i,] <- parallel.results[[i]]$zeta_icde # Regression
  
  # IIE #
  delta_iie1[i,] <- as.numeric(parallel.results[[i]]$delta_iie1) # delta estimate from the method 1
  zeta_iie0[i,] <- as.numeric(parallel.results[[i]]$zeta_iie0) # zeta estimate from the method 1
  
  delta_iie1_2[i,] <- as.numeric(parallel.results[[i]]$delta_iie1_2) # delta estimate from the method 2 
  zeta_iie0_2[i,] <- as.numeric(parallel.results[[i]]$zeta_iie0_2) # zeta estimate from the method 2

  # ICE #
  delta_ice1[i,] <- as.numeric(parallel.results[[i]]$delta_ice1) # delta estimate
  zeta_ice0[i,] <- as.numeric(parallel.results[[i]]$zeta_ice0) # zeta estimate 
}

ICDE_result <- cbind(zeta_icde.wm, zeta_icde)
colnames(ICDE_result) <- c("WM", "Reg")

IIE_result <- cbind(delta_iie1, zeta_iie0,
                    delta_iie1_2, zeta_iie0_2)

colnames(IIE_result) <- c("Est1_delta", "Est1_zeta", 
                          "Est2_delta", "Est2_zeta")

ICE_result <- cbind(delta_ice1, zeta_ice0)

colnames(ICE_result) <- c("Est1_delta", "Est1_zeta")

Est_result <- list(ICDE=ICDE_result, IIE=IIE_result, ICE=ICE_result)

names(true.value) <- c("ICDE", "delta_IIE", "zeta_IIE", "delta_ICE", "zeta_ICE")

# Bias #
bias_ICDE <- true.value["ICDE"] - colMeans(Est_result$ICDE)
bias_IIE.delta <- true.value["delta_IIE"] - colMeans(Est_result$IIE[,c(1,3)])
bias_IIE.zeta <- true.value["zeta_IIE"] - colMeans(Est_result$IIE[,c(2,4)])
bias_ICE.delta <- true.value["delta_ICE"] - colMeans(Est_result$ICE)[1]
bias_ICE.zeta <- true.value["zeta_ICE"] - colMeans(Est_result$ICE)[2]

Bias <- list(ICDE=bias_ICDE, IIE=c(bias_IIE.delta, bias_IIE.zeta), ICE=c(bias_ICE.delta, bias_ICE.zeta))

# RMSE $
rmse_ICDE <- sqrt(colMeans((true.value["ICDE"]- Est_result$ICDE)^2))
rmse_IIE.delta <- sqrt(colMeans((true.value["delta_IIE"] - Est_result$IIE[,c(1,3)])^2))
rmse_IIE.zeta <- sqrt(colMeans((true.value["zeta_IIE"] - Est_result$IIE[,c(2,4)])^2))
rmse_ICE.delta <- sqrt(mean((true.value["delta_ICE"] - Est_result$ICE[,1])^2))
rmse_ICE.zeta <- sqrt(mean((true.value["zeta_ICE"] - Est_result$ICE[,2])^2))

RMSE <- list(ICDE=rmse_ICDE, IIE=c(rmse_IIE.delta, rmse_IIE.zeta), ICE=c(rmse_ICE.delta, rmse_ICE.zeta))

Result <- list(Estimation=Est_result, True=true.value, Bias=Bias, RMSE=RMSE)

save_dir <- paste("Pop.Simulation(nobs=", nobs, "): ", Date_Time, ".RData", sep="") # Creating file name, contains date-time
save(Result, file = save_dir, version = 2) # Saving RData file, recommend using version 2