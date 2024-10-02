library(glmnet)
library(tidyverse)
ntrials <- 2000
n <- 500
p <- 50

all_res <- data.frame("setting"=NA, "method"=NA, "sel.ind" = NA, "cov"=NA, "pval"=NA)


### Note: We generate orthogonal Xs such that the target parameters are the same
### regardless of the selected model. This is a simplification that is sufficient to make
### the point about the marginal vs. conditional approach. 
set.seed(1)
X <- 10*pracma::randortho(n)[,1:p]
beta0 <- 0.6
beta1 <- -0.9
beta2 <- 2.1
beta3 <- -1.5

### Data generating mechanism for the non-global-null scenario. 
beta <- c(beta1, beta2, beta3, rep(0, p-3))
theta <- exp(beta0+beta1*X[,1]+beta2*X[,2]+beta3*X[,3])/
  (1+exp(beta0+beta1*X[,1]+beta2*X[,2]+beta3*X[,3])) 

### Let epsilon=0.8 such that the majority of the information is in the training set. 
epsilon=0.8

for (t in 1:ntrials) {
  if (t%%50==0) {print(t)}
  set.seed(t)
  
  #####################
  #### GLOBAL NULL ####
  #####################
  y <- rbinom(n, size=1, prob=beta0)
  
  # Step 1
  Z <- rbinom(n, size=1, prob=epsilon)
  y1fission <- y*(1-Z)+(1-y)*Z
  y2fission <- y
  # Step 2
  mod2 = cv.glmnet(X, y1fission, family = "binomial")
  selected.vars = which(coef(mod2, s = 'lambda.1se')!=0)[-1]-1
  
  if (length(selected.vars) > 0) {
    # Step 3, Example 2 method. 
    mod.fiss <- glm(y2fission~X[,selected.vars], family="binomial")
    pvals.leiner <- clubSandwich::coef_test(mod.fiss, vcov = "CR0", cluster = factor(1:n), 
                                            level = 0.95)[,6]
    counter=2
    for (var in selected.vars) {
      pval <-  pvals.leiner[counter]
      counter <- counter+1
      all_res <- all_res %>%  
        add_row(setting="globalnull", method="Example 2", 
                sel.ind = var, cov=pval > 0.05, pval = pval)
    }
    
    # Step 3, Example 2 method. 
    counter <- 2
    offsets <-  rep((1-epsilon)/epsilon, n)
    offsets[y1fission==0] <- epsilon/(1-epsilon)
    mod.fiss <- glm(y2fission~X[,selected.vars]+
                      offset(log(offsets)), family="binomial")
    for (var in selected.vars) {
      pval <- summary(mod.fiss)$coefficients[counter,4]
      counter <- counter+1
      all_res <- all_res %>%  add_row(setting="globalnull",
            method="Example 3", sel.ind = var, cov=pval > 0.05, pval = pval)
    }
  }
  

  #####################
  #### NON NULL ####
  #####################
  y <- rbinom(n, size=1, prob=theta)
  
  ### We skip" method 1 "Example 2" because we already know it does not work!!!! 
  
  #### "Example 3 method"
  
  # Step 1
  Z <- rbinom(n, size=1, prob=epsilon)
  y1fission <- y*(1-Z)+(1-y)*Z
  y2fission <- y
  
  # Step 2
  mod2 = cv.glmnet(X, y1fission, family = "binomial")
  selected.vars = which(coef(mod2, s = 'lambda.1se')!=0)[-1]-1

  # Step 3
  offsets <-  rep((1-epsilon)/epsilon, n)
  offsets[y1fission==0] <- epsilon/(1-epsilon)

  if (length(selected.vars) > 0) {
    mod.fiss <- glm(y2fission~X[,selected.vars]+
                      offset(log(offsets)), family="binomial")
    
    ## Standard inference. 
    counter <- 2
    for (var in selected.vars) {
      CIthin <- suppressMessages(confint(mod.fiss, counter))
      pvalthin <- summary(mod.fiss)$coefficients[counter,4]
      counter <- counter+1
      
      cov.thin <- CIthin[1] <= beta[var] & CIthin[2] >= beta[var]
      all_res <- all_res %>%  add_row(setting="nonnull", method="Example 3", sel.ind = var, cov=cov.thin, pval = pvalthin)
    }
  }
}

all_res <- all_res[-1,]

globalnullres <- all_res %>% filter(setting=="globalnull")

ggplot(data=globalnullres, aes(sample=pval, col=method))+
  geom_qq(distribution=qunif)+geom_abline(a=0,b=1)+theme_bw()+ylab("Sample quantiles")+xlab("Unif(0,1) quantiles")+
  theme(axis.title =element_text(size=14),
        plot.title=element_text(size=16),
        plot.subtitle=element_text(size=16),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))+
  labs(col="Method")+
  ggtitle("Uniform QQ plot of p-values", "under the global null")+
  scale_color_grey()
ggsave("global_null_qq.png", width=5, height=4)

ggplot(data=globalnullres, aes(sample=pval, col=method))+
  geom_qq(distribution=qunif)+geom_abline(a=0,b=1)+theme_bw()+ylab("Sample quantiles")+xlab("Unif(0,1) quantiles")+
  theme(axis.title =element_text(size=14),
        plot.title=element_text(size=16),
        plot.subtitle=element_text(size=16),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))+
  labs(col="Method")+
  ggtitle("Uniform QQ plot of p-values", "under the global null")+
  scale_color_grey()
ggsave("global_null_qq_bw.png", width=5, height=4)

all_res <- all_res %>% filter(setting != "globalnull")

### This is the correct entry for 
all_res %>% group_by(sel.ind, method) %>% summarize("prop"=n()/t) %>% 
  filter(sel.ind > 3) %>% group_by(method) %>%
  summarize(mean(prop))


### Row 4 column 2 is showing how often we selected at least one null variable.
### The above calculation shows  how often we selected a specific null variable, on average. 
all_res <- all_res %>% mutate("Variable" = ifelse(sel.ind > 3, "Null",  paste0("Beta", sel.ind)))
all_res %>% group_by(setting, method, Variable) %>% 
  summarize("Coverage:" = round(mean(cov),2), 
            "Number of Trials:"= n()/t, 
            "Proportion of times rejected:"=mean(pval < 0.05))

