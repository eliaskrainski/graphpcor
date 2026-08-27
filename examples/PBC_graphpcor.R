
library(INLA)
library(INLAjoint)
library(graphpcor)
library(ggplot2)
library(ggpubr)

## Followup of 312 randomised patients with primary biliary
##   cirrhosis, a rare autoimmune liver disease, at Mayo Clinic.
data(pbc2, package = "JM")

## covariates
xnames <- c("drug", "age", "sex")

## categorial variables
ycat <- c("ascites", "hematomgaly", "spiders", "edema")

## numeric variables
ynum <- c("serBilir", "serChol", "albumin", "alkaline",
          "SGOT", "platelets", "prothrombin", "histologic")


## reshape
library(dplyr)
library(tidyr)

pbc_long <- pbc2 |>
  select(id, year, drug, all_of(ynum)) |>
  pivot_longer(
    cols = all_of(ynum),
    names_to = "variable",
    values_to = "value"
  )

ggplot(pbc_long,
       aes(year, value, group = id)) +
    geom_line(aes(color = drug),
              alpha = 0.75, linewidth = 0.25) +
  geom_smooth(
      aes(group = drug),
      method = "loess",
      se = FALSE,
      colour = drug,
      linewidth = 1
  ) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_bw()

pbc2_diff <- pbc2 %>%
  group_by(id) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(
    across(
      all_of(ynum),
      ~ .x - .x[which(year == 0)[1]],
      .names = "{.col}"
    )
  ) %>%
    ungroup()
pbc2_diff_long <- pbc2_diff |>
  select(id, year, drug, all_of(ynum)) |>
  pivot_longer(
    cols = all_of(ynum),
    names_to = "variable",
    values_to = "value"
  )

ggplot(pbc2_diff_long,
       aes(year, value, group = id)) +
  geom_line(alpha = 0.15, linewidth = 0.25) +
  geom_smooth(
    aes(group = 1),
    method = "loess",
    se = FALSE,
    colour = "red",
    linewidth = 1
  ) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_bw()

ggplot(
  pbc2_diff_long,
  aes(
    x = year,
    y = value,
    group = id
  )
) +
  geom_hline(yintercept = 0, colour = "grey75") +
  geom_line(alpha = 0.15, linewidth = 0.25) +
  geom_smooth(
    aes(colour = drug),
    method = "loess",
    se = FALSE,
    linewidth = 1.1
  ) +
  facet_wrap(~variable, scales = "free_y") +
  labs(
    x = "Years",
    y = "Change from baseline",
    colour = "Drug"
  ) +
    theme_bw()


########################
########################

library(dplyr)
library(tidyr)
library(ggplot2)

## Numeric variables that were differenced
num_vars <- names(pbc2_diff)[sapply(pbc2_diff, is.numeric)]
num_vars <- setdiff(num_vars, c("id", "year"))
num_vars <- ynum

## Long format, keeping the treatment assignment
pbc_long_diff <- pbc2_diff %>%
  select(id, drug, year, all_of(num_vars)) %>%
  pivot_longer(
    cols = all_of(num_vars),
    names_to = "variable",
    values_to = "change"
  )

ggplot(
  pbc_long_diff,
  aes(
    x = year,
    y = change,
    group = id
  )
) +
  geom_hline(yintercept = 0, colour = "grey75") +
  geom_line(alpha = 0.15, linewidth = 0.25) +
  geom_smooth(
    aes(colour = drug),
    method = "loess",
    se = FALSE,
    linewidth = 1.1
  ) +
  facet_wrap(~variable, scales = "free_y") +
  labs(
    x = "Years",
    y = "Change from baseline",
    colour = "Drug"
  ) +
  theme_bw()

ggplot(
  pbc_long_diff,
  aes(
    x = year,
    y = change,
    group = id,
    colour = drug
  )
) +
  geom_hline(yintercept = 0, colour = "grey75") +
  geom_line(alpha = 0.10, linewidth = 0.25) +
  geom_smooth(
    aes(group = drug),
    method = "loess",
    se = FALSE,
    linewidth = 1.2
  ) +
  facet_wrap(~variable, scales = "free_y") +
  labs(
    x = "Years",
    y = "Change from baseline",
    colour = "Drug"
  ) +
    theme_bw()

## visualize 8 longitudinal vriables
gg0 <- ggplot(data = pbc2) +
    scale_y_log10()
ggarrange(
    gg0 + geom_line(aes(x = year, y = serBilir, group = id)),
    gg0 + geom_line(aes(x = year, y = serChol, group = id)),
    gg0 + geom_line(aes(x = year, y = albumin, group = id)),
    gg0 + geom_line(aes(x = year, y = alkaline, group = id)),
    gg0 + geom_line(aes(x = year, y = SGOT, group = id)),
    gg0 + geom_line(aes(x = year, y = platelets, group = id)),
    gg0 + geom_line(aes(x = year, y = prothrombin, group = id)),
    gg0 + geom_line(aes(x = year, y = histologic, group = id))
)

## how much of this correlation is
##  - 1: "intercept":
##    - 1a) data at the first time
##    - 1b) average data over time
##  - 2: "trajectory"
##    how the profile evolution correlates

## raw correlation
ycorr <- cor(pbc2[ynum], use = 'pair')
round(ycorr*100)

## correlation at first time
table(duplicated(pbc2$id), pbc2$year>0)
it1 <- pbc2$year==0
ycorr0 <- cor(pbc2[it1, ynum], use = 'pair')
round(ycorr0 * 100)

ycorr0 / ycorr 


## number of observations per individual: 1 to 16
table(table(pbc2$id)) ## mode: 44 id have 4 obs.
summary(as.data.frame(table(pbc2$id))$Freq) ## av. 6.234 obs. per id

## however, thera ar missing
ynames <- 

## extract some variable of interest without missing values
Longi <- na.omit(pbc2[, c("id", "years", "status","drug","age",
                          "sex","year","serBilir","SGOT", "albumin", "edema",
                          "platelets", "alkaline","spiders", "ascites")])

f1 <- function(x) x^2

### first prepare data with correct format using dataOnly and full sample
run0 <- joint(
    formLong = list(serBilir ~ (1 + year + f1(year)) * drug + (1 + year + f1(year)| id),
                    platelets ~ (1 + year + f1(year)) * drug + (1 + year + f1(year)| id)),
    dataLong = Longi, id = "id", timeVar = "year",
    family = c("lognormal", "poisson"),
    control=list(int.strategy = "eb", cfg = TRUE),
    corLong = TRUE, run = FALSE)
run.kd <- joint.run(run0)
vars.kd <- summary(run.kd, sdcor=T)$ReffList[[1]]

## graph-based structure model for partial correlations
G <- graphpcor(a1 ~ a2 + a3 + b1, ## 1st intercept with 1st slope and 2nd intercept
               b1 ~ b2 + b3) ## 2nd intercept with 2nd slope
plot(G)

## the correlation model
gmodel1 <- cgeneric(
  model = G,
  sigma.prior.reference = rep(5, dim(G)[1]),
  sigma.prior.probability = rep(0.2, dim(G)[1]),
  lambda = 3
)

## because of the indexing we need to do a "Kronecker":
## within individual correlation model (X) between individuals "iid" model
iidmodel <- cgeneric(
    model = "iid",
    n = 312, ## number of individuals
    param = c(1, 0.0)
)

## the model to actually use
gmodel <- kronecker(gmodel1, iidmodel)

## update the formula with the model
run.G <- run0
run0$.args$formula
run.G$.args$formula <-
    Yjoint ~ -1 + Intercept_L1 + year_L1 + f1year_L1 + drugDpenicil_L1 + 
    year.X.drugDpenicil_L1 + f1year.X.drugDpenicil_L1 + Intercept_L2 + 
    year_L2 + f1year_L2 + drugDpenicil_L2 + year.X.drugDpenicil_L2 + 
    f1year.X.drugDpenicil_L2 +
    f(IDIntercept_L1, WIntercept_L1, model = gmodel) +
    f(IDyear_L1, Wyear_L1, copy = "IDIntercept_L1") + 
    f(IDf1year_L1, Wf1year_L1, copy = "IDIntercept_L1") +
    f(IDIntercept_L2, WIntercept_L2, copy = "IDIntercept_L1") +
    f(IDyear_L2, Wyear_L2, copy = "IDIntercept_L1") +
    f(IDf1year_L2, Wf1year_L2, copy = "IDIntercept_L1")

##    Yjoint ~ -1 + Intercept_L1 + year_L1 + drugDpenicil_L1 + year.X.drugDpenicil_L1 +
  ##      Intercept_L2 + year_L2 + drugDpenicil_L2 + year.X.drugDpenicil_L2 +
    ##    f(IDIntercept_L1, WIntercept_L1, model = gmodel) +
      ##  f(IDyear_L1, Wyear_L1, copy = "IDIntercept_L1") +
        ##f(IDIntercept_L2, WIntercept_L2, copy = "IDIntercept_L1") +
        ##f(IDyear_L2, Wyear_L2, copy = "IDIntercept_L1")

run.G <- joint.run(run.G)

# save(JM_INLA_GRAPH, file="JM_INLA_GRAPH.RData")

vfit <- vcov(G, theta = run.G$mode$theta[-1])
round(cov2cor(vfit), 2)

hsamples <- inla.hyperpar.sample(n=1e4, result = run.G, intern = TRUE)
scsamples <- t(sapply(1:nrow(hsamples), function(i) {
    v <- vcov(G, theta = hsamples[i, -1])
    return(c(sqrt(diag(v)), cov2cor(v)[lower.tri(v)]))
}))

gsummary <- data.frame(
    mean = apply(scsamples, 2, mean),
    sd = apply(scsamples, 2, sd),
    q0.025 = apply(scsamples, 2, quantile, 0.025),
    q0.975 = apply(scsamples, 2, quantile, 0.975)
)

round(cbind(vars.kd[, 1:2], gsummary[, 1:2]), 4)

## "12", "13", "14", "15", "16",
##   "23", "24", "25", "26", 
##     "34", "35", "36"
##      "45", "46", "56"
ijc <- list(c(1,2), c(1,3), c(1,4), c(1,5), c(1,6),
            c(2,3), c(2,4), c(2,5), c(2,6),
            c(3,4), c(3,5), c(3,6), c(4,5), c(4,6), c(5,6))
nvars <- 6
ncors <- length(ijc)
npars <- nvars + ncors

par(mfrow = c(1, 2), mar = c(4,4,1,1), mgp = c(3,1.5,0), las = 1, bty = 'n')
plot(1:nvars-0.1, vars.kd[1:nvars, 1], pch = 19, axes = FALSE, log = "y",
     xlim = c(0.5, nvars+.5), ylim = range(vars.kd[1:nvars, 3:5], gsummary[1:nvars, 3:4]),
     xlab = '', ylab = expression(sigma))
axis(2)
axis(1, 1:nvars, as.expression(lapply(1:nvars, function(i) bquote(sigma[.(i)]))))
segments(1:nvars-0.1, vars.kd[1:nvars, 3], 1:nvars-0.1, vars.kd[1:nvars, 5])
points(1:nvars+0.1, gsummary[1:nvars, 1], col = 2, pch = 8)
segments(1:nvars+0.1, gsummary[1:nvars, 3], 1:nvars+0.1, gsummary[1:nvars, 4], col = 2)
plot(1:ncors, vars.kd[(nvars+1):npars, 1], pch = 19, axes = FALSE,
     ylim = range(vars.kd[(nvars+1):npars, 3:4]), xlim = c(-0.5, ncors+.5), 
     xlab = '', ylab = 'Correlation')
axis(2)
axis(1, 1:ncors, as.expression(lapply(ijc, function(ij)
    bquote(rho[.(ij[1])~.(ij[2])]))))
segments(1:ncors-0.1, vars.kd[(nvars+1):npars, 3],
   1:ncors-0.1, vars.kd[(nvars+1):npars, 5])
points(1:ncors+0.1,gsummary[(nvars+1):npars, 1], col = 2, pch = 8)
segments(1:ncors+0.1, gsummary[(nvars+1):npars, 3],
         1:ncors+0.1, gsummary[(nvars+1):npars, 4], col = 2)
abline(h=0)
legend("bottomleft", c("IIDKD", "GraphPCor"), bty = "n",
       lty = 1, col = 1:2, pch = c(19,8))
