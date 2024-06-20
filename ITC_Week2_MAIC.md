## A. Introduction

### A.1 Indirect Treatment Comparison (ITC)

Indirect Treatment Comparison (ITC) methods are crucial tools in health
technology assessment (HTA) and comparative effectiveness research,
especially when direct head-to-head clinical trial data are unavailable.
Head-to-head clinical trial data refers to information obtained from
clinical trials that directly compare two or more treatments against
each other within the same study population. This type of data is
considered the gold standard in clinical research for evaluating the
relative effectiveness and safety of different interventions. ITC allows
for the comparison of multiple treatments across different studies by
using a common comparator. The standard Bucher method is a simple form
of ITC that estimates the relative treatment effect by leveraging data
from studies that share a common comparator. For example, if Study 1
compares Treatment A to Control C and Study 2 compares Treatment B to
Control C, the indirect comparison of Treatment A and Treatment B can be
calculated as:

*l**o**g*(*H**R*<sub>*A*/*B*</sub>) = *l**o**g*(*H**R*<sub>*A*/*C*</sub>) − *l**o**g*(*H**R*<sub>*B*/*C*</sub>),

where *l**o**g*(*H**R*<sub>*A*/*B*</sub>) is the log hazard ratio
comparing Treatment 1 to Treatment 2,
*l**o**g*(*H**R*<sub>*A*/*C*</sub>) is the log hazard ratio of Treatment
1 versus Control, and *l**o**g*(*H**R*<sub>*B*/*C*</sub>) is the log
hazard ratio of Treatment 2 versus Control. This can be exponentiated to
obtain the HR:

$$
HR\_{A/B}=\frac{HR\_{A/C}}{HR\_{B/C}}.
$$

#### A.1.1 Anchored vs. Unanchored

Anchored ITC uses a common comparator to link the treatments from
different studies. This method relies on having a shared control or
placebo group in both studies, which acts as an anchor. Unanchored ITC
does not use a common comparator. Instead, it attempts to compare the
treatments directly across different studies without an explicit linking
through a shared control group. This approach is less common and more
challenging due to the increased potential for bias and confounding.

#### A.1.2 Individual Patient Data (IPD) and Aggregate Data (AgD)

In clinical research and meta-analyses, data can be presented in two
main forms: Individual Patient Data (IPD) and Aggregate Data (AgD).
Understanding the differences between these two types of data is crucial
for selecting appropriate methods for analysis and for interpreting
results correctly.

IPD provides detailed, patient-level information, including demographic
data (age, gender, etc.), clinical measurements (blood pressure,
biomarkers, etc.), treatment received, and outcomes. AgD provides
summary statistics, such as the mean age of participants, the proportion
of patients achieving a particular outcome, or hazard ratios for
treatment effects.

### A.2 Bucher method

The Bucher method is a widely used approach for indirect treatment
comparisons, especially in the absence of head-to-head clinical trial
data. It relies on the assumption that the trials being compared are
sufficiently similar in terms of their patient populations and study
designs. However, when there are differences in baseline covariates
(e.g., age, gender, disease severity) between the populations in Study 1
and Study 1, the Bucher method does not account for these differences,
potentially leading to biased estimates of the treatment effect. This
bias occurs because the method assumes that the baseline risks and
covariates are similar across the studies, which may not be the case in
reality.

Bucher Method formula for the indirect comparison is
*Δ̂*<sub>*A*/*B*</sub> = log (HR<sub>*A*/*B*</sub>) = *Δ̂*<sub>*A*/*C*</sub> − *Δ̂*<sub>*B*/*C*</sub>,

where *Δ̂*<sub>*A*/*C*</sub> is the estimated marginal treatment effect
of A vs C, and *Δ̂*<sub>*B*/*C*</sub> is the estimated marginal treatment
effect of B vs C.

Variance and Standard Error of the indirect effects is:
$$
\text{Var}(\hat{\Delta}\_{A/B}) = \text{Var}(\hat{\Delta}\_{A/C}) + \text{Var}(\hat{\Delta}\_{B/C}) \\
\text{SE}(\hat{\Delta}\_{A/B}) = \sqrt{\text{Var}(\hat{\Delta}\_{A/C}) + \text{Var}(\hat{\Delta}\_{B/C})}.
$$

Exponentiating to return to the HR scale and Confidence Intervals:
$$
\text{HR}\_{A/B} = \exp(\hat{\Delta}\_{A/B}) \\
\text{CI}\_{\text{lower}} = \exp(\hat{\Delta}\_{A/B} - 1.96 \cdot \text{SE}(\hat{\Delta}\_{A/B})) \\
\text{CI}\_{\text{upper}} = \exp(\hat{\Delta}\_{A/B} + 1.96 \cdot \text{SE}(\hat{\Delta}\_{A/B})).
$$

### A.3 Matching-Adjusted Indirect Comparison (MAIC)

The comparison using the Bucher method is only valid for a target
population where treatment effects are homogeneous, or where all
potential effect modifiers (for A vs. C and B vs. C) are equally
distributed across studies. Specifically, if covariate *X* modifies the
effect of Treatment A, causing treatment effect heterogeneity, and the
distribution of *X* in Study 1 differs from the distribution of that
factor in Study 2, we typically have
*Δ̂*<sub>*A*/*C*</sub><sup>(1)</sup> ≠ *Δ̂*<sub>*A*/*C*</sub><sup>(2)</sup>.
The difference between *Δ̂*<sub>*A*/*C*</sub><sup>(1)</sup> and
*Δ̂*<sub>*A*/*C*</sub><sup>(2)</sup> is influenced by the treatment
effect heterogeneity induced by the effect modifiers and the imbalance
(difference in means) of these effect modifiers across the study
populations. This heterogeneity is the primary motivation for using
Matching-Adjusted Indirect Comparison (MAIC) in the anchored setting.

In MAIC, IPD from the index study are weighted so that the moments of
selected covariates are balanced with respect to the published moments
of the competitor study. The weight *w*<sub>*i*</sub> for each
participant *i* in the index trial is estimated using a logistic
regression:

$$
\ln(w_i) = \ln\[w(z_i)\] = \ln\left\[\frac{\Pr(S = 2 \mid z_i)}{1 - \Pr(S = 2 \mid z_i)}\right\] = \alpha_0 + z_i \alpha_1,
$$

where *α*<sub>0</sub> is the model intercept and *α*<sub>1</sub> is a
vector of model coefficients. While most applications of weighting, such
as controlling for confounding in observational studies, construct
‘inverse probability’ weights for treatment assignment, MAIC uses ‘odds
weighting’ to model trial assignment. The weight *w*<sub>*i*</sub>
represents the conditional odds that an individual *i* with covariates
*z*<sub>*i*</sub>, selected as marginal effect modifiers, is enrolled in
the competitor study. Alternatively, the weight represents the inverse
conditional odds that the individual is enrolled in the index study.

The logistic regression parameters cannot be derived using conventional
methods such as maximum-likelihood estimation, due to unavailable IPD
for the competitor trial. Signorovitch et al. propose using a method of
moments instead to enforce covariate balance across studies. Prior to
balancing, the IPD covariates are centered on the means or proportions
published for the competitor trial. That is, we must find a solution to:

$$
\bar{x}\_{\text{agg}} \sum\_{i=1}^{n} \exp(\alpha_0 + z_i \alpha_1) = \sum\_{i=1}^{n} z_i \cdot \exp(\alpha_0 + z_i \alpha_1),
$$

*x̄*<sub>agg</sub> denotes a vector of means or proportions for the
covariates in Study 2. Without loss of generality, it can be assumed
that *x̄*<sub>agg</sub> = 0 (e.g we could transform baseline
characteristics in both trials by subtracting *x̄*<sub>agg</sub>) leaving
the equation

$$
0 = \sum\_{i=1}^{n} z_i \cdot \exp(\alpha_0 + z_i \alpha_1),
$$

The right hand side of this weight estimation is the first derivative of
the objective function:

$$
Q(\alpha_1)= \sum\_{i=1}^{n} \exp(\alpha_0 + z_i \alpha_1).
$$

The function *Q*(*α*<sub>1</sub>) is convex and can be minimized using
standard convex optimization algorithms. Provided that there is adequate
overlap, minimization yields the unique finite solution:
$\hat{\alpha_1}=\text{argmin}\[Q(\alpha_1)\]$. Feasible solutions do not
exist if all the values observed for a covariate in *z* are greater or
lesser than its corresponding element in *x̄*<sub>agg</sub>.

#### A.3.1. Model assumption

We briefly describe the assumptions required by MAIC and their
implications:

-   *Internal validity* of the effect estimates is ensured for both the
    index and competitor studies. This is certainly feasible in the case
    of RCTs, where randomization ensures exchangeability over treatment
    assignment on expectation. While internal validity is typically
    maintained in RCTs, it poses a more stringent condition for
    observational studies. It is assumed that there is no informative
    measurement error, missing data, or non-adherence, among other
    factors.
-   *Consistency under parallel studies* requires that there is only one
    well-defined version of each treatment, or that any variations in
    the versions of treatment are irrelevant.
-   *Conditional transportability (exchangeability)* of the marginal
    treatment effect for A versus C from the index to the competitor
    study. This means that trial assignment does not influence this
    measure when conditioned on *z*. Previous research has referred to
    this assumption as the conditional constancy of relative effects. It
    is plausible if *z* includes all covariates considered to modify the
    marginal treatment effect for A versus C, meaning there are no
    unmeasured effect modifiers.
-   *Sufficient overlap* requires that the ranges of the selected
    covariates in Study 1 encompass their respective moments in Study 2.
    Overlap violations can be either deterministic or random.
    Deterministic violations occur structurally due to non-overlapping
    trial target populations (eligibility criteria). Random violations
    occur empirically by chance, particularly when sample sizes are
    small. Consequently, overlap can be evaluated based on absolute
    sample sizes, with the Effective Sample Size (ESS) serving as a
    convenient one-number diagnostic.
-   *Correct specification* of the Study 2 covariate distribution is
    essential. Analysts can only approximate the joint distribution due
    to the unavailability of individual participant data (IPD) for the
    competitor study. Covariate correlations are seldom published for
    Study 2, making it impossible for MAIC to balance them. In such
    cases, it is assumed that these correlations are equal to those in
    the pseudo-sample created by weighting the IPD.

#### A.3.2 Effective sample size

The Effective Sample Size (ESS) is a concept used in various statistical
and methodological contexts to represent the number of independent
observations in a dataset, particularly after adjustments such as
reweighting or clustering. ESS is especially relevant in situations
where the original sample has been modified in some way that affects the
statistical independence or weight of the observations. When
observations are reweighted, the ESS can be calculated as follows:

$$
ESS=\\\frac{{(\sum\_{i=1}^{n}w_i)}^2}{\sum\_{i=1}^{n} w_i^2},
$$

where *n* is the original sample size and *w*<sub>*i*</sub> are the
weights assigned to each observation *i*.

### A.4 Effect modification

Effect modification occurs when the effect of a treatment on an outcome
varies according to the level of another variable, known as the effect
modifier. In the context of randomized controlled trials (RCTs) with
time-to-event (TTE) outcomes, effect modification implies that the
treatment effect on the time until the event of interest (e.g., survival
time, time to disease recurrence) differs across subgroups defined by
the effect modifier. An effect modifier is a variable that alters the
magnitude or direction of the treatment effect. Examples include age,
gender, baseline health status, genetic markers, and other demographic
or clinical characteristics.

## B. Data generating mechanism

We simulate two populations, each consisting of 1,000 subjects, for
Study 1 (Treatment A) and Study 2 (Treatment B), with 500 subjects
assigned to each treatment group. Using the overall covariate means and
ignoring sampling variability, we apply certain parametric assumptions
about the covariate distributions to simulate these populations. Since
the studies are assumed to be appropriately randomized, it makes no
difference whether we simulate the treatment arms jointly or separately,
using arm-specific covariate means. We assume that the covariates are
uncorrelated.

<span style="color:red"> **Here we could experiment around to explore
the imbalance of EM (Sex) between two studies.** </span>

``` r
source("//ktchpc/home/xuzhicha/Rscript/SimFunction.R")
set.seed(2025) # set random seed for reproducibility
N_1 <- 1000 # number of simulated subjects per population
N_2 <- 1000 # number of simulated subjects per population
Treatment_A <- rep(c(1, 0), each = N_1/2) # assume a 1:1 treatment allocation ratio in study 1
Treatment_B <- rep(c(1, 0), each = N_2/2) # assume a 1:1 treatment allocation ratio in study 2

# Covariates for population 1 (Trt A)
Age_1 <- rnorm(n = N_1, mean = 70, sd = 20) 
SBP_1 <- rnorm(n = N_1, mean = 130, sd = 20) 
Sex_1 <- rbinom(n = N_1, size = 1, prob = 0.5) 
Smoke_1 <- rbinom(n = N_1, size = 1, prob = 0.3) 

# Covariates for population 2 (Trt A)
Age_2 <- rnorm(n = N_2, mean = 65, sd = 10) 
SBP_2 <- rnorm(n = N_2, mean = 125, sd = 10) 
Sex_2 <- rbinom(n = N_2, size = 1, prob = 0.3) 
Smoke_2 <- rbinom(n = N_2, size = 1, prob = 0.2) 
```

For each patient, whatever the population, occurrence of the event of
interest is governed by a proportional hazard model with a constant
hazard:

*h* = *h*<sub>0</sub> ⋅ *e**x**p*(\[*β*<sub>*A*</sub>⋅Arm<sub>*A*</sub>+*β*<sub>*B*</sub>⋅Arm<sub>*B*</sub>\]+*β*<sub>1</sub>×Age+*β*<sub>2</sub>×SBP+*β*<sub>3</sub>×Sex+*β*<sub>4</sub>×Smoke+*β*<sub>*i**n**t*</sub>×Sex×(Arm<sub>*A*</sub>+Arm<sub>*B*</sub>)),

where *h*<sub>0</sub> denotes the baseline hazard function, assumed
constant; and *β*<sub>1</sub>, *β*<sub>2</sub>, *β*<sub>3</sub> and
*β*<sub>4</sub> represent the conditional effects of age, systolic blood
pressure, sex, and smoke, respectively. The true conditional treatment
effect (log hazard ratio) for A vs C (in both Study 1 and Study 1) is
*β*<sub>*A*</sub> = *l**o**g*(0.4) and the true conditional treatment
effect for B vs C in both study populations is
*β*<sub>*A*</sub> = *l**o**g*(0.8). Note that, due to the
noncollapsibility of the (log) hazard ratio, these conditional effects
are specific to the adjustment set of covariates used in the
outcome-generating process. Latent event times for both studies, without
accounting for censoring, follow an Gompertz distribution. We set the
shape parameter *α* and scale parameter *β*. Use the data-generating
process of Bender to simulate Gompertz-distributed survival times under
proportional hazards. We could also control the right-censoring rate
using `Censor_rate`.

<span style="color:red"> **Here we could experiment around to explore
the strength of EM (Sex) between two studies.** </span>

``` r
# TTE Outcome parameters setting
beta_A <- log(0.3) # conditional treatment effect (log HR) for A vs C in study 1 
beta_B <- log(0.5) # conditional treatment effect (log HR) for B vs C in study 2 
beta_1 <- 0.005 # conditional effect 1
beta_2 <- 0.03 # EM1
beta_3 <- 1.1 # conditional effect 3
beta_4 <- 1.2 # conditional effect 4
beta_int <- 0.3 # conditional effect 4

alpha <- 0.01
beta <- 0.01
Censor_rate <- 0.1 # rate of censoring distribution as set by AGS
```

Simulate survival times and event indicators under in Study 1 and Study
2.

``` r
# linear predictor 
LP_1 <- beta_A*Treatment_A + beta_1*Age_1 + beta_2*SBP_1 + beta_3*Sex_1 + beta_4*Smoke_1 + beta_int*Sex_1*Treatment_A  
SimDat_1 <- surv.sim_gompertz(N = N_1, LP = LP_1, alpha = alpha, beta = beta, Censor_rate = Censor_rate, Censor_type = "Exp") 
time_1 <- as.numeric(SimDat_1[, 1])
status_1 <- SimDat_1[, 2]

# simulate survival times and event indicators under A and C in S=1
LP_2 <- beta_B*Treatment_B + beta_1*Age_2 + beta_2*SBP_2 + beta_3*Sex_2 + beta_4*Smoke_2 + beta_int*Sex_2*Treatment_B 
SimDat_2 <- surv.sim_gompertz(N = N_2, LP = LP_2, alpha = alpha, beta = beta, Censor_rate = Censor_rate, Censor_type = "Exp") 
time_2 <- as.numeric(SimDat_2[, 1])
status_2 <- SimDat_2[, 2]

# Combine data into data frames
Study1_dat <- data.frame(SUBID = 1:N_1, 
                         Treatment = Treatment_A, 
                         Age = Age_1, SBP = SBP_1, 
                         Sex = Sex_1, Smoke = Smoke_1) %>%
  bind_cols(SimDat_1) %>%
  mutate(Center_Age = Age - mean(Age_2),
         Center_SBP = SBP - mean(SBP_2), 
         Center_Sex = Sex - mean(Sex_2), 
         Center_Smoke = Smoke - mean(Smoke_2))

Study2_dat <- data.frame(SUBID = (N_1 + 1):(N_1 + N_2), 
                         Treatment = Treatment_B, 
                         Age = Age_2, SBP = SBP_2, 
                         Sex = Sex_2, Smoke = Smoke_2) %>%
  bind_cols(SimDat_2)

# Variables to summarize for each study
vars <- c("Treatment", "Age", "SBP", "Sex", "Smoke", "time", "status")

# Summary statistics for Study 1
SummaryStats1 <- process_summary(Study1_dat, vars)
SummaryStats2 <- process_summary(Study2_dat, vars)

full_join(SummaryStats1, SummaryStats2, by = "Variable", suffix = c(".Study1", ".Study2")) %>%
  transmute(Variables = gsub("_.*", "", Variable), 
            Study1 = Statistics.Study1, 
            Study2 = Statistics.Study2)
```

    ## # A tibble: 7 × 3
    ##   Variables Study1          Study2         
    ##   <chr>     <chr>           <chr>          
    ## 1 Treatment 0.5 (0.5)       0.5 (0.5)      
    ## 2 Age       70.044 (19.88)  64.829 (9.94)  
    ## 3 SBP       129.144 (20.56) 124.595 (10.36)
    ## 4 Sex       0.509 (0.5)     0.311 (0.46)   
    ## 5 Smoke     0.324 (0.47)    0.194 (0.4)    
    ## 6 time      1.185 (1.81)    1.335 (1.64)   
    ## 7 status    0.869 (0.34)    0.856 (0.35)

### B.1 True marginal and conditional treatment effects

We calculate the true value of the marginal or population-average
treatment effect *Δ̂*<sub>*A*/*C*</sub><sup>(1)</sup> in the Study 1
population, and the true value of the marginal treatment effect
*Δ̂*<sub>*B*/*C*</sub><sup>(2)</sup> in the Study 2 population. The
marginal effect represents the expected difference in potential
outcomes, measured on the log hazard ratio scale, between a scenario
where all members of a given population receive the active treatment and
a scenario where all members of the population receive the common
comparator.

``` r
# Define the formulae
marginal_formula <- Surv(time, status) ~ Treatment
conditional_formula <- Surv(time, status) ~ Treatment + Age + SBP + Sex + Smoke + Sex * Treatment

# Fit the Cox proportional hazards models using as.formula
Marginal_AC_Model <- coxph(as.formula(marginal_formula), data = Study1_dat)
Marginal_BC_Model <- coxph(as.formula(marginal_formula), data = Study2_dat)

Cond_AC_Model <- coxph(as.formula(conditional_formula), data = Study1_dat)
Cond_BC_Model <- coxph(as.formula(conditional_formula), data = Study2_dat)

# Extract coefficients and calculate exponentiated values
coef_AC <- exp(coef(Marginal_AC_Model))
coef_BC <- exp(coef(Marginal_BC_Model))
coef_Cond_AC <- exp(coef(Cond_AC_Model))
coef_Cond_BC <- exp(coef(Cond_BC_Model))

# Create the summary dataframe
EffectDisplay <- data.frame(
  Coef = c("True_TreatmentA", "True_TreatmentB", "Age", "SBP", "Sex", "Smoke", "Sex_Trt"),
  Value = exp(c(beta_A, beta_B, beta_1, beta_2, beta_3, beta_4, beta_int)),
  Marginal_AC = c(coef_AC, rep(NA, 6)),
  Marginal_BC = c(NA, coef_BC, rep(NA, 5)),
  Conditional_AC = c(coef_Cond_AC[1], NA, coef_Cond_AC[2:6]),
  Conditional_BC = c(NA, coef_Cond_BC)
)

# Remove row names
rownames(EffectDisplay) <- NULL

# View the summary dataframe
print(EffectDisplay)
```

    ##              Coef    Value Marginal_AC Marginal_BC Conditional_AC
    ## 1 True_TreatmentA 0.300000   0.5345675          NA       0.301229
    ## 2 True_TreatmentB 0.500000          NA   0.6133679             NA
    ## 3             Age 1.005013          NA          NA       1.004629
    ## 4             SBP 1.030455          NA          NA       1.030125
    ## 5             Sex 3.004166          NA          NA       2.855746
    ## 6           Smoke 3.320117          NA          NA       3.487546
    ## 7         Sex_Trt 1.349859          NA          NA       1.465432
    ##   Conditional_BC
    ## 1             NA
    ## 2      0.5450184
    ## 3      1.0050435
    ## 4      1.0258524
    ## 5      2.8996974
    ## 6      3.3349322
    ## 7      1.2078218

## C. Bucher method

``` r
# The ratio of true HRs
True_HR <- exp(beta_A) / exp(beta_B)

ITC_Bucher_logHR <- coef(Marginal_AC_Model) - coef(Marginal_BC_Model)
ITC_Bucher_est <- exp(ITC_Bucher_logHR)
ITC_Bucher_bias <- as.numeric(ITC_Bucher_est - True_HR)
ITC_Bucher_se <- sqrt(coef(summary(Marginal_AC_Model))["Treatment", "se(coef)"]^2 +
                        coef(summary(Marginal_BC_Model))["Treatment", "se(coef)"]^2)
# # Compute the 95% confidence interval for the Bucher ITC HR ratio
ITC_Bucher_lowCI <- ITC_Bucher_logHR - qnorm(0.975) * ITC_Bucher_se
ITC_Bucher_upCI <- ITC_Bucher_logHR + qnorm(0.975) * ITC_Bucher_se

# # Print the 95% confidence interval for the Bucher ITC HR ratio
cat("95% CI for HR ratio (Bucher ITC):", round(ITC_Bucher_est, 4), 
    "[", round(exp(ITC_Bucher_lowCI), 4), ",", round(exp(ITC_Bucher_upCI), 4), "]\n")
```

    ## 95% CI for HR ratio (Bucher ITC): 0.8715 [ 0.7184 , 1.0574 ]

## D. MAIC

First, we try to match all the covariates, which are Age, SBP, Sex, and
Smoke. Then we run the diagnostic summary statistics for weights and get
the effective sample size (ESS).

``` r
cent_match_cov <- Study1_dat %>% 
  dplyr::select(Center_Age:Center_Smoke) %>% 
  colnames()
match_cov <- Study1_dat %>% 
  dplyr::select(Age:Smoke) %>% 
  colnames()

est_weights <- MAIC::estimate_weights(intervention_data = Study1_dat,
                                      matching_vars = cent_match_cov)
Diagnostics_1 <- wt_diagnostics(est_weights$analysis_data, vars = match_cov)
cat("The ESS is ", Diagnostics_1$ESS, "\n")
```

    ## The ESS is  707.3544

``` r
cat("The summary statistics for (rescaled) weights are: \n")
```

    ## The summary statistics for (rescaled) weights are:

``` r
Diagnostics_1$Summary_of_weights
```

    ##               type     mean        sd    median        min      max
    ## 1          Weights 0.831454 0.5350667 0.6895267 0.09973892 3.612190
    ## 2 Rescaled weights 1.000000 0.6435314 0.8293023 0.11995724 4.344425

``` r
weights_hist(est_weights$analysis_data, binwidth = 0.05)
```

![](ITC_Week2_MAIC_files/figure-markdown_github/unnamed-chunk-7-1.png)

Then we compare the baseline covariates before/after matching to check
the results of matching.

``` r
# Compute the weighted mean and standard deviation in one step
weighted_summary <- est_weights$analysis_data %>%
  summarise(across(c(Age, SBP, Sex, Smoke, Treatment, time, status), 
                   list(mean = ~ weighted.mean(., wt), 
                        sd = ~ weighted_sd(., wt))))

# Transform the summary to the desired format
SummaryStats3 <- weighted_summary %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "value") %>%
  separate(Variable, into = c("Variable", "Statistic"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = "Statistic", values_from = "value") %>%
  mutate(Statistics = paste0(round(mean, 3), " (", round(sd, 2), ")")) %>%
  select(Variable, Statistics)

# Combine the summaries
combined_summary <- full_join(SummaryStats1, SummaryStats2, by = "Variable", suffix = c(".Study1", ".Study2")) %>%
  full_join(SummaryStats3, by = "Variable") %>%
  transmute(Variables = gsub("_.*", "", Variable), 
            Study1 = Statistics.Study1, 
            Study2 = Statistics.Study2, 
            Study1_adj = Statistics)
print(combined_summary)
```

    ## # A tibble: 7 × 4
    ##   Variables Study1          Study2          Study1_adj     
    ##   <chr>     <chr>           <chr>           <chr>          
    ## 1 Treatment 0.5 (0.5)       0.5 (0.5)       0.49 (0.5)     
    ## 2 Age       70.044 (19.88)  64.829 (9.94)   64.829 (19.62) 
    ## 3 SBP       129.144 (20.56) 124.595 (10.36) 124.595 (20.83)
    ## 4 Sex       0.509 (0.5)     0.311 (0.46)    0.311 (0.46)   
    ## 5 Smoke     0.324 (0.47)    0.194 (0.4)     0.194 (0.4)    
    ## 6 time      1.185 (1.81)    1.335 (1.64)    1.691 (2.29)   
    ## 7 status    0.869 (0.34)    0.856 (0.35)    0.822 (0.38)

Marginal AC in Study 2

``` r
# fit weighted Cox proportional hazards model, robust=TRUE for robust sandwich variance 
ITC_MAIC <- coxph(Surv(time, status) ~ Treatment, 
                data = est_weights$analysis_data, 
                weights = wt, 
                robust = TRUE)

ITC_MAIC_est <- exp(coef(ITC_MAIC) - coef(Marginal_BC_Model))

# Create a dataframe for displaying the effects
EffectDisplay2 <- data.frame(
  Coef = c("True_TreatmentA", "True_TreatmentB", "Age", "SBP", "Sex", "Smoke", "Sex*Trt"),
  Value = exp(c(beta_A, beta_B, beta_1, beta_2, beta_3, beta_4, beta_int)),
  Marginal_AC = c(exp(coef(Marginal_AC_Model)), rep(NA, 6)),
  Marginal_AC_match = c(exp(coef(ITC_MAIC)), rep(NA, 6)),
  Marginal_BC = c(NA, exp(coef(Marginal_BC_Model)), rep(NA, 5)),
  Conditional_AC = c(exp(coef(Cond_AC_Model))[1], NA, exp(coef(Cond_AC_Model))[2:6]),
  Conditional_BC = c(NA, exp(coef(Cond_BC_Model)))
)

# Remove row names
rownames(EffectDisplay2) <- NULL

# View the summary dataframe
print(EffectDisplay2)
```

    ##              Coef    Value Marginal_AC Marginal_AC_match Marginal_BC
    ## 1 True_TreatmentA 0.300000   0.5345675         0.5207169          NA
    ## 2 True_TreatmentB 0.500000          NA                NA   0.6133679
    ## 3             Age 1.005013          NA                NA          NA
    ## 4             SBP 1.030455          NA                NA          NA
    ## 5             Sex 3.004166          NA                NA          NA
    ## 6           Smoke 3.320117          NA                NA          NA
    ## 7         Sex*Trt 1.349859          NA                NA          NA
    ##   Conditional_AC Conditional_BC
    ## 1       0.301229             NA
    ## 2             NA      0.5450184
    ## 3       1.004629      1.0050435
    ## 4       1.030125      1.0258524
    ## 5       2.855746      2.8996974
    ## 6       3.487546      3.3349322
    ## 7       1.465432      1.2078218

### Calculate the Standard Error

#### Robust Sandwich Variance

``` r
# ---- Sandwich Variance Estimation ----
library(clubSandwich)
```

    ## Registered S3 method overwritten by 'clubSandwich':
    ##   method    from    
    ##   bread.mlm sandwich

``` r
vmod <- clubSandwich::vcovCR(ITC_MAIC, cluster=est_weights$analysis_data$SUBID, type="CR2")
coef_res <- clubSandwich::conf_int(ITC_MAIC, vmod, coefs = "Treatment")

# MAIC method calculations
ITC_MAIC_logHR <- coef_res$beta - coef(Marginal_BC_Model)
ITC_MAIC_est <- as.numeric(exp(ITC_MAIC_logHR))
ITC_Bucher_bias <- as.numeric(ITC_MAIC_est - True_HR)

# Calculate standard errors
se_AC <- coef_res$SE
se_BC <- coef(summary(Marginal_BC_Model))["Treatment", "se(coef)"]
ITC_MAIC_se <- sqrt(se_AC^2 + se_BC^2)

# Compute the 95% confidence interval for the Bucher ITC HR ratio
ITC_MAIC_lowCI <- ITC_MAIC_logHR - qnorm(0.975) * ITC_MAIC_se
ITC_MAIC_upCI <- ITC_MAIC_logHR + qnorm(0.975) * ITC_MAIC_se

# Print the 95% confidence interval for the Bucher ITC HR ratio
cat("95% CI for HR ratio (Bucher ITC):", round(ITC_Bucher_est, 4), 
    "[", round(exp(ITC_Bucher_lowCI), 4), ",", round(exp(ITC_Bucher_upCI), 4), "]\n")
```

    ## 95% CI for HR ratio (Bucher ITC): 0.8715 [ 0.7184 , 1.0574 ]

``` r
cat("95% CI for HR ratio (MAIC ITC) using Sandwich Variance:", round(ITC_MAIC_est, 4), 
    "[", round(exp(ITC_MAIC_lowCI), 4), ",", round(exp(ITC_MAIC_upCI), 4), "]\n")
```

    ## 95% CI for HR ratio (MAIC ITC) using Sandwich Variance: 0.8489 [ 0.6035 , 1.1943 ]

#### Bootstrap-based SE

``` r
library(boot)
```

    ## 
    ## Attaching package: 'boot'

    ## The following object is masked from 'package:survival':
    ## 
    ##     aml

``` r
boot_results <- boot(data = Study1_dat, 
                     statistic = bootstrap_maic, 
                     R = 1000, 
                     comp_model = Marginal_BC_Model, 
                     cent_match_cov = cent_match_cov, 
                     match_cov = match_cov)
```

``` r
cat("95% CI for HR ratio (MAIC ITC) using Bootstrap:", round(median(boot_results$t), 4), 
    "[", round(quantile(boot_results$t, probs = 0.025), 4),",", round(quantile(boot_results$t, probs = 0.975), 4), "]\n")
```

    ## 95% CI for HR ratio (MAIC ITC) using Bootstrap: 0.8493 [ 0.6942 , 1.0163 ]

``` r
# boot.ci(boot.out = boot_results, index=1, type="perc")
as.data.frame(boot_results$t) %>% 
  rename(Sample = V1) %>% 
  ggplot(aes(x = Sample)) +
  geom_histogram(fill = "#00877c", color = "black", bins = 100)+
  theme_minimal(base_size = 16) +
  labs(
      x = "Boostrapped HR",
      y = "Frequency"
    )+
  geom_vline(xintercept = quantile(boot_results$t, probs = c(0.025, 0.975)),
             linewidth = 1, linetype = "dashed")+
  geom_vline(xintercept = True_HR,
             linewidth = 1.2, color = "#d47b22")+
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 16, face = "bold"),
      strip.background = element_rect(fill = "grey80", color = "black"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(linewidth = 0.5),
      panel.grid.minor = element_line(linewidth = 0.25),
      # plot.background = element_rect(color = "black", linewidth = 1.5)
    )
```

![](ITC_Week2_MAIC_files/figure-markdown_github/unnamed-chunk-12-1.png)

### K-M Plot

``` r
library(RColorBrewer)
library(survminer)
```

    ## Loading required package: ggpubr

    ## 
    ## Attaching package: 'survminer'

    ## The following object is masked from 'package:survival':
    ## 
    ##     myeloma

``` r
# Unweighted intervention data
KM_int <- survfit(formula = Surv(time, status) ~ Treatment,
                  data = Study1_dat,
                  type="kaplan-meier")

# Weighted intervention data
KM_int_wtd <- survfit(formula = Surv(time, status) ~ Treatment,
                      data = est_weights$analysis_data,
                      weights = wt,
                      type="kaplan-meier")
# Comparator data
KM_comp <- survfit(formula = Surv(time, status) ~ Treatment,
                   data = Study2_dat,
                   type="kaplan-meier")


# Extract survival summary statistics for each KM object
extract_summary <- function(survfit_obj, group_name) {
  summary_obj <- summary(survfit_obj)
  data.frame(
    time = summary_obj$time,
    surv = summary_obj$surv,
    strata = summary_obj$strata,
    group = group_name
  )
}

combined_data <- bind_rows(
  extract_summary(KM_int, "Study 1"),
  extract_summary(KM_int_wtd, "Study 1 Adjusted"),
  extract_summary(KM_comp, "Study 2")
)

# Determine linetype based on treatment status
combined_data <- combined_data %>%
  mutate(
    Treatment = ifelse(grepl("=1", strata), "Treatment", "Control")
  )
```

``` r
# Custom colors and linetypes
custom_linetypes <- c("Treatment" = "solid", "Control" = "dashed")

# Generate the combined plot using ggplot2
ggplot(combined_data, aes(x = time, y = surv, color = group, linetype = Treatment)) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Dark2")) +
  scale_linetype_manual(values = custom_linetypes) +
  labs(
    x = "Time",
    y = "Survival Probability",
    color = "Study",
    linetype = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.width = unit(1, "cm"),
    legend.box = "vertical", 
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(linewidth = 0.5),
    panel.grid.minor = element_line(linewidth = 0.25)
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = T),
    linetype = guide_legend(nrow = 1, byrow = T)
  )
```

![](ITC_Week2_MAIC_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
library(ggpubr)

KM_int <- survfit(formula = Surv(time, status) ~ Treatment + Sex,
                  data = Study1_dat,
                  type="kaplan-meier")

KM1 <- extract_summary(KM_int, "Study 1") %>%
  mutate(
    Treatment = ifelse(grepl("Treatment=1", strata), "Treatment", "Control"),
    Sex = ifelse(grepl("Sex=1", strata), "Male", "Female")
  ) %>% 
  ggplot(aes(x = time, y = surv, color = Sex, linetype = Treatment)) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Dark2")) +
  scale_linetype_manual(values = custom_linetypes) +
  labs(
    x = "Time",
    y = "Survival Probability",
    color = "Sex",
    linetype = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.width = unit(1, "cm"),
    legend.box = "vertical", 
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(linewidth = 0.5),
    panel.grid.minor = element_line(linewidth = 0.25)
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = T),
    linetype = guide_legend(nrow = 1, byrow = T)
  )


# Weighted intervention data
KM_int_wtd <- survfit(formula = Surv(time, status) ~ Treatment + Sex,
                      data = est_weights$analysis_data,
                      weights = wt,
                      type="kaplan-meier")

KM2 <- extract_summary(KM_int_wtd, "Study 1 Match") %>%
  mutate(
    Treatment = ifelse(grepl("Treatment=1", strata), "Treatment", "Control"),
    Sex = ifelse(grepl("Sex=1", strata), "Male", "Female")
  ) %>% 
  ggplot(aes(x = time, y = surv, color = Sex, linetype = Treatment)) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Dark2")) +
  scale_linetype_manual(values = custom_linetypes) +
  labs(
    x = "Time",
    y = "Survival Probability",
    color = "Sex",
    linetype = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.width = unit(1, "cm"),
    legend.box = "vertical", 
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(linewidth = 0.5),
    panel.grid.minor = element_line(linewidth = 0.25)
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = T),
    linetype = guide_legend(nrow = 1, byrow = T)
  )

# Comparator data
KM_comp <- survfit(formula = Surv(time, status) ~ Treatment + Sex,
                   data = Study2_dat,
                   type="kaplan-meier")

KM3 <- extract_summary(KM_comp, "Study 2") %>%
  mutate(
    Treatment = ifelse(grepl("Treatment=1", strata), "Treatment", "Control"),
    Sex = ifelse(grepl("Sex=1", strata), "Male", "Female")
  ) %>% 
  ggplot(aes(x = time, y = surv, color = Sex, linetype = Treatment)) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Dark2")) +
  scale_linetype_manual(values = custom_linetypes) +
  labs(
    x = "Time",
    y = "Survival Probability",
    color = "Sex",
    linetype = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.width = unit(1, "cm"),
    legend.box = "vertical", 
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(linewidth = 0.5),
    panel.grid.minor = element_line(linewidth = 0.25)
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = T),
    linetype = guide_legend(nrow = 1, byrow = T)
  )

ggpubr::ggarrange(KM1, KM2, KM3, common.legend = T, nrow = 1)
```

![](ITC_Week2_MAIC_files/figure-markdown_github/unnamed-chunk-15-1.png)

### The results varying strength and balance

#### Read in all the data

``` r
combined_data %>% # Calculate the mean of each column
  pivot_longer(cols = c("Bucher_bias", "MAIC_bias", "MAIC_boot_bias"), 
               names_to = "MethodBias", values_to = "Bias") %>%
  dplyr::filter(N_2 == 1000, True_HR == 0.5, beta_2 < 0.1) %>% 
  mutate(beta_int = paste0("Beta_int = ", round(beta_int, 2))) %>% 
  ggplot(aes(x = as.factor(Sex.perc.2), y = Bias, fill = MethodBias)) + 
  geom_boxplot() +
  facet_wrap(~ beta_int, nrow = 1)+
  scale_fill_manual(values = brewer.pal(n = 3, name = "Dark2")) +
  theme_minimal(base_size = 14, base_line_size = 0.5)+
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.width = unit(0.5, "cm"),
    legend.box = "vertical", 
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1),
    panel.grid.major = element_line(linewidth = 0.3),
    panel.grid.minor = element_line(linewidth = 0.1),
    strip.background = element_rect(fill = "grey90", color = "black"), 
    strip.text = element_text(size = 14, face = "bold")
  )+
  labs(
    x = "Sex Percentage in Study 2",
    y = "Bias",
    fill = "Method")
```

![](ITC_Week2_MAIC_files/figure-markdown_github/unnamed-chunk-17-1.png)
