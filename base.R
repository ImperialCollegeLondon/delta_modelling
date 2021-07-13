library(rstan)
library(matrixStats)
library(data.table)
library(EnvStats)
library(lubridate)
library(tidyverse)
library(abind)
library(here)
library(job)
library(rio)
library(optparse)
library(glue)
library(huxtable)
library(posterior)
option_list <- list(
  make_option(c("--sdate"),action="store", default="2020-03-01",help="When to start the model [default \"%default\"]"),
  make_option(c("--edate"),action="store", default="2021-05-25",help="End date for the model [default \"%default\"]"),
  make_option(c("--vdate"),action="store", default="2021-01-31",help="Start of variant date for the model [default \"%default\"]"),
  make_option(c("--jname"),action="store", default="base",help="Name of stanfile [default \"%default\"]"),
  make_option(c("--wpar"),action="store", type="integer", default=310,help="Wanning Rayleigh Par [default \"%default\"]"),
  make_option(c("--areaname"),action="store", default="Mumbai",help="Which area to run this for [default \"%default\"]"),
  make_option(c("--nchains"),action="store", type="integer", default=3,help="Number of Chains [default \"%default\"]"),
  make_option(c("--warmup"),action="store", type="integer", default=500,help="Number of warmup iterations [default \"%default\"]"),
  make_option(c("--iter"),action="store", type="integer", default=1000,help="Number of iterations [default \"%default\"]"),
  make_option(c("--thin"),action="store", type="integer", default=1,help="Amount of thinning of results [default \"%default\"]"),
  make_option(c("--seed"),action="store", type="integer", default=sample(1:1e5, 1),help="Seed for run [default \"%default\"]"),
  make_option(c("--delta"),action="store", type="double", default=0.95,help="Amount of adapt delta for sampling [default \"%default\"]"),
  make_option(c("--tree"),action="store", type="integer", default=12,help="Tree depth for MCMC [default \"%default\"]"),
  make_option(c("--ur"),action="store", type="double", default=1.0,help="Underreporting factor [default \"%default\"]"),
  make_option(c("--ifrMean"),action="store", type="double", default=0.3,help="IFR ~N(ifrMean,)  [default \"%default\"]"),
  make_option(c("--ifrSD"),action="store", type="double", default=0.02,help="IFR ~ N(,ifrSD) [default \"%default\"]")
)
opt <- parse_args(OptionParser(option_list=option_list))

start_date <- ymd(opt$sdate)
end_date <- ymd(opt$edate) 
job <- opt$jname
rayleigh_par <- opt$wpar
rel_IFR1 <- 1/100
T2_date <- ymd(opt$vdate)
ITER <- opt$iter
WARM <- opt$warmup
CORES <- opt$nchains
SEED <- opt$seed
DELTA <- opt$delta
TREE <- opt$tree
RANGE_TIME <- seq.Date(start_date, end_date,by = 'day')
selected_region <-  opt$areaname

if(selected_region=="Mumbai"){
  # get data for mumbai,, dat aexist only from 26th Apr 2020 so starting stuff 
  # from earlier and just inserting -1 to make sure I have seeded priperly
  # likelihood start s only from 1st date
  df<-
    import("https://api.covid19india.org/csv/latest/districts.csv") %>%
    filter(District == selected_region) %>%
    mutate(Date = ymd(Date),
           Deaths = Deceased - lag(Deceased),
           Cases = Confirmed - lag(Confirmed)) %>%
    mutate(Deaths = case_when(Deaths <= 0 ~ -1.0, TRUE ~ as.double(Deaths )),
           Cases = case_when(Cases <= 0 ~ -1.0, TRUE ~ as.double(Cases))) %>%
    select(Date, District, Deaths, Cases) %>%
    complete(Date = seq.Date(start_date, max(.$Date), by = "day"), nesting(District), fill = list(Deaths=-1,Cases=-1)) %>%
    rename(region=District)
}else{
  stop("Other regions not supported as of now")
}
sero_data <- import(here("data","sero_positivity.csv"))

df_pop <- import(here("data", "population.csv"))
serial.interval <- import(here("data","serial_interval.csv"))
onset_paras <-
  import(here("data","statelevel_OnsetDeathParams.csv")) %>%
  select(region,mean,cv)
pcr_and_sero <-import(here("data","pcr_and_sero_pos.csv"))
if(opt$areaname == "Mumbai")
{
  pcr_file = "pcr_genome_fraction_mumbai.csv"
} else {
  print("wrong location")
}

pcr_genome_fraction <-
  import(here("data",pcr_file)) %>%
  mutate(date = ymd(agg_dat),
         negative  = N_sampled - N_positive,
         positive = N_positive) %>%
  select(date, positive, negative) %>%
  filter(date >= T2_date, date <= end_date)
df <-
  df %>%
  left_join(df_pop) %>%
  arrange(Date) %>%
  filter(Date <= ymd(end_date))

#indices for likelihood
ll_idxs <- which(!df$Deaths==-1)
ll_len <- length(ll_idxs)

StanModel <- as.character(job)
d<-df
N2 <- length(RANGE_TIME)
dates <- list()
reported_cases <- list()
stan_data = list(M=length(selected_region),
                 N=NULL,
                 deaths=NULL,
                 f=NULL,
                 N0=6,
                 SI=NULL,
                 ll_len = ll_len,
                 ll_idxs = ll_idxs,
                 pop = NULL,
                 T2 = NULL,
                 phylo_N = NULL,
                 phylo_PSamples = NULL,
                 phylo_NSamples = NULL,
                 sero_N = NULL
)
deaths_by_country <-list()
deaths_by_country_combined <- list()
mean1 <- 5.1; cv1 <- 0.86;
x1 <- rgammaAlt(1e6,mean1,cv1)
mean2 <- onset_paras$mean[onset_paras$region == selected_region]
cv2 <- onset_paras$cv[onset_paras$region == selected_region]
x2 <- rgammaAlt(1e6,mean2,cv2)
ecdf.saved <- ecdf(x1+x2)
stan_data$pop <- c(stan_data$pop, df_pop[df_pop$region == selected_region,]$population)
stan_data$pop <- as.array(stan_data$pop)
stan_data$par <- rayleigh_par
dates[[selected_region]] <- d$Date
N <- length(d$Cases)
convolution1 <- function(u) (rel_IFR1 * ecdf.saved(u))
f <- rep(0,N2)
f[1] <- (convolution1(1.5) - convolution1(0))
for(i in 2:N2)
{
  f[i] <- (convolution1(i+.5) - convolution1(i-.5))
}

cases <- as.vector(as.numeric(d$Cases))
deaths <- as.vector(as.numeric(d$Deaths))
stan_data$N <- c(stan_data$N,N)
stan_data$f <- cbind(stan_data$f,f)
stan_data$deaths <- cbind(stan_data$deaths,deaths)
stan_data$cases <- cbind(stan_data$cases,cases)
stan_data$N2 <- N2

if(length(stan_data$N) == 1) {
  stan_data$N = as.array(stan_data$N)
}
stan_data$W <- ceiling(stan_data$N2/7)
stan_data$week_index <- matrix(1,stan_data$M,stan_data$N2)
for(state.i in 1:stan_data$M) {
  stan_data$week_index[state.i,] <- rep(2:(stan_data$W+1),each=7)[1:stan_data$N2]
  last_ar_week <- which(dates[[state.i]]==max(df$Date) -  7)
  stan_data$week_index[state.i,last_ar_week:ncol(stan_data$week_index)] <-
    stan_data$week_index[state.i,last_ar_week]
}
stan_data$AR_SD_MEAN = 0.2
stan_data$P <- 0
stan_data$SI <- serial.interval$fit[1:N2]
stan_data$T2 <- which(dates[[selected_region]]==T2_date)

sero_data_reg = sero_data[sero_data$region == selected_region,]
stan_data$sero_N  <- which(dates[[selected_region]] %in% sero_data_reg$date)
stan_data$sero_N_len <- length(stan_data$sero_N)
stan_data$sero_prev <- sero_data_reg$percentage

length_pcr_sero_data <- length(pcr_and_sero$p_PCR_positive)
padding <- stan_data$N2 - length_pcr_sero_data
stan_data$PCR_pos_prob <- as.matrix(c(pcr_and_sero$p_PCR_positive, rep(0, padding)))
stan_data$seroconv_cdf <- as.matrix(c(pcr_and_sero$cum_seropositive, rep(1, padding)))
stan_data$serorev_surv <- as.matrix(1 - pweibull(seq(1, N2), shape = 2.933, scale = 208))

stan_data$phylo_N  <- which(dates[[selected_region]] %in% pcr_genome_fraction$date)
stan_data$phylo_N_len <- length(stan_data$phylo_N)
stan_data$phylo_PSamples <- pcr_genome_fraction$positive
stan_data$phylo_NSamples <- pcr_genome_fraction$negative
stan_data$WI =  exp( - 0.5 * (1:stan_data$N2)^2 / (stan_data$par^2) )

stan_data$UR <- opt$ur
stan_data$ifrMean <- opt$ifrMean
stan_data$ifrSD <- opt$ifrSD
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m <- stan_model(here(paste0('stan_models/',StanModel,'.stan')))

fit <- sampling(m,
                data=stan_data,
                iter=ITER,
                warmup=WARM,
                chains=CORES,
                cores=CORES,
                seed=SEED,
                init = 0.01,
                control = list(adapt_delta = DELTA, max_treedepth = TREE))
out <- rstan::extract(fit)

if(!dir.exists(here("figures",job))){
  dir.create(here("figures",job), recursive = TRUE)
}
if(!dir.exists(here("results",job))){
  dir.create(here("results",job), recursive = TRUE)
}

saved_data = list(
  dates = dates,
  stan_data = stan_data,
  fit = fit,
  selected_region = selected_region
)

tr_median = round(100 * mean(out$R_difference - 1))
tr_lr = round(100 * quantile(out$R_difference - 1, probs = 0.25))
tr_ur =  round(100 * quantile(out$R_difference - 1, probs = 0.75))
ie_median = round(100 * mean(1 - out$cross))
ie_lr = round(100 * quantile(1 - out$cross, probs = 0.25))
ie_ur =  round(100 * quantile(1 - out$cross, probs = 0.75))
ifr1_median = signif(mean(out$ifr1), digits = 2)
ifr1_lr = signif( quantile(out$ifr1, probs = 0.25), digits = 2)
ifr1_ur = signif(quantile(out$ifr1, probs = 0.75), digits = 2)
ifr2_median = signif(mean(out$ifr2), digits = 2)
ifr2_lr = signif(quantile(out$ifr2, probs = 0.25), digits = 2)
ifr2_ur = signif(quantile(out$ifr2, probs = 0.75), digits = 2)

print_md(as_hux(tibble("Start Date" = opt$vdate,
              "Under reporting" = stan_data$UR,
              "Immune Escape" = glue("{ie_median}% ({ie_lr}% - {ie_ur}%)"),
              "Transmissibility Increase" = glue("{tr_median}% ({tr_lr}% - {tr_ur}%)"),
              "IFR1" = glue("{ifr1_median}% ({ifr1_lr}% - {ifr1_ur}%)"),
              "IFR2" = glue("{ifr2_median}% ({ifr2_lr}% - {ifr2_ur}%)")
)))


filename <- paste(opt$jname,gsub("-","_",opt$edate),gsub("-","_",opt$vdate),opt$wpar,opt$areaname,opt$iter,opt$ifrMean,opt$ifrSD,opt$ur,sep="_")
saveRDS(saved_data, here(paste0("results/",job,"/",filename,"_saved_data.rds")))



