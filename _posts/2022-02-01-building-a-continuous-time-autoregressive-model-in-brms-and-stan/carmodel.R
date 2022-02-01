
library("tidyverse")
library("brms")
library("rstan")

rel_path <- "_posts/2022-02-01-building-a-continuous-time-autoregressive-model-in-brms-and-stan/"

stan_seed <- 1256

# simulate a set of irregularly sampled AR(1) processes:

phi <- .8

withr::with_seed(stan_seed, {
  data <- tibble(
    x = 1:100,
    y1 = arima.sim(list(ar = phi), length(x)),
    y2 = arima.sim(list(ar = phi), length(x))
  ) %>% 
    pivot_longer(starts_with("y"), names_to = "g", values_to = "y")
  
  subset <- data %>% 
    slice_sample(prop = .6) %>% # change the proportion retained
    arrange(g, x) %>% 
    group_by(g) %>% 
    mutate(
      x_lag = lag(x),
      d_x = x - x_lag # spacing of observations
    ) %>% 
    ungroup()
})

# use brms to generate stan code/data for the model:

# scode <- brms::make_stancode(y ~ ar(time = x), data = subset)
sdata <- brms::make_standata(y ~ ar(time = x, gr = g), data = subset)

# modify stan code/data for CAR(1) model:
# (mods are all commented with the tag "CAR(1)")

sdata$s <- replace_na(subset$d_x, 0) # exponent (replace 1st element w/ 0)
# write_lines(scode, file = "scode.stan") # write brms-generated code, then modify
scode <- read_lines(file = paste0(rel_path, "scode.stan")) # read modified code 

# fit model:

# stanfit <- rstan::stan(
#   model_code = scode, 
#   data = sdata, 
#   cores = 4,
#   sample_file = "models/armodel", # output in csv format
#   seed = stan_seed
# )

# read csv to create stanfit object:

csvfiles <- list.files(
  path = rel_path,
  pattern = "armodel_\\d\\.csv",
  full.names = TRUE
)

stanfit <- rstan::read_stan_csv(csvfiles)

# feed the stan model back into brms (for postprocessing):

fit <- brm(y ~ ar(time = x, gr = g), data = subset, empty = TRUE)
fit$fit <- stanfit
fit <- rename_pars(fit)
summary(fit)

# plot data and model fit:

fitted_vals <- fitted(fit)

subset %>%
  bind_cols(as_tibble(fitted_vals)) %>% 
  ggplot(aes(x)) + 
  facet_wrap(vars(g)) + 
  geom_line(aes(y = y, col = "data")) + 
  geom_line(aes(y = Estimate, col = "model"))

