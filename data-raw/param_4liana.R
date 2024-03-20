require(data.table)
require(MCMCglmm)

# Load Chave et al. 2014 data
data <- fread("data-raw/424_biomass_lianas.csv")

data[, Dry.total.AGB..kg. := Dry.AGB.wood..kg. + Dry.weight.leaf..kg.]

# work on the data
data[, AGBlog := log(Dry.total.AGB..kg.)]
setnames(data, c("Trunk.diameter..cm."), c("D"))
data[, Comp := D ]
data <- data[!is.na(Comp), ]

# remove the trees where the heigth is less than 1.3 m
#data <- data[H >= 1.3, ]



# Model and MCMC
form <- as.formula("AGBlog~I(log(Comp))")
mod <- MCMCglmm(form, data = data, pr = T, nitt = 13001)

# extract the data from the MCMC
mod$VCV <- as.matrix(mod$VCV)
param_4liana <- data.frame(mod$Sol, sqrt(mod$VCV))


colnames(param_4liana) <- c("intercept", "logagbt", "sd")
usethis::use_data(param_4liana)
#BIOMASS::use_data(param_4liana, overwrite = TRUE)
