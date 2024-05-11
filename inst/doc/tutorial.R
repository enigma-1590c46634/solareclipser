## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(solareclipser)

## ----fake_inc_solareclipser, echo = FALSE, message = FALSE--------------------
#library(solareclipser)
#library(devtools)
#load_all("../")

## ----inc_solareclipser, eval = FALSE------------------------------------------
#  library(solareclipser)

## ----ver_solareclipser--------------------------------------------------------
packageVersion("solareclipser")

## ----inc_other----------------------------------------------------------------
library(plyr)

## ----minimal, cache = TRUE----------------------------------------------------
# load library
library(solareclipser)

# load data set
data(dat30)

# univariate polygenic model
mod1 <- solarPolygenic(trait1 ~ 1, dat30)

# bivariate polygenic model
mod2 <- solarPolygenic(trait1 + trait2 ~ 1, dat30,
  polygenic.options = '-testrhoe -testrhog')

# specify directory with IBD matrices and run linkage model
mibddir <- system.file('extdata', 'solarOutput',
  'solarMibdsCsv', package = 'solareclipser') 
link <- solarMultipoint(trait1 ~ 1, dat30,
  mibddir = mibddir, chr = 5)

# run association model in parallel
assoc <- solarAssoc(trait1 ~ 1, dat30, cores = 2,
  snpcovdata = genocovdat30, snpmap = mapdat30)

## ----dat30, cache = TRUE------------------------------------------------------
data(dat30)
str(dat30)
str(genocovdat30)
str(mapdat30)

## ----dat50, cache = TRUE------------------------------------------------------
data(dat50)
str(phenodata)
str(genodata)
str(genocovdata)

## ----multic, cache = TRUE-----------------------------------------------------
mdat <- loadMulticPhen()
str(mdat)

## ----ex, cache = TRUE---------------------------------------------------------
#gawdat <- loadExamplesPhen()
#str(gawdat)

## ----eval = FALSE-------------------------------------------------------------
#  dat.dir <- system.file("extdata", "solarOutput", package = "solareclipser")
#  
#  mdat <- readPhen(phen.file = file.path(dat.dir, "simulated.phen"), sep.phen = ",",
#    ped.file = file.path(dat.dir, "simulated.ped"), sep.ped = ",")

## ----kin_dat30, cache = TRUE--------------------------------------------------
plotKinship2(2*kin)

## ----kin_solar, cache = TRUE--------------------------------------------------
kin2.dat30 <- solarKinship2(dat30)
plotKinship2(kin2.dat30)

## ----sub_kin_plot, cache = TRUE-----------------------------------------------
plotKinship2(kin2.dat30[1:30, 1:30])

## ----hist_kin, cache = TRUE---------------------------------------------------
histKinship2(kin2.dat30)

## ----unique_k2, cahce = TRUE--------------------------------------------------
unique(as.vector(kin2.dat30))

## ----dat30_ped_2, cache = TRUE------------------------------------------------
plotPed(dat30, 2)

## ----exp_trait1---------------------------------------------------------------
dat30 <- mutate(dat30,
  exp_trait1 = exp(trait1))

## -----------------------------------------------------------------------------
availableTransforms()

## ----plot_transform, cache = TRUE---------------------------------------------
library(ggplot2)
library(gridExtra)

theme_set(theme_bw())

with(dat30, {
  grid.arrange(
    qplot(exp_trait1),
    qplot(transformTrait(exp_trait1, "out")),
    qplot(transformTrait(exp_trait1, "log")),
    qplot(transformTrait(exp_trait1, "inormal"))
  )
})

## ----M1, cache = TRUE---------------------------------------------------------
M1 <- solarPolygenic(trait1 ~ age + sex, dat30)
M1

## ----M1_sum, cache = TRUE-----------------------------------------------------
summary(M1)

## ----M2, cache = TRUE---------------------------------------------------------
M2 <- solarPolygenic(trait1 ~ age + sex, dat30, covtest = TRUE)
M2

## ----M2_cf, cache = TRUE------------------------------------------------------
M2$cf

## ----M3, cache = TRUE---------------------------------------------------------
M3 <- solarPolygenic(trait ~ 1, phenodata, kinship = kin)
M3

## ----M3_lf, cache = TRUE------------------------------------------------------
M3$lf

## ----B1, cache = TRUE---------------------------------------------------------
B1 <- solarPolygenic(trait1 + trait2 ~ 1, dat30)
B1

## ----B2, cache = TRUE---------------------------------------------------------
B2 <- solarPolygenic(trait1 + trait2 ~ 1, dat30, polygenic.options = "-testrhoe -testrhog")
B2

## ----B2_vcf, cache = TRUE-----------------------------------------------------
B2$vcf

## ----B3, cache = TRUE---------------------------------------------------------
B3 <- solarPolygenic(trait1 + trait2 ~ sex + age(trait2), dat30)
B3

## ----tail_mod-----------------------------------------------------------------
tail(B3$solar$files$model$poly.mod, 3)


## ----genodata, cache = TRUE---------------------------------------------------
A1 <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata, kinship = kin)
A1

## ----A1_summary---------------------------------------------------------------
summary(A1)

## ----A1_summary_2-------------------------------------------------------------
summary(A1, alpha = 1)

## ----genocovdata, cache = TRUE------------------------------------------------
A2 <- solarAssoc(trait ~ 1, phenodata, snpcovdata = genocovdata, kinship = kin)
A2

## ----genocovdata2, cache = TRUE-----------------------------------------------
genocovdata2 <- genocovdata
genocovdata2[genocovdata2 == 2] <- 1.9

A2 <- solarAssoc(trait ~ 1, phenodata, snpcovdata = genocovdata2, kinship = kin)
A2

## ----genocov, cache = TRUE----------------------------------------------------
dir <- package.file("extdata", "solarAssoc", package = "solareclipser")
genocov.files <- file.path(dir, "snp.genocov")
snplists.files <- file.path(dir, c("snp.geno-list1", "snp.geno-list2"))

A3 <- solarAssoc(trait ~ 1, phenodata, 
  genocov.files = genocov.files, snplists.files = snplists.files)
A3

## ----genocov2, cache = TRUE---------------------------------------------------
dir <- package.file("extdata", "solarAssoc", package = "solareclipser")
genocov.files <- file.path(dir, c("snp.genocov1", "snp.genocov2"))
snplists.files <- file.path(dir, c("snp.geno-list1", "snp.geno-list2"))

A4 <- solarAssoc(trait ~ 1, phenodata, 
  genocov.files = genocov.files, snplists.files = snplists.files)
A4

## ----assoc_map, cache = TRUE--------------------------------------------------
A5 <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata, snpmap = snpdata, kinship = kin)

## ----manh_plot, cache = TRUE--------------------------------------------------
plot(A5)

## ----manh_qq, cache = TRUE----------------------------------------------------
plot(A5, "qq")

## ----A6, cache = TRUE---------------------------------------------------------
A6 <- solarAssoc(trait ~ 1, phenodata, snpdata = genodata, kinship = kin, cores = 2)

## ----A16----------------------------------------------------------------------
A1$assoc$tprofile$cputime.sec
A6$assoc$tprofile$cputime.sec

## ----A7, cache = TRUE---------------------------------------------------------
A7 <- solarAssoc(trait1 ~ 1, dat30, 
  snpcovdata = genocovdat30, snpmap = mapdat30)

## ----A8, cache = TRUE---------------------------------------------------------
A8 <- solarAssoc(trait1 ~ 1, dat30, 
  snpcovdata = genocovdat30, snpmap = mapdat30,
  cores = 2)

## ----A_78---------------------------------------------------------------------
A7$assoc$tprofile$cputime.sec
A8$assoc$tprofile$cputime.sec

## ----data_multic_link---------------------------------------------------------
mibddir <- package.file("extdata", "solarOutput", "solarMibdsCsv", package = "solareclipser")  
list.files(mibddir)

## ----L1, cache = TRUE---------------------------------------------------------
L1 <- solarMultipoint(trait1 ~ 1, dat30, mibddir = mibddir)
L1

## ----L1_sum-------------------------------------------------------------------
summary(L1)

## ----plot_L1, cache = TRUE----------------------------------------------------
plot(L1)

## ----L2, cache = TRUE---------------------------------------------------------
L2 <- solarMultipoint(trait1 ~ 1, dat30, mibddir = mibddir, chr = 5, interval = 2, 
  multipoint.settings = "finemap off", multipoint.options = "3")
summary(L2)

## ----plot_L2, cache = TRUE----------------------------------------------------
plot(L2, pass = 2)

## ----L3, cache = TRUE---------------------------------------------------------
L3 <- solarMultipoint(trait1 + trait2 ~ 1, dat30, mibddir = mibddir, chr = 5, 
  interval = 2, multipoint.settings = "finemap off")
summary(L3)

## ----L4, cache = TRUE---------------------------------------------------------
L4 <- solarMultipoint(trait1 + trait2 ~ age + sex(trait1), dat30, 
  mibddir = mibddir, chr = 5, interval = 2, multipoint.settings = "finemap off")
summary(L4)

## ----L5, cache = TRUE---------------------------------------------------------
L5 <- solarMultipoint(trait1 ~ 1, dat30, mibddir = mibddir,
  chr = c(2, 5), cores = 2)

## ----L1_cpu-------------------------------------------------------------------
L1$multipoint$tprofile$cputime.sec
L5$multipoint$tprofile$cputime.sec

## ----session_info-------------------------------------------------------------
sessionInfo()

