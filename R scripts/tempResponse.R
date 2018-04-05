#===============================================================================
# Photosynthesis temperature response
#===============================================================================
# Prepared by Sean Michaletz (sean.michaletz@gmail.com), May 2017
# Revised by PFTC3-Peru group 2

#===============================================================================
# Contents
#===============================================================================

#  Part 1:  Photosynthesis temperature response curves (standard)
#  Part 2:  Arrhenius plots
#  Part 3:  Estimating activation energy
#  Part 4:  TPCFitting

#===============================================================================
# Introductory code: Load packages/libraries, set paths, load datasets
#===============================================================================
#--Clear memory
rm(list = ls(all = T))

# << DATASETS >>

# << PACKAGES >>
library(ggplot2)
library(polynom)
library(grid)
library(gridExtra)
library(plyr)
library(scales)
library(reshape2)
#library(nlsLoop)
library(magrittr)
library(tidyr)
library(dplyr)

#--Load photosynthesis data.
psyn <- read.csv("./workingData/dataMaster_Peru_area.csv", header = T) %>%
    filter(Taxon != "Paspallum") # eliminates the empty line

# << GLOBAL OPTIONS >>
#--Define colorblind palette (see http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/).
cbPalette <- c("#000000",
               "#E69F00",
               "#56B4E9",
               "#009E73",
               "#CC79A7",
               "#0072B2",
               "#D55E00",
               "#F0E442",
               "#999999")

# << CLEAN AND MERGE DATA >>
#--Add columns to identify individual temperature response curves in SUMO and RMBL data
psyn <- transform(psyn, curveID = as.integer(factor(Filename, unique(Filename))))

#--Add column for processID
psyn$rate <- "Net photosynthesis (umol m-2 s-1)"


# Calculate Boltzmann leaf temperature
psyn$invBT_eV <- 1 / (0.00008617 * (psyn$Tleaf + 273.15))



#===============================================================================
# Part 1: Photosynthesis temperature response curves (standard)
#=============================================================================

#--Plot all curves on a single set of axes.
ggplot(psyn, aes(x = Tleaf, y = Photo)) +
    stat_smooth(
        method = "lm",
        se = TRUE,
        fill = NA,
        formula = y ~ poly(x, 2, raw = TRUE),
        aes(colour = factor(curveID), linetype = factor(curveID)),
        data = psyn
    ) +
    geom_point(aes(color = as.factor(curveID)), data = psyn) +
    xlab(expression('Leaf temperature' ~ (degree * C))) +
    ylab(expression("Assimilation rate (" * mu ~ "mol" ~ m ^ -2 ~ s ^ -1 * ")")) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    ) +
    scale_colour_manual(values = rep(cbPalette, times = 22)) +
    # plot each 22 times to maximize solids (there are 9 colors)
    scale_linetype_manual(values = rep(c(1:6), each = 33))

#--Plot all curves on a single set of axes (curves only).
ggplot(psyn, aes(x = Tleaf, y = Photo)) + stat_smooth(
    method = "lm",
    se = TRUE,
    fill = NA,
    formula = y ~ poly(x, 2, raw = TRUE),
    aes(colour = factor(curveID)),
    data = psyn
) +
    xlab(expression('Leaf temperature' ~ (degree * C))) +
    ylab(expression("Assimilation rate (" * mu ~ "mol" ~ m ^ -2 ~ s ^ -1 * ")")) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    ) +
    scale_colour_manual(values = rep("gray60", times = 30))

#--Plot all curves on individual axes.
for (i in seq_along(unique(psyn$curveID))) {
    print(ggplot(
            subset(psyn, psyn$curveID == i),
            aes(x = Tleaf, y = Photo)
        ) +
            stat_smooth(
                method = "lm",
                se = TRUE,
                fill = NA,
                formula = y ~ poly(x, 2, raw = TRUE),
                aes(colour = factor(curveID)),
                data = subset(psyn, psyn$curveID == i)
            ) +
            geom_point(aes(color = as.factor(curveID)), data = subset(psyn, psyn$curveID == i)) +
            xlab(expression('Leaf temperature' ~ (degree * C))) +
            ylab(expression(
                "Assimilation rate (" * mu ~ "mol" ~ m ^ -2 ~ s ^ -1 * ")"
            )) +
            theme_bw(base_size = 12) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black")
            )
    )
}

#--Plot a single curve (e.g., curveID = 47).
ggplot(subset(psyn, psyn$curveID == "4"),
       aes(x = Tleaf, y = Photo)) +
    stat_smooth(
        method = "lm",
        se = TRUE,
        fill = NA,
        formula = y ~ poly(x, 2, raw = TRUE),
        aes(colour = factor(curveID))
    ) +
    geom_point(aes(color = as.factor(curveID))) +
    xlab(expression('Leaf temperature' ~ (degree * C))) +
    ylab(expression("Assimilation rate (" * mu ~ "mol" ~ m ^ -2 ~ s ^ -1 * ")")) +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    )

#=========================================================================================
# Part 2: Arrhenius plots
#=========================================================================================

#--Plot all curves on a single set of axes.----
ggplot(psyn, aes(x = invBT_eV, y = Photo)) +
    stat_smooth(
        method = "lm",
        se = TRUE,
        fill = NA,
        formula = y ~ poly(x, 2, raw = TRUE),
        aes(colour = factor(curveID), linetype = factor(curveID)),
        data = psyn
    ) +
    geom_point(aes(color = as.factor(curveID)), data = psyn) +
    xlab(expression(paste(
        'Leaf temperature ', '<1/', italic('kT'), '>', ' (',  eV ^ {
            -1
        }, ')'
    ))) +
    ylab(expression("Assimilation rate (" * mu ~ "mol" ~ m ^ -2 ~ s ^ -1 * ")")) +
    scale_y_continuous(
        trans = "log",
        breaks = trans_breaks("log", function(x)
            exp(x), n = 3),
        labels = trans_format("log", math_format(e ^ .x))
    ) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    ) +
    scale_colour_manual(values = rep(cbPalette, times = 22)) +
    # plot each 22 times to maximize solids (there are 9 colors)
    scale_linetype_manual(values = rep(c(1:6), each = 33))

#--Plot all curves on a single set of axes (curves only).----
ggplot(psyn, aes(x = invBT_eV, y = Photo)) +
    stat_smooth(
        method = "lm",
        se = TRUE,
        fill = NA,
        formula = y ~ poly(x, 2, raw = TRUE),
        aes(colour = factor(curveID)),
        data = psyn
    ) +
    xlab(expression(paste(
        'Leaf temperature ', '<1/', italic('kT'), '>', ' (',  eV ^ {
            -1
        }, ')'
    ))) +
    ylab(expression("Assimilation rate (" * mu ~ "mol" ~ m ^ -2 ~ s ^ -1 * ")")) +
    scale_y_continuous(
        trans = "log",
        breaks = trans_breaks("log", function(x)
            exp(x), n = 3),
        labels = trans_format("log", math_format(e ^ .x))
    ) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    ) +
    scale_colour_manual(values = rep("gray60", times = 197))

#--Plot curves 1 to 100 on individual axes.----
for (i in seq_along(1:25)) {
  print(
      ggplot(subset(psyn, psyn$curveID == i), aes(x = invBT_eV, y = Photo)) +
          stat_smooth(method = "lm", se = TRUE, fill = NA,
                      formula = y ~ poly(x, 2, raw = TRUE),
                      aes(colour = factor(curveID)),
                  data = subset(psyn, psyn$curveID == i)) +
          geom_point(aes(color = as.factor(curveID)), data = subset(psyn, psyn$curveID == i)) +
      xlab(expression(paste('Leaf temperature ', '<1/',italic('kT'),'>', ' (',  eV^{-1}, ')'))) +
      ylab(expression("Assimilation rate (" * mu ~ "mol" ~m^-2 ~s^-1 * ")")) +
      scale_y_continuous(trans = "log", breaks = trans_breaks("log", function(x) exp(x), n = 3),
                         labels = trans_format("log", math_format(e^.x))) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"))
  )
}


#--Plot a single curve (e.g., curveID = 47).----
ggplot(subset(psyn, psyn$curveID == "1"),
       aes(x = invBT_eV, y = Photo)) +
    stat_smooth(
        method = "lm",
        se = TRUE,
        fill = NA,
        formula = y ~ poly(x, 2, raw = TRUE),
        aes(colour = factor(curveID))
    ) +
    geom_point(aes(color = as.factor(curveID))) +
    xlab(expression(paste(
        'Leaf temperature ', '<1/', italic('kT'), '>', ' (',  eV ^ {
            -1
        }, ')'
    ))) +
    ylab(expression("Assimilation rate (" * mu ~ "mol" ~ m ^ -2 ~ s ^ -1 * ")")) +
    scale_y_continuous(
        trans = "log",
        breaks = trans_breaks("log", function(x)
            exp(x), n = 3),
        labels = trans_format("log", math_format(e ^ .x))
    ) +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    )


#==============================================================================
# Part 3: TPCFitting
#==============================================================================
#--Load TPCfitting code----
source("./R scripts/Michaletz_2018b.R") # Pawar code from Michaletz (2018)
#source("./code from Pawar & Dell/TPCFitting_v3_stm_fixed.R")


# Boltzmann
# Time execution
start.time <- Sys.time()
# Execute function
BA_results <-
    TPCFit(
        Data = psyn,
        trait = "Photo",
        ID = "curveID",
        temper = "Tleaf",
        species = "Taxon",
        traitName = "rate",
        PLOT = TRUE,
        OverPLOT = FALSE,
        Model = "Boltzmann",
        SchoolTpk = TRUE,
        rand.st = TRUE,
        n.rand = 100
    )
# Calculate time taken
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
# n.rand=100 took 1.49667 mins
median(BA_results$E_boltz, na.rm = T)
#AST: why is it zero?

# Schoolfield
# Time execution
start.time <- Sys.time()
# Execute function
SS_results <-
    TPCFit(
        Data = psyn,
        trait = "Photo",
        ID = "curveID",
        temper = "Tleaf",
        species = "Taxon",
        traitName = "rate",
        PLOT = TRUE,
        OverPLOT = FALSE,
        Model = "Schoolfield",
        SchoolTpk = TRUE,
        rand.st = TRUE,
        n.rand = 100
    )
# Calculate time taken
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
# n.rand=100 took 2.358905 mins
# n.rand=1000 took 20.06794 mins; this made no difference for median E (0.2918334 for 100 vs. 0.2918333 for 1000)

# Add in grouping variables from psyn dataframe
SS_results$Taxon <- psyn$Taxon[match(SS_results$id, psyn$curveID)]
SS_results$Genus <- psyn$Genus[match(SS_results$id, psyn$curveID)]
SS_results$Species <- psyn$Species[match(SS_results$id, psyn$curveID)]
SS_results$Pathway <- psyn$Pathway[match(SS_results$id, psyn$curveID)]
SS_results$Site <- psyn$Site[match(SS_results$id, psyn$curveID)]
SS_results$Treatment <- psyn$Treatment[match(SS_results$id, psyn$curveID)]
SS_results$Source <- psyn$Source[match(SS_results$id, psyn$curveID)]
SS_results$SourceCompiled <- psyn$SourceCompiled[match(SS_results$id, psyn$curveID)]
SS_results$SourcePrimary <- psyn$SourcePrimary[match(SS_results$id, psyn$curveID)]


# 95% CI for each leaf
SS_results$CIlow <- SS_results$E_sch - 1.96 * SS_results$E_SE_sch
SS_results$CIhi <- SS_results$E_sch + 1.96 * SS_results$E_SE_sch
SS_results$DiffFromPred <- ifelse("0.32" > SS_results$CIlow,
                                   ifelse("0.32" < SS_results$CIhi, 0, 1), 1) # 1=is different from 0.32, 0 is not different

# Calculate species mean value
sppMeanE <-
    ddply(
        SS_results,
        .(Taxon),
        summarize,
        E_sch = mean(E_sch, na.rm = T)
        #Q10 = mean(Q10, na.rm = T)
    )
# Calculate species median value
sppMedianE <- ddply(SS_results, .(Taxon), summarize,
                    E_sch = median(E_sch, na.rm = T)
                    #Q10=median(Q10, na.rm = T)
                    )

#--All samples: mean (95% CI) and median (95% CI)
# This is biased in favor of spp w/ many measured leaves
mean(sppMedianE$E_sch)
# Mean 95% CI (0.43 to 0.64); (0.36 to 0.59) for fixed
mean(SS_results$E_sch, na.rm = T) - 1.96 * (sd(SS_results$E_sch, na.rm = T) / sqrt(length(SS_results$E_sch)))
mean(SS_results2$E_sch, na.rm = T) + 1.96 * (sd(SS_results2$E_sch, na.rm = T) / sqrt(length(SS_results2$E_sch)))
# Median (0.30); 0.22 for fixed
median(SS_results$E_sch, na.rm = T)
# Median 95% CI (0.23 to 0.37); (0.17 to 0.27) for fixed see https://stats.stackexchange.com/questions/184516/why-is-the-95-ci-for-the-median-supposed-to-be-%C2%B11-57iqr-sqrtn
median(SS_results$E_sch, na.rm = T) - 1.57 * IQR(SS_results$E_sch, na.rm = T) / sqrt(length(SS_results$E_sch))
median(SS_results$E_sch) + 1.57 * IQR(SS_results$E_sch) / sqrt(length(SS_results$E_sch))

#--Species means: mean (95% CI) and median (95% CI)
# This is a way to avoid biasing towards highly-sampled spp
# Mean (0.74); 0.64 for fixed
mean(sppMeanE$E_sch)
# Mean 95% CI (0.57 to 0.91); (0.43 to 0.86) for fixed
mean(sppMeanE$E_sch) - 1.96 * (sd(sppMeanE$E_sch) / sqrt(length(sppMeanE$E_sch)))
mean(sppMeanE$E_sch) + 1.96 * (sd(sppMeanE$E_sch) / sqrt(length(sppMeanE$E_sch)))
# Median (0.72); 0.42 for fixed
median(sppMeanE$E_sch)
# Median 95% CI (0.55 to 0.89); (0.27 to 0.56) for fixed see https://stats.stackexchange.com/questions/184516/why-is-the-95-ci-for-the-median-supposed-to-be-%C2%B11-57iqr-sqrtn
median(sppMeanE$E_sch) - 1.57 * IQR(sppMeanE$E_sch) / sqrt(length(sppMeanE$E_sch))
median(sppMeanE$E_sch) + 1.57 * IQR(sppMeanE$E_sch) / sqrt(length(sppMeanE$E_sch))

#--Species medians: mean (95% CI) and median (95% CI)
# This is a way to avoid biasing towards highly-sampled spp
# Mean (0.63); 0.54 for fixed
mean(sppMedianE$E_sch)
# Mean 95% CI (0.46 to 0.80); (0.31 to 0.76) for fixed
mean(sppMedianE$E_sch)-1.96*(sd(sppMedianE$E_sch)/sqrt(length(sppMedianE$E_sch)))
mean(sppMedianE$E_sch)+1.96*(sd(sppMedianE$E_sch)/sqrt(length(sppMedianE$E_sch)))
# Median (0.47); 0.27 for fixed
median(sppMedianE$E_sch)
# 95% CI (0.32 to 0.62); (0.14 to 0.39) for fixed see https://stats.stackexchange.com/questions/184516/why-is-the-95-ci-for-the-median-supposed-to-be-%C2%B11-57iqr-sqrtn
median(sppMedianE$E_sch) - 1.57 * IQR(sppMedianE$E_sch, na.rm = T) / sqrt(length(sppMedianE$E_sch))
median(sppMedianE$E_sch) + 1.57 * IQR(sppMedianE$E_sch, na.rm = T) / sqrt(length(sppMedianE$E_sch))


#--Plots
#--Fig S1: Plot all fitted curves together
# First, use aggregate to get min and max Tleaf for each leaf, then convert result to dataframe
minMax <- aggregate(Tleaf ~ curveID, psyn, function(x) c(Tl_C_min = min(x), Tl_C_max = max(x)))
minMax <- cbind(minMax[-ncol(minMax)], minMax[[ncol(minMax)]])
# Define function to make sequence between min and max Tleaf_C (including max)
# see https://stackoverflow.com/questions/28419281/missing-last-sequence-in-seq-in-r; My function
# was modified to just add the max value "to" to the end of the sequence, so cases
# that already end on that will have a repeated value (which is OK)
seqlast <- function(from, to, by) {
    vec <- do.call(what = seq, args = list(from, to, by))
    return(c(vec, to))
}

# Make sequence of Tleaf between min and max for each leaf
SS_fits <-
    data.frame(
        id = rep(minMax$curveID, minMax$Tl_C_max - minMax$Tl_C_min + 2),
        #Tleaf_C=unlist(mapply(seq, minMax$Tl_C_min, minMax$Tl_C_max))); rm(minMax)
        Tleaf_C = unlist(mapply(
            seqlast, minMax$Tl_C_min, minMax$Tl_C_max, 1
        ))
    )
#rm(minMax)
# Merge Tleaf sequence with fitted parameter values from SS_results2 (this will omit bad curves)
SS_fits <- merge(SS_fits, SS_results, by = "id")
# Make Tleaf_K
SS_fits$Tleaf_K <- SS_fits$Tleaf_C + 273.15
# Now use Schoolfield function with TPCFit results to give fitted psyn value for each Tleaf for each curveID
SS_fits$Apred_umol_m2_s <-
    exp(
        Schoolfield(
            lnB0 = SS_fits$lnB0_sch,
            E = SS_fits$E_sch,
            E_D = SS_fits$E_D_sch,
            T_h = SS_fits$T_h_sch,
            temp = SS_fits$Tleaf_K,
            SchoolTpk = TRUE
        )
    )
# Plot modified Arrhenius plot (1/kT on bottom) [USE THIS]
FigS1b <-
    ggplot(SS_fits, aes(1 / (0.00008617 * (Tleaf_C + 273.15)), Apred_umol_m2_s, group =
                            factor(id))) +
    geom_line(alpha = 0.4) +
    xlab(expression(paste(
        'Leaf temperature ', '1/', italic('kT'), ' (',  eV ^ {
            -1
        }, ')'
    ))) +
    ylab(expression("Assimilation rate (" * mu ~ "mol" ~ m ^ -2 ~ s ^ -1 * ")")) +
    scale_x_continuous(sec.axis = sec_axis(trans = ~ (1 / (. * 0.00008617)) -
                                               273.15 , name = expression(paste(
                                                   "Leaf temperature (", degree, C, ")"
                                               )))) +
    scale_y_continuous(
        trans = "log",
        breaks = trans_breaks("log", function(x)
            exp(x), n = 3),
        labels = trans_format("log", math_format(e ^ .x))
    ) +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    )
FigS1b

# Plot a single curve with data and S-S fit (curveID=131)
FigS1a <-
    ggplot(subset(psyn, psyn$curveID == "2"), aes(x = 1 / (0.00008617 * (Tleaf + 273.15)), y = Photo)) +
    geom_point(
        shape = 21,
        size = 2.75,
        col = "black",
        fill = "white"
    ) +
    geom_line(
        data = subset(SS_fits, SS_fits$id == "131"),
        aes(x = 1 / (0.00008617 * (Tleaf_K)), y = Apred_umol_m2_s),
        alpha = 0.4
    ) +
    xlab(expression(paste(
        'Leaf temperature ', '1/', italic('kT'), ' (',  eV ^ {
            -1
        }, ')'
    ))) +
    ylab(expression("Assimilation rate (" * mu ~ "mol" ~ m ^ -2 ~ s ^ -1 * ")")) +
    scale_x_continuous(sec.axis = sec_axis(trans = ~ (1 / (. * 0.00008617)) -
                                               273.15 , name = expression(paste(
                                                   "Leaf temperature (", degree, C, ")")))) +
    scale_y_continuous(
        trans = "log",
        breaks = trans_breaks("log", function(x)
            exp(x), n = 3),
        labels = trans_format("log", math_format(e ^ .x))
    ) +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    )
FigS1a

#Plot E w/ 95% CI for all leaves
ggplot(SS_results, aes(x = Taxon, y = E_sch)) +
    geom_point(
        shape = 21,
        size = 2.75,
        col = "black",
        fill = "white"
    ) +
    geom_hline(aes(yintercept = 0.32)) +
    geom_errorbar(
        aes(
            ymin = E_sch - 1.96 * E_SE_sch,
            ymax = E_sch + 1.96 * E_SE_sch
        ),
        width = .2,
        position = position_dodge(.9)
    ) + scale_y_continuous(limits = c(0, 3.5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_bw()

# E boxplots by taxon

ggplot(SS_results, aes(x = Taxon, y = E_sch, fill = id_spp)) +
    geom_boxplot() +
    guides(fill = FALSE) +
    geom_hline(aes(yintercept = 0.32))  +
    theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.25
    ))

# Topt boxplots by taxon
ggplot(SS_results, aes(x = Taxon, y = (T_pk_sch - 273.15))) +
    geom_boxplot(fill = "grey") +
    guides(fill = FALSE) +
    #scale_y_continuous(limits = c(0,45)) +
    ylab(expression(paste(
        'Optimal temperature for photosynthesis (', degree, C, ')'
    ))) +
    xlab(expression("Taxon")) +
    theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.25
    )) +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    )

# Compare distributions across sites
ggplot(SS_results, aes(x = E_sch, fill = Site)) +
    geom_density(alpha = .3)

ggplot(subset(SS_results, !is.na(SS_results$Site)),
       aes(x = E_sch, fill = Site)) +
    geom_density(alpha = .3)

# thermoregulation
psyn
ggplot(psyn, aes(x = Tair, y = Tleaf)) +
    geom_point(aes(color = as.factor(Taxon)), data = psyn) +
    stat_smooth(
        method = "lm",
        se = TRUE,
        fill = NA,
        formula = y ~ poly(x, 2, raw = TRUE),
        aes(colour = factor(Taxon)),
        data = psyn
    ) +
    geom_abline(intercept = 0, slope = 1) +
    xlab(expression('Air temperature' ~ (degree * C))) +
    ylab(expression("Leaf temperature" ~ (degree * C))) +
    theme_bw(base_size = 12) +
    theme(
        #legend.position = "none",
        legend.title = element_text("Species"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    ) +
    #scale_colour_manual(values = rep(cbPalette, times = 3)) +
    # plot each 22 times to maximize solids (there are 9 colors)
    coord_equal()

##per site
ggplot(psyn, aes(x = Tair, y = Tleaf)) +
    geom_point(aes(color = as.factor(Site)), data = psyn) +
    stat_smooth(
        method = "lm",
        se = TRUE,
        fill = NA,
        formula = y ~ poly(x, 2, raw = TRUE),
        aes(colour = factor(Site)),
        data = psyn
    ) +
    geom_abline(intercept = 0, slope = 1) +
    xlab(expression('Air temperature' ~ (degree * C))) +
    ylab(expression("Leaf temperature" ~ (degree * C))) +
    theme_bw(base_size = 12) +
    theme(
        #legend.position = "none",
        legend.title = element_text("Species"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    ) +
    #scale_colour_manual(values = rep(cbPalette, times = 3)) +
    # plot each 22 times to maximize solids (there are 9 colors)
    coord_equal()

###por todo
psyn %>% dplyr::filter(Site %in% c("QUE", "WAY")) %>%
ggplot(aes(x = Tair, y = Tleaf)) +
    geom_point(aes(color = as.factor(Site))) +
    stat_smooth(
        method = "lm",
        se = TRUE,
        fill = NA,
        formula = y ~ poly(x, 2, raw = TRUE),
        aes(colour = factor(Site))
    ) +
    geom_abline(intercept = 0, slope = 1) +
    xlab(expression('Air temperature' ~ (degree * C))) +
    ylab(expression("Leaf temperature" ~ (degree * C))) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
    ) +
    #scale_colour_manual(values = rep(cbPalette, times = 3)) +
    # plot each 22 times to maximize solids (there are 9 colors)
    coord_equal() +
    facet_grid(Taxon ~ .)

