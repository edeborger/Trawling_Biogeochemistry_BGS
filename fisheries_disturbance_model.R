#-------------------------------------------------------------------------------
#                       BGS : TRAWLING DISTURBANCE MODEL                        
#                                EMIL DE BORGER                                 
#-------------------------------------------------------------------------------
#
# Function "perturbrange" is a wrapper around functions from the CNPDIA
# package available at R-Forge, which contains the actual model of sediment 
# biogeochemistry and auxilliary functions. This function simulates chronic 
# trawling in sediments as a combination of an erosion and mixing impact, which
# also affects faunal activity. The accompanying RMD file also explains 
# required input.
# See De Borger et al. 2021, Biogeosciences (DOI: ) for the reasoning behind the 
# model.
#
# Functions "modif" and "depCalc" are needed to describe the depletion and
# recovery of the species community.
#
# Function "extend" adapts input data of forcing environmental variables to the
# length of the simulation.
#
# Function boxCompare is a plotting method to compare the changes to modelled 
# process rates between two model runs. 
# 
# Function plotProfile is a plotting method to plot profiles of carbon fractions 
# or other nutrients
#
# CONTENTS:
# 1. function perturbRange
# 2. function modif
# 3. function depCalc
# 4. function extend
# 5. function boxCompare
# 6. function plotProfile
#-------------------------------------------------------------------------------

# 1. perturbRange
# Simulates sediment biogeochemistry subject to chronic bottomo trawling.

perturbRange <- function(mixdist, 
                         erosdist, 
                         freq         = 0:5, 
                         yrs          = 15,
                         mudContent   = 10, 
                         longevity    = 2.5, 
                         sedimentpars = list(),
                         ox           = NULL, 
                         nit          = NULL, 
                         amm          = NULL, 
                         pho          = NULL, 
                         temp         = NULL, 
                         carb         = NULL
                         ) {
  
  output <- NULL # Empty vector for output.

  # Reformatting nutrients depending on input type.
  # NULL input returns default values.
  # LIST input returns sinusoid function as specified.
  # VECTOR input extends values to length of simulation through repetition of annual data.

  # Reformat environmental forcings to 1 year.
  carb1 <- extend(carb, 1)
  ox1   <- extend(ox, 1)
  nit1  <- extend(nit, 1)
  amm1  <- extend(amm, 1)
  pho1  <- extend(pho, 1)
  temp1 <- extend(temp, 1)

  # Reformat environmental forcings to spinup run years.
  carbext <- extend(carb, yrs)
  oxext   <- extend(ox, yrs)
  nitext  <- extend(nit, yrs)
  ammext  <- extend(amm, yrs)
  phoext  <- extend(pho, yrs)
  tempext <- extend(temp, yrs)

  # Setup grid with N number of layers for L centimeters deep.
  gr <- setup.grid.1D(x.up = 0, dx.1 = 0.01, N = 100, L = 350)

  # Steady state solution to initialize dynamic solution later.
  stead <- CNPDIAsolve(parms       = sedimentpars,
                       Grid        = gr,
                       CfluxForc   = carb,
                       O2bwForc    = oxext,
                       NO3bwForc   = nitext,
                       NH3bwForc   = ammext,
                       PO4bwForc   = phoext,
                       temperature = tempext
                       )

  # FIRST A RUN WITHOUT TRAWLING DISTURBANCE (CNPDIAperturb).
  # Please see ?CNPDIAperturb for an explanation of the arguments.

  for (i in freq) {
    out2 <- NULL
    if (i == 0) {
      # Spinup run
      out <- CNPDIAperturb(yini         = stead$y,
                           Grid         = gr,
                           parms        = sedimentpars,
                           times        = 1:(365 * yrs),
                           CfluxForc    = carbext,
                           O2bwForc     = oxext,
                           NO3bwForc    = nitext,
                           NH3bwForc    = ammext,
                           PO4bwForc    = phoext,
                           temperature  = tempext,
                           perturbType  = c("erode"), 
                           perturbDepth = 0, 
                           perturbTimes = seq(0, (365 * yrs), by = 1)
                           )

      # Preparing output of spinup run for starting new output run
      out <- tail(out[, 2:1201], 1)
      out <- matrix(data = out, byrow = FALSE, nrow = 100, ncol = 12)

      # Output run
      out <- CNPDIAperturb(yini         = out,
                           Grid         = gr,
                           parms        = c(sedimentpars),
                           times        = 1:365,
                           CfluxForc    = carb1,
                           O2bwForc     = ox1,
                           NO3bwForc    = nit1,
                           NH3bwForc    = amm1,
                           PO4bwForc    = pho1,
                           temperature  = temp1,
                           perturbType  = c("erode"), 
                           perturbDepth = 0, 
                           perturbTimes = seq(0, 365, by = 1)
                           )

      # Calculate budgets of the output
      Cbudg <- CNPDIAbudgetC(out)
      Nbudg <- CNPDIAbudgetN(out)
      Obudg <- CNPDIAbudgetO2(out)
      Pbudg <- CNPDIAbudgetP(out)
      
      # Calculate annually averaged depth profiles of TOC, FDET, SDET and nutrients.
      Tocavg  <- colMeans(tail(out[, 1269:1368], 365))
      Fdetavg <- colMeans(tail(out[, 2:101]    , 365))
      Sdetavg <- colMeans(tail(out[, 102:201]  , 365))
      
      O2avg  <- colMeans(tail(out[, 202:301], 365))
      NO3avg <- colMeans(tail(out[, 302:401], 365))
      NO2avg <- colMeans(tail(out[, 402:501], 365))
      NH3avg <- colMeans(tail(out[, 502:601], 365))
      DICavg <- colMeans(tail(out[, 702:801], 365))
      
      # Bind output in a dataframe
      out <- data.frame("Events"          = i,
                        "Mixdepth"        = mixDepth,
                        "Erosdepth"       = eroDepth,
                        "Times"           = paste0(ts, collapse = "_"),
                        "Oxicmin"         = Cbudg[[2]][1] / 100,
                        "Denitrification" = Cbudg[[2]][2] / 100,
                        "Anoxicmin"       = Cbudg[[2]][3] / 100,
                        "NH3prod"         = Nbudg[[2]][1] / 100,
                        "Nitrification"   = ((Nbudg[[2]][5] * 1.5) + (Nbudg[[2]][6] * 0.5)) / 100,
                        "ODUox"           = Obudg[[2]][2] / 100,
                        "DICflux"         = Cbudg[[1]][1, 3] / 100,
                        "CBurial"         = Cbudg[[1]][2, 5] / 100,
                        "O2flux"          = Obudg[[1]][1, 1] / 100,
                        "Oxicrel"         = Cbudg[[2]][1] / Cbudg[[2]][4],
                        "Denitrificrel"   = Cbudg[[2]][2] / Cbudg[[2]][4],
                        "Anoxicrel"       = Cbudg[[2]][3] / Cbudg[[2]][4],
                        "DIPsur"          = Pbudg[[1]][1, 6] / 100,
                        "DIPbot"          = Pbudg[[1]][2, 6] / 100,
                        "TotalMin"        = Cbudg[[2]][4] / 100,
                        "NH3flux"         = Nbudg[[1]][1, 5] / 100,
                        "NO3flux"         = Nbudg[[1]][1, 3] / 100,
                        "DINflux"         = Nbudg[[1]][1, 7] / 100,
                        "pNremoved"       = Cbudg[[2]][2] * 0.8 / Nbudg[[2]][1],
                        rbind(Tocavg),
                        rbind(Fdetavg),
                        rbind(Sdetavg),
                        rbind(O2avg),
                        rbind(NO3avg),
                        rbind(NH3avg),
                        rbind(DICavg),
                        rbind(NO2avg)
      )
      
      out2 <- rbind(out2, out)
      print(paste(i, " rep. ", 1))

    # THEN THE RUNS WITH 1:N TRAWLING EVENTS PER YEAR
    } else {
      for (j in 1:length(mixdist)) { # For each of the 30 randomly generated mixing depths...

        # Randomization step to see where in the year the event(s) take(s) place.
        ts <- sort(sample(1:365, i))
        ts <- c(ts, ts + rep(seq(from = 365, to = 365 * yrs, by = 365), each = i))

        # Selecting penetration characteristics
        mixDepth <- mixdist[j]  # select mixing depth
        eroDepth <- erosdist[j] # select erosion layer

        # Calculate faunal activity alterations.
        depletion <- depCalc(erosion = eroDepth, mix = mixDepth, mud = mudContent)                            # calculate depletion factor of fauna
        modifBio  <- modif(times = 1:(365 * (yrs + 1)), impacts = ts, d = depletion, r = 1 / longevity / 365) # calculate faunal recovery and impacts over the year
        BiotLong  <- cbind(time = 1:(365 * (yrs)), data = sedimentpars[["biot"]] * modifBio[1:(365 * yrs)])   # apply to bioturbation coefficient in model
        BiotShort <- cbind(time = 1:365, data = sedimentpars[["biot"]] * tail(modifBio, 365))                 # apply to bioturbation coefficient in model

        # Running the simulation activey using CNPDIAperturb. Please check ?CNPDIAperturb to see how disturbances are implemented.
        # Spinup run
        out <- CNPDIAperturb(parms        = sedimentpars,
                             times        = 1:(365 * yrs),
                             yini         = stead$y,
                             Grid         = gr,
                             CfluxForc    = carbext,
                             O2bwForc     = oxext,
                             NO3bwForc    = nitext,
                             NH3bwForc    = ammext,
                             biotForc     = list(data = BiotLong),
                             PO4bwForc    = phoext,
                             temperature  = tempext,
                             perturbType  = c("erode", "mix"), 
                             perturbDepth = c(eroDepth, mixDepth),
                             perturbTimes = ts[1:(yrs * i)]
                             )

        # Preparing output of spinup run for starting new output run
        out <- tail(out[, 2:1201], 1)
        out <- matrix(data = out, byrow = FALSE, nrow = 100, ncol = 12)

        # Output run
        out <- CNPDIAperturb(parms        = sedimentpars,
                             times        = 1:365,
                             yini         = out,
                             Grid         = gr,
                             CfluxForc    = carb1,
                             O2bwForc     = ox1,
                             NO3bwForc    = nit1,
                             NH3bwForc    = amm1,
                             biotForc     = list(data = BiotShort),
                             PO4bwForc    = pho1,
                             temperature  = temp1,
                             perturbType  = c("erode", "mix"),
                             perturbDepth = c(eroDepth, mixDepth),
                             perturbTimes = ts[1:i]
                             )

        # Calculate budgets of the output
        Cbudg <- CNPDIAbudgetC(out)
        Nbudg <- CNPDIAbudgetN(out)
        Obudg <- CNPDIAbudgetO2(out)
        Pbudg <- CNPDIAbudgetP(out)

        # Calculate annually averaged depth profiles of TOC, FDET, SDET and nutrients.
        Tocavg  <- colMeans(tail(out[, 1269:1368], 365))
        Fdetavg <- colMeans(tail(out[, 2:101]    , 365))
        Sdetavg <- colMeans(tail(out[, 102:201]  , 365))

        O2avg  <- colMeans(tail(out[, 202:301], 365))
        NO3avg <- colMeans(tail(out[, 302:401], 365))
        NO2avg <- colMeans(tail(out[, 402:501], 365))
        NH3avg <- colMeans(tail(out[, 502:601], 365))
        DICavg <- colMeans(tail(out[, 702:801], 365))

        # Bind output in a dataframe
        out <- data.frame("Events"          = i,
                          "Mixdepth"        = mixDepth,
                          "Erosdepth"       = eroDepth,
                          "Times"           = paste0(ts, collapse = "_"),
                          "Oxicmin"         = Cbudg[[2]][1] / 100,
                          "Denitrification" = Cbudg[[2]][2] / 100,
                          "Anoxicmin"       = Cbudg[[2]][3] / 100,
                          "NH3prod"         = Nbudg[[2]][1] / 100,
                          "Nitrification"   = ((Nbudg[[2]][5] * 1.5) + (Nbudg[[2]][6] * 0.5)) / 100,
                          "ODUox"           = Obudg[[2]][2] / 100,
                          "DICflux"         = Cbudg[[1]][1, 3] / 100,
                          "CBurial"         = Cbudg[[1]][2, 5] / 100,
                          "O2flux"          = Obudg[[1]][1, 1] / 100,
                          "Oxicrel"         = Cbudg[[2]][1] / Cbudg[[2]][4],
                          "Denitrificrel"   = Cbudg[[2]][2] / Cbudg[[2]][4],
                          "Anoxicrel"       = Cbudg[[2]][3] / Cbudg[[2]][4],
                          "DIPsur"          = Pbudg[[1]][1, 6] / 100,
                          "DIPbot"          = Pbudg[[1]][2, 6] / 100,
                          "TotalMin"        = Cbudg[[2]][4] / 100,
                          "NH3flux"         = Nbudg[[1]][1, 5] / 100,
                          "NO3flux"         = Nbudg[[1]][1, 3] / 100,
                          "DINflux"         = Nbudg[[1]][1, 7] / 100,
                          "pNremoved"       = Cbudg[[2]][2] * 0.8 / Nbudg[[2]][1],
                          rbind(Tocavg),
                          rbind(Fdetavg),
                          rbind(Sdetavg),
                          rbind(O2avg),
                          rbind(NO3avg),
                          rbind(NH3avg),
                          rbind(DICavg),
                          rbind(NO2avg)
                          )
        
        out2 <- rbind(out2, out)
        print(paste(i, " rep. ", j)) # To keep track of where the function is while running.
      }
    }
    output <- rbind(output, out2)
  }
  return(output)
}

# 2. modif 
# Logistic growth modifier.
# This function simulates  logistic growth of a population with a certain 
# growth rate r, subject to impacts (at times impacts) which cull the population
# with a depletion factor d. This transience is scaled relative to the carrying
# capacity K.

modif <- function(times, impacts, d, r) {
  K <- 1
  out <- NULL
  out <- cbind(1, out)

  for (i in 2:(length(times) + 1)) {
    N <- out[i - 1]
    if (i %in% impacts) {
      newout <- pmax(0.1, N * pmax(0, (1 - d)))
    } else {
      newout <- N + r * N * (1 - N / K)
    }
    out <- c(out, newout)
  }
  return(out[-1])
}

# 3. depCalc
# Calculates depletion of species abundance (d) based on the penetration 
# depth and mud concent according to Sciberras et al., (2018).

depCalc <- function(erosion, mix, mud = 10) {
  pdepth <- erosion + mix
  d      <- pdepth * 3 + mud * 0.3
  d      <- d / 100
  return(round(d, 2))
}

# 4. Extend
# Handles environmental forcings dependent on input type.

extend <- function(var, yrs) {
  if (is.list(var)) {
    out <- var
  } else if (is.null(var)) {
    out <- NULL
  } else {
    out <- cbind(time = 1:(365 * yrs), data = rep(var, yrs))
    out <- list(data = out)
  }
  return(out)
}

# Plotting functions examples.

# 5. boxComparee
# Compare the output of two models for a given output variable.
# This function works for trawling frequencies of 0 - 5, other 
# frequencies require changes to the code.
# The whole function is a wrapper around the boxplot function.

boxCompare <- function(model1,
                       model2,
                       variable,
                       ylabz   = "",
                       ylimz   = c(),
                       plotpoints,
                       leg     = TRUE
                       ) {
  
  n <- 1:12
  n <- n+rep(seq(0, 2, by = 0.4), each = 2)
  l <- factor(c("Model1", "Model2"), levels = c("Model1", "Model2"))

  b <- boxplot(c(model1[,variable],model2[,variable])~rep(l, each = 5*30+1)+rep(model1[,"Events"], 2),
               xlab      = "",
               ylab      = "",
               main      = paste(variable),
               boxlty    = 0,
               whisklty  = 1,
               whisklwd  = 2,
               staplelwd = 2,
               staplewex = 0.25,
               boxwex    = 0.9,
               medcol    = c("blue", "red"),
               outcol    = NA,
               col       = c(alpha("blue", 0.1), alpha("red", 0.1)),
               axes      = FALSE, at = n)
  abline(v   = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5) + seq(0, 2, by  = 0.4),
         lty = 3,
         col = "gray75")
  axis(side    = 1,
       labels  = c("0", "1", "2", "3", "4", "5"),
       at      = c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5) + seq(0, 2, by = 0.4),
       lwd     = 2,
       cex.axis= 1.3)
  mtext(expression(paste("Frequency (y"^"-1", ")")),
        side  = 1,
        outer = FALSE,
        line  = 2.5,
        font  = 2,
        cex   = 1.0)
  axis(side = 2, lwd = 2, cex.axis= 1.3)
  mtext(ylabz, side = 2, line = 2.5, font = 2, cex = 1.1)
  box(bty = "L", lwd = 2)
  
  if(leg == TRUE){
    legend("bottom",
           c("Model 1", "Model 2"),
           title  = "Gear",
           bty    = "n",
           fill   = c(alpha("blue", 0.1), alpha("red", 0.1)),
           border = c(NA, NA),
           cex    = 1.4,
           horiz  = TRUE,
           pt.cex = 1.2)
    legend("bottom",
           c("Model 1", "Model 2"),
           title   = "Gear",
           bty     = "n",
           lty     = 1,
           col     = c("blue", "red"),
           seg.len = 0.8,
           cex     = 1.4,
           horiz   = TRUE)
  }
  
  if(plotpoints == 1){
    points(jitter(c(n[1], rep(n[c(3,5,7,9,11)] , each = 30)), 0.4),
           model1[[variable]],
           col = alpha("blue", 0.3),
           pch = 16,
           cex = 0.8)
    points(jitter(c(n[2], rep(n[c(4,6,8,10,12)], each = 30)), 0.4),
           model2[[variable]],
           col = alpha("red" , 0.3),
           pch = 16,
           cex = 0.8)
  }
}

# 6. plotProfile

`%notin%` <- Negate(`%in%`)

plotProfile <- function(model,
                        fraction, 
                        xlimz = c(),
                        ylimz = c(-3, 0),
                        depth,
                        colors,
                        plottitle = ""
                        ) {
  
  ## Extract carbon classes as a vector and calculate the mean and sd of the selected variable per frequency
  #### Tickler
  TOC  <- model[, 24:123]
  FDET <- model[, 124:223]
  SDET <- model[, 224:323]
  O2   <- model[, 324:423]
  NO3  <- model[, 424:523]
  NH3  <- model[, 524:623]
  DIC  <- model[, 624:723]
  NO2  <- model[, 724:823]
  TOC2 <- FDET + SDET        # The sum of FDET and SDET
  QUAL <- FDET / TOC2
  
  # Nutrients plotted normally
  if(fraction %notin% c("FDET", "SDET", "TOC2")) {
    sel     <- get(fraction)
    meansel <- apply(sel , 2, function(x) tapply(x, INDEX = model$Events, FUN = mean))
    xaxist  <- expression(paste("mmol m"^"-3"))
    
  # Carbon log-transformed to mol m-3
  } else{
    sel     <- get(fraction)
    meansel <- apply(sel , 2, function(x) tapply(x, INDEX = model$Events, FUN = mean))
    meansel <- meansel / 1000 # for mol m-3
    meansel <- meansel + 1
    xaxist  <- expression(paste("mol m"^"-3"))
  }
  
  if(fraction == "QUAL") {
    xaxist <- "Proportion of highly reactive OM"
    xlimz  <- c(0,1)
  }
  
  ## Plot the profile of the selected depth
  
  plot(meansel[1,],
       -depth,
       xlim = xlimz,
       ylim = ylimz,
       type = "l",
       lwd  = 2,
       col  = colors[1],
       axes = FALSE,
       xlab = "",
       ylab = "")
  axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  axis(side = 3, lwd = 2, cex.axis= 1.3)
  mtext(xaxist, side = 3, line = 2.5, font = 2, cex = 1.1)
  mtext(plottitle, side = 3, line = 5, font = 2, cex = 1.1)
  box(bty = "O", lwd = 2)
  
  lines(meansel[2,], -depth, lwd = 2, col = colors[2])
  lines(meansel[3,], -depth, lwd = 2, col = colors[3])
  lines(meansel[4,], -depth, lwd = 2, col = colors[4])
  lines(meansel[5,], -depth, lwd = 2, col = colors[5])
  lines(meansel[6,], -depth, lwd = 2, col = colors[6])
  
  legend("bottomright",
         paste(0:5),
         title = expression(paste("Frequency (y"^"-1", ")")),
         lwd   = 2,
         col   = colors,
         bty   = "n")
}
