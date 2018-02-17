source('Scripts/R/paths.R')

load(PATHS$HERV.RANGES.DATA)

get.herv.lengths <- function(herv.ranges) {
  return(end(herv.ranges)-start(herv.ranges))
}

hervS1.lengths <- get.herv.lengths(hervS1.ranges)

hervS2.lengths <- get.herv.lengths(hervS2.ranges)

hervS1.lengths[hervS1.lengths > 5000] <- 5000
