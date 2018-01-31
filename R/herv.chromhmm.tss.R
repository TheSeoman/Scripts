source('Scripts/R/paths.R')

load(PATHS$HERVS2.2KB.CHROMHMM.DATA)

tss.states <- c("1_TssA", "2_TssAFlnk", "10_TssBiv", "11_BivFlnk")

check.herv.for.tss <-
  function(herv.chromhmm.annotation.list, samples) {
    temp <-
      lapply(herv.chromhmm.annotation.list, function(herv.chromhmm.annotation) {
        temp2 <- list()
        for (herv in names(herv.chromhmm.annotation)) {
          for (annotation in herv.chromhmm.annotation[[herv]]) {
            temp3 <- list()
            if (annotation[4] %in% tss.states) {
              temp3[[length(temp3) + 1]] <- annotation
              cat(annotation, fill = T)
            }
            if (length(temp3) > 0) {
              temp2[[herv]] <- temp3
            }
          }
        }
      })
    return(temp)
  }

is.tss(annotation) {
  return(annotation[4] %in% )
}