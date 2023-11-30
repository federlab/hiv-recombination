#sample a start position
numsamp <- 100000

#draw a start position
starts <- floor(runif(numsamp, 0, 700))

#draw an insert size
inserts <- pmax(0, floor(runif(numsamp, 400, 700)) - 600)

simReads <- tibble((bind_cols(starts, starts + 300, starts + 300 + inserts, starts + 300 + inserts + 300)))
names(simReads) <- c("L_start", "L_end", "R_start", "R_end")

#Discard all reads that don't contain position 700
simReads <- simReads %>% filter((L_start < 700 & L_end > 700) | (R_start < 700 & R_end > 700))

overlap_helper <- function(x){
    return(between(1:1400, x[1], x[2]) | between(1:1400, x[3], x[4]))
}

overlap_boolean <- t(apply(simReads, 1, overlap_helper))

overlap_probability <- tibble(bind_cols(1:1400, apply(overlap_boolean, 2, mean)))
names(overlap_probability) <- c("pos", "p_sample")

empiricalDist <- overlap_probability %>% mutate(distFromCenter = abs(700 - pos)) %>% group_by(distFromCenter) %>%
    summarize(psamp = mean(p_sample))

empiricallyComputeProbOfSamplingGivenDeltaT <- function(disp){
    if(disp < 0 | disp > 700){return(0) }
    return(empiricalDist[disp + 1, ]$psamp    )
}

