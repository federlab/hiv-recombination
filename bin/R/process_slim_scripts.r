
process_slim_scripts <- function(dir, numSampled){

                                        #dirtmp <-  dir
genomeLength = 500

#Let's read in our hxb2 reference sequence. It's against this seq that all our dynamics will happen
hxb2_env <- strsplit(readLines(path_to_hxb2_env)[2], "")[[1]]

#We'll convert it to numerics, because it'll make swapping indices easier
dict <- paste(c(1, 2, 3, 4))
names(dict) <- c("A", "T", "C", "G")
hxb2_numeric <- as.numeric(str_replace_all(hxb2_env, dict))


res <- readLines(paste0(dir, "slim_output.txt"))

#Find where the mutational list starts - this has important info about mutational position 
fixedMutsStart <- which(grepl("^#OUT.*F", res)) + 2
segMutsStart <- which(grepl("^#OUT.*SS", res)) + 2
genStart <- which(res == "Genomes:") + 1

#we need ALL of the positions which have segregating mutations
mutInf_by_t <- foreach(i = 1:length(segMutsStart))%do%{
    mutInf <-matrix((unlist(lapply(strsplit(res[ (segMutsStart[i] ): (genStart[i] - 2)], " "), function(x){x[c(1, 4, 10)]} ) )), ncol = 3, byrow = TRUE)
    colnames(mutInf) <- c("mutnum", "mutpos", "nuc")
    tbl_df(mutInf) %>% mutate(mutnum = as.numeric(mutnum), mutpos = as.numeric(mutpos))
}

#Let's first make a dictionary of mutational numbers to genomic positions
numMuts_by_t <- lapply(mutInf_by_t, nrow)

#Let's also set up the sampled genomes at each timepoint
#Here, each row represents a genome with any mutations
#NOT observing a number in a row means its identity is its ancestor
genoInf_by_t <- foreach(i = 1:length(genStart))%do%{
    readInfo <- lapply(strsplit(res[ (genStart[i] ):(genStart[i]  + numSampled - 1)], " "), function(x){ as.numeric(x[-(1:2)])})
}

#Finally, we need to figure out which mutations are fixed at each timepoint
fixedMuts_by_t <- foreach(i = 1:length(fixedMutsStart))%do%{
    if(fixedMutsStart[i] < segMutsStart[i] - 3){
        mutInf <-matrix((unlist(lapply(strsplit(res[ (fixedMutsStart[i] ): (segMutsStart[i] - 3)], " "), function(x){x[c(1, 4, 10)]} ) )), ncol = 3, byrow = TRUE)
        colnames(mutInf) <- c("mutnum", "mutpos", "nuc")
        tbl_df(mutInf) %>% mutate(mutnum = as.numeric(mutnum), mutpos = as.numeric(mutpos))
    }else{
        #If no fixed mutations, return NULL
        return(NULL)
    }
}

allMutPositions <- sort(unique(unlist(lapply(mutInf_by_t, function(x){return(x$mutpos)}))))

#This will be the template datastructure
dat_template <- foreach(i = 1:(genomeLength))%do%{
    foreach(j = 1:(genomeLength))%do%{
        matrix(rep(0, 6 * 6), nrow =  6)
    }
}

#Let's make a map of our sampling timepoints
timepoints <- tbl_df(unique(as.numeric(gsub("#OUT: | [A-Z].*", "", res[which(grepl("^#OUT", res))])))) %>% 
mutate(ind = 1:n()) %>% select(ind, value)

write.table(timepoints, paste0(dir, "timepoint_info.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)

#Let's also grab the recombination rate:
recombination_rate <- gsub(".*\\(|\\).*", "", res[grep("initializeRecombinationRate", res)])

#For each of the timepoints
dats_by_tp <- foreach(tp = 1:length(fixedMutsStart))%do%{


    mutInf <- mutInf_by_t[[tp]] %>% 
        mutate(nuc_num = as.numeric(str_replace_all(nuc, dict)))
    numMuts <- numMuts_by_t[[tp]]
    genoInf <- genoInf_by_t[[tp]]
    fixedMuts <- fixedMuts_by_t[[tp]]

    #First, we update our hxb2_numeric reference based on fixedMuts
    hxb2_tp <- hxb2_numeric
    if(!is.null(fixedMuts[1])){
        hxb2_tp[fixedMuts$mutpos] <-  as.numeric(str_replace_all(fixedMuts$nuc, dict))
    }
    #Here - in some instances, we're replacing a position with itself, 
    #even though this should be a fixed mutation. Will need to look into this behavior


#Let's set up our data structure of nested lists
    dat <- dat_template

    convertToRead <- function(x, hxb2_tp, mutInf){

        subs <- mutInf %>% filter(is.element(mutnum, x))
        bitstring <- hxb2_tp
        if(length(subs$mutpos) != length(subs$nuc_num)){
            print(mutnum)
        }
        bitstring[subs$mutpos + 1] <- subs$nuc_num
        return(bitstring)
    }

    readsAsBinary <- lapply(genoInf, convertToRead, hxb2_tp, mutInf)
    readsAsBinaryMat <- (matrix(unlist(readsAsBinary), nrow = length(readsAsBinary), ncol = genomeLength , byrow = TRUE))
    colnames(readsAsBinaryMat) <- c(0:(genomeLength - 1))
    readsAsBinaryMat <- tbl_df(readsAsBinaryMat)


    #Just as a sanity check, let's see if the reads that we expect to be polymorphic are

    mutated_pos <- unique(mutInf$mutpos)

    captureNulls <- foreach(i = 1:length(allMutPositions))%do%{

        if(i %% 50 == 0){
            print(paste0("Timpoint ", tp, ", ", i, "/", length(allMutPositions)))
        }

        loc1 <- paste(allMutPositions[i])

        foreach(j = i:length(allMutPositions))%do%{

            loc2 <- paste(allMutPositions[j])

            if(i != j){

                relCols <- readsAsBinaryMat %>% select(loc1, loc2 )
                names(relCols) <- c("l1", "l2")

                toIncrement <- relCols %>% group_by(l1, l2) %>% 
                    summarize(n = n(), .groups = "drop") 

                dat[[allMutPositions[i] + 1]][[allMutPositions[j] + 1]][as.matrix(toIncrement[,c(1, 2)])] <- toIncrement$n
                dat[[allMutPositions[j] + 1]][[allMutPositions[i] + 1]][as.matrix(toIncrement[,c(1, 2)])] <- toIncrement$n

            }else if(i == j){
                relCols <- readsAsBinaryMat %>% select(loc1)
                names(relCols) <- c("l1")

                toIncrement <- relCols %>% group_by(l1) %>% 
                    summarize(n = n(), .groups = "drop") %>% 
                        mutate(l2 = l1) %>% select(l1, l2, n)
                dat[[allMutPositions[i] + 1]][[allMutPositions[i] + 1]][as.matrix(toIncrement[,c(1, 2)])] <- toIncrement$n
            }
        }
    }


    #I'm going to comment out this analysis for now because it runs pretty slowly

    if(0 == 1){
        print("Computing pairwise R^2")

                                        #Let's visualize the decay of r2 over time too
        muts <- sort(mutInf_by_t[[1]]$mutpos)

        r2 <- foreach(i = 1:(length(muts) - 1), .combine = "rbind", .errorhandling='pass')%do%{

            foreach(j = (i + 1):length(muts), .combine = "rbind", .errorhandling='pass')%do%{

                                        #            print(paste0(i, "-", j))
                d <- dat[[muts[i] + 1]][[muts[j] + 1]]

                maxRows <- apply(d, 1, sum)
                maxCols <- apply(d, 2, sum)

                topRows <- order(maxRows, 1:6, decreasing = TRUE)[1:2]
                topCols <- order(maxCols, 1:6, decreasing = TRUE)[1:2]

                                        #if polymorphic
                if(sum(c(maxRows[topRows], maxCols[topCols]) > 0) == 4){

                    toComp <- d[topRows, topCols]

                    pAB <- toComp[1,1]/sum(toComp)

                    pA <- (apply(toComp, 1, sum)/sum(toComp))[1]
                    pB <- (apply(toComp, 2, sum)/sum(toComp))[1]

                    return(c(muts[i], muts[j], (pAB - pA * pB)^2/(pA * pB * (1 - pA) * (1 - pB))))
                }else{
                    return(c(muts[i], muts[j], NA))
                }
            }
        }

        print("Plotting R^2 decay")

        colnames(r2) <- c("l1", "l2", "r2")
        r2 <- tbl_df(r2)
        r2 <- r2 %>% mutate(diff = abs(l2 - l1)) %>% arrange(diff) %>% filter(diff > 0) %>% unique() 
        binsize <- 10
        runningav <- r2 %>% mutate(bin = floor(diff/binsize) * binsize) %>% group_by(bin) %>% summarize(av = mean(r2, na.rm = TRUE)) %>% 
            mutate(bin = bin + binsize/2)
        r2_plot <- r2 %>% ggplot() + geom_point(aes(diff, r2), alpha = 0.2, size = 0.5) +
            geom_line(aes(x = bin, y = av), col = "red",  data = runningav, size = 1) + theme_classic() +
                ggtitle( paste0("Generation: ",(timepoints %>% slice(tp))$value, ", recombination rate = ", recombination_rate)) + 
                    labs(x = "Distance (bp)", y = "R^2")

        pdf(paste0(dir, "r2_plots/r2_", tp, ".pdf"), height = 5, width = 5)
        print(r2_plot)
        dev.off()

    }


#Roughly, I think this is working
#Next question: how do I print it so that it gets read in correctly by the Zanini reshape command?

                                        #I *think* the first 6 elements should be
    
#1) locus 1, A and locus 1 A, 
#2) locus 1, A and locus 2 A, 
#3) locus 1, A and locus 3 A, 
#... etc

    if(0 == 1){

#    print("Reformatting data (slowly)")
    #This is very slow, but I don't know if it really makes sense to optimize at this point
        data_1d <- foreach(nuc1 = 1:6, .combine = "c")%do%{
            foreach(nuc2 = 1:6, .combine = "c")%do%{
                print(paste0("nuc1: ", nuc1, ", nuc2: ", nuc2))
                foreach(locus1 = 1:genomeLength, .combine = "c")%do%{
                    foreach(locus2 = 1:genomeLength, .combine = "c")%do%{
                        dat[[locus1]][[locus2]][nuc1, nuc2]
                    }
                }
            }
        }
 
    }


    print("Reformatting data")

    #This is MUCH faster
    L <- genomeLength
    a <- c(36 * (0:(L * L - 1)) + 1 + 6*0, 
           36 * (0:(L * L - 1)) + 1 + 6*1, 
           36 * (0:(L * L - 1)) + 1 + 6*2, 
           36 * (0:(L * L - 1)) + 1 + 6*3, 
           36 * (0:(L * L - 1)) + 1 + 6*4, 
           36 * (0:(L * L - 1)) + 1 + 6*5)

    data_1d <- unlist(dat)[c(a, a+1, a+2, a+3, a+4, a+5)]
    


    print(paste0("Writing to ", paste0(dir, "slim_formatted_t",tp, ".txt")))

    write.table(data_1d, paste0(dir, "slim_formatted_t",tp, ".txt"),  
                row.names = FALSE, col.names = FALSE, quote = FALSE)

}
}



#Here's a little script intended to generate and run slim simulations
generate_slim_script <- function(mu = 1e-5, rho = 1.4e-5, Ne = 1000, samplingGenerations = c(1000, 1200, 1400, 1600), M = 1000){
 
    slimscript <- paste0('initialize() {

 initializeSLiMOptions(nucleotideBased=T); 
 defineConstant("L", initializeAncestralNucleotides("',path_to_hxb2_env,'" ));

 initializeMutationTypeNuc("m1", 0.5, "f", 0.0);

 initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(',mu,'));
 initializeGenomicElement(g1, 0, L-1);
 initializeRecombinationRate(',rho,');

                         
}
1 {   sim.addSubpop("p1",', Ne,'); }\n',  paste0(apply(expand.grid(c(" late() { sim.outputFixedMutations(); }", paste0(" late() { p1.outputSample(",M,"); }")), samplingGenerations)[, c(2:1)], 1, paste, collapse = ""), collapse = "\n"), "\n")


    return(slimscript)
}

setup_slim_dir <- function(base_directory, sub_directory, mu, rho, Ne, samplingGenerations, M){

    system(paste0("mkdir ../../data/simulated/", base_directory, "/", sub_directory))
    script <- generate_slim_script(mu, rho, Ne, samplingGenerations, M)
    write(script, paste0("../../data/simulated/", base_directory, "/", sub_directory, "/slim_script.txt" ))
    system(paste0("mkdir ../../data/simulated/", base_directory, "/", sub_directory, "/r2_plots"))
    
    
}

setup_slim_trials <- function(base_directory, params, samplingGenerations){

    system(paste0("mkdir ../../data/simulated/", base_directory))

    sto <- foreach(i = 1:nrow(params))%do%{
    
        par <- params[i,]
        dirname <- paste(apply(rbind(names(par), par), 2,paste0, collapse = ""), collapse = "_")
        setup_slim_dir(base_directory, sub_directory = dirname, mu = par$mu, rho = par$rho, Ne = par$Ne, samplingGenerations = samplingGenerations, M = par$M)
    }
}
                     

run_slim_trials <- function(base_directory){

    dirpath <- paste0("../../data/simulated/",base_directory, "/")
    foreach(run = list.files(dirpath))%do%{

        system( paste0("slim ", dirpath, run, "/slim_script.txt", " > ", dirpath, run, "/slim_output.txt"))
                
    }
}


process_slim_trials <- function(base_directory){

    dirpath <- paste0("../../data/simulated/",base_directory, "/")
    foreach(run = list.files(dirpath))%do%{

        process_slim_scripts(paste0(dirpath, run, "/"), M)
                
    }

}


