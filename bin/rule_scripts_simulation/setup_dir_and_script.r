library("optparse")

option_list = list(
  make_option(c("-d", "--directory"), type="character", default='tmp', 
              help="directory name for collection of runes", metavar="character"),
  make_option(c("-m", "--mu"), type="numeric", default="1e-06", 
              help="mutation rate [default= 1e-06]", metavar="character"),
  make_option(c("-r", "--rho"), type="numeric", default="1e-05", 
              help="recombination rate [default= 1e-05]", metavar="character"),
  make_option(c("-N", "--Ne"), type="numeric", default="1e04", 
              help="Effective population size [default= 1e04]", metavar="character"),
  make_option(c("-M", "--samplingdepth"), type="numeric", default="800", 
              help="Sampling depth [default= 800]", metavar="character"),
  make_option(c("-i", "--rep"), type="numeric", default=1, 
              help="rep", metavar="character"),
  make_option(c("-s", "--selection"), type="character", default='neutral', 
              help="selective regime [default = 'neutral', options = 'neutral, MPL']", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Input args
dirstub <- opt$directory
filepath <- paste0("/net/feder/vol1/project/hiv_recombination/data/simulated/", dirstub, "/")
mu <- opt$mu
rho <- opt$rho
Ne <- opt$Ne
M <- opt$samplingdepth
rep <- opt$rep
sel <- opt$selection

#Fixed args
path_to_hxb2_env <- "/net/feder/vol1/project/hiv_recombination/data/reference/hxb2/hxb2_env_full.fa.txt"
samplingGenerations <- 5*Ne + c(0, 150, 300, 450, 600)

trial_directory <- paste0("mu", mu, "_rho", rho, "_Ne", Ne, "_M", M, "_s", sel,"_rep", rep)

#If the date directory doesn't exist, make it
if(!file.exists(filepath)){
    system(paste0("mkdir ", filepath))
}

#create directory for the particular slim run
system(paste0("mkdir ", filepath, trial_directory))

#Write the slim file to run intself

if(sel == "neutral"){

slimscript <- paste0('initialize() {

 initializeSLiMOptions(nucleotideBased=T); 
 defineConstant("L", initializeAncestralNucleotides("',path_to_hxb2_env,'" ));

 initializeMutationTypeNuc("m1", 0.5, "f", 0.0);

 initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(',mu,'));
 initializeGenomicElement(g1, 0, L-1);
 initializeRecombinationRate(',rho,');
                        
}
1 {   sim.addSubpop("p1",', Ne,'); }\n',  paste0(apply(expand.grid(c(" late() { sim.outputFixedMutations(); }", paste0(" late() { p1.outputSample(",M,"); }")), samplingGenerations)[, c(2:1)], 1, paste, collapse = ""), collapse = "\n"), "\n")

}


if(sel == "MPL"){

#read in MPL simulation raw values

mpl.in <- read.table("/net/feder/vol1/project/hiv_recombination/data/reference/barton/total-selection.csv", sep = ",", header = TRUE)

#small_mpl_string <- paste(mpl.in[1:10, ]$s_MPL, collapse = ", ")
small_mpl_string <- paste(mpl.in$s_MPL, collapse = ", ")

slimscript <- paste0('initialize() {

 myscript = "return sample(c(', small_mpl_string, '), 1, replace = T);";

 initializeSLiMOptions(nucleotideBased=T); 
 defineConstant("L", initializeAncestralNucleotides("',path_to_hxb2_env,'" ));

 initializeMutationTypeNuc("m1", 0.5, "s", myscript);

 initializeGenomicElementType("g1", m1, 1, mmJukesCantor(',mu,'));
 initializeGenomicElement(g1, 0, L-1);
 initializeRecombinationRate(',rho,');
                         
}
1 {   sim.addSubpop("p1",', Ne,'); }\n',  paste0(apply(expand.grid(c(" late() { sim.outputFixedMutations(); }", paste0(" late() { p1.outputSample(",M,"); }")), samplingGenerations)[, c(2:1)], 1, paste, collapse = ""), collapse = "\n"), "\n")

}


write(slimscript, paste0(filepath, trial_directory, "/slim_script.txt" ))

