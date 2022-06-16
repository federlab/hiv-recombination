import os


def createVCF(indir, reference, outdir):
    """ Takes a directory of bam files and the reference sequence. Then creates
    VCF files and puts them in outdir.
    """
    for currFile in os.listdir(indir):
        #get the last entry so that if its bam.bai we can ignore it
        extension = (currFile.split("."))[-1]
        sample = (currFile.split("."))[0]
        if extension == "bam":
            os.system("bcftools mpileup --max-depth 30000 -f " + reference + " " + 
                    indir + currFile + 
                    " | bcftools call -mv -Ov -o " + outdir + 
                    sample + ".vcf")
    return 

