import sys
#for running on cluster
sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
import zaniniUtil as zu

# rule label_viral_loads:
#     input:
#      dataDir + "{par}_{frag}/linkage/d_ratio",
#      vlDataDir + "viralLoads/viralLoad_{par}.tsv"
#     output:
#      dataDir + "{par}_{frag}/linkage/d_ratio_labeled"
#     script:
#      ruleScriptDir + "label_viral_loads.py"

#This script uses the coCounts_arr and list of segregatingLoci
#It takes them and saves a dataframe with the R^2 and D statistics
d_ratioData = snakemake.input[0]
vlDataDir = snakemake.input[1]
outDir = snakemake.output[0]

print("Viral Load Directory is", file = sys.stderr)
print(vlDataDir, file = sys.stderr)

print("D Ratio Directory is", file = sys.stderr)
print(d_ratioData, file = sys.stderr)

print("Output Directory is", file = sys.stderr)
print(outDir, file = sys.stderr)

