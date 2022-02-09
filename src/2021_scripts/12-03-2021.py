import sys
# sys.path.append('/net/feder/vol1/home/evromero/2021_hiv-rec/bin')
#for running on desktop
sys.path.append('/Volumes/feder-vol1/home/evromero/2021_hiv-rec/bin')
import os
import r2Analysis as r2
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import zaniniUtil as zu
import statsmodels.api as sm
from lmfit import minimize, Parameters, fit_report


# #Today I am going to run a regression on estimated rho values vs viral load
# #cluster run
# data_in = '/net/feder/vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/rho/'
# plot_out = '/net/feder/vol1/home/evromero/2021_hiv-rec/results/zanini/researchReports/rho_reg/'

#desktop run
data_in = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/data/zanini/analysis/rho/'
plot_out = '/Volumes/feder-vol1/home/evromero/2021_hiv-rec/results/zanini/researchReports/rho_reg/'

#Open the relevant data file 
rho_data = pd.read_csv(data_in + "rho_predictions.csv", index_col= 0)

rho_data = rho_data[rho_data['n_snps']  > 100]

#Take the log of viral load
rho_data['log_viral_load'] = np.log10(rho_data['viral_load'])

#set our variables
X = rho_data['log_viral_load']
y = rho_data['rho']

#run the linear model
X2 = sm.add_constant(X)
est = sm.OLS(y, X2)
est2 = est.fit()
intercept = est2.params['const']
slope = est2.params['log_viral_load']

#get the values from the model
x_values =np.log10([1000, 10000, 100000, 1000000])

def make_line(slope, intercept, x):
    return slope * x + intercept
fitted_values = [make_line(slope, intercept, x)  for x in x_values ]
print(est2.summary())

#plot rho vs viral load color by participant
sns.set(rc={'figure.figsize':(15,5)}, font_scale = 2)
sns.scatterplot(x = 'log_viral_load', y = 'rho', hue = 'participant', data = rho_data, alpha = 0.5)
plt.plot(x_values, fitted_values, color = 'Black')
plt.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))
plt.xlabel("Log10 Viral Load [Virions/ml]")
plt.ylabel("Estimated Rho")
plt.ylim(-0.1,1.1)
plt.tight_layout()
plt.savefig(plot_out + "Rho_vs_ViralLoad_Participant")
plt.close()

