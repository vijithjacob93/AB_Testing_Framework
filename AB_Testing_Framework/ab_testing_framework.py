import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.stats as stats

class ab_testing_framework:
    def __init__(self, desired_mde, metric_baseline, alpha = 0.05, beta = 0.8, num_tails = 2, weekly_treatment_size = 0, weekly_control_size = 0, total_treatment_size = 0, total_control_size = 0, max_num_weeks = 24):
        self.max_num_weeks = max_num_weeks
        self.weekly_treatment_size = weekly_treatment_size
        self.weekly_control_size =  weekly_control_size
        self.total_treatment_size =  total_treatment_size
        self.total_control_size =  total_control_size
        self.num_tails = num_tails
        self.desired_mde = desired_mde
        self.metric_baseline = metric_baseline
        self.alpha = alpha
        self.beta = beta
    
    def treatment_size(self, mde_var, num_weeks):#.weekly_treatment_size, self.weekly_control_size, self.total_control_size, self.p1, self.z_alpha, self.z_beta, self.num_weeks):
        self.z_alpha = stats.norm.ppf(1-self.alpha/self.num_tails) 
        self.z_beta = stats.norm.ppf(self.beta) 
        q = 1-self.metric_baseline
        k = self.total_control_size/(self.weekly_treatment_size*num_weeks) if self.weekly_control_size == 0 else self.weekly_control_size/self.weekly_treatment_size
        return self.weekly_treatment_size*num_weeks - ( (k + 1)/k ) * ( self.metric_baseline * q * (self.z_alpha + self.z_beta)**2 / mde_var**2 )
    
    def mde_solve(self):
        self.mde_over_time = pd.DataFrame({'num_weeks':np.arange(1,self.max_num_weeks), 'mde':np.nan})
        self.mde_over_time['mde'] = self.mde_over_time['num_weeks'].apply(lambda x: abs(fsolve(self.treatment_size, self.desired_mde, args = (x)))[0])#.weekly_treatment_size, self.weekly_control_size, self.total_control_size, self.metric_baseline, self.z_alpha, self.z_beta, x)))[0])
        self.mde_over_time['desired_mde'] = self.desired_mde
        return self.mde_over_time

    def plot_mde_over_time(self):
        fig, ax = plt.subplots()
        ax.plot(self.mde_over_time['num_weeks'], self.mde_over_time['mde'])
        ax.axhline(self.mde_over_time.loc[0,'desired_mde'], linestyle = '--', color = 'k')
        ax.set_xlabel('Number of Weeks')
        ax.set_ylabel('Minimum Detectable Effect')
        fig.suptitle('Test Runtime Duration vs Minimum Detectable Effect')