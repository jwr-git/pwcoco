# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 14:13:15 2020

@author: Jamie
"""

class Phenotype:
    def __init__(self, filename, to_save, snpname, snpbp, snpp, snppc):
        try:
            self.df = pd.read_csv(filename, sep="\t", index_col=False)
        except:
            e = sys.exc_info()
            print("Error: {}, cannot plot this data.".format(e))
            
        self.snpname = snpname
        self.snpbp = snpbp
        self.snpp = snpp
        self.snppc = snppc
        
        self.to_save = to_save
            
    def format_data(self, ld_col, snp_col="SNP", chr_col="Chr",
                    bp_col="bp", pval_col="p"):
        # We do not need all the data
        self.chr = self.df['Chr'][0]
        self.df = self.df[['SNP', 'bp', 'p', 'pC', ld_col]]
        self.df.rename(columns={ld_col: "ld"}, inplace=True)
        self.df[['bp', 'p', 'pC', 'ld']] = self.df[['bp', 'p', 'pC', 'ld']].apply(pd.to_numeric)
        
        # Append the lead SNP to the df as this will not be in the file given
        self.df = self.df.append({'SNP': self.snpname,
                                  'bp': self.snpbp,
                                  'p': self.snpp,
                                  'pC': self.snppc,
                                  'ld': 1},
                                 ignore_index=True)
        
        # Convert data for plotting
        self.df['bp'] = self.df['bp'] / 1000000
        with np.errstate(divide='ignore'):
            self.df['logp'] = np.where(self.df['p'] > 0, -np.log10(self.df['p']), 0)
            self.df['logpC'] = np.where(self.df['pC'] > 0, -np.log10(self.df['pC']), 0)
        
        # Colour according to LD
        self.df['ld_plot'] = pd.cut(self.df['ld'], 
                                    bins=[-1, 0, 0.2, 0.4, 0.6, 0.8, 1],
                                    labels=[-1, 0, 0.2, 0.4, 0.6, 0.8])
                
    def plot_region(self):
        self.format_data(self.snpname)
        
        fig, ax = plt.subplots()
        
        ax.scatter(self.df['bp'], self.df['logpC'], marker='.')
                
        #ax = sns.scatterplot(data=self.df,
        #                     x='bp', y='logpC',
        #                     hue="ld_plot", 
        #                     palette=sns.color_palette("ch:s=.75,rot=-.57"))
        
        idx = self.df['logpC'].idxmax()
        ax.annotate(self.df['SNP'][idx],
                    xy=(self.df['bp'][idx], self.df['logpC'][idx]), xycoords='data',
                    xytext=(-2,2),
                    textcoords='offset points', ha='center', va='bottom')
                
        plt.ylim(bottom=0)
        plt.xlabel("Chromosomal position (chr {}, Mb)".format(self.chr))
        plt.ylabel(u"log\u2081\u2080(P)")
        #plt.show()
        plt.savefig(self.to_save, bbox_inches='tight')
        
if __name__ == "__main__":
    import sys
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    if len(sys.argv) != 7:
        sys.exit("Not enough arguments for plotting")

    ph = Phenotype(str(sys.argv[1]), str(sys.argv[2]),
                   str(sys.argv[3]), 
                   float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]))
    ph.plot_region()
