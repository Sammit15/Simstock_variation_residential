#script to collect results from multiple EP simulations

import os
import pandas as pd

ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
OUTPUTS_DIR = os.path.join(ROOT_DIR, 'ep_outputs')
#weather = 'GBR_London.Gatwick.037760_IWEC' #directory containing simulation results subdirectories
sim_dirs = [f.name for f in os.scandir(OUTPUTS_DIR) if f.is_dir() ]    #list of subdirectories
combined_results = pd.concat([pd.read_csv(OUTPUTS_DIR+'/'+f+'/eplusout.csv') for f in sim_dirs], sort=True)
combined_results['run']=sim_dirs
combined_results.set_index(['run'], inplace=True)
combined_results.to_csv('combined_results.csv')