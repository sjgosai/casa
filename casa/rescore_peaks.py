import numpy as np
import pandas as pd
import sys 
import csv 

def check_overlap(interval, array):
    height = array.shape[0]
    intervals = np.stack([np.tile(interval,(height,1)), array],axis=0)
    swaghook = (intervals[0,:,0] < intervals[1,:,0]).astype(int)
    return intervals[1-swaghook,np.arange(height),1] > intervals[swaghook,np.arange(height),0]

data = pd.read_table(sys.argv[1], header=None)
post = np.array([[float(k) for k in j.split(',')] for i,j in data[3].items()])
filt = check_overlap(np.array([-np.log(float(sys.argv[2])),np.log(float(sys.argv[2]))]),post) == False

data[4] = filt

data.to_csv(sys.argv[3],sep='\t', index=False, header=False)
