import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

benchmarks = pd.read_csv(sys.argv[1]).groupby(['name'])
for bench in benchmarks :
    # Skip groups that contain "/"
    if "/" in bench[0]:
        continue
    gbenchmarks = bench[1].groupby(['group']);        
    maxl = 0
    for gbench in gbenchmarks :                
        maxl = max(maxl, len(gbench[1]))
    for gbench in gbenchmarks :
        curl = len(gbench[1])
        if gbench[0] == 'core' : # core tests are discontinued
            crange = range(0, curl)
        else :
            crange = range(maxl-curl, maxl)
        plt.plot(range(maxl-curl, maxl), gbench[1]['cpu_time'], marker="o", label=gbench[0])
    plt.legend()
    plt.savefig(bench[0] + '.png', format='png')
    plt.cla()
