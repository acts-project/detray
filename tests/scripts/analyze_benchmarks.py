import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

benchmarks = pd.read_csv('benchmarks_history.csv').groupby(['name'])
for bench in benchmarks :
    gbenchmarks = bench[1].groupby(['group']);
    for gbench in gbenchmarks :                
        plt.plot(gbench[1]['commit'], gbench[1]['cpu_time'], marker="o", label=gbench[0])
    plt.legend()
    plt.savefig(bench[0] + '.png', format='png')
    plt.cla()
