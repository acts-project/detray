import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

benchmarks = pd.read_csv('benchmarks/benchmarks_history.csv')
unqiue_benchmarks = benchmarks.name.unique()

for ub in unqiue_benchmarks:
    benchmarks[benchmarks.name == ub].plot(x='commit',
                                           y='cpu_time',
                                           label='cpu time [ns]  ' + ub,
                                           style='-',
                                           marker="o",
                                           figsize=(10, 5))
    plt.savefig(ub + '.png', format='png')
