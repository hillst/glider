import matplotlib as mpl
import numpy as np
mpl.use('agg')
from viz.logo import plot_logo
import matplotlib.pyplot as plt


motif = np.asarray([[.25,.25,.25,.25],[.1,.1,.1,.7], [.3,.29,.4,0.01], [.97,0.01,0.01,0.01],  [0,0,0,1], [.25,.25,.25,.25], [.25,.25,.25,.25] , [.25,.25,.25,.25]])   
ax = plot_logo(motif)
import seaborn as sb
sb.despine(trim=True)
plt.savefig("test.png")
