LOWER_X = 2200
UPPER_X = 2500

import pyfits
import pylab as plt
import sys
import numpy as np
from scipy.ndimage.filters import median_filter

colours = ['r', 'b', 'g', 'y', 'k']
fl_std = []
for f_idx, fn in enumerate(sys.argv[1:]):
  f = pyfits.open(fn)
  sn = []
  fl = []
  fi = []
  for idx, row in enumerate(f[0].data):
    x = range(len(row))
    y = median_filter(row, size=100)
    y_norm = row/y
    y_99_dot_5_percentile = np.percentile(row, 99.5)
    y_norm_std = np.nanstd((row/y)[LOWER_X:UPPER_X])
    sn.append(np.nanmean(y_norm)/y_norm_std)
    fl.append(y_99_dot_5_percentile)
    fi.append(idx+1)
    if f_idx == 0:
      fl_std.append(y_99_dot_5_percentile)

  plt.plot(fl_std, sn, colours[f_idx] + '.', label=fn)
  plt.xlabel("99.5th percentile flux (nopt)")
  plt.ylabel("sn")
  plt.legend()

plt.show()

# as a f of fibre number
for f_idx, fn in enumerate(sys.argv[1:]):
  f = pyfits.open(fn)
  sn = []
  fl = []
  fi = []
  for idx, row in enumerate(f[0].data):
    x = range(len(row))
    y = median_filter(row, size=100)
    y_norm = row/y
    y_99_dot_5_percentile = np.percentile(row, 99.5)
    y_norm_std = np.nanstd((row/y)[LOWER_X:UPPER_X])
    sn.append(np.nanmean(y_norm)/y_norm_std)
    fl.append(y_99_dot_5_percentile)
    fi.append(idx+1)

  plt.plot(fi, sn, colours[f_idx] + '-', label=fn)
  plt.legend()

plt.ylabel("sn")
plt.xlabel("fibre number")
plt.show()

