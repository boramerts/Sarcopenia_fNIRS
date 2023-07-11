import numpy as np
import matplotlib.pyplot as plt
import bioread

import scipy.io
emg = scipy.io.loadmat('EMG_001.mat')

data = bioread.read_file('Acq/EMG_001.acq')

# plt.subplot(2,1,1)
# plt.plot(emg['data'][:,0])

# plt.subplot(2,1,2)
# plt.plot(emg['data'][:,1])

# plt.show()

print(data)