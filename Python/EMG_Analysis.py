import numpy as np
import matplotlib.pyplot as plt
import bioread

data = bioread.read_file('/Users/boramert/Desktop/Yüksek Lisans/Data_Results/Data/EMG_Data/Acq/EMG_001.acq')
emg = data.channels[0].data
dyno = data.channels[1].data
t = data.channels[0].time_index

plt.figure(dpi=1200)
plt.plot(t,emg, linewidth = 0.5)
plt.plot(t,dyno*0.01, linewidth = 0.5)
plt.show()