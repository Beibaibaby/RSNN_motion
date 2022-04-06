import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

mat = scipy.io.loadmat('motionEnergy.mat')
data= mat['motionEnergy']
time_step,contrast_=data.shape[0],data.shape[1]

low_contrast = data[:,0,:]
high_contrast = data[:,contrast_-1,:]
low_contrast_ratio= np.zeros(low_contrast[:,1].shape)
low_contrast_ratio[9:]=(low_contrast[9:,1])/(low_contrast[9:,0]+low_contrast[9:,1])
low_contrast_ratio[0:9]=low_contrast_ratio[9]
print(low_contrast_ratio)

def objective(x, a, b):
	return a + b**x

x=[0 + x*(low_contrast_ratio.shape[0]-0)/low_contrast_ratio.shape[0] for x in range(low_contrast_ratio.shape[0])]
x=np.asarray(x)
popt, _ = curve_fit(objective, x, low_contrast_ratio)
a, b = popt
mymodel = np.poly1d(np.polyfit(x, low_contrast_ratio,5))

fig, ax = plt.subplots()

ccc=['blue','red','green','brown','purple']
#plt.plot(time,smoothing(ts_transformed.T[0],50),'blue',label='pc1-'+str(pca.explained_variance_ratio_[0]),markersize=3)
ax.plot(low_contrast[:,0],ccc[1],label='Left',markersize=3)
ax.plot(low_contrast[:,1],ccc[0],label='Right',markersize=3)

#plt.title('Before stress: PC (smoothing =50)')
ax.set_xlabel('times',fontsize=14)
ax.set_ylabel('Energy',fontsize=14)

ax2=ax.twinx()
ax2.plot(low_contrast_ratio,ccc[2],label='Right_ratio',markersize=3,alpha=0.8)
ax2.plot(mymodel(x),ccc[3],label='Fitting Line',markersize=3)
ax2.set_ylabel("percent",fontsize=14)
fig.legend()
plt.title('Motion Energy-Low Contrast')
plt.show()


fig.savefig('./test.png')
plt.show()
