import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


cc=0.12
Je = 0.3103
Jot = -0.007
Jop = -0.046
contrast = 0.05
sum_plot=0
#3 decision with smoothing
cov=True
check=False
check_reverse=False
iter=1

mat = scipy.io.loadmat('motionEnergy.mat')
data= mat['motionEnergy']
time_step,contrast_=data.shape[0],data.shape[1]

low_contrast = data[:,0,:]
high_contrast = data[:,contrast_-1,:]
low_contrast_ratio= np.zeros(low_contrast[:,1].shape)
low_contrast_ratio[9:]=(low_contrast[9:,1])/(low_contrast[9:,0]+low_contrast[9:,1])
low_contrast_ratio[0:9]=low_contrast_ratio[9]

high_contrast_ratio= np.zeros(high_contrast[:,1].shape)
high_contrast_ratio[9:]=(high_contrast[9:,1])/(high_contrast[9:,0]+high_contrast[9:,1])
high_contrast_ratio[0:9]=high_contrast_ratio[9]


def objective(x, R_max, c_1,c_2):
	return R_max*(c_1**x/(c_2**x+c_1**x))+0.2

x=[0 + x*(low_contrast_ratio.shape[0]-0)/low_contrast_ratio.shape[0] for x in range(low_contrast_ratio.shape[0])]
x=np.asarray(x)
if contrast == 0.05:
    popt, _ = curve_fit(objective, x, low_contrast_ratio)
else:
    popt, _ = curve_fit(objective, x, high_contrast_ratio)
R_max, c_1,c_2 = popt

def convolve(x,covk):
    r=np.zeros(x.shape)
    for i in range(len(x)):
        xx=0
        for j in range(len(covk)):
            xx=x[i]


def input_function(R_max, c_1,c_2,cprime):
    ratio = objective([0 + x for x in range(350)], R_max, c_1,c_2)
    right = cprime**cc*ratio*0.8*60
    left = cprime**cc * (1 - ratio) * 0.8 * 60
    if cov==True:
      cov_c = np.asarray([0.1, 0.097, 0.089, 0.083, 0.08, 0.076, 0.088, 0.089, 0.092, 0.1])
    #left = cprime * (1 - ratio) * 0.9 * 60

      right = np.convolve(right, cov_c, 'same')

      for i in range(5):
          right[i] = right[5]
          right[len(right) - i - 1] = right[len(right) - 6]

      #for i in range(80):
      #    left[i] = left[0]

      #for i in range(100):
      #    left[60 + i] = left[200]

    up = cprime**cc*0.1*(np.zeros(left.shape)+1)*60
    down = cprime**cc*0.1*(np.zeros(left.shape)+1)*60
    #print(right)
    return left,right,up,down

hue = ['orange', 'red', 'blue', 'green']

def tuningc(degree, cprime):
    realdegree = 0
    r0 = 20
    r1 = 0
    r2 = 100
    sigma = 40
    v = r0 + cprime * (-r1 + r2 * math.e ** (-(degree - realdegree) ** 2 / (sigma ** 2)))
    return v

left,right,up,down =input_function(R_max, c_1,c_2,0.05)
#print(input_function(R_max,c_1,c_2,0.99))
plt.figure

fig, axs = plt.subplots(3)
fig.set_size_inches(10.5, 18.5)
axs[0].plot([0 + x for x in range(350)], left, color=hue[0], label='left')
axs[0].plot([0 + x for x in range(350)], up, color=hue[1],label='up')
axs[0].plot([0 + x for x in range(350)], right, color=hue[2] ,label='right')
axs[0].plot([0 + x for x in range(350)], down, color=hue[3] ,label='down')

axs[0].set(xlabel='x-label', ylabel='y-label')
# plt.ylim(0,20)
# plt.xlim(-600,750)
# plt.text(-125, 16, 'Threshold')
#plt.legend()
#plt.savefig('./figs/input.png')
#plt.show()

a = 270
b = 108
d = 0.154
gamma = 0.641
taus = 100 * 10 ** -3
tauampa = 2 * 10 ** -3
#Je = 0.3103
#Jot = -0.0185
#Jop = -0.0185


(J11,J12,J13,J14) = (Je,Jot,Jop,Jot)
(J21,J22,J23,J24) = (Jot,Je,Jot,Jop)
(J31,J32,J33,J34) = (Jop,Jot,Je,Jot)
(J41,J42,J43,J44) = (Jot,Jop,Jot,Je)


Jext = 5.2 * 10 ** -4
I0 = 0.3255
stdnoise = 0.5
mu0 = 30
dt = 0.1 * 10 ** -3

def objective(x, R_max, c_1,c_2):
	return R_max*(c_1**x/(c_2**x+c_1**x))+0.2

def H(xi):
    return (a * xi - b) / (1 - np.exp(-d * (a * xi - b)))


starttime = -0.5
endtime = 2
steps = int(abs(starttime - endtime) / dt)
time = np.linspace(starttime, endtime, steps)
time = time





def experiment(cprime):
    (H1, H2, S1, S2) = (np.zeros(steps + 1), np.zeros(steps + 1),
                        np.zeros(steps + 1), np.zeros(steps + 1))
    (H3, H4, S3, S4) = (np.zeros(steps + 1), np.zeros(steps + 1),
                        np.zeros(steps + 1), np.zeros(steps + 1))
    ini = 0
    S1[0] = ini
    S2[0] = ini
    S3[0] = ini
    S4[0] = ini

    H1[0] = 0
    H2[0] = 0
    H3[0] = 0
    H4[0] = 0

    (Inoise1, Inoise2, Inoise3, Inoise4) = (
    np.zeros(steps + 1), np.zeros(steps + 1), np.zeros(steps + 1), np.zeros(steps + 1))
    left,right,up,down= input_function(R_max, c_1, c_2, cprime)

    for (index, t) in enumerate(time):
        if index < 5000:
            mu1 = left[0]
            mu2 = up[0]
            mu3 = right[0]
            mu4 = down[0]
        elif (index>=5000 and index<8500 ):
            ik=int(index/10)-500
            mu1 = left[ik]
            mu2 = up[ik]
            mu3 = right[ik]
            mu4 = down[ik]
        else:
            mu1 = left[349]
            mu2 = up[349]
            mu3 = right[349]
            mu4 = down[349]
        Inoise1[index + 1] = 0

        Inoise2[index + 1] = 0
        Inoise3[index + 1] = 0
        Inoise4[index + 1] = 0
        if t > 0 :
                sks = 0

           #if index>=5000 and index<5500:

          #      x1 = J11 * S1[index] + J12 * S2[index] + J13 * S3[index] + J14 * S4[index] + I0 + Jext * mu1 \
           #          + Inoise1[index]

           #     x2 = x1
            #    x3 = J31 * S1[index] + J32 * S2[index] + J33 * S3[index] + J34 * S4[index] + I0 + Jext * mu3 \
            #         + Inoise3[index]
            #    x4 = x1




                x1 = J11 * S1[index] + J12 * S2[index] + J13 * S3[index] + J14 * S4[index] + I0 + Jext * mu1 \
                     + Inoise1[index]

                x2 = J21 * S1[index] + J22 * S2[index] + J23 * S3[index] + J24 * S4[index] + I0 + Jext * mu2 \
                     + Inoise2[index]

                x3 = J31 * S1[index] + J32 * S2[index] + J33 * S3[index] + J34 * S4[index] + I0 + Jext * mu3 \
                     + Inoise3[index]

                x4 = J41 * S1[index] + J42 * S2[index] + J43 * S3[index] + J44 * S4[index] + I0 + Jext * mu4 \
                     + Inoise4[index]


        else :
            x1 = J11 * S1[index] + J12 * S2[index] + J13 * S3[index] + J14 * S4[index] + I0

            x2 = J21 * S1[index] + J22 * S2[index] + J23 * S3[index] + J24 * S4[index] + I0

            x3 = J31 * S1[index] + J32 * S2[index] + J33 * S3[index] + J34 * S4[index] + I0

            x4 = J41 * S1[index] + J42 * S2[index] + J43 * S3[index] + J44 * S4[index] + I0


        H1[index + 1] = H(x1)
        H2[index + 1] = H(x2)
        H3[index + 1] = H(x3)
        H4[index + 1] = H(x4)

        S1[index + 1] = S1[index] + dt * (-S1[index] / taus + (1 \
                                                               - S1[index]) * gamma * H1[index])

        S2[index + 1] = S2[index] + dt * (-S2[index] / taus + (1 \
                                                               - S2[index]) * gamma * H2[index])

        S3[index + 1] = S3[index] + dt * (-S3[index] / taus + (1 \
                                                               - S3[index]) * gamma * H3[index])

        S4[index + 1] = S4[index] + dt * (-S4[index] / taus + (1 \
                                                               - S4[index]) * gamma * H4[index])
    print(S1[5002],S2[5002],S3[5002],S4[5002])
    return (H1[1:], H2[1:], H3[1:], H4[1:])


def slided(data):
    timestep = 5 * 10 ** -3
    slided_data = []
    for (index, value) in enumerate(data):
        if index % int(timestep / dt) == 0:
            slided_data.append(value)
    return slided_data


def smoothing(data):
    length = len(data)
    smoothed_data = np.zeros(length)
    width = 1000
    for i in range(length):
        if length - (i + 1) < width:
            smoothed_data[i] = np.average(data[i:])
        else:
            smoothed_data[i] = np.average(data[i:i + width])
    return smoothed_data



hue = ['orange', 'red', 'blue', 'green']


for cprime in [contrast]:
  if sum_plot==0:
    from scipy.special import softmax

    left_con = []
    up_con = []
    right_con = []
    down_con = []
    for i in range(iter):

        print(i)
        result = experiment(cprime)
        result = np.asarray(result)


        print(result.shape)
        result1=result[:,:4000]
        result2=result[:,4000:]

        for idd in range(4):
             result2[idd]=smoothing(result2[idd])
        result=np.concatenate((result1, result2), axis=1)
        if check==True:
            if contrast == 0.99:
                if result[0, 24000] > 5:
                    continue
        if check_reverse == True:
            if contrast == 0.99:
                if result[0, 24000] <=5:
                    continue

        print(result.shape)
        hue = ['orange', 'red', 'blue', 'green']
        time=time+0.1
        axs[1].set_xlim(-500, 500)
        if (i == 0):
            axs[1].plot(time * 1000, result[0], color=hue[0], label='left')
            axs[1].plot(time * 1000, result[1], color=hue[1], label='up')
            axs[1].plot(time * 1000, result[2], color=hue[2], label='right')
            axs[1].plot(time * 1000, result[3], color=hue[3], label='down')
        else:

            axs[1].plot(time * 1000, result[0], color=hue[0], label='left')
            axs[1].plot(time * 1000, result[1], color=hue[1], label='up')
            axs[1].plot(time * 1000, result[2], color=hue[2], label='right')
            axs[1].plot(time * 1000, result[3], color=hue[3], label='down')
        result_softmax = softmax(result, axis=0)

        left_con.append(result_softmax[0])
        up_con.append(result_softmax[1])
        right_con.append(result_softmax[2])
        down_con.append(result_softmax[3])

    axs[0].set(xlabel='Time(ms)', ylabel='Motion Energy')

    #plt.xlabel('Time(ms)')
    #plt.ylabel('Firing rate(Hz)')
    axs[1].set(xlabel='Time(ms)', ylabel='Firing Rate')
    axs[2].set_xlim(-500, 500)
    #plt.title('Contrast-' + str(cprime) + ' Iteration-' + str(iter) + ' Firing Rate')
    # plt.text(-125, 16, 'Threshold')
    axs[2].set(xlabel='Time(ms)', ylabel='Probability')







    left_con = np.asarray(left_con)
    left_mean = left_con.mean(axis=0)
    left_var = left_con.var(axis=0)

    right_con = np.asarray(right_con)
    right_mean = right_con.mean(axis=0)
    right_var = right_con.var(axis=0)

    up_con = np.asarray(up_con)
    up_mean = up_con.mean(axis=0)
    up_var = up_con.var(axis=0)

    down_con = np.asarray(down_con)
    down_mean = down_con.mean(axis=0)
    down_var = down_con.var(axis=0)
    axs[2].plot(time * 1000, left_mean, color=hue[0], label='left')
    axs[2].plot(time * 1000, up_mean, color=hue[1], label='up')
    axs[2].plot(time * 1000, right_mean, color=hue[2], label='right')
    axs[2].plot(time * 1000, down_mean, color=hue[3], label='down')
    axs[2].fill_between(time * 1000, left_mean - left_var, left_mean + left_var, color=hue[0], alpha=0.3)
    axs[2].fill_between(time * 1000, up_mean - up_var, up_mean + up_var, color=hue[1], alpha=0.3)
    axs[2].fill_between(time * 1000, right_mean - right_var, right_mean + right_var, color=hue[2], alpha=0.3)
    axs[2].fill_between(time * 1000, down_mean - down_var, down_mean + down_var, color=hue[3], alpha=0.3)
    # plt.plot(time * 1000, 15 * np.ones(steps))
    #.xlabel('Time(ms)')
    #plt.ylabel('Pro')
    #axs[2].set(xlabel='Time(ms)', ylabel='Time(ms)')
    #plt.title('Contrast-' + str(cprime) + ' Iteration-' + str(iter) + ' Decision Making(Smoothing)')
    # plt.text(-125, 16, 'Threshold')
    #plt.legend()
    #plt.savefig('./figs/psychophysics_' + str(cprime) + '_smoothing.png')
    #plt.show()
    plt.legend()
    plt.savefig('./figs/allthree_' + str(cprime) + '.png')
    plt.show()
  elif sum_plot== 1:
    for i in range(iter):
        print(i)
        result = experiment(cprime)
        result = np.asarray(result)
        if(i==0):
            sum=np.zeros((4,20000))
        else:
            sum= sum +result
        hue = ['orange', 'red', 'blue', 'green']


    sum=sum/iter
    plt.plot(time * 1000, smoothing(sum[0]), color=hue[0], label='left')
    plt.plot(time * 1000, smoothing(sum[1]), color=hue[1], label='up')
    plt.plot(time * 1000, smoothing(sum[2]), color=hue[2], label='right')
    plt.plot(time * 1000, smoothing(sum[3]), color=hue[3], label='down')
    plt.xlabel('Time(ms)')
    plt.ylabel('Firing rate(Hz)')
    # plt.ylim(0,20)
    # plt.xlim(-600,750)
    # plt.text(-125, 16, 'Threshold')
    plt.legend()
    plt.savefig('./figs/low.png')
    plt.show()
  elif sum_plot == 2:
      from scipy.special import softmax
      for i in range(iter):
          print(i)
          result = experiment(cprime)
          result = np.asarray(result)
          result = softmax(result, axis=0)
          hue = ['orange', 'red', 'blue', 'green']
          if (i == 0):
              plt.plot(time * 1000, smoothing(result[0]), color=hue[0], label='left')
              plt.plot(time * 1000, smoothing(result[1]), color=hue[1], label='up')
              plt.plot(time * 1000, smoothing(result[2]), color=hue[2], label='right')
              plt.plot(time * 1000, smoothing(result[3]), color=hue[3], label='down')
          else:

              plt.plot(time * 1000, smoothing(result[0]), color=hue[0])
              plt.plot(time * 1000, smoothing(result[1]), color=hue[1])
              plt.plot(time * 1000, smoothing(result[2]), color=hue[2])
              plt.plot(time * 1000, smoothing(result[3]), color=hue[3])
#plt.plot(time * 1000, 15 * np.ones(steps))
      plt.xlabel('Time(ms)')
      plt.ylabel('Pro')
      # plt.text(-125, 16, 'Threshold')
      plt.legend()
      plt.savefig('./figs/low_probablity.png')
      plt.show()
  elif sum_plot == 3:
      from scipy.special import softmax
      left_con = []
      up_con=[]
      right_con=[]
      down_con=[]

      for i in range(iter):
          print(i)
          result = experiment(cprime)
          result = np.asarray(result)
          result = softmax(result, axis=0)
          hue = ['orange', 'red', 'blue', 'green']
          left_con.append(result[0])
          up_con.append(result[1])
          right_con.append(result[2])
          down_con.append(result[3])

      left_con=np.asarray(left_con)
      left_mean=left_con.mean(axis=0)
      left_var = left_con.var(axis=0)

      right_con=np.asarray(right_con)
      right_mean=right_con.mean(axis=0)
      right_var = right_con.var(axis=0)

      up_con = np.asarray(up_con)
      up_mean = up_con.mean(axis=0)
      up_var = up_con.var(axis=0)

      down_con = np.asarray(down_con)
      down_mean = down_con.mean(axis=0)
      down_var = down_con.var(axis=0)
      plt.plot(time * 1000, left_mean, color=hue[0], label='left')
      plt.plot(time * 1000, up_mean, color=hue[1], label='up')
      plt.plot(time * 1000, right_mean, color=hue[2], label='right')
      plt.plot(time * 1000, down_mean, color=hue[3], label='down')
      plt.fill_between(time * 1000, left_mean-left_var, left_mean+left_var, color=hue[0],alpha=0.3)
      plt.fill_between(time * 1000, up_mean - up_var, up_mean + up_var, color=hue[1], alpha=0.3)
      plt.fill_between(time * 1000, right_mean - right_var, right_mean + right_var, color=hue[2], alpha=0.3)
      plt.fill_between(time * 1000, down_mean - down_var, down_mean + down_var, color=hue[3], alpha=0.3)
      # plt.plot(time * 1000, 15 * np.ones(steps))
      plt.xlabel('Time(ms)')
      plt.ylabel('Pro')
      plt.title('Contrast-'+str(cprime)+' Iteration-'+str(iter)+' Decision Making')
      # plt.text(-125, 16, 'Threshold')
      plt.legend()
      plt.savefig('./figs/psychophysics_'+str(cprime)+'.png')
      plt.show()
  elif sum_plot == 4:
      from scipy.special import softmax
      left_con = []
      up_con=[]
      right_con=[]
      down_con=[]

      for i in range(iter):
          print(i)
          result = experiment(cprime)
          result = np.asarray(result)
          result = softmax(result, axis=0)
          hue = ['orange', 'red', 'blue', 'green']
          left_con.append(smoothing(result[0]))
          up_con.append(smoothing(result[1]))
          right_con.append(smoothing(result[2]))
          down_con.append(smoothing(result[3]))

      left_con=np.asarray(left_con)
      left_mean=left_con.mean(axis=0)
      left_var = left_con.var(axis=0)

      right_con=np.asarray(right_con)
      right_mean=right_con.mean(axis=0)
      right_var = right_con.var(axis=0)

      up_con = np.asarray(up_con)
      up_mean = up_con.mean(axis=0)
      up_var = up_con.var(axis=0)

      down_con = np.asarray(down_con)
      down_mean = down_con.mean(axis=0)
      down_var = down_con.var(axis=0)
      plt.plot(time * 1000, left_mean, color=hue[0], label='left')
      plt.plot(time * 1000, up_mean, color=hue[1], label='up')
      plt.plot(time * 1000, right_mean, color=hue[2], label='right')
      plt.plot(time * 1000, down_mean, color=hue[3], label='down')
      plt.fill_between(time * 1000, left_mean-left_var, left_mean+left_var, color=hue[0],alpha=0.3)
      plt.fill_between(time * 1000, up_mean - up_var, up_mean + up_var, color=hue[1], alpha=0.3)
      plt.fill_between(time * 1000, right_mean - right_var, right_mean + right_var, color=hue[2], alpha=0.3)
      plt.fill_between(time * 1000, down_mean - down_var, down_mean + down_var, color=hue[3], alpha=0.3)
      # plt.plot(time * 1000, 15 * np.ones(steps))
      plt.xlabel('Time(ms)')
      plt.ylabel('Pro')
      plt.title('Contrast-'+str(cprime)+' Iteration-'+str(iter)+' Decision Making(Smoothing)')
      # plt.text(-125, 16, 'Threshold')
      plt.legend()
      plt.savefig('./figs/psychophysics_'+str(cprime)+'_smoothing.png')
      plt.show()