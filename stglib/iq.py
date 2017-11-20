from __future__ import division, print_function
import scipy.io as spio
import matplotlib.pyplot as plt

iqmat = spio.loadmat('/Volumes/Backstaff/field/1097/1097IQ_20170815_124439.mat')


# 1e6
# sample time is is microseconds since the beginning of the decade.
# for example, 556122600000000.000000 is 17.63453196347032 years
print('%f' % iqmat['FlowData_SampleTime'][7])

year = floor(iqmat['FlowData_SampleTime'][7]

print(iqmat['Data_Units'])
print(iqmat['FlowData_SampleTime'][1] - iqmat['FlowData_SampleTime'][0])
# %%
for k in iqmat:
    # if 'FlowData' in k:
    print(k)
for k in iqmat:
    if 'Profile' in k:
        print(k)
# %%

print(iqmat.keys())

# %%
plt.figure(figsize=(18,8))
plt.plot(iqmat['FlowData_Flow'])
# plt.ylim(-200,200)

plt.show()
