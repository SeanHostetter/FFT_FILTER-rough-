#-------------CONFIGURATION--------------------------------------------------
file_head = "coupon_board_into_scope_256x25sec_6GHz256x1sk_DownSampled.csv"
write_FFT = False
write_FFT_Clean = False
write_Data_Clean = False
graph_256 = True
Filter = "Next20"       #Next20, Zero, Average, none
#-------------CONFIGURATION--------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['agg.path.chunksize'] = 10000
plt.rcParams['figure.figsize'] = [15, 7]
plt.rcParams.update({'font.size' : 18})
file_name = "C:\\Users\\brian\\Desktop\\DataConversionThing\\DownSampledData\\" + file_head

def graph(dv):                      #function to graph passed in csv
    print(dv)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    x = dv[200000:350000,0]
    y = dv[200000:350000,1]
    ax.scatter(x,y)
    plt.title('Data')
    plt.ylabel('CV')
    plt.xlabel('DV')
    plt.show()


def file_len(fname):                 #function to calculate file length(for graph calibration)
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


import csv

##DELETE#LINE#1#####################
a_file = open(file_name, "r")
lines = a_file.readlines()
a_file.close()
del lines[0]
######################################
#t = np.arange(0,1.62760417e-8,dt)
dv = np.loadtxt(lines, delimiter=',')
length = file_len(file_name)
#graph(dv)

#CREATING VARIABLES FOR FFT
f = dv[0:,1]                           #output values
t = dv[0:,0]                           #time values
dt = 1/(6000000000 * 1024)             #time step

#FFT
n = len(t)
fhat = np.fft.fft(f,n)                      #compute FFT
PSD = fhat * np.conj(fhat) / n              #Power Spectrum
freq = (1/(dt*n)) * np.arange(n)            #create vector of powers at each frequency
L = np.arange(1,np.floor(n/2),dtype='int')  #only plot first half

#declaring plotting variables

PSD_clone = PSD
range_bot = 6.1e9 - 0.2e9
range_top = 6.1e9 + 0.2e9

tot_domain = freq[1],freq[-1]
tot_range = PSD[1],PSD[-1]
print("domain: ")
print(tot_domain)
print("range: ")
print(tot_range)
middle = 4e12
cutoff = 1.75e4
whole = False

if (whole):
  sv = np.where(freq > freq[1])[0][0]          #display entire FFT graph
  ev = np.where(freq > freq[-2])[0][0]         #^^^
  height = 1e4                                 #^^^
else:
  sv = np.where(freq > 6.1e9-10e9)[0][0]        #zoom in on 6ghz signal
  ev = np.where(freq > 6.1e9+10e9)[0][0]        #^^^
  height = 1e9                                #^^^

print("Starding index:", sv)
print("Ending index: ", ev)


def graph_fft(range_bot, range_top, freq, PSD, freq_clean, PSD_clean):
    fig, axs = plt.subplots(1, 1)

    plt.sca(axs)

    if (whole):
        plt.plot(freq[L], PSD[L], color='k', linewidth=0.5, label='Noisy')  # plot entire thing(first half)
        plt.plot(freq_clean[L], PSD_clean[L], color='c', linewidth=0.5, label='Clean')  # plot entire thing(first half)
        plt.xlim(freq[L[0]], freq[L[-1]])  # ^^^
    else:
        plt.plot(freq[sv:ev], PSD[sv:ev], color='k', linewidth=0.5, label='Noisy')  # plot selection
        plt.plot(freq_clean[sv:ev], PSD_clean[sv:ev], color='c', linewidth=0.5, label='Clean')  # plot entire thing(first half)
        plt.xlim(freq[sv], freq[ev])  # ^^^

    plt.vlines(range_bot,0,1e11,colors='k',linestyles='solid',label='')     #vertical lines to reference frequency ranges
    plt.vlines(range_top,0,1e11,colors='k',linestyles='solid',label='')

    plt.ylim(0, height)
    plt.legend()
    plt.show()

# FILTERING OUT 6ghz AND 4 HARMONICS----------------------------------------------


start_index = np.where(freq > range_bot)
end_index = np.where(freq > range_top)

freq_clean = freq.copy()
PSD_clean = PSD.copy()

get_rid_of_base = start_index[0][0] + np.argmax(PSD[start_index[0][0]:end_index[0][0]])
print(get_rid_of_base)

mult=2
# ---"Smart"-Filtering------------------------------------------------------------
if Filter == "Next20":
    for i in range(4):  # each loop filters a multiple of ~6ghz
        for j in range(20):
            get_rid_of = start_index[0][0] + np.argmax(PSD_clean[start_index[0][0]:end_index[0][0]])
            PSD_clean[get_rid_of] = PSD_clean[get_rid_of + 20]  # replaces local maxes with nearby data to even it out
            fhat[get_rid_of] = fhat[get_rid_of + 20]
            PSD_clean[-get_rid_of] = PSD_clean[-get_rid_of - 20]  # replaces local maxes with nearby data to even it out
            fhat[-get_rid_of] = fhat[-get_rid_of - 20]
            print(PSD_clean[get_rid_of])
        new_target = (mult * freq[get_rid_of_base])  # go after next harmonic
        start_index = np.where(freq > (new_target - 0.2e9))
        end_index = np.where(freq > (new_target + 0.2e9))
        mult *= 2
# ---Zero-Signal-Filtering--------------------------------------------------------
elif Filter == "Zero":
    for i in range(4):                   #each loop filters a multiple of ~6ghz
        get_rid_of = start_index[0][0] + np.argmax(PSD[start_index[0][0]:end_index[0][0]])
        PSD_clean[get_rid_of-20:get_rid_of+20] = 0                                     #replaces local maxes with minimum value in spectrum
        fhat[get_rid_of-20:get_rid_of+20] = fhat[354694]
        PSD_clean[-get_rid_of-20:get_rid_of+20] = 0                                    #replaces local maxes with minimum value in spectrum
        fhat[-get_rid_of-20:get_rid_of+20] = fhat[354694]
        print(get_rid_of)
        new_target = (mult*freq[get_rid_of_base])             #go after next harmonic
        start_index = np.where(freq > (new_target-0.2e9))
        end_index = np.where(freq > (new_target+0.2e9))
        mult*=2
# ---Average-Filtering------------------------------------------------------------
elif Filter == "Average":
    for i in range(4):  # each loop filters a multiple of ~6ghz
        for j in range(20):
            get_rid_of = start_index[0][0] + np.argmax(PSD_clean[start_index[0][0]:end_index[0][0]])
            PSD_clean[get_rid_of] = np.mean(PSD_clean[get_rid_of + 20:get_rid_of + 40])  # replaces local maxes with nearby data to even it out
            fhat[get_rid_of] = np.mean(fhat[get_rid_of + 20:get_rid_of + 40])
            PSD_clean[-get_rid_of] = np.mean(PSD_clean[-get_rid_of - 40:-get_rid_of - 20])  # replaces local maxes with nearby data to even it out
            fhat[-get_rid_of] = np.mean(fhat[-get_rid_of - 40:-get_rid_of - 20])
            print(PSD_clean[get_rid_of])
        new_target = (mult * freq[get_rid_of_base])  # go after next harmonic
        start_index = np.where(freq > (new_target - 0.2e9))
        end_index = np.where(freq > (new_target + 0.2e9))
        mult *= 2
# --------------------------------------------------------------------------------

ffilt = np.fft.ifft(fhat)  # reverse the fft, returns ffilt which is the filtered data

graph_fft(range_bot=range_bot,range_top=range_top, freq=freq,PSD=PSD,freq_clean=freq_clean,PSD_clean=PSD_clean)

#DISPLAYING THE FILTERED DATA

fig,axs = plt.subplots(1,1)

plt.sca(axs)
plt.plot(t,f,color='k',linewidth = 0.1,label='Filtered')
plt.plot(t,ffilt,color='c',linewidth = 0.1,label='Filtered')
if (graph_256):
    plt.xlim(t[10000],t[100000])
else:
    plt.xlim(t[100000], t[300000])
plt.legend()
plt.show()


if write_FFT:
    print("writing FFT to file...")
    fft_file = file_name[:file_name.find('.csv')] + '_FFT_' + file_name[file_name.find('.csv'):]
    a = [freq,PSD]
    with open(fft_file,"w+") as fft_csv:
        csvWriter = csv.writer(fft_csv,delimiter=',')
        csvWriter.writerows(a)

if write_FFT_Clean:
    print("writing filtered FFT to file...")
    fft_clean_file = file_name[:file_name.find('.csv')] + '_FFT_CLEAN_' + file_name[file_name.find('.csv'):]
    a = [freq_clean,PSD_clean]
    with open(fft_file,"w+") as fft_csv:
        csvWriter = csv.writer(fft_csv,delimiter=',')
        csvWriter.writerows(a)

if write_FFT:
    print("writing Filtered Dataset to file...")
    clean_file = file_name[:file_name.find('.csv')] + '_CLEAN_' + file_name[file_name.find('.csv'):]
    a = [t,ffilt]
    with open(fft_file,"w+") as fft_csv:
        csvWriter = csv.writer(fft_csv,delimiter=',')
        csvWriter.writerows(a)