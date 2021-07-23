""" se definen varias funciones"""
import numpy as np
from numpy.core.numeric import False_
from scipy.signal import filtfilt, butter, lfilter, find_peaks, detrend, hilbert
from scipy.integrate import cumtrapz, trapezoid
from scipy.interpolate import interp1d
from newclasses import *
from datetime import date
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import dialog as dlg
from PyQt5 import QtCore, QtGui, QtWidgets
import copy
#def median_filter:

def remove_dup(a):
    # quita elementos repetidos en un vector

    #Argumentos:
    #   a:<dict> lista con elementos repetidos
    r = list(set(tuple(a)))
    #Return:
    #   r:<list> lista sin elementos repertidos
    return r
    
def datetoday(datev):
    #convierte un vetor con fechas a un vector con enteros, cada entero son los días desde cada fecha a la fecha más pequeña

    #Argumentos:
    #   datev: <list> vector con fechas (elementos date)
    m = min(datev)
    r = [(k-m).days for k in datev]
    #Return:
    #   r: <list> vector con enteros (días)
    return r


def sortbyrow(a,f):
    #esta función no sirve
    b = [[] for k in a]
    r = copy.deepcopy(a[f])
    order = []
    for k in r:
        order.append(np.argmin(r))
        r[np.argmin(r)]=max(r)
    print(order)
    return [[k[j] for j in order] for k in a]

def strtodate(s,f=False):
    #convierte fecha en formato "AA-MM-DD" a formato "DD/MM/AA" o a objeto date

    #Argumentos:
    #   s: <str> fecha "AA-MM-DD"
    #   f: <bool> define si se devuelve "DD/MM/AA" o un objeto date
    a = s.split('/')
    if f:
        r = '/'.join(list(reversed(str(date(int(a[2]),int(a[1]),int(a[0]))).split('-'))))
        #Return:
        #   r:<str> fecha en formato "DD/MM/AA"
        return r
    else:
        r=date(int(a[2]),int(a[1]),int(a[0]))
        #Return:
        #   r:<str> fecha en objeto date
        return r


def calculateAge(birthDate):
    #calcula la edad actual con base en una fecha de nacimiento

    #Argumentos:
    #   birthDate: <datetime.date> objeto date con la fecha de nacimiento
	days_in_year = 365.2425	
	age = int((date.today() - birthDate).days / days_in_year)
    #Return
    #   age: <int> edad en años
	return age

def bandpass_filter(data, f1, f2, fs, order):
    #aplica un filtro pasabanda (pasabajas en serie con pasaaltas)

    #Argumentos:
    #   data: <numpy.ndarray> señal que se filtra
    #   f1: <float> frecuencia de corte menor
    #   f2: <float> frecuencia de corte mayor
    #   fs: <float> frecuencia de muestreo
    #   order: <int> orden del filtro
    filtlow = lowpass_filter(data, f2, fs, order)
    filthigh = highpass_filter(filtlow, f1, fs, order)
    #Return:
    #   filthigh: <numpy.ndarray> señal filtrada
    return filthigh

def lowpass_filter(data, f, fs, order):
    #aplica un filtro butterworth pasabajas
    
    #Argumentos:
    #   data: <numpy.ndarray> señal que se filtra
    #   f: <float> frecuencia de corte
    #   fs: <float> frecuencia de muestreo
    #   order: <int> orden del filtro
    b, a = butter(order, f / (fs/2), btype='lowpass')
    r = filtfilt(b, a, data)
    #Return:
    #   r: <numpy.ndarray> señal filtrada
    return r

def highpass_filter(data, f, fs, order):
    #aplica un filtro butterworth pasaaltas
    
    #Argumentos:
    #   data: <numpy.ndarray> señal que se filtra
    #   f: <float> frecuencia de corte
    #   fs: <float> frecuencia de muestreo
    #   order: <int> orden del filtro
    b, a = butter(order, f / (fs/2), btype='highpass')
    r = filtfilt(b, a, data)
    #Return:
    #   r: <numpy.ndarray> señal filtrada
    return r

def get_avg(data):
    #calcula la media de un vector

    #Argumentos:
    #   data: <numpy.ndarray> o <list> arreglo con valores
    r = float(np.mean(data))
    #Return:
    #   r: <float> promedio
    return r

def get_vpp(data):
    #calcula valor pico-pico de una señal

    #Argumentos:
    #   data: <numpy.ndarray> o <list> arreglo con valores
    r = float(abs(np.max(data) - np.min(data)))
    #Return:
    #   r: <float> valor pico-pico
    return r

def fft(data, fs): 
    #calcula la transformada rápida de fourier de una señal

    #Argumentos:
    #   data: <numpy.ndarray>  señal
    #   fs: <float> frecuencia de muestreo
    mag = abs(np.fft.fft(data))
    freq = np.linspace(0, fs / 2, int(len(data) / 2))
    #devuelve la fft  [frecuencia, magnitud]
    f = freq[1:int(len(data) / 2)].copy()
    m = mag[1:int(len(data) / 2)].copy()
    #Return:
    #   f: <numpy.ndarray> vector de frecuencias de la transformada
    #   m: <numpy.ndarray> magnitudes de la transformada con respecto a f
    return [f, m]

def segment_ppg(rc, irc, fs, time, fnirs=True):
    #Segmenta PPG o fNIRS, por picos sistólicos

    #Argumentos
    #   rc: <numpy.ndarray> señal de corriente en el rojo
    #   irc: <numpy.ndarray> señal de corriente en el infrarrojo
    #   fs: <float> frecuencia de muestreo
    filtlow = lowpass_filter(irc, 5, fs, 3)
    filtered = highpass_filter(filtlow, 0.5, fs, 3)
    if fnirs:
        #si es fnirs se invierte para obtener atenuación y se encuentran los cruces por cero de negativo a positivo
        locs = np.array(get_zeros(-filtered,'positive'))
    else:
        #si es PPG se encuentran los cruces por cero de negativo a positivo
        locs = np.array(get_zeros(filtered,'positive'))
    segmentsI = [(irc[locs[i]:locs[i+1]]) for i in range(0,len(locs)-1)]
    segmentsR = [(rc[locs[i]:locs[i+1]]) for i in range(0,len(locs)-1)] 
    times =[time[locs[i]:locs[i+1]] for i in range(0,len(locs)-1)] 
    if fnirs:
        #si es fnirs se invirte para calcular con atenuación
        segmentsF = [(-filtered[locs[i]:locs[i+1]]) for i in range(0,len(locs)-1)]
    else:
        segmentsF = [(filtered[locs[i]:locs[i+1]]) for i in range(0,len(locs)-1)]
    print(len(segmentsF))
    if fnirs:
        segmentsF.append(-filtered[locs[-1]:])
    else:
        segmentsF.append(filtered[locs[-1]:])
    print(len(segmentsF))
    locs = [fstmax(segmentsF[i])+locs[i] for i in range(0,len(segmentsF)) if type(fstmax(segmentsF[i]))!=None]
    segmentsI = [(irc[locs[i]:locs[i+1]]) for i in range(0,len(locs)-1)]
    segmentsR = [(rc[locs[i]:locs[i+1]]) for i in range(0,len(locs)-1)] 
    times =[time[locs[i]:locs[i+1]] for i in range(0,len(locs)-1)]
    #Return:
    #  times: <list> arreglo con los segmentos de tiempo que corresponden a los segmentos de señal
    #  segmentsR: <list> arreglo con los segmentos en el Rojo
    #  segmentsIR: <list> arreglo con los segmentos en el Infrarrojo
    #  locs: <list> arreglo con los índices donde ocurren los picos sistólicos  
    return times, segmentsR, segmentsI, locs

def fstmax(s):
    try:
        return min([i for i in range(1,len(s)-1) if s[i]>s[i-1] and s[i]>s[i+1]])
    except:
        return None

def integral(data,fs):
    #calcula la integral de una señal

    #Argumentos:
    #   data: <numpy.ndarray> señal para integrar
    #   fs: <float> frecuencia de muestreo
    time = (1/fs)*np.arange(len(data))
    i = cumtrapz(data, time, initial=0)
    #Return:
    #   i: <numpy.ndarray> señal integrada
    return i

def derivative(data):
    #calcula la derivada de una señal

    #Argumentos:
    #   data: <numpy.ndarray> señal para derivar
    d =np.array(data[1:])-np.array(data[:-1])
    #Return:
    #   d: <numpy.ndarray> señal derivada
    return d

def adj_mean(data):
    #calcula el promedio entre elementos adyacentes de un arreglo

    #Argumentos:
    #   data: <numpy.ndarray> arreglo
    r = (data[1:]+data[:-1])/2
    #Return:
    #   r: <numpy.ndarray> arreglo con los promedios sucesivos
    return r

def get_zeros(data,dir):
    #localiza los ceros de una señal

    #Argumentos: 
    #   data: <numpy.ndarray> o <list> señal para localizar cruces por cero
    #   dir: <str>, si es 'positive' localiza los cruces de - a +, si es 'negative' loccaliza los de + a -, en otro caso localiza ambos 
    if (dir=='positive'):
        locs = [k for k in range(1,len(data)) if ((data[k]>0) and (data[k-1]<0))]
    elif (dir=='negative'):
        locs = [k for k in range(1,len(data)) if ((data[k]<0) and (data[k-1]>0))]
    else:
        locs = [k for k in range(1,len(data)) if (((data[k]>0) and (data[k-1]<0)) or ((data[k]<0) and (data[k-1]>0)))]
    #Return:
    #   locs: <list> arreglo con los índices donde ocurren los cruces por cero
    return locs

def get_peaks(data,fs,threshold):
    #localiza los picos de una señal

    #Argumentos:
    #   data: <numpy.ndarray> o <list> señal para encontrar picos
    #   fs: <float> frecuencia de muestreo
    #   threshold: <float> umbral para contar el pico como válido
    locs = [k for k in range(1,len(data)-1) if (((data[k]>data[k-1]) and (data[k]>data[k+1])) and (data[k]>threshold))]
    v_locs = verify_locs(locs,fs)
    peaks = [data[i] for i in v_locs]
    #Return
    #   v_locs: <list> arreglo con los índices donde ocurren picos válidos
    #   peaks: <list> arreglo con los picos (señal evaluada en v_locs)
    return v_locs, peaks

def interpolate(x, y, fs):
    #interpola una señal

    #Argumentos:
    #   x: <list> valores en el eje x para interpolar
    #   y: <list> valores en y que corresponden a x
    #   fs: <float> frecuencia de muestreo deseada
    f = interp1d(x, y, kind='cubic')
    x2 = np.arange(start=x[0], stop=x[-1], step=1/fs)
    y2 = f(x2)
    #Return:
    #   x2: <numpy.ndarray> vector en x interpolado
    #   y2: <list> vector en y interpolado 
    return x2, y2    

def verify_locs(locs, fs):
    locs1 = [k for k in locs]
    Mmin = 0.3*fs
    e=1
    while (e<len(locs1)):
        if ((locs1[e]-locs1[e-1])<Mmin):
            d = [locs1.pop(e)]
        e = e+1
    return locs1

def verify_locs2(locs, fs):
    locs1 = [k for k in locs]
    Mmin = 1.5*fs
    e=1
    while (e<len(locs1)):
        if ((locs1[e]-locs1[e-1])<Mmin):
            d = [locs1.pop(e)]
        e = e+1
    return locs1

def bandpower(freq, psd, f1, f2):
    #calcula la potencia en la banda de frecuencia f1-f2 sobre la psd (densidad espectral de potencia)

    #Argumentos
    #   freq: <numpy.ndarray> vector de frecuencias
    #   psd: <numpy.ndarray> densidad espectral de potencia con respecto a freq
    #   f1: <float> inicio de la banda de frecuencia 
    #   f2: <float> fin de la banda de frecuencia
    r = float(trapezoid(psd[np.argmin(abs(f1-freq)):np.argmin(abs(f2-freq))], x=freq[np.argmin(abs(f1-freq)):np.argmin(abs(f2-freq))]))
    #Return:
    #   r: <float> potencia en la banda f1-f2 
    return trapezoid(psd[np.argmin(abs(f1-freq)):np.argmin(abs(f2-freq))], x=freq[np.argmin(abs(f1-freq)):np.argmin(abs(f2-freq))])
    
def min_max_normalize(data):
    #normaliza con mínimo-máximo

    #Argumentos:
    #   data: <numpy.ndarray> datos para normalizar con mínimo-máximo
    n = (data-min(data))/(max(data)-min(data))
    #Return:
    #   n: <numpy.ndarray> datos normalizados
    return n

def get_hrt_eda(scr, fs):
    #calcula el tiempo de recuperación media en un segmento que contiene una respuesta de conductancia de la piel SCR

    #Argumentos:
    #   scr: <numpy.ndarray> segmento con una SCR de EDA
    #   fs: <float> frecuencia de muestreo
    amp = max(scr)-scr[0]
    r = abs(abs(scr[np.argmax(scr):]-max(scr))-amp/2)
    r2 = (abs(scr[np.argmax(scr):]-max(scr))-amp/2)
    z = [k+0.5 for k in range(0,len(r2)-1) if r2[k]<0 and r2[k+1]>=0]
    if len(z)>0:
        hrt = float((1/fs)*z[0])
        print("hrt",hrt)
    else:
        hrt=None
    #Return:
    #   hrt: <float> o None, tiempo de recuperación media 
    return hrt

def rms(data):
    #calcula el valor RMS de una señal

    #Argumentos:
    #   data: <list> o <numpy.ndarray> señal para calcular el RMS
    rms = float(np.sqrt(np.mean(np.array(data)**2)))
    #Return:
    #   rms: <float> valor rms
    return rms

def readtxt(dirfilename):
    #lee un archivo

    #Argumentos:
    #   dirfilename: <str> nombre del archivo
    f = open(dirfilename,'r')
    lin = [k[0:-1] for k in f]
    f.close()
    #Return:
    #   lin: <list> líneas del archivo
    return lin

def det_assessment(real, predicted, fs=100, w=0.05):
    sw = math.floor(fs*w)
    FN = 0
    distance = [[abs(k-i) for i in predicted if (abs(k-i)<=sw)] for k in real]
    VP = sum([1 for k in distance if len(k)>0])
    FP = len(predicted)-VP
    FN = len(real)-VP
    return VP, FP, FN

def plot31(x,y,xlbl,ylbl,xt,txt,coor=[[0,0],[0,0],[0,0]]):
    # crea una gráfica de 3x1
    fig, [ax1,ax2,ax3] = plt.subplots(3,1,sharex=True)
    ax1.plot(x[0],y[0],'#1A535C')
    ax1.set_ylabel(ylbl[0],labelpad=1)
    ax2.plot(x[1],y[1],'#59CAC2')
    ax2.set_ylabel(ylbl[1],labelpad=1)
    ax3.plot(x[2],y[2],'#FFA044')
    ax3.set_ylabel(ylbl[2],labelpad=1)
    ax3.set_xlabel(xlbl,labelpad=1)
    ax1.grid(b=True,axis='both',which='both')
    ax2.grid(b=True,axis='both',which='both')
    ax3.grid(b=True,axis='both',which='both')
    ax3.set_xlim(xt[0],xt[1])
    ax1.text(coor[0][0],coor[0][1],txt[0],fontfamily='monospace')
    ax2.text(coor[1][0],coor[1][1],txt[1],fontfamily='monospace')
    ax3.text(coor[2][0],coor[2][1],txt[2],fontfamily='monospace')
    return [fig,ax1,ax2,ax3]

def plot21(x,y,xlbl,ylbl,xt,txt):
    #crea una gráfica de 2x1
    fig, [ax2,ax3] = plt.subplots(2,1,sharex=True)
    ax2.plot(x[0],y[0],'#59CAC2')
    ax2.set_ylabel(ylbl[0],labelpad=1)
    ax3.plot(x[1],y[1],'#FFA044')
    ax3.set_ylabel(ylbl[1],labelpad=1)
    ax3.set_xlabel(xlbl,labelpad=1)
    ax2.grid(b=True,axis='both',which='both')
    ax3.grid(b=True,axis='both',which='both')
    ax3.set_xlim(xt[0],xt[1])
    ax2.text(180,0,txt[0],fontfamily='monospace')
    ax3.text(200,0.04,txt[1],fontfamily='monospace')
    return [fig,ax2,ax3]

def plots(x,y,xlbl,ylbl,xt,txtbox="",coor=[0.6,30],title="",color='#1A535C'):
    #crea una gráfica
    fig, ax1 = plt.subplots(1)
    ax1.plot(x,y,color)
    ax1.set_xlabel(xlbl,labelpad=1)
    ax1.set_ylabel(ylbl,labelpad=1)
    ax1.set_xlim(xt[0],xt[1])
    ax1.set_title(title)
    ax1.grid(b=True,axis='both',which='both')
    if txtbox!="":
        ax1.text(coor[0],coor[1],txtbox)
    return [fig, ax1]

def plot2(t1,x1,t2,x2,xt,tag1,tag2,txt,coor=[235,0.98]):
    #crea 2 gráficas
    fig, [ax1, ax2] = plt.subplots(2,1,sharex=True)
    ax1.plot(t1,x1,'#1A535C')
    ax1.set_ylabel(tag1)
    ax1.grid(b=True,axis='both',which='both')
    ax2.plot(t2,x2,'#59CAC2')
    ax2.set_ylabel(tag2)
    ax2.set_xlabel("tiempo (s)",labelpad=1)
    ax2.grid(b=True,axis='both',which='both')
    ax2.set_xlim(xt[0],xt[1])
    ax1.text(coor[0],coor[1],txt)
    return [fig, ax1, ax2]   

def plotlocs(t,signal,locs,tag,units,xt):
    #Crea una gráfica con eñaladores en locs
    tlocs = [t[k] for k in locs]
    ylocs = [signal[k] for k in locs]
    fig, ax1 = plt.subplots(1)
    ax1.plot(t,signal,'#1A535C')
    ax1.plot(tlocs,ylocs,'#FFA044', linestyle='', marker=mpl.markers.CARETDOWN)
    ax1.set_xlabel("tiempo (s)",labelpad=0.1)
    ax1.set_ylabel(tag+" ("+units+")")
    ax1.set_xlim(xt[0],xt[1])
    ax1.grid(b=True,axis='both',which='both')
    return [fig, ax1]

def plotIC(t,I,C,tagI,tagC,units,x):
    #grafica
    fig, [ax1,ax2] = plt.subplots(2,1,sharex=True)
    ax1.plot(t,I,'#1A535C')
    ax1.set_ylabel(tagI)
    ax1.grid(b=True,axis='both',which='both')
    ax2.plot(t,C,'#59CAC2')
    ax2.set_ylabel(tagC)
    ax2.set_xlabel("tiempo (s)",labelpad=1)
    ax2.grid(b=True,axis='both',which='both')
    ax2.set_xlim(x[0],x[1])
    return [fig, ax1, ax2]

def plotCF_fft(t,C,F,tag,units,xt,xf,fs,dcr=1):
    #gráfica
    fig, [[ax3,ax4],[ax5,ax6]] = plt.subplots(2,2)
    [freqsC, magC] = fft(C,fs)
    [freqsF, magF] = fft(F,fs)
    magC = magC/dcr
    magF = magF/dcr
    ax3.plot(t,C,'#FFA044')
    ax3.set_ylabel(tag+" ("+units+")")
    ax3.set_xlim(xt[0],xt[1])
    ax4.plot(freqsC,magC,'#FFA044')
    ax4.set_xlim(xf[0],xf[1])
    ax4.set_ylabel("|FFT| (×10^"+str(int(round(1/math.log(10,dcr),0)))+")")
    ax5.plot(t,F,'#1A535C')
    ax5.set_ylabel(tag+" ("+units+")")
    ax5.set_xlim(xt[0],xt[1])
    ax5.set_xlabel("tiempo (s)")
    ax6.plot(freqsF,magF,"#1A535C")
    ax6.set_xlabel("frecuencia (Hz)")
    ax6.set_ylabel("|FFT| (×10^"+str(int(round(1/math.log(10,dcr),0)))+")")
    ax6.set_xlim(xf[0],xf[1])
    ax3.grid(b=True,axis='both',which='both')
    ax4.grid(b=True,axis='both',which='both')
    ax5.grid(b=True,axis='both',which='both')
    ax6.grid(b=True,axis='both',which='both')
    return [fig, ax3, ax4, ax5, ax6]

def plotCF_fft2(t,C,F,tag,units,xt,xf,fs,dcr=[1,1]):
    #gráfica
    fig, [[ax3,ax4],[ax5,ax6]] = plt.subplots(2,2)
    [freqsC, magC] = fft(C,fs)
    [freqsF, magF] = fft(F,fs)
    magC = magC/dcr[0]
    magF = magF/dcr[1]
    ax3.plot(t,C,'#FFA044')
    ax3.set_ylabel(tag+" (RAW)")
    ax3.set_xlim(xt[0],xt[1])
    ax4.plot(freqsC,magC,'#FFA044')
    ax4.set_xlim(xf[0],xf[1])
    ax4.set_ylabel("|FFT| (×10^"+str(int(round(1/math.log(10,dcr[0]),0)))+")")
    ax5.plot(t,F,'#1A535C')
    ax5.set_ylabel(tag+" ("+units+")")
    ax5.set_xlim(xt[0],xt[1])
    ax5.set_xlabel("tiempo (s)")
    ax6.plot(freqsF,magF,"#1A535C")
    ax6.set_xlabel("frecuencia (Hz)")
    ax6.set_ylabel("|FFT| (×10^"+str(int(round(1/math.log(10,dcr[1]),0)))+")")
    ax6.set_xlim(xf[0],xf[1])
    ax3.grid(b=True,axis='both',which='both')
    ax4.grid(b=True,axis='both',which='both')
    ax5.grid(b=True,axis='both',which='both')
    ax6.grid(b=True,axis='both',which='both')
    return [fig, ax3, ax4, ax5, ax6]

def plot22(x,y,xlbl,ylbl,xt):
    #gráfica
    fig, [[ax3,ax4],[ax5,ax6]] = plt.subplots(2,2)
    ax3.plot(x[0],y[0],'#FFA044')
    ax3.set_ylabel(ylbl[0])
    ax3.set_xlim(xt[0],xt[1])
    ax4.plot(x[1],y[1],'#1A535C')
    ax4.set_xlim(xt[0],xt[1])
    ax5.plot(x[2],y[2],'#FFA044')
    ax5.set_ylabel(ylbl[1])
    ax5.set_xlim(xt[0],xt[1])
    ax5.set_xlabel(xlbl[0])
    ax6.plot(x[3],y[3],"#1A535C")
    ax6.set_xlabel(xlbl[1])
    ax6.set_ylabel(ylbl[1])
    ax6.set_xlim(xt[0],xt[1])
    ax3.grid(b=True,axis='both',which='both')
    ax4.grid(b=True,axis='both',which='both')
    ax5.grid(b=True,axis='both',which='both')
    ax6.grid(b=True,axis='both',which='both')
    return [fig, ax3, ax4, ax5, ax6]

def showdialog(size,icon="🛈",msg="",toast=True):
    #muestra cuadro de diálogo para la interfaz
    dlog = QtWidgets.QDialog()
    uidlog = dlg.Ui_Dialog()
    uidlog.setupUi(dlog,icon,msg,toast,size)
    if toast:
        dlog.exec()
    else:
        dlog.exec()

def concat(x,y,sort=False):
    #concatena dos arreglos

    #Argumentos:
    #   x: <list> arreglo 1
    #   y: <list> arreglo 2
    #   sort: <bool> si es True se ordena el arreglo resultante
    z = []
    [z.append(k) for k in x]
    [z.append(k) for k in y]
    if sort:
        z.sort()
    #Return:
    #   z: <list> concatenación de x e y
    return z
    
def print_dict(d, key=False):
    if key!=False:
        print(key,':{')
    [print(k,':',type(d[k])) if type(d[k])!=type({}) else print_dict(d[k], key=k) for k in d]
    if key!=False:
        print('}')

def reshape(m):
    #crea un arreglo plano de uno multidimensional

    #Argumentos:
    #   m: <list> arreglo multidimensional
    flatten_list = lambda y:[x for a in y for x in flatten_list(a)] if type(y) is list else [y]
    flat = flatten_list(m)
    #Return;
    #   flat: <list> areglo unidimensional con todos los elementos del arreglo inicial
    return flat

