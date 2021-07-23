#-------------Objeto para procesar fotopletismografía (PPG) y espectroscopia de infrarrojo cercano (fNIRS)

from functions import *
import biosppy.signals as bio
class ppg():
    def __init__(self,r,ir,n,fs,convert=False,fnirs=True):
        # se crea un objeto para procesar PPG o fNIRS

        #Argumentos:
        #   r: <numpy.ndarray> o <list>  señal de corriente en el rojo
        #   ir: <numpy.ndarray> o <list>  señal de corriente en el infrarrojo
        #   n: <int> resolución en bits
        #   fs: <float> frecuencia de muestreo
        #   convert: <bool> si es True se aplca conversión a % de expansipon torácica
        #   fnirs: <bool> si es True se procesa como fNIRS, si no se procesa como PPG
        self.fs = fs
        if convert==True:
            self.rc_uA = lowpass_filter(((0.15*r)/(2**n)),30,fs,3)
            self.irc_uA = lowpass_filter(((0.15*ir)/(2**n)),30,fs,3)
        else:
            self.rc_uA = lowpass_filter(r,15,fs,3)
            self.irc_uA = lowpass_filter(ir,15,fs,3)               
        self.time = (1/self.fs)*np.arange(len(self.rc_uA))
        self.fnirs = fnirs
    
    def get_filtered(self):
        #Return:
        #   self.rc_uA : <numpy.ndarray> señal de corriente en el rojo filtrada
        #   seld.irc_uA: <numpy.ndarray> señal de corriente en el infrarrojo filtrada
        return [self.rc_uA, self.irc_uA]

    def get_spo2(self,sfs=14):
        #calcula saturación de O2

        #Argumentos:
        #   sfs: <int> nueva frecuencia de muestreo para la señal de saturación interpolada.
        times, segmentsR, segmentsI, self.locs = segment_ppg(self.rc_uA, self.irc_uA, self.fs, self.time, fnirs=self.fnirs)
        self.filtered = bandpass_filter(-self.irc_uA,0.5,5,self.fs,3)
        VppR = np.array([get_vpp(np.array(k)) for k in segmentsR])
        VppI = np.array([get_vpp(np.array(k)) for k in segmentsI])
        VavgR = np.array([get_avg(np.array(k)) for k in segmentsR])
        VavgI = np.array([get_avg(np.array(k)) for k in segmentsI])
        R = (VppR*VavgI)/(VavgR*VppI)
        Sat = 110-(25*R)
        tsat = np.array([np.mean(times[k]) for k in range(0,len(times))])
        intt, intS = interpolate(tsat, Sat, sfs)
        filS = lowpass_filter(intS,0.2,sfs,2)
        if max(Sat)>100:
            maxS = 100
        else:
            maxS = float(round(max(Sat),1))
        minS = float(round(min(Sat),1))
        meanS = float(round(np.mean(Sat),1))

        #Return:
        #   fils: <numpy.ndarray> señal de saturación interpolada a la frecuencia sfs y filtrada
        #   intt: <numpy.ndarray> vector de tiempo para la saturación interpolada
        #   Sat: <numpy.ndarray> vector con la serie de tiempo de la saturación
        #   tsat: <numpy.ndarray> vector de tiempo de la serie de tiempo de la saturación
        #   minS: <float> saturación mínima
        #   meanS: <float> saturación media
        #   maxS: <maxS> saturación máxima
        return filS, intt, Sat, tsat, minS, meanS, maxS
    
    def get_HRV(self):
        #calcula variabilidad de la frecuencia cardiaca desde los picos sistólicos
        if self.fnirs:
            filtered = bandpass_filter(-self.irc_uA,0.5,5,self.fs,3)
        else:
            filtered = bandpass_filter(self.irc_uA,0.5,5,self.fs,3)
        times, segmentsR, segmentsI, self.locs = segment_ppg(self.rc_uA, self.irc_uA, self.fs, self.time, fnirs=self.fnirs)
        rrintervals = (1/self.fs)*derivative(np.array(self.locs))
        rrintervals_times = (1/self.fs)*adj_mean(np.array(self.locs))
        inter_rrt, inter_rri = interpolate(rrintervals_times,rrintervals,14)
        detr_rri = detrend(inter_rri,type='linear',overwrite_data=False)
        fft_intervals = fft(detr_rri,14)
        [Xfft, Yfft] = [fft_intervals[0],fft_intervals[1]]
        VLF = bandpower(Xfft,Yfft,0,0.04)
        LF = bandpower(Xfft,Yfft,0.04,0.15)
        HF = bandpower(Xfft,Yfft,0.15,0.4)
        LF_HF = LF/HF
        freqbpm = (1/rrintervals)*(60)
        fMax = max(freqbpm)
        fMean = np.mean(freqbpm)
        fMin = min(freqbpm)
        #Return:
        #   filtered: <numpy.ndarray> señal filtrada
        #   self.locs: <list> índices donde ocurren los picos sistólicos
        #   rrintervals: <numpy.ndarray> serie de tiempo de los periodos cardiacos
        #   rrintervals_times: <numpy.ndarray> vector de tiempo de la serie de tiempo de los periodos cardiacos
        #   inter_rrt: <numpy.ndarray> vector de tiempo interpolado de los periodos cardiacos
        #   detr_rri: <numpy.ndarray> periodos cardiacos interpolados y sin tendencia lineal
        #   Xfft: <numpy.ndarray> eje x de la densidad espectral de potencia de los periodos cardiacos
        #   Yfft: <numpy.ndarray> eje y de la densidad espectral de potencia de los periodos cardiacos
        #   VLF: <float> potencia en la banda de muy baja frecuencia (0-0.04Hz) 
        #   LF: <float> potencia en la banda de baja frecuencia (0.04-0.15Hz)
        #   HF: <float> potencia en la banda de alta frecuencia (0.15-0.4Hz)
        #   LF_HF: <float> división LF/HF
        #   fMax: <float> frecuencia cardiaca máxima
        #   fMean: <float> frecuencia cardiaca media
        #   fMin: <float> frecuencia cardiaca mínima
        return filtered, self.locs, rrintervals, rrintervals_times, inter_rrt, detr_rri, Xfft, Yfft, VLF, LF, HF, LF_HF, fMax, fMean, fMin
