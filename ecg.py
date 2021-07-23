#--------------------------------objeto para procesar electrocardiografía (ECG)
from functions import *
import biosppy.signals as bio


class ecg():
    def __init__(self,v,n,fs,convert=False):
        #se crea el objeto para procesar ECG

        #Argumentos
        #   v: <numpy.ndarray> señal de voltaje ECG
        #   n: <int> resolución en bits
        #   fs: <float> frecuencia de muestreo
        #   convert: <bool> si es True se aplca conversión a mV
        try:
            v = np.array([k[0] for k in v])
        except:
            pass
        self.fs = fs
        self.time = (1/fs)*np.arange(len(v))
        self.ecg_mV = (((((v/(2**n))-0.5)*3.3)/(1100))*1000)
        if convert==True:
            self.filtered_ecg = bandpass_filter((((((v/(2**n))-0.5)*3.3)/(1100))*1000), 2.5, 40, fs, 3)
        else:
            self.filtered_ecg = bandpass_filter(v, 2.5, 40, fs, 3)
        ts, self.filtered_ecg_QRS, rpeaks, templates_ts, templates, heart_rate_ts, heart_rate = bio.ecg.ecg(signal=self.filtered_ecg,sampling_rate=fs,show=False)

    def get_mV(self):
        #Return: 
        #   self.time: <numpy.ndarray> vector de tiempo de la señal
        #   self.ecg_mV: <numpy.ndarray> vector con el ECG convertido a mV (sin filtrar)
        #   self.filtered_ecg: <numpy.ndarray> vector con el ECG filtrado
        return self.time, self.ecg_mV, self.filtered_ecg

    def get_rpeaks(self):
        #Return:
        #   rpeaks_locs: <list> vector con los índices donde occuren picos R en el ECG
        ts, filtered, rpeaks, templates_ts, templates, heart_rate_ts, heart_rate = bio.ecg.ecg(signal=self.filtered_ecg,sampling_rate=self.fs,show=False)
        rpeaks_locs = verify_locs(rpeaks,self.fs)
        return rpeaks_locs

    def get_FC(self):
        VLF, LF, HF, LF_HF, inter_rrt, detr_rri, rrintervals_times, rrintervals, arr, rpeaks_peak, Xfft, Yfft = self.get_HRV()
        freqbpm = (1/rrintervals)*(60)
        M = float(max(freqbpm))
        m = float(min(freqbpm))
        m_ = float(np.mean(freqbpm))
        #Return: 
        #   M: <float> Frecuencia cardiaca máxima
        #   m: <float> Frecuencia cardiaca mínima
        #   m_: <float> Frecuencia cardiaca media
        return M, m, m_


    def get_HRV(self):
        #calcula HRV
        rpeaks_locs = self.get_rpeaks()
        rrintervals = (1/self.fs)*derivative(np.array(rpeaks_locs))
        rrintervals_times = (1/self.fs)*adj_mean(np.array(rpeaks_locs))
        inter_rrt, inter_rri = interpolate(rrintervals_times,rrintervals,14)
        detr_rri = detrend(inter_rri,type='linear',overwrite_data=False)
        rpeaks_peak = [self.filtered_ecg_QRS[i] for i in rpeaks_locs]
        fft_intervals = fft(detr_rri,14)
        Xfft, Yfft = interpolate(fft_intervals[0],fft_intervals[1],3000)
        VLF = bandpower(Xfft,Yfft,0,0.04)
        LF = bandpower(Xfft,Yfft,0.04,0.15)
        HF = bandpower(Xfft,Yfft,0.15,0.4)
        LF_HF = LF/HF
        a = np.array([self.time[i] for i in rpeaks_locs])
        #Return: 
        #   VLF: <float> potencia en muy baja frecuencia de HRV (0-0.04 Hz), 
        #   LF: <float> potencia en baja frecuencia de HRV (0.04-0.15Hz), 
        #   HF: <float> potencia en alta frecuencia de HRV (0.15-0.4Hz), 
        #   LF_HF: <float> índice LF/HF, 
        #   inter_rrt: <numpy.ndarray> vector de tiempo de los periodos cardiacos interpolado, 
        #   detr_rri: <numpy.ndarray> periodos cardiacos interpolados y sin tendencia lineal,
        #   rrintervals_times: <numpy.ndarray>  vector de tiempo de la serie de tiempo de los periodos cardiacos
        #   rrintervals: <numpy.ndarray> serie de tiempo de los periodos  cardiacos
        #   a: <numpy.ndarray> tiempos en que ocurren los picos R
        #   rpeaks_peak: <list> picos evaluados en los índices donde se localizaron
        #   fft_intervals[0]: <numpy.ndarray> vector de frecuencia de la densidad espectral de potencia de la HRV
        #   fft_intervals[1]: <numpy.ndarray> vector de magnitud de la densidad espectral de potencia de la HRV
        return VLF, LF, HF, LF_HF, inter_rrt, detr_rri, rrintervals_times, rrintervals, a, rpeaks_peak, fft_intervals[0], fft_intervals[1]