#-------------Objeto para procesar electromiografía (EMG)

from functions import *
import biosppy.signals as bio

class emg():
    def __init__(self, v, n, fs, convert=False):
        # se crea un objeto para procesar respiración

        #Argumentos:
        #   v: <numpy.ndarray> o <list>  señal de voltage de EMG
        #   n: <int> resolución en bits
        #   fs: <float> frecuencia de muestreo
        #   convert: <bool> si es True se aplca conversión a mV
        self.fs = fs
        self.time = (1/fs)*np.arange(len(v))
        if convert==True:
            self.emg_mV = (((((v/(2**n))-0.5)*3.3)/(1009))*1000)
        else:
            self.emg_mV = v
        
        #self.filtered_emg = highpass_filter(self.emg_mV,20,fs,3)
        self.filtered_emg = bandpass_filter(self.emg_mV,20,400,fs,3)
        #rectified = abs(self.filtered_emg)
        #self.emg_envelope = envelope(rectified)

    def get_mV(self):
        #Return:
        #   self.time: <numpy.ndarray> vector de tiempo correspondiente a la señal
        #   self.emg_mV: <numpy.ndarray> EMG convertido a mV sin filtrar
        #   self.filtered_emg: <numpy.ndarray> EMG filtrado
        return self.time, self.emg_mV, self.filtered_emg

    def get_meanrms(self,w,mvc=True):
        #calcula el rms medio

        #Argumentos:
        #   w: <float> ancho de la ventana  en segundos para segmentar y calccular RMS
        #   mvc: <bool> o <float> si es booleano no se obtiene el porcentaje relativo a la máxima contracción voluntaria
        #                           si es float, se halla el porcentaje relativo a mvc que seria el RMS en la máxima contracción voluntaria
        s = int(round(w*self.fs,0))
        segments = [self.emg_mV[i:i+s] for i in range(0,len(self.emg_mV)-1,s)]
        times = w*np.arange(len(segments))+w/2
        vrms = [np.sqrt(np.mean(np.array(k)**2)) for k in segments]  
        if type(mvc)==type(True):
            meanrms = float(np.mean(vrms))
            percrms = meanrms
        else:
            meanrms = float(np.mean(vrms))
            percrms = meanrms/mvc
        #Return:
        #   meanrms: <float> valor rms promedio en todos los segmentos de tamaño de la ventana w
        #   vrms:   <list> valor rms de cada segmento
        #   percrms: <float> porcentaje con respecto a la máxima contracción voluntaria
        #   times: <numpy.ndarray> vectores de tiempo de los segmentos
        return meanrms, vrms, percrms, times

        


