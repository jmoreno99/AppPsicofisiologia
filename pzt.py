#-------------Objeto para procesar respiración (RESP)
from functions import *
import biosppy.signals as bio

class pzt():
    def __init__(self, p, n, fs, convert=False):
        # se crea un objeto para procesar respiración

        #Argumentos:
        #   t: <numpy.ndarray> o <list>  señal de respiración
        #   n: <int> resolución en bits
        #   fs: <float> frecuencia de muestreo
        #   convert: <bool> si es True se aplca conversión a % de expansipon torácica
        try:
            p = np.array([k[0] for k in p])
        except:
            pass
        self.fs = fs
        self.time = (1/fs)*np.arange(len(p))
        if convert==True:
            self.pzt_p = ((p/(2**n))-0.5)*100
        else:
            self.pzt_p = p
        self.filtered_pzt = bandpass_filter(self.pzt_p,0.1,0.5,self.fs,3)
        self.locs, d = find_peaks(self.filtered_pzt)

   
    def get_filtered(self):
        #Return:
        # self.filtered_pzt: <numpy.ndarray> señal de respiración filtrada
        return self.filtered_pzt
    
    def get_converted(self):
        #Return:
        # self.filtered_pzt: <numpy.ndarray> señal de respiración convertida a % sin filtrar
        return self.pzt_p

    def getfreq(self):
        #Extrae frecuencia respiratoria

        rrintervals = (1/self.fs)*derivative(np.array(self.locs))
        rrintervals_times = (1/self.fs)*adj_mean(np.array(self.locs))
        freq_rpm = 60/rrintervals
        maxR  = float(max(freq_rpm))
        minR = float(min(freq_rpm))
        meanR = float(np.mean(freq_rpm))
        #Return:
        #   maxR: <float> frecuencia respiratoria máxima
        #   minR: <float> frecuencia respiratoria máxima
        #   meanR: <float> frecuencia respiratoria máxima
        #   freq_rpm: <numpy.ndarray> serie de tiempo de la frecuencia respiratoria
        #   rrintervals_times: <numpy.ndarray> tiempos de la serie de tiempo de la frecuencia respiratoria
        return maxR, minR, meanR, freq_rpm, rrintervals_times

    def getexp(self):
        #esto función no se usa
        maxP = max(self.filtered_pzt)
        minP = min(self.filtered_pzt)
        meanP = np.mean(self.filtered_pzt)
        return maxP, minP, meanP

    def get_rpeaks(self):
        #Return:
        #   self.locs: <list> arreglo con los índices donde ocurren las respiraciones
        return self.locs