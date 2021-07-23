#-------------Objeto para procesar temperatura (TEMP)
from functions import *
import biosppy.signals as bio

class temp():
    def __init__(self, t, n, fs,convert=False):
        #se crea un objeto para procesar temperatura
        
        #Argumentos:
        #   t: <numpy.ndarray> o <list>  señal de temperatura
        #   n: <int> resolución en bits
        #   fs: <float> frecuencia de muestreo
        #   convert: <bool> si es True se aplca conversión a °C
        self.fs = fs
        self.time = (1/fs)*np.arange(len(t))
        if convert==True:
            ntcv = np.array((t*3.3)/(2**n))
            ntcr = (10000*ntcv)/(3.3-ntcv)
            a0 = 0.00112764514
            a1 = 0.000234282709
            a2 = 0.0000000877303013
            temp_kelvin =  1/(a0+(a1*np.log(ntcr))+(a2*((np.log(ntcr))**3)))
            temp_celcius = temp_kelvin-273.15
            temp_farenheit = (temp_celcius*9/5)+32 
        else:
            temp_celcius = t
            temp_kelvin = temp_celcius+273.15
        self.temp_K = temp_kelvin
        self.temp_C = temp_celcius
        self.temp_F = temp_farenheit
        self.filtered_K = lowpass_filter(self.temp_K,0.5,self.fs,3)
        self.filtered_C = lowpass_filter(self.temp_C,0.5,self.fs,3)
        self.filtered_F = lowpass_filter(self.temp_F,0.5,self.fs,3)

    def get_C(self):
        # extrae la temperatura máxima, mínima y media en °C
        maxT = float(max(self.filtered_C))
        minT = float(min(self.filtered_C))
        meanT = float(np.mean(self.filtered_C))
        #Return:
        #   maxT: <float> temperatura máxima en °C
        #   minT: <float> temperatura mínima en °C
        #   medT: <float> temperatura media en °C
        return maxT, minT, meanT
    
    def get_F(self):
        # extrae la temperatura máxima, mínima y media en °F
        maxT = max(self.filtered_F)
        minT = min(self.filtered_F)
        meanT = np.mean(self.filtered_F)
        #Return:
        #   maxT: <float> temperatura máxima en °F
        #   minT: <float> temperatura mínima en °F
        #   medT: <float> temperatura media en °F
        return maxT, minT, meanT



