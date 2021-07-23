#-------------Objeto para procesar actividad electrodérmica (EDA)


from functions import *
import biosppy.signals as bio
from pyEDA.main import *
class eda():
    def __init__(self, c, n, fs, convert=False):
        # se crea un objeto para procesar actividad electrodérmica

        #Argumentos:
        #   c: <numpy.ndarray> o <list>  señal de conductancia de la piel
        #   n: <int> resolución en bits
        #   fs: <float> frecuencia de muestreo
        #   convert: <bool> si es True se aplca conversión a uS
        self.fs = fs
        self.time = (1/fs)*np.arange(len(c))
        if convert==True:
            self.eda_uS = ((c/(2**n))*3.3)/0.132
        else:
            self.eda_uS = c
        self.prefiltered_eda = lowpass_filter(self.eda_uS, 0.37, fs, 3)
        m, self.wd, filtered2 = process_statistical(self.prefiltered_eda, use_scipy=True, sample_rate=self.fs, new_sample_rate=40, segment_width=600, segment_overlap=0)
        #ver si se puede procesar sin normalizar 
        self.filtered_eda = filtered2[0][int(40):-1]
        self.time40Hz = (1/40)*np.arange(len(self.filtered_eda))
        self.filtered_phasic_eda = self.wd['filtered_phasic_gsr'][0] [int(self.fs):-1]
        self.phasic_eda = self.wd['phasic_gsr'][0][int(40):-1]
        self.tonic_eda = self.wd['tonic_gsr'][0][int(40):-1]
        locs = [i for i in range(1,len(self.phasic_eda)-1) if self.phasic_eda[i]>self.phasic_eda[i-1] and self.phasic_eda[i]>self.phasic_eda[i+1]]
        locs_min = [i for i in range(1,len(self.phasic_eda)-1) if self.phasic_eda[i]<self.phasic_eda[i-1] and self.phasic_eda[i]<self.phasic_eda[i+1]]
        self.nscr = 0
        if len(locs) == len(locs_min):
            if locs[0]<locs_min[0]:
                locs_min.insert(0,0)
                amps = [self.phasic_eda[locs[i]]-self.phasic_eda[locs_min[i]] for i in range(0,len(locs))]
                self.nscr = len([locs[i] for i in range(0,len(amps)) if amps[i]>0.01])
                locs_min = locs_min[1:]
                amps = [self.phasic_eda[locs[i+1]]-self.phasic_eda[locs_min[i]] for i in range(0,len(locs)-1)]
                locs = [locs[i+1] for i in range(0,len(amps)-1) if amps[i]>0.01]    
                locs_min = [locs_min[i] for i in range(0,len(amps)) if amps[i]>0.01]
            else:
                amps = [self.phasic_eda[locs[i]]-self.phasic_eda[locs_min[i]] for i in range(0,len(locs))]
                locs = [locs[i] for i in range(0,len(amps)) if amps[i]>0.01]
                self.nscr = len(locs)
                locs_min = [locs_min[i] for i in range(0,len(amps)) if amps[i]>0.01]
        elif len(locs)>len(locs_min):
            locs_min.insert(0,0)
            amps = [self.phasic_eda[locs[i]]-self.phasic_eda[locs_min[i]] for i in range(0,len(locs_min))]
            self.nscr = len([locs[i] for i in range(0,len(amps)) if amps[i]>0.01])
            locs_min = locs_min[1:]
            amps = [self.phasic_eda[locs[i+1]]-self.phasic_eda[locs_min[i]] for i in range(0,len(locs)-1)]
            locs = [locs[i+1] for i in range(0,len(amps)-1) if amps[i]>0.01]
            locs_min = [locs_min[i] for i in range(0,len(amps)) if amps[i]>0.01]
        elif len(locs)<len(locs_min):
            amps = [self.phasic_eda[locs[i]]-self.phasic_eda[locs_min[i]] for i in range(0,len(locs))]
            locs = [locs[i] for i in range(0,len(amps)) if amps[i]>0.01]
            self.nscr = len(locs)
            locs_min = [locs_min[i] for i in range(0,len(amps)) if amps[i]>0.01]
        else:
            pass
        print('locs_min',locs_min)
        print('locs',locs)
        self.scr_freq = (self.nscr/((1/40)*len(self.phasic_eda)))*(60)
        if len(locs)==0:
            self.nscr = 0
        if self.nscr>0:
            if locs[-1]>locs_min[-1]:
                locs_min.append(len(self.phasic_eda)-1)
            segments = [self.phasic_eda[locs_min[i]:locs_min[i+1]] for i in range(0,len(locs_min)-1)]
            segments_times = [self.time40Hz[locs_min[i]:locs_min[i+1]] for i in range(0,len(locs_min)-1)]
            self.aux4 = segments_times
            self.aux3 = segments
            tonic_segments = [self.tonic_eda[locs_min[i]:locs_min[i+1]] for i in range(0,len(locs_min)-1)]
            self.aux5 = tonic_segments
            scls = [np.mean(k) for k in tonic_segments]
            rtimes = (1/40)*np.array([np.argmax(k) for k in segments])
            amplitudes = [max(k)-k[0] for k in segments]
            hlfrt = [get_hrt_eda(k,40) for k in segments]
            self.aux1 = locs
            self.aux2 = locs_min
            self.aux4 = segments_times
            self.aux3 = segments
            self.aux6 = scls
            self.aux8 = amplitudes
            self.aux7 = rtimes
            self.aux9 = hlfrt
            self.ext_scrs = [scls,amplitudes,rtimes,hlfrt]           
         
    def get_filtered(self):
        #Return:
        #   self.filtered_eda: <numpy.ndarray> actividad electrodérmica filtrada
        return self.filtered_eda
    
    def get_filtered_fsi(self):
        t, y = interpolate(self.time40Hz,self.filtered_eda,self.fs)
        return y
    
    def get_converted(self):
        #Return:
        #   self.filtered_eda: <numpy.ndarray> actividad electrodérmica en uS sin filtrar   
        return self.eda_uS
    
    def get_meanresp(self):
        # obitene las características promedio entre todas las respuestas SCR
        if self.nscr==0:
            scl, hlfrt, amp, riset = None, None, None, None
        elif self.nscr==1:
            amp = self.ext_scrs[1][0]
            riset = self.ext_scrs[2][0]
            scl = self.ext_scrs[0][0]
            hlfrt = self.ext_scrs[3][0]
            #para los 4 valores toca tomar el de la posición 0 porque están así: [valor] y se necesita así: valor
        else:
            try:
                hlfrt = float(np.mean([k for k in self.ext_scrs[3] if k!=None]))
            except:
                hlfrt = None
            amp = float(np.mean(self.ext_scrs[1]))
            riset = float(np.mean(self.ext_scrs[2]))
            scl = float(np.mean(self.ext_scrs[0]))
        #Return
        #   scl: <float> o  None, conductancia basal promedio o None en caso que no se encuentren SCRs
        #   amp: <float> o None, amplitud promedio o None en caso que no se encuentren SCRs
        #   riset: <float> o None, tiempo de elevación promedio o None en caso que no se encuentren SCRs
        #   scl: <float> o None, tiempo de recuperación media promedio, None si la señal no alcanza para calcularlo
        return scl, amp, riset, hlfrt

    def get_maxresp(self):
        # obtiene las características de la SCR con mayor amplitud
        if self.nscr==0:
           scl, hlfrt, amp, riset = None, None, None, None
        elif self.nscr==1:
            amp = self.ext_scrs[1][0]
            riset = self.ext_scrs[2][0]
            scl = self.ext_scrs[0][0]
            hlfrt = self.ext_scrs[3][0]
        else:
            scls, amps, risets, hlfrts = self.ext_scrs
            i = np.argmax(amps)
            scl, amp, riset, hlfrt = float(scls[i]), float(amps[i]), float(risets[i]), float(hlfrts[i])
        #Return:
        #   scl: <float> o  None, conductancia basal  o None en caso que no se encuentren SCRs
        #   amp: <float> o None, amplitud o None en caso que no se encuentren SCRs
        #   riset: <float> o None, tiempo de elevación o None en caso que no se encuentren SCRs
        #   scl: <float> o None, tiempo de recuperación media, None si la señal no alcanza para calcularlo
        return scl, amp, riset, hlfrt

    def get_parameters(self):
        #Return:
        #   self.time: <numpy.ndarray> vector de tiempo de la señal
        #   self.time40Hz: <numpy.ndarray> vector de tiempo de la señal interpolada a 40Hz
        #   self.filtered_eda: <numpy.eda> actividad electrodérmica filtrada
        #   self.wd: <dict> diccionario con los datos que procesa la librería pyEDA
        #   self.filtered_phasic_eda: <numpy.ndarray> actividad fásica extraída de EDA filtrada y con frecuencia de muestreo 40Hz
        #   self.tonic_eda: <numpy.ndarray> actividad tónica extraída de EDA y con frecuencia de muestreo de 40 Hz
        #   self.ext_scrs: <list> [scls,amplitudes,rtimes,hlfrt]
        return self.time, self.time40Hz, self.filtered_eda, self.wd, self.filtered_phasic_eda, self.tonic_eda, self.ext_scrs
