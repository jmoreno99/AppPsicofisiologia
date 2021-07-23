#Archivo para abrir, cerrar archivos y manipular datos

import json

import h5py
from functions import *
from ecg import *
from emg import *
from eda import *
from temp import *
from pzt import *
from ppg import *
import os

def readdata(filename=''):
    # lee archivo .json

    # Argumentos:
    #     filename: <str> nombre de archivo .json para leer
    try:
        json_file = open(filename,'r')
        data = json.load(json_file)
        json_file.close()
    except FileNotFoundError:
        print('archivo'+filename+' no encontrado')
    # Return:
    #     data: <dict> diccionario con los datos del .json
    return data

def writedata(data,filename='data.json'):
    # escribe archivo .json

    #Argumentos:
    #   data: <dict> dccionario para guardar en .json
    #   filename: <str> nombre de archivo .json para guardar el diccionario
    with open(filename,'w+') as json_file:
        json.dump(data,json_file)
        json_file.close()
    
def adduser(id,nombre="",f_nacimiento="",sexo="",peso_kg="",estatura_m="",escolaridad="",antecedentes="",semiologia="",historia_clinica=""):
    # crea un .json y una entrada en tabla.csv, para un paciente nuevo
    
    #Argumentos:
        # id: <str> identificador del paciente
        # f_nacimiento: <str> fecha de necimiento del paciente
        # sexo: <str> sexo del paciente
        # peso_kg: <str> peso del paciente
        # estatura_m: <str> estatura del paciente
        # escolaridad: <str> escolaridad del paciente
        # antecedentes: <str> antecedentes del paciente
        # semiologia: <str> semiología
        # historia_clinica: <str> historia clínica del paciente
    x = len([k.split('.')[0] for k in os.listdir('data')])
    p = {'nombre':nombre,'f_nacimiento':f_nacimiento,'sexo':sexo,'peso_kg':peso_kg,'estatura_m':estatura_m,'escolaridad':escolaridad,'antecedentes':antecedentes,'semiologia':semiologia,'historia_clinica':historia_clinica,'eventos':[]}
    writedata(p,'data/'+id+'.json')
    f = open('table.csv','a')
    f.write(str(id)+';'+str(nombre)+';'+str(f_nacimiento)+';'+str(sexo)+';'+str(antecedentes)+';'+str(peso_kg)+';'+str(estatura_m)+';'+str(escolaridad)+'\n')
    f.close()
    if len([k.split('.')[0] for k in os.listdir('data')])>x:
        # Return:
        #   True si se guradó correctamente
        #   False si no
        return True
    else:
        return False
        
def open_table():
    #abre el archivo tabla.csv donde están los datos personales de los pacientes
 
    f = open('table.csv','r')
    l = [k.replace('\n','').split(';') for k in f]
    print('l inicial',l)
    if len(l)<2:
        l=[]
    else:
        l = [[k[j] if j!=2 else f_toage(k[j]) for j in range(0,len(k))] for k in l[1:]]
        print('lfinal',l)
    f.close()
    
    #Return:
    #    l: <list> arreglo 2D con los datos de table.csv
    return l

def remove_pat(idn):
    #elimina el .json y la entrada en tabla.csv, del paciente 

    #Argumentos:
    #   idn: <str> ID del paciente para eliminar
    os.remove('data/'+str(idn)+'.json')
    f = open('table.csv','r')
    r = [k for k in f]
    l = [k.replace('\n','').split(';') for k in r]
    f.close()
    l = [r[j] if l[j][0]!=idn else [] for j in range(0,len(l))]
    print(l)
    f = open('table.csv','w')
    [f.write(k) for k in l if k!=[]]
    f.close()

def f_toage(f):
    # convierte fecha de nacimiento a edad 

    #Argumentos:
    #   f: <str> fecha en formato "DD/MM/AA"
    r = str(calculateAge(date(int(f.split('/')[2]),int(f.split('/')[1]),int(f.split('/')[0]))))
    #Return
    #   r: <str> edad en cadena
    return r

def get_patdat(file):
    # lee el .json del paciente y devuelve el diccionario

    #Argumentos:
    #   file:<str> nombre del archivo .json del paciente
    d = readdata(filename='data/'+file)
    #Return
    #   d: <dict> diccionario con la información del .json del paciente
    return d

def save_pat(id,pat):
    #guarda el paciente en el respectivo .json

    #Argumentos:
    #   id: <str> ID del paciente para guardar
    #   pat: <dict> diccionario con la información del paciente para guardar
    writedata(pat,'data/'+id+'.json')

""" def append_patdat(idn,jsobj,multiple=False):
    print('entra a append patdat')
    p = readdata(filename='data/'+idn+'.json')
    if multiple:
        [p['eventos'].append(k) for k in jsobj]
    else:
        p['eventos'].append(jsobj)
    writedata(p,'data/'+idn+'.json') """

def edit_patdat(idn,tag,val):
    #esta función no se usa
    p = readdata(filename='data/'+idn+'.json')
    p[tag]=val
    writedata(p,'data/'+idn+'.json')

def search_dctlst(lst,key,val):
    return [k for k in lst if k[key]==val]


def previewsignals(cdir):
    #crea una lista con diccionarios (uno por señal), de las señales guardadas en un archivo txt
    # cada diccionario de la forma: {'med':<medida>, 'filt':{'fs':<frecuencia de muestreo>, 'signal':<señal filtrada>, 'unid':<unidades>}}

    #Argumentos:
    #   cdir: <str> dirección de archivo de adquisición con señales 
    typef = cdir.split(".")[1]
    if typef=="h5":
        f = h5py.File(cdir,'r')
        dev = str(f.keys())[16:-3]
        #cambiar!!!!
        raise ValueError
    elif typef=="txt":
        f = open(cdir,'r')
        d = [k for k in f]
        dt1 = [k.replace('\n','').split('\t') for k in d if not('#' in k)]
        headdct = json.loads(d[1].replace('#','').replace('\n',''))
        dev = str(headdct.keys()).split("'")[1]
        fs = headdct[dev]['sampling rate']
        sig_names =  convertsignames(headdct[dev]['sensor'])
        datelist = list(reversed(headdct[dev]['date'].split('-')))
        dateacq =  datelist[0]+"/"+datelist[1]+"/"+datelist[2]
        ns =  headdct[dev]['resolution']
        n2 = ns[len(ns)-len(sig_names):len(ns)]
        dt2 = [k[len(ns)-len(sig_names):len(ns)] for k in dt1]
        dt2 = np.transpose(np.array([[int(i) for i in j] for j in dt2]))
        f.close()
        obj = []
        i = 0
        q = True
        for name in sig_names:
            if 'PPG' in name:
                ppgobj = ppg(dt2[i],dt2[i],n2[i],fs,convert=True,fnirs=False)
                fIR = ppgobj.irc_uA
                jsobj = {'med':'PPG', 'filt':{'fs':fs, 'signal':fIR, 'unid':'(µA)'}}
                obj.append(jsobj)
            elif 'EDA' in name:
                c = dt2[i]
                nr = n2[i]
                C = ((c/(2**nr))*3.3)/0.132
                F = lowpass_filter(C, 0.37, fs, 3)
                jsobj = {'med':'EDA', 'filt':{'fs':fs, 'signal':F, 'unid':'(µS)'}}
                obj.append(jsobj)
            elif 'RESP' in name:
                respobj = pzt(dt2[i],n2[i],fs,convert=True)
                F = respobj.filtered_pzt
                jsobj = {'med':'RESP', 'filt':{'fs':fs, 'signal':F, 'unid':'(%)'}}
                obj.append(jsobj)
            elif 'EMG' in name:
                emgobj = emg(dt2[i],n2[i],fs,convert=True)
                F = emgobj.filtered_emg
                jsobj = {'med':'EMG', 'filt':{'fs':fs, 'signal':F, 'unid':'(mV)'}}
                obj.append(jsobj)
            elif 'TEMP' in name:
                tempobj = temp(dt2[i], n2[i], fs, convert=True)
                F = tempobj.filtered_C
                jsobj = {'med':'TEMP', 'filt':{'fs':fs, 'signal':F, 'unid':'(°C)'}}
                obj.append(jsobj)
            elif 'ECG' in name:
                ecgobj = ecg(dt2[i], n2[i],fs,convert=True)
                F = ecgobj.filtered_ecg
                jsobj = {'med':'ECG', 'filt':{'fs':fs, 'signal':F, 'unid':'(mV)'}}
                obj.append(jsobj)
            elif 'fNIRS' in name:
                if q:
                    channels = headdct[dev]['channels']
                    indR = channels.index(9)
                    indIR = channels.index(10)
                    R = dt2[indR]
                    IR = dt2[indIR]
                    ppgobj = ppg(R,IR,n2[i],fs,convert=True, fnirs=True)
                    fR = ppgobj.rc_uA
                    fIR =  ppgobj.irc_uA
                    filS, intt, Sat, tsat, minS, meanS, maxS = ppgobj.get_spo2()
                    jsobj = {'med':'fNIRS-rojo', 'filt':{'fs':fs, 'signal':fR, 'unid':'(µA)'}}
                    obj.append(jsobj)
                    jsobj = {'med':'fNIRS-infrarrojo', 'filt':{'fs':fs, 'signal':fIR, 'unid':'(µA)'}}
                    obj.append(jsobj)
                    #jsobj = {'med':'SpO2', 'filt':{'t':tsat,'signal':Sat,'unid':'(%)'}}
                    
                    q = False
            else:
                pass
            i = i+1
        #Return:
        #   obj: <list> lista con un dicionario por señal
        return obj
    else:
        raise ValueError

def extractsignals(cdir, arr_seg):
    #crea una lista de diccionarios con un diccionario(evento) por segmento de la adquisición

    #Argumentos:
    #   cdir: <str> dirección del archivo con adquisición para segmentar
    #   arr_seg: <list> arreglo 2D con la información para segmentar
    f = open(cdir,'r')
    d = [k for k in f]
    dt1 = [k.replace('\n','').split('\t') for k in d if not('#' in k)]
    headdct = json.loads(d[1].replace('#','').replace('\n',''))
    dev = str(headdct.keys()).split("'")[1]
    fs = headdct[dev]['sampling rate']
    sig_names0 =  convertsignames(headdct[dev]['sensor'])
    datelist = list(reversed(headdct[dev]['date'].split('-')))
    dateacq =  datelist[0]+"/"+datelist[1]+"/"+datelist[2]
    ns =  headdct[dev]['resolution']
    n0 = ns[len(ns)-len(sig_names0):len(ns)]
    dt0 = [k[len(ns)-len(sig_names0):len(ns)] for k in dt1]
    dt0 = np.transpose(np.array([[int(i) for i in j] for j in dt0]))
    f.close()
    obj = []
    q = True
    if arr_seg != []:
        samples = [int(round(k[0]*fs,0)) for k in arr_seg]
        samples.append(len(dt0[0])-1)
        print('samples:',samples)
        dt2 = []
        sig_names = []
        tags = []
        n2 = []
        for k in range(0,len(dt0)):
            for s in range(0, len(samples)-1):
                if (arr_seg[s][1].lower()!='eliminar'):
                    dt2.append(dt0[k][samples[s]:samples[s+1]])
                    sig_names.append(sig_names0[k])
                    tags.append(arr_seg[s][1])
                    n2.append(n0[k])
                else:
                    pass
    else:
        sig_names = sig_names0
        dt2 = dt0
        n2 = n0
        tags = ['' for k in sig_names]
    print('segmentos: ',[len(k) for k in dt2])
    print('nombres: ',sig_names)
    print('tags: ',tags)
    print('len(dt2)',len(dt2))
    i = 0
    for name in sig_names:
        if 'PPG' in name:
            ppgobj = ppg(dt2[i],dt2[i],n2[i],fs,convert=True,fnirs=False)
            fIR = ppgobj.irc_uA
            filtered, locs, x1, t1, t2, x2, Xfft, Yfft, VLF, LF, HF, LF_HF, maxima, media, minima = ppgobj.get_HRV() 
            maxima = int(round(maxima,0))
            media = int(round(media,0))
            minima = int(round(minima,0))
            VLF = float(round(VLF,2))
            LF = float(round(LF,2))
            HF = float(round(HF,2))
            LF_HF = float(round(LF_HF,2))
            jsobj = {'fecha':dateacq, 'tag':tags[i], 'med':'PPG', 'datos':{'raw':{'fs':fs,'n':n2[i],'signal':[int(k) for k in dt2[i]]}, '~filt':{'fs':fs, 'signal':[float(k) for k in fIR], 'unid':'(µA)'}, '~HRV@':{'t':[float(k) for k in t1], 'signal':[float(k) for k in x1], 'unid':'(s)'},'extr':{'FC$max.bpm':maxima,'FC$med.bpm':media,'FC$min.bpm':minima, 'VLF.s²':VLF, 'LF.s²':LF, 'HF.s²':HF, 'LF/HF':LF_HF}}}
            obj.append(jsobj)
        elif 'EDA' in name:
            edaobj = eda(dt2[i],n2[i],fs,convert=True)
            C = edaobj.eda_uS
            F = edaobj.filtered_eda
            T =  edaobj.tonic_eda
            PH = edaobj.phasic_eda
            if edaobj.nscr>0:
                SCLm, Am, RTm, T1_2m = edaobj.get_meanresp()
                SCLM, AM, RTM, T1_2M = edaobj.get_maxresp()
                SCLm = float(round(SCLm,2))
                Am = float(round(Am,2))
                RTm = float(round(RTm,2))
                if T1_2m!=None:
                    T1_2m = float(round(T1_2m,2))
                else:
                    T1_2m = 'N'
                SCLM = float(round(SCLM,2))
                AM = float(round(AM,2))
                RTM = float(round(RTM,2))
                FREQ = float(round(edaobj.scr_freq,2))
                if T1_2M!=None:
                    T1_2M = float(round(T1_2M,2))
                else:
                    T1_2M = 'N'
                jsobj = {'fecha':dateacq, 'tag':tags[i], 'med':'EDA', 'datos':{'raw':{'fs':fs,'n':n2[i],'signal':[int(k) for k in dt2[i]]}, '~filt':{'fs':40, 'signal':[float(k) for k in F], 'unid':'(µS)'}, '~T':{'fs':40, 'signal':[float(k) for k in T], 'unid':'(µS)'}, '~PH':{'fs':40, 'signal':[float(k) for k in PH], 'unid':'(µS)'}, 'extr':{'SCL$med.µS':SCLm, 'A$med.µS':Am, 'RT$med.s':RTm, 'T1/2$med.s':T1_2m,'SCL$max.µS':SCLM, 'A$max.µS':AM, 'RT$max.s':RTM, 'T1/2$max.s':T1_2M, 'FREQ.scr/min':FREQ}}}
                obj.append(jsobj)
            else:
                jsobj = {'fecha':dateacq, 'tag':tags[i], 'med':'EDA', 'datos':{'raw':{'fs':fs,'n':n2[i],'signal':[int(k) for k in dt2[i]]}, '~filt':{'fs':40, 'signal':[float(k) for k in F], 'unid':'(µS)'}, '~T':{'fs':40, 'signal':[float(k) for k in T], 'unid':'(µS)'}, '~PH':{'fs':40, 'signal':[float(k) for k in PH], 'unid':'(µS)'}, 'extr':{}}}
                obj.append(jsobj)
        elif 'RESP' in name:
            respobj = pzt(dt2[i],n2[i],fs,convert=True)
            F = respobj.filtered_pzt
            maxR, minR, meanR, freqrpm, t = respobj.getfreq()
            maxR = int(round(maxR,0))
            minR = int(round(minR,0))
            meanR = int(round(meanR,0))
            jsobj = {'fecha':dateacq, 'tag':tags[i], 'med':'RESP', 'datos':{'raw':{'fs':fs,'n':n2[i],'signal':[int(k) for k in dt2[i]]}, '~filt':{'fs':fs, 'signal':[float(k) for k in F], 'unid':'(%)'}, '~FR@':{'t':[float(k) for k in t], 'signal':[float(k) for k in freqrpm], 'unid':'(rpm)'},'extr':{'FR$max.rpm':maxR, 'FR$med.rpm':meanR, 'FR$min.rpm':minR}}}
            obj.append(jsobj)
        elif 'EMG' in name:
            emgobj = emg(dt2[i],n2[i],fs,convert=True)
            F = emgobj.filtered_emg
            meanrms, vrms, percrms, times = emgobj.get_meanrms(0.5)
            meanrms = float(round(meanrms,3))
            jsobj = {'fecha':dateacq, 'tag':tags[i], 'med':'EMG', 'datos':{'raw':{'fs':fs,'n':n2[i],'signal':[int(k) for k in dt2[i]]}, '~filt':{'fs':fs, 'signal':[float(k) for k in F], 'unid':'(mV)'}, 'extr':{'RMS$med.mV':meanrms}}}
            obj.append(jsobj)
        elif 'TEMP' in name:
            tempobj = temp(dt2[i], n2[i], fs, convert=True)
            C = tempobj.temp_C
            F = tempobj.filtered_C
            Fa = tempobj.filtered_F
            maxC, minC, meanC = tempobj.get_C()
            maxF, minF, meanF = tempobj.get_F()
            maxC = float(round(maxC,2))
            minC = float(round(minC,2))
            meanC = float(round(meanC,2))
            maxF = float(round(maxF,2))
            minF = float(round(minF,2))
            meanF = float(round(meanF,2))
            jsobj = {'fecha':dateacq, 'tag':tags[i], 'med':'TEMP', 'datos':{'raw':{'fs':fs,'n':n2[i],'signal':[int(k) for k in dt2[i]]}, '~filt':{'fs':fs, 'signal':[float(k) for k in F], 'unid':'(°C)'}, 'extr':{'$max.°C':maxC, '$min.°C':minC, '$med.°C':meanC, '$max.°F':maxF, '$min.°F':minF, '$med.°F':meanF}}}
            obj.append(jsobj)
        elif 'ECG' in name:
            ecgobj = ecg(dt2[i],10,fs,convert=True)
            C = ecgobj.ecg_mV
            F = ecgobj.filtered_ecg
            t = (1/fs)*np.arange(len(dt2[i]))
            locs = ecgobj.get_rpeaks()
            VLF, LF, HF, LF_HF, t2, x2, t1, x1, arr, rpeaks_peak, Xfft, Yfft = ecgobj.get_HRV()
            maxima,minima,media = ecgobj.get_FC()
            maxima = int(round(maxima,0))
            media = int(round(media,0))
            minima = int(round(minima,0))
            VLF = float(round(VLF,2))
            LF = float(round(LF,2))
            HF = float(round(HF,2))
            LF_HF = float(round(LF_HF,2))
            jsobj = {'fecha':dateacq, 'tag':tags[i], 'med':'ECG', 'datos':{'raw':{'fs':fs,'n':n2[i],'signal':[int(k) for k in dt2[i]]}, '~filt':{'fs':fs, 'signal':[float(k) for k in F], 'unid':'(mV)'}, '~HRV@':{'t':[float(k) for k in t1], 'signal':[float(k) for k in x1], 'unid':'(s)'},'extr':{'FC$max.bpm':maxima,'FC$med.bpm':media,'FC$min.bpm':minima, 'VLF.s²':VLF, 'LF.s²':LF, 'HF.s²':HF, 'LF/HF':LF_HF}}}
            obj.append(jsobj)
        elif 'fNIRS' in name:
            if q:
                channels = headdct[dev]['channels']
                indR = channels.index(9)
                indIR = channels.index(10)
                R = dt2[indR]
                IR = dt2[indIR]
                ppgobj = ppg(R,IR,n2[i],fs,convert=True, fnirs=True)
                fR = ppgobj.rc_uA
                fIR = ppgobj.irc_uA
                filS, intt, Sat, tsat, minS, meanS, maxS = ppgobj.get_spo2()
                jsobj = {'fecha':dateacq, 'tag':tags[i], 'med':'fNIRS', 'datos':{'rawR':{'fs':fs,'n':n2[i],'signal':[int(k) for k in R]}, 'rawIR':{'fs':fs,'n':n2[i],'signal':[int(k) for k in IR]}, '~filtR':{'fs':fs, 'signal':[float(k) for k in fR], 'unid':'(µA)'}, '~filtIR':{'fs':fs, 'signal':[float(k) for k in fIR], 'unid':'(µA)'}, '~SpO2s@':{'t':[float(k) for k in tsat],'signal':[float(k) for k in Sat],'unid':'(%)'},'extr':{'SpO2$max.%':maxS, 'SpO2$med.%':meanS, 'SpO2$min.%':minS}}}
                obj.append(jsobj)
                q = False
        else:
            pass
        i = i+1
    print('obj: ', len(obj))
    #Return:
    #   obj: <list> arreglo con los diccionarios ()
    return obj
    
def tagtodesc(med,tag,c=True):
    tags = {'MRC': {'Descr': 'MRC-Medical Research Council(Fuerza Muscular)', 'mus': 'Músculo', 'sc': 'Score'}, 'EEP': {'Descr': 'EEP-Escala de estrés percibido', 'sc': 'Score'}, 'FSS': {'Descr': 'FSS-Severidad de la Fatiga', 'sc': 'Score'}, 'SF-36': {'Descr': 'SF-36-Instrumento SF-36', 'PF': 'Función física', 'EF': 'Vitalidad', 'LPP': 'Rol físico', 'SF': 'Función social', 'PA': 'Dolor corporal', 'LEP': 'Rol emocional', 'GH': 'Salud general', 'EWB': 'Salud mental'}, 'BDI': {'Descr': 'BDI-Test de Depresión de Beck', 'sc': 'Score'}, 'BAI': {'Descr': 'BAI-Test de Ansiedad de Beck', 'sc': 'Score'}, 'FCSRT': {'Descr': 'FCSRT-Free and Cued Selective Reminding Test', 'RLE1': 'Recuerdo libre ensayo 1', 'RLT': 'Recuerdo libre total', 'RT': 'Recuerdo total'}, 'DIGIT': {'Descr': 'DIGIT SPAN', 'DIR': 'Orden directo', 'INV': 'Orden inverso'}, 'STROOP': {'Descr': 'STROOP-Color word test', 'SP': 'Stroop palabra', 'SC': 'Stroop color', 'II': 'Índice de interferencia', 'PC': 'Stroop palabra-color'}, 'TOL': {'Descr': 'TOL-Torre de Londres', 'TC': 'Puntaje total correcto', 'TM': 'Puntaje total movimiento', 'TE': 'Tiempo total de ejecución', 'TR': 'Tiempo total de resolución', 'TI': 'Tiempo total de inicio'}, 'VS': {'Descr': 'Signos Vitales', 'FC': 'Frecuencia cardiaca', 'FR': 'Frecuencia respiratoria', 'SpO2': 'Saturación de O2', 'S': 'Presión arterial sistólica', 'D': 'Presión arterial diastólica', 'TEMP': 'Temperatura'}, 'ACQ': {'HRV': 'Variabilidad de la frecuencia cardiaca', 'FC': 'Frecuencia cardiaca', 'max': 'Máxima', 'med': 'Media', 'min': 'Mínima', 'VLF': 'Muy baja frecuencia(0-0.04Hz)', 'LF': 'Baja frecuencia(0.04-0.15)', 'HF': 'Alta frecuencia(0.15-0.4Hz)', 'LF/HF': 'Índice LF/HF', 'SCL': 'Conductancia Basal', 'A': 'Amplitud', 'RT': 'Tiempo de elevación', 'T1/2': 'tiempo de recuperación media', 'FREQ': 'SCR/min', 'SCR': 'Respuestas de conductancia de la piel', 'T': 'EDA tónico', 'PH': 'EDA Fásico', 'FR': 'Frecuencia Respiratoria', 'RMS': 'Valor RMS', 'SpO2': 'Saturación de O₂'}}
    if c:
        print('med',med,' tag',tag)
        return tags[med][tag]
    else:
        return tag

def convertsignames(sig_names):
    return [k.replace('HR','PPG (µA)').replace('EDABITREV','EDA (µS)').replace('RESPBIT','RESP(%)').replace('EMGBITREV','EMG(mV)').replace('TMP','TEMP(°C)').replace('ECGBIT','ECG(mV)').replace('hSpO2','fNIRS(µA)') for k in sig_names]

