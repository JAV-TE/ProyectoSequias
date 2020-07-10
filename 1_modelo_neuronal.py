###
#Este programa realiza el relleno con datos esternos 3 variables extenas de
#princenton y 1 o mas var de estaciones de referencia completas para rellenar 
#laa demas restantes
#####

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys


#ruta variable externas entrada----------------------------------------
file1 = '/home/Physics/work/1980_2014PRCP/princenton/reg22_-17.375_-66.125.txt'
file2 = '/home/Physics/work/1980_2014PRCP/series_fill_reg/reg22_fill/CP2775_prcp.scv'
#file3 = '/home/Physics/work/1980_2014PRCP/series_fill_reg/reg22_fill/S2548_prcp.scv'
#file4 = '/home/Physics/work/1980_2014PRCP/series_fill_reg/reg22_fill/P2300_prcp.scv' 


#---archivo que sera rellenado------------------------------------
f_int = '/home/Physics/work/1980_2014PRCP/series_fill_reg/reg22_fill/S2548_prcp.scv'
#ruta salias de resultados----------------------------------------
out =   '/home/Physics/work/1980_2014PRCP/series_fill_reg/reg22_fill/salida/'
#-----------------------------------------------------------------

#obtener nombre ref2------------------------------
base=os.path.basename(file2)
os.path.splitext(base)
nam_ref = str(os.path.splitext(base)[0])
nam_ref = nam_ref[:-5]

#---------------------------------------------------------------
#obtener nombre file3
#base3=os.path.basename(file3)
#os.path.splitext(base3)
#nam_ref3 = str(os.path.splitext(base3)[0])
#nam_ref3 = nam_ref3[:-5]

#get name file4
#base4=os.path.basename(file4)
#os.path.splitext(base4)
#nam_ref4 = str(os.path.splitext(base4)[0])
#nam_ref4 = nam_ref4[:-5]
#---------------------------------------------------------------

#obtener nombre f_int obj
base1=os.path.basename(f_int)
os.path.splitext(base1)
nam_obj = str(os.path.splitext(base1)[0])
nam_obj = nam_obj[:-5]
#------------------------------------------

#read var ext as df
df = pd.read_csv(file1)
df.replace(-999.0, np.NaN, inplace = True)
#crear data frame con fechas#############################
dates = pd.date_range(start='1980', end='2015', freq="M")
#data frame  df 
df.set_index(dates, inplace = True)
df.index.names = ['Tiempo']

#####read file2 y file3
est = pd.read_csv(file2 , header = None, delimiter = '\t')
est.replace(-99.9, np.NaN, inplace  =True)
est.head()
df[nam_ref] = est[1].values

#-----------------------------------------------------------------
#est3 = pd.read_csv(file3 , header = None, delimiter = '\t')
#est3.replace(-99.9, np.NaN, inplace  =True)
#est3.head()
#df[nam_ref3] = est3[1].values

#est4 = pd.read_csv(file4 , header = None, delimiter = '\t')
#est4.replace(-99.9, np.NaN, inplace  =True)
#est4.head()
#df[nam_ref4] = est4[1].values
#-----------------------------------------------------------------

list_var = list(df)
print('lISTA DATA FRAME COLUMNS X[]')
print(list_var)
#####read file obj
obj = pd.read_csv(f_int, header = None, delimiter = '\t')
obj.replace(-99.9, np.NaN, inplace  =True)
df[nam_obj] = obj[1].values



################################################################33
corr = df.corr(method = 'pearson')
plt.figure(figsize=(10,6))
sns.heatmap(corr,annot=True, fmt=".2f",cmap = 'RdBu',linewidth= 0.5)
plt.title('Correlaciones_Prcp[mm]_'+nam_obj)
plt.savefig(out+nam_obj+'_corr.pdf')
plt.close()

###################################################################


#CREAMOS VAR X[] y Y[] y modelamos
df1 = df.dropna()  #sin Nas para todos
X_train = np.asanyarray(df1[list_var]).astype('float32')
Y_train = np.asanyarray(df1[[nam_obj]]).astype('float32')


#argumrnto minimo sin estandarizar
ar_min  =np.sort(Y_train, axis = None)
a_min = ar_min[0]

#Ahora normalizamos-------------------------------------------------
from sklearn.preprocessing import MinMaxScaler


#Rescaling (Min-Max normalization) formule:
#x'= (x-min(x))/(max(x)-min(x))

scaler_x = MinMaxScaler(feature_range=(0, 1))
scaler_x = scaler_x.fit(X_train)
x_train = scaler_x.transform(X_train)

scaler_y = MinMaxScaler(feature_range=(0, 1))
scaler_y = scaler_y.fit(Y_train)      #esta config servira para llevar a datos origami
y_train = scaler_y.transform(Y_train) 

#argumento minimo estadarizado
ar_minest  =np.sort(y_train, axis = None)
a_minest = ar_minest[0]


#######################SKLEARN-KERAS##PROCESING###################################
#######################-----------------------####################################
from sklearn import linear_model
regr = linear_model.LinearRegression()
regr.fit (x_train,y_train)

#%%%%%%%%%%%%%%%%%%%%%%GUARDAR datos en txt%%%%%%%%%%%%%%
#abrir cadena{						%
orig_stdout = sys.stdout                                #
f = open(out + nam_obj+'_info.txt', 'w')          # 
sys.stdout = f  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print('--------------------Coefients SKlearn----------------------------')
print('Coefients:', regr.coef_)
print('-----------------------------------------------------------------')
###################################
y_hat = regr.predict(x_train)     #
#e#################################
y_hat_sk =np.reshape(y_hat ,(len(y_hat)))   #predecidos 1D

y_true = np.reshape(y_train,(len(y_train))) #verdaderos 1D

print('Suma errores residuales: %.5f'
       % np.mean((y_hat_sk - y_true)**2))
print('Score Variance : %.2f'  % regr.score(x_train,y_train))
print('------------------------CORR./y_hat-y_true/----------------------')
print(np.corrcoef(y_true, y_hat_sk))
print('-----------------------------------------------------------------')

#PLOT RESULTADOS#---###---###----#################
plt.figure(figsize=(9,4))
plt.plot(df1.index,y_true, color = 'black',label = 'Valores verdaderos', linewidth=0.4)
plt.plot(df1.index,y_hat_sk,label='Predichos RLM',color = 'red',linewidth=0.5)
plt.xlabel('Tiempo[meses]')
plt.ylabel('Precipitación[mm]')
#plt.title('keras, model')
plt.legend()
plt.savefig(out+nam_obj+'_RML.pdf')
plt.close()

#######################------keras---####################################
nx_var = len(list_var) #numero de X[]vars

import keras
from keras.models import Sequential
from keras.layers import Dense

model = Sequential()

#Aquí necesitas especificar las épocas (epochs); estos son el número de iteraciones para que el #proceso de entrenamiento se ejecute a través del conjunto de datos y el tamaño del grupo, que es el #número de instancias que se evalúan antes de una actualización de peso.

model.add(Dense(27, activation='relu', input_shape=(nx_var,)))
model.add(Dense(27, activation='linear'))
model.add(Dense(27, activation='linear'))
model.add(Dense(1,activation='relu'))
model.summary()

#relu, linear , linear, relu

model.compile(loss='mean_squared_error',
                optimizer='adam',
                metrics=["mean_absolute_error"])
model.fit(x_train,y_train,epochs = 21,batch_size=32, verbose = 0)

#predecimos con los valores modelados
#################################
y_pred = model.predict(x_train) #
#################################

y_predic = np.reshape(y_pred,(len(y_pred))) #verdaderos 1D

#evaluacion metrics
result = model.evaluate(x_train,y_train)
print(model.metrics_names)
print(result)
print('----------------Correlacion--------')
print(np.corrcoef(y_true, y_predic))
print('-----------------------------------')

#PLOT RESULTS#---###---###----#################
plt.figure(figsize=(9,4))
plt.plot(df1.index,y_true, color = 'blue',label = 'Valores verdaderos', linewidth=0.4)
plt.plot(df1.index,y_predic,label='Predichos DNN',color = 'red',linewidth=0.5)
plt.xlabel('Tiempo[meses]')
plt.ylabel('Precipitación[mm]')
plt.legend()
plt.savefig(out+nam_obj+'_SP.pdf')
plt.close()

##########DROP NAN values solo en X[]########
df2 = df[list_var]
df2 = df2.dropna() #aunque no hay Nas pero quizas en un futuro si

X_obj  = np.asanyarray(df2.astype('float32'))
obj_x = scaler_x.transform(X_obj) 

#predescimos para datos solo x[]
y_fin = model.predict(obj_x)   #
################################
ni = len(y_fin)
y_final = np.zeros([ni,1],dtype='f')


for i in range(ni):
    if y_fin[i] < 0:
        y_final[i,0] = a_minest
    else :
        y_final[i,0] = y_fin[i]

#transformacion inversa
final_y  = scaler_y.inverse_transform(y_final)
fin_y =np.reshape(final_y ,(len(final_y))) #como array¿

######################################################
####aqui falta un lugar para interar concaTENAR datos originales#

df2[nam_obj+'_pred'] = fin_y

print(df.corr())

original = df[[nam_obj]]
predict  = df2[[nam_obj+'_pred']]

#UNir los data frames y ori y y_prd
fin = pd.concat([original,predict], axis = 1)

#cremos un solo vector concatenando
n = len(fin[nam_obj])
#obtener valores
ori_data = fin[nam_obj].values
pred_data  = fin[nam_obj+'_pred'].values

y_fill = np.zeros([len(ori_data)])
for i in range(len(ori_data)):
    if np.isnan(ori_data[i]):
        y_fill[i] = pred_data[i]
    else :
        y_fill[i] = ori_data[i]


#Almacenamos en el dataframe fin el concatenado
df[nam_obj+'_fill'] =np.around(y_fill,3)
#PLOT DEL RESULTADOs
plt.figure(figsize=(12,12))
plt.plot(df[nam_obj+'_fill'], label = 'Rellenado' ,linewidth=0.3, color = 'red')
plt.plot(df[nam_obj], label = 'Original' ,linewidth=0.6, color = 'black')
plt.legend()
plt.savefig(out+ nam_obj+'_fill.pdf')


#guardar rellenado en csv
df.replace(np.NaN,-99.9,inplace = True)
df[nam_obj+'_fill'].to_csv(out+nam_obj+'_prcp.scv', index = True, sep = '\t',header = None)
#-------------------------------------------------------------------------
sys.stdout = orig_stdout
f.close()
print('------------------------------------------------------------------------')
print('SE HA RELLENADO EXITOSAMENTE LA ESTACION:'+nam_obj)
