import pandas as pd
import numpy as np
from numpy import unravel_index
from sklearn.cluster import KMeans
from collections import Counter, defaultdict
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# load combine data - 4 years
data = pd.read_csv('CombinedData.csv') 
yeardays= 1439

df=data.copy()
df['datetime']=pd.date_range(start='2015-01-01 00:00:00', end='2018-12-09 23:00:00', freq='60min')
df=df.set_index('datetime')
df.head()

# find peak electricity and gas day(s)
print('electricity peak day index =', df['Elec'].idxmax())
gas_cols_Domestic=['EA_Domestic', 'EM_Domestic', 'NE_Domestic', 'NO_Domestic', 'NT_Domestic', 'NW_Domestic', 'SC_Domestic', 'SE_Domestic', 'SO_Domestic', 'SW_Domestic', 'WM_Domestic', 'WN_Domestic', 'WS_Domestic']
print('gas_Domestic peak day index =', df[gas_cols_Domestic].sum(axis=1).idxmax())
gas_cols_Services=['EA_Services', 'EM_Services', 'NE_Services', 'NO_Services', 'NT_Services', 'NW_Services', 'SC_Services', 'SE_Services', 'SO_Services', 'SW_Services', 'WM_Services', 'WN_Services', 'WS_Services']
print('gas_Services peak day index =', df[gas_cols_Services].sum(axis=1).idxmax())
gas_cols_Industrial = ['EA_Industrial', 'EM_Industrial', 'NE_Industrial', 'NO_Industrial', 'NT_Industrial', 'NW_Industrial', 'SC_Industrial', 'SE_Industrial', 'SO_Industrial', 'SW_Industrial', 'WM_Industrial', 'WN_Industrial', 'WS_Industrial']
print('gas_Industrial peak day index =', df[gas_cols_Services].sum(axis=1).idxmax())


df_nopeak=df.copy()
df_nopeak=df_nopeak.loc[df_nopeak.index.date!=df['Elec'].idxmax().date()]
#df_nopeak.describe()

date_peak = pd.to_datetime(df['Elec'].idxmax()).date()

# take only one peak day out from clustering-electricity peak day

days = [data[i].values.reshape(yeardays,24) for i in data] # yeardays days x 24 hours array of the days in a year
peakELEC=unravel_index(data['Elec'].values.reshape(yeardays,24).argmax(), data['Elec'].values.reshape(yeardays,24).shape)[0] #
DaysNoPeak = np.delete(days,peakELEC,axis=1) #Days without peak electricity day

## no need to take out second highest heat demand day
# peakHEAT=unravel_index(np.asarray(DaysNoPeak[7:20]).sum(axis=0).argmax(), np.asarray(DaysNoPeak[7:20]).sum(axis=0).shape)[0] -1 
# DaysNoPeak = np.delete(DaysNoPeak,peakHEAT,axis=1)#Days without peak electricity and gas day

DaysMax = [DaysNoPeak[i].max() for i in range(0,len(days))] #max non peak day for each attribute 
DaysNormalised= [DaysNoPeak[i]/DaysMax[i] for i in range(0,len(DaysNoPeak))] #Normalised attributes for all the remaining days (0,1)

print('# of days in DaysNoPeak/DaysNormalised =', len(DaysNormalised[0]))

data_new = np.concatenate((DaysNormalised), axis=1) #array with rows each day and columns each hour of each attribute (yeardays, 24*70)
days_new = np.concatenate((days),axis=1) #Some as above but not normalised 

print(data_new.shape, days_new.shape)

# set number of clusters
# Calling the Kmeans Clustering 
nClusters= 5

### original ###
# kmeans = KMeans(nClusters,init='k-means++', max_iter=2000).fit(data_new)

# somehow I need to set random seed INSIDE KMeans to ensure reproducibility # 
kmeans = KMeans(nClusters,init='k-means++', max_iter=2000, random_state=40).fit(data_new)

centroids=kmeans.cluster_centers_ #centroids of each cluster: Array with (NClusters, 24*70)

#sillouette measure the less the better :)
print('inertia=', kmeans.inertia_)

#Using the "labels_" attribute the population of the days can be counted  
nDays = Counter(kmeans.labels_) #how many days in each cluster
z={i: np.where(kmeans.labels_ == i)[0] for i in range(kmeans.n_clusters)} #Dictionary of the days within each cluster 
nDays = sorted(nDays.items()) #sort nDays by key, i.e. cluster 

weights = [];
[weights.append(nDays[v][1]) for v in range(0,len(nDays))]

DaysCentroids = [centroids[:,i*24:24*(1+i)] for i in range(0,len(DaysNoPeak))] # Array (70,30,24)

centroids_means = np.vstack([DaysCentroids[i].mean(axis=1) for i in range(0, len(DaysNoPeak))])
items_mean = np.vstack([DaysNormalised[i].mean(axis=1) for i in range(0, len(DaysNoPeak))])

def representative_day(centroid_means,z, nClusters):
    yy=[]
    for clu in range(0,nClusters):
        centroids_means.transpose()[clu]
        # focus on gas columns
        X = centroids_means.transpose()[clu][7:20]-items_mean.transpose()[z[clu]][:,7:20]
        #X = centroids_means.transpose()[clu] -items_mean.transpose()[z[clu]] 
        #yy.append(z[clu][np.abs(X.sum(axis=1)).argmin()])
        yy.append(z[clu][np.where(X.sum(axis=1) > 0, X.sum(axis=1), np.inf).argmin()])
    return yy
rep_days=representative_day(centroids_means,z, nClusters)

# add peak day(s)
rep_days.insert(0,peakELEC) #insert peak elec day 
weights.insert(0,1) # append weight 1 for Peak elec Day
Clusters= pd.DataFrame({"Cluster_iD":range(nClusters+1),"Repres_Days":rep_days,"Weights": pd.Series(weights)})

Clusters["Cluster_iD"]=range(nClusters+1)
DaysFinal= [pd.DataFrame(days_new[rep_days,i*24:24*(1+i)], 
                         index=range(len(rep_days)),columns=range(24)) for i in range(0,len(DaysNormalised))]

peak_ids=[0] # only one peak

M_profiles=[]
for iD in Clusters['Cluster_iD']:
    rep_day=Clusters.loc[Clusters['Cluster_iD']==iD, 'Repres_Days'].values[0]
    if iD in peak_ids:
        _M=df.iloc[rep_day*24:(rep_day+1)*24]
    else:
        _M=df_nopeak.iloc[rep_day*24:(rep_day+1)*24]
    M_profiles.append(_M)

df_ClusteredProfiles = pd.DataFrame(
    index=tuple(zip(np.array([np.repeat(i, 24) for i in range(nClusters+1)]).ravel(),
                             np.tile(range(24),nClusters+1))), 
    columns=df.columns,
    data= pd.concat(M_profiles, axis=0).values)

writer = pd.ExcelWriter(str(20152018)+'Cluster_'+str(nClusters)+'.xlsx')
Clusters.to_excel(writer, 'NDaysCluster', index=False)
df_ClusteredProfiles.to_excel(writer, sheet_name='ClusteredProfiles', index=True)
writer.save()

# load 2035 data
data_2035 = pd.read_csv('Hydrogen_demand_2050.csv') 
yeardays= 1439

df = data_2035.copy()
df['datetime']=pd.date_range(start='2015-01-01 00:00:00', end='2018-12-09 23:00:00', freq='60min')
df=df.set_index('datetime')
#df = df[~(df.index.date == date_peak)]
#df.head()

df_nopeak=df.copy()
df_nopeak=df_nopeak.loc[df_nopeak.index.date!=date_peak ]

#Extract peak day data

df_peakday = df.loc[df.index.date == date_peak]

columns = ['EA', 'EM', 'NE', 'NO', 'NT', 'NW', 'SC', 'SE', 'SO', 'SW', 'WM', 'WN', 'WS']
matrices = {}
df_peakday_data = {}
    
for column in columns:
    peakday_data_column = df_peakday[column].values
    new_rows = [peakday_data_column[i:i+24] for i in range(0, len(peakday_data_column), 24)]
    matrix = pd.DataFrame(new_rows)
    df_peakday_data[column] = matrix     

peakday = pd.DataFrame({'Cluster_ID': [1]*24, 'Row_Number': list(range(1, 25))})

for column in columns:
    # 从 df_peakday 中提取列数据
    peakday_data_column = df_peakday[column].values

    # 将每个 column 的数据添加到 combined_df DataFrame 中
    peakday[column] = peakday_data_column

# save to Excel
peakday.to_csv('peakday_data.csv', index=False)

days_2035 = [data_2035[i].values.reshape(yeardays,24) for i in data_2035] # yeardays days x 24 hours array of the days in a year
DaysNoPeak_2035 = np.delete(days_2035,peakELEC,axis=1) #Days without peak electricity day



for cluster_index, date_indices in z.items():

    cluster_data_2035 = pd.DataFrame()
    

    for date_index in date_indices:
        date_data = df_nopeak.iloc[date_index*24 : (date_index+1)*24]
        cluster_data_2035 = pd.concat([cluster_data_2035, date_data])
    

    filename = f"cluster_{cluster_index+2}.csv"
    

    cluster_data_2035.to_csv(filename, index=False)

acculated_clusters = {}  # Using a dictionary to store DataFrames
df_result_ave= pd.DataFrame() 
df_para_result_neg= pd.DataFrame()      
df_para_result_pos= pd.DataFrame()
combined_ave_df = pd.DataFrame()
combined_para_result_neg = pd.DataFrame()
combined_para_result_pos = pd.DataFrame()

for id in range(2, cluster_index+3): #从2开始排序
    filename = f"cluster_{id}.csv"  # Corrected file name string
    acculated_clusters[id] = pd.read_csv(filename)  # Storing DataFrame in a dictionary
    data_cluster = acculated_clusters[id]
    columns = ['EA', 'EM', 'NE', 'NO', 'NT', 'NW', 'SC', 'SE', 'SO', 'SW', 'WM', 'WN', 'WS']
    matrices = {}

    for column in columns:
        column_data = data_cluster[column].values
        new_rows = [column_data[i:i+24] for i in range(0, len(column_data), 24)]
        matrix = pd.DataFrame(new_rows)
        matrices[column] = matrix 
    
    for column in columns:
        new_variable_name = f"Gas_{column}"
        globals()[new_variable_name] = matrices[column]
        
    Quantile_results = {}
    u_aveS= {}
    for column in columns:
        globals()[f"u_ave_{column}"] = np.mean(globals()[f"Gas_{column}"], axis=0) #Calculate the average of each column
        #u_aveS[column] = np.mean(globals()[f"Gas_{column}"], axis=0)
        globals()[f"Num_{column}"] = globals()[f"Gas_{column}"].shape[0] 
        globals()[f"X_0_{column}"] = globals()[f"Gas_{column}"] - np.tile(globals()[f"u_ave_{column}"], (globals()[f"Num_{column}"], 1))
        globals()[f"S_{column}"]  = (1 / (globals()[f"Num_{column}"] - 1)) * globals()[f"X_0_{column}"].T @ globals()[f"X_0_{column}"]  # 协方差矩阵
        globals()[f"eigenvalues_{column}"], globals()[f"eigenvectors_{column}"] = np.linalg.eig(globals()[f"S_{column}"])
        globals()[f"Eigv_{column}"] = np.transpose(globals()[f"eigenvectors_{column}"])  # Each row represents an eigenvector
        globals()[f"T_{column}"] = globals()[f"Eigv_{column}"] @ globals()[f"X_0_{column}"].T    
        df_result_ave[column] = globals()[f"u_ave_{column}"]   
        
    df_result_ave['Cluster_ID'] = id
    df_result_ave['Row_Number'] = np.arange(1, 25)  # 生成从 1 到 24 的序列   
    df_result_ave = df_result_ave[['Cluster_ID', 'Row_Number'] + columns]
    combined_ave_df = pd.concat([combined_ave_df, df_result_ave], axis=0)
    
    #Kernel density estimation, computes the probability density function    
    alpha = 0.01
    kde={}
    x={}
    cdf={}
    for column in columns:
        num_attributes_column = globals()[f"T_{column}"].shape[0]
        globals()[f"num_attributes_{column}"] = num_attributes_column
        globals()[f"lower_quantiles_{column}"] = np.zeros((num_attributes_column,1))
        globals()[f"upper_quantiles_{column}"] = np.zeros((num_attributes_column,1))  
        for i in range(num_attributes_column):
            kde[column] = gaussian_kde(globals()[f"T_{column}"].iloc[i]) # 为每一行属性计算核密度函数
            x[column] = np.linspace(np.min(globals()[f"T_{column}"].iloc[i]), np.max(globals()[f"T_{column}"].iloc[i]), 10000)   
            cdf[column] = np.zeros_like(x[column])     
       #Calculate the CDF(cumulative density function) of the probability density function
            for j, vall in enumerate(x[column]):
                cdf[column][j] = kde[column].integrate_box_1d(-np.inf, vall)     
       # 找到满足条件的最小数据点 lower bound
            for j, vall in enumerate(x[column]):
                if cdf[column][j] >= alpha:
                    globals()[f"lower_quantiles_{column}"][i] = vall
                    break
            else:
                    globals()[f"lower_quantiles_{column}"][i] = np.min(globals()[f"T_{column}"].iloc[i]) 
                    
            for j, vall in enumerate(x[column]):
                if cdf[column][j] >= 1- alpha:
                    globals()[f"upper_quantiles_{column}"][i] = vall
                    break
            else:
                    globals()[f"upper_quantiles_{column}"][i] = np.max(globals()[f"T_{column}"].iloc[i])        
    
    #Calculate xi 
    for column in columns:
        num_attributes = globals()[f"num_attributes_{column}"]  #  num_attributes_{column} 
        e_column = np.ones((num_attributes, 1))  
        eigenvectors_att = globals()[f"eigenvectors_{column}"]
        lower_quantiles_att = globals()[f"lower_quantiles_{column}"]
        up_quantiles_att = globals()[f"upper_quantiles_{column}"]
        PZ_lb_att = np.multiply (eigenvectors_att,  np.ones((num_attributes, 1))@lower_quantiles_att.T) 
        PZ_ub_att = np.multiply (eigenvectors_att,  np.ones((num_attributes, 1))@up_quantiles_att.T)
        PZ_att =  np.hstack((PZ_lb_att, PZ_ub_att)) #24*48 matrix
        PZ_att_trans = np.transpose(PZ_att) #48*24 matrix
        PZ_lb_att_trans =np.transpose(PZ_lb_att) # 24*24 matrix
        PZ_ub_att_trans =np.transpose(PZ_ub_att)
        Colum_PZ_att = PZ_att_trans.reshape(-1, order='F') # stack 1152*1
        Colum_PZ_lb_att = PZ_lb_att_trans.reshape(-1, order='F')
        Colum_PZ_ub_att = PZ_ub_att_trans.reshape(-1, order='F')
        df_para_result_neg[column] = Colum_PZ_lb_att
        df_para_result_pos[column] = Colum_PZ_ub_att
        
     
    df_para_result_neg['Cluster_ID'] = id  # Cluster ID 
    df_para_result_neg['Sub_Row_Number'] = np.tile(np.arange(1, 25), 24)  #  1  24     
    df_para_result_neg['Row_Number'] = np.repeat(np.arange(1, 25), 24)  # 1  24 
    df_para_result_neg = df_para_result_neg[['Cluster_ID','Row_Number','Sub_Row_Number'] + columns]
    df_para_result_pos['Cluster_ID'] = id
    df_para_result_pos['Sub_Row_Number'] = np.tile(np.arange(1, 25), 24)       
    df_para_result_pos['Row_Number'] = np.repeat(np.arange(1, 25), 24)  
    df_para_result_pos = df_para_result_pos[['Cluster_ID','Row_Number','Sub_Row_Number'] + columns] 
    combined_para_result_pos =  pd.concat([combined_para_result_pos, df_para_result_pos], axis=0)
    combined_para_result_neg =  pd.concat([combined_para_result_neg, df_para_result_neg], axis=0)   
    
combined_para_result_neg.to_csv('BigQneg_2050.csv', index=False)
combined_para_result_pos.to_csv('BigQpos_2050.csv', index=False)  
combined_ave_df.to_csv('ave_clusters.csv', index=False)

# read data
peakday_data = pd.read_csv('peakday_data.csv')
ave_clusters = pd.read_csv('ave_clusters.csv')

# Concatenate data by row
final_average_data = pd.concat([peakday_data, ave_clusters], axis=0)

# save to Excel
final_average_data.to_csv('final_average_data_2050.csv', index=False)