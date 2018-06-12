#!/usr/bin/env python3
#!/usr/bin/python3
# Evaluation of number of clusters for different timesteps in a xyz file.

from fnc import *
#from gui import *
ppo = 9
box = 30
conc = 0.05
frameskip = 25
flag_hist = 0
flag_ave = 0
flag_pbc = 1
interval = 250000
RHO 			= 3
FRAME_COLLECTION 	= 500
PLURCHAIN		= 15
timecounter 		= []
CheckCluster 		= []
check_status 		= 0

import sklearn
import numpy as np
from sklearn import cluster
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import errno
import subprocess
import os
import math as mth
from math import sqrt
import traceback
import sys
import pylab as pl
from scipy.spatial.distance import pdist, squareform
import scipy
from collections import Counter
sys.tracebacklimit = None 
filename	= "video.xyz"
times 		= [TimeStep([],[],[],[],[],[]) for _ in range(5000)]    
print ("I am collecting coordinates")
xBeads 	= box*box*box*RHO*conc/PLURCHAIN
Beads 	= xBeads * ppo
coordinates = []

xyz = open(filename,"r")
coordinate_counter = 0
for line in xyz:
        try:
            atom, x, y, z 	= line.split()
            index 		= int (coordinate_counter / int(Beads))
            times[index].add_atom(atom, x, y, z)
            coordinate_counter += 1
        except ValueError:      
            pass
xyz.close()
box_norm = 0
clusterSize = []
clusterGyrRad = []
print ("I am creating clusters and plotting ...")
for j in range (1000,2000,frameskip):
	
	coordinates 		= np.column_stack((times[j].xCoord,times[j].yCoord,times[j].zCoord))
	coordinates_norm 	= coordinates/box
	distance_matrix = [[0 for x in range(len(coordinates[0]))]for y in range (len(coordinates[0]))]
	if flag_pbc == 1:	
		box_norm = box / box
		distance_matrix = pbc_distance_matrix(coordinates, box)
		db = DBSCAN(eps =2.2, metric='precomputed').fit(distance_matrix)
		core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
		core_samples_mask[db.core_sample_indices_] = True
		labels = db.labels_	

	else:
		box_norm = box / box
		db = DBSCAN(eps=2.2, metric='euclidean', algorithm='auto', p=None).fit(coordinates)			
		core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
		core_samples_mask[db.core_sample_indices_] = True
		labels = db.labels_
	current_clusters = len(set(labels))
	n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
	import matplotlib.pyplot as plt
	fig = plt.figure()

	if flag_hist == 1:
		ax = fig.add_subplot(121, projection ='3d')	
		bx = fig.add_subplot(122)	
	else:	
		ax = fig.add_subplot(111, projection='3d')
	unique_labels = set(labels)
	colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
	for k, col in zip(unique_labels, colors):
    		if k == -1:
 
        		col = 'k'

    		class_member_mask = (labels == k)
    		xy = coordinates[class_member_mask & core_samples_mask]
		
    		ax.plot(xy[:, 0], xy[:, 1],xy[:,2],'o', markerfacecolor=col,
			markeredgecolor='k', markersize=14)
    		#print ("################ Inizio Gyration #####################")		
    		#print ("XCoord: ")
    		#print (xy[:,0])
    		if (len(xy[:,0]) != 0) :
    			xCm = (sum(xy[:,0])/len(xy[:,0]))
    			xDistCm = [x - xCm for x in xy[:,0]]
    			xDistCmSquared = [i ** 2 for i in xDistCm]
    		
    		#print (X_dist_cm)	
    		#print (X_dist_cm_squared)
    		#print ("X_cm: ")
    		#print (X_cm)
    		#print ("YCoord: ")
    		#print (xy[:,1])
    			yCm = (sum(xy[:,1])/len(xy[:,1]))
    			yDistCm = [y - yCm for y in xy[:,1]]
    			yDistCmSquared = [i ** 2 for i in yDistCm]
    		#print ("Y_cm: ")
    		#print (Y_cm)
    		#print ("ZCoord: ")		
    		#print (xy[:,2])
    			zCm = (sum(xy[:,2])/len(xy[:,2]))
    			zDistCm = [z - zCm for z in xy[:,2]]
    			zDistCmSquared = [i ** 2 for i in zDistCm]
    			xyzDistCmSquared = [x+y+z for x,y,z in zip(xDistCmSquared, yDistCmSquared,zDistCmSquared)]
    			#xyzDistCmSquaredCorrected = [x*3/5 for x in xyzDistCmSquared]
    			xyzDistCmSquaredAv = (sum(xyzDistCmSquared[:]) / len(xyzDistCmSquared[:]))
    			gyrationRadius = sqrt(xyzDistCmSquaredAv)
    			#xyzDistCmRooted = [sqrt(x) for x in xyzDistCmSquared] 
    			#gyrationRadius = (sum(xyzDistCmRooted)/len(xyzDistCmRooted))
    			if (len(xy[:,0]) != 0) and (gyrationRadius > 1):		
    				clusterSize.append(len(xy[:,0]))	    			
    				clusterGyrRad.append(gyrationRadius)
			
    			#print (clusterSize)
    			#print (clusterGyrRad) 
    		#		print (gyrationRadius)
    		#print (X_dist_cm_squared)
    		#print (Y_dist_cm_squared)
    		#print (Z_dist_cm_squared)
    		#print (XYZ_dist_cm_squared)
    		#print ("Y_cm: ")
    		#print (Y_cm)
    		#print ("################ Fine Gyration #####################")
    		times[j].beadHisto.append(len(xy[:,0]))
    		if times[j].beadHisto != 0:							
    			times[j].chainHisto.append(round((len(xy[:,0]) / ppo)))	
    		ax.plot(xy[:, 0], xy[:, 1],xy[:,2],'o', markerfacecolor=col,
			markeredgecolor='k', markersize=6)
    		check_status = check_status + (len(xy[:,0]))
	if flag_hist == 1:	
		normalized_bin = len(times[j].chainHisto)
		bx.hist(times[j].chainHisto, histtype='stepfilled')
	CheckCluster.append(n_clusters_)
	timecounter.append(j*FRAME_COLLECTION)
	plt.ylim((0,30))
	plt.xlim((0,30))
	#plt.zlim((0,30))		
	plt.title('Estimated number of clusters: %d' % n_clusters_)
	plt.savefig('timestep_'+str((1+j)*FRAME_COLLECTION)+'.png')
	plt.close()
	if (j >= (index-frameskip)):	
		fig_2 = plt.figure()
		cx = fig_2.add_subplot(111)
		cx.plot(timecounter,CheckCluster,'o',markersize=6,ls='--')
		plt.title('#Cluster Vs Time')
		plt.xlabel('Timestep')
		plt.ylabel('#Cluster')		
		plt.savefig('cluster_'+str(j*500)+'.png')
		plt.close()
	fig_3 = plt.figure()
	gx = fig_3.add_subplot(111)
	gx.plot(clusterSize,clusterGyrRad,'o',markersize=5)
	plt.title('Gyr Rad')
	plt.xlabel('Size')
	plt.ylabel('GyrRad')		
	plt.savefig('gyrRadCluster_'+str(j*500)+'.png')
	plt.close()
from scipy.optimize import curve_fit


clusterSize = np.asarray([mth.log(n) for n in clusterSize])

clusterGyrRad = np.asarray([mth.log(t) for t in clusterGyrRad])

from pylab import * 

m,b = polyfit(clusterSize,clusterGyrRad,1)


print(len(clusterSize))
print(len(clusterGyrRad))
print (m)
print (b)
fig_4 = plt.figure()
plt.plot(clusterSize,m*clusterSize+b,'--k')
plt.plot(clusterSize,clusterGyrRad,'o',markersize=5)
plt.title('Gyr Rad Log')
plt.xlabel('Size')
plt.ylabel('GyrRad')		
plt.savefig('LogGyrRad.png')
plt.close()
print ("End of the Cluster Analysis")
print ("Starting Histogram Averaging")
if (flag_ave == 1 ):
	average_histograms_in_time(times, 500000, 990000,frameskip,FRAME_COLLECTION)
#	average_histograms_in_time(times, 500000, 700000,frameskip,FRAME_COLLECTION)
#	average_histograms_in_time(times, 850000, 1050000,frameskip,FRAME_COLLECTION)
print ("Check Your Results in the 'Results' Folder")



