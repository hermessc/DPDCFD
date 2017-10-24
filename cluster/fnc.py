#Class: 	TimeStep  
#Comment:	It contains positions for different PPO atoms at different timesteps)
#Param:
#		@number: the "label of the atoms"
#		@xCoord: x position of one atom
#		@yCoord: y position of one atom 
#		@zCoord: z position of one atom
#		@beadHisto: number of beads per cluster
#		@chainHisto: number of chains per cluster
class TimeStep(object):
	number = []
	xCoord = []
	yCoord = []
	zCoord = []
	beadHisto = []
	chainHisto = []

#Function: 
#Comment:	Initialize the object. It assigns the lists to the object
	def __init__(self,number,xCoord,yCoord,zCoord,beadHisto,chainHisto):
		self.number = number
		self.xCoord = xCoord
		self.yCoord = yCoord
		self.zCoord = zCoord
		self.beadHisto = beadHisto
		self.chainHisto = chainHisto

	
#Function:
#Comment:	#Function to assign the values read from the file to the objects
#Param:
#	@self: reference object
#	@atoms: label of the atom
#	@xCoord: x coordinate
#	@yCoord: y coordinate
#	@zCoord: z coordinate
	def add_atom(self,atoms,mxCoord,myCoord,mzCoord):
		self.number.append(atoms)
		self.xCoord.append(float(mxCoord))
		self.yCoord.append(float(myCoord))
		self.zCoord.append(float(mzCoord))

#Function:
#Comment:	It takes a list and the number of clusters and create a small matrix which counts all the occurrencies of a value into the list
#Param:
#	@histo: list containing all the cluster sizes
#	@size: number of clusters
#Return:
	#@matrix: the matrix counting the occurrencies
def AverageManager(histo,size,matrix):											
	columns = 2
	for cluster_size in histo:
		for mini_row in range(0,size):
			if (matrix[mini_row][0] != cluster_size) and (matrix[mini_row][0] == 0):			
				matrix[mini_row][0] = cluster_size
				matrix[mini_row][1] += 1
				break
			elif (matrix[mini_row][0] == cluster_size):
				matrix[mini_row][1] += 1;
				break
	return (matrix)


def sortingHistograms(normalized_histogram_x, normalized_histogram_y):
	ordered_normalized_histogram_y = []
	ordered_normalized_histogram_x = sorted(normalized_histogram_x)
	for x in range(0,len(ordered_normalized_histogram_x),1):
			for y in range(0, len(normalized_histogram_x),1):
				if (normalized_histogram_x[y] == ordered_normalized_histogram_x[x]):
					ordered_normalized_histogram_y.append(normalized_histogram_y[y])
	return ordered_normalized_histogram_y

def reBuildingHistograms(ordered_normalized_histogram_x, ordered_normalized_histogram_y):
	reBuiltHisto = []
	for x_element in range(0,len(ordered_normalized_histogram_x),1):
		bin_size = int(round(ordered_normalized_histogram_y[x_element] * 10)) 
		for y_element in range(0, bin_size, 1):
			reBuiltHisto.append(ordered_normalized_histogram_x[x_element])
	return reBuiltHisto


def pbc_distance_matrix(coordinates,box):
	import numpy as np
	xdistance = 0.0
	ydistance = 0.0
	zdistance = 0.0
	distance = 0.0
	matrix_dimension = len(coordinates)-1	
	sq_dist = 0.0
	j=0
	distance_matrix = [[0 for x in range(matrix_dimension)]for y in range (matrix_dimension)]
	for i in range (0,matrix_dimension,1):
		j=0
		if i != j:
			while j < i:
				xdistance = coordinates[i][0] - coordinates[j][0]
				if xdistance > 0.5*box:
					xdistance -=box
				ydistance = coordinates[i][1] - coordinates[j][1]
				if xdistance > 0.5*box:
					ydistance -=box
				zdistance = coordinates[i][2] - coordinates[j][2]
				if xdistance > 0.5*box:
					zdistance -=box
				sq_dist = xdistance**2 + ydistance**2 + zdistance**2
				distance = np.sqrt(sq_dist)	
				distance_matrix[i][j] = distance
				distance_matrix[j][i] = distance
				j+=1
	return distance_matrix
def average_histograms_in_time(times, starting_interval, ending_interval, frameskip, frame_collection):
	AverageMatrix 	= [[0 for x in range(2)]for y in range (500)]    	
	#print ("inside and outside")
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	clusterRowOne 		= []
	clusterRowTwo 		= []
	average_histogram 	= []
	average_histogram_x 	= []
	average_histogram_y 	= []
	compatibility_histogram = []
	ending_interval 	= int(ending_interval/frame_collection)
	corrected_interval 	= int (starting_interval/frame_collection)
	print (ending_interval)
	print (corrected_interval)
	print ((ending_interval-(corrected_interval)))
	for counter in range(corrected_interval,ending_interval+1,1):
		histogram 	= times[counter].chainHisto
		size 		= len(histogram)
		tempMatrix 	= AverageManager(histogram,size,AverageMatrix)
	print (counter)
	for x in range (0,size,1):
		if(AverageMatrix[x][0] != 0):		
			clusterRowOne.append(AverageMatrix[x][0])
			clusterRowTwo.append(AverageMatrix[x][1])	
	for x in range(0,len(clusterRowTwo),1):
		average_histogram_y.append(clusterRowOne[x]*clusterRowTwo[x])
		average_histogram_x.append(clusterRowOne[x])
	for x in range(0,len(clusterRowTwo),1):
		for y in range(0,int(clusterRowTwo[x]),1):
			average_histogram.append(clusterRowOne[x])
	print (AverageMatrix)	
	print("End")
	fig_3 	 		= plt.figure()
	dx 	 		= fig_3.add_subplot(111)
	min_hist 		= 0
	max_hist 		= 0
	hist_bin 		= 0
	min_hist 		= min(average_histogram)
	max_hist 		= max(average_histogram)
	hist_bin 		= max(10,int(max_hist - min_hist))	
	from scipy.stats import norm
	import matplotlib
	(mu,sigma) 		= norm.fit(average_histogram)
	n,bins,patches 		= dx.hist(average_histogram,normed=1,facecolor='green',alpha=0.75,bins='auto')
	y 			= mlab.normpdf( bins, mu, sigma)
	l 			= dx.plot(bins,y,'r--',linewidth=2)
	plt.title('Histogram: cluster size')
	plt.xlabel('Dimension')
	dx.set_xlim(xmin=0)
	plt.ylabel('Frequency')
	plt.grid(True)	
	plt.savefig('average_dist_'+str(ending_interval*frame_collection)+'.png')
	plt.close()
	fig_4 			= plt.figure()
	ex 			= fig_4.add_subplot(111)
	ordered_normalized_histogram_x = []
	ordered_normalized_histogram_y = []
	max_y 			= max(average_histogram_y)
	normalized_histogram_x  = average_histogram_x 
	normalized_histogram_y  = [ y / max_y for y in average_histogram_y]
	ordered_normalized_histogram_x = sorted(normalized_histogram_x)
	ordered_normalized_histogram_y = sortingHistograms(normalized_histogram_x,normalized_histogram_y)
	compatibility_histogram = reBuildingHistograms(ordered_normalized_histogram_x, ordered_normalized_histogram_y)	
	ex.plot(ordered_normalized_histogram_x,ordered_normalized_histogram_y, color='r', marker = 'o', alpha = 0.75)
	plt.title('Normalized Plot')
	plt.xlabel('Size')
	plt.ylabel('Frequency')
	plt.grid(True)	
	plt.savefig('normalized_plot_'+str(ending_interval*frame_collection)+'.png')
	plt.close()
	fig_5 = plt.figure()
	fx = fig_5.add_subplot(111)
	min_comp_hist 		= 0
	min_comp_hist 		= min(compatibility_histogram)
	max_comp_hist 		= max(compatibility_histogram)
	bin_comp_hist 		= max(10,int(max_comp_hist - min_comp_hist))
	fx.hist(compatibility_histogram,bins='auto',facecolor='green',alpha=0.75)
	plt.title('Normalized Hisrogram')
	plt.xlabel('Mimmo')
	plt.ylabel('Alessio Domenico Lavino')
	plt.grid(True)	
	plt.savefig('normalized_dist_'+str(ending_interval*frame_collection)+'.png')
	plt.close()

