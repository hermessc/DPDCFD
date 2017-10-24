#GUI
#press function: Cancel = quit, Clusterize = Get values and run the algorithm, Create Video = make a video.
from appJar import gui
def press(btn):
	global ppo, box, conc, video_frequency, frameskip, flag, flag_ave, interval
	ppo = 9   		#number of PPO beads per chain
	box = 30		#Box dimension
	conc = 0.15		#concentration of surfactant in water
	switcher = {"1x":1, "2x":2, "5x":5, "10x":10, "25x":25} #declare switch cases
	flag = 0
	if btn =="Cancel":
		app.stop()     	#close the app
	elif btn =="Create Video": #make a video
		try:	 		#check if the cluster algorithm has been executed
			my_file = Path("/home/hermes/Desktop/cluster/example/timestep_0.png")
			if my_file.is_file():
				subprocess.run(["avconv", "-i", "timestep_%d.png", "-r", "25", "-c:v" ,"libx264", "-crf", "20" , "-pix_fmt", "yuv420p", "img.mov"]) #Create a video.mov
				subprocess.run(["convert", "img.mov","img.mp4"]) #convert into .mp4
				subprocess.run(["mv", "img.mp4", "video_result"])
				app.stop() #close the GUI
				print (quit()) #close python
			else:
				print ("Error! Run the cluster algorithm before") #Cluster algorithm must be executed before! (TODO: add a popup window)            			
				app.stop()	#close the gui
				print (quit()) 	#quit python
		except Exception as e:   		#remove exceptions
				if e.errno == errno.EACCES:  # 13
					print(quit()) #close python	
	elif btn == "Cheat!":
		ppo = 9
		box = 30
		conc = 0.05
		frameskip = 10
		flag = 0
		flag_ave = 1
		interval = 250000
		app.stop()	
	else: 				
		ppo = int(app.getEntry("PPO"))   #read PPO val
		box = int(app.getEntry("box"))	 #read box val
		conc = float(app.getEntry("conc"))	#read conc val	
		video_frequency = app.getSpinBox("Video Frequency") #take the value from video_frequency spinbox
		frameskip = int(switcher.get(video_frequency,"1x")) #convert string into int
		flag = int(app.getCheckBox("Histogram"))
		app.stop()			#close the GUI
def ave(btn):
	global interval, flag_ave
	flag_ave = 0
	if btn == "Save":
		interval = int(app.getEntry("num_int"))
		flag_ave = 1



app = gui("Hermes' Hut", "900x600")	#gui title
app.setBg("#FCF3CF")   			#gui background
app.setFont("18", "Comic Sans")		#gui font
app.startTabbedFrame("TabbedFrame")
app.setTabbedFrameTabExpand("TabbedFrame", expand=True)

			#TAB 1#
app.startTab("Data")
app.addLabel("title", "HermesClusterManager",0,0,2) 	#gui Main title
app.addLabel("PPO", "PPO: ",1,0)			#row 1, col 0, label
app.addEntry("PPO",1,1)					#row 1, col 1, editText:@+id/PPO
app.setEntryDefault("PPO", "PPO beads per chain")	#hint for PPO
app.addLabel("box","Box: ",2,0)				#row 2, col 0, label
app.addEntry("box",2,1)					#row 2, col 1, editText:@+id/box
app.setEntryDefault("box", "box size, eg: 20")		#hint for box
app.addLabel("conc","Concentration: ",3,0)		#row 3, col 0, label
app.addEntry("conc",3,1)				#row 3, col 1, editText:@+id/conc
app.setEntryDefault("conc", "Pluronic concentration, eg: 0.15")	#hint for conc
 #Change the frequency for videos
app.addCheckBox("Histogram",4,1)
app.setCheckBox("Histogram",False, True)
app.addButtons(["Clusterize","Cheat!"],press,5,0,2)
app.enableEnter(press) #assign "press enter to function "press"
app.setEntryFocus("PPO")	#the first entry will be PPO
app.setGuiPadding(4,4)		#Padding of the gui

app.stopTab()


			#TAB 2#
app.startTab("Post-Proc")
app.addLabel("num_int","#Intervals",1,0)
app.addEntry("num_int",1,1)
app.setEntryDefault("num_int","Average on last #intervals")	#hint for number of intervals
app.addButtons(["Save"],ave,1,2)
app.addLabelSpinBox("Video Frequency",["1x","2x","5x","10x","25x"],2,1)
app.setSpinBox("Video Frequency","1x",callFunction=True)		#set the spinbox and set 1x as default
app.addButtons(["ClusterMe","Create Video","Cancel"],press,3,0,3)	#define the 3 buttons (names), assign a function (press), row 4, col 0 , tab 2)
app.stopTab()
app.stopTabbedFrame()


app.go()			#launch
						#--------  END GUI --------#

