import numpy as np
import matplotlib.pyplot as plt
import os
import array
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import time


def fun(x, y):
  return x**2 + y


fig = plt.figure()
ax = fig.gca(projection='3d')
#ax = Axes3D(fig)
step = 0.1;


numFrames=40
for  i in range(1,numFrames):	
	fileParam = open('param','r');
	fileParam.readline()
	size_x=int(fileParam.readline()) 
	fileParam.readline()
	size_y=int(fileParam.readline()) 
	fileParam.readline()
	x_domains=int(fileParam.readline())
	fileParam.readline()
	y_domains=int(fileParam.readline()) 
	fileParam.readline()
	MaxStep=int(fileParam.readline()) 
	fileParam.readline()
	dt=float(fileParam.readline())  
	fileParam.readline()
	cnv_tol=float(fileParam.readline())   
	fileParam.close();
	#[X,Y]=meshgrid(0:sizex+1,0:sizey+1)
	x=np.linspace(0, size_x, size_x+1)
	y=np.linspace(0, size_y, size_y+1)
	X, Y = np.meshgrid(x, y)
	zfile='outputPar'+str(1)+'.dat'
	fin = open(zfile, 'rb')
	data=[]
	#print fin.readline()
	for j in range(0, size_y+1):
      		data.append(map(float, fin.readline().split()[:size_x+1]))	


	print len(data[0])
	fin.close()
	
	ax.plot_surface(X, Y, data)
	plt.show()
	time.sleep(0.5)

