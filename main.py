import matplotlib.pyplot as plt
import scipy
from scipy import special
import numpy as np
import harmo_sphe as hs
import math

"""
Paramètres
"""
m=0
l=2
h=0.1
x=np.arange(-1,1+h,h)

"""
Calcul des harmoniques sphériques 
"""
#Ma fonction 
n=np.size(x)
y1=np.zeros(n)
y2=np.zeros(n)
i=0
for k in x:
	Y=hs.legendre(l,k)
	Y2=scipy.special.lpmn(m,l,k)
	y1[i]=Y[l,m]
	y2[i]=Y2[0][m,l]
	i+=1
#Python 
#y2=scipy.special.sph_harm(m, l, 0, x)
"""
Plot 
"""
#Mon truc
plt.plot(x, y1, color='red', linewidth = 3)
#Python 
plt.plot(x, y2, color='blue', linestyle='dashed', linewidth = 3) 
#Légendes et affichages
plt.title('Harmoniques sphériques')
plt.grid()
plt.show()
