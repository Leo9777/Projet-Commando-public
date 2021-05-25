#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

def wavelength_to_rgb(wavelength, gamma=0.8):
    ''' taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an 
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range
    '''
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 750:
        A = 1.
    else:
        A=0.5
    if wavelength < 380:
        wavelength = 380.
    if wavelength >750:
        wavelength = 750.
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R,G,B,A)


# In[6]:


def make_color_map ( wavelength ):
    """ Return a LinearSegmentedColormap going from transparent to full intensity
    for a wavelength given in nanometer .
    wavelength : float (in nm)
    """
    R, G, B, A = wavelength_to_rgb ( wavelength )
    colors = [(R,G,B,c) for c in np. linspace (1 ,0 ,100)]
    return matplotlib.colors.LinearSegmentedColormap . from_list ("mycmap ", colors )


# In[46]:


from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

def sinc(x):
    return np.sinc(x/np.pi)

def I(I0,d,a,X,lamb,f,N,theta):
    c=2*np.pi/((f*lamb))
    y=[]
    for x in X:
        if x!=0:
            y.append(I0*(np.sin(c*N*(x-f*np.sin(theta))*d/2)/(N*np.sin(c*(x-f*np.sin(theta))*d/2)))**2*(sinc(c*(x-f*np.sin(theta))*a/2))**2)
        else:
            y.append(y[len(y-1)])
    return y

#Paramètres
I0=1  #Intensité initiale
d=50*10**(-6) #Distance entre les fentes
a=10*10**(-6) #Taille des fentes
lamb=570*10**(-9) #longueur d'onde
f=1 #focale
N=3 #Nombres de fentes
theta=0 #Angle d'incidence de la source

#Taille de l'affichage de la figure de diffraction (n=1 : largeur tâche centrale)
n=2
ordre=(n*lamb*f/a)#+2*f*np.sin(theta)


# Création et affichage du graphique de l'intensité en fonction de la position sur l'écran                
X=np.linspace(-ordre,ordre,1000)
I_graph=I(I0,d,a,X,lamb,f,N,theta)
plt.figure()
plt.xlabel("Position sur l'écran (m)")
plt.ylabel("Intensité lumineuse")
plt.plot(X,I_graph)
plt.show()


#Création et affichage de la figure d'interférence
y2,x2=np.meshgrid(I_graph,X)
wv = lamb*10**9# nm
R, G, B, A = wavelength_to_rgb (wv)
colors = [(R,G,B,a) for a in np. linspace (0 ,1 ,100) ]
comap=matplotlib.colors.LinearSegmentedColormap.from_list ("mycmap ", colors )
#Axes
fig, ax = plt.subplots()
ax.imshow(y2,cmap=comap, extent=[-ordre,ordre,-ordre,ordre])
ax.yaxis.set_visible(False)
plt.xlabel("Position sur l'écran (m)")


#Définition de la Colorbar
colors = [(0, 0, 0), (R, G, B)]
cm = LinearSegmentedColormap.from_list("Custom", colors, N=100)
cmappable = ScalarMappable(norm=Normalize(0,1), cmap=cm)
plt.colorbar(cmappable)

plt. gca (). set_facecolor ("0") # add a black background
plt. show()

