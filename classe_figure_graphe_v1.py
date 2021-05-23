#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


def make_color_map ( wavelength ):
    """ Return a LinearSegmentedColormap going from transparent to full intensity
    for a wavelength given in nanometer .
    wavelength : float (in nm)
    """
    R, G, B, A = wavelength_to_rgb ( wavelength )
    colors = [(R,G,B,c) for c in np. linspace (0 ,1 ,100)]
    return matplotlib.colors.LinearSegmentedColormap . from_list ("mycmap ", colors )


# In[5]:


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
            y.append(I0)
    return y

class interference:
    
    def __init__(self,a,d,lamb,f,N,theta):
        self.a=a
        self.d=d
        self.lamb=lamb
        self.f=f
        self.N=N
        self.theta=theta
        n=2
        I0=1
        self.I0=I0
        ordre=(n*self.lamb*self.f/self.a)
        self.ordre=ordre
        self.wv = self.lamb*10**9
        
    def graphe_intensite(self):
        X=np.linspace(-self.ordre,self.ordre,1000)
        In=I(self.I0,self.d,self.a,X,self.lamb,self.f,self.N,self.theta)
        return X,In
    
    def figure(self):
        X,I=self.graphe_intensite()
        y,x=np.meshgrid(I,X)
        comap=make_color_map(self.wv)
        extent=[-self.ordre,self.ordre,-self.ordre,self.ordre]
        return y, comap, extent
    
    def colorbar(self):
        R, G, B, A = wavelength_to_rgb (self.wv)
        colors = [(0, 0, 0), (R, G, B)]
        cm = LinearSegmentedColormap.from_list("Custom", colors, N=100)
        cmappable = ScalarMappable(norm=Normalize(0,1), cmap=cm)
        return cmappable
    
    
    
#Exemple de l'utilisation de la classe  
inte=interference(10*10**(-6),50*10**(-6),550*10**(-9),1,2,0 )

#Récupération des données nécessaires à l'affichage du graphique
X,Ig=inte.graphe_intensite()
#Récupération des données nécessaires à l'affichage de la figure
y,comap,extent=inte.figure()

#Affichage du graphique
plt.figure()
plt.xlabel("Position sur l'écran (m)")
plt.ylabel("Intensité lumineuse")
plt.plot(X,Ig)
plt.show()

#Affichage de la figure 
fig, ax = plt.subplots()
ax.imshow(y,cmap=comap, extent=extent)
ax.yaxis.set_visible(False)
plt.xlabel("Position sur l'écran (m)")

plt.colorbar(inte.colorbar())#Affichage de la colorbar grâce à la classe

plt. gca (). set_facecolor ("0")
plt. show()

