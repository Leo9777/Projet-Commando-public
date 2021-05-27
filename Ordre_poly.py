#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 23 14:55:06 2021
@author: mwa
"""

import numpy as np
from tkinter import *
from PIL import Image, ImageTk
import matplotlib.pyplot as pp
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


def make_color_map ( wavelength ):
    """ Return a LinearSegmentedColormap going from transparent to full intensity
    for a wavelength given in nanometer .
    wavelength : float (in nm)
    """
    R, G, B, A = wavelength_to_rgb ( wavelength )
    colors = [(R,G,B,c) for c in np. linspace (0 ,1 ,100)]
    return matplotlib.colors.LinearSegmentedColormap . from_list ("mycmap ", colors )

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
    
    def __init__(self):
        self.a=0
        self.d=0
        self.f=0
        self.N=0
        self.theta=0
            
    def graphe_intensite(self, lamb,ordre_poly=0, I0=1):
        ordre = 2*lamb*self.f/self.a
        if ordre_poly!=0:
            ordre=ordre_poly
        X=np.linspace(-ordre,ordre,1000)
        In=I(I0,self.d,self.a,X,lamb,self.f,self.N,self.theta)
        return X,In
    
    def figure(self, lamb, I0=1):
        ordre = 2*lamb*self.f/self.a
        X,I=self.graphe_intensite(lamb, I0)
        y,x=np.meshgrid(I,X)
        comap=make_color_map(self.wv)
        extent=[-ordre,ordre,-ordre,ordre]
        return y, comap, extent


    def ordrepoly(self,lamb_tab):
        ordre=[]
        for lamb in lamb_tab:
            ordre.append(2*lamb*self.f/self.a)
        ordrepoly=np.mean(ordre)

        return ordrepoly
        

    
    def colorbar(self, wv):
        R, G, B, A = wavelength_to_rgb (wv)
        colors = [(0, 0, 0), (R, G, B)]
        cm = LinearSegmentedColormap.from_list("Custom", colors, N=100)
        cmappable = ScalarMappable(norm=Normalize(0,1), cmap=cm)
        return cmappable


class Images():
    def __init__(self):
        self.sequence = 0
        self.interferences = interference()
        self.lamb = []
        self.fig_graph = pp.figure(1)
        self.fig_fig = pp.figure(2)

    def update_var(self):
        self.interferences.a = float(scale_taille_fentes.get()*10**(-6))
        self.interferences.d = float(scale_distance_fentes.get()*10**(-6))
        self.interferences.f = float(scale_focale.get()/100)
        self.interferences.N = int(scale_nombre_fentes.get())
        self.interferences.theta = float(scale_angle.get())

    def trace_graphe(self, parent,lamb,ordre_poly=0, largeur=0):
        X,Ig = self.interferences.graphe_intensite(lamb,ordre_poly)
        R, G, B, A = wavelength_to_rgb(lamb*10**9)
        pp.figure(1)
        pp.plot(X, Ig,color=[R,G,B])
        pp.xlabel("Position sur l'écran (m)")
        pp.ylabel("Intensité lumineuse (cd)")
        pp.savefig("fig1.png")
        
        new_img = ImageTk.PhotoImage(Image.open("fig1.png").resize((432,288), Image.ANTIALIAS))
        img = Label(parent, image=new_img)
        img.image = new_img
        img.grid(column=0, row=0)
        
    def trace_figures(self, parent, lamb=0):
        y,comap,extent = self.interferences.figure(lamb)
        
        fig, ax = pp.subplots()
        ax.imshow(y,cmap=comap, extent=extent)
        ax.yaxis.set_visible(False)
        pp.xlabel("Position sur l'écran (m)")
        pp.colorbar(self.interferences.colorbar(lamb*10**9))#Affichage de la colorbar grâce à la classe
        pp. gca (). set_facecolor ("0")
        pp.savefig("fig2.png")
        
        new_img = ImageTk.PhotoImage(Image.open("fig2.png").resize((432,288), Image.ANTIALIAS))
        img = Label(parent, image=new_img)
        img.image = new_img
        img.grid(column=0, row=0)

class Interface_onde():
    def __init__(self):
        self.curseur_lambda = Scale(frame_HD, orient=HORIZONTAL, label="Longueur d'onde", length=130, from_=400, to=800)
        
    def param_onde(self):
        bouton.grid(column=0, row=10)
        bouton_supp.grid(column=0, row=11)
        self.curseur_lambda.grid(column=0, row=8)
        bouton_onde.configure(text="Valider", command=self.add_onde)
        

    def add_onde(self):
        if affichage.sequence == 0:
            lamb=int(self.curseur_lambda.get())*10**(-9)
            bouton_onde.configure(text='Nouvelle onde', command = self.param_onde)
            affichage.lamb.append(lamb)
            self.curseur_lambda.grid_forget()
            
            longueur_onde = 1
            
        else:
            lamb=int(self.curseur_lambda.get())*10**(-9)
            bouton_onde.configure(text='Nouvelle onde', command = self.param_onde)
            affichage.update_var()
            affichage.lamb.append(lamb)
            simulation()
            print(affichage.lamb)
            self.curseur_lambda.grid_forget()
        
  
def simulation():
    if len(affichage.lamb) == 0:
        print("Veuillez choisir")
    
    else:
        if affichage.sequence == 0:
            affichage.sequence = 1
            fenetre.geometry("580x728")
        
        affichage.update_var()
        ordre_poly=affichage.interferences.ordrepoly(affichage.lamb)
        for x in affichage.lamb:
            #affichage.trace_figures(frame_HG, x)
            affichage.trace_graphe(frame_BG, x,ordre_poly)
    
def effacer():
        affichage.lamb = []
        print(affichage.lamb)
        affichage.fig_graph.clf()

first_button_active = 1

fenetre = Tk()
fenetre.geometry("180x728")
fenetre.resizable(width=0, height=0)

frame_HG = Frame(fenetre)
frame_HD = Frame(fenetre)
frame_BG = Frame(fenetre)
frame_lambdas = Frame(frame_HD)
frame_HG.grid(column=0, row=0)
frame_HD.grid(column=1, row=0)
frame_BG.grid(column=0, row=1)

affichage = Images()
interface_onde = Interface_onde()


#scale_lambda = Scale(frame_HD, orient=HORIZONTAL, label="Longueur d'onde", length=130, from_=400, to=800)
#scale_lambda.grid(column=0, row=0)

scale_nombre_fentes = Scale(frame_HD, orient=HORIZONTAL, label="Nombre de fentes", length=130, from_=1, to=10)
scale_nombre_fentes.grid(column=0, row=1)

scale_taille_fentes = Scale(frame_HD, orient=HORIZONTAL, label="Taille des fentes", length=130, from_=1, to=1000)
scale_taille_fentes.grid(column=0, row=2)

scale_distance_fentes = Scale(frame_HD, orient=HORIZONTAL, label="Distance interfentes", length=130, from_=10, to=1000, resolution=10)
scale_distance_fentes.grid(column=0, row=3)

scale_angle = Scale(frame_HD, orient=HORIZONTAL, label="Angle", length=130)
scale_angle.grid(column=0, row=4)

scale_focale = Scale(frame_HD, orient=HORIZONTAL, label="Distance focale [cm]", length=130, from_=40, to=160)
scale_focale.grid(column=0, row=5)

txt_contraste = Label(frame_HD, text="Contraste: []")
txt_contraste.grid(column=0, row=6)

bouton_onde= Button(frame_HD, text='Ajouter une onde', command = interface_onde.param_onde)
bouton_onde.grid(column=0, row=7)

bouton = Button(frame_HD, text='Nouvelle simulation', command = simulation)
bouton.grid(column=0, row=8)

bouton_supp = Button(frame_HD, text='Effacer mémoire', command = effacer)
bouton.grid(column=0, row=9)



fenetre.mainloop()
