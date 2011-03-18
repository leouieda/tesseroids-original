"""
Module with functions that return rotation and reflection matrices.
Uses Numpy matrix objects.
All functions recieve angles in decimal degrees.
"""

import numpy as np
import math

# The ROTAION and REFLECTION matrices
# ##############################################################################
def R1(angle):
    """
    Returns in a matrix object the rotation matrix around the x axis.
    angle should be in degrees!
    """
    
    d2r = math.pi/180.0
    
    return np.matrix([[1.,0.,0.],[0.,np.cos(d2r*angle),np.sin(d2r*angle)],[0.,-np.sin(d2r*angle),np.cos(d2r*angle)]])

def R2(angle):
    """
    Returns in a matrix object the rotation matrix around the y axis.
    angle should be in degrees!
    """

    d2r = math.pi/180.0
    
    return np.matrix([[np.cos(d2r*angle),0.,-np.sin(d2r*angle)],[0.,1.,0.],[np.sin(d2r*angle),0.,np.cos(d2r*angle)]])

def R3(angle):
    """
    Returns in a matrix object the rotation matrix around the z axis.
    angle should be in degrees!
    """

    d2r = math.pi/180.0
    
    return np.matrix([[np.cos(d2r*angle),np.sin(d2r*angle),0.],[-np.sin(d2r*angle),np.cos(d2r*angle),0.],[0.,0.,1.]])

def P1():
    """
    Returns in a matrix object the reflection matrix of the x axis.
    """
    return np.matrix([[-1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    
def P2():
    """
    Returns in a matrix object the reflection matrix of the y axis.
    """
    return np.matrix([[1.,0.,0.],[0.,-1.,0.],[0.,0.,1.]])

def P3():
    """
    Returns in a matrix object the reflection matrix of the z axis.
    """
    return np.matrix([[1.,0.,0.],[0.,1.,0.],[0.,0.,-1.]])
# ##############################################################################

def L2G(x, y, z, lonO, latO, hO):
    """
    Converte x,y,z no sistema local para coordenadas globais esfericas
    """

    # Pega o raio medio da Terra
    R = 6378137.0

    d2r = math.pi/180.0
    r2d = 180.0/math.pi

    eL = np.matrix([x, y, z]).T

    # Converte a coordenada da origem para cartesiano
    eO = [math.cos(d2r*latO)*math.cos(d2r*lonO), \
          math.cos(d2r*latO)*math.sin(d2r*lonO), \
          math.sin(d2r*latO)]
    # Cria um vetor coluna
    eO = np.matrix(eO).T
    # Multiplica pelo raio da Terra para dar a dimensao correta
    eO = eO*(R + hO)

    # Faz a matriz de tranformacao
    l2g = R3(180. - lonO)*R2(90. - latO)*P2()

    # Alinha os dois sistemas e desloca a origem para eO
    eG = l2g*eL + eO

    # Converte pra esfericas
    eG = eG.T.tolist()[0]
    xG = eG[0]
    yG = eG[1]
    zG = eG[2]    
    r = math.sqrt(xG**2 + yG**2 + zG**2)
    lat = r2d*math.asin(zG/r)
    lon = r2d*math.atan2(yG, xG)
    h = r - R
    
    # Retorna uma lista com as coordenadas
    return [lon, lat, h]

    
def G2L(lonG, latG, hG, lonO, latO, hO):
    """
    Converte coordenadas globais esfericas para coordenadas cartesianas locais.
    (lonG,latG,hG) sao as coordendas globais que serao mudadas. hG eh a altitude
    acima do raio medio da Terra. (lonO,latO,hO) sao as coordendas da origem do
    sistema local.
    """

    # Pega o raio medio da Terra
    R = 6378137.0

    d2r = math.pi/180.0

    # Converte a coordenada da origem para cartesiano
    eO = [math.cos(d2r*latO)*math.cos(d2r*lonO), \
          math.cos(d2r*latO)*math.sin(d2r*lonO), \
          math.sin(d2r*latO)]
    # Cria um vetor coluna
    eO = np.matrix(eO).T
    # Multiplica pelo raio da Terra para dar a dimensao correta
    eO = eO*(R + hO)

    # Converte a coord global de esferico para cartesiano
    eG = [math.cos(d2r*latG)*math.cos(d2r*lonG), \
          math.cos(d2r*latG)*math.sin(d2r*lonG), \
          math.sin(d2r*latG)]
    # Cria um vetor coluna
    eG = np.matrix(eG).T
    # Multiplica pelo raio da Terra para dar a dimensao correta
    eG = eG*(R + hG)

    # Faz a matriz de tranformacao
    g2l = P2()*R2(latO - 90.)*R3(lonO - 180.)

    # Alinha os dois sistemas e desloca a origem para eO
    eL = g2l*(eG - eO)

    # Retorna uma lista com as coordenadas
    return eL.T.tolist()[0]

if __name__ == '__main__':
    eL = G2L(0, 45, 200, 0, 35, 0)
    print eL
    eG = L2G(eL[0], eL[1], eL[2], 0, 35, 0)
    print eG
