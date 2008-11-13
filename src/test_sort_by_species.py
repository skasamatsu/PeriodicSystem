import random
import numpy
import PeriodicSystem
import os
import copy
import math
import shutil

def getUniqueCombinations(alist, numb):
   for i in range(len(alist)):
      if numb==1:
         #templist=[alist[i]]
         yield [alist[i]]
      else:
         for combo in getUniqueCombinations(alist[i+1:], numb-1):
            #print templist
            templist=[alist[i]]
            templist.extend(combo)
            yield templist

#for combo in getUniqueCombinations([3,2,23,45,44],0):
#    print combo


m=3
n=3
az=5.12 # ZrO2 lattice constant


### Construct ZrO2 fcc (111) surfaces
# Distance between nearest neighbors on the same layer
dz = az/numpy.sqrt(2)
# Distance between zirconium layer and oxygen layer
dz_o = az*numpy.sqrt(3)/12.0

# Distance between oxygen layer and oxygen layer
do_o = 2.0*dz_o
# vector generators:
uz = numpy.array([1.0, 0.0, 0.0])*dz
vz = numpy.array([0.5, numpy.sqrt(3)/2.0, 0.0])*dz
duz = numpy.array([0.5, 1.0/2.0/numpy.sqrt(3), 0.0])*dz

ZA = numpy.reshape(numpy.zeros(3*m*m), (m*m,3))*1.0

for j in range(m):
    for i in range(m):
        ZA[i+m*j] = i*uz*1.0 + j*vz*1.0
     
ZA1 = ZA.copy()
ZB1 = ZA + 1.0*duz
ZB2 = ZA + 1.0*duz
ZC1 = ZA + 2.0*duz
OA1 = ZA.copy()
OA2 = ZA.copy()
OB1 = OA1 + 1.0*duz
OB2 = OA1 + 1.0*duz
OC1 = OA1 + 2.0*duz
OC2 = OA1 + 2.0*duz

# Set z coordinates for each layer
ZB1[:,2] += 0.0
OC1[:,2] += dz_o
OB1[:,2] += dz_o + do_o
ZC1[:,2] += dz_o + do_o + dz_o
OA1[:,2] += dz_o + do_o + dz_o + dz_o
OC2[:,2] += dz_o + do_o + dz_o + dz_o + do_o
ZA1[:,2] += dz_o + do_o + dz_o + dz_o + do_o + dz_o
OB2[:,2] += dz_o + do_o + dz_o + dz_o + do_o + dz_o + dz_o
OA2[:,2] += dz_o + do_o + dz_o + dz_o + do_o + dz_o + dz_o + do_o

# size of supercell lattice
zlength = dz_o + do_o + dz_o + dz_o + do_o + dz_o + dz_o + do_o + dz_o 

#zlength_0=round(zlength_0,1)
#getcontext().prec=5
#zlength_0 = Decimal(str(zlength_0))
xylength = az/numpy.sqrt(2)*m
# zlength_0=25.7

# Add atom number specification for Siesta
#zr = numpy.resize(numpy.array([1]),(m*m,1))
#o = numpy.resize(numpy.array([2]),(m*m,1))
#ZA1 = numpy.concatenate((ZA1,zr),1)
#ZB1 = numpy.concatenate((ZB1,zr),1)
#ZB2 = numpy.concatenate((ZB2,zr),1)
#ZC1 = numpy.concatenate((ZC1,zr),1)
#OA1 = numpy.concatenate((OA1,o),1)
#OA2 = numpy.concatenate((OA2,o),1)
#OB1 = numpy.concatenate((OB1,o),1)
#OB2 = numpy.concatenate((OB2,o),1)
#OC1 = numpy.concatenate((OC1,o),1)
#OC2 = numpy.concatenate((OC2,o),1)
# numpy.concatenate each layer to one numpy.array

scattering_positions = numpy.concatenate((ZB1, OC1, OB1, ZC1, OA1, OC2, ZA1, OB2, OA2))
lattice_vectors = numpy.array(
                              [[xylength,0.,0.],
                               [xylength*math.cos(math.pi/3.),xylength*math.sin(math.pi/3.),0.],
                               [0.,0.,zlength]]
                                                            )
#testvectors=numpy.array([[xylength,0.,0.],[xylength*math.cos(math.pi/3.),xylength*math.sin(math.pi/3.,0.),4],[3,4,5]])
species = tuple((['Zr']*9+['O']*18)*3)
ZrO2_111=PeriodicSystem.PeriodicSystem(species,scattering_positions,lattice_vectors)
ZrO2_111.sortBySpecies(['Zr','O'])
ZrO2_111.toXYZ('ZrO2.xyz')
ZrO2_111.toblockAtomicCoordinatesAndSpeciesInAng(['Zr','O'],'siesta.dat')