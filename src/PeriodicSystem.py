'''Updated March 7, 2009'''
from numpy import *
from numpy.linalg import inv
import math

from elements import *

class PeriodicSystem:
    '''A class for working with periodic structures'''
    
    def __init__(self, specietuple, coordinates, translation_vectors):
        '''Return a PeriodicSystem object
        Inputs:
         specielist - a Python tuple of species in the periodic system
         coordinates - a Numpy array (len(species) by 3) of the cartesian 
         coordinates of each specie in specielist (in Angstroms only at the 
         moment)
         translation_vectors - a numpy array of 3 translation vectors 
         defining the unit cell
         '''
        self.species = list(specietuple)
        self.coordinates = coordinates
        self.transvec = translation_vectors
        self.numAtoms = len(self.species)

    def expandxf(self, factor):
        '''multiplies the size of the PeriodicSystem by "factor" in 3 
        directions defined by self.transvec'''
        N = len(self.species)
        Nexpanded = N * (factor ** 3)
        speciesExpanded = self.species * (factor ** 3)
        coordinatesExpanded = reshape(zeros(Nexpanded * 3, float),
                                      (Nexpanded, 3))
        indexExpanded = 0
        for a in range(factor):
            for b in range(factor):
                for c in range(factor):
                    for atom in range(N):
                        coordinatesExpanded[indexExpanded] = \
                            self.coordinates[atom] + \
                            a * self.transvec[0] + \
                            b * self.transvec[1] + \
                            c * self.transvec[2]
                        indexExpanded += 1
        self.species = speciesExpanded
        assert len(speciesExpanded) == Nexpanded
        self.coordinates = coordinatesExpanded
        self.transvec = self.transvec * factor
        self.numAtoms = len(speciesExpanded)
        
    def calculatedN(self, dr, alpha, beta):
        '''Counts the number of beta atoms at distance r from each alpha atom 
        and calculates the average dN(r)(can be used for calculating 2-body 
        distribution function g(r)). Returns a list of values of r and 
        another list of values of dN corresponding to the list of r
        Inputs:
          dr - the thickness of the shell for calculating dN(r)
          alpha - type of atom
          beta - type of atom
          '''

        rmax = min(sqrt(dot(self.transvec[0], self.transve[0])), 
                   sqrt(dot(self.transvec[1], self.transve[1])), 
                   sqrt(dot(self.transvec[2], self.transve[2]))) / 2.0
        N = len(self.species)
        Nexpanded = N * 27
        speciesExpanded = self.species * 27
        coordinatesExpanded = reshape(zeros(Nexpanded * 3, float),
                                      (Nexpanded, 3))
        indexExpanded = 0
        for a in range(-1, 2):
            for b in range(-1, 2):
                for c in range(-1, 2):
                    for atom in range(N):
                        coordinatesExpanded[indexExpanded] = \
                            self.coordinates[atom] + \
                            a * self.transvec[0] + \
                            b * self.transvec[1] + \
                            c * self.transvec[2]
                        indexExpanded += 1
        rij = reshape(zeros(N * Nexpanded, float), (N, Nexpanded))
        for i in range(N):
            for j in range(Nexpanded):
                rij[i, j] = sqrt(sum((self.coordinates[i] - 
                                      coordinatesExpanded[j]) **2 ))
        dNi = reshape(zeros(N * int(rmax/dr), float), (N, (int(rmax / dr))))
        for i in range(N):
            if self.species[i] == alpha:
                for k in range(int(rmax/dr)):
                    r=(k+0.5)*dr
                    for j in range(Nexpanded):
                        if (speciesExpanded[j] == beta and 
                            rij[i, j] > r - dr / 2 and 
                            rij[i, j] <= r + dr / 2): 
                            dNi[i, k]+=1
        dN = 0.0
        for i in range(N):
            dN += dNi[i]
        return [(k + 0.5) * dr for k in range(int(rmax / dr))], dN

    def addAtom(self, index, specie, coordinate):
        '''Adds an atom of type specie at the specified coordinates and index
        Input
          index - the index number of the atom to be added
          specie - the specie of the atom (ex/Ni,Ar,...)
          coordinate - 3-dimensional numpy array of the cartesian coordinates 
          of the new atom
          '''
        newspecielist = self.species[0:index]
        newspecielist.append(specie)
        newspecielist.extend(self.species[index:])
        self.species = newspecielist
        self.coordinates = insert(self.coordinates, index, coordinate, axis=0)
        self.numAtoms = len(newspecielist)

    def deleteAtom(self, index):
        '''Deletes an atom of the specified index'''
        del self.species[index]
        self.coordinates = delete(self.coordinates,index,axis=0)
        self.numAtoms = len(self.species)

    def countAtom(self, alpha):
        '''Counts the number of alpha atoms and returns an integer'''
        Nalpha = 0
        for i in range(self.numAtoms):
            if self.species[i] == alpha: Nalpha += 1
        return Nalpha
    
    def make_neighbor_list(self, alpha, beta, distance1, distance2):
        '''Counts the number of beta atoms at distance1 or 2 +-0.0000001 
        from alpha atoms and returns the result in a dictionary 
        in the following form:
           {index_of_alpha0:{1:[list of index_of_beta at distance1],
           2[list of index_of_beta at distance2]},index_of_alpha1:{...}}'''
        species = self.species
        N = len(self.species)
        Nexpanded = N*27
        speciesExpanded = self.species*27
        index_of_original_cell = []
        coordinatesExpanded = reshape(zeros(Nexpanded * 3, float),
                                      (Nexpanded, 3))
        indexExpanded = 0
        neighbor_dict = {}
        delta = 1e-7
        for a in range(-1, 2):
            for b in range(-1, 2):
                for c in range(-1, 2):
                    for atom in range(N):
                       coordinatesExpanded[indexExpanded] = \
                           self.coordinates[atom] + \
                           a * self.transvec[0] + \
                           b * self.transvec[1] + \
                           c * self.transvec[2]
                       index_of_original_cell.append(atom)
                       indexExpanded += 1
        for i in range(N):
            if species[i] != alpha: 
                continue
            neighbor_dict[i] = {1:[],2:[]}
            for j in range(Nexpanded):
                if speciesExpanded[j] != beta: 
                    continue
                if sqrt(sum((self.coordinates[i] - 
                             coordinatesExpanded[j]) ** 2)) > \
                             distance1 - delta and \
                             sqrt(sum((self.coordinates[i] - 
                                       coordinatesExpanded[j]) ** 2)) < \
                                       distance1 + delta:
                    neighbor_dict[i][1].append(index_of_original_cell[j])
            for j in range(Nexpanded):
                if speciesExpanded[j] != beta: 
                    continue
                if sqrt(sum((self.coordinates[i] - 
                             coordinatesExpanded[j]) **2)) > \
                             distance2-delta and \
                             sqrt(sum((self.coordinates[i] - 
                                       coordinatesExpanded[j]) ** 2)) < \
                                       distance2 + delta:
                    neighbor_dict[i][2].append(index_of_original_cell[j])
        return neighbor_dict

    def toXYZ(self,xyzfile):
        '''Creates an xyz file to be read by JMOL,etc.
        Input
          xyzfile - file name of the xyz file'''
        output=open(xyzfile,'w')
        output.write(str(self.numAtoms)+'\n\n')
        for i in range(self.numAtoms):
            if len(self.species[i]) == 2: output.write(self.species[i] + '  ')
            if len(self.species[i]) == 1: output.write(self.species[i] + '   ')
            if len(self.species[i]) == 3: output.write(self.species[i] + ' ')
            output.write('%12.6f %12.6f %12.6f \n'% tuple(self.coordinates[i]))
        output.close()
        
    def toblockAtomicCoordinatesAndSpeciesInAng(self,species,datfile):
        '''Creates a text file containing %blockAtomicCoordinatesAndSpecies 
        in Angstroms to be used by SIESTA
        Input
          species - a list of the name of the species in the same order 
          as in the fdf input for SIESTA
          datfile - name of the output file'''
        output = open(datfile,'w')
        for i in range(self.numAtoms):
            for j in range(len(species)):
                if self.species[i] == species[j]:
                    output.write('%f \t %f \t %f \t' 
                                 %tuple(self.coordinates[i]) + str(j+1)+'\n')
                    break
        output.close()

    def toblockLatticeVectorsInAng(self,datfile):
        '''Creates a text file containing %blockLatticeVectors in
        Angstroms to be used by SIESTA
        Input
          datfile - name of the output file'''
        output = open(datfile,'w')
        for i in range(3):
            output.write('%f \t %f \t %f \t' %tuple(self.transvec[i]) +'\n')
        output.close()
    
    def makeAtomIndexDict(self):
        '''Create a dictionary of indexes with the name of atom species
        as keys'''
        index_dict = {}
        for i in self.species:
            index_dict[i] = []
        for i in range(self.numAtoms):
            index_dict[self.species[i]].append(i)
        return index_dict

    def sortBySpecies(self, orderlist):
        sorted_coor_tmp = zeros([self.numAtoms,3])
        species_list_tmp = []
        sorted_index = []
        for i in range(len(orderlist)): sorted_index.append([])
        for i in range(self.numAtoms):
            assert self.species[i] in orderlist       
            for j in range(len(orderlist)):
                if self.species[i] == orderlist[j]: sorted_index[j].extend([i])
        indexlist = []
        for i in range(len(sorted_index)):
            indexlist.extend(sorted_index[i])      
        for i in range(self.numAtoms):
            sorted_coor_tmp[i] = self.coordinates[indexlist[i]]
            species_list_tmp.extend([self.species[indexlist[i]]])
        self.coordinates = sorted_coor_tmp
        self.species = species_list_tmp

    def toPOSCAR(self, orderlist, system_name, d="c"):
        '''Creates POSCAR file for VASP'''
        self.sortBySpecies(orderlist)
        fi=open('POSCAR','w')
        fi.write(system_name+' \n'+ 
                 '   1.00 \n' +
                 ' %10.4f  %10.4f  %10.4f \n' % tuple(self.transvec[0]) +  
                 ' %10.4f  %10.4f  %10.4f \n' % tuple(self.transvec[1]) +  
                 ' %10.4f  %10.4f  %10.4f \n' % tuple(self.transvec[2])
                 )
        for i in orderlist:
            fi.write(' '+str(self.countAtom(i)))
        if d=="c":
            fi.write('\nCartesian\n')
            for i in range(self.numAtoms):
                fi.write(' %12.6f  %12.6f  %12.6f \n' % tuple(self.coordinates[i]))
            fi.close()
        elif d=="d":
            fi.write('\nDirect\n')
            invLmat = inv(self.transvec)
            for i in range(self.numAtoms):
                fi.write(' %12.11f  %12.11f  %12.11f \n' 
                         % tuple(dot(self.coordinates[i],invLmat)))
        else:
            print "Please specify c or d"
            fi.close()
    def toXSF(self, xsffi):
        '''Create XSF file of the PeriodicSystem
        Input
        xsffi -  the name of the xsf file
        '''
        file = open(xsffi, 'w')
        file.write('CRYSTAL\n' + 
                   'PRIMVEC\n' +
                   ' %12.6f  %12.6f  %12.6f \n' % tuple(self.transvec[0]) +  
                   ' %12.6f  %12.6f  %12.6f \n' % tuple(self.transvec[1]) +  
                   ' %12.6f  %12.6f  %12.6f \n' % tuple(self.transvec[2]) +
                   'PRIMCOORD\n' + 
                   str(self.numAtoms) + ' 1\n'
                   )
        for i in range(self.numAtoms):
            file.write(str(Name2Number(self.species[i])) + '    ' + 
                       '%12.6f  %12.6f  %12.6f \n' % tuple(self.coordinates[i])
                       )


#####################################################################
# Functions for creating PeriodicSystem from different file formats #
#####################################################################

def xyz2PeriodicSystem(xyzfile, translation_vectors):
    '''Creates a PeriodicSystem object from xyz file
    Input
      xyzfile - name of the input xyz file
      translation_vectors - translation vectors defining the unit cell
      '''
    file = open(xyzfile)
    file.readline()
    file.readline()
    lines = file.readlines()
    N = len(lines)
    file.close()
    specielist = []
    coordinates = reshape(zeros(N * 3, float), (N, 3))
    for i in range(N):
        lines[i] = lines[i].split()
        specielist.append(lines[i][0])
        coordinates[i] = [float(lines[i][j]) for j in [1,2,3]]
    species = tuple(specielist)
    return PeriodicSystem(species, coordinates, translation_vectors)

def POSCAR2PeriodicSystem(POSCAR, orderlist):
    '''Creates a PeriodicSystem object from VASP POSCAR file
    Input
      POSCAR - location of the input POSCAR file (can use a name other 
      than POSCAR
      orderlist - list of species in the same order as in the POTCAR file'''
    file = open(POSCAR)

    # Read the comment line.
    file.readline()  
    # Lattice constant, lattice vectors
    latcons = float(file.readline())  
    transvec = reshape(zeros(9, float), (3,3))
    for i in range(3):
        vec = file.readline().split()
        transvec[i] = [(float(vec[j]) * latcons) for j in [0,1,2]]

    # Count the number of atoms for each specie from POSCAR
    # Creation of species list
    nspecie_list = file.readline().split()  
    species = []
    assert len(nspecie_list) == len(orderlist)
    for i in range(len(nspecie_list)):
        nspecie_list[i] = int(nspecie_list[i])
    for i in range(len(orderlist)):
        species.extend([orderlist[i]]*nspecie_list[i])
    species = tuple(species)
    natoms = sum(nspecie_list)

    DirCar = file.readline().lstrip()
    if DirCar[0] == 'S' or DirCar[0] == 's':  # Selective dynamics?
        DirCar = file.readline().lstrip()  # Read direct or cartesian
    if DirCar[0] != 'D' and DirCar[0] != 'd' and \
            DirCar[0] != 'C' and DirCar[0] != 'c':
        raise IOError, 'Error in input in specification of cartesian' \
            'or direct'
    counter = 0
    coor = reshape(zeros(natoms * 3, float), (natoms, 3))
    while counter < natoms:
        line = file.readline().lstrip()
        line = line.split()
        for j in range(3):
            line[j] = float(line[j])
        if DirCar[0] == 'C' or DirCar[0] == 'c':
            coor[counter] = array(line[0:3]) * latcons
        else:
            coor[counter] = transvec[0] * line[0] + \
                transvec[1] * line[1] + \
                transvec[2] * line[2]
        counter += 1
    return PeriodicSystem(species, coor, transvec)

def XSF2PeriodicSystem(xsffile):
    '''Create a PeriodicSystem from XCrysDen XSF file'''
    file = open(xsffile)
