from numpy import *
import math
class PeriodicSystem:
    '''A class for working with periodic structures'''
    
    def __init__(self,specietuple,coordinates,translation_vectors):
        '''Return a PeriodicSystem object
        Inputs:
         specielist - a Python tuple of species in the periodic system
         coordinates - a Numpy array (len(species) by 3) of the cartesian coordinates of each specie in specielist (in Angstroms only at the moment)
         translation_vectors - a numpy array of 3 translation vectors defining the unit cell
         '''
        self.species=list(specietuple)
        self.coordinates=coordinates
        self.translation_vectors=translation_vectors
        self.numAtoms=len(self.species)

    def expandx2(self):
        '''doubles the size of the PeriodicSystem in 3 directions defined by self.translation_vectors'''
        N=len(self.species)
        Nexpanded=N*8
        speciesExpanded=self.species*8
        coordinatesExpanded=reshape(zeros(Nexpanded*3,float),(Nexpanded,3))
        indexExpanded=0
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    for atom in range(N):
                       coordinatesExpanded[indexExpanded]=self.coordinates[atom]+a*self.translation_vectors[0]+b*self.translation_vectors[1]+c*self.translation_vectors[2]
                       indexExpanded+=1
        self.species=speciesExpanded
        self.coordinates=coordinatesExpanded
        self.translation_vectors=self.translation_vectors*3
        self.numAtoms=len(speciesExpanded)
        
    def expandx3(self):
        '''same as expandx2 but triples the size of the system'''
        N=len(self.species)
        Nexpanded=N*27
        speciesExpanded=self.species*27
        coordinatesExpanded=reshape(zeros(Nexpanded*3,float),(Nexpanded,3))
        indexExpanded=0
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    for atom in range(N):
                       coordinatesExpanded[indexExpanded]=self.coordinates[atom]+a*self.translation_vectors[0]+b*self.translation_vectors[1]+c*self.translation_vectors[2]
                       indexExpanded+=1
        self.species=speciesExpanded
        self.coordinates=coordinatesExpanded
        self.translation_vectors=self.translation_vectors*3
        self.numAtoms=len(speciesExpanded)
        
    def calculatedN(self,dr,alpha,beta):
        '''Counts the number of beta atoms at distance r from each alpha atom and calculates the average dN(r)(can be used for calculating 2-body distribution function g(r)). Returns a list of values of r and another list of values of dN corresponding to the list of r
        Inputs:
          dr - the thickness of the shell for calculating dN(r)
          alpha - type of atom
          beta - type of atom
          '''
        veca=self.translation_vectors[0]
        vecb=self.translation_vectors[1]
        vecc=self.translation_vectors[2]
        rmax=min(sqrt(dot(veca,veca)),sqrt(dot(vecb,vecb)),sqrt(dot(vecc,vecc)))/2.0
        N=len(self.species)
        Nexpanded=N*27
        speciesExpanded=self.species*27
        coordinatesExpanded=reshape(zeros(Nexpanded*3,float),(Nexpanded,3))
        indexExpanded=0
        for a in range(-1,2):
            for b in range(-1,2):
                for c in range(-1,2):
                    for atom in range(N):
                       coordinatesExpanded[indexExpanded]=self.coordinates[atom]+a*self.translation_vectors[0]+b*self.translation_vectors[1]+c*self.translation_vectors[2]
                       indexExpanded+=1
        rij=reshape(zeros(N*Nexpanded,float),(N,Nexpanded))
        for i in range(N):
            for j in range(Nexpanded):
                rij[i,j]=sqrt(sum((self.coordinates[i]-coordinatesExpanded[j])**2))
        dNi=reshape(zeros(N*int(rmax/dr),float),(N,(int(rmax/dr))))
        for i in range(N):
            if self.species[i]==alpha:
                for k in range(int(rmax/dr)):
                    r=(k+0.5)*dr
                    for j in range(Nexpanded):
                        if speciesExpanded[j]==beta and rij[i,j]>r-dr/2 and rij[i,j]<=r+dr/2: dNi[i,k]+=1
        dN=0.0
        for i in range(N):
            dN+=dNi[i]
        return [(k+0.5)*dr for k in range(int(rmax/dr))],dN

    def addAtom(self,index,specie,coordinate):
        '''Adds an atom of type specie at the specified coordinates and index
        Input
          index - the index number of the atom to be added
          specie - the specie of the atom (ex/Ni,Ar,...)
          coordinate - 3-dimensional numpy array of the cartesian coordinates of the new atom
          '''
        newspecielist=self.species[0:index]
        newspecielist.append(specie)
        newspecielist.extend(self.species[index:])
        self.species=newspecielist
        self.coordinates=insert(self.coordinates,index,coordinate,axis=0)
        self.numAtoms=len(newspecielist)

    def deleteAtom(self,index):
        '''Deletes an atom of the specified index'''
        del self.species[index]
        self.coordinates=delete(self.coordinates,index,axis=0)
        self.numAtoms=len(self.species)

    def countAtom(self,alpha):
        '''Counts the number of alpha atoms and returns an integer'''
        Nalpha=0
        for i in range(self.numAtoms):
            if self.species[i]==alpha: Nalpha+=1
        return Nalpha
    
    def make_neighbor_list(self,alpha,beta,distance1,distance2):
        '''Counts the number of beta atoms at distance1 or 2 +-0.0000001 from alpha atoms and returns the result in a dictionary in the following form:
           {index_of_alpha0:{1:[list of index_of_beta at distance1],2[list of index_of_beta at distance2]},index_of_alpha1:{...}}'''
        species=self.species
        veca=self.translation_vectors[0]
        vecb=self.translation_vectors[1]
        vecc=self.translation_vectors[2]
        N=len(self.species)
        Nexpanded=N*27
        speciesExpanded=self.species*27
        index_of_original_cell=[]
        coordinatesExpanded=reshape(zeros(Nexpanded*3,float),(Nexpanded,3))
        indexExpanded=0
        neighbor_dict={}
        for a in range(-1,2):
            for b in range(-1,2):
                for c in range(-1,2):
                    for atom in range(N):
                       coordinatesExpanded[indexExpanded]=self.coordinates[atom]+a*veca+b*vecb+c*vecc
                       index_of_original_cell.append(atom)
                       indexExpanded+=1
        #print alpha,' - ',beta,' neighboring index'
        for i in range(N):
            if species[i]!=alpha: continue
            neighbor_dict[i]={1:[],2:[]}
            #print '1st neighbor'
            for j in range(Nexpanded):
                if speciesExpanded[j]!=beta: continue
                if sqrt(sum((self.coordinates[i]-coordinatesExpanded[j])**2))>distance1-0.0000001 and sqrt(sum((self.coordinates[i]-coordinatesExpanded[j])**2))<distance1+0.0000001:
                    #print index_of_original_cell[j]
                    neighbor_dict[i][1].append(index_of_original_cell[j])
            #print '2nd neighbor'
            for j in range(Nexpanded):
                if speciesExpanded[j]!=beta: continue
                if sqrt(sum((self.coordinates[i]-coordinatesExpanded[j])**2))>distance2-0.0000001 and sqrt(sum((self.coordinates[i]-coordinatesExpanded[j])**2))<distance2+0.0000001:
                    #print index_of_original_cell[j]
                    neighbor_dict[i][2].append(index_of_original_cell[j])
        return neighbor_dict
    def toXYZ(self,xyzfile):
        '''Creates an xyz file to be read by JMOL,etc.
        Input
          xyzfile - file name of the xyz file'''
        output=open(xyzfile,'w')
        output.write(str(self.numAtoms)+'\n\n')
        for i in range(self.numAtoms):
            if len(self.species[i])==2:output.write(self.species[i]+'  ')
            if len(self.species[i])==1:output.write(self.species[i]+'   ')
            if len(self.species[i])==3:output.write(self.species[i]+' ')
            for j in range(3): output.write('%12.6f'% self.coordinates[i,j])
            output.write('\n')
        output.close()
        
    def toblockAtomicCoordinatesAndSpeciesInAng(self,species,datfile):
        '''Creates a text file containing %blockAtomicCoordinatesAndSpecies in Angstroms to be used by SIESTA
        Input
          species - a list of the name of the species in the same order as in the fdf input for SIESTA
          datfile - name of the output file'''
        output=open(datfile,'w')
        for i in range(self.numAtoms):
            for j in range(len(species)):
                if self.species[i]==species[j]:
                    output.write('%f \t %f \t %f \t' %tuple(self.coordinates[i]) + str(j+1)+'\n')
                    break
        output.close()
    

def xyz2PeriodicSystem(xyzfile,translation_vectors):
    '''Creates a PeriodicSystem object from xyz file
    Input
      xyzfile - name of the input xyz file
      translation_vectors - translation vectors defining the unit cell
      '''
    file=open(xyzfile)
    file.readline()
    file.readline()
    lines=file.readlines()
    N=len(lines)
    file.close()
    specielist=[]
    coordinates=reshape(zeros(N*3,float),(N,3))
    for i in range(N):
        lines[i]=lines[i].split()
        specielist.append(lines[i][0])
        coordinates[i]=[float(lines[i][j]) for j in [1,2,3]]
    return PeriodicSystem(specielist,coordinates,translation_vectors)

        
        
    
