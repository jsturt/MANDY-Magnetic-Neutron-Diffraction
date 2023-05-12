# import plotly.graph_objects as go
# import plotly.express as px
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cmath
import math
import os
import pathlib
from dataclasses import dataclass,field
from crystals import Crystal
from gemmi import cif


@dataclass
class site:
    """
    Encapsulates the information about a magnetic site. Defaults given for Fe2.
    
    Parameters
    ----------
    moment : list
        'classical' description of the moment direction and magnitude.
    l_s    : tuple
        l and s quantum numbers for the site, used in magnetic form factor calculation.
    ion_name   : string
        name of the desired ion, used in the magnetic form factor calculation.
    label  : string
        the internal label used to describe this site, useful for distinguishing between inequivalent sites of the same element e.g. fe2a, fe6h.
    Returns
    -------
    self : site
        dataclass representing a magnetic site.
    """

    moment: list
    l_s: tuple = (2,2)
    ion_name: str = 'Fe2'
    label: str = 'Fe2a'


class mandyCrystal:
    def __init__(self,cifpath:str,sites:list):
        self.cifPath = cifpath
        self.sites = sites
        #self.n = np.array([1,1,1])
    
    def build(self):    
        """
        Builds the crystal for use in calculations.
        
        Parameters
        ----------
        cifPath : STRING
            File path of the CIF file to build a crystal from.
    
        Returns
        -------
        pos_df : PANDAS DATAFRAME
            dataframe containing unit cell positions indexed by site name.
        mom_df : PANDAS DATAFRAME
            dataframe containing the moments associated with each type of site.
        len : NUMPY ARRAY (3,1)
            lattice parameters of the crystal.
        ang : NUMPY ARRAY (3,1)
            lattice angles of the crystal.
        """
        print('Building Crystal')
        print('----------------------------------------------------------------------------------------------')
       
        #Read the cif file and save it in a block, standard procedure when working with gemmi
        doc = cif.read_file(self.cifPath)
        block = doc.sole_block()
        # Read CIF file to facilitate extraction of lattice parameters
        crystalData = Crystal.from_cif(self.cifPath)
        crystalArray = np.array(crystalData)
        # Storing lattice parameters 
        self.a1, self.a2, self.a3 = crystalData.lattice_vectors
        self.b1, self.b2, self.b3 = crystalData.reciprocal_vectors
        
        reciprocalLengths = np.array(crystalData.reciprocal.lattice_parameters[0:3])
        reciprocalAngles = np.array(crystalData.reciprocal.lattice_parameters[3:6])
        

        # Produce list of lists containing the information of each magnetic site
        momentData = []
        for site in self.sites:
            print(site)
            momentData.append([site.label,site.moment[0],site.moment[1],site.moment[2]])

        ### Construct the pandas dataframes 
        # name the corresponding columns and also, make element to be the index
        mom_df = pd.DataFrame(momentData, columns=['Element','m1','m2','m3'])
        mom_df = mom_df.set_index('Element')
        # Recover the index labels by mom_df.index.tolist()
        print(mom_df)
        
        pos_df = pd.DataFrame(crystalArray,columns=['Element','x','y','z'])
        pos_df = pos_df.set_index('Element')
        print('----------------------------------------------------------------------------------------------')
        
        # Mainly work with crystals package here to take the coordinates with their atomic number and just change 9 -> 'F', 28 -> 'Ni', etc
        csvFilePath = pathlib.Path(__file__).parent.joinpath("atomic_data.csv")
        atomicData = pd.read_csv(csvFilePath)
        atomicData = atomicData.set_index('AtomicNumber')
        
        #Create this empty list that will contain all the symbols and later we will set this as an index for df2
        symbForPos = []
        for index, row2 in atomicData.iterrows():
            #print(index)
            for index2, rows in pos_df.iterrows():
                #print(j[0])
                if index == index2:
                    #create a new column that contains the symbol
                    symbForPos.append(row2['Symbol'])
                    
        #Print to see all the coresponding symbols and finally atach them to the df and make them an index, exacly as planned!      
        pos_df['symbols'] = symbForPos 
        pos_df = pos_df.set_index('symbols')
        print(pos_df)
        
        # Produce labels for the sites
        # Try to find a site name file in the current working directory (cwd)
    
        #Create an emty list where the new names will be stored
        newIndex = []
        
        cifFileName = self.cifPath[self.cifPath.rfind('/')+1:-4]
        filename = '{}SiteNames.dat'.format( cifFileName ) # File name based upon the CIF file.
        filepath = os.path.join(os.getcwd(),filename)   # Join the file name and cwd in an operating-system independent manner
        
        try:
            file = open(filepath,'r' )
            for line in file:
                newIndex.append(line.strip(f'\n'))
            file.close()
            
        except FileNotFoundError:    
            for ind, j in pos_df.iterrows():  # Iterate over indices and coordinates called j, ind = index
                coordinates = [j['x'], j['y'], j['z']]
            
                print("moment label for : "+str(coordinates))
                user_input = input()
                newIndex.append(user_input)
                #Save the coords in a list 
               
            file = open('{}SiteNames.dat'.format( cifFileName  ),'w' )
            [file.write(line + '\n') for line in newIndex]
            file.close()
        
        print('----------------------------------------------------------------------------------------------')
        
        # Convert fractional coordinates to cartesian coordinates

        pos_df.reset_index(inplace=True,drop=True)
        for i,row in enumerate(pos_df.itertuples()):
             cartesianArray = np.array((np.array(pos_df.loc[i].x).dot(self.a1)) + (np.array(pos_df.loc[i].y).dot(self.a2)) + (np.array(pos_df.loc[i].z).dot(self.a3)))
             pos_df.at[i,'x'] = cartesianArray[0]
             pos_df.at[i,'y'] = cartesianArray[1]
             pos_df.at[i,'z'] = cartesianArray[2]
            
        pos_df['site_name'] = newIndex   #create a new column called newIndex
        pos_df = pos_df.set_index('site_name') #Now change its index to newIndex  
        
        print(pos_df)  #I have created new dataframce with the corresponding names that I can use 
        #All that need to be done is add Fe2a and Fe6h to the magnetic data and you can combine them
        print('----------------------------------------------------------------------------------------------')
        
        self.indices = newIndex
        self.pos_df = pos_df
        self.mom_df = mom_df
        self.reciprocalLengths = reciprocalLengths
        self.reciprocalAngles = reciprocalAngles
        
   
    def createModulation(self, q, u, v, m1=1, m2=1, n=None):
        """
        default value of n set to ceiling function of 1 / q
        default direction is in c
        Parameters
        -------
        q
            wavevector of the modulation
        u,v
            unit vectors specifying the plane of the modulation
        m1,m2
            envelope controls of the modulation
        n
            numpy array of the unit cells

        Returns
        -------
        None.

       """
        # convert from fractional coords
        q = q[0]*self.b1 + q[1]*self.b2 + q[2]*self.b3

        # n is a [3x1] numpy array, want to set each component to 1/q
        tempN = []
        if(n is None):
            for i in range(3):
                tempN += [None]                         # Pad out the list since not using append here.
                if(q[i]>=0.01):
                    tempN[i] = math.ceil(1/q[i])     # Set value based upon the propagation vector, prevent fractional unit cells from being produced with ceiling function.
                    print('Using n{}={}'.format(['x','y','z'][i],tempN[i]))
                else:                                   # Prevent automatically generating 100s of unit cells
                    tempN[i] = 1
                    print('q < 0.01 ({0}-dir), setting n{0}=1 to prevent generation of 100s of unit cells. If this is a mistake please specify the kwarg "n" ([3,1] numpy array).'.format(*['x','y','z'][i]))
            self.n = np.array(tempN)
        else:
            self.n = n

        values = [ np.array(row3[1:4]) for row3 in self.pos_df.itertuples() ]    

        # Add unit cells in a,b,c directions
        values += [row + self.a1*(ix+1) for ix in range(self.n[0]-1) for row in values if self.n[0]>1] 
        values += [row + self.a2*(iy+1) for iy in range(self.n[1]-1) for row in values if self.n[1]>1] 
        values += [row + self.a3*(iz+1) for iz in range(self.n[2]-1) for row in values if self.n[2]>1] 
 
         
        reIndex = pd.Index(self.indices*np.prod(self.n), name = 'site_name')      # self.indices*n tiles the indices n times to allign with the supercell

        self.pos_df = pd.DataFrame(values, columns=['x','y','z'], index = reIndex)  # Updating positions dataframe to contain the supercell
        
        
        # Assigning the correct moments to each position in the supercell
        momentList = []
        for row in self.pos_df.itertuples():
            m = np.array( self.mom_df.loc[ row[0] ] )  # Recover moment from dataframe
            moment = (m) * (m1*u + 1j*m2*v) * cmath.exp( complex(0, -(np.dot(q,row[1:4]))) ) ### moment * modulation | removing 2pi factor
            # moment = np.linalg.norm(m) * (m1*u + 1j*m2*v) * cmath.exp( complex(0, -(np.dot(q,row[1:4]))) ) ### removing 2pi factor
            # moment = np.linalg.norm(m) * (m1*u + 1j*m2*v) * cmath.exp( complex(0, -2*cmath.pi * (np.dot(q,row[1:4]))) )
            momentList.append(moment.real)
         
        self.mom_df = pd.DataFrame(momentList, columns=['m1','m2','m3'],index=self.pos_df.index.copy())  # Updating positions dataframe to contain the supercell, copy over the indices from pos_df  




    #def plotCrystal(self, moment_scale = 0.6, aspect_ratio = dict(x=1, y=1, z=1)):
    #    scatter_points = px.scatter_3d(self.pos_df.reset_index(),x='x',y='y',z='z',color='site_name',size_max=2)
    #    scatter_points.update_layout(
    #                        scene_aspectmode='manual',
    #                        scene_aspectratio=aspect_ratio)
    #    scatter_points.add_trace(go.Cone(
    #        x = np.array(self.pos_df.x),
    #        y = np.array(self.pos_df.y),
    #        z = np.array(self.pos_df.z),
    #        u = np.array(self.mom_df.m1),
    #        v = np.array(self.mom_df.m2),
    #        w = np.array(self.mom_df.m3),
    #        sizemode = "absolute",
    #        sizeref = moment_scale,
    #        anchor = "tail",
    #        colorscale = "blues",
    #        hoverinfo = 'skip'))
    #    scatter_points.update_traces(showscale=False, selector=dict(type='cone'))
    #    scatter_points.show()

        
    def plotCrystal(self, moment_scale = 1):
            x = np.array(self.pos_df.x)
            y = np.array(self.pos_df.y)
            z = np.array(self.pos_df.z)
            u = np.array(self.mom_df.m1)
            v = np.array(self.mom_df.m2)
            w = np.array(self.mom_df.m3)

            fig = plt.figure()
            ax = plt.axes(projection='3d')
            ax.scatter(x,y,z)
            ax.quiver(x,y,z,u,v,w,length=moment_scale,pivot='middle')
            ax.set_aspect('equal')
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
            plt.show()
