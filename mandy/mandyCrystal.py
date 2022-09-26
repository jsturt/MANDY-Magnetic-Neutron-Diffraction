import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import math
import os
import pathlib
from dataclasses import dataclass,field
from crystals import Crystal
from gemmi import cif


@dataclass
class site:
    """
    Encapsulates the information about a magnetic site.
    Default values are an example of 'Fe2'
    """
    moment: list
    l_s: tuple = (2,2)
    name: str = 'Fe2'
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
        print(crystalData.reciprocal_vectors)
        
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
        
        filename = '{}SiteNames.dat'.format( block.find_value('_chemical_formula_sum')) # File name based upon the CIF file.
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
               
            file = open('{}SiteNames.dat'.format( block.find_value('_chemical_formula_sum')),'w' )
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
        # Extract data for given sites as follows:
        # new_pos.loc[['Fe6h']]
        
        print(pos_df)  #I have created new dataframce with the corresponding names that I can use 
        #All that need to be done is add Fe2a and Fe6h to the magnetic data and you can combine them
        print('----------------------------------------------------------------------------------------------')
        
        self.indices = newIndex
        self.pos_df = pos_df
        self.mom_df = mom_df
        self.reciprocalLengths = reciprocalLengths
        self.reciprocalAngles = reciprocalAngles
        
    def createSDW(self, qSDW, n=None):
        """
        default value of n set to ceiling function of 1 / q
        default direction is in c

        Returns
        -------
        None.

       """
        
        # n is a [3x1] numpy array, want to set each component to 1/q
        tempN = []
        if(n is None):
            for i in range(3):
                tempN += [None]                         # Pad out the list since not using append here.
                if(qSDW[i]>=0.01):
                    tempN[i] = math.ceil(1/qSDW[i])     # Set value based upon the propagation vector, prevent fractional unit cells from being produced with ceiling function.
                    print('Using n{}={}'.format(['x','y','z'][i],tempN[i]))
                else:                                   # Prevent automatically generating 100s of unit cells
                    tempN[i] = 1
                    print('qSDW < 0.01 ({0}-dir), setting n{0}=1 to prevent generation of 100s of unit cells. If this is a mistake please specify the kwarg "n" ([3,1] numpy array).'.format(*['x','y','z'][i]))
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
            moment = np.array( self.mom_df.loc[ row[0] ] )#[2]  # Recover moment from dataframe
            momentList.append( moment * math.cos(2 * math.pi * np.array(row[1:4]).dot(qSDW) ) ) # Implement SDW, current_moment * cos(2pi * R.Q)
           
        self.mom_df = pd.DataFrame(momentList, columns=['m1','m2','m3'],index=self.pos_df.index.copy())  # Updating positions dataframe to contain the supercell, copy over the indices from pos_df  



    def plotCrystal(self, moment_scale = 0.6, aspect_ratio = dict(x=1, y=1, z=1)):
        scatter_points = px.scatter_3d(self.pos_df.reset_index(),x='x',y='y',z='z',color='site_name',size_max=2)
        scatter_points.update_layout(
                            scene_aspectmode='manual',
                            scene_aspectratio=aspect_ratio)
        scatter_points.add_trace(go.Cone(
            x = np.array(self.pos_df.x),
            y = np.array(self.pos_df.y),
            z = np.array(self.pos_df.z),
            u = np.array(self.mom_df.m1),
            v = np.array(self.mom_df.m2),
            w = np.array(self.mom_df.m3),
            sizemode = "absolute",
            sizeref = moment_scale,
            anchor = "tail",
            colorscale = "blues",
            hoverinfo = 'skip'))
        scatter_points.update_traces(showscale=False, selector=dict(type='cone'))
        scatter_points.show()

        