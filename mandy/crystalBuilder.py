import pandas as pd
import numpy as np
import os
import pathlib
from crystals import Crystal
from gemmi import cif

"""
Sub-routine designed to take a CIF file and produce the information needed for the diffraction calculation.
"""

def build(cifPath):    
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
    doc = cif.read_file(cifPath)
    block = doc.sole_block()
    # Read CIF file to facilitate extraction of lattice parameters
    crystalData = Crystal.from_cif(cifPath)
    # Storing calculated reciprocal lattice parameters 
    reciprocalLengths = np.array(crystalData.reciprocal.lattice_parameters[0:3])
    reciprocalAngles = np.array(crystalData.reciprocal.lattice_parameters[3:6])
    
    #create an empty list that will store all the data in this form: ['Ni', m1, m2, m3, 'F', m1, m2 m3]
    data = []
    for row in block.find_mmcif_category('_atom_site_moment.'):
        row = list(row)
        #print(row)
        for value in row:
            try:
                data.append(float(value))
            except ValueError:
                data.append(value)
        
    finalData = []
    i = 0
    while i < len(data):
        finalData.append(data[i:i+4])
        i+=4
    
    #Here we just create a pandas df, the code is standard and makes everything easy to acces later on
    #name the corresponding columns and also, make element to be the index
    mom_df = pd.DataFrame(finalData, columns=['Element','m1','m2','m3'])
    mom_df = mom_df.set_index('Element')
    # Recover the index labels by mom_df.index.tolist()
    print(mom_df)

    crystalArray = np.array(crystalData)
    
    #make it into and array that will make access and edits easier
    pos_df = pd.DataFrame(crystalArray,columns=['Element','x','y','z'])
    pos_df = pos_df.set_index('Element')
    print('----------------------------------------------------------------------------------------------')
    
    #mainly work with crystals package here to take the coordinates with their atomic number and just change 9 -> 'F', 28 -> 'Ni', etc
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
    
   
    pos_df['site_name'] = newIndex   #create a new column called newIndex
    pos_df = pos_df.set_index('site_name') #Now change its index to newIndex
    
    # Extract data for given sites as follows:
    # new_pos.loc[['Fe6h']]
    
    print(pos_df)  #I have created new dataframce with the corresponding names that I can use 
    #All that need to be done is add Fe2a and Fe6h to the magnetic data and you can combine them
    print('----------------------------------------------------------------------------------------------')
    
    return pos_df, mom_df, reciprocalLengths, reciprocalAngles



