import mandyCrystal as md
import requestNIST as rN
import requests

crys = md.mandyCrystal(r'C:\Users\John\Documents\Uni\2022_Internship\Python_Code\MANDY-Magnetic-Neutron-Diffraction-off-axis-moments\NbFe2 Example\NbFe2_mp-568901_primitive.cif')
crys.build()

moments = [row[:] for row in crys.mom_df.itertuples()]
# rN.find_L_S('jeff')
# crys.createSDW(qSDW)
name = 'jef'

url = 'https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra={}%2B&submit=Retrieve+Data&units=1&format=3&order=0&ion_charge_out=on&el_name_out=on&level_out=on&e_out=0'.format('jef') # Formatted URL

try:
    x = requests.get(url) # Send request to NIST database
    if('Unrecognized token.' in x.text):
        raise NameError('Unrecognised Token "{}".'.format(name))
        
    x = x.text.splitlines()[1].split('\t') # Remove formatting
    x = [elem.replace('"','') for elem in x] # Remove quotation marks
    x = list(filter(None,x))[:-1] # Remove empty elements and ionisation energy collumn
    print('Retrieving term symbol: ' + str(x)) # Print finding : |Ion Charge|Element|Ground-state level|
    
    x = x[-1] # Leave only the term symbol
    termSymbol = x[:x.find('<')] # Strip out the J component
    
except requests.exceptions.RequestException as e:
    print('Could not reach NIST database :\n',e)
    termSymbol = '3P'
    print('Using L,S = 1')
    
except NameError as e:
    print(e)
    termSymbol = '3P'
    print('Using L,S = 1')
    