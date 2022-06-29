import requests

term_letters = 'SPDFGHIKLMNOQRTUVWXYZ'

def find_L_S(name):
    """
    Parameters
    ----------
    name : STRING
        Name of the Ion, e.g. 'Fe2' or 'Nb0'.

    Returns
    -------
    L : INT
        Total angular momentum quantum number.
    S : FLOAT
        Total spin quantum number.
    """

    url = 'https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra={}%2B&submit=Retrieve+Data&units=1&format=3&order=0&ion_charge_out=on&el_name_out=on&level_out=on&e_out=0'.format(name) # Formatted URL
    
    try:
        x = requests.get(url) # Send request to NIST database
        
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


    # Interpreting the term symbol
    S = (float(termSymbol[0]) - 1) / 2
    L = term_letters.index(termSymbol[1])
    return L,S

