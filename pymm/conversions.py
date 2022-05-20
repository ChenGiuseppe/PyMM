''' Module that collects some often occurring units convertions (in particular
    to and from a.u.)'''

# Dictionaries to help converting between different formats
Z2mass = {'1': 1.00784,
          '6': 12.0107,
          '7': 14.0067,
          '8': 15.999,
          '16': 32.065}

mass2symbol = {1: 'H',
               12: 'C',
               14: 'N',
               16: 'O',
               32: 'S'}

symbol2mass = {v: k for k, v in mass2symbol.items()}

# Units convertion
au2eV = 27.2113961
Debye2au = 0.393456
Bohr2Ang = 0.529177249

#
k_B = 8.617333*10**(-5) # eV/K
R = 8.314 # J/(K*mol)

if __name__ == '__main__':
    pass
