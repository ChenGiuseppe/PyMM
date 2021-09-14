''' Module that collects some often occurring units convertions (in particular
    to and from a.u.)'''

# Dictionaries to help converting between different formats
Z2mass = {'1': 1.0, '6': 12.0, '7': 14.0, '8': 16.0, '16': 32.0}
mass2symbol = {1.0: 'H', 12.0: 'C', 14.0: 'N', 16.0: 'O', 32.0: 'S'}
symbol2mass = {v: k for k, v in mass2symbol.items()}

# Units convertion
au2eV = 27.2113961
Debye2au = 0.393456
Bohr2Ang = 0.529177249

if __name__ == '__main__':
    print(mass2symbol[32.0])
