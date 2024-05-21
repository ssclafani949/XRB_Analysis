import numpy as np
import astropy.io.fits as fits

source_list_1  = np.genfromtxt("./swift_select.txt",dtype=str,delimiter=',')  #swift sources
names_1        = [source_list_1[i][0].replace(' ','') for i in range(len(source_list_1))]
dec_1          = source_list_1[:,-2]
ra_1           = source_list_1[:,-1]

names_1 = np.array(names_1)
dec_1 = np.array(dec_1)
ra_1 = np.array(ra_1)

missing = ['4U 0352+309', '2S 1417-624', 'RX J0440.9+4431', '4U 2206+543', '2S 1553-542', '2S 1845-024', 'EXO 1722-363', '4U 1223-624', 'IGR J16320-4751', 'OAO 1657-415', '4U 1119-603', '4U 1907+09', 'GRO J1750-27', '4U 1258-61', '4U 1538-52', '4U 1700-37', '2S 0114+650', '1H 2138+579', '4U 0900-40', 'XTE J1739-302', '4U 2030+40','GS 1354-64', '4U 1659-487', '4U 1630-47', 'H1743-322', 'GRO J1719-24', 'Swift J1753.5-0127', 'GS 2023+338', '3A 1516-569', '4U 1608-52', '4U 1642-45', '4U 1735-444', '4U 1705-44', '3A 1702-363', '4U 1728-34', 'MXB 1730-335', '4U1820-30', 'KS 1741-293', '4U 1758-25', 'HETE J1900.1-2455', 'AX J1824.5-2451', '3A 1728-247', 'GS 1826-238', '4U 1811-17', 'H 1617-155', '4U 1813-14', '4U1908+005', '4U 0614+091', '3A 1954+319', '2A 1655+353', '4U 2142+38','CXOGCJ174540.0-290031', 'GRS 1741.9-2853', 'SAX J1748.9-2021', '3A1909+048', '4U1956+35']

hdul = fits.open('/home/ssclafani/Downloads/BAT_catalog.fits')
new_cat = hdul[1]
new_cat_ar = np.array(new_cat.data)
new_decs = np.array(new_cat_ar['DEC_OBJ'], dtype = 'float')
new_ras = np.array(new_cat_ar['RA_OBJ'], dtype = 'float')
new_names = np.array(new_cat_ar['NAME'], dtype = 'str')

replacement_names = ['GS 1843-02', 'IGR J18245-2452'] #missing from .fit file need to add manually 
print(names_1)

for m in missing:
    print('------------------------------')
    print('Missing Name: {}'.format(m))
    m = m.replace(' ','') 
    mask = (names_1 == m) 
    ra = ra_1[mask]
    dec = dec_1[mask]
    print(ra, dec)
    idx = np.where(np.isclose(new_decs ,float(dec[0]), atol=.1))
    #print(idx, len(idx[0]))
    if len(idx[0]) == 1:
        if np.isclose(new_ras[idx], float(ra[0]), atol=.1):
            print('Match Found:')
            print(new_ras[idx], new_decs[idx])
            print(new_names[idx])
            replacement_names.append(new_names[idx][0].strip())
    elif len(idx[0]) == 0:
        print(m)
        print('NO MATCH!!!!!!')
    else:
        found = False
        for i, id in enumerate(idx[0]):
            if np.isclose(new_ras[id], float(ra[0]), atol=.1):
                print('Match Found:')
                print(new_ras[id], new_decs[id])
                print(new_names[id])
                found = True
                replacement_names.append(str(new_names[id]).strip())
        if found == False:
            print(m)
            print('Did not find counterpart')

print(len(missing))
print(len(replacement_names))
print(replacement_names)
np.save('replacement_names.npy', list(zip(replacement_names, missing)))
