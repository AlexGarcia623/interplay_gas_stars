'''Module used to load EAGLE data from SQL query

'''
# EAGLE Subhalo Query
import numpy as np
import sys
import h5py
from os import path, mkdir

WHICH_SIM = "EAGLE".upper()

# https://github.com/kyleaoman/eagleSqlTools
EAGLE_SQL_TOOLS = '/home/alexgarcia/github/eagleSqlTools' 

#########################################################################
################# YOU WILL PROBABLY HAVE TO CHANGE THIS #################
savedir = '../Data/%s/' %WHICH_SIM
#########################################################################

SAVE_DATA = True

if not (path.exists(savedir)):
    mkdir( savedir )
    print('New directory %s added' %savedir)

sys.path.insert(1,EAGLE_SQL_TOOLS)
import eagleSqlTools as sql

# con = sql.connect( ENTER_YOUR_USER_NAME, password=ENTER_YOUR_PASSWORD )
if !con:
    raise Exception('You need to make your own EAGLE username/password!: https://virgodb.dur.ac.uk/')

if (WHICH_SIM == 'EAGLE'):
    SIM = 'RefL0100n1504'

snaps = [28,19,15,12,10,8,6,5,4,3,2]

for snap in snaps:
    
    currentDir = savedir + 'snap%s/' %snap
    
    if not (path.exists(currentDir)):
        mkdir(currentDir)
        print('New directory %s added' %currentDir)
    
    myQuery = '''SELECT \
        SF_Metallicity as Zgas,\
        Stars_Metallicity as Zstar,\
        StarFormationRate as SFR,\
        MassType_Star as Stellar_Mass,\
        MassType_Gas as Gas_Mass,\
        HalfMassRad_Gas as R_gas,\
        HalfMassRad_Star as R_star\
    FROM \
        RefL0100N1504_SubHalo as SH\
    WHERE \
        SnapNum = %s \
        and SH.SubGroupNumber = 0 
        and SH.StarFormationRate > 0.0\
        and SH.MassType_Star > 1E8\
        and SH.SubGroupNumber = 0''' %(snap) 
    
    print('Starting Query... snapshot %s' %snap)
    myData = sql.execute_query(con, myQuery)
        
    if SAVE_DATA:
        print('\tSaving Data')
        np.save(currentDir+'Zgas'        , np.array(myData['Zgas'][:]) )
        np.save(currentDir+'Zstar'       , np.array(myData['Zstar'][:]) )
        np.save(currentDir+'SFR'         , np.array(myData['SFR'][:]) )
        np.save(currentDir+'Stellar_Mass', np.array(myData['Stellar_Mass'][:]) )
        np.save(currentDir+'Gas_Mass'    , np.array(myData['Gas_Mass'][:]))
        np.save(currentDir+'R_gas'       , np.array(myData['R_gas'][:]))
        np.save(currentDir+'R_star'      , np.array(myData['R_star'][:]))
    
    print('Query Complete')
    print('\n')
    