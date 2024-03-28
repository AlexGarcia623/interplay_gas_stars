'''
This file is used to create Figure 5 of "Interplay of Stellar
and Gas-Phase Metallicities: Unveiling Insights for Stellar 
Feedback Modeling with Illustris, IllustrisTNG, and EAGLE"

Paper: https://ui.adsabs.harvard.edu/abs/2024MNRAS.tmp..787G/abstract

Code written by: Alex Garcia, 2023-24
'''
# Import from this library
from getAlpha import plot

SAVEDIR = '../Figures (pdf)/'

## see getAlpha.plot
plot('EAGLE',SAVEDIR + 'Figure5',redshift=2)