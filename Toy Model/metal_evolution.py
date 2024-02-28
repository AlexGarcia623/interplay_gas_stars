import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from sklearn.linear_model import LinearRegression

class metals():
        
    def __init__(self, Zstar0, Zgas0, sSFR0, dt):
        # Initial values
        self.Zstar0  = Zstar0
        self.Zgas0   = Zgas0
        # Current values
        self.Zstar   = self.Zstar0
        self.Zgas    = self.Zgas0
        # sSFR
        self.sSFR    = sSFR0
        # time counters
        self.t       = 0
        self.T       = -dt
        # initialize current target
        self.targetZ = None
        # set timestep size
        self.deltaZ  = None
    
        # Constants
        self.dt = dt
                
        # Helpers
        self.nSteps = 1e8
        
        # Initialize these arrays with arbitrarily large number of elements
        self.Zstars  = np.zeros(int(self.nSteps))
        self.Zgasses = np.zeros(int(self.nSteps))
        self.ts      = np.zeros(int(self.nSteps))
        
        # Set initial value of arrays to initial values
        self.Zstars [0] = self.Zstar
        self.Zgasses[0] = self.Zgas
        self.ts     [0] = self.t     
        
    def doIntegration( self ):
        while (self.t < self.T):            
            if (self.index == self.nSteps):
                print('Index larger than nSteps!!!')
                break
                
            self.Zstar += self.sSFR*(self.Zgas - self.Zstar)*self.dt
            
            self.Zgas += self.deltaZ
            
            self.t += self.dt
            
            self.Zstars [self.index] = self.Zstar
            self.Zgasses[self.index] = self.Zgas
            self.ts     [self.index] = self.t 
            
            self.index += 1
        
    def run( self, ZgasPertb, coherenceTime ):
                
        self.ZgasPertb = ZgasPertb
        
        self.index = 0
        for index, newZgasTarget in enumerate(self.ZgasPertb):
            self.T += coherenceTime[index]
            
            self.targetZ = newZgasTarget
            
            if (self.T - self.t)/self.dt < 1:
                # If coherence time is less than dt modify slightly
                # This is a bug fix, this happens ~1-5/10,000 times
                self.deltaZ = self.targetZ - self.Zgas
            else:
                self.deltaZ = self.dt*(self.targetZ - self.Zgas) / (self.T - self.t)
                
            self.doIntegration()
                        
        mask = self.ts > 0
        
        # Remove last step (is an extra 0)
        self.Zstars  = self.Zstars [mask]
        self.Zgasses = self.Zgasses[mask]
        self.ts      = self.ts     [mask]
        
    def plot(self):
        
        plt.hist2d( self.Zgasses, self.Zstars, bins=(100,100), norm=LogNorm() )

        plt.axhline(0,color='k')
        plt.axvline(0,color='k')

        x = self.Zgasses
        y = self.Zstars
        
        lm = LinearRegression(fit_intercept = False)
        
        lm.fit(x.reshape(-1, 1), y)
        X = np.arange( np.min(x), np.max(x), 0.1)
        preds = lm.predict( X.reshape(-1,1) )

        a,b = lm.coef_, lm.intercept_
        plt.plot(X, preds, color='k',lw=2.0)
        plt.plot(X, preds, color='w',label=f"Slope: {a}")
        plt.legend()
        
        plt.show()
        
    def time_series(self):
        
        plt.plot( self.ts, self.Zgasses )
        plt.plot( self.ts, self.Zstars )
        
        plt.axhline(0,color='k')
        
        plt.show()
        
        
class metals_MZR():
        
    def __init__(self, Zstar0, Zgas0, sSFR0, dt, MZR_g_init, MZR_s_init, dZsdt, dZgdt):
        # Initial values
        self.Zstar0 = Zstar0
        self.Zgas0 = Zgas0
        self.MZR_g_init = MZR_g_init
        self.MZR_s_init = MZR_s_init
        self.dZsdt = dZsdt
        self.dZgdt = dZgdt
        # Current values
        self.Zstar = self.Zstar0
        self.Zgas = self.Zgas0
        # sSFR
        self.sSFR = sSFR0
        # time counters
        self.t = 0
        self.T = 0
        # initialize current target
        self.targetZ = None
        # set timestep size
        self.deltaZ  = None
    
        # Constants
        self.dt = dt
        
        # Helpers
        self.nSteps = 1e8
        
        # Initialize these arrays with arbitrarily large number of elements
        self.Zstars  = np.zeros(int(self.nSteps))
        self.Zgasses = np.zeros(int(self.nSteps))
        self.ts      = np.zeros(int(self.nSteps))
        # Set initial value of arrays to initial values
        self.Zstars [0] = self.Zstar
        self.Zgasses[0] = self.Zgas
        self.ts     [0] = self.t
        
    def doIntegration( self ):
        while (self.t < self.T):            
            if (self.index == self.nSteps):
                print('Index larger than nSteps!!!')
                break
                
            newT = self.t + self.dt
            
            MZR_diff_term = self.sSFR * (self.MZR_s_init - self.MZR_g_init) * self.dt
            delta_term = self.sSFR*(self.Zgas - self.Zstar)*self.dt 
            MZR_evo_term = self.dZsdt * self.dt
            
            self.MZR_s_init += self.dZsdt * self.dt
            self.MZR_g_init += self.dZgdt * self.dt
            
            newZstar = self.Zstar + (MZR_diff_term + delta_term + MZR_evo_term)
            
            newZgas = self.Zgas + self.deltaZ
                        
            self.Zstars [self.index] = newZstar
            self.Zgasses[self.index] = newZgas 
            self.ts     [self.index] = newT
                        
            self.Zstar = newZstar 
            self.Zgas  = newZgas 
            self.t     = newT
            
            self.index += 1
        
    def run( self, ZgasPertb, coherenceTime ):
                
        self.ZgasPertb = ZgasPertb
        
        self.index = 0
        for index, newZgasTarget in enumerate(self.ZgasPertb):
            self.T += coherenceTime[index]
                        
            self.targetZ = newZgasTarget - self.MZR_g_init
            
            if (self.T - self.t)/self.dt < 1:
                # If coherence time is less than dt modify slightly
                # This is a bug fix, this happens ~1-5/10,000 times
                self.deltaZ = self.targetZ - self.Zgas
            else:
                self.deltaZ = self.dt*(self.targetZ - self.Zgas) / (self.T - self.t)
            self.doIntegration()
                        
        mask = self.ts > 0
        
        # Remove last step (is an extra 0)
        self.Zstars  = self.Zstars [mask]
        self.Zgasses = self.Zgasses[mask]
        self.ts      = self.ts     [mask]
    
    def plot(self):
        
        plt.hist2d( self.Zgasses, self.Zstars, bins=(100,100), norm=LogNorm() )

        plt.axhline(0,color='k')
        plt.axvline(0,color='k')

        x = self.Zgasses
        y = self.Zstars
        
        lm = LinearRegression(fit_intercept = False)
        
        lm.fit(x.reshape(-1, 1), y)
        X = np.arange( np.min(x), np.max(x), 0.1)
        preds = lm.predict( X.reshape(-1,1) )

        a,b = lm.coef_, lm.intercept_
        plt.plot(X, preds, color='k',lw=2.0)
        plt.plot(X, preds, color='w',label=f"Slope: {a}")
        plt.legend()
        
        plt.show()
        
    def time_series(self):
        
        plt.plot( self.ts, self.Zgasses )
        plt.plot( self.ts, self.Zstars )
        
        plt.axhline(0,color='k')
        
        plt.show()
        
def run_n_times(sSFR,dt,MZR_g_init,MZR_s_init,dZsdt,dZgdt,n=2):

    all_Zstar = []
    all_Zgas = []

    for _ in range(n):

        num_Zgas_jumps = 50
        Zpertb = 0.12 * np.random.randn( num_Zgas_jumps ) 

        coherenceTimes = np.random.exponential( scale=0.1*1/sSFR, size=num_Zgas_jumps )

        this_run = metals_MZR(0.12*np.random.randn(),0.12*np.random.randn(),
                              sSFR,dt,MZR_g_init,MZR_s_init,dZsdt,dZgdt)
        this_run.run(Zpertb, coherenceTimes)
        
        all_Zstar += list(this_run.Zstars)
        all_Zgas += list(this_run.Zgasses)
        
    plt.hist2d(all_Zgas, all_Zstar, bins=(100,100), norm=LogNorm() )
    
    
    plt.axhline(0,color='k')
    plt.axvline(0,color='k')

    x = np.array(all_Zstar)
    y = np.array(all_Zgas)

    lm = LinearRegression(fit_intercept = False)

    lm.fit(x.reshape(-1, 1), y)
    X = np.arange( np.min(x), np.max(x), 0.1)
    preds = lm.predict( X.reshape(-1,1) )

    a,b = lm.coef_, lm.intercept_
    plt.plot(X, preds, color='k',lw=2.0)
    plt.plot(X, preds, color='w',label=f"Slope: {a}")
    plt.legend()

    plt.show()