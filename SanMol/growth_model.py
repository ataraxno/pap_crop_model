import numpy as np

class CropTemplate:
    def __init__(self):
        ## Monitoring ##
        self.dat = 0
        self.thermal_time = 0
        self.phot = 0
        self.resp = 0
        self.dwl = 0 # dry weight of leaves (g m^-2)
        self.dws = 0
        self.dwr = 0
        self.dwf = 0
        
        ## LAI ##
        self.lai_max = 6
        self.lai_ini = 0.1
        self.lai = self.lai_ini
        self.Cp_lai = 0.00018 # Tuning parameter of Gompertz function
        
        ## Photosnythesis ##
        self.Vc_max = 200
        self.gamma_star = None # calculated in photosynthesis function
        self.Kc = 300
        self.O = 21000
        self.Ko = 300
        self.eq = 0.063
        self.alpha = 0.76 # in paper; 0.953 = empricial value
        self.R = 8.314 # Universal gas constant
    
        ## Maintenance respiration ##
        self.respiration_rate = None # calculated in maintenance_respiration function
        self.refer_temp = 25
        
        ## Dry weight ##
        self.fc = 0.29 # A conversion factor for g CH2O to g DW
        self.ff = None # fractions of the total daily dry weight to fruits; calculated in life_cycle function
        self.fr = None # fractions of the total daily dry weight to roots; calculated in life_cycle function
        self.vegitative_stage = 30 # days of vegitative stage
        
    
    def LAI(self, thermal_time): # Daily
        self.lai = self.lai_max*np.exp((-self.lai_max/self.lai_ini)*np.exp(-self.Cp_lai*thermal_time))
        return self.lai
    
    
    def photosynthesis(self, temp_K, If, C): # Hourly
        self.gamma_star = 42.75*np.exp(37830*(temp_K-298)/(298*self.R*temp_K))
        Vc = self.Vc_max*((C - self.gamma_star)/(C + self.Kc*(1 + self.O/self.Ko)))
        Vq = self.eq*self.alpha*If*((C - self.gamma_star)/(C + 2*self.gamma_star))
        Vs = self.Vc_max/2
                
        photo_leaf = min(Vc, Vq, Vs)
        photo_crop = photo_leaf*self.lai
        return photo_crop
    
    
    def maintenance_respiration(self, temp): # Hourly
        self.respiration_rate = self.dwl*0.03 + self.dws*0.015 + self.dwr*0.015 + self.dwf*0.01
        main_resp = self.respiration_rate*(2**((temp - self.refer_temp)/10))
        return main_resp
            
    
    def fraction_sinusoidal(self, DAT, a, b, c, d):
        fraction = a + b*np.sin(2*np.pi*(DAT + c)/d)
        return fraction
    
        
    def life_cycle(self, temp, rh, co2, rad): # Environment have to be hourly input.
        self.thermal_time += (temp - 10).sum() # input will be considered as numpy.array
        intercellular_C = co2*(1 - 1/(6.698*rh)) # m = 6.698 from empirical data (NOT IN THE PAPER)
        temp_K = temp + 273.15 # Kelvin temperature
        lai = self.LAI(self.thermal_time)
        
        phot = 0
        resp = 0
        for _ in range(24): # 24 hours
            phot += self.photosynthesis(temp_K[_], rad[_], intercellular_C[_])
            resp += self.maintenance_respiration(temp[_])
        
        self.phot = phot
        self.resp = resp
        daily_tot_dw = self.fc*(phot - resp)
        if daily_tot_dw < 0:
            daily_tot_dw = 0
            
        if self.dat < self.vegitative_stage:
            self.ff = 0
            self.fr = 0.35
            self.fl = 0.31
            self.fs = 0.34
        else:
            self.ff = self.fraction_sinusoidal(self.dat, 0.48, 0.17, 60, 79)
            self.fr = self.fraction_sinusoidal(self.dat, 0.17, 0.17, 25, 79)
            self.fl = 0.045*self.ff + 0.133 # Linear relation from fig. 5
            self.fs = 1 - self.ff - self.fr - self.fl
        self.dwf += self.ff*daily_tot_dw
        self.dwr += self.fr*daily_tot_dw
        self.dwl += self.fl*daily_tot_dw
        self.dws += self.fs*daily_tot_dw
        
        self.dat += 1
        return lai, phot, resp
        
    def status(self):
        print(f'<{self.dat} DAT>')
        print()
        print('# Assimilation #')
        print(f'Photosynthesis: {self.phot} / Respiration: {self.resp}')
        print()
        print('# Dry weights #')
        print(f'leaves: {self.dwl} / stems: {self.dws} / roots: {self.dwr} / fruits: {self.dwf}')
        print()
        print('# LAI #')
        print(f'LAI: {self.lai}')
    
    
class SweetPepper(CropTemplate):
    def __init__(self):
        super(SweetPepper, self).__init__()
        ## Monitoring ##
        self.dat = 0
        self.thermal_time = 0
        self.phot = 0
        self.resp = 0
        self.dwl = 0 # dry weight of leaves (g m^-2)
        self.dws = 0
        self.dwr = 0
        self.dwf = 0
        
        ## LAI ##
        self.lai_max = 6
        self.lai_ini = 0.1
        self.lai = self.lai_ini
        self.Cp_lai = 0.00018 # Tuning parameter of Gompertz function
        
        ## Photosnythesis ##
        self.Vc_max = 200
        self.gamma_star = None # calculated in photosynthesis function
        self.Kc = 300
        self.O = 21000
        self.Ko = 300
        self.eq = 0.063
        self.alpha = 0.76 # in paper; 0.953 = empricial value
        self.R = 8.314 # Universal gas constant
    
        ## Maintenance respiration ##
        self.respiration_rate = None # calculated in maintenance_respiration function
        self.refer_temp = 25
        
        ## Dry weight ##
        self.fc = 0.29 # A conversion factor for g CH2O to g DW
        self.ff = None # fractions of the total daily dry weight to fruits; calculated in life_cycle function
        self.fr = None # fractions of the total daily dry weight to roots; calculated in life_cycle function
        self.vegitative_stage = 30 # days of vegitative stage