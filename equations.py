"""
Implementazione delle equazioni del modello calcinatore
"""

import math
import numpy as np
from parameters import ModelParameters

class CalcinerModelEquation:
    """Implementazione delle equazioni del modello calcinatore"""
    
    def __init__(self):
        self.params = ModelParameters()
    
    def fe_barker(self, T_K, P_atm):
        """
        Equilibrio di Barker (Eq. 11): f_e(T,P) = 10^(7.079 - 8308/T) / P.
        
        Args:
            T_K (float): Temperatura [K]
            P_atm (float): Pressione [atm]
            
        Returns:
            float: Frazione di equilibrio CO2
        """
        return (10.0 ** (7.079 - 8308.0 / T_K)) / P_atm
    
    def kc_arrhenius(self, T_K):
        """
        Equazione di Arrhenius per k_c(T).
        
        Args:
            T_K (float): Temperatura [K]
            
        Returns:
            float: Costante cinetica [m³/(mol·s)]
        """
        return self.params.kc0 * math.exp(-self.params.Ea_J / (self.params.R_SI * T_K))
    
    def t_c_star(self, T_K, P_atm, fCO2, kc_T, Xcarb, eps=1e-12):
        """
        Tempo per completare la calcinazione per particella (Eq. 7).
        
        Args:
            T_K (float): Temperatura [K]
            P_atm (float): Pressione [atm]
            fCO2 (float): Frazione molare CO2
            kc_T (float): Costante cinetica alla temperatura T
            Xcarb (float): CaCO3/Ca alla particella
            eps (float): Tolleranza numerica
            
        Returns:
            float: Tempo caratteristico [s]
        """
        C_T = P_atm / (self.params.R_ATM * T_K)  # mol/m³
        C_eq = self.fe_barker(T_K, P_atm) * C_T
        C_co2 = fCO2 * C_T
        driving = C_eq - C_co2
        if driving <= eps or kc_T <= eps:
            return float('inf')
        return 3.0 * Xcarb / (kc_T * driving)
    
    def sigma_residence(self, NCa, FCa, F0, eps=1e-12):
        """
        Tempo medio di residenza dei solidi nel calcinatore (Eq. 10).
        
        Args:
            NCa (float): Inventario solidi [mol/m²]
            FCa (float): Flusso Ca [mol/(m²·s)]
            F0 (float): Make-up fresco [mol/(m²·s)]
            eps (float): Tolleranza numerica
            
        Returns:
            float: Tempo di residenza medio [s]
        """
        denom = FCa + F0
        if denom <= eps:
            return float('inf')
        return NCa / denom
    
    def f_active_fraction(self, t_star, sigma_s):
        """
        Frazione attiva (Eq. 9) sotto solidi perfettamente miscelati.
        
        Args:
            t_star (float): Tempo caratteristico [s]
            sigma_s (float): Tempo di residenza [s]
            
        Returns:
            float: Frazione attiva [-]
        """
        if not np.isfinite(t_star) or not np.isfinite(sigma_s) or sigma_s <= 0.0:
            return 0.0
        return 1.0 - math.exp(-t_star / sigma_s)
    
    def E_calc_from_fa(self, fa, eps=1e-15):
        """
        Efficienza del calcinatore (Eq. 17).
        
        Args:
            fa (float): Frazione attiva [-]
            eps (float): Tolleranza numerica
            
        Returns:
            float: Efficienza calcinatore [-]
        """
        fa = min(max(fa, 0.0), 1.0 - 1e-15)
        denom = math.log(1.0 / (1.0 - fa))
        if denom <= eps:
            return 0.0
        return fa / denom
    
    def Xcalc_out(self, Xcarb_in, E_calc):
        """
        CaCO3/Ca in uscita dai solidi (Eq. 1).
        
        Args:
            Xcarb_in (float): CaCO3/Ca in ingresso [mol/mol]
            E_calc (float): Efficienza calcinatore [-]
            
        Returns:
            float: CaCO3/Ca in uscita [mol/mol]
        """
        E_calc = min(max(E_calc, 0.0), 1.0)
        return Xcarb_in * (1.0 - E_calc)
    
    def r_calc_avg(self, kc_T, T_K, P_atm, fCO2):
        """
        Velocità media di calcinazione (Eq. 8) — utile per interpretazione.
        
        Args:
            kc_T (float): Costante cinetica alla temperatura T
            T_K (float): Temperatura [K]
            P_atm (float): Pressione [atm]
            fCO2 (float): Frazione molare CO2
            
        Returns:
            float: Velocità media di calcinazione
        """
        C_T = P_atm / (self.params.R_ATM * T_K)
        C_eq = self.fe_barker(T_K, P_atm) * C_T
        C_co2 = fCO2 * C_T
        driving = max(C_eq - C_co2, 0.0)
        return kc_T * driving / 3.0
