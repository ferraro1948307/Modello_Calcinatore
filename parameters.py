"""
Parametri del modello Calcinatore basati sul paper di Ortiz et al.
"""

import os

class ModelParameters:
    """Parametri del modello estratti dal paper di Ortiz et al."""
    
    def __init__(self):
        # Costanti dei gas
        self.R_SI = 8.314462618     # costante dei gas universale J/(mol·K)
        self.R_ATM = 8.205733e-5    # costante dei gas universale m³·atm/(mol·K)
        
        # Parametri cinetici (dal paper)
        self.kc0 = 2050.0           # fattore pre-esponenziale (m³/(mol·s))
        self.Ea_J = 112000.0        # energia di attivazione (J/mol)
        
        # Parametri operativi di default (Figura 3 baseline)
        self.T_default = 1183.0     # temperatura di default (K)
        self.P_default = 1.0        # pressione di default (atm)
        self.fCO2_default = 0.80    # frazione molare CO2 di default
        self.NCa_default = 10000.0  # inventario solidi di default (mol/m²)
        self.FCa_default = 45.0     # flusso Ca di default (mol/(m²·s))
        self.F0_default = 2.25      # make-up fresco di default (mol/(m²·s))
        self.Xcarb_default = 0.20   # CaCO3/Ca all'ingresso di default
        
        # Parametri per Figura 3
        self.fig3_temperatures = [1163.0, 1168.0, 1173.0, 1183.0, 1193.0]  # temperature per Fig. 3
        self.fig3_N_range = (3000.0, 15000.0, 25)                          # range NCa per Fig. 3
        
        # Scenari per Figura 4
        self.fig4_scenarios = [
            {"name": "low make-up (F0/FCa=0.008, Xcarb=0.1)", "FCa": 90.0, "F0": 0.72, "Xcarb": 0.10},
            {"name": "high make-up (F0/FCa=0.12, Xcarb=0.3)", "FCa": 30.0, "F0": 3.60, "Xcarb": 0.30},
        ]
        self.fig4_T_range = (1160.0, 1195.0, 71)                          # range temperatura per Fig. 4
        
        # Parametri per Figura 5
        self.fig5_fCO2_list = [0.70, 0.80, 0.90, 1.00]                    # liste fCO2 per Fig. 5
        self.fig5_T_range = (1160.0, 1195.0, 71)                          # range temperatura per Fig. 5
        
        # Parametri per Figura 6
        self.fig6_T_range = (1165.0, 1185.0, 121)                         # range temperatura per Fig. 6
        self.fig6_L_range = (6.0, 11.0, 101)                              # range carico per Fig. 6
        self.fig6_ratio_F0_FCa = 0.05                                      # rapporto F0/FCa per Fig. 6
        self.fig6_efficiency_levels = [0.95, 0.96, 0.97]                  # livelli efficienza per Fig. 6

def ensure_dir(path: str):
    """Utility per creare directory se non esistono"""
    os.makedirs(path, exist_ok=True)
    return path
