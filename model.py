"""
Modello per il calcolo dell'efficienza di cattura di CO2 — Equazioni 8-23 da Ortiz et al. (2015)
"""

from dataclasses import dataclass
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from parameters import ModelParameters, ensure_dir
from equations import CalcinerModelEquation

class CalcinerModel:
    """Modello per il calcolo dell'efficienza di cattura di CO2 — Equazioni 8-23 da Ortiz et al. (2015)"""
    
    def __init__(self):
        self.params = ModelParameters()
        self.equations = CalcinerModelEquation()

@dataclass
class Inputs:
    T_K: float = None
    P_atm: float = None
    fCO2: float = None
    NCa: float = None
    FCa: float = None
    F0: float = None
    Xcarb: float = None
    
    def __post_init__(self):
        """Inizializza i valori di default se non forniti"""
        params = ModelParameters()
        if self.T_K is None:
            self.T_K = params.T_default
        if self.P_atm is None:
            self.P_atm = params.P_default
        if self.fCO2 is None:
            self.fCO2 = params.fCO2_default
        if self.NCa is None:
            self.NCa = params.NCa_default
        if self.FCa is None:
            self.FCa = params.FCa_default
        if self.F0 is None:
            self.F0 = params.F0_default
        if self.Xcarb is None:
            self.Xcarb = params.Xcarb_default

@dataclass
class Outputs:
    feq: float
    kcT: float
    t_star: float
    sigma_s: float
    fa: float
    Ecalc: float
    Xcalc: float
    rcalc: float

def run_model(inp: Inputs) -> Outputs:
    """Esegue il modello calcinatore con i parametri di input specificati"""
    equations = CalcinerModelEquation()
    
    feq = equations.fe_barker(inp.T_K, inp.P_atm)
    kcT = equations.kc_arrhenius(inp.T_K)
    tstar = equations.t_c_star(inp.T_K, inp.P_atm, inp.fCO2, kcT, inp.Xcarb)
    sigma_s = equations.sigma_residence(inp.NCa, inp.FCa, inp.F0)
    fa = equations.f_active_fraction(tstar, sigma_s)
    E = equations.E_calc_from_fa(fa)
    X_out = equations.Xcalc_out(inp.Xcarb, E)
    r_avg = equations.r_calc_avg(kcT, inp.T_K, inp.P_atm, inp.fCO2)
    return Outputs(feq, kcT, tstar, sigma_s, fa, E, X_out, r_avg)

# -------------- Figure generators (save PNG + CSV) --------------

def fig3(output_dir: str):
    """
    Figura 3: Efficienza vs N_Ca per diverse temperature
    Reproduce il grafico del paper di Ortiz et al. con stile professionale
    """
    ensure_dir(output_dir)
    params = ModelParameters()
    equations = CalcinerModelEquation()
    
    P = params.P_default
    fCO2 = params.fCO2_default
    Xcarb = params.Xcarb_default
    FCa = params.FCa_default
    F0 = params.F0_default
    T_list = params.fig3_temperatures
    N_range = params.fig3_N_range
    N_vals = np.linspace(N_range[0], N_range[1], N_range[2])

    # Calcolo dei dati
    rows = []
    plt.style.use('default')  # Reset dello stile
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Colori distintivi per ogni temperatura
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    markers = ['o', 's', '^', 'D', 'v']
    
    for i, T in enumerate(T_list):
        E_vals = []
        for N in N_vals:
            kT = equations.kc_arrhenius(T)
            t  = equations.t_c_star(T, P, fCO2, kT, Xcarb)
            sig= equations.sigma_residence(N, FCa, F0)
            fa = equations.f_active_fraction(t, sig)
            E  = equations.E_calc_from_fa(fa)
            rows.append({"T_K":T, "NCa_mol_per_m2":N, "Ecalc":E, "t_star_s":t, "sigma_s":sig})
            E_vals.append(E)
        
        # Plot con stile migliorato
        ax.plot(N_vals, E_vals, 
                color=colors[i % len(colors)], 
                marker=markers[i % len(markers)],
                markersize=4,
                linewidth=2,
                markevery=3,  # Mostra marker ogni 3 punti
                label=f"{int(T)} K")

    # Personalizzazione degli assi
    ax.set_xlabel("N$_{Ca}$ [mol/m²]", fontsize=12, fontweight='bold')
    ax.set_ylabel("E$_{calc}$ [-]", fontsize=12, fontweight='bold')
    ax.set_title("Figura 3: Efficienza del Calcinatore vs Inventario Solidi\n" + 
                f"(f$_{{CO2}}$ = {fCO2}, F$_0$/F$_{{Ca}}$ = {F0/FCa:.3f}, Xcarb = {Xcarb})", 
                fontsize=13, fontweight='bold', pad=20)
    
    # Griglia e leggenda
    ax.grid(True, which="major", linestyle="-", alpha=0.3)
    ax.grid(True, which="minor", linestyle=":", alpha=0.2)
    ax.legend(title="Temperatura", fontsize=10, title_fontsize=11, 
              loc='lower right', frameon=True, fancybox=True, shadow=True)
    
    # Range degli assi per migliore visualizzazione
    ax.set_xlim(N_vals[0], N_vals[-1])
    ax.set_ylim(0.82, 1.02)
    
    # Personalizzazione tick
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    
    plt.tight_layout()
    
    # Salvataggio con alta qualità
    import os
    plt.savefig(os.path.join(output_dir, "fig3_E_vs_NCa.png"), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, "fig3_E_vs_NCa.pdf"), 
                bbox_inches='tight', facecolor='white')  # Anche in PDF
    plt.close()
    
    # Salvataggio dati
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, "fig3_E_vs_NCa.csv"), index=False)
    print(f"Figura 3 salvata in {output_dir} (PNG e PDF)")

def fig4(output_dir: str):
    """
    Figura 4: Efficienza vs Temperatura per diversi scenari di make-up
    Confronta scenari con diversi rapporti F0/FCa e Xcarb
    """
    ensure_dir(output_dir)
    params = ModelParameters()
    equations = CalcinerModelEquation()
    
    P = params.P_default
    NCa = params.NCa_default
    fCO2 = params.fCO2_default
    scenarios = params.fig4_scenarios
    T_range = params.fig4_T_range
    T_grid = np.linspace(T_range[0], T_range[1], T_range[2])

    # Calcolo dei dati
    rows = []
    plt.style.use('default')  # Reset dello stile
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Colori e stili distintivi per ogni scenario
    colors = ['#e74c3c', '#3498db']  # Rosso e blu
    linestyles = ['-', '--']
    markers = ['o', 's']
    
    for i, sc in enumerate(scenarios):
        E_vals = []
        for T in T_grid:
            kT = equations.kc_arrhenius(T)
            t  = equations.t_c_star(T, P, fCO2, kT, sc["Xcarb"])
            sig= equations.sigma_residence(NCa, sc["FCa"], sc["F0"])
            fa = equations.f_active_fraction(t, sig)
            E  = equations.E_calc_from_fa(fa)
            rows.append({"T_K":T, "scenario":sc["name"], "Ecalc":E, "t_star_s":t, "sigma_s":sig})
            E_vals.append(E)
        
        # Etichette più pulite
        if "low" in sc["name"].lower():
            label = f"Basso make-up (F₀/F_Ca = {sc['F0']/sc['FCa']:.3f})"
        else:
            label = f"Alto make-up (F₀/F_Ca = {sc['F0']/sc['FCa']:.2f})"
        
        # Plot con stile migliorato
        ax.plot(T_grid, E_vals, 
                color=colors[i], 
                linestyle=linestyles[i],
                marker=markers[i],
                markersize=4,
                linewidth=2.5,
                markevery=8,  # Mostra marker ogni 8 punti
                label=label)

    # Personalizzazione degli assi
    ax.set_xlabel("Temperatura [K]", fontsize=12, fontweight='bold')
    ax.set_ylabel("E$_{calc}$ [-]", fontsize=12, fontweight='bold')
    ax.set_title("Figura 4: Efficienza del Calcinatore vs Temperatura\n" + 
                f"(N_Ca = {NCa:.0f} mol/m², f_CO₂ = {fCO2})", 
                fontsize=13, fontweight='bold', pad=20)
    
    # Griglia e leggenda
    ax.grid(True, which="major", linestyle="-", alpha=0.3)
    ax.grid(True, which="minor", linestyle=":", alpha=0.2)
    ax.legend(fontsize=10, loc='lower right', frameon=True, fancybox=True, shadow=True)
    
    # Range degli assi per migliore visualizzazione
    ax.set_xlim(T_grid[0], T_grid[-1])
    ax.set_ylim(0.8, 1.0)
    
    # Personalizzazione tick
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    
    plt.tight_layout()
    
    # Salvataggio con alta qualità
    import os
    plt.savefig(os.path.join(output_dir, "fig4_E_vs_T_makeup.png"), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, "fig4_E_vs_T_makeup.pdf"), 
                bbox_inches='tight', facecolor='white')
    plt.close()
    
    # Salvataggio dati
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, "fig4_E_vs_T_makeup.csv"), index=False)
    print(f"Figura 4 salvata in {output_dir} (PNG e PDF)")

def fig5(output_dir: str):
    """
    Figura 5: Efficienza vs Temperatura per diverse frazioni molari di CO2
    Mostra l'effetto della composizione del gas sulla performance del calcinatore
    """
    ensure_dir(output_dir)
    params = ModelParameters()
    equations = CalcinerModelEquation()
    
    P = params.P_default
    NCa = params.NCa_default
    Xcarb = params.Xcarb_default
    FCa = params.FCa_default
    F0 = params.F0_default
    fCO2_list = params.fig5_fCO2_list
    T_range = params.fig5_T_range
    T_grid = np.linspace(T_range[0], T_range[1], T_range[2])

    # Calcolo dei dati
    rows = []
    plt.style.use('default')  # Reset dello stile
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Colori dal viola al rosso (scala di concentrazione CO2)
    colors = ['#9b59b6', '#3498db', '#2ecc71', '#f39c12']  # Viola, blu, verde, arancione
    markers = ['o', 's', '^', 'D']
    linestyles = ['-', '-', '-', '-']
    
    for i, fco2 in enumerate(fCO2_list):
        E_vals = []
        for T in T_grid:
            kT = equations.kc_arrhenius(T)
            t  = equations.t_c_star(T, P, fco2, kT, Xcarb)
            sig= equations.sigma_residence(NCa, FCa, F0)
            fa = equations.f_active_fraction(t, sig)
            E  = equations.E_calc_from_fa(fa)
            rows.append({"T_K":T, "fCO2":fco2, "Ecalc":E, "t_star_s":t, "sigma_s":sig})
            E_vals.append(E)
        
        # Plot con stile migliorato
        ax.plot(T_grid, E_vals, 
                color=colors[i], 
                linestyle=linestyles[i],
                marker=markers[i],
                markersize=4,
                linewidth=2.5,
                markevery=8,  # Mostra marker ogni 8 punti
                label=f"f_CO₂ = {fco2:.2f}")

    # Personalizzazione degli assi
    ax.set_xlabel("Temperatura [K]", fontsize=12, fontweight='bold')
    ax.set_ylabel("E$_{calc}$ [-]", fontsize=12, fontweight='bold')
    ax.set_title("Figura 5: Efficienza del Calcinatore vs Temperatura\n" + 
                f"(N_Ca = {NCa:.0f} mol/m², F₀/F_Ca = {F0/FCa:.3f}, X_carb = {Xcarb})", 
                fontsize=13, fontweight='bold', pad=20)
    
    # Griglia e leggenda
    ax.grid(True, which="major", linestyle="-", alpha=0.3)
    ax.grid(True, which="minor", linestyle=":", alpha=0.2)
    ax.legend(title="Frazione molare CO₂", fontsize=10, title_fontsize=11,
              loc='lower right', frameon=True, fancybox=True, shadow=True)
    
    # Range degli assi per migliore visualizzazione
    ax.set_xlim(T_grid[0], T_grid[-1])
    ax.set_ylim(0.88, 1.0)
    
    # Personalizzazione tick
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    
    plt.tight_layout()
    
    # Salvataggio con alta qualità
    import os
    plt.savefig(os.path.join(output_dir, "fig5_E_vs_T_fCO2.png"), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, "fig5_E_vs_T_fCO2.pdf"), 
                bbox_inches='tight', facecolor='white')
    plt.close()
    
    # Salvataggio dati
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, "fig5_E_vs_T_fCO2.csv"), index=False)
    print(f"Figura 5 salvata in {output_dir} (PNG e PDF)")

def fig6(output_dir: str):
    """
    Figura 6: Efficienza del calcinatore in funzione di NCa e FCa·Xcarb per diverse temperature
    Replica il formato della figura del paper con linee iso-efficienza
    """
    ensure_dir(output_dir)
    params = ModelParameters()
    equations = CalcinerModelEquation()
    
    P = params.P_default
    fCO2 = params.fCO2_default
    Xcarb = params.Xcarb_default
    ratio = params.fig6_ratio_F0_FCa  # F0/FCa = 0.05
    
    # Temperature da testare (seleziono alcune da fig3)
    T_list = [1168.0, 1173.0, 1183.0]  # K
    
    # Range per gli assi
    FCa_Xcarb_vals = np.linspace(6.0, 11.0, 50)  # FCa·Xcarb [mol CO2/(m²·s)]
    NCa_vals = np.linspace(2500.0, 22500.0, 200)  # NCa [mol/m²]
    
    # Livelli di iso-efficienza
    efficiency_levels = [0.95, 0.97]
    
    # Creazione del grafico
    plt.style.use('default')
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Colori e stili per diverse temperature
    colors = ['black', 'gray', 'darkgray']
    linestyles = ['-', '-', '-']
    linewidths = [2.0, 2.0, 2.0]
    
    rows = []
    
    # Per ogni temperatura
    for i, T in enumerate(T_list):
        print(f"Calcolando per T = {T} K...")
        
        # Per ogni livello di efficienza
        for E_target in efficiency_levels:
            NCa_isoline = []
            FCa_Xcarb_isoline = []
            
            # Per ogni valore di FCa·Xcarb, trova NCa che dà l'efficienza target
            for FCa_Xcarb in FCa_Xcarb_vals:
                FCa = FCa_Xcarb / Xcarb  # FCa = (FCa·Xcarb) / Xcarb
                F0 = ratio * FCa
                
                # Trova NCa che dà efficienza E_target
                for NCa in NCa_vals:
                    kT = equations.kc_arrhenius(T)
                    t_star = equations.t_c_star(T, P, fCO2, kT, Xcarb)
                    sigma_s = equations.sigma_residence(NCa, FCa, F0)
                    fa = equations.f_active_fraction(t_star, sigma_s)
                    E = equations.E_calc_from_fa(fa)
                    
                    # Salva tutti i dati per CSV
                    rows.append({
                        "T_K": T, "NCa_mol_per_m2": NCa, "FCa_Xcarb": FCa_Xcarb,
                        "FCa": FCa, "F0": F0, "Ecalc": E
                    })
                    
                    # Se siamo vicini all'efficienza target, aggiungi il punto
                    if abs(E - E_target) < 0.005:  # Tolleranza di ±0.5%
                        NCa_isoline.append(NCa)
                        FCa_Xcarb_isoline.append(FCa_Xcarb)
                        break
            
            # Plot della linea iso-efficienza
            if len(NCa_isoline) > 3:  # Solo se abbiamo abbastanza punti
                label = f"E$_{{calc}}$ = {E_target:.2f}" if i == 0 else ""  # Label solo per la prima temperatura
                ax.plot(FCa_Xcarb_isoline, NCa_isoline, 
                       color=colors[i], linestyle=linestyles[i], 
                       linewidth=linewidths[i], label=label)
        
        # Aggiungi annotazione della temperatura
        # Trova un punto sulla curva E=0.97 per posizionare l'etichetta
        for FCa_Xcarb in [7.5, 9.0, 10.5]:  # Prova diverse posizioni
            FCa = FCa_Xcarb / Xcarb
            F0 = ratio * FCa
            for NCa in np.linspace(8000, 18000, 20):
                kT = equations.kc_arrhenius(T)
                t_star = equations.t_c_star(T, P, fCO2, kT, Xcarb)
                sigma_s = equations.sigma_residence(NCa, FCa, F0)
                fa = equations.f_active_fraction(t_star, sigma_s)
                E = equations.E_calc_from_fa(fa)
                
                if abs(E - 0.97) < 0.01:  # Vicino alla linea E=0.97
                    ax.annotate(f'{int(T)} K', xy=(FCa_Xcarb, NCa), 
                               xytext=(FCa_Xcarb+0.3, NCa+1000),
                               fontsize=12, fontweight='bold',
                               arrowprops=dict(arrowstyle='->', color=colors[i], lw=1.5))
                    break
            else:
                continue
            break
    
    # Personalizzazione degli assi
    ax.set_xlabel("F$_{Ca}$ · X$_{carb}$ [mol CO₂/(m²·s)]", fontsize=14, fontweight='bold')
    ax.set_ylabel("N$_{Ca}$ [mol/m²]", fontsize=14, fontweight='bold')
    ax.set_title("Fig. 6. Calciner efficiency as a function of solid inventory (N$_{Ca}$) and CO₂ captured in\n" +
                f"the carbonator (F$_{{Ca}}$ · X$_{{carb}}$) at different temperatures (X$_{{carb}}$ = {Xcarb}, f$_{{CO2}}$ = {fCO2})",
                fontsize=12, fontweight='bold', pad=20)
    
    # Griglia e leggenda
    ax.grid(True, which="major", linestyle="-", alpha=0.3)
    ax.grid(True, which="minor", linestyle=":", alpha=0.2)
    ax.legend(fontsize=12, loc='upper left', frameon=True, fancybox=True, shadow=True)
    
    # Range degli assi
    ax.set_xlim(6, 11)
    ax.set_ylim(2500, 22500)
    
    # Personalizzazione tick
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=10)
    
    plt.tight_layout()
    
    # Salvataggio con alta qualità
    import os
    plt.savefig(os.path.join(output_dir, "fig6_isoefficiency.png"), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, "fig6_isoefficiency.pdf"), 
                bbox_inches='tight', facecolor='white')
    plt.close()
    
    # Salvataggio dati
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, "fig6_isoefficiency.csv"), index=False)
    print(f"Figura 6 salvata in {output_dir} (PNG e PDF)")
