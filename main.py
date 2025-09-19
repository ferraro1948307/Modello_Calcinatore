# main.py — CLI entrypoint for the modular calciner model
# Usage examples (from project root, with venv activated):
#   python main.py --all
#   python main.py --fig 3
#   python main.py --custom --T 1183 --fCO2 0.8 --NCa 10000 --FCa 45 --F0 2.25 --Xcarb 0.2

import argparse
import os
import pandas as pd

from parameters import ModelParameters, ensure_dir
from model import Inputs, run_model, fig3, fig4, fig5, fig6

def main():
    ap = argparse.ArgumentParser(description="Calciner model (modular) — Figures 3–6 + custom run.")
    ap.add_argument("--out", default="outputs", help="Output directory (default: ./outputs)")
    ap.add_argument("--fig", choices=["3","4","5","6"], help="Generate only a specific figure.")
    ap.add_argument("--all", action="store_true", help="Generate all figures (3,4,5,6).")
    ap.add_argument("--custom", action="store_true", help="Run a custom single-point calculation.")
    params = ModelParameters()
    ap.add_argument("--T", type=float, default=params.T_default, help="Temperature [K] for custom run.")
    ap.add_argument("--P", type=float, default=params.P_default, help="Pressure [atm] for custom run.")
    ap.add_argument("--fCO2", type=float, default=params.fCO2_default, help="CO2 mole fraction for custom run.")
    ap.add_argument("--NCa", type=float, default=params.NCa_default, help="Solids inventory [mol/m^2] for custom run.")
    ap.add_argument("--FCa", type=float, default=params.FCa_default, help="Ca flow [mol/(m^2·s)] for custom run.")
    ap.add_argument("--F0", type=float, default=params.F0_default, help="Fresh make-up [mol/(m^2·s)] for custom run.")
    ap.add_argument("--Xcarb", type=float, default=params.Xcarb_default, help="CaCO3/Ca at calciner inlet [mol/mol].")

    args = ap.parse_args()
    outdir = ensure_dir(args.out)

    if args.all:
        fig3(outdir); fig4(outdir); fig5(outdir); fig6(outdir)
        print(f"All figures saved under: {outdir}")
        return

    if args.fig == "3":
        fig3(outdir); print(f"Fig. 3 saved under: {outdir}"); return
    if args.fig == "4":
        fig4(outdir); print(f"Fig. 4 saved under: {outdir}"); return
    if args.fig == "5":
        fig5(outdir); print(f"Fig. 5 saved under: {outdir}"); return
    if args.fig == "6":
        fig6(outdir); print(f"Fig. 6 saved under: {outdir}"); return

    if args.custom:
        inp = Inputs(T_K=args.T, P_atm=args.P, fCO2=args.fCO2, NCa=args.NCa,
                     FCa=args.FCa, F0=args.F0, Xcarb=args.Xcarb)
        out = run_model(inp)
        df = pd.DataFrame([{
            "T_K": inp.T_K, "P_atm": inp.P_atm, "fCO2": inp.fCO2,
            "NCa": inp.NCa, "FCa": inp.FCa, "F0": inp.F0, "Xcarb": inp.Xcarb,
            "feq": out.feq, "kcT": out.kcT, "t_star_s": out.t_star, "sigma_s_s": out.sigma_s,
            "fa": out.fa, "Ecalc": out.Ecalc, "Xcalc_out": out.Xcalc, "r_calc_avg": out.rcalc
        }])
        csv_path = os.path.join(outdir, "custom_run.csv")
        df.to_csv(csv_path, index=False)
        print("Custom run saved at:", csv_path)
        return

    ap.print_help()

if __name__ == "__main__":
    main()
