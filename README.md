# Calcinatore Model (Paper Replication)

Questo progetto replica le figure chiave (Fig. 3–6) del paper sul calcinatore.
Il modello implementa le equazioni principali per il calcolo dell'efficienza del calcinatore.

## Struttura del progetto
- `parameters.py` — costanti e parametri fissi
- `equations.py` — tutte le equazioni fondamentali (Barker, Arrhenius, t*, σ, f_a, E_calc, ecc.)
- `model.py` — esecuzione del modello e generatori delle figure
- `main.py` — interfaccia a linea di comando (CLI)
- `requirements.txt` — librerie necessarie

## Installazione
1. Creare un ambiente virtuale (opzionale ma consigliato):
   ```bash
   python -m venv .venv
   source .venv/bin/activate   # macOS/Linux
   .venv\Scripts\activate    # Windows
   ```

2. Installare le dipendenze:
   ```bash
   pip install -r requirements.txt
   ```

## Utilizzo
Esempi da terminale, nella cartella del progetto:

- Tutte le figure (3–6):
  ```bash
  python main.py --all
  ```

- Una figura specifica:
  ```bash
  python main.py --fig 3
  python main.py --fig 4
  python main.py --fig 5
  python main.py --fig 6
  ```

- Run personalizzato:
  ```bash
  python main.py --custom --T 1183 --fCO2 0.8 --NCa 10000 --FCa 45 --F0 2.25 --Xcarb 0.2
  ```

- Specificare cartella output:
  ```bash
  python main.py --all --out risultati
  ```

## Output
Le figure vengono salvate in formato `.png` e i dati in `.csv` nella cartella `outputs/` (default) o in quella specificata con `--out`.
