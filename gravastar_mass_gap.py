#!/usr/bin/env python3
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass
"""
Gravastar mass-gap population model.
STATUS: STUB (R4-M2). Full population synthesis in companion paper.
Reference: Manuscript §5; App D.2; Table 9 row 4.
"""
import numpy as np

c, G, M_sun, PI = 2.998e8, 6.674e-11, 1.989e30, np.pi

def f_schw(M): return c**3 / (2*PI*6.29*G*M*M_sun)
def f_ve(M):   return PI*c**3 / (32*G*M*M_sun)

def main():
    print("Gravastar mass-gap population — STUB (R4-M2)")
    print(f"{'M[Msol]':>7} {'f_Schw[Hz]':>11} {'f_VE[Hz]':>10} {'Ratio':>7} {'Band':>10}")
    print("-"*50)
    for M in [2.5,3,3.5,4,5,7,10,15,30]:
        fs, fv = f_schw(M), f_ve(M)
        band = "LIGO" if 10<=fv<=5000 else "LISA" if fv<10 else "CE/ET"
        gap = " <GAP" if 2.5<=M<=5 else ""
        print(f"{M:>7.1f} {fs:>11.0f} {fv:>10.0f} {fv/fs:>7.3f} {band:>10}{gap}")
    print(f"\nf_VE/f_Schw = {PI**2/(16*0.37367):.3f} (mass-independent)")
    print("WARNING: Needs: progenitor f(M), binary rates, NR waveforms")

if __name__ == "__main__":
    main()
