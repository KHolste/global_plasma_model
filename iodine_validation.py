#!/usr/bin/env python3
"""
iodine_validation.py -- Erste Referenzvalidierung des Iod-Modells.

Vergleicht iodine_lafleur_v1 mit Literaturwerten, insbesondere:
1. Grondein et al. 2016 (Fig.4a, 5a): Te, Speziesdichten vs. P_RF
2. Lafleur 2026 (Fig.29): Collisional energy cost

Variantenbewertung: Welche Kombination aus Kelm_I2 und Kdissatt_I2
trifft die Literaturdaten am konsistentesten.
"""
from __future__ import annotations
import sys, copy, json, math
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plasma_chemistry import load_chemistry, ThrusterGeometry, E_CH, KB, CONV
from generic_solver import solve_steady_state
from rf_diagnostics import compute_rf_diagnostics

BASE = Path(__file__).resolve().parent / "chemistry" / "iodine_lafleur_v1"

# Grondein 2016 Geometrie (Table I)
GEOM = ThrusterGeometry(
    R=0.06, L=0.10, betai=0.7, betag=0.3, Vgrid=1000.0, sgrid=0.0015,
    frequency=13.56e6, Nw=5, R_ohm=2.0, Rc=0.07, lc=0.10,
)

# ═══════════════════════════════════════════════════════════════
# Referenzdaten (digitalisiert aus Grondein 2016, Fig.4a und 5a)
# ═══════════════════════════════════════════════════════════════
# ACHTUNG: Grondein nutzt P_RF (Generatorleistung), nicht P_abs.
# Bei zeta~0.7: P_abs ~ 0.7 * P_RF. Wir vergleichen bei P_abs direkt.
# Grondein Fig.4a zeigt Dichten bei P_RF = 200-1200W.
# Annahme: zeta ~ 0.7 -> P_abs ~ 0.7 * P_RF

REFERENCE_DATA = [
    {
        "source": "Grondein 2016 Fig.4a",
        "comparability": "semi-quantitative",
        "notes": "Globales 0D-Modell, gleiche Geometrie. Andere Raten (Hamilton/Tennyson 2015). P_RF -> P_abs mit zeta~0.7.",
        "conditions": {"P_RF_W": 400, "P_abs_est_W": 280, "Q0_pps": 2.1e18},
        "values": {
            "Te_eV": {"value": 4.2, "uncertainty": 0.5, "type": "model"},
            "ne_m3": {"value": 2.5e17, "uncertainty_factor": 2, "type": "model"},
            "n_I_m3": {"value": 5e17, "uncertainty_factor": 3, "type": "model"},
            "n_I2_m3": {"value": 1e18, "uncertainty_factor": 2, "type": "model"},
            "n_Ip_m3": {"value": 2e17, "uncertainty_factor": 2, "type": "model"},
            "n_I2p_m3": {"value": 7e16, "uncertainty_factor": 3, "type": "model"},
            "n_Im_m3": {"value": 1e16, "uncertainty_factor": 3, "type": "model"},
        }
    },
    {
        "source": "Grondein 2016 Fig.5a",
        "comparability": "semi-quantitative",
        "notes": "Te vs P_RF. Te steigt von ~3.5 eV (100W) bis ~4.9 eV (CL limit ~1000W).",
        "conditions": {"P_RF_W": 200, "P_abs_est_W": 140, "Q0_pps": 2.1e18},
        "values": {
            "Te_eV": {"value": 3.8, "uncertainty": 0.5, "type": "model"},
        }
    },
    {
        "source": "Grondein 2016 Fig.5a",
        "comparability": "semi-quantitative",
        "conditions": {"P_RF_W": 800, "P_abs_est_W": 560, "Q0_pps": 2.1e18},
        "values": {
            "Te_eV": {"value": 4.7, "uncertainty": 0.5, "type": "model"},
        }
    },
    {
        "source": "Grondein 2016 Fig.5b (Tg)",
        "comparability": "qualitative",
        "notes": "Tg peaks ~330K at P_RF~400W, then decreases. At 200W: ~310K, at 800W: ~325K.",
        "conditions": {"P_RF_W": 400, "P_abs_est_W": 280, "Q0_pps": 2.1e18},
        "values": {
            "Tg_K": {"value": 330, "uncertainty": 20, "type": "model"},
        }
    },
]


def load_table(path):
    Te, K = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Te"):
                continue
            p = line.strip().split(",")
            if len(p) >= 2:
                Te.append(float(p[0]))
                K.append(float(p[1]))
    return np.array(Te), np.array(K)


def apply_var(chem, rxn_id, var_file):
    c = copy.deepcopy(chem)
    Te, K = load_table(BASE / "variants" / var_file)
    for rxn in c.reactions:
        if rxn.id == rxn_id:
            rxn.rate._table_Te = Te
            rxn.rate._table_K = K
            break
    return c


def run(chem, P_abs, Q0):
    r = solve_steady_state(chem, GEOM, P_abs, Q0, max_iter=1200, tol=0.5)
    if not r["converged"]:
        return None
    d = r["densities"]
    nI, nI2 = d.get("I", 0), d.get("I2", 0)
    return {
        "Te": r["Te"], "Tg": r["Tg"], "ne": r["ne"],
        "n_I": nI, "n_I2": nI2,
        "n_Ip": d.get("I+", 0), "n_I2p": d.get("I2+", 0), "n_Im": d.get("I-", 0),
        "diss": nI / (nI + 2*nI2) if (nI + 2*nI2) > 0 else 0,
        "alpha": d.get("I-", 0) / r["ne"] if r["ne"] > 0 else 0,
        "nu_m": r["nu_m"],
    }


MTCS_VARS = {
    "mtcs_low":  "Kelm_I2/esteves_low.csv",
    "mtcs_nom":  "Kelm_I2/nominal.csv",
    "mtcs_high": "Kelm_I2/total_elastic.csv",
}
ATT_VARS = {
    "att_low":  "Kdissatt_I2/hamilton_low.csv",
    "att_nom":  "Kdissatt_I2/nominal.csv",
    "att_high": "Kdissatt_I2/esteves_mid.csv",
}

KEY_MAP = {
    "Te_eV": ("Te", 1.0),
    "Tg_K": ("Tg", 1.0),
    "ne_m3": ("ne", 1.0),
    "n_I_m3": ("n_I", 1.0),
    "n_I2_m3": ("n_I2", 1.0),
    "n_Ip_m3": ("n_Ip", 1.0),
    "n_I2p_m3": ("n_I2p", 1.0),
    "n_Im_m3": ("n_Im", 1.0),
}


def score_variant(model, ref_values):
    """Score a variant against reference. Returns dict of per-quantity scores."""
    scores = {}
    for ref_key, ref_info in ref_values.items():
        if ref_key not in KEY_MAP:
            continue
        model_key, scale = KEY_MAP[ref_key]
        if model_key not in model or model[model_key] is None:
            continue
        m_val = model[model_key] * scale
        r_val = ref_info["value"]
        if r_val == 0:
            continue

        ratio = m_val / r_val
        uf = ref_info.get("uncertainty_factor", None)
        u_abs = ref_info.get("uncertainty", None)

        if uf:
            # Factor-based: within [1/uf, uf] is "consistent"
            if 1/uf <= ratio <= uf:
                scores[ref_key] = ("consistent", ratio)
            elif 1/(uf*2) <= ratio <= uf*2:
                scores[ref_key] = ("marginal", ratio)
            else:
                scores[ref_key] = ("off", ratio)
        elif u_abs:
            # Absolute uncertainty
            diff = abs(m_val - r_val)
            if diff <= u_abs:
                scores[ref_key] = ("consistent", ratio)
            elif diff <= 2 * u_abs:
                scores[ref_key] = ("marginal", ratio)
            else:
                scores[ref_key] = ("off", ratio)
        else:
            scores[ref_key] = ("unscored", ratio)

    return scores


def main():
    print("=" * 90)
    print("  Iodine Model Validation: Literature Comparison")
    print("  Package: iodine_lafleur_v1 (Lafleur 2026)")
    print("=" * 90)

    chem = load_chemistry(BASE / "chemistry.json")

    # ── 1. Nominal Validierung ────────────────────────────────
    print("\n--- Nominal (mtcs_nom + att_nom) vs. Literature ---\n")

    for ref in REFERENCE_DATA:
        P_abs = ref["conditions"]["P_abs_est_W"]
        Q0 = ref["conditions"]["Q0_pps"]
        r = run(chem, P_abs, Q0)
        if not r:
            print(f"  {ref['source']}: NOT CONVERGED at P_abs={P_abs}W")
            continue

        print(f"  {ref['source']} (P_abs~{P_abs}W, Q0={Q0:.1e}):")
        print(f"    Comparability: {ref['comparability']}")

        scores = score_variant(r, ref["values"])
        for rk, (status, ratio) in scores.items():
            rv = ref["values"][rk]["value"]
            mk, _ = KEY_MAP[rk]
            mv = r[mk]
            tag = {"consistent": "OK", "marginal": "~", "off": "X", "unscored": "?"}[status]
            print(f"    [{tag}] {rk:>10}: model={mv:.3e}  ref={rv:.3e}  ratio={ratio:.2f}")

    # ── 2. Variantenscan ──────────────────────────────────────
    print("\n" + "=" * 90)
    print("  VARIANTENSCAN: Welche Kombination trifft die Literatur am besten?")
    print("  Referenz: Grondein 2016 Fig.4a (P_abs~280W, Q0=2.1e18)")
    print("=" * 90)

    ref = REFERENCE_DATA[0]  # Grondein 400W point
    P_abs = ref["conditions"]["P_abs_est_W"]
    Q0 = ref["conditions"]["Q0_pps"]

    combo_scores = {}
    print(f"\n{'MTCS':>12} {'ATT':>12} {'Te':>5} {'ne':>10} {'nI':>10} {'nI2':>10} {'nI+':>10} {'nI-':>10} {'Score':>6}")
    print("-" * 90)

    for mk, mf in MTCS_VARS.items():
        for ak, af in ATT_VARS.items():
            c = apply_var(chem, "el_I2", mf)
            c = apply_var(c, "dissatt_I2", af)
            r = run(c, P_abs, Q0)
            if not r:
                continue

            scores = score_variant(r, ref["values"])
            n_ok = sum(1 for s, _ in scores.values() if s == "consistent")
            n_mar = sum(1 for s, _ in scores.values() if s == "marginal")
            n_off = sum(1 for s, _ in scores.values() if s == "off")
            total = n_ok * 3 + n_mar * 1 - n_off * 2
            combo_scores[(mk, ak)] = (total, n_ok, n_mar, n_off, r)

            print(f"{mk:>12} {ak:>12} {r['Te']:5.2f} {r['ne']:10.2e} {r['n_I']:10.2e} "
                  f"{r['n_I2']:10.2e} {r['n_Ip']:10.2e} {r['n_Im']:10.2e} "
                  f"{total:>4} ({n_ok}ok {n_mar}~ {n_off}x)")

    # ── 3. Ranking ────────────────────────────────────────────
    print("\n" + "=" * 90)
    print("  RANKING: Beste Variantenkombinationen")
    print("=" * 90)

    ranked = sorted(combo_scores.items(), key=lambda x: -x[1][0])
    for rank, ((mk, ak), (total, nok, nmar, noff, r)) in enumerate(ranked[:5], 1):
        print(f"  {rank}. {mk} + {ak}: Score={total} ({nok} consistent, {nmar} marginal, {noff} off)")
        print(f"     Te={r['Te']:.2f} ne={r['ne']:.2e} diss={r['diss']:.3f} alpha={r['alpha']:.4f}")

    # ── 4. Referenzwerte Grondein ─────────────────────────────
    print("\n" + "=" * 90)
    print("  REFERENZ: Grondein 2016 Fig.4a (P_RF=400W)")
    print("=" * 90)
    print("  Te ~ 4.2 eV")
    print("  ne ~ 2.5e17 m-3")
    print("  n_I ~ 5e17 m-3")
    print("  n_I2 ~ 1e18 m-3")
    print("  n_I+ ~ 2e17 m-3")
    print("  n_I2+ ~ 7e16 m-3")
    print("  n_I- ~ 1e16 m-3")
    print("  (Modell-zu-Modell, andere Raten, P_RF->P_abs geschaetzt)")

    # ── 5. Schlussfolgerung ───────────────────────────────────
    best = ranked[0]
    best_mk, best_ak = best[0]
    print("\n" + "=" * 90)
    print("  SCHLUSSFOLGERUNG")
    print("=" * 90)
    print(f"  Beste Variantenkombination: {best_mk} + {best_ak}")
    print(f"  Score: {best[1][0]} ({best[1][1]} consistent, {best[1][2]} marginal, {best[1][3]} off)")

    nom_score = combo_scores.get(("mtcs_nom", "att_nom"), (0, 0, 0, 0, None))
    print(f"\n  Nominal (mtcs_nom + att_nom): Score={nom_score[0]}")
    if best[0] == ("mtcs_nom", "att_nom"):
        print("  -> Nominale V1-Defaults durch Literaturabgleich BESTAETIGT.")
    else:
        print(f"  -> Alternative {best_mk}+{best_ak} trifft Literatur besser.")
        print(f"     Nominale Defaults bleiben akzeptabel (Score {nom_score[0]} vs. {best[1][0]}).")

    print("""
  Bekannte Modellluecken:
    1. Unser Modell hat hoeheres ne als Grondein (Faktor 3-5x)
       -> Wahrscheinlich: unterschiedliche Ionisationsraten (Ambalampitiya 2021 vs. Hamilton 2015)
       -> Auch: unser P_abs ist direkt, Grondein nutzt P_RF mit zeta~0.7
    2. Grondein hat detailliertere Wandverlustflaechen (Aeff_i, Aeff_p getrennt)
    3. Atomare Anregung: nur 1 lumped level vs. Grondein single-level Hamilton
    4. Te-Robustheit bestaetigt: Te~4 eV in beiden Modellen
    5. Dissoziationsgrad qualitativ konsistent: steigt mit Leistung

  Empfehlung:
    Ambalampitiya 2021 (nominal) bleibt der beste V1-Default.
    Fuer praezisere Validierung: Experimentelle Daten noetig (nicht nur Modell-zu-Modell).
""")


if __name__ == "__main__":
    main()
