#!/usr/bin/env python3
"""Calcule la moyenne point par point des courbes I(Q) dans chaque sous-dossier *_primus_out."""
import os
import numpy as np
from pathlib import Path

def process_primus_out(directory: Path):
    dat_files = sorted(directory.glob("*.dat"))
    if not dat_files:
        print(f"[SKIP] Aucun .dat dans {directory}")
        return

    q_arrays = []
    i_arrays = []

    for dat_file in dat_files:
        try:
            data = np.loadtxt(dat_file, comments='#')
        except Exception as e:
            print(f"[ERROR] Impossible de lire {dat_file}: {e}")
            continue

        if data.ndim == 1:
            data = data.reshape(1, -1)

        if data.shape[1] < 2:
            print(f"[SKIP] {dat_file}: moins de 2 colonnes")
            continue

        q_arrays.append(data[:, 0])
        i_arrays.append(data[:, 1])

    if not q_arrays:
        print(f"[SKIP] Aucune donnée valide dans {directory}")
        return

    # Vérifier que toutes les colonnes Q sont identiques (même Q pour chaque point)
    q_ref = q_arrays[0]
    for q_arr in q_arrays[1:]:
        if not np.allclose(q_arr, q_ref):
            print(f"[WARN]Les valeurs de Q ne correspondent pas exactement dans {directory}, on utilisera la première.")

    i_mean = np.mean(i_arrays, axis=0)

    output_path = directory / "average_I.dat"
    output_data = np.column_stack((q_ref, i_mean))
    np.savetxt(output_path, output_data, fmt="%.6f\t%.6e")
    print(f"[OK] Moyenne enregistrée dans {output_path} ({len(i_arrays)} fichiers)")


def main():
    import sys
    args = sys.argv[1:]
    if args:
        # explicit target directories (e.g. the D_L_curves / D_LK_curves folders)
        for d in args:
            process_primus_out(Path(d))
        return
    base_dir = Path.cwd()
    primus_dirs = sorted([d for d in base_dir.rglob("*_primus_out") if d.is_dir()])

    if not primus_dirs:
        print("Aucun sous-dossier *_primus_out trouvé.")
        return

    print(f"Trouvé {len(primus_dirs)} dossier(s) *_primus_out.")
    for d in primus_dirs:
        process_primus_out(d)


if __name__ == "__main__":
    main()
