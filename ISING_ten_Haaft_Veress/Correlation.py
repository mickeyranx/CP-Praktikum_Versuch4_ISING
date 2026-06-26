import numpy as np
import puwr as uwerr
from pathlib import Path

# Einstellungen
L = 32
beta = 0.600
mh_values = [1, 2, 3, 4, 5]

input_dir = Path("project/output_files/testing_therm_params")
output_file = input_dir / f"uwerr_results_L={L}_beta={beta:.3f}.txt"

# Ergebniszeilen sammeln
results = []

# Header für Ausgabedatei
header = (
    "MH "
    "mean_E err_E tau_E dtau_E cost_E "
    "mean_M err_M tau_M dtau_M cost_M "
    "mean_absM err_absM tau_absM dtau_absM cost_absM "
    "mean_t_sweep\n"
)

for mh in mh_values:
    filename = input_dir / f"testing_metro_L={L}_beta={beta:.3f}_MH={mh}_start=true.txt"

    if not filename.exists():
        print(f"Datei nicht gefunden: {filename}")
        continue

    data = np.loadtxt(filename, skiprows=1)

    # Spalten aus deiner Datei
    # Annahme:
    # data[:, 1] = E
    # data[:, 2] = M
    # data[:, 3] = |M|
    # data[:, 4] = t_sweep
    E = data[:, 1]
    M = data[:, 2]
    absM = data[:, 3]
    t_sweeps = data[:, 4]

    mean_t_sweep = np.mean(t_sweeps)

    # py-uwerr Format:
    # data[observable][replica][measurement]
    data_uwerr = [
        [E],
        [M],
        [absM],
    ]

    mean_E, err_E, tau_E, dtau_E = uwerr.tauint(data_uwerr, 0)
    mean_M, err_M, tau_M, dtau_M = uwerr.tauint(data_uwerr, 1)
    mean_absM, err_absM, tau_absM, dtau_absM = uwerr.tauint(data_uwerr, 2)

    cost_E = tau_E * mean_t_sweep
    cost_M = tau_M * mean_t_sweep
    cost_absM = tau_absM * mean_t_sweep

    results.append([
        mh,
        mean_E, err_E, tau_E, dtau_E, cost_E,
        mean_M, err_M, tau_M, dtau_M, cost_M,
        mean_absM, err_absM, tau_absM, dtau_absM, cost_absM,
        mean_t_sweep,
    ])

    print(f"\nMH = {mh}")
    print(f"E:    mean = {mean_E:.8f}, err = {err_E:.8f}, tau = {tau_E:.4f}, cost = {cost_E:.8e}")
    print(f"M:    mean = {mean_M:.8f}, err = {err_M:.8f}, tau = {tau_M:.4f}, cost = {cost_M:.8e}")
    print(f"|M|:  mean = {mean_absM:.8f}, err = {err_absM:.8f}, tau = {tau_absM:.4f}, cost = {cost_absM:.8e}")
    print(f"mean t_sweep = {mean_t_sweep:.8e}")

# Ergebnisse speichern
with open(output_file, "w") as f:
    f.write(header)

    for row in results:
        f.write(
            f"{int(row[0])} "
            f"{row[1]:.12e} {row[2]:.12e} {row[3]:.12e} {row[4]:.12e} {row[5]:.12e} "
            f"{row[6]:.12e} {row[7]:.12e} {row[8]:.12e} {row[9]:.12e} {row[10]:.12e} "
            f"{row[11]:.12e} {row[12]:.12e} {row[13]:.12e} {row[14]:.12e} {row[15]:.12e} "
            f"{row[16]:.12e}\n"
        )

print(f"\nErgebnisse gespeichert in:\n{output_file}")