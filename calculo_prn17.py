import math
import re
import csv
from pathlib import Path

import matplotlib.pyplot as plt


# ============================================================
# Datos del grupo
# Último dígito suma DNI = 0 -> G17, t = 311400 s
# ============================================================

RINEX_FILE = "Modified_RINEX_I2.rnx"

PRN = 17
T0 = 311400.0          # tiempo GPS pedido [s]
DT = 1800.0            # paso de 30 min [s]
HOURS = 8              # intervalo total para la traza

OUT_SV = "SV_PRN17.txt"
OUT_CSV = "PRN17_resultados.csv"
OUT_PNG = "PRN17_ground_track.png"


# ============================================================
# Constantes
# ============================================================

MU = 3.986005e14              # constante gravitacional GPS [m^3/s^2]
OMEGA_E = 7.2921151467e-5     # velocidad angular Tierra [rad/s]

# Elipsoide WGS84
A_WGS84 = 6378137.0
F_WGS84 = 1 / 298.257223563
E2_WGS84 = F_WGS84 * (2 - F_WGS84)


# Columnas que pide el enunciado para el fichero SV_PRNnum.txt.
# Las dos primeras son week y toe.
SV_COLUMNS = [
    "week", "toe",
    "sqrtA", "e", "M0", "dn", "omega",
    "Omega0", "i0", "Omegadot", "IDOT",
    "Cuc", "Cus", "Crc", "Crs", "Cic", "Cis",
    "af0", "af1", "af2",
    "IODE", "IODC", "TGD", "SV_health", "SV_accuracy"
]


# ============================================================
# Lectura del RINEX
# ============================================================

def nums(line):
    pattern = r"[-+]?\d+\.\d+(?:[DEde][-+]?\d+)?|[-+]?\d+(?:[DEde][-+]?\d+)?"
    return [float(x.replace("D", "E").replace("d", "E")) for x in re.findall(pattern, line)]


def read_rinex_nav(path):
    lines = Path(path).read_text(errors="ignore").splitlines()

    start = next(i for i, line in enumerate(lines) if "END OF HEADER" in line) + 1

    ephs = []
    i = start

    while i < len(lines):
        line0 = lines[i]

        # Me quedo solo con satélites GPS: G01, G02, ...
        if not line0.strip() or not line0.startswith("G"):
            i += 1
            continue

        block = lines[i:i + 8]
        if len(block) < 8:
            break

        prn = int(line0[1:3])
        clock = nums(line0[23:])
        rows = [nums(line) for line in block[1:]]

        ephs.append({
            "prn": prn,
            "epoch": line0[4:23].strip(),

            # Parámetros de reloj
            "af0": clock[0],
            "af1": clock[1],
            "af2": clock[2],

            # Línea 1 del bloque
            "IODE": rows[0][0],
            "Crs": rows[0][1],
            "dn": rows[0][2],
            "M0": rows[0][3],

            # Línea 2
            "Cuc": rows[1][0],
            "e": rows[1][1],
            "Cus": rows[1][2],
            "sqrtA": rows[1][3],

            # Línea 3
            "toe": rows[2][0],
            "Cic": rows[2][1],
            "Omega0": rows[2][2],
            "Cis": rows[2][3],

            # Línea 4
            "i0": rows[3][0],
            "Crc": rows[3][1],
            "omega": rows[3][2],
            "Omegadot": rows[3][3],

            # Línea 5
            "IDOT": rows[4][0],
            "week": rows[4][2],

            # Línea 6
            "SV_accuracy": rows[5][0] if len(rows[5]) > 0 else "",
            "SV_health": rows[5][1] if len(rows[5]) > 1 else "",
            "TGD": rows[5][2] if len(rows[5]) > 2 else "",
            "IODC": rows[5][3] if len(rows[5]) > 3 else "",
        })

        i += 8

    return ephs


def write_sv_txt(ephs, filename):
    with open(filename, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SV_COLUMNS, delimiter="\t")
        writer.writeheader()

        for eph in ephs:
            writer.writerow({col: eph.get(col, "") for col in SV_COLUMNS})


def read_sv_txt(filename):
    ephs = []

    with open(filename, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            eph = {}
            for key, value in row.items():
                eph[key] = None if value == "" else float(value)
            ephs.append(eph)

    return ephs


# ============================================================
# Selección de efemérides
# ============================================================

def gps_time_diff(t, toe):
    tk = t - toe

    if tk > 302400:
        tk -= 604800
    elif tk < -302400:
        tk += 604800

    return tk


def nearest_eph(ephs, t):
    return min(ephs, key=lambda eph: abs(gps_time_diff(t, eph["toe"])))


def select_from_sv_file(filename, t):
    ephs = read_sv_txt(filename)
    eph = nearest_eph(ephs, t)
    tk = gps_time_diff(t, eph["toe"])

    return eph["toe"], eph, tk


# ============================================================
# Kepler
# ============================================================

def kepler(M, e, tol=1e-12, max_iter=50):
    E = M

    for _ in range(max_iter):
        f = E - e * math.sin(E) - M
        fp = 1 - e * math.cos(E)

        dE = -f / fp
        E += dE

        if abs(dE) < tol:
            break

    return E


def validate_kepler():
    M = 1.0
    e = 0.01
    E = kepler(M, e)

    residual = E - e * math.sin(E) - M

    print("\nVALIDACIÓN KEPLER")
    print(f"M = {M}")
    print(f"e = {e}")
    print(f"E = {E:.12f}")
    print(f"residual = {residual:.3e}")


# ============================================================
# Algoritmo GPS: efemérides -> ECEF
# ============================================================

def satpos(eph, t):
    tk = gps_time_diff(t, eph["toe"])

    # Semieje mayor y movimiento medio
    A = eph["sqrtA"] ** 2
    n0 = math.sqrt(MU / A**3)
    n = n0 + eph["dn"]

    # Anomalía media y anomalía excéntrica
    M = eph["M0"] + n * tk
    E = kepler(M, eph["e"])

    # Anomalía verdadera
    nu = math.atan2(
        math.sqrt(1 - eph["e"]**2) * math.sin(E),
        math.cos(E) - eph["e"]
    )

    # Argumento de latitud sin corregir
    phi = nu + eph["omega"]

    # Correcciones armónicas
    du = eph["Cuc"] * math.cos(2 * phi) + eph["Cus"] * math.sin(2 * phi)
    dr = eph["Crc"] * math.cos(2 * phi) + eph["Crs"] * math.sin(2 * phi)
    di = eph["Cic"] * math.cos(2 * phi) + eph["Cis"] * math.sin(2 * phi)

    # Argumento de latitud, radio e inclinación corregidos
    u = phi + du
    r = A * (1 - eph["e"] * math.cos(E)) + dr
    inc = eph["i0"] + eph["IDOT"] * tk + di

    # Posición en el plano orbital
    xp = r * math.cos(u)
    yp = r * math.sin(u)

    # Longitud corregida del nodo ascendente
    Omega = eph["Omega0"] + (eph["Omegadot"] - OMEGA_E) * tk - OMEGA_E * eph["toe"]

    # Transformación final a ECEF
    x = xp * math.cos(Omega) - yp * math.cos(inc) * math.sin(Omega)
    y = xp * math.sin(Omega) + yp * math.cos(inc) * math.cos(Omega)
    z = yp * math.sin(inc)

    return x, y, z, tk


# ============================================================
# ECEF -> geodésicas WGS84
# ============================================================

def ecef_to_geodetic(x, y, z, tol=1e-14, max_iter=10):
    lon = math.atan2(y, x)
    p = math.hypot(x, y)

    # Estimación inicial de latitud
    lat = math.atan2(z, p * (1 - E2_WGS84))

    for _ in range(max_iter):
        N = A_WGS84 / math.sqrt(1 - E2_WGS84 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N

        lat_new = math.atan2(z, p * (1 - E2_WGS84 * N / (N + h)))

        if abs(lat_new - lat) < tol:
            lat = lat_new
            break

        lat = lat_new

    N = A_WGS84 / math.sqrt(1 - E2_WGS84 * math.sin(lat) ** 2)
    h = p / math.cos(lat) - N

    return math.degrees(lat), math.degrees(lon), h


def validate_geodetic():
    lat, lon, h = ecef_to_geodetic(A_WGS84, 0.0, 0.0)

    print("\nVALIDACIÓN ECEF -> GEODÉSICAS")
    print(f"lat = {lat:.12f} deg")
    print(f"lon = {lon:.12f} deg")
    print(f"h = {h:.6f} m")


# ============================================================
# Resultados y traza
# ============================================================

def compute_interval(ephs, t0, hours, dt):
    rows = []

    t = t0 - hours * 3600 / 2
    t_end = t0 + hours * 3600 / 2

    while t <= t_end + 1e-9:
        eph = nearest_eph(ephs, t)

        x, y, z, tk = satpos(eph, t)
        lat, lon, h = ecef_to_geodetic(x, y, z)

        rows.append({
            "t_GPS_s": round(t, 3),
            "toe_s": round(eph["toe"], 3),
            "tk_s": round(tk, 3),
            "x_ECEF_m": round(x, 3),
            "y_ECEF_m": round(y, 3),
            "z_ECEF_m": round(z, 3),
            "lat_deg": round(lat, 6),
            "lon_deg": round(lon, 6),
            "h_m": round(h, 3),
        })

        t += dt

    return rows


def write_results_csv(rows, filename):
    with open(filename, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def plot_ground_track(rows, filename):
    lons = [row["lon_deg"] for row in rows]
    lats = [row["lat_deg"] for row in rows]

    plt.figure(figsize=(10, 5))
    plt.plot(lons, lats, marker="o")
    plt.scatter([lons[len(lons)//2]], [lats[len(lats)//2]], marker="x", s=80)

    plt.title(f"Ground track G{PRN:02d} - intervalo de {HOURS} h")
    plt.xlabel("Longitud [deg]")
    plt.ylabel("Latitud [deg]")

    plt.xlim(-180, 180)
    plt.ylim(-90, 90)
    plt.grid(True)

    plt.savefig(filename, dpi=200, bbox_inches="tight")
    plt.close()


# ============================================================
# Programa principal
# ============================================================

def main():
    # 1) Leo RINEX y filtro mi satélite
    ephs_all = read_rinex_nav(RINEX_FILE)
    ephs_prn = [eph for eph in ephs_all if eph["prn"] == PRN]

    if not ephs_prn:
        raise ValueError(f"No hay efemérides para G{PRN:02d}")

    print(f"Efemérides encontradas para G{PRN:02d}: {len(ephs_prn)}")

    # 2) Genero SV_PRN17.txt como tabla
    write_sv_txt(ephs_prn, OUT_SV)

    # 3) Uso el propio fichero SV_PRN17.txt para seleccionar la efeméride
    toe, eph0, tk = select_from_sv_file(OUT_SV, T0)

    # 4) Calculo posición principal
    x, y, z, tk = satpos(eph0, T0)
    lat, lon, h = ecef_to_geodetic(x, y, z)

    print("\nRESULTADO PRINCIPAL")
    print(f"PRN = G{PRN:02d}")
    print(f"t = {T0:.0f} s")
    print(f"toe = {toe:.0f} s")
    print(f"tk = {tk:.0f} s")
    print(f"x = {x:.3f} m")
    print(f"y = {y:.3f} m")
    print(f"z = {z:.3f} m")
    print(f"lat = {lat:.6f} deg")
    print(f"lon = {lon:.6f} deg")
    print(f"h = {h:.3f} m")

    # 5) Validaciones pedidas
    validate_kepler()
    validate_geodetic()

    # 6) Tabla y ground track de 8 horas
    ephs_from_file = read_sv_txt(OUT_SV)
    rows = compute_interval(ephs_from_file, T0, HOURS, DT)

    write_results_csv(rows, OUT_CSV)
    plot_ground_track(rows, OUT_PNG)

    print("\nArchivos generados:")
    print(f"- {OUT_SV}")
    print(f"- {OUT_CSV}")
    print(f"- {OUT_PNG}")


if __name__ == "__main__":
    main()
