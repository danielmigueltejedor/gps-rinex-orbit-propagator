import math
import re
import csv
from pathlib import Path
import matplotlib.pyplot as plt

# Caso del grupo
RINEX_FILE = "Modified_RINEX_I2.rnx"
PRN = 17
T0 = 311400.0
DT = 1800.0
HOURS = 8

OUT_SV = "SV_PRN17.txt"
OUT_CSV = "PRN17_resultados.csv"
OUT_PNG = "PRN17_ground_track.png"

# Constantes GPS y WGS84
MU = 3.986005e14
OMEGA_E = 7.2921151467e-5

A_WGS84 = 6378137.0
F_WGS84 = 1 / 298.257223563
E2_WGS84 = F_WGS84 * (2 - F_WGS84)

SV_COLUMNS = [
    "week", "toe", "sqrtA", "e", "M0", "dn", "omega",
    "Omega0", "i0", "Omegadot", "IDOT",
    "Cuc", "Cus", "Crc", "Crs", "Cic", "Cis",
    "af0", "af1", "af2", "IODE", "IODC", "TGD",
    "SV_health", "SV_accuracy"
]


def nums(line):
    """Extrae números de una línea del RINEX."""
    pattern = r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?|[-+]?\d+(?:[Ee][-+]?\d+)?"
    return [float(x) for x in re.findall(pattern, line)]


def read_rinex_nav(path):
    """Lee el RINEX y devuelve las efemérides GPS."""
    lines = Path(path).read_text(errors="ignore").splitlines()
    start = next(i for i, line in enumerate(lines) if "END OF HEADER" in line) + 1

    ephs = []
    i = start

    while i < len(lines):
        line0 = lines[i]

        if not line0.strip() or not line0.startswith("G"):
            i += 1
            continue

        block = lines[i:i + 8]
        prn = int(line0[1:3])
        clock = nums(line0[23:])
        r = [nums(line) for line in block[1:]]

        ephs.append({
            "prn": prn,
            "epoch": line0[4:23].strip(),
            "af0": clock[0], "af1": clock[1], "af2": clock[2],
            "IODE": r[0][0], "Crs": r[0][1], "dn": r[0][2], "M0": r[0][3],
            "Cuc": r[1][0], "e": r[1][1], "Cus": r[1][2], "sqrtA": r[1][3],
            "toe": r[2][0], "Cic": r[2][1], "Omega0": r[2][2], "Cis": r[2][3],
            "i0": r[3][0], "Crc": r[3][1], "omega": r[3][2], "Omegadot": r[3][3],
            "IDOT": r[4][0], "week": r[4][2],
            "SV_accuracy": r[5][0] if len(r[5]) > 0 else "",
            "SV_health": r[5][1] if len(r[5]) > 1 else "",
            "TGD": r[5][2] if len(r[5]) > 2 else "",
            "IODC": r[5][3] if len(r[5]) > 3 else "",
        })

        i += 8

    return ephs


def write_sv_txt(ephs, filename):
    """Genera SV_PRN17.txt como tabla."""
    with open(filename, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SV_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows({col: eph.get(col, "") for col in SV_COLUMNS} for eph in ephs)


def read_sv_txt(filename):
    """Lee el fichero SV_PRN17.txt."""
    with open(filename, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [
            {key: None if value == "" else float(value) for key, value in row.items()}
            for row in reader
        ]


def gps_time_diff(t, toe):
    """Calcula tk ajustando el cambio de semana GPS si hace falta."""
    tk = t - toe
    if tk > 302400:
        tk -= 604800
    elif tk < -302400:
        tk += 604800
    return tk


def nearest_eph(ephs, t):
    """Selecciona la efeméride con menor |t - toe|."""
    return min(ephs, key=lambda eph: abs(gps_time_diff(t, eph["toe"])))


def select_from_sv_file(filename, t):
    """Función pedida: entrada SV_PRN17.txt y t; salida toe, fila y tk."""
    eph = nearest_eph(read_sv_txt(filename), t)
    return eph["toe"], eph, gps_time_diff(t, eph["toe"])


def kepler(M, e, tol=1e-12, max_iter=50):
    """Resuelve E - e sin(E) = M por Newton-Raphson."""
    E = M

    for _ in range(max_iter):
        dE = -(E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E += dE

        if abs(dE) < tol:
            break

    return E


def satpos(eph, t):
    """Calcula la posición ECEF del satélite."""
    tk = gps_time_diff(t, eph["toe"])

    A = eph["sqrtA"] ** 2
    n = math.sqrt(MU / A**3) + eph["dn"]
    M = eph["M0"] + n * tk
    E = kepler(M, eph["e"])

    nu = math.atan2(
        math.sqrt(1 - eph["e"]**2) * math.sin(E),
        math.cos(E) - eph["e"]
    )

    phi = nu + eph["omega"]

    du = eph["Cuc"] * math.cos(2 * phi) + eph["Cus"] * math.sin(2 * phi)
    dr = eph["Crc"] * math.cos(2 * phi) + eph["Crs"] * math.sin(2 * phi)
    di = eph["Cic"] * math.cos(2 * phi) + eph["Cis"] * math.sin(2 * phi)

    u = phi + du
    r = A * (1 - eph["e"] * math.cos(E)) + dr
    inc = eph["i0"] + eph["IDOT"] * tk + di

    xp = r * math.cos(u)
    yp = r * math.sin(u)

    Omega = eph["Omega0"] + (eph["Omegadot"] - OMEGA_E) * tk - OMEGA_E * eph["toe"]

    x = xp * math.cos(Omega) - yp * math.cos(inc) * math.sin(Omega)
    y = xp * math.sin(Omega) + yp * math.cos(inc) * math.cos(Omega)
    z = yp * math.sin(inc)

    return x, y, z, tk


def ecef_to_geodetic(x, y, z, tol=1e-14):
    """Convierte ECEF a latitud, longitud y altura WGS84."""
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1 - E2_WGS84))

    for _ in range(10):
        N = A_WGS84 / math.sqrt(1 - E2_WGS84 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N
        new_lat = math.atan2(z, p * (1 - E2_WGS84 * N / (N + h)))

        if abs(new_lat - lat) < tol:
            lat = new_lat
            break

        lat = new_lat

    N = A_WGS84 / math.sqrt(1 - E2_WGS84 * math.sin(lat) ** 2)
    h = p / math.cos(lat) - N

    return math.degrees(lat), math.degrees(lon), h


def validate_kepler():
    """Validación de la función de Kepler."""
    M, e = 1.0, 0.01
    E = kepler(M, e)
    residual = E - e * math.sin(E) - M

    print("\nVALIDACIÓN KEPLER")
    print(f"M = {M}")
    print(f"e = {e}")
    print(f"E = {E:.12f}")
    print(f"residual = {residual:.3e}")


def validate_geodetic():
    """Validación ECEF -> geodésicas con un punto conocido."""
    lat, lon, h = ecef_to_geodetic(A_WGS84, 0.0, 0.0)

    print("\nVALIDACIÓN ECEF -> GEODÉSICAS")
    print(f"lat = {lat:.12f} deg")
    print(f"lon = {lon:.12f} deg")
    print(f"h = {h:.6f} m")


def compute_interval(ephs):
    """Calcula posiciones cada 30 minutos durante 8 horas."""
    rows = []
    t = T0 - HOURS * 3600 / 2
    t_end = T0 + HOURS * 3600 / 2

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

        t += DT

    return rows


def write_results_csv(rows):
    """Guarda la tabla de resultados."""
    with open(OUT_CSV, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def plot_ground_track(rows):
    """Representa la traza terrestre."""
    lons = [row["lon_deg"] for row in rows]
    lats = [row["lat_deg"] for row in rows]

    plt.figure(figsize=(10, 5))
    plt.plot(lons, lats, marker="o")
    plt.scatter(lons[len(lons)//2], lats[len(lats)//2], marker="x", s=80)

    plt.title(f"Ground track G{PRN:02d} - intervalo de {HOURS} h")
    plt.xlabel("Longitud [deg]")
    plt.ylabel("Latitud [deg]")
    plt.xlim(-180, 180)
    plt.ylim(-90, 90)
    plt.grid(True)

    plt.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    # 1) Leo el RINEX y filtro el satélite del grupo
    ephs = [e for e in read_rinex_nav(RINEX_FILE) if e["prn"] == PRN]

    if not ephs:
        raise ValueError(f"No hay efemérides para G{PRN:02d}")

    print(f"Efemérides encontradas para G{PRN:02d}: {len(ephs)}")

    # 2) Genero el fichero pedido SV_PRN17.txt
    write_sv_txt(ephs, OUT_SV)

    # 3) Selecciono la efeméride usando el propio fichero generado
    toe, eph0, tk = select_from_sv_file(OUT_SV, T0)

    # 4) Resultado principal
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

    # 5) Validaciones
    validate_kepler()
    validate_geodetic()

    # 6) Tabla del intervalo y traza
    rows = compute_interval(read_sv_txt(OUT_SV))
    write_results_csv(rows)
    plot_ground_track(rows)

    print("\nArchivos generados:")
    print(f"- {OUT_SV}")
    print(f"- {OUT_CSV}")
    print(f"- {OUT_PNG}")


if __name__ == "__main__":
    main()
