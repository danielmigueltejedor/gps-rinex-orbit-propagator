import math
import re
import csv
from pathlib import Path

# =========================
# Configuración
# =========================

RINEX_FILE = "Modified_RINEX_I2.rnx"
PRN = 17
T0 = 311400.0          # tiempo GPS objetivo [s]
DT = 1800.0            # paso 30 min [s]
HOURS = 8              # intervalo total centrado en T0

OUT_TXT = "SV_PRN17.txt"
OUT_CSV = "PRN17_resultados.csv"

# Constantes GPS
MU = 3.986005e14
OMEGA_E = 7.2921151467e-5

# WGS84
A_WGS84 = 6378137.0
F_WGS84 = 1 / 298.257223563
E2_WGS84 = F_WGS84 * (2 - F_WGS84)


# =========================
# Lectura RINEX
# =========================

def nums(line):
    """Extrae números de una línea RINEX, aceptando exponentes D o E."""
    pattern = r"[-+]?\d+\.\d+(?:[DEde][-+]?\d+)?|[-+]?\d+(?:[DEde][-+]?\d+)?"
    return [float(x.replace("D", "E").replace("d", "E")) for x in re.findall(pattern, line)]


def read_rinex_nav(path):
    """Lee un RINEX de navegación y devuelve efemérides GPS."""
    lines = Path(path).read_text(errors="ignore").splitlines()

    start = next(i for i, l in enumerate(lines) if "END OF HEADER" in l) + 1
    ephs = []
    i = start

    while i < len(lines):
        l0 = lines[i]

        if not l0.strip() or not l0.startswith("G"):
            i += 1
            continue

        block = lines[i:i + 8]
        if len(block) < 8:
            break

        prn = int(l0[1:3])
        clk = nums(l0[23:])
        v = [nums(b) for b in block[1:]]

        if len(clk) < 3 or any(len(row) < 4 for row in v[:4]) or len(v[4]) < 3:
            raise ValueError(f"Bloque RINEX mal leído para {l0[:3]}:\n" + "\n".join(block))

        ephs.append({
            "prn": prn,
            "epoch": l0[4:23].strip(),
            "af0": clk[0],
            "af1": clk[1],
            "af2": clk[2],

            "IODE": v[0][0],
            "Crs": v[0][1],
            "dn": v[0][2],
            "M0": v[0][3],

            "Cuc": v[1][0],
            "e": v[1][1],
            "Cus": v[1][2],
            "sqrtA": v[1][3],

            "toe": v[2][0],
            "Cic": v[2][1],
            "Omega0": v[2][2],
            "Cis": v[2][3],

            "i0": v[3][0],
            "Crc": v[3][1],
            "omega": v[3][2],
            "Omegadot": v[3][3],

            "IDOT": v[4][0],
            "week": v[4][2],

            "SV_accuracy": v[5][0] if len(v[5]) > 0 else None,
            "SV_health": v[5][1] if len(v[5]) > 1 else None,
            "TGD": v[5][2] if len(v[5]) > 2 else None,
            "IODC": v[5][3] if len(v[5]) > 3 else None,
        })

        i += 8

    return ephs


# =========================
# Algoritmo GPS
# =========================

def gps_time_diff(t, toe):
    """Diferencia temporal ajustada a media semana GPS."""
    tk = t - toe
    if tk > 302400:
        tk -= 604800
    if tk < -302400:
        tk += 604800
    return tk


def nearest_eph(ephs, t):
    """Selecciona la efeméride con toe más cercano al tiempo t."""
    return min(ephs, key=lambda e: abs(gps_time_diff(t, e["toe"])))


def kepler(M, e, tol=1e-12):
    """Resuelve E - e sin(E) = M por Newton-Raphson."""
    E = M
    for _ in range(50):
        dE = -(E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E


def satpos(eph, t):
    """Calcula posición ECEF del satélite a partir de efemérides GPS."""
    tk = gps_time_diff(t, eph["toe"])

    A = eph["sqrtA"] ** 2
    n0 = math.sqrt(MU / A**3)
    n = n0 + eph["dn"]

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


def ecef_to_geodetic(x, y, z):
    """Conversión aproximada ECEF -> latitud, longitud y altura WGS84."""
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1 - E2_WGS84))

    for _ in range(10):
        N = A_WGS84 / math.sqrt(1 - E2_WGS84 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N
        lat_new = math.atan2(z, p * (1 - E2_WGS84 * N / (N + h)))
        if abs(lat_new - lat) < 1e-14:
            lat = lat_new
            break
        lat = lat_new

    N = A_WGS84 / math.sqrt(1 - E2_WGS84 * math.sin(lat) ** 2)
    h = p / math.cos(lat) - N

    return math.degrees(lat), math.degrees(lon), h


# =========================
# Programa principal
# =========================

def main():
    ephs_all = read_rinex_nav(RINEX_FILE)
    ephs = [e for e in ephs_all if e["prn"] == PRN]

    if not ephs:
        raise ValueError(f"No hay efemérides para G{PRN:02d}")

    print(f"Efemérides encontradas para G{PRN:02d}: {len(ephs)}")

    # Guardar efemérides extraídas
    with open(OUT_TXT, "w", encoding="utf-8") as f:
        for e in ephs:
            f.write(f"G{PRN:02d} epoch={e['epoch']} toe={e['toe']} week={e['week']}\n")
            f.write(str(e) + "\n\n")

    # Resultado principal
    eph0 = nearest_eph(ephs, T0)
    x, y, z, tk = satpos(eph0, T0)
    lat, lon, h = ecef_to_geodetic(x, y, z)

    print("\nRESULTADO PRINCIPAL")
    print(f"PRN: G{PRN:02d}")
    print(f"t = {T0:.0f} s")
    print(f"toe = {eph0['toe']:.0f} s")
    print(f"tk = {tk:.0f} s")
    print(f"x = {x:.3f} m")
    print(f"y = {y:.3f} m")
    print(f"z = {z:.3f} m")
    print(f"lat = {lat:.6f}º")
    print(f"lon = {lon:.6f}º")
    print(f"h = {h:.3f} m")

    # Tabla de 8 h centrada en T0
    rows = []
    t_ini = T0 - HOURS * 3600 / 2
    t_fin = T0 + HOURS * 3600 / 2

    t = t_ini
    while t <= t_fin:
        e = nearest_eph(ephs, t)
        x, y, z, tk = satpos(e, t)
        lat, lon, h = ecef_to_geodetic(x, y, z)

        rows.append({
            "t_GPS_s": round(t, 3),
            "toe_s": round(e["toe"], 3),
            "tk_s": round(tk, 3),
            "x_ECEF_m": round(x, 3),
            "y_ECEF_m": round(y, 3),
            "z_ECEF_m": round(z, 3),
            "lat_deg": round(lat, 6),
            "lon_deg": round(lon, 6),
            "h_m": round(h, 3),
        })

        t += DT

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nArchivos generados:")
    print(f"- {OUT_TXT}")
    print(f"- {OUT_CSV}")


if __name__ == "__main__":
    main()
