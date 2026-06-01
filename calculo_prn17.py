import math
import csv
from pathlib import Path


# ============================================================
# CONFIGURACIÓN DEL GRUPO
# ============================================================

RINEX_FILE = "Modified_RINEX_I2.rnx"

PRN_OBJETIVO = 17
T_OBJETIVO = 311400.0          # segundos de la semana GPS
PASO = 1800.0                  # 30 minutos
INTERVALO_HORAS = 8            # intervalo total centrado en t

OUTPUT_SV = "SV_PRN17.txt"
OUTPUT_CSV = "PRN17_resultados.csv"


# ============================================================
# CONSTANTES GPS
# ============================================================

MU = 3.986005e14               # m^3/s^2
OMEGA_E = 7.2921151467e-5      # rad/s
PI = math.pi

# WGS84 para conversión ECEF -> geodéticas
WGS84_A = 6378137.0
WGS84_F = 1 / 298.257223563
WGS84_E2 = WGS84_F * (2 - WGS84_F)


# ============================================================
# UTILIDADES
# ============================================================

def rinex_float(value: str) -> float:
    """
    Convierte números en formato RINEX.
    RINEX puede usar D en vez de E para exponentes.
    """
    value = value.replace("D", "E").replace("d", "E").strip()
    return float(value)


def split_rinex_values(line: str):
    """
    Extrae valores numéricos de una línea RINEX de navegación.
    En RINEX 3 los números suelen ocupar campos de 19 caracteres.
    """
    values = []
    for i in range(0, len(line), 19):
        chunk = line[i:i + 19].strip()
        if chunk:
            try:
                values.append(rinex_float(chunk))
            except ValueError:
                pass
    return values


def normalize_time(tk: float) -> float:
    """
    Ajusta tk al rango recomendado por GPS:
    -302400 <= tk <= 302400
    """
    half_week = 302400.0
    week = 604800.0

    if tk > half_week:
        tk -= week
    elif tk < -half_week:
        tk += week

    return tk


# ============================================================
# LECTURA DEL RINEX DE NAVEGACIÓN
# ============================================================

def read_gps_rinex_nav(filename: str):
    """
    Lee un archivo RINEX de navegación y devuelve una lista de efemérides GPS.

    El diccionario generado contiene los parámetros orbitales principales:
    - sqrtA, e, M0, Delta_n, omega, Omega0, i0, Omega_dot, IDOT
    - Cuc, Cus, Crc, Crs, Cic, Cis
    - toe, af0, af1, af2, etc.
    """

    path = Path(filename)

    if not path.exists():
        raise FileNotFoundError(f"No se encontró el archivo: {filename}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    # Saltar cabecera
    start_index = 0
    for i, line in enumerate(lines):
        if "END OF HEADER" in line:
            start_index = i + 1
            break

    ephs = []
    i = start_index

    while i < len(lines):
        line0 = lines[i]

        if not line0.strip():
            i += 1
            continue

        # En RINEX 3 los satélites GPS empiezan por Gxx
        sat_id = line0[0:3].strip()

        if not sat_id.startswith("G"):
            i += 1
            continue

        try:
            prn = int(sat_id[1:])
        except ValueError:
            i += 1
            continue

        # Cada bloque GPS tiene 8 líneas: línea inicial + 7 líneas
        block = lines[i:i + 8]

        if len(block) < 8:
            break

        # Línea 0: fecha + reloj
        # Formato aproximado:
        # G17 2022 12 21 14 00 00 af0 af1 af2
        year = int(line0[4:8])
        month = int(line0[9:11])
        day = int(line0[12:14])
        hour = int(line0[15:17])
        minute = int(line0[18:20])
        second = float(line0[21:23])

        clock_values = split_rinex_values(line0[23:])
        af0, af1, af2 = clock_values[0:3]

        vals1 = split_rinex_values(block[1])
        vals2 = split_rinex_values(block[2])
        vals3 = split_rinex_values(block[3])
        vals4 = split_rinex_values(block[4])
        vals5 = split_rinex_values(block[5])
        vals6 = split_rinex_values(block[6])
        vals7 = split_rinex_values(block[7])

        eph = {
            "prn": prn,

            "year": year,
            "month": month,
            "day": day,
            "hour": hour,
            "minute": minute,
            "second": second,

            "af0": af0,
            "af1": af1,
            "af2": af2,

            # Línea 1
            "IODE": vals1[0],
            "Crs": vals1[1],
            "Delta_n": vals1[2],
            "M0": vals1[3],

            # Línea 2
            "Cuc": vals2[0],
            "e": vals2[1],
            "Cus": vals2[2],
            "sqrtA": vals2[3],

            # Línea 3
            "toe": vals3[0],
            "Cic": vals3[1],
            "Omega0": vals3[2],
            "Cis": vals3[3],

            # Línea 4
            "i0": vals4[0],
            "Crc": vals4[1],
            "omega": vals4[2],
            "Omega_dot": vals4[3],

            # Línea 5
            "IDOT": vals5[0],
            "codes_L2": vals5[1] if len(vals5) > 1 else None,
            "GPS_week": vals5[2] if len(vals5) > 2 else None,
            "L2P_flag": vals5[3] if len(vals5) > 3 else None,

            # Línea 6
            "SV_accuracy": vals6[0] if len(vals6) > 0 else None,
            "SV_health": vals6[1] if len(vals6) > 1 else None,
            "TGD": vals6[2] if len(vals6) > 2 else None,
            "IODC": vals6[3] if len(vals6) > 3 else None,

            # Línea 7
            "transmission_time": vals7[0] if len(vals7) > 0 else None,
            "fit_interval": vals7[1] if len(vals7) > 1 else None,
        }

        ephs.append(eph)
        i += 8

    return ephs


# ============================================================
# SELECCIÓN DE EFEMÉRIDES
# ============================================================

def filter_prn(ephs, prn: int):
    return [eph for eph in ephs if eph["prn"] == prn]


def select_nearest_ephemeris(ephs_prn, t: float):
    """
    Selecciona la efeméride cuyo toe esté más cerca del instante t.
    """
    if not ephs_prn:
        raise ValueError("No hay efemérides para el PRN indicado.")

    return min(ephs_prn, key=lambda eph: abs(normalize_time(t - eph["toe"])))


# ============================================================
# ALGORITMO GPS
# ============================================================

def solve_kepler(M: float, e: float, tol=1e-12, max_iter=50):
    """
    Resuelve la ecuación de Kepler:
        E - e sin(E) = M

    usando Newton-Raphson.
    """
    E = M

    for _ in range(max_iter):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)

        dE = -f / fp
        E += dE

        if abs(dE) < tol:
            break

    return E


def compute_satellite_position(eph, t: float):
    """
    Calcula la posición ECEF del satélite GPS en el instante t.

    Devuelve:
    - x, y, z en metros
    - diccionario con variables intermedias
    """

    toe = eph["toe"]
    tk = normalize_time(t - toe)

    # Semieje mayor
    A = eph["sqrtA"] ** 2

    # Movimiento medio calculado
    n0 = math.sqrt(MU / A**3)

    # Movimiento medio corregido
    n = n0 + eph["Delta_n"]

    # Anomalía media
    Mk = eph["M0"] + n * tk

    # Anomalía excéntrica
    Ek = solve_kepler(Mk, eph["e"])

    # Anomalía verdadera
    vk = math.atan2(
        math.sqrt(1.0 - eph["e"]**2) * math.sin(Ek),
        math.cos(Ek) - eph["e"]
    )

    # Argumento de latitud sin corregir
    phik = vk + eph["omega"]

    # Correcciones armónicas
    du = eph["Cuc"] * math.cos(2.0 * phik) + eph["Cus"] * math.sin(2.0 * phik)
    dr = eph["Crc"] * math.cos(2.0 * phik) + eph["Crs"] * math.sin(2.0 * phik)
    di = eph["Cic"] * math.cos(2.0 * phik) + eph["Cis"] * math.sin(2.0 * phik)

    # Parámetros corregidos
    u = phik + du
    r = A * (1.0 - eph["e"] * math.cos(Ek)) + dr
    i = eph["i0"] + eph["IDOT"] * tk + di

    # Coordenadas en el plano orbital
    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    # Longitud corregida del nodo ascendente
    Omega_k = (
        eph["Omega0"]
        + (eph["Omega_dot"] - OMEGA_E) * tk
        - OMEGA_E * toe
    )

    # Transformación a ECEF
    x = x_orb * math.cos(Omega_k) - y_orb * math.cos(i) * math.sin(Omega_k)
    y = x_orb * math.sin(Omega_k) + y_orb * math.cos(i) * math.cos(Omega_k)
    z = y_orb * math.sin(i)

    debug = {
        "toe": toe,
        "tk": tk,
        "A": A,
        "n0": n0,
        "n": n,
        "Mk": Mk,
        "Ek": Ek,
        "vk": vk,
        "phik": phik,
        "du": du,
        "dr": dr,
        "di": di,
        "u": u,
        "r": r,
        "i": i,
        "x_orb": x_orb,
        "y_orb": y_orb,
        "Omega_k": Omega_k,
    }

    return x, y, z, debug


# ============================================================
# CONVERSIÓN ECEF -> GEODÉTICAS
# ============================================================

def ecef_to_geodetic(x: float, y: float, z: float):
    """
    Convierte coordenadas ECEF a latitud, longitud y altura usando WGS84.

    Devuelve:
    - latitud en grados
    - longitud en grados
    - altura en metros
    """

    lon = math.atan2(y, x)

    p = math.sqrt(x**2 + y**2)
    lat = math.atan2(z, p * (1.0 - WGS84_E2))

    for _ in range(10):
        N = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat)**2)
        h = p / math.cos(lat) - N
        lat_new = math.atan2(z, p * (1.0 - WGS84_E2 * N / (N + h)))

        if abs(lat_new - lat) < 1e-14:
            lat = lat_new
            break

        lat = lat_new

    N = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat)**2)
    h = p / math.cos(lat) - N

    lat_deg = math.degrees(lat)
    lon_deg = math.degrees(lon)

    return lat_deg, lon_deg, h


# ============================================================
# SALIDAS
# ============================================================

def write_sv_file(ephs_prn, filename: str):
    """
    Guarda en TXT las efemérides extraídas del PRN objetivo.
    """

    with open(filename, "w", encoding="utf-8") as f:
        f.write(f"Efemérides extraídas para PRN{PRN_OBJETIVO:02d}\n")
        f.write("=" * 70 + "\n\n")

        for eph in ephs_prn:
            f.write(
                f"PRN{eph['prn']:02d} "
                f"{eph['year']:04d}-{eph['month']:02d}-{eph['day']:02d} "
                f"{eph['hour']:02d}:{eph['minute']:02d}:{int(eph['second']):02d} "
                f"toe={eph['toe']:.3f} "
                f"GPS_week={eph['GPS_week']}\n"
            )

            campos = [
                "af0", "af1", "af2",
                "IODE", "Crs", "Delta_n", "M0",
                "Cuc", "e", "Cus", "sqrtA",
                "toe", "Cic", "Omega0", "Cis",
                "i0", "Crc", "omega", "Omega_dot",
                "IDOT", "GPS_week", "SV_accuracy",
                "SV_health", "TGD", "IODC",
            ]

            for campo in campos:
                f.write(f"  {campo:12s}: {eph[campo]}\n")

            f.write("\n")


def write_results_csv(rows, filename: str):
    """
    Guarda la tabla de resultados en CSV.
    """

    fieldnames = [
        "t_GPS_s",
        "toe_usado_s",
        "tk_s",
        "x_ECEF_m",
        "y_ECEF_m",
        "z_ECEF_m",
        "lat_deg",
        "lon_deg",
        "h_m",
    ]

    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def print_main_result(eph, x, y, z, lat, lon, h, debug):
    """
    Muestra por pantalla el resultado principal.
    """

    print()
    print("=" * 70)
    print("RESULTADO PRINCIPAL")
    print("=" * 70)
    print(f"PRN objetivo:       PRN{PRN_OBJETIVO:02d}")
    print(f"t objetivo:         {T_OBJETIVO:.3f} s")
    print(f"toe seleccionado:   {eph['toe']:.3f} s")
    print(f"tk:                 {debug['tk']:.3f} s")
    print()
    print("Coordenadas ECEF:")
    print(f"x = {x: .3f} m")
    print(f"y = {y: .3f} m")
    print(f"z = {z: .3f} m")
    print()
    print("Coordenadas geodéticas aproximadas:")
    print(f"lat = {lat: .6f} grados")
    print(f"lon = {lon: .6f} grados")
    print(f"h   = {h: .3f} m")
    print("=" * 70)
    print()


# ============================================================
# PROGRAMA PRINCIPAL
# ============================================================

def main():
    # 1. Leer RINEX
    ephs = read_gps_rinex_nav(RINEX_FILE)

    # 2. Filtrar PRN objetivo
    ephs_prn = filter_prn(ephs, PRN_OBJETIVO)

    if not ephs_prn:
        raise ValueError(f"No se encontraron efemérides para PRN{PRN_OBJETIVO:02d}")

    # 3. Guardar archivo SV_PRN17.txt
    write_sv_file(ephs_prn, OUTPUT_SV)

    # 4. Resultado para el instante objetivo
    eph_main = select_nearest_ephemeris(ephs_prn, T_OBJETIVO)
    x, y, z, debug = compute_satellite_position(eph_main, T_OBJETIVO)
    lat, lon, h = ecef_to_geodetic(x, y, z)

    print_main_result(eph_main, x, y, z, lat, lon, h, debug)

    # 5. Tabla para 8 horas centradas en t
    half_interval = INTERVALO_HORAS * 3600.0 / 2.0
    t_start = T_OBJETIVO - half_interval
    t_end = T_OBJETIVO + half_interval

    rows = []

    t = t_start
    while t <= t_end + 1e-9:
        eph = select_nearest_ephemeris(ephs_prn, t)
        x, y, z, dbg = compute_satellite_position(eph, t)
        lat, lon, h = ecef_to_geodetic(x, y, z)

        rows.append({
            "t_GPS_s": round(t, 3),
            "toe_usado_s": round(eph["toe"], 3),
            "tk_s": round(dbg["tk"], 3),
            "x_ECEF_m": round(x, 3),
            "y_ECEF_m": round(y, 3),
            "z_ECEF_m": round(z, 3),
            "lat_deg": round(lat, 6),
            "lon_deg": round(lon, 6),
            "h_m": round(h, 3),
        })

        t += PASO

    write_results_csv(rows, OUTPUT_CSV)

    # 6. Imprimir tabla
    print("TABLA DE RESULTADOS")
    print("-" * 120)
    print(
        f"{'t GPS (s)':>10} "
        f"{'toe (s)':>10} "
        f"{'tk (s)':>10} "
        f"{'x (m)':>15} "
        f"{'y (m)':>15} "
        f"{'z (m)':>15} "
        f"{'lat (deg)':>12} "
        f"{'lon (deg)':>12}"
    )
    print("-" * 120)

    for row in rows:
        print(
            f"{row['t_GPS_s']:10.0f} "
            f"{row['toe_usado_s']:10.0f} "
            f"{row['tk_s']:10.0f} "
            f"{row['x_ECEF_m']:15.3f} "
            f"{row['y_ECEF_m']:15.3f} "
            f"{row['z_ECEF_m']:15.3f} "
            f"{row['lat_deg']:12.6f} "
            f"{row['lon_deg']:12.6f}"
        )

    print("-" * 120)
    print(f"\nArchivo generado: {OUTPUT_SV}")
    print(f"Archivo generado: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
