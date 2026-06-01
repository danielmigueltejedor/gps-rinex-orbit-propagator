# GPS RINEX Orbit Propagator

Python tool for computing the ECEF position of a GPS satellite from broadcast ephemerides contained in a RINEX navigation file.

This project was developed as part of an academic assignment focused on GPS satellite positioning using RINEX navigation data and the algorithm described in the GPS interface specification.

## Overview

The program reads a RINEX navigation file, extracts the ephemerides of a selected GPS satellite, chooses the ephemeris closest to the desired GPS time, and computes the satellite position in the Earth-Centered Earth-Fixed coordinate system.

For the current assignment case, the selected configuration is:

```text
Satellite: PRN17
GPS time: 311400 s
Group: 0
```

The program also generates a table of satellite positions every 30 minutes over an 8-hour interval centered on the target time.

## Features

- Reads GPS navigation data from a RINEX file.
- Filters broadcast ephemerides by PRN.
- Selects the closest ephemeris using the `toe` parameter.
- Solves Kepler's equation using the Newton-Raphson method.
- Applies GPS broadcast ephemeris corrections.
- Computes satellite coordinates in the ECEF reference frame.
- Converts ECEF coordinates to approximate geodetic latitude, longitude and height.
- Exports extracted ephemerides to a text file.
- Exports computed results to a CSV file.

## Repository Structure

```text
.
├── calculo_prn17.py
├── Modified_RINEX_I2.rnx
├── SV_PRN17.txt
├── PRN17_resultados.csv
└── README.md
```

## Requirements

The project only uses the Python standard library.

No external dependencies are required.

Tested with:

```text
Python 3.10+
```

## Usage

Place the RINEX navigation file in the same directory as the Python script.

The default expected input file is:

```text
Modified_RINEX_I2.rnx
```

Run the program with:

```bash
python calculo_prn17.py
```

## Configuration

The main configuration parameters are defined at the beginning of the script:

```python
RINEX_FILE = "Modified_RINEX_I2.rnx"

PRN_OBJETIVO = 17
T_OBJETIVO = 311400.0
PASO = 1800.0
INTERVALO_HORAS = 8
```

These values can be modified to compute the position of another satellite or another GPS time.

## Output Files

The program generates two output files:

### `SV_PRN17.txt`

Contains all extracted broadcast ephemerides for the selected satellite.

### `PRN17_resultados.csv`

Contains the computed satellite position table, including:

```text
t_GPS_s
toe_usado_s
tk_s
x_ECEF_m
y_ECEF_m
z_ECEF_m
lat_deg
lon_deg
h_m
```

## Main Result

For the assignment case:

```text
PRN17
t = 311400 s
toe = 309600 s
tk = 1800 s
```

The computed ECEF coordinates are approximately:

```text
x = -3780438.789 m
y = 14605294.608 m
z = 22300010.756 m
```

## Algorithm Summary

The implemented algorithm follows the standard GPS broadcast ephemeris procedure:

1. Read the RINEX navigation file.
2. Extract the ephemerides of the selected PRN.
3. Select the ephemeris with the closest `toe`.
4. Compute the time difference `tk`.
5. Compute the semi-major axis from `sqrtA`.
6. Compute the corrected mean motion.
7. Compute the mean anomaly.
8. Solve Kepler's equation for the eccentric anomaly.
9. Compute the true anomaly.
10. Apply harmonic corrections to the argument of latitude, orbital radius and inclination.
11. Compute the satellite position in the orbital plane.
12. Transform the coordinates to the ECEF frame.
13. Convert ECEF coordinates to approximate geodetic coordinates.

## Constants Used

```python
MU = 3.986005e14
OMEGA_E = 7.2921151467e-5
```

Where:

- `MU` is the Earth's universal gravitational parameter used by GPS.
- `OMEGA_E` is the Earth's rotation rate.

## Notes

This project is intended for educational use and for validating the implementation of the GPS broadcast ephemeris algorithm. It is not intended for operational navigation.

## Author

Daniel Miguel Tejedor

## License

Academic project. No specific license assigned.
