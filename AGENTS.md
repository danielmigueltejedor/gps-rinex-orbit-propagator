# Agent instructions — gps-rinex-orbit-propagator

## Repository-specific deterministic quality contract (2026-10-26)

**Project:** GPS RINEX navigation ephemeris orbit propagation Python script.

**Engineering requirements:** Preserve time, units, Kepler solver accuracy and WGS84 coordinate transforms; compare sampled known ephemerides without silently adjusting tolerances.

**Verification commands:** python -m py_compile calculo_prn17.py; verify Modified_RINEX_I2.rnx exists; run scripted numerical fixture checks when added.

Before completion, run relevant checks and report actual outputs, including any failures. Add regression tests for changed behavior, edge cases and error paths. Never weaken/remove/skip tests, suppress failures or lower CI thresholds to force a pass. Changes to `QUALITY_GATES.md`, workflows or these rules need human review. Configure the `quality` status as a required check in branch protection. If tooling is unavailable, report it instead of claiming checks passed. Preserve existing project instructions.
