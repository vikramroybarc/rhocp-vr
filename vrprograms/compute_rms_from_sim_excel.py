#!/usr/bin/env python3
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# ---------- Config: sim column positions (0-based) ----------
SIM_X_IDX = 1   # 2nd column in sim sheet
SIM_Y_IDX = 13  # 14th column in sim sheet
# ------------------------------------------------------------

# Regex for sheet names like "873K", "773K"
SHEET_TEMP_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*[Kk]\s*$")

def read_table_robust(path: Path) -> pd.DataFrame:
    if path is None:
        return None
    for enc in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            return pd.read_csv(path, sep=None, engine="python", encoding=enc)
        except Exception:
            continue
    return pd.read_csv(path)

def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")

def find_exp_file(exp_dir: Path, temp_str: str) -> Path | None:
    exact = exp_dir / f"EXP_{temp_str}K.txt"
    if exact.exists():
        return exact
    # case-insensitive fallback
    for p in exp_dir.glob(f"*{temp_str}K*.txt"):
        if p.name.lower() == f"exp_{temp_str.lower()}k.txt":
            return p
    for p in exp_dir.glob(f"*{temp_str}K*.txt"):
        if p.name.lower().startswith("exp_"):
            return p
    return None

def compute_error_df(exp_df: pd.DataFrame, sim_df: pd.DataFrame) -> pd.DataFrame | None:
    """
    Build error table aligned on experimental X:
      X_exp, Exp_Y, Sim_Y_interp, Error = (Exp_Y - Sim_Y_interp) / Exp_Y
    """
    if exp_df is None or sim_df is None:
        return None
    if exp_df.shape[1] < 2 or sim_df.shape[1] <= max(SIM_X_IDX, SIM_Y_IDX):
        return None

    x_exp = coerce_numeric(exp_df.iloc[:, 0])
    y_exp = coerce_numeric(exp_df.iloc[:, 1])
    x_sim = coerce_numeric(sim_df.iloc[:, SIM_X_IDX])
    y_sim = coerce_numeric(sim_df.iloc[:, SIM_Y_IDX])

    # Drop NaNs
    exp_mask = ~(x_exp.isna() | y_exp.isna())
    x_exp, y_exp = x_exp[exp_mask], y_exp[exp_mask]
    sim_mask = ~(x_sim.isna() | y_sim.isna())
    x_sim, y_sim = x_sim[sim_mask], y_sim[sim_mask]

    if len(x_exp) == 0 or len(x_sim) < 2:
        return None

    # Sort sim by X for interpolation
    order = np.argsort(x_sim.values)
    xs, ys = x_sim.values[order], y_sim.values[order]

    sim_interp = np.interp(x_exp.values, xs, ys)

    with np.errstate(divide="ignore", invalid="ignore"):
        err = (y_exp.values - sim_interp) / y_exp.values  # NaN where Exp_Y == 0

    return pd.DataFrame({
        "X_exp": x_exp.values,
        "Exp_Y": y_exp.values,
        "Sim_Y_interp": sim_interp,
        "Error": err
    })

def rms_error_from_series(err_series: pd.Series) -> float | None:
    vals = pd.to_numeric(err_series, errors="coerce").dropna().values
    if vals.size == 0:
        return None
    return float(np.sqrt(np.mean(vals**2)))

def main():
    if len(sys.argv) < 3:
        print(
            "Usage:\n"
            "  python compute_rms_from_sim_excel.py <SIM_EXCEL> <EXP_DIR>\n\n"
            "Example:\n"
            "  python compute_rms_from_sim_excel.py ./combined_3e-3.xlsx ./experiment\n"
        )
        sys.exit(1)

    sim_excel = Path(sys.argv[1]).expanduser().resolve()
    exp_dir   = Path(sys.argv[2]).expanduser().resolve()

    if not sim_excel.exists():
        print(f"[ERR] Simulation Excel not found: {sim_excel}")
        sys.exit(2)
    if not exp_dir.exists():
        print(f"[ERR] EXP_DIR not found: {exp_dir}")
        sys.exit(2)

    # Gather temperature-named sheets
    xls = pd.ExcelFile(sim_excel)
    temp_sheets = []
    for sheet in xls.sheet_names:
        m = SHEET_TEMP_RE.match(sheet)
        if m:
            temp_sheets.append((sheet, m.group(1)))  # (sheet_name, temp_str)

    summary_rows = []
    all_errors = []  # collect all per-point errors across temps

    # Append/replace Summary in the same file (requires openpyxl)
    with pd.ExcelWriter(sim_excel, mode="a", if_sheet_exists="replace", engine="openpyxl") as writer:
        for sheet_name, temp_str in temp_sheets:
            print(f"[INFO] Processing sheet: {sheet_name}")

            # Read simulation sheet
            sim_df = pd.read_excel(sim_excel, sheet_name=sheet_name, header=0)

            sim_ok = sim_df.shape[1] > SIM_Y_IDX
            if sim_ok:
                sim_df.iloc[:, SIM_X_IDX] = coerce_numeric(sim_df.iloc[:, SIM_X_IDX])
                sim_df.iloc[:, SIM_Y_IDX] = coerce_numeric(sim_df.iloc[:, SIM_Y_IDX])
            else:
                print(f"       [WARN] '{sheet_name}' has only {sim_df.shape[1]} columns; need >= {SIM_Y_IDX+1}")

            # Read experiment
            exp_path = find_exp_file(exp_dir, temp_str)
            exp_df = read_table_robust(exp_path) if exp_path else None
            if exp_df is not None and exp_df.shape[1] >= 2:
                exp_df.iloc[:, 0] = coerce_numeric(exp_df.iloc[:, 0])
                exp_df.iloc[:, 1] = coerce_numeric(exp_df.iloc[:, 1])
            else:
                exp_df = None
                print(f"       [INFO] No EXP file for {temp_str}K")

            # Compute error + RMS
            err_df = compute_error_df(exp_df, sim_df) if (sim_ok and exp_df is not None) else None
            if err_df is not None:
                # Accumulate all errors (numeric only)
                err_vals = pd.to_numeric(err_df["Error"], errors="coerce").dropna()
                all_errors.append(err_vals)

                temp_rms = rms_error_from_series(err_vals)
                n_pts = int(err_vals.shape[0])
            else:
                temp_rms = np.nan
                n_pts = 0

            # Summary row
            try:
                temp_val = float(temp_str)
            except Exception:
                temp_val = temp_str
            summary_rows.append({
                "Temperature": temp_val,
                "RMS_Error": temp_rms,
                "N_points": n_pts
            })

        # Compute overall RMS across all temperatures
        if len(all_errors) > 0:
            all_concat = pd.concat(all_errors, ignore_index=True)
            overall_rms = rms_error_from_series(all_concat)
            overall_n = int(all_concat.shape[0])
        else:
            overall_rms = np.nan
            overall_n = 0

        # Build Summary DataFrame (add final ALL row)
        summary_df = pd.DataFrame(summary_rows)

        # Sort by numeric Temperature if possible (leaves 'ALL' for end)
        try:
            summary_df["_t"] = pd.to_numeric(summary_df["Temperature"], errors="coerce")
            summary_df = summary_df.sort_values("_t").drop(columns="_t")
        except Exception:
            pass

        # Append overall row
        summary_df = pd.concat([
            summary_df,
            pd.DataFrame([{"Temperature": "ALL", "RMS_Error": overall_rms, "N_points": overall_n}])
        ], ignore_index=True)

        # Write/replace Summary sheet
        summary_df.to_excel(writer, index=False, sheet_name="Summary")

    print(f"[OK] Summary (with overall RMS) added/updated in: {sim_excel}")

if __name__ == "__main__":
    main()
