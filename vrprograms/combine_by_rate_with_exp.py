#!/usr/bin/env python3
import re
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# ----------------------- Engine selection -----------------------
def pick_excel_engine():
    try:
        import xlsxwriter  # noqa: F401
        return "xlsxwriter"
    except Exception:
        try:
            import openpyxl  # noqa: F401
            return "openpyxl"
        except Exception:
            raise RuntimeError(
                "Neither 'xlsxwriter' nor 'openpyxl' is installed.\n"
                "Install one of them:\n"
                "  pip install xlsxwriter   # (recommended for charts)\n"
                "or\n"
                "  pip install openpyxl"
            )

ENGINE = pick_excel_engine()
HAS_XLSXWRITER = (ENGINE == "xlsxwriter")

# ----------------------- Naming patterns ------------------------
# Simulation files: out_<temp>K_<rate>.csv  (e.g., out_873K_3e-3.csv)
SIM_RE = re.compile(r"^out_(\d+(?:\.\d+)?)K_([0-9.eE+\-]+)\.csv$")

# ----------------------- I/O helpers ----------------------------
def read_table_robust(path: Path) -> pd.DataFrame:
    """Read CSV/TXT with auto delimiter detection and robust encodings."""
    for enc in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            return pd.read_csv(path, sep=None, engine="python", encoding=enc)
        except Exception:
            continue
    return pd.read_csv(path)  # let pandas raise

def coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")

def find_exp_file(exp_dir: Path, temp_str: str) -> Path | None:
    """EXP files are named exactly: EXP_<temp>K.txt (e.g., EXP_873K.txt)."""
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

# ----------------------- Error metric ---------------------------
def compute_error(exp_df: pd.DataFrame, sim_df: pd.DataFrame,
                  sim_x_idx: int, sim_y_idx: int) -> pd.DataFrame | None:
    """
    Error(X) = (EXP_Y(X) - SIM_Y_interp(X)) / EXP_Y(X)
    exp_df: X=col0, Y=col1
    sim_df: X=sim_x_idx, Y=sim_y_idx (interp at exp X)
    Returns DataFrame [X_exp, Exp_Y, Sim_Y_interp, Error] or None.
    """
    if exp_df is None or sim_df is None:
        return None
    if exp_df.shape[1] < 2 or sim_df.shape[1] <= max(sim_x_idx, sim_y_idx):
        return None

    x_exp = coerce_numeric(exp_df.iloc[:, 0])
    y_exp = coerce_numeric(exp_df.iloc[:, 1])
    x_sim = coerce_numeric(sim_df.iloc[:, sim_x_idx])
    y_sim = coerce_numeric(sim_df.iloc[:, sim_y_idx])

    # drop NaNs
    exp_mask = ~(x_exp.isna() | y_exp.isna())
    x_exp, y_exp = x_exp[exp_mask], y_exp[exp_mask]
    sim_mask = ~(x_sim.isna() | y_sim.isna())
    x_sim, y_sim = x_sim[sim_mask], y_sim[sim_mask]

    if len(x_exp) == 0 or len(x_sim) < 2:
        return None

    # sort sim by x for interpolation
    order = np.argsort(x_sim.values)
    x_sim_sorted = x_sim.values[order]
    y_sim_sorted = y_sim.values[order]

    # interpolate sim_y at exp_x
    sim_y_interp = np.interp(x_exp.values, x_sim_sorted, y_sim_sorted)

    # avoid division by zero (gives NaN where Exp_Y == 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        error = (y_exp.values - sim_y_interp) / y_exp.values

    return pd.DataFrame({
        "X_exp": x_exp.values,
        "Exp_Y": y_exp.values,
        "Sim_Y_interp": sim_y_interp,
        "Error": error
    })

def rms_error(err_df: pd.DataFrame) -> float | None:
    if err_df is None or "Error" not in err_df.columns:
        return None
    vals = pd.to_numeric(err_df["Error"], errors="coerce").dropna().values
    if vals.size == 0:
        return None
    return float(np.sqrt(np.mean(vals**2)))

# ----------------------- Main -----------------------------------
def main():
    if len(sys.argv) < 4:
        print("Usage:\n  python combine_by_rate_with_exp.py <INPUT_DIR> <RATE> <EXP_DIR>")
        sys.exit(1)

    input_dir = Path(sys.argv[1]).expanduser().resolve()
    rate_arg  = sys.argv[2].strip()
    exp_dir   = Path(sys.argv[3]).expanduser().resolve()

    if not input_dir.exists():
        print(f"[ERR] INPUT_DIR not found: {input_dir}"); sys.exit(2)
    if not exp_dir.exists():
        print(f"[ERR] EXP_DIR not found: {exp_dir}"); sys.exit(2)

    # Collect simulation files for the given rate
    sim_files = []
    for p in list(input_dir.glob("*.csv")) + list(input_dir.glob("*.CSV")):
        m = SIM_RE.match(p.name)
        if not m:
            continue
        temp_str, rate_str = m.group(1), m.group(2)
        try:
            if float(rate_str) == float(rate_arg):
                sim_files.append((float(temp_str), temp_str, p))
        except Exception:
            if rate_str == rate_arg:
                sim_files.append((float(temp_str) if temp_str.replace('.', '', 1).isdigit() else 1e99, temp_str, p))

    if not sim_files:
        print(f"[WARN] No simulation files for rate '{rate_arg}' in {input_dir}")
        print("       Expected: out_<temp>K_<rate>.csv")
        sys.exit(0)

    # Sort by temperature
    sim_files.sort(key=lambda t: t[0])

    # Sim X/Y indices (0-based): X=2nd col => 1; Y=14th col => 13
    SIM_X_IDX = 1
    SIM_Y_IDX = 13

    out_xlsx = input_dir / f"combined_{rate_arg}.xlsx"
    print(f"[INFO] Writing workbook: {out_xlsx}")
    if not HAS_XLSXWRITER:
        print("[WARN] xlsxwriter not available -> charts will be skipped (data only).")

    summary_rows = []          # per-temp rows for Summary
    all_error_series = []      # collect all Error values across temps

    with pd.ExcelWriter(out_xlsx, engine=ENGINE) as writer:
        for temp_num, temp_str, sim_path in sim_files:
            sheet = f"{temp_str}K"[:31]
            print(f"[READ] Sim: {sim_path.name}  -> sheet '{sheet}'")

            # ----- Read simulation -----
            sim_df = read_table_robust(sim_path)
            sim_ok_for_metrics = sim_df.shape[1] > SIM_Y_IDX
            if sim_ok_for_metrics:
                sim_df.iloc[:, SIM_X_IDX] = coerce_numeric(sim_df.iloc[:, SIM_X_IDX])
                sim_df.iloc[:, SIM_Y_IDX] = coerce_numeric(sim_df.iloc[:, SIM_Y_IDX])
            else:
                print(f"       [WARN] Simulation file has only {sim_df.shape[1]} columns; "
                      f"need >= {SIM_Y_IDX+1} (for 14th col).")

            # ----- Read experiment for same temperature -----
            exp_path = find_exp_file(exp_dir, temp_str)
            exp_df = None
            if exp_path:
                print(f"[READ] Exp: {exp_path.name}")
                exp_df = read_table_robust(exp_path)
                if exp_df.shape[1] >= 1:
                    exp_df.iloc[:, 0] = coerce_numeric(exp_df.iloc[:, 0])  # X
                if exp_df.shape[1] >= 2:
                    exp_df.iloc[:, 1] = coerce_numeric(exp_df.iloc[:, 1])  # Y
            else:
                print(f"       [INFO] No experimental TXT found for {temp_str}K in {exp_dir} "
                      f"(expected 'EXP_{temp_str}K.txt').")

            # ----- Write tables to sheet -----
            sim_start_row, sim_start_col = 0, 0
            sim_df.to_excel(writer, index=False, sheet_name=sheet,
                            startrow=sim_start_row, startcol=sim_start_col)

            exp_start_row, exp_start_col = None, None
            if exp_df is not None:
                exp_start_row, exp_start_col = 0, 7  # column H
                exp_df.to_excel(writer, index=False, sheet_name=sheet,
                                startrow=exp_start_row, startcol=exp_start_col)

            # ----- Error metrics (and RMS) -----
            err_df = None
            temp_rms = None
            n_pts = 0
            if exp_df is not None and sim_ok_for_metrics:
                err_df = compute_error(exp_df, sim_df, SIM_X_IDX, SIM_Y_IDX)
                if err_df is not None:
                    # write error table starting at column O
                    err_start_row, err_start_col = 0, 14
                    err_df.to_excel(writer, index=False, sheet_name=sheet,
                                    startrow=err_start_row, startcol=err_start_col)
                    temp_rms = rms_error(err_df)
                    # collect numeric error values for global RMS
                    errs_num = pd.to_numeric(err_df["Error"], errors="coerce").dropna()
                    n_pts = errs_num.shape[0]
                    if n_pts > 0:
                        all_error_series.append(errs_num)

            # Collect for Summary
            try:
                temp_val = float(temp_str)
            except Exception:
                temp_val = np.nan
            summary_rows.append({
                "Temperature": temp_val if not np.isnan(temp_val) else temp_str,
                "RMS_Error": temp_rms if temp_rms is not None else np.nan,
                "N_points": int(n_pts)
            })

            # ----- Chart (xlsxwriter only) -----
            if HAS_XLSXWRITER:
                workbook  = writer.book
                worksheet = writer.sheets[sheet]
                chart = workbook.add_chart({"type": "scatter", "subtype": "straight_with_markers"})

                # Simulation series: X = col index 1, Y = col index 13 (if present)
                sim_rows = len(sim_df)
                if sim_rows >= 2 and sim_ok_for_metrics:
                    chart.add_series({
                        "name":       "Simulation (col2 vs col14)",
                        "categories": [sheet, sim_start_row + 1, sim_start_col + SIM_X_IDX,
                                              sim_start_row + sim_rows, sim_start_col + SIM_X_IDX],
                        "values":     [sheet, sim_start_row + 1, sim_start_col + SIM_Y_IDX,
                                              sim_start_row + sim_rows, sim_start_col + SIM_Y_IDX],
                        "marker":     {"type": "circle"},
                    })
                else:
                    print(f"       [WARN] Skipping simulation series for {sheet} (missing columns/rows).")

                # Experiment series: X = col 0, Y = col 1 (if present)
                if exp_df is not None and exp_df.shape[1] >= 2:
                    exp_rows = len(exp_df)
                    chart.add_series({
                        "name":       "Experiment",
                        "categories": [sheet, 1, exp_start_col + 0, exp_rows, exp_start_col + 0],
                        "values":     [sheet, 1, exp_start_col + 1, exp_rows, exp_start_col + 1],
                        "marker":     {"type": "square"},
                    })

                chart.set_title({"name": f"{temp_str}K  @  rate {rate_arg}"})
                chart.set_x_axis({"name": "X"})
                chart.set_y_axis({"name": "Y"})
                chart.set_legend({"position": "bottom"})

                # place chart below the longer of the two data blocks
                exp_len = len(exp_df) if exp_df is not None else 0
                chart_row = max(sim_rows, exp_len) + 3
                worksheet.insert_chart(chart_row, 0, chart, {"x_scale": 1.2, "y_scale": 1.2})

        # ----- Write Summary sheet (with ALL-data RMS) -----
        summary_df = pd.DataFrame(summary_rows)
        # sort numerically by Temperature if possible
        try:
            summary_df["_t"] = pd.to_numeric(summary_df["Temperature"], errors="coerce")
            summary_df = summary_df.sort_values("_t").drop(columns="_t")
        except Exception:
            pass

        # overall RMS across all temps
        if all_error_series:
            all_concat = pd.concat(all_error_series, ignore_index=True)
            overall_rms = float(np.sqrt(np.mean(all_concat.values**2))) if all_concat.size else np.nan
            overall_n = int(all_concat.shape[0])
        else:
            overall_rms = np.nan
            overall_n = 0

        # append final ALL row
        summary_df = pd.concat([
            summary_df,
            pd.DataFrame([{"Temperature": "ALL", "RMS_Error": overall_rms, "N_points": overall_n}])
        ], ignore_index=True)

        summary_df.to_excel(writer, index=False, sheet_name="Summary")

    print(f"[OK] Saved workbook with sheets per temperature + Summary (with ALL-data RMS): {out_xlsx}")
    if not HAS_XLSXWRITER:
        print("[NOTE] Install xlsxwriter to get charts:\n  pip install xlsxwriter")

if __name__ == "__main__":
    main()
