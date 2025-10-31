#!/usr/bin/env python3
import re
import sys
import math
from pathlib import Path

import pandas as pd

# ---------------- Config ----------------
INPUT_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
OUTPUT_XLSX = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("combined_creep.xlsx")
SHEET_SUFFIX = "64g"  # -> "<temp>K_<stress>MPa_64g"
# ---------------------------------------

def pick_excel_engine():
    """Prefer xlsxwriter; fall back to openpyxl."""
    try:
        import xlsxwriter  # noqa: F401
        return "xlsxwriter"
    except Exception:
        try:
            import openpyxl  # noqa: F401
            return "openpyxl"
        except Exception:
            raise RuntimeError(
                "Neither 'xlsxwriter' nor 'openpyxl' is installed. "
                "Install one: pip install xlsxwriter  (or)  pip install openpyxl"
            )

ENGINE = pick_excel_engine()

# Regex to extract temperature & stress; tolerates case and optional spaces
NAME_RE = re.compile(
    r"^out_creep_(\d+(?:\.\d+)?)\s*K_(\d+(?:\.\d+)?)\s*M[Pp]a\.csv$",
    re.IGNORECASE,
)

def pick_x(n_rows: int) -> int:
    if n_rows < 20:
        return 2
    elif n_rows < 30:
        return 5
    elif n_rows <= 40:
        return 8
    elif n_rows <= 60:
        return 10   
    elif n_rows <= 100:
        return 20     
    else:
        return 25

def read_csv_robust(path: Path) -> pd.DataFrame:
    for enc in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            df = pd.read_csv(path, sep=None, engine="python", encoding=enc)
            return df
        except Exception:
            continue
    return pd.read_csv(path)  # let pandas raise

def to_numeric_first_two(df: pd.DataFrame) -> pd.DataFrame:
    if df.shape[1] >= 1:
        df.iloc[:, 0] = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    if df.shape[1] >= 2:
        df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    return df

def last_segment_ratio(df: pd.DataFrame) -> float:
    if df.shape[1] < 2:
        return math.nan

    clean = df[[df.columns[0], df.columns[1]]].dropna()
    n = len(clean)
    if n < 2:
        return math.nan

    x = min(pick_x(n), n - 1)
    A = clean.iloc[:, 0]
    B = clean.iloc[:, 1]
    num = float(B.iloc[-1]) - float(B.iloc[-1 - x])
    den = float(A.iloc[-1]) - float(A.iloc[-1 - x])
    if den == 0:
        return math.nan
    return num / den

def main():
    candidates = list(INPUT_DIR.glob("*.csv")) + list(INPUT_DIR.glob("*.CSV"))
    files = [p for p in candidates if NAME_RE.match(p.name)]

    if not files:
        print(f"[WARN] No matching files found in {INPUT_DIR.resolve()}")
        with pd.ExcelWriter(OUTPUT_XLSX, engine=ENGINE) as writer:
            pd.DataFrame(columns=["Temperature", "Stress", "(B_last - B_last-x)/(A_last - A_last-x)"]).to_excel(
                writer, index=False, sheet_name="Summary"
            )
        print(f"[OK] Wrote empty workbook with Summary to: {OUTPUT_XLSX.resolve()}")
        return

    def key_fn(p: Path):
        m = NAME_RE.match(p.name)
        return (float(m.group(1)), float(m.group(2)))

    files.sort(key=key_fn)

    summary_rows = []
    with pd.ExcelWriter(OUTPUT_XLSX, engine=ENGINE) as writer:
        for path in files:
            m = NAME_RE.match(path.name)
            temp_str, stress_str = m.group(1), m.group(2)

            def clean_num_str(s: str) -> str:
                v = float(s)
                return s if ("." in s and not s.endswith(".0")) else str(int(v))

            temp_clean = clean_num_str(temp_str)
            stress_clean = clean_num_str(stress_str)

            sheet_name = f"{temp_clean}K_{stress_clean}MPa_{SHEET_SUFFIX}"
            sheet_name = sheet_name[:31]

            print(f"[READ] {path.name}  ->  sheet '{sheet_name}'")
            df = read_csv_robust(path)

            if not df.empty:
                df = to_numeric_first_two(df)
            else:
                print(f"       [WARN] {path.name} parsed as empty. Writing empty sheet.")

            df.to_excel(writer, index=False, sheet_name=sheet_name)

            ratio = last_segment_ratio(df)
            summary_rows.append({
                "Temperature": temp_clean,
                "Stress": stress_clean,
                "(B_last - B_last-x)/(A_last - A_last-x)": ratio
            })

        summary_df = pd.DataFrame(summary_rows)
        try:
            summary_df["_t"] = summary_df["Temperature"].astype(float)
            summary_df["_s"] = summary_df["Stress"].astype(float)
            summary_df = summary_df.sort_values(["_t", "_s"]).drop(columns=["_t", "_s"])
        except Exception:
            pass
        summary_df.to_excel(writer, index=False, sheet_name="Summary")

    print(f"[OK] Wrote Excel workbook to: {OUTPUT_XLSX.resolve()}")
    print(f"[INFO] Engine used: {ENGINE}")

if __name__ == "__main__":
    main()
