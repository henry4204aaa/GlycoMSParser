
def compare_against_manual(
    df_filtered: pd.DataFrame,
    manual_excel_path: str,
    manual_sheet: str = "MSlist",
    keep_top_n: int = DEFAULTS["keep_top_n"]
):
    xls = pd.ExcelFile(manual_excel_path)
    if manual_sheet not in xls.sheet_names:
        raise ValueError(f"Sheet '{manual_sheet}' not found in {manual_excel_path}.")
    manual = xls.parse(manual_sheet)

    if "MS2scan_no" not in manual.columns or "Structure" not in manual.columns:
        raise ValueError("Manual sheet must contain 'MS2scan_no' and 'Structure' columns.")

    manual2 = manual[["MS2scan_no", "Structure"]].dropna(subset=["MS2scan_no"]).copy()
    manual2["MS2scan_no"] = pd.to_numeric(manual2["MS2scan_no"], errors="coerce").astype("Int64")
    manual2 = manual2.dropna(subset=["MS2scan_no"]).drop_duplicates(subset=["MS2scan_no"], keep="first")

    if df_filtered is None or len(df_filtered) == 0:
        joined = manual2.copy()
        joined["composition"] = pd.NA
        joined["ion score"] = pd.NA
        joined["ion hit count"] = pd.NA
        joined["ppm_error"] = pd.NA
        joined["match"] = False
        tp = 0
        fp = 0
        fn = int(len(joined))
        precision = 0.0
        recall = 0.0
        f1 = 0.0
        metrics = dict(TP=tp, FP=fp, FN=fn, precision=precision, recall=recall, f1=f1)
        return joined, metrics

    pseudo2 = df_filtered[["MS2scan_no", "composition", "ion score", "ion hit count", "ppm_error"]].copy()
    pseudo2["MS2scan_no"] = pd.to_numeric(pseudo2["MS2scan_no"], errors="coerce").astype("Int64")

    joined = manual2.merge(pseudo2, on="MS2scan_no", how="left", suffixes=("_manual", "_pseudo"))
    joined["match"] = (joined["Structure"].astype(str).str.strip() == joined["composition"].astype(str).str.strip())

    tp = int((joined["match"] == True).sum())
    fp = int((joined["match"] == False).sum())
    fn = int(joined["composition"].isna().sum())

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    metrics = dict(TP=tp, FP=fp, FN=fn, precision=precision, recall=recall, f1=f1)
    return joined, metrics
