import sys, os, platform, traceback, importlib.util

def have(mod):
    return importlib.util.find_spec(mod) is not None

def show():
    print("=== GlycoMSP environment check ===")
    print(f"Python: {sys.version.split()[0]}  ({sys.executable})")
    print(f"Platform: {platform.platform()}")
    for name in ("tkinter","numpy","pandas","joblib","sklearn"):
        spec = importlib.util.find_spec(name)
        if not spec:
            print(f"MISS {name}")
            continue
        try:
            m = __import__(name)
            v = getattr(m,"__version__","present")
            print(f"OK   {name} {v}")
        except Exception as e:
            print(f"MISS {name} (import failed: {e})")
    if os.name == "nt":
        spec = importlib.util.find_spec("pymsfilereader")
        print(("OK   " if spec else "MISS ") + "pymsfilereader  (RAW conversion only)")
    if sys.version_info >= (3,13):
        print("\n[Notice] Python 3.13: some packages may lack wheels yet. "
              "Prefer 3.11/3.12 if installs fail.")

if __name__ == "__main__":
    try:
        show()
    finally:
        try:
            if os.name == "nt" and not sys.stdin.isatty():
                input("\nPress Enter to close…")
        except Exception:
            pass