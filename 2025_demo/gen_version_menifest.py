import os
import re
import argparse
import configparser
from datetime import datetime

version = "0.4"
last_update = 20260412

VERSION_RE = re.compile(
    r'^\s*version\s*=\s*["\']?([A-Za-z0-9._-]+)["\']?',
    re.MULTILINE
)

LAST_UPDATE_RE = re.compile(
    r'^\s*last_update\s*=\s*["\']?([0-9./_-]+)["\']?',
    re.MULTILINE
)

EXCLUDE_DIRS = {
    ".venv", "venv", "__pycache__", ".git", ".idea", ".pytest_cache"
}

EXCLUDE_FILES = {
    "gen_version_menifest.py",
    "gen_version_manifest.py",
}

def extract_version_info(filepath):
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()
    except UnicodeDecodeError:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()

    version_match = VERSION_RE.search(content)
    update_match = LAST_UPDATE_RE.search(content)

    found_version = version_match.group(1) if version_match else None
    found_last_update = update_match.group(1) if update_match else None
    return found_version, found_last_update


def iter_python_files(directory=".", recursive=False):
    directory = os.path.abspath(directory)

    if recursive:
        for root, dirs, files in os.walk(directory):
            dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS]
            for filename in files:
                if not filename.endswith(".py"):
                    continue
                if filename in EXCLUDE_FILES:
                    continue
                yield os.path.join(root, filename)
    else:
        for filename in os.listdir(directory):
            filepath = os.path.join(directory, filename)
            if not os.path.isfile(filepath):
                continue
            if not filename.endswith(".py"):
                continue
            if filename in EXCLUDE_FILES:
                continue
            yield filepath


def preview_scan(directory=".", recursive=False):
    files = list(iter_python_files(directory, recursive))
    print(f"Scan root: {os.path.abspath(directory)}")
    print(f"Include child folders: {recursive}")
    print(f"Python files found: {len(files)}")
    print("-" * 50)

    for path in files[:30]:
        print(os.path.relpath(path, directory))

    if len(files) > 30:
        print(f"... and {len(files) - 30} more")

    return files


def generate_manifest(directory=".", recursive=False, ask_confirm=True):
    directory = os.path.abspath(directory)

    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")

    files = preview_scan(directory, recursive)

    if ask_confirm:
        answer = input("\nProceed to generate manifest from these files? (y/n): ").strip().lower()
        if answer != "y":
            print("Cancelled.")
            return

    config = configparser.ConfigParser()
    now = datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    config["version_check"] = {
        "checked_on": now,
        "scan_root": directory,
        "recursive": str(recursive),
    }

    skipped = []

    for filepath in files:
        found_version, found_last_update = extract_version_info(filepath)

        if found_version is None or found_last_update is None:
            skipped.append(os.path.relpath(filepath, directory))
            continue

        section_name = os.path.relpath(filepath, directory).replace("\\", "/").replace(".py", "")
        config[section_name] = {
            "version": str(found_version),
            "last_update": str(found_last_update),
        }

    outpath = os.path.join(directory, "project_version_manifest.ini")
    with open(outpath, "w", encoding="utf-8") as configfile:
        config.write(configfile)

    print(f"\n✅ Version manifest written to: {outpath}")

    if skipped:
        print("\nSkipped files missing version/last_update:")
        for s in skipped:
            print(f" - {s}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate project version manifest from Python files."
    )
    parser.add_argument(
        "--path",
        default=".",
        help="Exact folder to scan. Default: current folder"
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Include child folders"
    )
    parser.add_argument(
        "--no-confirm",
        action="store_true",
        help="Skip preview confirmation"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    generate_manifest(
        directory=args.path,
        recursive=args.recursive,
        ask_confirm=not args.no_confirm
    )


#python .\2025_demo\gen_version_menifest.py --path ".\2025_demo"