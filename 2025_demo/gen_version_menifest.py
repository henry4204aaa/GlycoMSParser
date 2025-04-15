import os
import re
import configparser
from datetime import datetime
version = 0.2
last_update = 20250404

def extract_version_info(filepath):
    versions = None
    last_updates = None

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if re.match(r'\s*version\s*=', line):
                versions = re.findall(r'["\']?([\w\.\-]+)["\']?', line)[1]
            elif re.match(r'\s*last_update\s*=', line):
                last_updates = re.findall(r'["\']?([\w\-/]+)["\']?', line)[1]
            if versions and last_updates:
                break
    return versions, last_updates

def generate_manifest(directory=".", recursive=False): #directory = "."
    config = configparser.ConfigParser()
    
    # Add a [version_check] section for timestamp
    now = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
    config['version_check'] = {"checked_on": now}

    # Walk through files
    if recursive:
        for root, _, files in os.walk(directory):
            for filename in files:
                if filename.endswith(".py") and filename != "gen_version_manifest.py":
                    filepath = os.path.join(root, filename)
                    version, last_update = extract_version_info(filepath)
                    if version and last_update:
                        section_name = os.path.relpath(filepath, directory).replace("\\", "/").replace(".py", "")
                        config[section_name] = {
                            "version": version,
                            "last_update": last_update
                        }
    else:
        for filename in os.listdir(directory):
            if filename.endswith(".py") and filename != "gen_version_manifest.py":
                filepath = os.path.join(directory, filename)
                version, last_update = extract_version_info(filepath)
                if version and last_update:
                    module_name = os.path.splitext(filename)[0]
                    config[module_name] = {
                        "version": version,
                        "last_update": last_update
                    }

    with open(f"{directory}\project_version_manifest.ini", "w", encoding="utf-8") as configfile: 
        #changed "project_version_manifest.ini"  to f"{directory}\..." 
        config.write(configfile)
    print("✅ Version manifest written to project_version_manifest.ini in current folder")

if __name__ == "__main__":
    generate_manifest(recursive=False)  # Change to True when you're ready

#version 0.2 changelog
#change directory settings to make the script grab current folder correctly
#version 0.1 changelog
#basic version grabbing function deployed 