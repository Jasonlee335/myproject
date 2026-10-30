"""Verify public code and aggregates without reading source datasets."""

import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "FILE_MANIFEST.json").read_text(encoding="utf-8"))
    failures = []
    for entry in manifest:
        path = root / entry["file"]
        if not path.is_file():
            failures.append(entry["file"] + ": missing")
        elif hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            failures.append(entry["file"] + ": checksum mismatch")
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"Verified {len(manifest)} deposited files.")


if __name__ == "__main__":
    main()
