"""Prepare the raw-SMILES provenance record needed by source-cohort analysis.

The existing cached representations are compared with a fresh reconstruction.
No model is fitted and the cache arrays are not replaced. The resulting
CACHE_VALIDATION.json preserves the distinction between current validation
and the unknown generation environment of a historical cache.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from _common import CACHE, load_study_data
from _cache_validation import audit_cache_from_raw


def main():
    audit_cache_from_raw(load_study_data(), CACHE)
    print(CACHE / "CACHE_VALIDATION.json")


if __name__ == "__main__":
    main()
