# Descriptors Formulation Script.py
# Code written by Ryan Cheung in March 2026 on behalf of Lancaster University working as part of a 4th year industrial group project. 
# The project is titled: "Input feature generation to build machine learning-based quantitative structure activity relationships (QSARs) for formulated products."
# Copyright (c) 2026 Ryan Cheung. Licensed under the MIT License, see LICENSE for more details. See CITATION.cff for information on how to cite this project. GitHub Repository:

# This script reads a table of compounds from CSV or Excel, resolves each entry through PubChem/OPSIN, computes configurable molecular descriptors and selected external properties, and writes the results to an Excel spreadsheet.
# Usage example:
# python "Descriptors Formulation Script.py" --input Input.csv --output features.xlsx

# -----------------------------
# Imports needed for synonyms and config (before any optional imports)
from __future__ import annotations
import re
from typing import Dict, Optional, List, Tuple, Any

SYNONYMS: Dict[str, str] = {
    # Add more as you discover failures
}

# Default configuration for which descriptors to compute. You can modify this or provide your own JSON config file.
# See the properties_config.json file for an example.
DEFAULT_STRUCTURE_CONFIG = {
    "rdkit": {
        "enabled": True,
        "families": {
            "estate": False,
            "qed_physchem": False,
            "charge": False,
            "fingerprint_density": False,
            "topological_indices": False,
            "connectivity_kappa": False,
            "surface_polarity": False,
            "vsa_estate": False,
            "counts_basic": False,
            "ring_counts": False,
            "hydrogen_bonding": False,
            "functional_groups": False,
        },
        "include_descriptors": [
            "MolWt",
            "ExactMolWt",
            "MolLogP",
            "MolMR",
            "TPSA",
            "FractionCSP3",
            "HeavyAtomCount",
            "NumHAcceptors",
            "NumHDonors",
            "NumHeteroatoms",
            "NumRotatableBonds",
            "RingCount",
            "NumAromaticRings",
            "NumAliphaticRings",
            "MaxPartialCharge",
            "MinPartialCharge",
            "BalabanJ",
            "BertzCT",
            "HallKierAlpha",
            "Kappa1",
            "Kappa2",
            "Kappa3"
        ],
        "exclude_descriptors": [],
        "fingerprints": {
            "morgan": {
                "enabled": False,
                "radius": 2,
                "nbits": 2048,
                "format": "bitstring",
            },
            "maccs": {
                "enabled": False,
                "format": "bitstring",
            }
        }
    },
    "mordred": {
        "enabled": False,
        "use_all": False,
        "families": {
            "constitutional": False,
            "adjacency_matrix": False,
            "autocorrelation": False,
            "topological": False,
            "charge": False,
            "surface": False,
            "ring": False,
            "fragment": False,
            "estate": False,
            "geometry_3d": False
        },
        "include_descriptors": [],
        "exclude_descriptors": []
    },
    "external_properties": {
        "enabled": True,
        "pubchem_rest": {
            "enabled": True,
            "properties": []
        },
        "pubchem_pug_view": {
            "enabled": True,
            "extract_boiling_point": True,
            "extract_vapour_pressure": True,
            "extract_density": True,
            "extract_water_solubility": True
        },
        "thermo": {
            "enabled": True,
            "temperature_K": 298.15,
            "pressure_Pa": 101325
        }
    }
}

# This part of the code defines the descriptor families for RDKit and Mordred, which are used to group related descriptors together. 
# This allows users to enable or disable entire families of descriptors in the configuration, and the code will automatically include all descriptors in that family. The `descriptor_in_family` function checks if a given descriptor name belongs to a specified family, which is useful for dynamically determining which descriptors to calculate based on the configuration.
# This is a Work In Progress! It won't be 100% correct, and some properties (especially mordred) will spill into different categories, but for now it will do, I think. Also bear in mind my chemistry knowledge isn't the best so trying to categorise about 1800 descriptors is very hard, so feel free to change it up!
RDKIT_FAMILIES = {
    "estate": {
        "MaxEStateIndex", "MinEStateIndex", "MaxAbsEStateIndex", "MinAbsEStateIndex"
    },
    "qed_physchem": {
        "qed", "MolWt", "HeavyAtomMolWt", "ExactMolWt",
        "MolLogP", "MolMR", "FractionCSP3"
    },
    "charge": {
        "MaxPartialCharge", "MinPartialCharge",
        "MaxAbsPartialCharge", "MinAbsPartialCharge"
    },
    "fingerprint_density": {
        "FpDensityMorgan1", "FpDensityMorgan2", "FpDensityMorgan3"
    },
    "topological_indices": {
        "BalabanJ", "BertzCT", "Ipc"
    },
    "connectivity_kappa": {
        "Chi0", "Chi0n", "Chi0v",
        "Chi1", "Chi1n", "Chi1v",
        "Chi2n", "Chi2v",
        "Chi3n", "Chi3v",
        "Chi4n", "Chi4v",
        "HallKierAlpha",
        "Kappa1", "Kappa2", "Kappa3"
    },
    "surface_polarity": {
        "LabuteASA", "TPSA"
    },
    "counts_basic": {
        "HeavyAtomCount", "NHOHCount", "NOCount",
        "NumValenceElectrons", "NumRadicalElectrons",
        "NumHeteroatoms", "NumRotatableBonds"
    },
    "ring_counts": {
        "NumAliphaticCarbocycles", "NumAliphaticHeterocycles", "NumAliphaticRings",
        "NumAromaticCarbocycles", "NumAromaticHeterocycles", "NumAromaticRings",
        "NumSaturatedCarbocycles", "NumSaturatedHeterocycles", "NumSaturatedRings",
        "RingCount"
    },
    "hydrogen_bonding": {
        "NumHAcceptors", "NumHDonors"
    }
}

MORDRED_FAMILIES = {
    "constitutional": {
        "ABCIndex",
        "AcidBase",
        "Aromatic",
        "AtomCount",
        "BondCount",
        "CarbonTypes",
        "Constitutional",
        "Framework",
        "FragmentComplexity",
        "Weight",
        "McGowanVolume",
        "VdwVolumeABC",
        "Lipinski",
        "HydrogenBond",
        "RotatableBond",
        "RingCount",
        "TopoPSA",
        "SLogP",
        "LogS",
    },
    "adjacency_matrix": {
        "AdjacencyMatrix",
        "DistanceMatrix",
        "DetourMatrix",
        "BaryszMatrix",
    },
    "autocorrelation": {
        "Autocorrelation",
        "MoeType",
        "MoRSE",
    },
    "topological": {
        "BalabanJ",
        "BertzCT",
        "Chi",
        "InformationContent",
        "KappaShapeIndex",
        "PathCount",
        "TopologicalCharge",
        "TopologicalIndex",
        "VertexAdjacencyInformation",
        "WalkCount",
        "WienerIndex",
        "ZagrebIndex",
        "MolecularDistanceEdge",
        "MolecularId",
        "EccentricConnectivityIndex",
        "ExtendedTopochemicalAtom",
    },
    "charge": {
        "BCUT",
        "Polarizability",
    },
    "surface": {
        "CPSA",
        "PBF",
        "SurfaceArea",
    },
    "ring": {
        "RingCount",
    },
    "fragment": {
        "FragmentComplexity",
        "Framework",
    },
    "estate": {
        "EState",
    },
    "geometry_3d": {
        "GeometricalIndex",
        "GravitationalIndex",
        "MomentOfInertia",
    }
}

# Necessary imports for the main code (after config)
import argparse
import json
import time
import os
import uuid
import random
from dataclasses import dataclass, asdict
from pathlib import Path

import pandas as pd
import requests

from rdkit import Chem
from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.Chem import rdFingerprintGenerator

import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

# Optional: Progress bar
try:
    from tqdm import tqdm
    HAS_TQDM = True
except Exception:
    HAS_TQDM = False

# Optional: mordred 
try:
    from mordred import Calculator, descriptors as mordred_descriptors
    HAS_MORDRED = True
except Exception:
    HAS_MORDRED = False

# Optional: thermo 
try:
    from thermo.chemical import Chemical
    HAS_THERMO = True
except Exception:
    HAS_THERMO = False


# Dataclass to hold resolved molecule information
@dataclass
class ResolvedMolecule:
    name: str
    smiles: Optional[str] = None
    cid: Optional[int] = None
    iupac_name: Optional[str] = None
    molecular_formula: Optional[str] = None
    pubchem_xlogp: Optional[float] = None
    pubchem_molwt: Optional[float] = None
    resolution_source: Optional[str] = None
    resolution_note: Optional[str] = None

# Cleaning and normalisation functions for input names and lookup candidates
def clean_text(s: Any) -> str:
    if s is None:
        return ""
    try:
        if pd.isna(s):
            return ""
    except Exception:
        pass

    s = str(s).strip()
    s = re.sub(r"\s+", " ", s) # Collapse multiple whitespaces to single space
    s = s.replace("\u2013", "-").replace("\u2014", "-") # normalise dashes
    return s

# Clean lookup names by removing common non-structural words and normalising whitespace.
def clean_lookup_name(name: str) -> str:
    n = clean_text(name)
    # Remove major/minor that sometimes appear in trade names (like linalool major)
    n = re.sub(r"\bmajor\b", "", n, flags=re.IGNORECASE).strip() 
    n = re.sub(r"\bminor\b", "", n, flags=re.IGNORECASE).strip()
    n = re.sub(r"\s+", " ", n).strip()
    return n

# Pre-clean synonyms keys for faster lookup
SYNONYMS = {clean_lookup_name(k).lower(): v for k, v in SYNONYMS.items()}

# Stereochemical prefixes to remove for lookup variants
_STEREO_PREFIXES = [
    r"^\(\s*[RrSs]\s*\)\s*[-–—]?\s*",     # (R)-, (S)-
    r"^\(\s*[EeZz]\s*\)\s*[-–—]?\s*",     # (E)-, (Z)-
    r"^\(\s*\+\s*\)\s*[-–—]?\s*",         # (+)-
    r"^\(\s*-\s*\)\s*[-–—]?\s*",          # (-)-
    r"^\+\s*", r"^-\s*",                  # + / - at start
    r"^(?:d|l|dl)\s*[-–—]\s*",            # d-, l-, dl-
    r"^(?:cis|trans)\s*[-–—]\s*",         # cis-, trans-
]

def strip_stereo_prefixes(s: str) -> str:
    out = s
    for pat in _STEREO_PREFIXES:
        out2 = re.sub(pat, "", out, flags=re.IGNORECASE)
        out = out2
    return out.strip()

# Normalise common Greek letters to words for better matching (e.g., α- -> alpha-)
def normalise_greek_words(s: str) -> str:
    # optional: helps when people type "alpha-" vs "α-"
    rep = {
        "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta" # etc, but these are the most common ones
    }
    out = s
    for k, v in rep.items():
        out = out.replace(k, v)
    return out

# Generate lookup variants by cleaning the name, removing stereo prefixes, normalising Greek letters, and applying synonyms. This helps to increase the chances of successful resolution through PubChem/OPSIN.
# However there needs to be a balance as we want the variants to be similar enough to the original name to still match the intended compound, and not explode into too many unrelated candidates.
def basic_variants(raw: str) -> List[str]:
    s = clean_lookup_name(raw)
    if not s:
        return []

    s = s.replace("®", "").replace("™", "").strip()
    s = s.rstrip(" -–—")  # fix trailing "-" that happens sometimes

    variants = [] # Keep track of variants to avoid duplicates
    seen = set()

    def add(x: str):
        x = clean_lookup_name(x)
        if x and x.lower() not in seen:
            seen.add(x.lower())
            variants.append(x)

    add(s)
    add(normalise_greek_words(s))
    # Remove stereo prefixes (but keep the original first!)
    add(strip_stereo_prefixes(s))
    # Also try removing stereo prefixes after greek normalisation
    add(strip_stereo_prefixes(normalise_greek_words(s)))

    return variants

# Here I decided to split candidates according to strong/weak seperators
# For example, "Tripal/Ligustral 1" should split into "Tripal" and "Ligustral 1", but "2,6-dimethylpyrazine" should not split on the comma.
_SPLIT_PAT_STRONG = re.compile(r"\s*(?:/|;|&)\s*", flags=re.IGNORECASE)
_SPLIT_PAT_WEAK   = re.compile(r"\s*(?:/|;|\bor\b|\band\b|&)\s*", flags=re.IGNORECASE)

# Function that does this
def split_candidates(v: str) -> List[str]:
    pat = _SPLIT_PAT_STRONG if ("(" in v or ")" in v) else _SPLIT_PAT_WEAK
    parts = [p.strip() for p in pat.split(v) if p.strip()]

    out: List[str] = []
    for p in parts:
        if " + " in p:
            out.extend([q.strip() for q in p.split(" + ") if q.strip()])
        else:
            out.append(p)

    return out

# Generate lookup candidates function
def generate_lookup_candidates(raw: str) -> List[str]:
    variants = basic_variants(raw)
    if not variants:
        return [] # No variants

    candidates: List[str] = []
    seen = set()

    def add(x: str):
        x = clean_lookup_name(x)
        if x and x.lower() not in seen:
            seen.add(x.lower()) # Track
            candidates.append(x) 

    for v in variants:
        # If it looks IUPAC-ish (commas + digits), DO NOT split it. Example: "1,2,3,4-tetrahydronaphthalene"
        if "," in v and re.search(r"\d", v):
            add(v)
            continue

        # Otherwise split only on true “alternatives” separators (/, and, or, &, +, ;)
        parts = split_candidates(v)
        if not parts:
            continue

        for p in parts:
            add(p)
            # Try dropping a trailing grade number: "Ligustral 1" -> "Ligustral"
            p2 = re.sub(r"\s+\d{1,3}$", "", p).strip()
            if p2 != p:
                add(p2)

    return candidates

# Make sure the input dataframe has standardised columns for name, smiles, and cas, and clean the name cells. This allows the rest of the code to rely on these standard columns existing and being clean.
def normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [clean_text(c).lower() for c in df.columns]

    rename_map = {}
    # Note that this is as forgiving as possible. For example, if you have both "name" and "compound" columns, they will both be renamed to "name" and the code will use the first one it sees. It’s better to have some duplicates than to miss a column due to strict matching.
    for c in df.columns:
        if c in {"name", "prm", "compound", "chemical", "molecule"}:
            rename_map[c] = "name"
        if c in {"smiles", "canonical_smiles", "isomeric_smiles"}:
            rename_map[c] = "smiles"
        if c in {"cas", "cas#", "cas_number", "cas number"}:
            rename_map[c] = "cas"

    df = df.rename(columns=rename_map)

    if "name" not in df.columns:
        raise ValueError("Input must contain a name column (e.g., Name/PRM/compound/chemical).")
    if "smiles" not in df.columns:
        df["smiles"] = ""
    if "cas" not in df.columns:
        df["cas"] = ""
    df["name"] = df["name"].apply(clean_text)
    return df

# Normalise cache keys- lowercase
def norm_cache_key(s: str) -> str:
    return clean_lookup_name(s).lower()

# Current UTC time in ISO format for timestamping in cache
def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()

# Treat as failure record if SMILES is missing/empty
def is_failure_record(d: dict) -> bool:
    return not d.get("smiles")

# HTTP status codes
# 429: Too Many Requests (Slow down!)
# 500: Internal Server Error (Someone's getting fired)
# 502: Bad Gateway (Check https://www.cloudflarestatus.com/)
# 503: Service Unavailable (Check your internet connection)
# 504: Gateway Timeout (Don't have all day, do you?)
_RETRY_STATUSES = {429, 500, 502, 503, 504} # Any more?

def get_json_retry(session: requests.Session, url: str, timeout_s: float, tries: int = 5): # 5 tries is enough
    for attempt in range(tries):
        r = session.get(url, timeout=timeout_s)
        if r.status_code == 200:
            return r.json(), 200 # 200 means it worked!
        if r.status_code in _RETRY_STATUSES:
            retry_after = r.headers.get("Retry-After") # Sometimes servers tell you how long to wait, respect that
            if retry_after:
                try:
                    sleep_s = float(retry_after)
                except ValueError:
                    sleep_s = 0.5 * (2 ** attempt) # Exponential backoff if header is malformed
            else:
                sleep_s = 0.5 * (2 ** attempt) # Exponential backoff if no Retry-After header
            sleep_s *= (0.75 + 0.5 * random.random())  # Add jitter to avoid a thundering herd
            time.sleep(min(sleep_s, 8.0))
            continue
        return None, r.status_code # For non-retryable errors, return None and the status code
    return None, 503 # After retries, give up and return 503

# MoleculeResolver: Handles names to SMILES using PubChem and OPSIN. This is the core of the name resolution process.
class MoleculeResolver:
    def __init__(
        self,
        cache_path: Path = Path("smiles_cache.json"), # Cache file, can change the location or name if you want
        sleep_s: float = 0.12,
        timeout_s: float = 20.0,
        user_agent: str = "fragrance-features/1.0 (name2smiles; rdkit)", # Server identification for HTTP requests- be nice to the servers!
        max_candidates: int = 6,
        failure_ttl_days: int = 7, # six seven
    ):
        self.cache_path = cache_path
        self.sleep_s = float(sleep_s)
        self.timeout_s = float(timeout_s)
        self.headers = {"User-Agent": user_agent}

        self.max_candidates = int(max_candidates)
        self.failure_ttl_days = int(failure_ttl_days)

        self.cache: Dict[str, Dict] = self._load_cache()
        self._tls = threading.local() # Threading to make it quicker, while being mindful of rate limits etc
        self._cache_lock = threading.Lock()  # Avoid file corruption when using threads
        self._pubchem_sem = threading.Semaphore(2) # 2 threads is apparently nice to PubChem's servers
        self._opsin_sem = threading.Semaphore(2) # Again, 2 threads for OPSIN 

        self._cache_dirty = False # Track if cache has unsaved changes
        self._cache_writes = 0
        self.cache_flush_every = 50  # Can be adjusted based on how often you want to write to disk vs risk losing data

    # Cache handling functions with locking for thread safety
    def _cache_get(self, key: str) -> Optional[Dict]:
        with self._cache_lock:
            value = self.cache.get(key)
            return dict(value) if value is not None else None # Return a copy to avoid accidental modifications outside the lock

    # Set cache with locking and mark as dirty for later flushing to disk
    def _cache_set(self, key: str, value: Dict) -> None:
        with self._cache_lock:
            self.cache[key] = dict(value)

    # Get a JSON snapshot of the cache for saving to disk, with locking to ensure consistency
    def _cache_snapshot_json(self) -> str:
        with self._cache_lock:
            return json.dumps(self.cache, indent=2, ensure_ascii=False)

    # Mark cache as dirty and flush to disk if we've reached the flush threshold
    def _mark_cache_dirty(self):
        self._cache_dirty = True
        self._cache_writes += 1
        if self._cache_writes % self.cache_flush_every == 0: # Flush to disk every N writes to avoid frequent disk I/O
            self._cache_dirty = False

    # Check if a failed cache entry is old enough to retry (i.e., if the TTL has expired)
    def _get_session(self) -> requests.Session:
        s = getattr(self._tls, "session", None)
        if s is None:
            s = requests.Session()
            s.headers.update(self.headers)
            self._tls.session = s
        return s

    # Determine if a cache entry that is a failure record is old enough to retry based on the failure timestamp and the configured TTL
    def _load_cache(self) -> Dict[str, Dict]:
        if self.cache_path.exists():
            try:
                return json.loads(self.cache_path.read_text(encoding="utf-8")) # UTF-8 is important for some chemical names with special characters (e.g., δ-cadinol)
            except Exception: 
                return {}
        return {}

    # Save the cache by writing to a temporary file and then replacing the old cache file
    def _save_cache(self) -> None:
        data = self._cache_snapshot_json()

        tmp = self.cache_path.with_name(self.cache_path.name + f".{uuid.uuid4().hex}.tmp")
        tmp.write_text(data, encoding="utf-8")

        for attempt in range(10): # Retry in case of syncing issues (OneDrive cough)
            try:
                os.replace(tmp, self.cache_path)
                return
            except PermissionError:
                time.sleep(0.05 * (attempt + 1)) # Wait a bit longer each time
        # If we fail after retries, try to clean up the temp file and raise an error with instructions for the user
        try:
            if tmp.exists():
                tmp.unlink()
        except Exception:
            pass

        raise PermissionError(
            f"Failed to replace cache file after retries: {self.cache_path}. "
            f"Close any programs using it and/or move cache outside OneDrive." # Sometimes I don't like OneDrive
        )

    # Main resolution function for PubChem and OPSIN
    def resolve(self, name: str) -> ResolvedMolecule: 
        original = clean_lookup_name(name)
        if original.endswith("-"):
            repaired = original.rstrip("-").strip()
            if repaired and repaired != original:
                original = repaired
        if not original:
            return ResolvedMolecule(name=name, resolution_note="Empty name")
        # Generate lookup candidates from the original name
        candidates = generate_lookup_candidates(original) or [original]

        expanded: List[str] = [] 
        seen: set[str] = set() # Track seen candidates to avoid duplicates

        for c in candidates: # Candidate lookup
            c_key = norm_cache_key(c)
            if c_key not in seen:
                expanded.append(c)
                seen.add(c_key)

            syn = SYNONYMS.get(c_key) # Synonym lookup
            if syn:
                syn_key = norm_cache_key(syn)
                if syn_key not in seen:
                    expanded.append(syn)
                    seen.add(syn_key)

        # Limit the number of candidates to avoid excessive lookups

        best_fail: ResolvedMolecule | None = None # Keep track of best failure record
        tried: List[str] = [] # Store candidate attempts here

        for lookup in expanded: # Lookup loop with caching and retries
            tried.append(lookup)
            cache_key = norm_cache_key(lookup)

            cached_dict = self._cache_get(cache_key)
            if cached_dict:
                try:
                    cached = ResolvedMolecule(**{
                        k: v for k, v in cached_dict.items() 
                        if k in ResolvedMolecule.__annotations__
                    })
                    if cached.smiles:
                        cached.name = name
                        return cached

                    # Only retry failed cache entries if TTL expired
                    if not self._cache_allows_retry(cached_dict):
                        continue
                except Exception:
                    pass

            mol = self._resolve_pubchem(lookup) # Try PubChem first
            if mol.smiles:
                mol.name = name
                mol.resolution_note = f"Matched via candidate: {lookup} | " + (mol.resolution_note or "") # Keep track of which candidate worked
                self._cache_set(cache_key, asdict(mol))
                self._mark_cache_dirty() # Mark cache as dirty for later flushing to disk
                time.sleep(self.sleep_s)
                return mol

            mol_opsin = self._resolve_opsin(lookup) # Now try OPSIN as a fallback as it can handle IUPAC names
            if mol_opsin.smiles:
                mol_opsin.name = name
                mol_opsin.resolution_note = f"Matched via candidate: {lookup} | " + (mol_opsin.resolution_note or "")
                self._cache_set(cache_key, asdict(mol_opsin))
                self._mark_cache_dirty()
                time.sleep(self.sleep_s)
                return mol_opsin # Same as above

            best_fail = mol

            # Cache failures with timestamp so you don’t repeatedly hammer PubChem
            fail_record = asdict(best_fail)
            fail_record["failed_at"] = utc_now_iso()

            # Don't cache transient errors as hard failures
            note = (best_fail.resolution_note or "").lower()
            if "http 503" in note or "http 429" in note:
                # keep in-memory only (or cache with very short TTL)
                time.sleep(self.sleep_s)
                continue

            self._cache_set(cache_key, fail_record)
            self._mark_cache_dirty()

            time.sleep(self.sleep_s)

        if best_fail is None:
            best_fail = ResolvedMolecule(name=name, resolution_note="No candidates") # Technically shouldn't happen

        best_fail.name = name
        best_fail.resolution_note = f"Failed candidates tried: {tried} | " + (best_fail.resolution_note or "") # Keep track of which candidates were tried
        return best_fail

    # Resolve a name through PubChem's REST API, with retries and fallbacks. This is the main function that tries to get the SMILES and other properties from PubChem, and if it fails, it will try OPSIN as a fallback.
    def _resolve_pubchem(self, name: str) -> ResolvedMolecule:
        props = [
            "CanonicalSMILES",
            "IsomericSMILES",
            "IUPACName",
            "MolecularFormula",
            "MolecularWeight",
            "XLogP", 
        ] # You can add more but be careful of running into other issues

        def fetch(url: str):
            # Semaphore protects PubChem concurrency
            with self._pubchem_sem:
                session = self._get_session()
                return get_json_retry(session, url, self.timeout_s)

        url1 = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" # PubChem name lookup endpoint
            + requests.utils.quote(name)
            + "/property/"
            + ",".join(props)
            + "/JSON" # Request the properties in JSON format for easier parsing
        )

        try:
            data, status = fetch(url1) # First try: name -> properties (including SMILES if we're lucky)

            if status != 200 or data is None: # If we get a non-200 status or no data, return a failure record with the status code in the note
                return ResolvedMolecule(
                    name=name,
                    resolution_source="pubchem",
                    resolution_note=f"HTTP {status}",
                ) 

            # Name -> properties (including SMILES)
            props_list = data.get("PropertyTable", {}).get("Properties", []) # Take the properties list from the response
            if not props_list:
                return ResolvedMolecule(
                    name=name,
                    resolution_source="pubchem",
                    resolution_note="No properties returned", # If we get a 200 but no properties, that’s also a failure
                )

            p0 = props_list[0] # Take the first match (there can be multiple, but we’ll just take the first one for simplicity)
            cid = p0.get("CID")
            smiles = p0.get("CanonicalSMILES") or p0.get("IsomericSMILES") # Either canonical or isomeric SMILES, whichever is available

            resolution_source = "pubchem"
            resolution_note = "OK (name->props)" # :)

            # If SMILES is missing, try some fallbacks before giving up:
            if (not smiles) and cid:
                url2 = (
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" # PubChem CID lookup endpoint
                    + str(cid)
                    + "/property/CanonicalSMILES,IsomericSMILES/JSON" # In JSON
                )
                data2, status2 = fetch(url2) # CID -> properties (including SMILES)
                if status2 == 200 and data2:
                    props_list2 = data2.get("PropertyTable", {}).get("Properties", []) # Take the properties list from the second response
                    if props_list2:
                        p2 = props_list2[0]
                        smiles = p2.get("CanonicalSMILES") or p2.get("IsomericSMILES") # Again eitehr canonical or isomeric SMILES

            # If still missing, CID -> SDF -> RDKit -> SMILES
            if (not smiles) and cid:
                sdf_url = (
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" # PubChem CID lookup for SDF format, which sometimes has better sanitisation than the property endpoint
                    + str(cid)
                    + "/SDF" # SDF format
                )
                try:
                    with self._pubchem_sem:
                        session = self._get_session()
                        r = session.get(sdf_url, timeout=self.timeout_s)
                    if r.status_code == 200 and r.text:
                        mol = Chem.MolFromMolBlock(r.text, sanitize=True, removeHs=True) # Use RDKit to parse the SDF, which can sometimes succeed where PubChem's own SMILES generation fails. Remove hydrogens for better sanitisation.
                        if mol is not None:
                            smiles = Chem.MolToSmiles(mol, isomericSmiles=True) # Isomeric SMILES for better specificity
                except Exception:
                    pass

            # If STILL missing, use OPSIN fallback via IUPAC
            if (not smiles) and p0.get("IUPACName"):
                mol_opsin_iupac = self._resolve_opsin(p0["IUPACName"]) # Try OPSIN with the IUPAC name from PubChem
                if mol_opsin_iupac.smiles:
                    smiles = mol_opsin_iupac.smiles
                    resolution_source = "pubchem+opsin(iupac)" # Note that we got the IUPAC name from PubChem, but had to use OPSIN to resolve it to SMILES, so it's a hybrid resolution
                    resolution_note = "OK (OPSIN via IUPAC)" # :')

            return ResolvedMolecule( # Finally, return the resolved molecule with all the properties we got from PubChem, and the resolution source/note for transparency
                name=name,
                smiles=smiles,
                cid=cid,
                iupac_name=p0.get("IUPACName"),
                molecular_formula=p0.get("MolecularFormula"),
                pubchem_molwt=p0.get("MolecularWeight"),
                pubchem_xlogp=p0.get("XLogP"),
                resolution_source=resolution_source,
                resolution_note=resolution_note,
            )

        except Exception as e:
            return ResolvedMolecule(
                name=name,
                resolution_source="pubchem",
                resolution_note=f"Exception: {type(e).__name__}: {e}", # Catch any exceptions and return a failure record
            )
        
    def _resolve_opsin(self, name: str) -> ResolvedMolecule: # OPSIN is a great fallback for IUPAC names that PubChem might struggle with, and it can sometimes resolve names that PubChem can't find at all. It’s not perfect, but it’s worth trying if PubChem fails to give us a SMILES
        url = "https://opsin.ch.cam.ac.uk/opsin/" + requests.utils.quote(name) + ".json" # OPSIN name lookup endpoint in JSON 
        try:
            with self._opsin_sem:
                session = self._get_session()
                r = session.get(url, timeout=self.timeout_s)

            if r.status_code != 200:
                return ResolvedMolecule(
                    name=name,
                    resolution_source="opsin",
                    resolution_note=f"HTTP {r.status_code}", # Same as with PubChem
                )

            data = r.json()
            smiles = data.get("smiles") # OPSIN returns SMILES directly if it can resolve the name, but it doesn't give us as much additional info as PubChem

            return ResolvedMolecule(
                name=name,
                smiles=smiles,
                resolution_source="opsin", # Note that we don't get properties like CID or IUPAC name from OPSIN, so those fields will be None
                resolution_note="OK", # I guess :/
            )

        except Exception as e:
            return ResolvedMolecule(
                name=name,
                resolution_source="opsin",
                resolution_note=f"Exception: {type(e).__name__}: {e}", # Same as before
            )
        
    def _cache_allows_retry(self, cached_dict: dict) -> bool: # Determine if a failed cache entry is old enough to retry based on the failure timestamp and the configured TTL
        if not is_failure_record(cached_dict):
            return True

        failed_at = cached_dict.get("failed_at")
        if not failed_at:
            return True

        try:
            t = datetime.fromisoformat(failed_at.replace("Z", "+00:00")) # Parse the timestamp (Z is apparently +00:00)
        except Exception:
            return True

        age_days = (datetime.now(timezone.utc) - t).total_seconds() / 86400.0 # Calculate the age of the failure record in days
        return age_days >= self.failure_ttl_days    

# Descriptor calculation functions for RDKit and Mordred with family handling
def descriptor_in_family(desc_name: str, family_name: str) -> bool:
    if family_name == "functional_groups":
        return desc_name.startswith("fr_") # I don't wanna be French

    if family_name == "vsa_estate": # Turns out it's EState, not Estate, which means electrotopological state Van Der Waals surface area descriptors, and has nothing to do with housing :(
        prefixes = ("PEOE_VSA", "SMR_VSA", "SlogP_VSA", "EState_VSA", "VSA_EState") 
        return desc_name.startswith(prefixes)

    return desc_name in RDKIT_FAMILIES.get(family_name, set())

# But where are those good old fashioned family values? On which we used to rely! Lucky there's a family guy!A

# Build the set of enabled RDKit descriptors based on the configuration
def build_enabled_rdkit_descriptor_set(rdkit_cfg: dict) -> set[str]:
    enabled = set()

    families = rdkit_cfg.get("families", {}) # Get the families configuration
    for family_name, is_on in families.items():
        if not is_on:
            continue

        if family_name == "functional_groups":
            enabled.add("__FUNCTIONAL_GROUPS__")
        elif family_name == "vsa_estate":
            enabled.add("__VSA_ESTATE__")
        else:
            enabled.update(RDKIT_FAMILIES.get(family_name, set()))

    # Make sure to apply includes/excludes after families, so that they can override family selections if needed
    include_descriptors = set(rdkit_cfg.get("include_descriptors", []))
    exclude_descriptors = set(rdkit_cfg.get("exclude_descriptors", []))

    enabled.update(include_descriptors)
    enabled.difference_update(exclude_descriptors)

    return enabled

# Same but with mordred
def build_enabled_mordred_module_set(mordred_cfg: dict) -> set[str]:
    enabled = set()
    families = mordred_cfg.get("families", {})
    for family_name, is_on in families.items():
        if is_on:
            enabled.update(m.lower() for m in MORDRED_FAMILIES.get(family_name, set())) # mordred is a bit different you see
    return enabled

# DescriptorCalculator: Computes RDKit and mordred descriptors, basically 
class DescriptorCalculator:
    # Initialisation stuff
    def __init__(self, use_mordred: bool = False, structure_cfg: Optional[dict] = None): # You can choose to use mordred descriptors or not, and also pass in a configuration for which descriptors to calculate
        self.cfg = structure_cfg or DEFAULT_STRUCTURE_CONFIG # Hey if you don't provide a config I'm just gonna use the one I put as default

        # Get the RDKit and Mordred configurations from the main config
        rdkit_cfg = self.cfg.get("rdkit", {}) 
        mordred_cfg = self.cfg.get("mordred", {})

        # Determine RDKit descriptors
        self.rdkit_enabled = rdkit_cfg.get("enabled", True)
        self.enabled_rdkit_descriptors = build_enabled_rdkit_descriptor_set(rdkit_cfg)
        self.enabled_rdkit_descriptors_lower = {d.lower() for d in self.enabled_rdkit_descriptors}

        # Determine Mordred descriptors
        self.mordred_enabled = mordred_cfg.get("enabled", False)
        self.mordred_use_all = mordred_cfg.get("use_all", False)
        self.enabled_mordred_modules = build_enabled_mordred_module_set(mordred_cfg)
        self.mordred_include = {d.lower() for d in mordred_cfg.get("include_descriptors", [])}
        self.mordred_exclude = {d.lower() for d in mordred_cfg.get("exclude_descriptors", [])}

        self.use_mordred = (use_mordred or self.mordred_enabled) and HAS_MORDRED # Only use mordred if it's enabled in the config and we have it installed
        self.mordred_calc = None

        if self.use_mordred:
            self.mordred_calc = Calculator(mordred_descriptors, ignore_3D=True) # Be my guest if you want 3D desciptors, but I'm not implementing 3D conformer generation. I have 2 weeks to do this, remember?

        self.rdkit_descriptor_fns = [] # List of (name, function) tuples for the RDKit descriptors we want to calculate
        for desc_name, fn in Descriptors._descList:
            if self._rdkit_descriptor_is_enabled(desc_name):
                self.rdkit_descriptor_fns.append((desc_name, fn)) 

    # Is this thing on?
    def _rdkit_descriptor_is_enabled(self, desc_name: str) -> bool: 
        if not self.rdkit_enabled:
            return False # If no, then no (wow)

        if desc_name.lower() in self.enabled_rdkit_descriptors_lower: 
            return True # First check if the descriptor is explicitly included (case-insensitive)

        if "__FUNCTIONAL_GROUPS__" in self.enabled_rdkit_descriptors and desc_name.startswith("fr_"):
            return True # Then check if it's a functional group descriptor and the family is enabled

        if "__VSA_ESTATE__" in self.enabled_rdkit_descriptors:
            prefixes = ("PEOE_VSA", "SMR_VSA", "SlogP_VSA", "EState_VSA", "VSA_EState")
            if desc_name.startswith(prefixes):
                return True # Then check if it's a VSA/EState descriptor and the family is enabled

        return False # Otherwise, it's not enabled           

    # Calculate the selected RDKit descriptors for a given molecule and return them in a dictionary.
    def rdkit_selected_descriptors(self, mol: Chem.Mol) -> Dict[str, float]: 
        out: Dict[str, float] = {}
        # Loop through the selected RDKit descriptors and calculate their values for the given molecule
        for desc_name, fn in self.rdkit_descriptor_fns:
            try:
                value = fn(mol)
                out[f"rd_{desc_name}"] = float(value) if value is not None else float("nan") # Convert to float and handle None values as NaN 
            except Exception:
                out[f"rd_{desc_name}"] = float("nan") # If your molecule is weird, then you deserve NaN values
        return out
    
    # Determine if a given Mordred descriptor should be calculated based on the configuration
    def _mordred_descriptor_is_enabled(self, descriptor_obj) -> bool:
        if not self.use_mordred:
            return False

        desc_name = str(descriptor_obj)
        desc_name_l = desc_name.lower()

        module_name = descriptor_obj.__class__.__module__.split(".")[-1].lower() # Get family/module name (last part of the module path)

        if desc_name_l in self.mordred_exclude:
            return False # Nuh uh

        if desc_name_l in self.mordred_include:
            return True # Yes please

        if self.mordred_use_all:
            return True # Yes to everything, except for the ones I said no to

        return module_name in self.enabled_mordred_modules

    # Get Morgan's fingerprint 
    def _morgan_fp(self, mol: Chem.Mol) -> Dict[str, object]:
        rdkit_cfg = self.cfg.get("rdkit", {})
        fp_cfg = rdkit_cfg.get("fingerprints", {}).get("morgan", {}) # Extract Morgan and their fingerprint

        if not fp_cfg.get("enabled", False):
            return {}

        radius = int(fp_cfg.get("radius", 2)) # Take their fingerprint away from them and give it a radius of 2
        nbits = int(fp_cfg.get("nbits", 2048)) # And make the fingerprint 2048 bits long
        fmt = str(fp_cfg.get("format", "bitstring")).lower() # And put it in lowercase bitstring format

        gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits) 
        bv = gen.GetFingerprint(mol)

        if fmt == "expanded": # Secret option to make your Excel spreadsheets another 2048 columns wide!
            return {f"fp_morgan_{i}": int(bv.GetBit(i)) for i in range(nbits)} 

        return {"fp_morgan_bitstring": bv.ToBitString()} # But most people don't want that, so the bitstring will be compact
    # On a real note though: I chose radius 2 and 2048 bits because it's common, you can change if you feel like it
    
    # Molecular ACCess System keys (MACCS keys), used to determine molecule substrcture in 166 bits
    # Each bit means something, check here: https://github.com/rdkit/rdkit/blob/a9d6642f7dd732ee05814937f0279a05a7485754/rdkit/Chem/MACCSkeys.py#L39
    def _maccs(self, mol: Chem.Mol) -> Dict[str, object]:
        rdkit_cfg = self.cfg.get("rdkit", {})
        fp_cfg = rdkit_cfg.get("fingerprints", {}).get("maccs", {}) # Extract MACCS and its fingerprint

        if not fp_cfg.get("enabled", False):
            return {} # -

        fmt = str(fp_cfg.get("format", "bitstring")).lower()
        bv = MACCSkeys.GenMACCSKeys(mol) # Generate the MACCS keys fingerprint for the molecule using RDKit's built-in function

        if fmt == "expanded": # If you really don't want to write code to extract the bit(s) you want for analysis and can only write code to take certain columns out, I guess
            return {f"fp_maccs_{i}": int(bv.GetBit(i)) for i in range(bv.GetNumBits())} 

        return {"fp_maccs_bitstring": bv.ToBitString()}

    # Mordred descriptors
    def mordred_descriptors(self, mol: Chem.Mol) -> Dict[str, float]: 
        if not self.use_mordred or self.mordred_calc is None:
            return {} 

        out: Dict[str, float] = {} # Calculate selected Mordred descriptors and return them here

        try:
            vals = self.mordred_calc(mol) # Calculate all Mordred descriptors for the molecule 
            for descriptor_obj, value in vals.items(): # Loop through the calculated descriptors and their values
                if not self._mordred_descriptor_is_enabled(descriptor_obj):
                    continue 

                key = f"md_{descriptor_obj}" # Turn descriptor name into string and prefix with "md_"
                # Some Mordred descriptors clash with RDKit ones
                try:
                    out[key] = float(value) # Convert to float
                except Exception:
                    out[key] = float("nan") # Otherwise NaN

        except Exception:
            out["md_status_failed"] = 1.0 # womp womp

        return out

    # This is the code that actually takes the SMILES string and calculates all the descriptors we want
    def featurise_smiles(self, smiles: str) -> Tuple[Optional[Dict[str, object]], str]:
        mol = Chem.MolFromSmiles(smiles) # Parse the SMILES string into an RDKit molecule object
        if mol is None:
            return None, "RDKit failed to parse SMILES" # :(

        feats: Dict[str, object] = {} # This will hold all the calculated features

        if self.rdkit_enabled:
            feats.update(self.rdkit_selected_descriptors(mol)) # Calculate the selected RDKit descriptors and add them to the features dictionary

        # If Morgan or MACCS fingerprints are enabled, calculate and add them to the features as well
        feats.update(self._morgan_fp(mol))
        feats.update(self._maccs(mol))
        feats.update(self.mordred_descriptors(mol))

        return feats, "OK" # OK

# Halfway through this project I was asked to add more descriptors, but RDKit and Mordred doesn't have them
# So this class is basically stapled on so we can get some classic properties, like boiling point, density, etc. These will always be useful for basically anything
class ExternalPropertyEnricher:
    def __init__(self, cfg: Optional[dict] = None, headers: Optional[dict] = None, timeout_s: float = 20.0): # You can pass in a configuration for which external properties to get
        self.cfg = cfg or {} # But if you don't provide one, I'll just use an empty one
        self.headers = headers or {} # See above
        self.timeout_s = float(timeout_s)

    def enrich(self, resolved: ResolvedMolecule) -> Dict[str, Any]: # Add the external properties here and store them
        out: Dict[str, Any] = {}

        ext_cfg = self.cfg.get("external_properties", {})
        if not ext_cfg.get("enabled", False):
            return out

        # PubChem REST API properties (XlogP, MolWt, etc.)
        rest_cfg = ext_cfg.get("pubchem_rest", {})
        if rest_cfg.get("enabled", False):
            props = rest_cfg.get("properties", [])
            if props:
                out.update(
                    pubchem_enrich_by_cid( # Then pass them on
                        cid=resolved.cid,
                        properties=props,
                        headers=self.headers,
                        timeout_s=self.timeout_s,
                    )
                )

        # PubChem PUG View properties (boiling point, density, vapour pressure, water solubility)
        pv_cfg = ext_cfg.get("pubchem_pug_view", {})
        if pv_cfg.get("enabled", False):
            out.update(
                extract_pubchem_pugview_properties(
                    cid=resolved.cid,
                    headers=self.headers,
                    timeout_s=self.timeout_s,
                    extract_boiling_point=pv_cfg.get("extract_boiling_point", True),
                    extract_vapour_pressure=pv_cfg.get("extract_vapour_pressure", True),
                    extract_density=pv_cfg.get("extract_density", True),
                    extract_water_solubility=pv_cfg.get("extract_water_solubility", True),
                )
            )

        # Thermo properties from the Chemical module (boiling point, density, vapour pressure, viscosity, heat of vapourisation)
        thermo_cfg = ext_cfg.get("thermo", {})
        if thermo_cfg.get("enabled", False):
            out.update(
                self._thermo_properties(
                    resolved=resolved,
                    T=float(thermo_cfg.get("temperature_K", 298.15)), # STP by default but you can change them
                    P=float(thermo_cfg.get("pressure_Pa", 101325)),
                )
            )
        return out
    # Get properties from thermo
    def _thermo_properties(self, resolved: ResolvedMolecule, T: float = 298.15, P: float = 101325) -> Dict[str, Any]:
        out: Dict[str, Any] = {}

        if not HAS_THERMO:
            out["prop_thermo_status"] = "thermo_not_installed" # pip install thermo :)
            return out

        # Prefer CAS, then name, then SMILES-ish identifiers only if needed
        candidate_ids: List[str] = []
        if resolved.name:
            candidate_ids.append(resolved.name)
        if resolved.iupac_name and resolved.iupac_name not in candidate_ids:
            candidate_ids.append(resolved.iupac_name)

        chem = None 
        last_err = None

        for ident in candidate_ids:
            try:
                chem = Chemical(ident, T=T, P=P) # Try to create a Chemical object with the identifier, temperature, and pressure
                out["prop_thermo_identifier_used"] = ident # Store which identifier worked 
                break
            except Exception as e:
                last_err = e 

        if chem is None:
            out["prop_thermo_status"] = f"lookup_failed_{type(last_err).__name__}" if last_err else "lookup_failed" # :(
            return out

        try:
            # What kind of property is this?
            # This is the boiling point in Kelvin, which is the temperature where the chemical transitions from liquid to gas. 
                out["prop_thermo_boiling_point_K"] = float(chem.Tb)
                out["prop_thermo_boiling_point_C"] = float(chem.Tb - 273.15) # Conversion: Kelvin to Celsius
        except Exception:
            pass

        try:
            # What kind of property is this?
            # This is the density in kg/m^3, which is the mass per unit volume of the chemical. Useful to know how heavy something is.
            if getattr(chem, "rho", None) is not None:
                rho = float(chem.rho)
                out["prop_thermo_density_kg_m3"] = rho
                out["prop_thermo_density_g_ml"] = rho / 1000.0 # Conversion: kg/m^3 to g/mL
        except Exception:
            pass

        try:
            # What kind of property is this?
            # This is the vapour pressure in Pascals, which is the pressure exerted by the chemical in its gaseous state when in equilibrium.
            Psat = getattr(chem, "Psat", None)
            if Psat is not None:
                psat_pa = float(Psat)
                out["prop_thermo_vapour_pressure_Pa"] = psat_pa
                out["prop_thermo_vapour_pressure_mmHg"] = psat_pa / 133.322368 # Conversion: Pascals to mmHg
                out["prop_thermo_vapour_pressure_T_K"] = T # Store temperature as vapour pressure is temperature dependent
        except Exception:
            pass

        try:
            # What kind of property is this?
            # This is the viscosity in Pascal-seconds, which is the resistance of the chemical to flow.
            mu = getattr(chem, "mu", None)
            if mu is not None:
                out["prop_thermo_viscosity_Pa_s"] = float(mu)
                out["prop_thermo_viscosity_cP"] = float(mu) * 1000.0 # Conversion: Pa.s to cP (centipoise)
                out["prop_thermo_viscosity_T_K"] = T # Store temperature as viscosity is temperature dependent
        except Exception:
            pass

        try:
            # What kind of property is this?
            # This is the heat of vapourisation in J/mol. It’s the amount of energy required to convert one mole of the chemical from liquid to gas at constant temperature and pressure.
            Hvapm = getattr(chem, "Hvapm", None)
            if Hvapm is not None:
                out["prop_thermo_heat_vapourisation_J_mol"] = float(Hvapm)
                out["prop_thermo_heat_vapourisation_kJ_mol"] = float(Hvapm) / 1000.0 # Conversion: J/mol to kJ/mol
        except Exception:
            pass

        out["prop_thermo_status"] = "ok" # :)
        return out

# Normalise these properties
# Once you've experienced computational chemistry, you'll never stop wanting to beat whoever decided not to standardise these databases with base SI units to death with your bare hands
def normalise_external_properties(block: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}

    # 1) Boiling point (in Celsius, more intuitive). Priority: thermo > PubChem PUG View > missing
    thermo_bp = block.get("prop_thermo_boiling_point_C")
    if thermo_bp is not None:
        out["prop_boiling_point_C"] = float(thermo_bp)
    else:
        raw = block.get("prop_pugview_boiling_point_raw")
        if raw:
            raw = str(raw)
            # Ok so this looks like a mess, but basically we're trying to find some semblance of a number followed by a unit that has a C in it
            # We're dealing with people who write 37 degrees Celsius as: 37 °C, 37 C, 37 degrees C, 37 degC, 37C, 37.0 °C, etc.
            m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*°?\s*C", raw, re.I) 
            if m:
                out["prop_boiling_point_C"] = float(m.group(1))
            else:
                m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*K", raw, re.I) # Same comment above but for Kelvin
                if m:
                    out["prop_boiling_point_C"] = float(m.group(1)) - 273.15

    # 2) Density [kg/m^3]. Priority: thermo > PubChem PUG View > missing
    thermo_rho = block.get("prop_thermo_density_kg_m3")
    if thermo_rho is not None:
        out["prop_density_kg_m3"] = float(thermo_rho)
    else:
        raw = block.get("prop_pugview_density_raw")
        if raw:
            raw = str(raw)

            # g/cm3 -> kg/m3
            m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*g\s*/?\s*cm3", raw, re.I)
            if m:
                out["prop_density_kg_m3"] = float(m.group(1)) * 1000.0
            else:
                # kg/m3 already
                m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*kg\s*/?\s*m3", raw, re.I)
                if m:
                    out["prop_density_kg_m3"] = float(m.group(1))
                else:
                    # specific gravity is dimensionless, approximate as * 1000 kg/m3
                    m = re.search(r"(?:specific gravity|sp\.?\s*gr\.?)\s*[:=]?\s*([-+]?[0-9]*\.?[0-9]+)", raw, re.I)
                    if m:
                        out["prop_density_kg_m3"] = float(m.group(1)) * 1000.0

    
    # 3) Vapour pressure [Pa]. Priority: thermo > PubChem PUG View > missing
    thermo_vp = block.get("prop_thermo_vapour_pressure_Pa")
    if thermo_vp is not None:
        out["prop_vapour_pressure_Pa"] = float(thermo_vp)
    else:
        raw = block.get("prop_pugview_vapour_pressure_raw")
        if raw:
            raw = str(raw)
            m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*kPa\b", raw, re.I) # KPa
            if m:
                out["prop_vapour_pressure_Pa"] = float(m.group(1)) * 1000.0 # KPa to Pa
            else:
                m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*Pa\b", raw, re.I) # Pa
                if m:
                    out["prop_vapour_pressure_Pa"] = float(m.group(1)) # Already in Pa
                else:
                    m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*mmHg\b", raw, re.I) # mmHg
                    if m:
                        out["prop_vapour_pressure_Pa"] = float(m.group(1)) * 133.322 # mmHg to Pa

    # 4) Viscosity [Pa.s]. # Priority: thermo only (I gave up with PubChem PUG View, values are just too inconsistent and hard to parse)
    thermo_mu = block.get("prop_thermo_viscosity_Pa_s")
    if thermo_mu is not None:
        out["prop_viscosity_Pa_s"] = float(thermo_mu)

    # 5) Heat of vapourisation [kJ/mol]. # Priority: thermo only
    thermo_hvap = block.get("prop_thermo_heat_vapourisation_kJ_mol")
    if thermo_hvap is not None:
        out["prop_heat_vapourisation_kJ_mol"] = float(thermo_hvap)

    # 6) Water solubility [g/L]. # Priority: PubChem PUG View only for now
    raw = block.get("prop_pugview_water_solubility_raw")
    if raw:
        raw = str(raw)

        # g/L
        m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*g\s*/?\s*L", raw, re.I)
        if m:
            out["prop_water_solubility_g_L"] = float(m.group(1))
        else:
            # mg/L -> g/L
            m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*mg\s*/?\s*L", raw, re.I)
            if m:
                out["prop_water_solubility_g_L"] = float(m.group(1)) / 1000.0
            else:
                # mg/mL -> g/L
                m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*mg\s*/?\s*mL", raw, re.I)
                if m:
                    out["prop_water_solubility_g_L"] = float(m.group(1))
                else:
                    # ug/mL or µg/mL -> g/L
                    m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*(?:ug|µg)\s*/?\s*mL", raw, re.I)
                    if m:
                        out["prop_water_solubility_g_L"] = float(m.group(1)) / 1000.0
                    else:
                        # g/100 mL -> g/L
                        m = re.search(r"([-+]?[0-9]*\.?[0-9]+)\s*g\s*/?\s*100\s*mL", raw, re.I)
                        if m:
                            out["prop_water_solubility_g_L"] = float(m.group(1)) * 10.0

    return out # It's because of the awful inconsistency that 4, 5, and 6 have lower priority and can be very sparse in datasets. 

# THE main Pipeline! 
# This is where we take the input DataFrame, resolve the molecules, calculate descriptors, enrich with external properties, and return the final DataFrame with all the features
class FeaturePipeline:
    def __init__( 
        self,
        resolver,
        calculator,
        prefer_input_smiles=True, # If True, use the input SMILES, then fall back to name and CAS
        properties_cfg: Optional[dict] = None, # Configuration for which properties to enrich from external sources
        external_enricher: Optional[ExternalPropertyEnricher] = None, # Optional external property enricher 
    ):
        # Resolver, calculator, enricher config
        self.resolver = resolver 
        self.calculator = calculator
        self.prefer_input_smiles = prefer_input_smiles
        self.properties_cfg = properties_cfg or {"enabled": False}
        self.external_enricher = external_enricher

    def run(self, df_in: pd.DataFrame) -> pd.DataFrame:
        if "name" not in df_in.columns:
            raise ValueError("Input must contain a 'name' column.") # We need at least a name column to have something to resolve
        if "smiles" not in df_in.columns:
            df_in = df_in.copy()
            df_in["smiles"] = "" # If no SMILES column, add an empty one so we can handle it uniformly
        if "cas" not in df_in.columns:
            df_in = df_in.copy()
            df_in["cas"] = "" # Same for CAS numbers, we can handle them uniformly as well

        df_in = df_in.reset_index(drop=True) # Ensure the index is a simple range from 0 to n-1 for easier parallel processing

        # Process one row of the input DataFrame: resolve the molecule, calculate descriptors, enrich with external properties, and return a dictionary of all the features for that row
        def process_one(idx: int) -> tuple[int, Dict[str, Any]]:
            try:
                row = df_in.loc[idx]
                # Clean input text
                name = clean_text(row.get("name", ""))
                cas = clean_text(row.get("cas", ""))
                provided_smiles = clean_text(row.get("smiles", ""))
                resolved: Optional[ResolvedMolecule] = None

                if self.prefer_input_smiles and provided_smiles:
                    mol_test = Chem.MolFromSmiles(provided_smiles) # Use SMILES if valid

                    if mol_test is not None:
                        resolved = ResolvedMolecule(
                            name=name,
                            smiles=provided_smiles,
                            resolution_source="input",
                            resolution_note="User-provided SMILES",
                        )
                    elif cas: # If SMILES is invalid, fall back to CAS
                        resolved = self.resolver.resolve(cas)
                        resolved.name = name
                        resolved.resolution_note = (
                            f"Input SMILES invalid; from CAS: {cas} | "
                            + (resolved.resolution_note or "")
                        )
                    elif name: # If no CAS, fall back to name
                        resolved = self.resolver.resolve(name)
                        resolved.name = name
                        resolved.resolution_note = (
                            "Input SMILES invalid; " + (resolved.resolution_note or "")
                        )
                    else: # This means you messed up the input
                        resolved = ResolvedMolecule(
                            name="",
                            smiles=None,
                            resolution_source="input",
                            resolution_note="Input SMILES invalid and no fallback identifier available",
                        )
                # If not preferring input SMILES, try CAS first
                elif cas: 
                    resolved = self.resolver.resolve(cas)
                    resolved.name = name
                    resolved.resolution_note = f"From CAS: {cas} | " + (resolved.resolution_note or "")

                elif name: # Then try the name
                    resolved = self.resolver.resolve(name)

                else: # Otherwise, you really messed up the input. What did you give me?
                    resolved = ResolvedMolecule(
                        name="",
                        smiles=None,
                        resolution_source="none",
                        resolution_note="No valid name, CAS, or SMILES provided",
                    )
                # Build the output dictionary
                pubchem_block: Dict[str, Any] = asdict(resolved)
                rdkit_block: Dict[str, Any] = {}
                mordred_block: Dict[str, Any] = {}
                external_block: Dict[str, Any] = {}

                if resolved.smiles:
                    feats, status = self.calculator.featurise_smiles(resolved.smiles) # Calculate descriptors and get status
                    pubchem_block["feature_status"] = status 

                    # Distribute the calculated features into the appropriate blocks based on their prefixes
                    if feats:
                        for k, v in feats.items():
                            if k.startswith("rd_") or k.startswith("fp_"): # rd_ = RDKit descriptor, fp_ = fingerprint
                                rdkit_block[k] = v
                            elif k.startswith("md_") or k.startswith("mordred_"): # md_ = Mordred descriptor
                                mordred_block[k] = v
                            else:
                                rdkit_block[k] = v
                else:
                    pubchem_block["feature_status"] = "No SMILES" # :(

                # Enrich with external properties if enabled and we have a resolved molecule with a CID
                if self.external_enricher is not None: 
                    raw_external = self.external_enricher.enrich(resolved)
                    external_block.update(normalise_external_properties(raw_external))

                # Update the base dictionary with all the blocks of features
                base: Dict[str, Any] = {}
                base.update(pubchem_block)
                base.update(rdkit_block)
                base.update(mordred_block)
                base.update(external_block)

                return idx, base

            # Something went wrong during processing this row- catch the exception, return dictionary with error info
            except Exception as e:
                row = df_in.loc[idx]
                return idx, {
                    "name": clean_text(row.get("name", "")),
                    "smiles": clean_text(row.get("smiles", "")),
                    "feature_status": f"row_failed: {type(e).__name__}",
                    "resolution_note": f"Row processing exception: {e}",
                }

        results: List[Optional[Dict[str, Any]]] = [None] * len(df_in) # Pre-allocate results list for parallel processing

        # Parallel processing is difficult, but I tried my best here, given these datasets will be used for ML
        with ThreadPoolExecutor(max_workers=4) as ex: # I've set 4 here, adjust depending on how potato your computer is
            futures = [ex.submit(process_one, i) for i in range(len(df_in))]

            iterator = as_completed(futures)
            if HAS_TQDM: # Progress bar, so the user has something to look at while they wait around
                iterator = tqdm(
                    iterator,
                    total=len(futures),
                    desc="Processing entries",
                    unit="entry"
                )

            for f in iterator:
                idx, out = f.result()
                results[idx] = out

        return pd.DataFrame(results)

# Utility functions and helpers
# All that stuff we did above? It needs to be put in an Excel spreadsheet for the user, so this stuff makes it look nice at least
def read_input(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in [".xlsx", ".xls"]: # If it's an Excel file, read with pandas' read_excel
        return pd.read_excel(path)
    return pd.read_csv(path) # Otherwise, assume it's a CSV file and read with read_csv

# Standardisation of column labels for better readability in output tables and Excel sheets
COLUMN_LABEL_OVERRIDES = {
    "name": "Name",
    "smiles": "SMILES",
    "cid": "CID",
    "iupac_name": "IUPAC Name",
    "molecular_formula": "Molecular Formula",
    "pubchem_xlogp": "PubChem XLogP",
    "pubchem_molwt": "PubChem Molecular Weight",
    "resolution_source": "Resolution Source",
    "resolution_note": "Resolution Note",
    "feature_status": "Feature Status",

    "prop_boiling_point_C": "Boiling Point (°C)",
    "prop_density_kg_m3": "Density (kg/m³)",
    "prop_vapour_pressure_Pa": "Vapour Pressure (Pa)",
    "prop_viscosity_Pa_s": "Viscosity (Pa·s)",
    "prop_heat_vapourisation_kJ_mol": "Heat of Vapourisation (kJ/mol)",
    "prop_water_solubility_g_L": "Water Solubility (g/L)",
}

# Build the summary table for the overview
def build_summary_tables(df: pd.DataFrame, runtime_s: Optional[float] = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    n_rows = len(df)
    if "feature_status" in df.columns:
        resolved_mask = df["feature_status"] == "OK" # Consider a row resolved if its feature_status is "OK"
        resolved_rows = int(resolved_mask.sum()) # Count how many rows are resolved
        unresolved_rows = int((~resolved_mask).sum()) # Count how many rows are unresolved (not "OK")
    else:
        resolved_rows = n_rows
        unresolved_rows = 0 # If no feature_status column, assume all rows are resolved

    total_cols = len(df.columns) # Total number of columns in the DataFrame
    # Categorise columns based on prefixes to determine how many are from RDKit, Mordred, external properties, or just metadata
    rdkit_cols = [c for c in df.columns if c.startswith("rd_") or c.startswith("fp_")] 
    mordred_cols = [c for c in df.columns if c.startswith("md_") or c.startswith("mordred_")]
    external_cols = [c for c in df.columns if c.startswith("prop_") and not c.startswith("prop_pubchem_")]
    metadata_cols = [c for c in df.columns if c not in rdkit_cols + mordred_cols + external_cols]

    # Build this information DataFrame with the useful info
    overview_rows = [
        ("Number of Input Rows", n_rows),
        ("Resolved Rows", resolved_rows),
        ("Unresolved Rows", unresolved_rows),
        ("Resolution Rate (%)", round(100.0 * resolved_rows / n_rows, 2) if n_rows else 0.0),
        ("Total Columns", total_cols),
        ("Metadata Columns", len(metadata_cols)),
        ("RDKit Columns", len(rdkit_cols)),
        ("Mordred Columns", len(mordred_cols)),
        ("External Property Columns", len(external_cols)),
    ]
    # Runtime calculation
    if runtime_s is not None:
        minutes = int(runtime_s // 60)
        seconds = int(runtime_s % 60)

        runtime_str = f"{minutes} m {seconds} s"

        overview_rows.extend([
            ("Runtime", runtime_str),
            ("Seconds per Row", round(runtime_s / n_rows, 4) if n_rows else None),
        ])

    overview_df = pd.DataFrame(overview_rows, columns=["Metric", "Value"])
    
    # Build the coverage table for each column
    coverage_rows = []
    for col in df.columns:
        non_missing = int(df[col].notna().sum())
        missing = int(df[col].isna().sum())
        availability_pct = round(100.0 * non_missing / n_rows, 2) if n_rows else 0.0 # Availability calculation

        # Determine the group for this column based on its prefix
        if col in rdkit_cols:
            group = "rdkit"
        elif col in mordred_cols:
            group = "mordred"
        elif col in external_cols:
            group = "external"
        else:
            group = "metadata"

        # Append the coverage information for this column to the list of coverage rows
        coverage_rows.append({
            "column": col,
            "group": group,
            "non_missing": non_missing,
            "missing": missing,
            "availability_percent": availability_pct,
        })

    # Create DataFrame
    coverage_df = pd.DataFrame(coverage_rows)
    coverage_df = coverage_df.sort_values(
        by=["group", "availability_percent", "column"],
        ascending=[True, False, True]
    ).reset_index(drop=True)

    # Then apply the pretty column names for the coverage table
    coverage_df = coverage_df.rename(columns={
        "column": "Column",
        "group": "Group",
        "non_missing": "Non-Missing",
        "missing": "Missing",
        "availability_percent": "Availability (%)",
    })

    return overview_df, coverage_df

# Function that writes the Excel spreadsheets and splits the sheets up according to the feature types
def write_excel(df: pd.DataFrame, out_path: Path, runtime_s: Optional[float] = None) -> None:
    unresolved = df[df["feature_status"] != "OK"].copy() # Separate out unresolved rows for a dedicated sheet

    # Keep raw/internal names in exported data sheets
    df_pretty = df.rename(columns=COLUMN_LABEL_OVERRIDES).copy()
    unresolved_pretty = unresolved.rename(columns=COLUMN_LABEL_OVERRIDES).copy()

    # For the main sheet, keep all columns
    pubchem_cols = [
        "name", "smiles", "cid", "iupac_name", "molecular_formula",
        "pubchem_xlogp", "pubchem_molwt", "resolution_source",
        "resolution_note", "prop_pubchem_MolecularWeight",
        "prop_pubchem_XLogP", "prop_pubchem_MolecularFormula",
        "prop_pubchem_IUPACName", "prop_pubchem_status", "feature_status"
    ]
    pubchem_cols = [c for c in pubchem_cols if c in df.columns] # Only include PubChem columns that are actually present in the DataFrame
    rdkit_cols = [c for c in df.columns if c.startswith("rd_") or c.startswith("fp_")] # Same but with RDKit 
    mordred_cols = [c for c in df.columns if c.startswith("md_") or c.startswith("mordred_")] # Same but with Mordred
    external_cols = [c for c in df.columns if c.startswith("prop_") and not c.startswith("prop_pubchem_")] # Same but with external properties (excluding PubChem REST API properties which are already in pubchem_cols)

    overview_df, coverage_df = build_summary_tables(df, runtime_s=runtime_s) # Build the summary tables for the overview sheet

    # For the other sheets, split them up by feature type and apply pretty column names
    pubchem_rdkit_pretty = df[pubchem_cols + rdkit_cols].rename(columns=COLUMN_LABEL_OVERRIDES).copy() if (pubchem_cols + rdkit_cols) else None
    pubchem_mordred_pretty = df[pubchem_cols + mordred_cols].rename(columns=COLUMN_LABEL_OVERRIDES).copy() if (pubchem_cols + mordred_cols) else None
    pubchem_external_pretty = df[pubchem_cols + external_cols].rename(columns=COLUMN_LABEL_OVERRIDES).copy() if (pubchem_cols + external_cols) else None

    # Use openpyxl to write the Excel file with multiple sheets and auto-fit column widths
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        # Summary
        overview_df.to_excel(writer, index=False, sheet_name="Summary", startrow=0)
        coverage_df.to_excel(writer, index=False, sheet_name="Summary", startrow=len(overview_df) + 3)

        # Main sheets
        df_pretty.to_excel(writer, index=False, sheet_name="All_Features")
        unresolved_pretty.to_excel(writer, index=False, sheet_name="Unresolved")
        if pubchem_rdkit_pretty is not None:
            pubchem_rdkit_pretty.to_excel(writer, index=False, sheet_name="PubChem_RDKit")
        if pubchem_mordred_pretty is not None:
            pubchem_mordred_pretty.to_excel(writer, index=False, sheet_name="PubChem_Mordred")
        if pubchem_external_pretty is not None:
            pubchem_external_pretty.to_excel(writer, index=False, sheet_name="PubChem_External")

        # Auto-fit column widths. Should work fine but depending on how you use this it may be problematic
        for sheet in writer.book.worksheets:
            for col_cells in sheet.columns:
                max_len = 0
                col_letter = col_cells[0].column_letter
                for cell in col_cells:
                    try:
                        if cell.value is not None:
                            max_len = max(max_len, len(str(cell.value)))
                    except Exception:
                        pass
                sheet.column_dimensions[col_letter].width = min(max_len + 2, 40) # Cap max width to 40 for readability

# Loading properties config from a JSON file
def load_properties_config(path: Path) -> dict:
    if not path or not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}

# Enrichment functions for external properties
def pubchem_enrich_by_cid(
    cid: Optional[int],
    properties: List[str],
    headers: dict,
    timeout_s: float,
    session: Optional[requests.Session] = None,
) -> Dict[str, object]:
    out: Dict[str, object] = {}
    if not cid:
        out["prop_pubchem_status"] = "no_cid" # If we don't have a CID, we can't query PubChem, so we return a status indicating that the CID is missing and skip the API call
        return out

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" # Base URL for PubChem PUG REST API to retrieve compound information by CID
        + str(cid)
        + "/property/"
        + ",".join(properties)
        + "/JSON" 
    )

    try:
        sess = session or requests.Session() # Use provided session or create a new one for making the API request
        if headers:
            sess.headers.update(headers)

        data, status = get_json_retry(sess, url, timeout_s) # Make the API request to PubChem and get the JSON response along with the HTTP status code

        if status != 200 or data is None: # API call failed, return status indicating the failure
            out["prop_pubchem_status"] = f"http_{status}"
            return out

        props_list = data.get("PropertyTable", {}).get("Properties", [])
        if not props_list:
            out["prop_pubchem_status"] = "no_properties" # API call succeeded but no properties were found in the response
            return out

        p0 = props_list[0] # We expect only one result since we're querying by CID, so we take the first item in the Properties list
        for k in properties:
            out[f"prop_pubchem_{k}"] = p0.get(k)

        out["prop_pubchem_status"] = "ok" # ok
        return out

    except Exception as e:
        out["prop_pubchem_status"] = f"exception_{type(e).__name__}" # Raise exception, something went wrong during API call
        return out

def pubchem_pugview_record( # Retries full PUG View record for a given CID
    cid: Optional[int],
    headers: dict,
    timeout_s: float,
    session: Optional[requests.Session] = None,
) -> Dict[str, Any]:
    if not cid:
        return {}

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON" # URL for PubChem PUG View API 

    try:
        sess = session or requests.Session() # Use provided session or create a new one for making the API request
        if headers:
            sess.headers.update(headers)

        data, status = get_json_retry(sess, url, timeout_s) # Make the API request to PubChem PUG View and get the JSON response along with the HTTP status code

        if status != 200 or data is None: # As above
            return {"pug_view_status": f"http_{status}"}

        return data

    except Exception as e:
        return {"pug_view_status": f"exception_{type(e).__name__}"}


def _iter_pugview_sections(node: Any): # Recursively iterate through the nested structure of a PubChem PUG View record to find all sections that contain a TOCHeading
    if isinstance(node, dict):
        if "TOCHeading" in node: # As above
            yield node
        for v in node.values():
            yield from _iter_pugview_sections(v)
    elif isinstance(node, list):
        for item in node:
            yield from _iter_pugview_sections(item)

# Extract the string values
def _extract_string_values_from_information(info_list: list) -> List[str]:
    out: List[str] = []

    for info in info_list or []:
        value = info.get("Value", {})

        # Common string forms
        if "StringWithMarkup" in value:
            for item in value["StringWithMarkup"]:
                s = item.get("String")
                if s:
                    out.append(s)

        if "String" in value:
            s = value.get("String")
            if s:
                out.append(s)

        # Sometimes numbers appear in Number / NumberUnit / etc.
        if "Number" in value:
            num = value.get("Number")
            if num is not None:
                unit = value.get("Unit", "")
                out.append(f"{num} {unit}".strip())

    return out

# Find the first section
def _find_first_section_texts_by_heading(pug_view_json: dict, headings: List[str]) -> List[str]:
    headings_l = {h.lower() for h in headings}

    for sec in _iter_pugview_sections(pug_view_json):
        heading = str(sec.get("TOCHeading", "")).strip().lower()
        if heading in headings_l:
            return _extract_string_values_from_information(sec.get("Information", []))

    return []

# Then all the properties from PUG View are taken
def extract_pubchem_pugview_properties(
    cid: Optional[int],
    headers: dict,
    timeout_s: float,
    extract_boiling_point: bool = True,
    extract_vapour_pressure: bool = True,
    extract_density: bool = True,
    extract_water_solubility: bool = True,
) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    if not cid:
        out["prop_pugview_status"] = "no_cid" # If no CID
        return out

    record = pubchem_pugview_record(cid=cid, headers=headers, timeout_s=timeout_s)
    if not record or "Record" not in record:
        out["prop_pugview_status"] = record.get("pug_view_status", "no_record") if isinstance(record, dict) else "no_record"
        return out

    # Extract boiling point, vapour pressure, density/specific gravity, solubility/water solubility
    if extract_boiling_point:
        vals = _find_first_section_texts_by_heading(record, ["Boiling Point"])
        if vals:
            out["prop_pugview_boiling_point_raw"] = " | ".join(vals)

    if extract_vapour_pressure:
        vals = _find_first_section_texts_by_heading(record, ["Vapour Pressure"])
        if vals:
            out["prop_pugview_vapour_pressure_raw"] = " | ".join(vals)

    if extract_density:
        vals = _find_first_section_texts_by_heading(record, ["Density", "Specific Gravity"])
        if vals:
            out["prop_pugview_density_raw"] = " | ".join(vals)

    if extract_water_solubility:
        vals = _find_first_section_texts_by_heading(record, ["Solubility", "Water Solubility"])
        if vals:
            out["prop_pugview_water_solubility_raw"] = " | ".join(vals)

    out["prop_pugview_status"] = "ok"
    return out

# Recursively merge two dictionaries, this is used to combine the default descriptor configuration with a user configuration file.
def deep_update(base: dict, updates: dict) -> dict:
    result = dict(base)
    for k, v in updates.items():
        if isinstance(v, dict) and isinstance(result.get(k), dict):
            result[k] = deep_update(result[k], v)
        else:
            result[k] = v
    return result

def main() -> None:
    t0 = time.perf_counter() # Begin time
    # Parse arguments into strings
    ap = argparse.ArgumentParser(
        description="Resolve fragrance/chemical names to SMILES and generate descriptors."
    )
    # CLI arguments. I like to work in the terminal, so these arguments can really help out
    ap.add_argument("--input", type=str, required=True, help="Input CSV/XLSX with columns: name, optional smiles")
    ap.add_argument("--output", type=str, default="features.xlsx", help="Output XLSX path")
    ap.add_argument("--cache", type=str, default="smiles_cache.json", help="Local SMILES cache JSON")
    ap.add_argument("--use-mordred", action="store_true", help="Add Mordred 2D descriptors (large feature set)")
    ap.add_argument("--no-prefer-input-smiles", action="store_true", help="Ignore provided smiles column and always resolve by name")
    ap.add_argument("--properties-config", type=str, default="properties_config.json", help="Optional JSON config for database enrichment")
    ap.add_argument("--clear-cache", action="store_true", help="Delete the SMILES cache before running")
    args = ap.parse_args()

    df_in = read_input(Path(args.input))
    cache_path = Path(args.cache)

    if args.clear_cache and cache_path.exists():
        cache_path.unlink()
        print("SMILES cache cleared.") # Clear cache every time you run

    # This part of the code runs the pipelines
    resolver = MoleculeResolver(cache_path=cache_path)
    props_cfg_raw = load_properties_config(Path(args.properties_config))
    props_cfg = deep_update(DEFAULT_STRUCTURE_CONFIG, props_cfg_raw)

    calc = DescriptorCalculator(use_mordred=args.use_mordred, structure_cfg=props_cfg)

    external_enricher = ExternalPropertyEnricher(
    cfg=props_cfg,
    headers=resolver.headers,
    timeout_s=resolver.timeout_s,
    )

    pipe = FeaturePipeline(
    resolver=resolver,
    calculator=calc,
    prefer_input_smiles=(not args.no_prefer_input_smiles),
    properties_cfg=props_cfg,
    external_enricher=external_enricher,
    )
    df_in = normalise_columns(df_in)

    df_out = pipe.run(df_in)

    if getattr(resolver, "_cache_dirty", False):
        resolver._save_cache()

    dt = time.perf_counter() - t0 # End timer here

    write_excel(df_out, Path(args.output), runtime_s=dt) # Add to Excel file with timing info

    # Print stuff to terminal indicating that it's done
    print(f"Saved: {args.output}")
    print(f"Rows: {len(df_out)} | Columns: {len(df_out.columns)}")
    print(f"Total runtime: {dt:.2f} s ({dt/60:.2f} min)")
    if len(df_out) > 0:
        print(f"Avg time per row: {dt/len(df_out):.3f} s")
    if args.use_mordred and not HAS_MORDRED:
        print("Note: --use-mordred was set but mordred import failed; only RDKit descriptors were produced.")

# Run main
if __name__ == "__main__":
    main()

# Code and comments written by Ryan. 
# If you are trying to debug this code, good luck and may the odds ever be in your favour.
