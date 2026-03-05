"""
BRENDA Database Query for GlgA (Glycogen Synthase) Kinetic Parameters
EC 2.4.1.21 - ADP-glucose:1,4-alpha-D-glucan 4-alpha-D-glucosyltransferase

Queries BRENDA SOAP API (https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl)
for Km and kcat values for GlgA, prioritising Lactobacillus acidophilus.

SOAP method signatures (from WSDL introspection):
  getKmValue(email, password, ecNumber, organism, kmValue,
             kmValueMaximum, substrate, commentary, ligandStructureId, literature)
  getTurnoverNumber(email, password, ecNumber, organism, turnoverNumber,
                    turnoverNumberMaximum, substrate, commentary,
                    ligandStructureId, literature)
  getKcatKmValue(email, password, ecNumber, organism, kcatKmValue,
                 kcatKmValueMaximum, substrate, commentary,
                 ligandStructureId, literature)
  getOrganism(email, password, ecNumber, organism, sequenceCode,
              commentary, literature, textmining)

Password must be SHA-256 hashed before passing to the API.
"""

import hashlib
import re
import os
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Credential loading
# ---------------------------------------------------------------------------

def load_credentials():
    """Load BRENDA credentials from environment variables or .env files."""
    email = (os.environ.get("BRENDA_EMAIL")
             or os.environ.get("BRENDA_EMIAL"))   # legacy typo variant
    password = os.environ.get("BRENDA_PASSWORD")

    env_search_paths = [
        Path(__file__).parent / ".env",
        Path.home() / ".claude" / "skills" / "brenda-database" / ".env",
        Path.home() / ".claude" / ".env",
    ]
    for env_path in env_search_paths:
        if env_path.exists():
            with open(env_path) as fh:
                for line in fh:
                    line = line.strip()
                    if "=" not in line or line.startswith("#"):
                        continue
                    key, _, val = line.partition("=")
                    val = val.strip().strip("\"'")
                    if key in ("BRENDA_EMAIL", "BRENDA_EMIAL") and not email:
                        email = val
                    elif key == "BRENDA_PASSWORD" and not password:
                        password = val
            if email and password:
                print(f"  Credentials loaded from: {env_path}")
                break

    return email, password


def hash_password(password: str) -> str:
    """Return SHA-256 hex digest of the password (BRENDA requirement)."""
    return hashlib.sha256(password.encode()).hexdigest()


# ---------------------------------------------------------------------------
# BRENDA SOAP client builder
# ---------------------------------------------------------------------------

def build_client():
    """Connect to BRENDA SOAP API and return zeep client."""
    from zeep import Client, Settings
    from zeep.transports import Transport
    import requests

    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    session = requests.Session()
    transport = Transport(session=session, timeout=60)
    settings = Settings(strict=False, xml_huge_tree=True)
    client = Client(wsdl, transport=transport, settings=settings)
    return client


# ---------------------------------------------------------------------------
# Data-parsing helpers
# ---------------------------------------------------------------------------

def parse_entry(obj) -> dict:
    """
    BRENDA SOAP returns either a string like
      'organism*Escherichia coli#substrate*glucose#kmValue*0.12#...'
    or a zeep object with direct attributes.  Handle both.
    """
    if obj is None:
        return {}

    # If it has __dict__ attributes (zeep object)
    if hasattr(obj, "__dict__"):
        raw = {k: v for k, v in vars(obj).items() if not k.startswith("_")}
        if raw:
            return raw

    # Fallback: string parsing
    text = str(obj)
    result = {}
    for part in text.split("#"):
        if "*" in part:
            key, _, value = part.partition("*")
            result[key.strip()] = value.strip()
    return result


def to_float(s) -> float | None:
    """Extract the first float from a string or numeric value."""
    if s is None:
        return None
    if isinstance(s, (int, float)):
        return float(s)
    m = re.search(r"(\d+\.?\d*(?:[eE][+-]?\d+)?)", str(s))
    return float(m.group(1)) if m else None


def extract_ph(text: str) -> float | None:
    if not text:
        return None
    m = re.search(r"pH\s*([0-9]+\.?[0-9]*)", text)
    return float(m.group(1)) if m else None


def extract_temp(text: str) -> float | None:
    if not text:
        return None
    m = re.search(r"(\d+)\s*[°\u00b0]?C\b", text)
    return float(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# Query functions (using correct individual-parameter calling convention)
# ---------------------------------------------------------------------------

WILDCARD = "*"


def query_km(client, email, pw_hash, ec, organism=WILDCARD, substrate=WILDCARD):
    """Retrieve Km values from BRENDA."""
    try:
        result = client.service.getKmValue(
            email, pw_hash,
            ec, organism, WILDCARD, WILDCARD,
            substrate, WILDCARD, WILDCARD, WILDCARD
        )
        if result is None:
            return []
        return list(result) if hasattr(result, "__iter__") and not isinstance(result, str) else [result]
    except Exception as exc:
        print(f"    [getKmValue] Error: {exc}")
        return []


def query_kcat(client, email, pw_hash, ec, organism=WILDCARD):
    """Retrieve turnover number (kcat) values from BRENDA."""
    try:
        result = client.service.getTurnoverNumber(
            email, pw_hash,
            ec, organism, WILDCARD, WILDCARD,
            WILDCARD, WILDCARD, WILDCARD, WILDCARD
        )
        if result is None:
            return []
        return list(result) if hasattr(result, "__iter__") and not isinstance(result, str) else [result]
    except Exception as exc:
        print(f"    [getTurnoverNumber] Error: {exc}")
        return []


def query_kcat_km(client, email, pw_hash, ec, organism=WILDCARD):
    """Retrieve kcat/Km catalytic efficiency values from BRENDA."""
    try:
        result = client.service.getKcatKmValue(
            email, pw_hash,
            ec, organism, WILDCARD, WILDCARD,
            WILDCARD, WILDCARD, WILDCARD, WILDCARD
        )
        if result is None:
            return []
        return list(result) if hasattr(result, "__iter__") and not isinstance(result, str) else [result]
    except Exception as exc:
        print(f"    [getKcatKmValue] Error: {exc}")
        return []


def query_organisms(client, email, pw_hash, ec):
    """Retrieve organisms with data for the given EC number."""
    try:
        result = client.service.getOrganism(
            email, pw_hash,
            ec, WILDCARD, WILDCARD, WILDCARD, WILDCARD, WILDCARD
        )
        if result is None:
            return []
        return list(result) if hasattr(result, "__iter__") and not isinstance(result, str) else [result]
    except Exception as exc:
        print(f"    [getOrganism] Error: {exc}")
        return []


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

LINE = "=" * 72


def section(title: str):
    print(f"\n{LINE}")
    print(f"  {title}")
    print(LINE)


def display_km(entries: list, header="Km results"):
    if not entries:
        print("  (no data returned)")
        return

    parsed_all = []
    for raw in entries:
        p = parse_entry(raw)
        if p:
            p["_km_num"] = to_float(p.get("kmValue") or p.get("kmvalue"))
            p["_ph"]     = extract_ph(str(p.get("commentary", "")))
            p["_temp"]   = extract_temp(str(p.get("commentary", "")))
            parsed_all.append(p)

    print(f"  Total entries : {len(parsed_all)}\n")

    # Group by organism
    by_org: dict = {}
    for p in parsed_all:
        org = str(p.get("organism", "unknown"))
        by_org.setdefault(org, []).append(p)

    for org in sorted(by_org):
        lacto_flag = " <-- LACTOBACILLUS" if "Lactobacillus" in org or "Limosilactobacillus" in org else ""
        print(f"  Organism: {org}{lacto_flag}")
        for p in by_org[org]:
            km_val  = p.get("kmValue") or p.get("kmvalue", "N/A")
            km_max  = p.get("kmValueMaximum") or p.get("kmvalueMaximum", "")
            substr  = p.get("substrate", "N/A")
            comment = str(p.get("commentary", ""))
            lit     = str(p.get("literature", ""))

            km_display = str(km_val)
            if km_max and str(km_max).strip():
                km_display += f" - {km_max}"

            print(f"    Substrate   : {substr}")
            print(f"    Km (mM)     : {km_display}")
            if p["_ph"] is not None:
                print(f"    pH          : {p['_ph']}")
            if p["_temp"] is not None:
                print(f"    Temperature : {p['_temp']} C")
            if comment.strip():
                print(f"    Conditions  : {comment[:130]}")
            if lit.strip():
                print(f"    Reference   : {lit[:90]}")
            print()


def display_kcat(entries: list):
    if not entries:
        print("  (no data returned)")
        return

    parsed_all = []
    for raw in entries:
        p = parse_entry(raw)
        if p:
            p["_kcat_num"] = to_float(p.get("turnoverNumber") or p.get("turnovernumber"))
            p["_ph"]       = extract_ph(str(p.get("commentary", "")))
            p["_temp"]     = extract_temp(str(p.get("commentary", "")))
            parsed_all.append(p)

    print(f"  Total entries : {len(parsed_all)}\n")

    by_org: dict = {}
    for p in parsed_all:
        org = str(p.get("organism", "unknown"))
        by_org.setdefault(org, []).append(p)

    for org in sorted(by_org):
        lacto_flag = " <-- LACTOBACILLUS" if "Lactobacillus" in org or "Limosilactobacillus" in org else ""
        print(f"  Organism: {org}{lacto_flag}")
        for p in by_org[org]:
            kcat_val  = p.get("turnoverNumber") or p.get("turnovernumber", "N/A")
            kcat_max  = p.get("turnoverNumberMaximum") or p.get("turnovernumberMaximum", "")
            substr    = p.get("substrate", "N/A")
            comment   = str(p.get("commentary", ""))
            lit       = str(p.get("literature", ""))

            kcat_display = str(kcat_val)
            if kcat_max and str(kcat_max).strip():
                kcat_display += f" - {kcat_max}"

            print(f"    Substrate   : {substr}")
            print(f"    kcat (1/s)  : {kcat_display}")
            if p["_ph"] is not None:
                print(f"    pH          : {p['_ph']}")
            if p["_temp"] is not None:
                print(f"    Temperature : {p['_temp']} C")
            if comment.strip():
                print(f"    Conditions  : {comment[:130]}")
            if lit.strip():
                print(f"    Reference   : {lit[:90]}")
            print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(LINE)
    print("  BRENDA Kinetic Parameter Query: GlgA (Glycogen Synthase)")
    print("  EC 2.4.1.21  |  ADP-glucose:1,4-alpha-D-glucan 4-alpha-glucosyltransferase")
    print("  Primary target: Lactobacillus acidophilus  |  Fallback: all organisms")
    print(LINE)

    # -----------------------------------------------------------------------
    # Credentials
    # -----------------------------------------------------------------------
    print("\n[1] Loading BRENDA credentials ...")
    email, password = load_credentials()

    if not email or not password:
        print("\n  WARNING: No BRENDA credentials found.")
        print("  Looked for:")
        print("    - Env vars BRENDA_EMAIL / BRENDA_PASSWORD")
        print("    - .env file in skill directory or ~/.claude/skills/brenda-database/")
        print()
        print("  To obtain credentials:")
        print("    1. Register at https://www.brenda-enzymes.org/")
        print("    2. Set BRENDA_EMAIL and BRENDA_PASSWORD environment variables")
        print("       OR create ~/.claude/skills/brenda-database/.env with those values")
        print()
        print("  Continuing with demo credentials to show API call structure ...")
        email    = "user@example.com"
        password = "demo_password_not_valid"
    else:
        print(f"  Email loaded: {email}")

    pw_hash = hash_password(password)
    print(f"  Password SHA-256: {pw_hash[:20]}...  (full hash passed to API)")

    # -----------------------------------------------------------------------
    # Connect
    # -----------------------------------------------------------------------
    print("\n[2] Connecting to BRENDA SOAP API ...")
    print("  WSDL: https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl")
    try:
        client = build_client()
        print("  Connection successful.")
    except Exception as e:
        print(f"  FATAL: Could not build SOAP client: {e}")
        sys.exit(1)

    EC_PRIMARY = "2.4.1.21"
    EC_ALT     = "2.4.1.11"
    TARGET_ORG = "Lactobacillus acidophilus"

    # -----------------------------------------------------------------------
    # Query 1: Km for Lactobacillus acidophilus
    # -----------------------------------------------------------------------
    section(f"Query 1 -- Km values: {TARGET_ORG} / EC {EC_PRIMARY}")
    print(f"  Calling getKmValue(ecNumber={EC_PRIMARY}, organism={TARGET_ORG}) ...")
    km_lacto = query_km(client, email, pw_hash, EC_PRIMARY, organism=TARGET_ORG)
    time.sleep(1)
    display_km(km_lacto)

    # -----------------------------------------------------------------------
    # Query 2: Km for all organisms
    # -----------------------------------------------------------------------
    section(f"Query 2 -- Km values: ALL organisms / EC {EC_PRIMARY}")
    print(f"  Calling getKmValue(ecNumber={EC_PRIMARY}, organism=*) ...")
    km_all = query_km(client, email, pw_hash, EC_PRIMARY)
    time.sleep(1)
    display_km(km_all)

    # Filter Lactobacillus within all-organism results
    lacto_subset = []
    for raw in km_all:
        if "Lactobacillus" in str(raw) or "Limosilactobacillus" in str(raw):
            p = parse_entry(raw)
            if p:
                lacto_subset.append(p)
    if lacto_subset:
        print(f"\n  Lactobacillus entries within all-organism Km results ({len(lacto_subset)}):\n")
        for p in lacto_subset:
            print(f"    Organism  : {p.get('organism', 'N/A')}")
            print(f"    Substrate : {p.get('substrate', 'N/A')}")
            print(f"    Km (mM)   : {p.get('kmValue', 'N/A')}")
            print(f"    Conditions: {str(p.get('commentary', ''))[:100]}")
            print()
    else:
        print("\n  No Lactobacillus-genus entries found in all-organism Km results.")

    # -----------------------------------------------------------------------
    # Query 3: kcat (turnover number) for all organisms
    # -----------------------------------------------------------------------
    section(f"Query 3 -- kcat (turnover number): ALL organisms / EC {EC_PRIMARY}")
    print(f"  Calling getTurnoverNumber(ecNumber={EC_PRIMARY}, organism=*) ...")
    kcat_all = query_kcat(client, email, pw_hash, EC_PRIMARY)
    time.sleep(1)
    display_kcat(kcat_all)

    # -----------------------------------------------------------------------
    # Query 4: kcat/Km catalytic efficiency
    # -----------------------------------------------------------------------
    section(f"Query 4 -- kcat/Km (catalytic efficiency): EC {EC_PRIMARY}")
    print(f"  Calling getKcatKmValue(ecNumber={EC_PRIMARY}, organism=*) ...")
    eff_all = query_kcat_km(client, email, pw_hash, EC_PRIMARY)
    time.sleep(1)
    if not eff_all:
        print("  (no data returned)")
    else:
        print(f"  Total entries : {len(eff_all)}\n")
        for raw in eff_all:
            p = parse_entry(raw)
            org    = p.get("organism", "N/A")
            substr = p.get("substrate", "N/A")
            val    = p.get("kcatKmValue", p.get("kcatkmvalue", "N/A"))
            val_mx = p.get("kcatKmValueMaximum", "")
            comment = str(p.get("commentary", ""))
            display_val = str(val) + (f" - {val_mx}" if val_mx and str(val_mx).strip() else "")
            print(f"  Organism     : {org}")
            print(f"  Substrate    : {substr}")
            print(f"  kcat/Km (mM-1 s-1): {display_val}")
            if comment.strip():
                print(f"  Conditions   : {comment[:110]}")
            print()

    # -----------------------------------------------------------------------
    # Query 5: Organism list
    # -----------------------------------------------------------------------
    section(f"Query 5 -- All organisms with EC {EC_PRIMARY} in BRENDA")
    print(f"  Calling getOrganism(ecNumber={EC_PRIMARY}) ...")
    org_entries = query_organisms(client, email, pw_hash, EC_PRIMARY)
    time.sleep(1)
    if not org_entries:
        print("  (no organism data returned)")
    else:
        org_names = set()
        for raw in org_entries:
            p = parse_entry(raw)
            name = p.get("organism") or str(raw)
            if name:
                org_names.add(str(name))
        print(f"  {len(org_names)} distinct organisms:\n")
        for name in sorted(org_names):
            flag = ""
            if "Lactobacillus" in name or "Limosilactobacillus" in name or "Lacticaseibacillus" in name:
                flag = " <-- TARGET GENUS"
            print(f"    {name}{flag}")

    # -----------------------------------------------------------------------
    # Fallback: EC 2.4.1.11 if primary yielded nothing
    # -----------------------------------------------------------------------
    if not km_all and not kcat_all:
        section(f"FALLBACK -- EC {EC_ALT} (UDP-glucose-dependent glycogen/starch synthase)")
        print(f"  EC {EC_PRIMARY} returned no results; attempting EC {EC_ALT} ...")

        km_alt = query_km(client, email, pw_hash, EC_ALT)
        time.sleep(1)
        section(f"Km values: ALL organisms / EC {EC_ALT}")
        display_km(km_alt)

        kcat_alt = query_kcat(client, email, pw_hash, EC_ALT)
        time.sleep(1)
        section(f"kcat values: ALL organisms / EC {EC_ALT}")
        display_kcat(kcat_alt)

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    section("SUMMARY")
    km_n    = len(km_all)
    km_l_n  = len(km_lacto)
    kcat_n  = len(kcat_all)
    eff_n   = len(eff_all)

    print(f"  Enzyme       : GlgA -- glycogen synthase (ADP-glucose donor)")
    print(f"  EC number    : {EC_PRIMARY}")
    print(f"  Target org   : {TARGET_ORG}")
    print()
    print(f"  Km entries returned (L. acidophilus only) : {km_l_n}")
    print(f"  Km entries returned (all organisms)       : {km_n}")
    print(f"  kcat entries returned (all organisms)     : {kcat_n}")
    print(f"  kcat/Km entries returned                  : {eff_n}")

    # Compute numeric Km stats
    km_nums = []
    for raw in km_all:
        p = parse_entry(raw)
        v = to_float(p.get("kmValue") or p.get("kmvalue"))
        if v is not None:
            km_nums.append(v)

    if km_nums:
        avg = sum(km_nums) / len(km_nums)
        print(f"\n  Km numeric statistics across all organisms (mM):")
        print(f"    n       : {len(km_nums)}")
        print(f"    min     : {min(km_nums):.4f}")
        print(f"    max     : {max(km_nums):.4f}")
        print(f"    mean    : {avg:.4f}")

    kcat_nums = []
    for raw in kcat_all:
        p = parse_entry(raw)
        v = to_float(p.get("turnoverNumber") or p.get("turnovernumber"))
        if v is not None:
            kcat_nums.append(v)

    if kcat_nums:
        avg_k = sum(kcat_nums) / len(kcat_nums)
        print(f"\n  kcat numeric statistics across all organisms (s-1):")
        print(f"    n       : {len(kcat_nums)}")
        print(f"    min     : {min(kcat_nums):.4f}")
        print(f"    max     : {max(kcat_nums):.4f}")
        print(f"    mean    : {avg_k:.4f}")

    print()
    print("  Done.")
    print(LINE)


if __name__ == "__main__":
    main()
