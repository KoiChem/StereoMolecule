#!/usr/bin/env python3
"""Build and verify the offline 3D molecule bundle used by StereoMolecule.

The application reads this bundle before attempting any network request, so
every record must be a verified PubChem 3D conformer.  A build failure leaves
the previously published bundle untouched rather than silently publishing a
2D fallback.
"""

from __future__ import annotations

import json
import math
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "assets" / "default-molecules.json"
BUNDLE_VERSION = 2
PUBCHEM_REQUEST_INTERVAL_SECONDS = 1.1

# Pin each default to the CID selected for the screen saver.  This prevents a
# synonym lookup from changing the displayed compound when the bundle is rebuilt.
# The third value allows intentionally planar structures such as porphine.
DEFAULT_COMPOUNDS = (
    ("L-cysteine", "5862", False),
    ("nanokid", "11353257", False),
    ("adamantane", "9238", False),
    ("caffeine", "2519", False),
    ("cholesterol", "5997", False),
    ("perfluorohexanesulfonic acid", "67734", False),
    ("morphine", "5288826", False),
    ("strychnine", "441071", False),
    ("sucrose", "5988", False),
    ("adenosine triphosphate", "5957", False),
    ("porphine", "66868", True),
    ("penicillin v", "6869", False),
    ("hexahelicene", "98863", False),
    ("mirex", "16945", False),
    ("tetrodotoxin", "11174599", False),
    ("capsaicin", "1548943", False),
    ("1-methylsilatrane", "16797", False),
    ("reserpine", "5770", False),
    ("ergotamine", "8223", False),
    ("doxorubicin", "31703", False),
    ("tetracycline", "54675776", False),
    ("cryptand-222", "72801", False),
    ("astaxanthin", "5281224", False),
    ("octanitrocubane", "11762357", False),
    ("dodecahedrane", "123218", False),
    ("aconitine", "245005", False),
    ("Fullerene-C36", "101760755", False),
)

ELEMENTS = [
    None, "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
]


def normalized(value: str) -> str:
    return " ".join(value.strip().casefold().split())


def hill_formula(elements: list[str]) -> str:
    counts = Counter(elements)
    ordered = []
    if "C" in counts:
        ordered.append("C")
    if "H" in counts:
        ordered.append("H")
    ordered.extend(sorted(element for element in counts if element not in {"C", "H"}))
    return "".join(element + (str(counts[element]) if counts[element] != 1 else "") for element in ordered)


def get_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "StereoMolecule verified 3D bundle builder"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8")


def get_pubchem_3d(cid: str) -> dict:
    encoded_cid = urllib.parse.quote(cid, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{encoded_cid}/JSON?record_type=3d"
    return json.loads(get_text(url))


def parse_pubchem_3d(payload: dict, expected_cid: str) -> dict:
    compound = payload.get("PC_Compounds", [None])[0]
    if not compound:
        raise ValueError("PubChem record is missing")
    received_cid = str(compound.get("id", {}).get("id", {}).get("cid", ""))
    if received_cid != expected_cid:
        raise ValueError(f"CID mismatch: expected {expected_cid}, received {received_cid or 'none'}")

    atom_ids = compound.get("atoms", {}).get("aid")
    atomic_numbers = compound.get("atoms", {}).get("element")
    coordinate_set = next(
        (
            coords for coords in compound.get("coords", [])
            if coords.get("aid") and coords.get("conformers")
            and all(key in coords["conformers"][0] for key in ("x", "y", "z"))
        ),
        None,
    )
    if not isinstance(atom_ids, list) or not isinstance(atomic_numbers, list) or not coordinate_set:
        raise ValueError("PubChem 3D coordinate set is missing")

    conformer = coordinate_set["conformers"][0]
    coordinate_ids = coordinate_set["aid"]
    if not (len(coordinate_ids) == len(conformer["x"]) == len(conformer["y"]) == len(conformer["z"])):
        raise ValueError("PubChem 3D coordinate array lengths do not match")

    coordinates = {
        aid: {"x": conformer["x"][index], "y": conformer["y"][index], "z": conformer["z"][index]}
        for index, aid in enumerate(coordinate_ids)
    }
    atom_index = {aid: index for index, aid in enumerate(atom_ids)}
    atoms = []
    for index, aid in enumerate(atom_ids):
        coordinate = coordinates.get(aid)
        element = ELEMENTS[atomic_numbers[index]] if atomic_numbers[index] < len(ELEMENTS) else None
        if not coordinate or not element:
            raise ValueError("PubChem 3D coordinates do not cover every atom")
        atoms.append({**coordinate, "element": element})

    bonds = [
        {"a": atom_index[aid1], "b": atom_index[aid2], "order": order}
        for aid1, aid2, order in zip(
            compound.get("bonds", {}).get("aid1", []),
            compound.get("bonds", {}).get("aid2", []),
            compound.get("bonds", {}).get("order", []),
        )
        if aid1 in atom_index and aid2 in atom_index
    ]
    return {"atoms": atoms, "bonds": bonds, "formula": hill_formula([atom["element"] for atom in atoms])}


def spatial_shape_ratio(atoms: list[dict]) -> float:
    """Return a scale-normalized 3D covariance determinant.

    A strictly planar coordinate set has a zero determinant regardless of how
    its plane is rotated.  The ratio is intentionally only used to reject a
    numerical zero; it is not a chemistry or conformer-energy assessment.
    """
    count = len(atoms)
    if count < 4:
        return 0.0
    center = {axis: sum(atom[axis] for atom in atoms) / count for axis in ("x", "y", "z")}
    covariance = {
        (left, right): sum((atom[left] - center[left]) * (atom[right] - center[right]) for atom in atoms) / count
        for left in ("x", "y", "z") for right in ("x", "y", "z")
    }
    determinant = (
        covariance[("x", "x")] * (covariance[("y", "y")] * covariance[("z", "z")] - covariance[("y", "z")] * covariance[("z", "y")])
        - covariance[("x", "y")] * (covariance[("y", "x")] * covariance[("z", "z")] - covariance[("y", "z")] * covariance[("z", "x")])
        + covariance[("x", "z")] * (covariance[("y", "x")] * covariance[("z", "y")] - covariance[("y", "y")] * covariance[("z", "x")])
    )
    scale = (covariance[("x", "x")] + covariance[("y", "y")] + covariance[("z", "z")]) / 3
    return abs(determinant) / (scale ** 3) if scale > 0 else 0.0


def validate_molecule(name: str, molecule: dict, allow_planar: bool) -> None:
    atoms = molecule.get("atoms", [])
    bonds = molecule.get("bonds", [])
    if not atoms or not bonds or not molecule.get("formula"):
        raise ValueError(f"{name}: atoms, bonds, or formula are missing")
    if any(
        not atom.get("element") or not all(math.isfinite(atom.get(axis, math.nan)) for axis in ("x", "y", "z"))
        for atom in atoms
    ):
        raise ValueError(f"{name}: coordinates are incomplete")
    if any(not (0 <= bond["a"] < len(atoms) and 0 <= bond["b"] < len(atoms) and bond["a"] != bond["b"]) for bond in bonds):
        raise ValueError(f"{name}: bond indices are invalid")
    if not allow_planar and spatial_shape_ratio(atoms) <= 1e-12:
        raise ValueError(f"{name}: rejected a planar coordinate set for a non-planar candidate")


def atomic_write_bundle(molecules: dict) -> None:
    OUTPUT.parent.mkdir(exist_ok=True)
    temporary = OUTPUT.with_name(f".{OUTPUT.name}.tmp")
    temporary.write_text(
        json.dumps({"version": BUNDLE_VERSION, "molecules": molecules}, ensure_ascii=False, separators=(",", ":")),
        encoding="utf-8",
    )
    temporary.replace(OUTPUT)


def verify_bundle() -> dict:
    bundle = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if bundle.get("version") != BUNDLE_VERSION:
        raise ValueError(f"Expected bundle version {BUNDLE_VERSION}, received {bundle.get('version')}")
    molecules = bundle.get("molecules", {})
    expected = {normalized(name): (cid, allow_planar) for name, cid, allow_planar in DEFAULT_COMPOUNDS}
    if set(molecules) != set(expected):
        missing = ", ".join(sorted(set(expected) - set(molecules)))
        extra = ", ".join(sorted(set(molecules) - set(expected)))
        raise ValueError(f"Bundle must contain exactly the default candidates (missing: {missing or 'none'}; extra: {extra or 'none'})")
    for key, (cid, allow_planar) in expected.items():
        molecule = molecules[key]
        if molecule.get("cid") != cid:
            raise ValueError(f"{molecule.get('name', key)}: expected CID {cid}, received {molecule.get('cid')}")
        validate_molecule(molecule.get("name", key), molecule, allow_planar)
    return molecules


def main() -> None:
    if sys.argv[1:] == ["--verify"]:
        molecules = verify_bundle()
        print(f"Verified {OUTPUT.relative_to(ROOT)} with {len(molecules)} molecules.")
        return

    requested = {normalized(name) for name in sys.argv[1:]}
    known = {normalized(name) for name, _, _ in DEFAULT_COMPOUNDS}
    if requested and not requested <= known:
        unknown = ", ".join(sorted(requested - known))
        raise SystemExit(f"Unknown default compound: {unknown}")

    if requested and OUTPUT.exists():
        molecules = json.loads(OUTPUT.read_text(encoding="utf-8")).get("molecules", {})
    else:
        molecules = {}

    selected = [entry for entry in DEFAULT_COMPOUNDS if not requested or normalized(entry[0]) in requested]
    for index, (name, cid, allow_planar) in enumerate(selected):
        print(f"PubChem 3D: {name} (CID {cid})", flush=True)
        molecule = parse_pubchem_3d(get_pubchem_3d(cid), cid)
        validate_molecule(name, molecule, allow_planar)
        molecules[normalized(name)] = {"name": name, "cid": cid, **molecule}
        if index < len(selected) - 1:
            time.sleep(PUBCHEM_REQUEST_INTERVAL_SECONDS)

    atomic_write_bundle(molecules)
    verified = verify_bundle()
    print(f"Wrote and verified {OUTPUT.relative_to(ROOT)} with {len(verified)} molecules.")


if __name__ == "__main__":
    main()
