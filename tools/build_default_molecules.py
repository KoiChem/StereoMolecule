#!/usr/bin/env python3
"""Build the offline molecule bundle used by StereoMolecule.

The app consumes a small, renderer-ready JSON representation instead of
requesting PubChem every time the default screen-saver list advances.  NCI/CIR
is used for the ordinary names; the one CID-only default uses PubChem's JSON
record endpoint because its name is deliberately not resolved by name.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "assets" / "default-molecules.json"
DEFAULT_NAMES = [
    "L-cysteine",
    "nanokid",
    "adamantane",
    "caffeine",
    "cholesterol",
    "perfluorohexanesulfonic acid",
    "morphine",
    "strychnine",
    "sucrose",
    "adenosine triphosphate",
    "porphine",
    "penicillin v",
    "hexahelicene",
    "mirex",
    "tetrodotoxin",
    "capsaicin",
    "1-methylsilatrane",
    "reserpine",
    "ergotamine",
    "doxorubicin",
    "tetracycline",
    "cryptand-222",
    "astaxanthin",
    "octanitrocubane",
    "dodecahedrane",
    "aconitine",
]
CID_ONLY = {"Fullerene-C36": "101760755"}


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


def parse_sdf(sdf: str) -> dict:
    lines = sdf.splitlines()
    if len(lines) < 4:
        raise ValueError("SDF has no counts line")
    atom_count = int(lines[3][0:3].strip())
    bond_count = int(lines[3][3:6].strip())
    atoms = []
    for line in lines[4 : 4 + atom_count]:
        atoms.append({
            "x": float(line[0:10].strip()),
            "y": float(line[10:20].strip()),
            "z": float(line[20:30].strip()),
            "element": line[31:34].strip(),
        })
    bonds = []
    for line in lines[4 + atom_count : 4 + atom_count + bond_count]:
        bonds.append({
            "a": int(line[0:3].strip()) - 1,
            "b": int(line[3:6].strip()) - 1,
            "order": int(line[6:9].strip()),
        })
    return {"atoms": atoms, "bonds": bonds, "formula": hill_formula([atom["element"] for atom in atoms])}


ELEMENTS = [
    None, "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
]


def parse_pcjson(payload: dict) -> dict:
    compound = payload["PC_Compounds"][0]
    atoms = compound["atoms"]
    coords = compound["coords"][0]
    conformer = coords["conformers"][0]
    by_aid = {}
    for index, aid in enumerate(coords["aid"]):
        by_aid[aid] = {
            "x": conformer["x"][index],
            "y": conformer["y"][index],
            "z": conformer["z"][index],
            "element": ELEMENTS[atoms["element"][index]],
        }
    ordered_aids = atoms["aid"]
    atom_index = {aid: index for index, aid in enumerate(ordered_aids)}
    result_atoms = [by_aid[aid] for aid in ordered_aids]
    bonds = [
        {"a": atom_index[a], "b": atom_index[b], "order": order}
        for a, b, order in zip(compound["bonds"]["aid1"], compound["bonds"]["aid2"], compound["bonds"]["order"])
    ]
    return {"atoms": result_atoms, "bonds": bonds, "formula": hill_formula([atom["element"] for atom in result_atoms])}


def get_text(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "StereoMolecule offline bundle builder"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return response.read().decode("utf-8")


def cir_molecule(name: str) -> dict:
    identifier = urllib.parse.quote(name, safe="")
    url = f"https://cactus.nci.nih.gov/chemical/structure/{identifier}/file?format=sdf&get3d=True"
    return parse_sdf(get_text(url))


def pubchem_molecule(namespace: str, identifier: str) -> dict:
    encoded_identifier = urllib.parse.quote(identifier, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{namespace}/{encoded_identifier}/JSON?record_type=3d"
    return parse_pcjson(json.loads(get_text(url)))


def main() -> None:
    selected_names = set(DEFAULT_NAMES + list(CID_ONLY))
    if len(sys.argv) > 1:
        selected_names = {name for name in sys.argv[1:]}
    if OUTPUT.exists():
        molecules = json.loads(OUTPUT.read_text(encoding="utf-8")).get("molecules", {})
    else:
        molecules = {}

    def save() -> None:
        OUTPUT.parent.mkdir(exist_ok=True)
        OUTPUT.write_text(json.dumps({"version": 1, "molecules": molecules}, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")

    for name in DEFAULT_NAMES:
        if name not in selected_names:
            continue
        print(f"CIR: {name}", flush=True)
        try:
            molecule = cir_molecule(name)
        except Exception as error:
            print(f"  CIR unavailable ({error}); trying PubChem JSON", flush=True)
            molecule = pubchem_molecule("name", name)
            time.sleep(1.1)
        molecules[normalized(name)] = {"name": name, **molecule}
        save()
        time.sleep(0.35)
    for name, cid in CID_ONLY.items():
        if name not in selected_names:
            continue
        print(f"PubChem CID: {name}", flush=True)
        molecules[normalized(name)] = {"name": name, "cid": cid, **pubchem_molecule("cid", cid)}
        save()
    print(f"Wrote {OUTPUT.relative_to(ROOT)} with {len(molecules)} molecules.")


if __name__ == "__main__":
    main()
