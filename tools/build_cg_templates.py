"""Build the curated display templates from the unmodified literature ITPs.

Run with the repository Python; no GROMACS installation or network is required.
This is a deliberately limited reader for the archived files, not a grompp clone.
Only default (FLEXIBLE undefined) bonds and constraints become display edges.
Virtual-site construction records are retained separately, never made into bonds.
"""
from pathlib import Path
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "data/cg_templates/paloncyova_2023_martini3"
SOURCE = PACKAGE / "source"


def read_molecules(path):
    molecules = {}
    current = None
    section = ""
    active = True
    stack = []
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("#"):
            parts = line.split()
            if parts[0] in ("#ifdef", "#ifndef") and parts[1:] == ["FLEXIBLE"]:
                stack.append((active, parts[0] == "#ifndef"))
                active = active and stack[-1][1]
            elif parts == ["#else"] and stack:
                parent, condition = stack[-1]
                stack[-1] = (parent, not condition)
                active = parent and not condition
            elif parts == ["#endif"] and stack:
                active = stack.pop()[0]
            else:
                raise ValueError(f"Unsupported directive at {path.name}:{lineno}: {line}")
            continue
        if not active:
            continue
        match = re.fullmatch(r"\[\s*(\w+)\s*\]", line)
        if match:
            section = match[1].lower()
            continue
        cols = line.split()
        if section == "moleculetype":
            current = {"molecule": cols[0], "source_file": path.name,
                       "atoms": [], "bonds": [], "constraints": [], "virtual_sites": []}
            if cols[0] in molecules:
                raise ValueError("Duplicate molecule " + cols[0])
            molecules[cols[0]] = current
            section = ""
        elif current is not None:
            if section == "atoms":
                if int(cols[0]) != len(current["atoms"]) + 1:
                    raise ValueError("Atom numbering must be consecutive")
                current["atoms"].append({"id": int(cols[0]), "type": cols[1],
                    "resnr": int(cols[2]), "resn": cols[3], "name": cols[4],
                    "charge": float(cols[6]),
                    "mass": float(cols[7]) if len(cols) > 7 else None})
            elif section in ("bonds", "constraints"):
                # Store zero-based molecular indices, not global object indices.
                current[section].append(sorted([int(cols[0]) - 1, int(cols[1]) - 1]))
            elif section.startswith("virtual_sites"):
                current["virtual_sites"].append({"section": section, "fields": cols})
    if stack:
        raise ValueError("Unclosed conditional in " + path.name)
    for mol in molecules.values():
        n = len(mol["atoms"])
        for pair in mol["bonds"] + mol["constraints"]:
            if not 0 <= pair[0] < pair[1] < n:
                raise ValueError("Invalid edge in " + mol["molecule"])
    return molecules


def build():
    specs = [
        ("ALC_0159_CG_final.itp", "AL15", ["AL15", "AL1"], "ALC-0159", "al15"),
        ("ALC_0315_charged_CG_final.itp", "AL31", ["AL31", "AL3"], "ALC-0315（带电源模板）", "al31"),
        ("martini_v3.0.0_phospholipids_v1.itp", "DPPC", ["DPPC", "DPP"], "辅助磷脂（源名称 DPPC）", "dppc"),
        ("martini_v3.0_sterols_beta.itp", "CHOL", ["CHOL", "CHO"], "胆固醇", "chol"),
        ("RNA.itp", "RNA", ["U"], "poly-U RNA 10-mer", "rna"),
    ]
    templates = []
    for filename, molecule, aliases, label, component in specs:
        mol = read_molecules(SOURCE / filename)[molecule]
        mol.update({"id": molecule, "label": label, "resn_aliases": aliases,
                    "component": component, "scope": "chain" if molecule == "RNA" else "residue"})
        # GROMACS editconf rotates four-letter names: EC10 -> 0EC1.
        mol["atom_aliases"] = {
            atom["name"]: [atom["name"][-1] + atom["name"][:-1]]
            for atom in mol["atoms"] if molecule == "AL15" and re.fullmatch(r"EC\d\d", atom["name"])
        }
        templates.append(mol)
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted(SOURCE.iterdir()) if p.is_file()}
    result = {
        "schema_version": 1,
        "id": "paloncyova_2023_martini3",
        "label": "Paloncýová 2023 · 定制 Martini 3.0.0",
        "forcefield": "Martini 3.0.0, author modified SN3r-water interactions",
        "citation": "Paloncýová et al., J. Phys. Chem. B 2023, 127, 1158–1166",
        "doi": "https://doi.org/10.1021/acs.jpcb.2c07671",
        "source_url": "https://zenodo.org/records/7456781",
        "archive_subdirectory": "data_CG/topologies",
        "source_license": "CC-BY-4.0",
        "retrieved_date": "2026-11-14",
        "preprocessor_defines": [],
        "connectivity_sections": ["bonds", "constraints"],
        "limitations": "Template name matching does not identify a PDB's force-field version or charge state. Virtual-site construction records are not chemical bonds. No all-atom reconstruction.",
        "source_sha256": files,
        "templates": templates,
    }
    (PACKAGE / "templates.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print("Saved", len(templates), "templates:",
          [(t["id"], len(t["atoms"]), len(t["bonds"]), len(t["constraints"])) for t in templates])


if __name__ == "__main__":
    build()
