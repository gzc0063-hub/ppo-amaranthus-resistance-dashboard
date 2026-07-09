from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class CAAtom:
    chain: str
    resseq: int
    icode: str
    resname: str
    coord: np.ndarray

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)


@dataclass(frozen=True)
class MutationFocus:
    kind: str
    label: str
    position: int | None
    ref_residues: list[int]
    comp_residues: list[int]
    approximate: bool


@dataclass(frozen=True)
class AlignmentResult:
    transformed_pdb_text: str | None
    global_rmsd: float | None
    local_rmsd: float | None
    matched_residue_count: int
    local_residue_count: int
    error: str | None = None


def parse_mutation_focus(mutation: str) -> MutationFocus:
    normalized = mutation.strip().replace("_", " ")
    if normalized.upper() == "WT":
        return MutationFocus("wt", "WT", None, [], [], False)

    deletion_match = re.fullmatch(r"Delta\s*G?(\d+)", normalized, flags=re.IGNORECASE)
    if deletion_match:
        position = int(deletion_match.group(1))
        return MutationFocus(
            kind="deletion",
            label=f"Deletion site: {mutation}",
            position=position,
            ref_residues=[position],
            comp_residues=list(range(position - 8, position + 9)),
            approximate=True,
        )

    substitution_match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", normalized, flags=re.IGNORECASE)
    if substitution_match:
        position = int(substitution_match.group(2))
        return MutationFocus(
            kind="substitution",
            label=f"Mutation site: {mutation}",
            position=position,
            ref_residues=[position],
            comp_residues=[position],
            approximate=False,
        )

    return MutationFocus("unknown", f"Mutation site: {mutation}", None, [], [], False)


def parse_ca_atoms(pdb_text: str) -> dict[tuple[str, int, str], CAAtom]:
    atoms: dict[tuple[str, int, str], CAAtom] = {}
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            chain = line[21].strip() or "_"
            resseq = int(line[22:26])
            icode = line[26].strip()
            resname = line[17:20].strip()
            coord = np.array(
                [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ],
                dtype=float,
            )
        except ValueError:
            continue
        atom = CAAtom(chain=chain, resseq=resseq, icode=icode, resname=resname, coord=coord)
        atoms[atom.key] = atom
    return atoms


def common_ca_keys(
    reference_atoms: dict[tuple[str, int, str], CAAtom],
    comparison_atoms: dict[tuple[str, int, str], CAAtom],
) -> list[tuple[str, int, str]]:
    return sorted(
        set(reference_atoms).intersection(comparison_atoms),
        key=lambda key: (key[0], key[1], key[2]),
    )


def _kabsch_transform(reference_coords: np.ndarray, comparison_coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ref_centroid = reference_coords.mean(axis=0)
    comp_centroid = comparison_coords.mean(axis=0)
    ref_centered = reference_coords - ref_centroid
    comp_centered = comparison_coords - comp_centroid
    covariance = comp_centered.T @ ref_centered
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    rotation = u_matrix @ vt_matrix
    if np.linalg.det(rotation) < 0:
        u_matrix[:, -1] *= -1
        rotation = u_matrix @ vt_matrix
    translation = ref_centroid - comp_centroid @ rotation
    return rotation, translation


def _transform_point(coord: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return coord @ rotation + translation


def transform_pdb_text(pdb_text: str, rotation: np.ndarray, translation: np.ndarray) -> str:
    transformed_lines: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                coord = np.array(
                    [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ],
                    dtype=float,
                )
            except ValueError:
                transformed_lines.append(line)
                continue
            new_coord = _transform_point(coord, rotation, translation)
            line = (
                f"{line[:30]}"
                f"{new_coord[0]:8.3f}{new_coord[1]:8.3f}{new_coord[2]:8.3f}"
                f"{line[54:]}"
            )
        transformed_lines.append(line)
    return "\n".join(transformed_lines) + "\n"


def _rmsd(reference_coords: np.ndarray, comparison_coords: np.ndarray) -> float | None:
    if len(reference_coords) == 0 or len(reference_coords) != len(comparison_coords):
        return None
    diff = reference_coords - comparison_coords
    return float(math.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def residue_window(position: int | None, radius: int = 8) -> list[int]:
    if position is None:
        return []
    return list(range(position - radius, position + radius + 1))


def keys_in_residue_numbers(
    keys: Iterable[tuple[str, int, str]],
    residue_numbers: Iterable[int],
) -> list[tuple[str, int, str]]:
    residue_set = set(residue_numbers)
    return [key for key in keys if key[1] in residue_set]


def align_pdbs_and_calculate_rmsd(
    reference_pdb_text: str,
    comparison_pdb_text: str,
    mutation: str,
    local_radius: int = 8,
) -> AlignmentResult:
    reference_atoms = parse_ca_atoms(reference_pdb_text)
    comparison_atoms = parse_ca_atoms(comparison_pdb_text)
    common_keys = common_ca_keys(reference_atoms, comparison_atoms)
    if len(common_keys) < 3:
        return AlignmentResult(None, None, None, len(common_keys), 0, "fewer than 3 matched C-alpha atoms")

    reference_coords = np.array([reference_atoms[key].coord for key in common_keys])
    comparison_coords = np.array([comparison_atoms[key].coord for key in common_keys])
    rotation, translation = _kabsch_transform(reference_coords, comparison_coords)
    transformed_common = np.array([
        _transform_point(comparison_atoms[key].coord, rotation, translation) for key in common_keys
    ])
    global_rmsd = _rmsd(reference_coords, transformed_common)

    focus = parse_mutation_focus(mutation)
    local_keys: list[tuple[str, int, str]] = []
    if focus.position is not None:
        local_keys = keys_in_residue_numbers(common_keys, residue_window(focus.position, local_radius))
    if local_keys:
        local_reference_coords = np.array([reference_atoms[key].coord for key in local_keys])
        local_comparison_coords = np.array([
            _transform_point(comparison_atoms[key].coord, rotation, translation) for key in local_keys
        ])
        local_rmsd = _rmsd(local_reference_coords, local_comparison_coords)
    else:
        local_rmsd = None

    return AlignmentResult(
        transformed_pdb_text=transform_pdb_text(comparison_pdb_text, rotation, translation),
        global_rmsd=global_rmsd,
        local_rmsd=local_rmsd,
        matched_residue_count=len(common_keys),
        local_residue_count=len(local_keys),
    )
