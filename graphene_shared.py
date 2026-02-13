import sys
import time
import uuid
from pathlib import Path

import numpy as np


class Profiler:
    def __init__(self, enabled):
        self.enabled = enabled
        self.last = time.perf_counter()

    def mark(self, label):
        if self.enabled:
            now = time.perf_counter()
            print(f"[profile] {label}: {now - self.last:.3f}s", file=sys.stderr)
            self.last = now


def make_temp_path(target_path):
    return target_path.parent / f".{target_path.stem}.{uuid.uuid4().hex}{target_path.suffix}"


def resolve_output_paths(output, output_dir):
    if output_dir is None:
        output_prefix = Path(output)
    else:
        output_prefix = Path(output_dir) / Path(output).name

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    gro_path = output_prefix.with_suffix(".gro")
    itp_path = output_prefix.with_suffix(".itp")
    return gro_path, itp_path


def warn_gro_atom_index_wrap(n_atoms, quiet):
    if n_atoms > 99999 and not quiet:
        print(
            "Warning: GRO atom index field is limited to 5 digits; indices wrap above 99999.",
            file=sys.stderr,
        )


def build_spatial_index(xy_positions, cell_size=2.60):
    cell_coordinates = np.floor(xy_positions / cell_size).astype(np.int32)
    spatial_index = {}
    for atom_index, (cx, cy) in enumerate(cell_coordinates):
        key = (int(cx), int(cy))
        spatial_index.setdefault(key, []).append(atom_index)
    return spatial_index, cell_coordinates


def collect_hexagons_from_virtual_sites(
    universe,
    virtual_site_indices,
    neighbor_distance_min=2.55,
    neighbor_distance_max=2.57,
    cell_size=2.60,
):
    xy_positions = universe.atoms.positions[:, :2]
    spatial_index, cell_coordinates = build_spatial_index(xy_positions, cell_size=cell_size)
    min_distance_sq = neighbor_distance_min ** 2
    max_distance_sq = neighbor_distance_max ** 2

    hexagon_indices = []
    for virtual_site_index in virtual_site_indices:
        center = xy_positions[virtual_site_index]
        cell_x, cell_y = cell_coordinates[virtual_site_index]
        neighbors = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for candidate_index in spatial_index.get((int(cell_x + dx), int(cell_y + dy)), []):
                    if candidate_index == virtual_site_index:
                        continue
                    delta = xy_positions[candidate_index] - center
                    distance_sq = delta[0] * delta[0] + delta[1] * delta[1]
                    if min_distance_sq <= distance_sq <= max_distance_sq:
                        neighbors.append(candidate_index)

        unique_neighbors = sorted(set(neighbors))
        if len(unique_neighbors) != 6:
            raise RuntimeError(
                f"Expected 6 neighbors for virtual site {virtual_site_index}, found {len(unique_neighbors)}."
            )
        unique_neighbors.append(virtual_site_index)
        hexagon_indices.append(unique_neighbors)
    return hexagon_indices


def unique_bonds_from_template(hexagons, template):
    unique = set()
    result = []
    for hexa in hexagons:
        for i, j in template:
            a, b = hexa[i], hexa[j]
            if a > b:
                a, b = b, a
            key = (a, b)
            if key not in unique:
                unique.add(key)
                result.append([a, b])
    return result


def build_row_groups(universe, y_coordinates, atol=1e-3):
    y_values = universe.atoms.positions[:, 1]
    groups = [universe.atoms[np.isclose(y_values, y, atol=atol, rtol=0)] for y in y_coordinates]
    if any(group.n_atoms == 0 for group in groups):
        raise RuntimeError("Failed to map one or more atoms to y-row groups.")
    return groups
