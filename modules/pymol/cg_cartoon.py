"""Smooth display geometry for bonded coarse-grained nucleic acid backbones."""

from collections import defaultdict
from math import dist

from pymol import cgo
from pymol.coarse_grained import NA_BB_CAP


BACKBONE_NAMES = frozenset(('BB1', 'BB2', 'BB3'))
TUBE_RADIUS = 0.42
BASE_RADIUS = 0.18


def _residue(atom):
    return atom.segi, atom.chain, atom.resi, atom.resn


def _index_pair(i, j):
    return (i, j) if i < j else (j, i)


def _paths(graph):
    """Trace maximal runs, splitting at termini, breaks, and junctions."""
    traversed = set()
    # Open paths first, then any cycles. Each edge is visited once.
    for start in sorted(graph, key=lambda i: (len(graph[i]) == 2, i)):
        if not graph[start]:
            yield [start]
        for following in sorted(graph[start]):
            path = [start]
            current = start
            while _index_pair(current, following) not in traversed:
                traversed.add(_index_pair(current, following))
                path.append(following)
                if len(graph[following]) != 2:
                    break
                current, following = following, next(i for i in graph[following] if i != current)
            if len(path) > 1:
                yield path


def _draw_tube(result, points, colors):
    """Interpolate through the bead centers with local bounded tangents."""
    from numpy import array

    points = array(points)
    colors = array(colors)
    tangents = [points[1] - points[0]]
    for prev, point, following in zip(points, points[1:], points[2:]):
        tangent = (following - prev) * 0.5
        limit = min(dist(prev, point), dist(point, following))
        length = dist((0, 0, 0), tangent)
        tangents.append(tangent * (limit / length) if length > limit else tangent)
    tangents.append(points[-1] - points[-2])
    if (points[0] == points[-1]).all():
        # Give both ends of a closed run the same tangent at the seam.
        seam = 0.5 * (points[1] - points[-2])
        tangents[0] = tangents[-1] = seam

    for i, (a, b) in enumerate(zip(points, points[1:])):
        left, right = tangents[i:i + 2]
        rgb0, rgb1 = colors[i:i + 2]
        last_point, last_rgb = a, rgb0
        for step in range(1, 7):
            t = step / 6
            point = ((2*t**3 - 3*t**2 + 1)*a + (t**3 - 2*t**2 + t)*left
                     + (-2*t**3 + 3*t**2)*b + (t**3 - t**2)*right)
            rgb = (1 - t)*rgb0 + t*rgb1
            result.extend([cgo.SAUSAGE, *last_point, *point, TUBE_RADIUS, *last_rgb, *rgb])
            last_point, last_rgb = point, rgb


def geometry(model, cg_indices, colors):
    """Create one frame from the selected atoms and their existing bond graph.

    The distance cap only discards discontinuities from periodic wrapping.
    It does not attach neighboring residues or components.
    """
    atoms = model.atom
    names = [atom.name.strip().upper() for atom in atoms]
    selected = {i for i, atom in enumerate(atoms) if atom.index in cg_indices}
    backbone = {i: set() for i in selected if names[i] in BACKBONE_NAMES}
    base = {i for i in selected if names[i].startswith('SC')}
    base_graph = defaultdict(set)
    attachments = defaultdict(set)
    has_backbone_bond = False
    for bond in model.bond:
        i, j = bond.index
        if i not in selected or j not in selected:
            continue
        a, b = atoms[i], atoms[j]
        if (a.segi, a.chain) != (b.segi, b.chain) or not 0.01 < dist(a.coord, b.coord) <= NA_BB_CAP:
            continue
        same_residue = _residue(a) == _residue(b)
        if i in backbone and j in backbone:
            if same_residue or {names[i], names[j]} == {'BB1', 'BB3'}:
                backbone[i].add(j)
                backbone[j].add(i)
                has_backbone_bond = True
        elif same_residue:
            if i in base and j in base:
                base_graph[i].add(j)
                base_graph[j].add(i)
            elif i in base and j in backbone:
                attachments[j].add(i)
            elif j in base and i in backbone:
                attachments[i].add(j)

    def rgb(i):
        return colors[atoms[i].index]

    result = []
    for path in _paths(backbone):
        if len(path) == 1:
            i = path[0]
            result.extend([cgo.COLOR, *rgb(i), cgo.SPHERE, *atoms[i].coord, TUBE_RADIUS])
        else:
            _draw_tube(result, [atoms[i].coord for i in path], [rgb(i) for i in path])

    for anchor, seeds in sorted(attachments.items()):
        component = set(seeds)
        pending = list(seeds)
        while pending:
            neighbors = base_graph[pending.pop()] - component
            component.update(neighbors)
            pending.extend(neighbors)
        if len(component) < 2:
            continue
        origin = atoms[anchor].coord
        central = [sum(values) / len(component) for values in
                   zip(*(atoms[i].coord for i in sorted(component)))]
        if 0.4 < dist(origin, central) <= NA_BB_CAP:
            result.extend([cgo.SAUSAGE, *origin, *central, BASE_RADIUS,
                           *rgb(anchor), *rgb(min(seeds))])
    return result, has_backbone_bond


def display_states(cmd, selection, cg_indices):
    """Use the representation copy's coordinates and colors for each state."""
    color_indices = {}
    cmd.iterate(selection, 'indices[index] = color', space={'indices': color_indices})
    palette = {index: cmd.get_color_tuple(index) for index in set(color_indices.values()) if index >= 0}
    colors = {i: palette.get(color) or (0.7, 0.7, 0.7) for i, color in color_indices.items()}
    states = []
    has_backbone_bond = False
    for state in range(1, cmd.count_states(selection) + 1):
        data, connected = geometry(cmd.get_model(selection, state), cg_indices, colors)
        has_backbone_bond |= connected
        states.append(data)
    return states, has_backbone_bond
