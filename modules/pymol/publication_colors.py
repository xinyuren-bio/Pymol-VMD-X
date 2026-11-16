"""Color overrides for publication figures, without GUI dependencies."""

import math


THEME_COLORS = {
    "journal": {"protein": "slate", "ligand": "tv_orange"},
    "illustration": {"protein": "teal", "ligand": "salmon"},
    "surface": {"protein": "lightblue", "ligand": "tv_orange"},
}


def normalize_color(value, cmd):
    """Return an RGB triple or None; reject invalid input before scene mutation."""
    if value is None:
        return None
    if isinstance(value, str):
        try:
            value = cmd.get_color_tuple(value)
        except Exception as exc:
            raise ValueError("无法识别颜色，请选择有效的颜色。") from exc
    try:
        if len(value) != 3:
            raise ValueError()
        rgb = tuple(float(component) for component in value)
        if not all(math.isfinite(v) and 0 <= v <= 1 for v in rgb):
            raise ValueError()
    except (TypeError, ValueError) as exc:
        raise ValueError("颜色必须是三个 0 到 1 之间的 RGB 数值。") from exc
    return rgb


def color_literal(rgb):
    """Use true colors without modifying any shared named-color definitions."""
    return "0x" + "".join("%02x" % round(component * 255) for component in rgb)


def apply_color_overrides(cmd, protein_selection, ligand_selection,
                          protein_color=None, ligand_color=None,
                          pocket_selection=None):
    from pymol import util

    if protein_color is not None:
        cmd.color(color_literal(protein_color), protein_selection)
        if pocket_selection:
            util.cnc(pocket_selection, _self=cmd)
    if ligand_color is not None and ligand_selection:
        cmd.color(color_literal(ligand_color),
                  "(%s) and elem C" % ligand_selection)
        util.cnc(ligand_selection, _self=cmd)
