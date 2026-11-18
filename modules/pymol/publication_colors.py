"""Color overrides for publication figures, without GUI dependencies."""

import math
from colorsys import hsv_to_rgb


LIGAND_CARBON_COLORS = {
    "journal": "tv_orange",
    "illustration": "salmon",
    "surface": "tv_orange",
}
ELEMENT_CARBON_COLOR = "0x465268"

CHAIN_PALETTES = {
    "journal": ("486F94", "B16C58", "47877D", "8C71A0", "AC873F", "6D8054", "AD6384", "586476"),
    "illustration": ("139C9B", "EE7955", "7358BE", "427AD4", "E2A92D", "76A743", "D45189", "6D5760"),
    "surface": ("618EB7", "BA7F79", "589B90", "9684B1", "B39555", "829E65", "BA829C", "6B859B"),
}
_CHAIN_FALLBACK_SV = {
    "journal": (0.4, 0.63),
    "illustration": (0.72, 0.85),
    "surface": (0.42, 0.72),
}


def chain_color(theme, index):
    palette = CHAIN_PALETTES[theme]
    if index < len(palette):
        code = palette[index]
        return tuple(int(code[i:i + 2], 16) / 255 for i in (0, 2, 4))
    return hsv_to_rgb((index * 0.61803398875) % 1, *_CHAIN_FALLBACK_SV[theme])


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


def apply_colors(cmd, theme, protein_selection, chains, ligand_selection, *,
                 protein_mode, ligand_mode, protein_color, ligand_color, chain_colors):
    if chains:
        if protein_mode == "single":
            rgb = protein_color if protein_color is not None else chain_color(theme, 0)
            cmd.color(color_literal(rgb), protein_selection)
        else:
            # A chain ID can be blank or contain selection wildcards. Match the
            # actual atom fields so neither can recolor an unrelated chain.
            colors = {
                chain.key: cmd.get_color_index(color_literal(
                    chain_colors.get(chain.key, chain_color(theme, i))))
                for i, chain in enumerate(chains)
            }
            cmd.alter(protein_selection, "color = colors.get((model, segi, chain), color)",
                      space={"colors": colors})
            cmd.recolor(protein_selection)

    if ligand_selection:
        cmd.color("atomic", ligand_selection)
        carbon = ELEMENT_CARBON_COLOR
        if ligand_mode == "carbon":
            carbon = (color_literal(ligand_color) if ligand_color is not None
                      else LIGAND_CARBON_COLORS[theme])
        cmd.color(carbon, f"({ligand_selection}) and elem C")
