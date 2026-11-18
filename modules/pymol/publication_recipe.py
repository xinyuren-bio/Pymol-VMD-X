"""Portable figure recipes. JSON contains settings and identifiers, not coordinates."""

import json
import math
from pathlib import Path

from . import publication as pub
from .publication_colors import color_literal, normalize_color


APPEARANCE_SETTINGS = (
    "bg_rgb", "bg_rgb_top", "bg_rgb_bottom", "bg_gradient", "ambient", "direct",
    "reflect", "light_count", "shininess", "specular", "specular_intensity",
    "ambient_occlusion_mode", "ray_trace_mode", "ray_trace_color", "ray_trace_gain",
    "ray_shadows", "depth_cue", "fog_start", "ray_trace_fog", "ray_trace_fog_start",
    "ray_orthoscopic",
) + ("light",) + tuple(f"light{i}" for i in range(2, 10))


def read_settings(cmd, names):
    settings = {}
    for name in names:
        kind, values = cmd.get_setting_tuple(name)
        value = values if kind == 4 else values[0]
        if kind == 5 and value >= 0:
            value = color_literal(cmd.get_color_tuple(value))
        settings[name] = value
    return settings


def ligand_id(ligand):
    return ligand.object, ligand.segi, ligand.chain, ligand.resi, ligand.resn


def save_recipe(path, *, view, appearance, theme, style, colors, chains, ligands, width, height,
                dpi, transparent):
    overrides = dict(colors)
    chain_colors = overrides.pop("chain_colors")
    recipe = {
        "format": "pymol-vmd-x-figure",
        "version": 1,
        "view": list(view),
        "appearance": appearance,
        "theme": theme,
        "style": style,
        "colors": overrides,
        "chains": [list(chain.key) for chain in chains],
        "chain_colors": [{"chain": list(key), "rgb": rgb} for key, rgb in sorted(chain_colors.items())],
        "ligands": [list(ligand_id(ligand)) for ligand in ligands],
        "png": {"width": width, "height": height, "dpi": dpi, "transparent": transparent},
    }
    Path(path).write_text(json.dumps(recipe, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8")


def _vector(value, length):
    if not isinstance(value, list) or len(value) != length:
        raise ValueError("图稿中的数值向量长度不正确。")
    if any(type(x) not in (float, int) or abs(x) > 1e9 or not math.isfinite(x) for x in value):
        raise ValueError("图稿中的数值超出范围。")
    return tuple(value)


def _key(value, length):
    if not isinstance(value, list) or len(value) != length or not all(isinstance(x, str) for x in value):
        raise ValueError("图稿中的链或配体标识不完整。")
    return tuple(value)


def _appearance(values, cmd):
    if values is None:
        return None
    if not isinstance(values, dict) or not values.keys() <= set(APPEARANCE_SETTINGS):
        raise ValueError("图稿包含不支持的光照设置。")
    for name, value in values.items():
        kind = cmd.get_setting_tuple(name)[0]
        if kind == 4:
            _vector(value, 3)
        elif kind == 5 and isinstance(value, str):
            if len(value) != 8 or not value.startswith("0x"):
                raise ValueError("图稿的颜色必须为 RGB 十六进制值。")
            int(value[2:], 16)
        elif kind == 3:
            _vector([value], 1)
        elif type(value) is not int or not -10 <= value <= 100:
            raise ValueError(f"图稿的设置值无效：{name}")
    return values


def _match(key, available):
    if key in available:
        return key
    candidates = [other for other in available if other[1:] == key[1:]]
    return candidates[0] if len(candidates) == 1 else None


def load_recipe(path, cmd, scene):
    try:
        recipe = json.loads(Path(path).read_text(encoding="utf-8"))
        if recipe["format"] != "pymol-vmd-x-figure" or type(recipe["version"]) is not int or recipe["version"] != 1:
            raise ValueError("无法识别图稿格式或版本。")
        theme = pub.validate_theme(recipe["theme"])
        style = recipe["style"]
        if style not in (pub.STYLE_OVERVIEW, pub.STYLE_POCKET):
            raise ValueError("无法识别图稿的布局。")
        view = _vector(recipe["view"], 18)
        appearance = _appearance(recipe["appearance"], cmd)
        colors = recipe["colors"]
        if colors["protein_mode"] not in ("chain", "single") or colors["ligand_mode"] not in ("element", "carbon"):
            raise ValueError("无法识别图稿的着色方式。")
        colors = {
            "protein_mode": colors["protein_mode"],
            "ligand_mode": colors["ligand_mode"],
            "protein_color": normalize_color(_vector(colors["protein_color"], 3), cmd)
            if colors["protein_color"] is not None else None,
            "ligand_color": normalize_color(_vector(colors["ligand_color"], 3), cmd)
            if colors["ligand_color"] is not None else None,
            "chain_colors": {},
        }
        png = recipe["png"]
        for key, low, high in (("width", 100, 8000), ("height", 100, 8000), ("dpi", 72, 600)):
            if type(png[key]) is not int or not low <= png[key] <= high:
                raise ValueError("图稿的 PNG 尺寸或分辨率超出范围。")
        if type(png["transparent"]) is not bool:
            raise ValueError("图稿的透明背景设置无效。")

        unmatched = []
        chains = {chain.key for chain in scene.protein_chains}
        missing_chains = {_key(entry, 3) for entry in recipe["chains"]}
        missing_chains = {key for key in missing_chains if _match(key, chains) is None}
        for entry in recipe["chain_colors"]:
            key = _key(entry["chain"], 3)
            color = normalize_color(_vector(entry["rgb"], 3), cmd)
            match = _match(key, chains)
            if match is None or match in colors["chain_colors"]:
                missing_chains.add(key)
            else:
                colors["chain_colors"][match] = color
        unmatched.extend("链 " + "/".join(key) for key in sorted(missing_chains))
        ligands = {ligand_id(ligand): ligand for ligand in scene.ligands}
        chosen = {}
        for entry in recipe["ligands"]:
            key = _key(entry, 5)
            match = _match(key, ligands)
            if match is None:
                unmatched.append("配体 " + "/".join(key))
            else:
                chosen[match] = ligands[match]
    except (TypeError, KeyError, AttributeError) as exc:
        raise ValueError("图稿文件缺少必需字段，或字段类型不正确。") from exc

    return {
        "theme": theme, "style": style, "view": view, "appearance": appearance,
        "colors": colors, "png": png, "ligands": list(chosen.values()), "unmatched": unmatched,
    }
