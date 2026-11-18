"""Compare two protein selections on single-state display copies."""

from dataclasses import dataclass

import numpy as np

from .analysis_selection import atom_selections
from .coordinate_fit import rigid_fit


@dataclass(frozen=True)
class Comparison:
    reference: str
    mobile: str
    rmsd: float
    matched: int
    visible_sources: tuple

    @property
    def objects(self):
        return self.reference, self.mobile


def compare(cmd, reference, mobile, reference_state, mobile_state, reference_scope, mobile_scope,
            *, syntax="pymol"):
    names = cmd.get_names("objects")
    for obj, state in ((reference, reference_state), (mobile, mobile_state)):
        if obj not in names or cmd.get_type(obj) != "object:molecule":
            raise ValueError("对照对象已不存在，请刷新列表。")
        if not 1 <= state <= cmd.count_states(obj):
            raise ValueError(f"{obj} 的指定坐标状态不可用。")

    ref_copy = cmd.get_unused_name("compare_reference")
    mob_copy = cmd.get_unused_name("compare_mobile")
    default_scope = "protein" if syntax == "vmd" else "polymer.protein"
    try:
        for output, source, state, scope, color in (
            (ref_copy, reference, reference_state, reference_scope, "0x497CA8"),
            (mob_copy, mobile, mobile_state, mobile_scope, "0xD48550"),
        ):
            expressions = (scope.strip() or default_scope,)
            with atom_selections(cmd, source, expressions, syntax=syntax, state=state) as [(selection, count)]:
                if not count:
                    raise ValueError(f"{source} 的范围没有匹配原子。")
                cmd.create(output, selection, source_state=state, target_state=1, zoom=0)
            cmd.set("state", 1, output)
            cmd.hide("everything", output)
            cmd.show("cartoon", f"{output} and polymer.protein")
            cmd.show("sticks", f"{output} and not polymer.protein")
            cmd.color(color, output)

        ref_fit, mob_fit = (f'{name} and polymer.protein and name CA and alt ""+A'
                            for name in (ref_copy, mob_copy))
        if min(cmd.count_atoms(ref_fit), cmd.count_atoms(mob_fit)) < 3:
            raise ValueError("参考和比较范围各需至少 3 个蛋白 Cα 原子。")
        mapping = cmd.get_unused_name("_compare_pairs")
        cmd.feedback("push")
        try:
            # PyMOL supplies the sequence pairing. The rigid fit below computes
            # the transform and RMSD from those pairs.
            cmd.feedback("disable", "matrix", "details")
            cmd.align(mob_fit, ref_fit, mobile_state=1, target_state=1, cycles=0,
                      transform=0, object=mapping)
            pairs = [dict(pair) for pair in cmd.get_raw_alignment(mapping)]
        finally:
            cmd.feedback("pop")
            cmd.delete(mapping)
        pairs = [pair for pair in pairs if mob_copy in pair and ref_copy in pair]
        if len(pairs) < 3:
            raise ValueError("匹配 Cα 少于 3 个，无法完成结构对齐。")
        ref_coords = dict(zip((index for _, index in cmd.index(ref_fit)), cmd.get_coords(ref_fit)))
        mob_coords = dict(zip((index for _, index in cmd.index(mob_fit)), cmd.get_coords(mob_fit)))
        rotation, translation, rmsd = rigid_fit(
            [mob_coords[pair[mob_copy]] for pair in pairs],
            [ref_coords[pair[ref_copy]] for pair in pairs],
        )
        transform = np.eye(4)
        transform[:3, :3] = rotation.T
        transform[:3, 3] = translation
        cmd.transform_object(mob_copy, transform.ravel(), state=1, homogenous=1)
    except Exception:
        for obj in (ref_copy, mob_copy):
            if obj in cmd.get_names("objects"):
                cmd.delete(obj)
        raise
    enabled = set(cmd.get_names("objects", enabled_only=1))
    visible = tuple(dict.fromkeys(obj for obj in (reference, mobile) if obj in enabled))
    return Comparison(ref_copy, mob_copy, rmsd, len(pairs), visible)


def show_overlay(cmd, result):
    for obj in result.visible_sources:
        if obj in cmd.get_names("objects"):
            cmd.disable(obj)
    for obj in result.objects:
        cmd.enable(obj)
    cmd.zoom(" or ".join(result.objects), buffer=4, state=1, animate=0)


def clear_comparison(cmd, result):
    names = set(cmd.get_names("objects"))
    for obj in result.objects:
        if obj in names:
            cmd.delete(obj)
    for obj in result.visible_sources:
        if obj in names:
            cmd.enable(obj)


def render_pair(cmd, result, left, right, *, width=600, height=600):
    """Render both panels at the same orientation, center, and scale."""
    from pymol2 import PyMOL

    if not set(result.objects) <= set(cmd.get_names("objects")):
        raise ValueError("对照副本已变化，请重新对齐。")
    session = cmd.get_session(" ".join(result.objects), partial=1, cache=0)
    with PyMOL() as engine:
        c = engine.cmd
        c.set_session(session, partial=1)
        c.set_view(cmd.get_view(), animate=0)
        c.bg_color("white")
        c.set("depth_cue", 0)
        c.set("orthoscopic", 1)
        c.set("ray_shadows", 0)
        c.zoom("all", buffer=4, state=1, animate=0)
        common_view = c.get_view()
        for show, hide, filename in ((result.reference, result.mobile, left),
                                     (result.mobile, result.reference, right)):
            c.enable(show)
            c.disable(hide)
            c.set_view(common_view, animate=0)
            c.png(str(filename), width, height, ray=1)
