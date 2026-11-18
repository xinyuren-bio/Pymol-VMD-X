"""Distances and least-squares RMSD over the loaded coordinate states."""

import csv
from dataclasses import dataclass

import numpy as np

from .analysis_selection import atom_selections
from .coordinate_fit import rigid_fit


@dataclass(frozen=True)
class TrajectoryTrace:
    object: str
    mode: str
    frames: tuple
    values: tuple
    reference: int
    selection_a: str
    selection_b: str
    selection_syntax: str = "pymol"

    def write_csv(self, path):
        with open(path, "w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(("state", f"{self.mode}_angstrom", "object", "reference_state",
                             "selection_a", "selection_b", "selection_syntax"))
            for frame, value in zip(self.frames, self.values):
                writer.writerow((frame, f"{value:.9g}", self.object, self.reference or "",
                                 self.selection_a, self.selection_b, self.selection_syntax))


def analyze(cmd, obj, mode, start, end, step, selection_a, selection_b="", reference=1,
            *, syntax="pymol"):
    if obj not in cmd.get_names("objects") or cmd.get_type(obj) != "object:molecule":
        raise ValueError("轨迹对象已不存在，请刷新。")
    states = cmd.count_states(obj)
    if not 1 <= start <= end <= states or step < 1:
        raise ValueError("帧范围或步长无效。")
    if mode not in {"distance", "rmsd"}:
        raise ValueError("无法识别轨迹分析方式。")
    if not selection_a.strip() or (mode == "distance" and not selection_b.strip()):
        raise ValueError("请输入要测量的原子选择。")
    if mode == "rmsd" and not 1 <= reference <= states:
        raise ValueError("RMSD 参考帧不存在。")

    frames = tuple(range(start, end + 1, step))
    expressions = (selection_a, selection_b) if mode == "distance" else (selection_a,)
    selection_state = start if mode == "distance" else reference

    def coords(selection, frame):
        name, count = selection
        values = cmd.get_coords(name, state=frame)
        if values is None or len(values) != count or not np.isfinite(values).all():
            raise ValueError(f"第 {frame} 帧缺少固定选区的坐标，请核对原子是否贯穿所有帧。")
        return np.asarray(values, dtype=np.float64)

    with atom_selections(cmd, obj, expressions, syntax=syntax, state=selection_state) as selections:
        a = selections[0]
        if mode == "distance":
            b = selections[1]
            if a[1] != 1 or b[1] != 1:
                raise ValueError(f"距离需要两个单原子选择，当前分别匹配 {a[1]} 和 {b[1]} 个。")
            if cmd.count_atoms(f"{a[0]} and {b[0]}"):
                raise ValueError("请选择两个不同的原子。")
            values = tuple(float(np.linalg.norm(coords(a, frame) - coords(b, frame)))
                           for frame in frames)
        else:
            if a[1] < 3:
                raise ValueError("RMSD 需选择至少 3 个原子。")
            target = coords(a, reference)
            values = tuple(rigid_fit(coords(a, frame), target)[2] for frame in frames)
    return TrajectoryTrace(obj, mode, frames, values, reference if mode == "rmsd" else 0,
                           selection_a, selection_b if mode == "distance" else "", syntax)
