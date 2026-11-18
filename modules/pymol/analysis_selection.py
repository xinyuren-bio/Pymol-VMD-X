"""Resolve user selections to fixed atom sets for analysis."""

from contextlib import contextmanager

from .publication import _quote_sel


@contextmanager
def atom_selections(cmd, obj, expressions, *, syntax="pymol", state=1):
    """Yield (selection name, atom count) pairs and clean up temporary names."""
    if syntax not in ("pymol", "vmd"):
        raise ValueError("请选择 VMD 或 PyMOL 选择语法。")
    if syntax == "vmd":
        from .visual_representations import evaluate_selection, load_atom_records

        records = load_atom_records(obj, cmd, state=state)
        indices = [evaluate_selection(expression, records)[0] for expression in expressions]
    else:
        indices = None

    names = []
    selected = []
    try:
        for i, expression in enumerate(expressions):
            name = cmd.get_unused_name("_analysis_atoms")
            names.append(name)
            if indices is None:
                count = cmd.select(name, f"({_quote_sel(obj)}) and ({expression})", enable=0)
            elif indices[i]:
                count = cmd.select_list(name, obj, indices[i], mode="index")
                cmd.disable(name)
            else:
                count = cmd.select(name, "none", enable=0)
            selected.append((name, count))
        yield selected
    finally:
        if names:
            cmd.delete(" ".join(names))
