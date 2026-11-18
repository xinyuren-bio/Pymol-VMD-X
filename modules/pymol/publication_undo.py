"""One-step undo of figure display changes without reloading the session."""

from uuid import uuid4

from .publication_recipe import read_settings


BOND_SETTINGS = ("stick_radius", "stick_color", "stick_transparency")
ATOM_KEY = "(model, rank, segi, chain, resi, name, alt)"
ATOM_REPS = ("spheres", "sticks", "surface", "labels", "nb_spheres", "cartoon", "ribbon",
             "lines", "dots", "mesh", "nonbonded", "ellipsoid")


class FigureUndo:
    def __init__(self, cmd, objects, scope, global_names, local_names):
        self.cmd = cmd
        self.selection = "_pubfig_undo_" + uuid4().hex
        self.global_values = read_settings(cmd, global_names)
        self.view = cmd.get_view()
        indices = {cmd.setting._get_index(name): name for name in local_names}
        self.atoms = {}
        self.bonds = []
        self.objects = {}
        self.reps = {}
        for obj in objects:
            self.objects[obj] = [
                (state, indices[i], value)
                for state in range(cmd.count_states(obj) + 1)
                for i, _, value in (cmd.get_object_settings(obj, state) or [])
                if i in indices
            ]

        atom_keys = {}

        def collect(model, index, key, color, settings):
            atom_keys[model, index] = key
            self.atoms[key] = color, {i: settings[i] for i in settings if i in indices}

        cmd.iterate(scope, f"collect(model, index, {ATOM_KEY}, color, s)", space={"collect": collect})
        for setting in BOND_SETTINGS:
            self.bonds.extend(
                (setting, atom_keys[model, a], atom_keys[model, b], value)
                for model, bonds in cmd.get_bond(setting, scope)
                for a, b, value in bonds if value is not None
            )
        self.local_names = local_names
        self.indices = indices
        try:
            cmd.select(self.selection, scope, enable=0, quiet=1)
            # Selections track the same atoms even when their indices shift.
            for rep in ATOM_REPS:
                name = f"{self.selection}_{rep}"
                if cmd.select(name, f"{self.selection} and rep {rep}", enable=0, quiet=1):
                    self.reps[rep] = name
                else:
                    cmd.delete(name)
        except Exception:
            self.discard()
            raise

    def discard(self):
        names = set(self.cmd.get_names("selections"))
        for name in (self.selection, *self.reps.values()):
            if name in names:
                self.cmd.delete(name)

    def restore(self):
        cmd = self.cmd
        if self.selection not in cmd.get_names("selections"):
            self.discard()
            return False
        try:
            remaining = set(cmd.get_object_list(self.selection))
            for obj in remaining & self.objects.keys():
                for name in self.local_names:
                    cmd.unset(name, obj, updates=0)
                for state, name, value in self.objects[obj]:
                    if state <= cmd.count_states(obj):
                        cmd.set(name, value, obj, state=state, updates=0)

            if remaining:
                atom_indices = {}

                def restore_atom(index, key, color, settings):
                    atom_indices[key] = index
                    old_color, old = self.atoms.get(key, (color, {}))
                    for i in (set(settings) & self.indices.keys()) - old.keys():
                        del settings[i]
                    for i, value in old.items():
                        settings[i] = value
                    return old_color

                cmd.alter(self.selection, f"color = restore_atom(index, {ATOM_KEY}, color, s)",
                          space={"restore_atom": restore_atom})
                for setting in BOND_SETTINGS:
                    cmd.unset_bond(setting, self.selection, updates=0)
                for setting, a, b, value in self.bonds:
                    if a in atom_indices and b in atom_indices:
                        cmd.set_bond(setting, value,
                                     f"({self.selection} and model {a[0]} and index {atom_indices[a]})",
                                     f"({self.selection} and model {b[0]} and index {atom_indices[b]})",
                                     updates=0)
                cmd.hide("everything", self.selection)
                available = set(cmd.get_names("selections"))
                for rep, selection in self.reps.items():
                    if selection in available:
                        cmd.show(rep, selection)
                cmd.rebuild(self.selection)
            for name, value in self.global_values.items():
                cmd.set(name, value)
            cmd.set_view(self.view, animate=0)
            return True
        finally:
            self.discard()
