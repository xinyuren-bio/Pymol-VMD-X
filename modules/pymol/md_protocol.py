"""Portable GROMACS stages for the exported ff14SB/TIP3P simulation bundle."""

STAGE_LABELS = {'em': '能量最小化', 'nvt': 'NVT 平衡', 'npt': 'NPT 平衡', 'md': '生产模拟'}

COMMON = """\
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
vdwtype = Cut-off
DispCorr = EnerPres
pbc = xyz
"""
EM = """\
integrator = steep
nsteps = 50000
emtol = 1000
emstep = 0.01
nstenergy = 100
""" + COMMON

DYNAMICS = """\
integrator = md
dt = 0.002
constraints = h-bonds
constraint-algorithm = LINCS
tcoupl = V-rescale
tc-grps = System
tau-t = 1.0
ref-t = {temperature}
nstlog = 1000
nstenergy = 1000
nstxout-compressed = 5000
""" + COMMON

PRESSURE = """\
pcoupl = C-rescale
pcoupltype = isotropic
tau-p = 5.0
ref-p = 1.0
compressibility = 4.5e-5
"""


def stage_inputs(temperature, production_steps):
    settings = DYNAMICS.format(temperature=temperature)
    return {
        'em': EM,
        'nvt': settings + f"""\
nsteps = 250000
define = -DPOSRES
pcoupl = no
continuation = no
gen-vel = yes
gen-temp = {temperature}
gen-seed = 20261116
""",
        'npt': settings + PRESSURE + """\
nsteps = 500000
define = -DPOSRES
refcoord-scaling = com
continuation = yes
gen-vel = no
""",
        'md': settings + PRESSURE + f"""\
nsteps = {production_steps}
continuation = yes
gen-vel = no
""",
    }


_RUN_BASE = r"""#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"

gmx_bin="${GMX_BIN:-gmx}"
case "${1:---gpu}" in
    --gpu) nonbonded=gpu ;;
    --cpu) nonbonded=cpu ;;
    *) echo "Usage: bash run_md.sh [--gpu|--cpu]" >&2; exit 2 ;;
esac
command -v "$gmx_bin" >/dev/null
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
mkdir .run-lock || { echo "Another run is active, or .run-lock needs review." >&2; exit 2; }
trap 'rmdir .run-lock' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
"$gmx_bin" --version | tee runtime-version.txt

if [[ "$nonbonded" == gpu ]] && grep -Eiq '^GPU support:[[:space:]]*(disabled|none|no)[[:space:]]*$' runtime-version.txt; then
    echo "This GROMACS executable has no GPU support. Set GMX_BIN to a GPU build." >&2
    exit 2
fi

run_stage() {
    local stage="$1" input="$2" accelerator="$3"
    shift 3
    if [[ -f "$stage.done" ]]; then
        [[ -s "$stage.gro" && -s "$stage.tpr" && -s "$stage.log" ]]
        if [[ "$stage" != em* ]]; then [[ -s "$stage.cpt" ]]; fi
        return
    fi
    echo "=== $stage ==="
    if [[ ! -s "$stage.tpr" ]]; then
        "$gmx_bin" grompp -f "mdp/$stage.mdp" -c "$input" -p "${topology:-system.top}" -n index.ndx \
            -po "${stage}-mdout.mdp" -pp "${stage}-processed.top" -o "$stage.tpr" "$@"
    fi
    local run=(mdrun -deffnm "$stage" -nb "$accelerator")
    if [[ -s "$stage.cpt" ]]; then
        run+=(-cpi "$stage.cpt" -append)
    elif [[ -e "$stage.log" || -e "$stage.edr" || -e "$stage.gro" ]]; then
        echo "$stage already has outputs and no resumable checkpoint. Review those files first." >&2
        return 1
    fi
    "$gmx_bin" "${run[@]}"
    [[ -s "$stage.gro" && -s "$stage.log" ]]
    if grep -Eiq 'LINCS WARNING|SETTLE.*(WARNING|error)|non-finite|(Potential Energy|Maximum force)[[:space:]]*=[[:space:]]*[-+]?(nan|inf)' "$stage.log"; then
        echo "$stage reported a stability warning. Review $stage.log." >&2
        return 1
    fi
    if [[ "$stage" == em* ]]; then
        if ! grep -Eq 'converged to Fmax <' "$stage.log"; then
            echo "EM did not reach the force threshold. Review the input geometry and em.log." >&2
            return 1
        fi
    else
        [[ -s "$stage.cpt" ]]
        local end_step final_step
        end_step=$(awk -F= '
            $1 ~ /^nsteps[[:space:]]*$/ {steps = $2 + 0}
            $1 ~ /^init-step[[:space:]]*$/ {initial = $2 + 0}
            END {if (steps <= 0) exit 1; printf "%.0f\n", steps + initial}
        ' "$stage-mdout.mdp")
        final_step=$("$gmx_bin" dump -cp "$stage.cpt" 2>/dev/null |
            awk '$1 == "step" && $2 == "=" && !found {print $3; found = 1} END {if (!found) exit 1}')
        if [[ "$final_step" != "$end_step" ]]; then
            echo "$stage checkpoint is at step $final_step of $end_step. Re-run to continue." >&2
            return 1
        fi
    fi
    touch "$stage.done"
}

"""

RUN_SCRIPT = _RUN_BASE + r"""
run_stage em system.gro cpu -r system.gro
run_stage nvt em.gro "$nonbonded" -r em.gro
run_stage npt nvt.gro "$nonbonded" -t nvt.cpt -r em.gro
run_stage md npt.gro "$nonbonded" -t npt.cpt -r em.gro
echo "Finished. See md.log and the production files in this directory."
"""


def charmm_gui_run_script(stages):
    rows = ['topology=topol.top', 'input=step3_input.gro', 'reference=step3_input.gro']
    previous = None
    for stage in stages:
        offload = 'cpu' if stage.startswith('em') else '"$nonbonded"'
        command = f'run_stage {stage} "$input" {offload} -r "$reference"'
        if previous is not None and not previous.startswith('em'):
            command += f' -t {previous}.cpt'
        rows.extend((command, f'input={stage}.gro'))
        previous = stage
    rows.append('echo "Finished. See md.log and the production files in this directory."')
    return _RUN_BASE + '\n'.join(rows) + '\n'
