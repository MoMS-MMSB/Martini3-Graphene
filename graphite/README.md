# Graphite Slab Builder



The script automates the construction of a graphitic slab by stacking multiple graphene layers along the z direction in these steps :

1. Build periodic graphene (`martini3-graphene-periodic.py`)
2. Set the one-layer box height to `0.382 nm`
3. Stack layers along `z` to obtain graphite

## Usage

From the repository root:

```bash
python graphite/build_graphite_slab.py -x 10 -y 10 --layers 5 -o graphite
```

Useful options:

- `--spacing 0.382`: interlayer spacing in nm (default `0.382`)
- `--graphene-prefix graphene`: prefix for single-layer files
- `--posre-fc 1000`: force constant for generated `posre_GRA.itp`
- `--output-dir DIR`: write all files into `DIR`
- `--force`: overwrite existing outputs

## Generated Files

With default names, the script writes:

- `graphene.gro` (single periodic layer)
- `graphene.itp` (single-layer topology, with `POSRES` include block)
- `graphene_1layer.gro` (single layer with `z = spacing`)
- `graphite.gro` (stacked slab)
- `graphite.pdb` (stacked slab for visualization)
- `posre_GRA.itp` (position restraints for one graphene layer)
- `topol.top` (default topology with Martini + graphene includes and `GRA <layers>`)

To use the stacked slab in a system topology, include `graphene.itp` and set:

```text
[ molecules ]
GRA    <layers>
```

where `<layers>` is the value passed to `--layers` (default `5`).
