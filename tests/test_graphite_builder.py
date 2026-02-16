import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "graphite" / "build_graphite_slab.py"
TIMEOUT_SECONDS = 180


def run_builder(*extra_args):
    cmd = [sys.executable, str(SCRIPT_PATH), *extra_args]
    return subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=TIMEOUT_SECONDS,
    )


def read_gro_atom_count(gro_path):
    lines = gro_path.read_text().splitlines()
    return int(lines[1].strip())


class GraphiteBuilderTests(unittest.TestCase):
    def test_generates_graphitic_surface_outputs(self):
        layers = 3
        spacing = 0.382

        with tempfile.TemporaryDirectory() as tmp_dir:
            result = run_builder(
                "-x",
                "3",
                "-y",
                "3",
                "--layers",
                str(layers),
                "--spacing",
                str(spacing),
                "-o",
                "graphite_test",
                "--output-dir",
                tmp_dir,
                "--quiet",
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)

            graphene_gro = Path(tmp_dir) / "graphene.gro"
            graphene_itp = Path(tmp_dir) / "graphene.itp"
            graphene_1layer = Path(tmp_dir) / "graphene_1layer.gro"
            graphite_gro = Path(tmp_dir) / "graphite_test.gro"
            graphite_pdb = Path(tmp_dir) / "graphite_test.pdb"
            posre_itp = Path(tmp_dir) / "posre_GRA.itp"
            topol_top = Path(tmp_dir) / "topol.top"

            for path in (
                graphene_gro,
                graphene_itp,
                graphene_1layer,
                graphite_gro,
                graphite_pdb,
                posre_itp,
                topol_top,
            ):
                self.assertTrue(path.exists(), msg=f"Missing {path}")

            single_layer_atoms = read_gro_atom_count(graphene_gro)
            stacked_atoms = read_gro_atom_count(graphite_gro)
            self.assertEqual(stacked_atoms, single_layer_atoms * layers)

            final_box = graphite_gro.read_text().splitlines()[-1].split()
            self.assertEqual(len(final_box), 3)
            self.assertAlmostEqual(float(final_box[2]), layers * spacing, places=3)

            itp_text = graphene_itp.read_text()
            self.assertIn('#include "posre_GRA.itp"', itp_text)
            self.assertIn("#ifdef POSRES", itp_text)

            posre_entries = [
                line
                for line in posre_itp.read_text().splitlines()
                if line.strip() and not line.strip().startswith(";") and not line.strip().startswith("[")
            ]
            self.assertEqual(len(posre_entries), single_layer_atoms)

            topol_text = topol_top.read_text()
            self.assertIn('#include "martini_v3.0.0.itp"', topol_text)
            self.assertIn('#include "graphene.itp"', topol_text)
            self.assertIn(f"GRA        {layers}", topol_text)

    def test_overwrite_requires_force(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            first = run_builder("-x", "3", "-y", "3", "--output-dir", tmp_dir, "--quiet")
            self.assertEqual(first.returncode, 0, msg=first.stderr)

            second = run_builder("-x", "3", "-y", "3", "--output-dir", tmp_dir, "--quiet")
            self.assertNotEqual(second.returncode, 0)
            self.assertIn("Refusing to overwrite existing output file", second.stderr)

            forced = run_builder(
                "-x",
                "3",
                "-y",
                "3",
                "--output-dir",
                tmp_dir,
                "--force",
                "--quiet",
            )
            self.assertEqual(forced.returncode, 0, msg=forced.stderr)


if __name__ == "__main__":
    unittest.main()
