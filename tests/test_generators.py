import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = (
    ("non_periodic", REPO_ROOT / "Non-Periodic" / "martini3-graphene-topology.py"),
    ("periodic", REPO_ROOT / "Periodic" / "martini3-graphene-periodic.py"),
)
BASE_ARGS = ("-x", "3", "-y", "3", "-z", "5")
TIMEOUT_SECONDS = 120


def run_script(script_path, *extra_args):
    cmd = [sys.executable, str(script_path), *BASE_ARGS, *extra_args]
    return subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        timeout=TIMEOUT_SECONDS,
    )


def parse_itp_section_counts(itp_path):
    counts = {}
    section = None
    for raw_line in itp_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line
            counts.setdefault(section, 0)
            continue
        if section is not None:
            counts[section] += 1
    return counts


class GeneratorRegressionTests(unittest.TestCase):
    def test_generates_outputs_for_both_scripts(self):
        required_sections = (
            "[ atoms ]",
            "[ bonds ]",
            "[ angles ]",
            "[ dihedrals ]",
            "[ virtual_sitesn ]",
            "[ exclusions ]",
        )
        for label, script_path in SCRIPTS:
            with self.subTest(script=label):
                with tempfile.TemporaryDirectory() as tmp_dir:
                    result = run_script(
                        script_path,
                        "-o",
                        "graphene_test",
                        "--output-dir",
                        tmp_dir,
                    )
                    self.assertEqual(result.returncode, 0, msg=result.stderr)
                    gro_path = Path(tmp_dir) / "graphene_test.gro"
                    itp_path = Path(tmp_dir) / "graphene_test.itp"
                    self.assertTrue(gro_path.exists(), msg=f"Missing {gro_path}")
                    self.assertTrue(itp_path.exists(), msg=f"Missing {itp_path}")

                    counts = parse_itp_section_counts(itp_path)
                    for section in required_sections:
                        self.assertIn(section, counts, msg=f"Missing section {section}")
                        self.assertGreater(counts[section], 0, msg=f"Section {section} is empty")

    def test_overwrite_requires_force(self):
        for label, script_path in SCRIPTS:
            with self.subTest(script=label):
                with tempfile.TemporaryDirectory() as tmp_dir:
                    first = run_script(
                        script_path,
                        "-o",
                        "graphene_test",
                        "--output-dir",
                        tmp_dir,
                    )
                    self.assertEqual(first.returncode, 0, msg=first.stderr)

                    second = run_script(
                        script_path,
                        "-o",
                        "graphene_test",
                        "--output-dir",
                        tmp_dir,
                    )
                    self.assertNotEqual(second.returncode, 0)
                    self.assertIn("Refusing to overwrite existing output file", second.stderr)

                    forced = run_script(
                        script_path,
                        "-o",
                        "graphene_test",
                        "--output-dir",
                        tmp_dir,
                        "--force",
                    )
                    self.assertEqual(forced.returncode, 0, msg=forced.stderr)

    def test_profile_outputs_timings(self):
        _, script_path = SCRIPTS[1]
        with tempfile.TemporaryDirectory() as tmp_dir:
            result = run_script(
                script_path,
                "-o",
                "graphene_test",
                "--output-dir",
                tmp_dir,
                "--profile",
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertIn("[profile]", result.stderr)


if __name__ == "__main__":
    unittest.main()
