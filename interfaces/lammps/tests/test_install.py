"""Installer checks without downloading or compiling LAMMPS."""
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("lammps_install", HERE / "install.py")
INSTALLER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(INSTALLER)
REPO = HERE.parents[1]


class InstallTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "src").mkdir()
        (self.root / "cmake/Modules/Packages").mkdir(parents=True)
        self.version("22 Jul 2025", "Update 6")
        (self.root / "cmake/CMakeLists.txt").write_text(
            "set(STANDARD_PACKAGES\n  KSPACE\n  PYTHON\n)\n"
            "foreach(PKG_WITH_INCL KSPACE PYTHON\n  ML-HDNNP)\nendforeach()\n"
            "foreach(PKG_WITH_INCL OPENMP GPU)\nendforeach()\n")
        (self.root / "src/atom.cpp").write_text(
            "#ifdef LMP_GPU\n  if (userbinsize == 0.0) {\n  }\n#endif\n" * 2)

    def version(self, release, update):
        (self.root / "src/version.h").write_text(
            f'#define LAMMPS_VERSION "{release}"\n#define LAMMPS_UPDATE "{update}"\n')

    def snapshot(self):
        return {str(p.relative_to(self.root)): p.read_bytes()
                for p in self.root.rglob("*") if p.is_file()}

    def test_supported_releases_and_idempotent_install(self):
        for release, update in sorted(INSTALLER.SUPPORTED):
            with self.subTest(release=release):
                self.version(release, update)
                self.assertEqual(INSTALLER.install(self.root), f"{release} {update}")
                first = self.snapshot()
                INSTALLER.install(self.root)
                self.assertEqual(first, self.snapshot())
                for name in ("accelnet.h", "accelnet_target.h"):
                    self.assertEqual((self.root / "src/ACCELNET" / name).read_bytes(),
                                     (REPO / "AccelNetPredictor/include" / name).read_bytes())
                cmake = (self.root / "cmake/CMakeLists.txt").read_text()
                self.assertEqual(cmake.count("ACCELNET"), 2)
                self.assertEqual((self.root / "src/atom.cpp").read_text().count("!domain->triclinic"), 2)
                self.assertTrue((self.root / "lib/gpu/lal_accelnet_ext.cpp").is_file())

    def test_unsupported_version_does_not_modify_tree(self):
        for release, update in (("2 Sep 2026", ""), ("22 Jul 2025", "Update 5")):
            with self.subTest(release=release, update=update):
                self.version(release, update)
                before = self.snapshot()
                with self.assertRaisesRegex(ValueError, "stable_22Jul2025_update6"):
                    INSTALLER.install(self.root)
                self.assertEqual(before, self.snapshot())

    def test_bad_cmake_does_not_partially_install(self):
        path = self.root / "cmake/CMakeLists.txt"
        path.write_text(path.read_text().replace("PKG_WITH_INCL KSPACE", "OTHER_LIST KSPACE"))
        before = self.snapshot()
        with self.assertRaisesRegex(ValueError, "registration"):
            INSTALLER.install(self.root)
        self.assertEqual(before, self.snapshot())

    def test_bad_sorting_code_does_not_partially_install(self):
        (self.root / "src/atom.cpp").write_text("// Unrecognized implementation\n")
        before = self.snapshot()
        with self.assertRaisesRegex(ValueError, "atom-sorting"):
            INSTALLER.install(self.root)
        self.assertEqual(before, self.snapshot())

    def test_legacy_entry_point(self):
        self.version("29 Aug 2024", "Update 4")
        result = subprocess.run([sys.executable, str(HERE / "29Aug2024/install.py"),
                                 str(self.root)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("29 Aug 2024 Update 4", result.stdout)


if __name__ == "__main__":
    unittest.main()
