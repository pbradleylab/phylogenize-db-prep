"""Exercise database downloads with local archives and a fake wget, without network access."""

import gzip
import importlib.util
import io
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
DATABASES = {
    "eggnog.db": b"annotation database\n",
    "eggnog_proteins.dmnd": b"diamond database\n",
    "eggnog.taxa.db": b"taxonomy database\n",
    "eggnog.taxa.db.traverse.pkl": b"taxonomy traversal\n",
}


@unittest.skipUnless(importlib.util.find_spec("eggnogmapper"), "eggNOG-mapper is required")
class EggnogDownloadTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="eggnog downloads ' ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.data_dir = self.directory / "data"
        self.data_dir.mkdir()
        self.archives = self.directory / "archives"
        self.archives.mkdir()
        self.requests = self.directory / "requests.txt"
        for name in ["eggnog.db", "eggnog_proteins.dmnd"]:
            (self.archives / (name + ".gz")).write_bytes(gzip.compress(DATABASES[name]))
        with tarfile.open(self.archives / "eggnog.taxa.tar.gz", "w:gz") as archive:
            for name in ["eggnog.taxa.db", "eggnog.taxa.db.traverse.pkl"]:
                entry = tarfile.TarInfo(name)
                entry.size = len(DATABASES[name])
                archive.addfile(entry, io.BytesIO(DATABASES[name]))
        self.bin_dir = self.directory / "bin"
        self.bin_dir.mkdir()
        wget = self.bin_dir / "wget"
        wget.write_text(f"#!{sys.executable}\n" + '''\
import os
from pathlib import Path
import shutil
import sys

url = sys.argv[-1]
with open(os.environ["EGGNOG_TEST_REQUESTS"], "a") as handle:
    handle.write(url + "\\n")
if os.environ.get("EGGNOG_TEST_FAIL"):
    sys.exit(8)
output = sys.argv[sys.argv.index("-O") + 1]
shutil.copyfile(Path(os.environ["EGGNOG_TEST_ARCHIVES"]) / url.rsplit("/", 1)[-1], output)
''')
        wget.chmod(0o755)

    def run_download(self, *options, fail=False):
        result = subprocess.run([
            sys.executable, "-B", str(REPOSITORY / "workflow/scripts/download_eggnog_data.py"),
            "-y", "--data_dir", str(self.data_dir), *options,
        ], env={
            **os.environ, "PATH": str(self.bin_dir) + os.pathsep + os.environ["PATH"],
            "EGGNOG_TEST_REQUESTS": str(self.requests), "EGGNOG_TEST_ARCHIVES": str(self.archives),
            "EGGNOG_TEST_FAIL": "1" if fail else "",
        }, capture_output=True, text=True, timeout=30)
        return result

    def requested_archives(self):
        if not self.requests.exists():
            return []
        return [url.rsplit("/", 1)[-1] for url in self.requests.read_text().splitlines()]

    def cache_databases(self):
        for name, contents in DATABASES.items():
            (self.data_dir / name).write_bytes(contents)

    def test_default_downloads_all_required_databases(self):
        result = self.run_download()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertCountEqual(self.requested_archives(), [
            "eggnog.db.gz", "eggnog.taxa.tar.gz", "eggnog_proteins.dmnd.gz",
        ])
        for name, contents in DATABASES.items():
            self.assertEqual((self.data_dir / name).read_bytes(), contents)

    def test_missing_diamond_download_preserves_other_databases(self):
        self.cache_databases()
        (self.data_dir / "eggnog_proteins.dmnd").unlink()
        (self.data_dir / "eggnog_proteins.dmnd.gz").touch()
        cached = {name: (self.data_dir / name).stat().st_mtime_ns for name in DATABASES
                  if name != "eggnog_proteins.dmnd"}
        result = self.run_download("--db", "diamond")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.requested_archives(), ["eggnog_proteins.dmnd.gz"])
        for name, contents in DATABASES.items():
            self.assertEqual((self.data_dir / name).read_bytes(), contents)
        for name, mtime in cached.items():
            self.assertEqual((self.data_dir / name).stat().st_mtime_ns, mtime)

    def test_complete_databases_are_reused(self):
        self.cache_databases()
        result = self.run_download()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.requested_archives(), [])

    def test_missing_taxonomy_traversal_is_restored(self):
        self.cache_databases()
        (self.data_dir / "eggnog.taxa.db.traverse.pkl").unlink()
        result = self.run_download("--db", "taxa")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.requested_archives(), ["eggnog.taxa.tar.gz"])
        self.assertEqual((self.data_dir / "eggnog.taxa.db.traverse.pkl").read_bytes(),
                         DATABASES["eggnog.taxa.db.traverse.pkl"])

    def test_empty_diamond_is_replaced(self):
        (self.data_dir / "eggnog_proteins.dmnd").touch()
        result = self.run_download("--db", "diamond")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual((self.data_dir / "eggnog_proteins.dmnd").read_bytes(),
                         DATABASES["eggnog_proteins.dmnd"])

    def test_download_failure_exits_with_an_error(self):
        result = self.run_download("--db", "diamond", fail=True)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse((self.data_dir / "eggnog_proteins.dmnd").exists())
        self.assertNotIn("Finished.", result.stdout)

    def test_decompression_failure_exits_with_an_error(self):
        (self.archives / "eggnog_proteins.dmnd.gz").write_bytes(b"invalid gzip archive")
        result = self.run_download("--db", "diamond")
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotIn("Finished.", result.stdout)


if __name__ == "__main__":
    unittest.main()
