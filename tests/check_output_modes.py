#!/usr/bin/env python3
"""Check compression-mode preflight through the public pair command."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "tests" / "data" / "pair"
BIN = os.environ.get("TEVOX_BIN", str(REPO / "tevox"))
TABLES = (
    "evidence", "observation_scores", "candidates", "candidate_features",
    "candidate_contexts", "decisions", "edges", "relations", "solver",
    "loci", "instances", "states", "summary", "synteny.blocks",
    "synteny.anchors", "contexts", "te_contexts",
)


def run_pair(prefix, compressed):
    args = [
        BIN, "pair", "--genome-a", "A", "--fasta-a", str(DATA / "A.fa"),
        "--te-a", str(DATA / "A.gff3"), "--genome-b", "B",
        "--fasta-b", str(DATA / "B.fa"), "--te-b", str(DATA / "B.gff3"),
        "--paf", str(DATA / "A_B.paf"), "--output", str(prefix),
    ]
    if compressed:
        args.append("--gzip-output")
    return subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, check=False)


def snapshot(directory):
    entries = {}
    for path in sorted(directory.rglob("*")):
        st = path.lstat()
        if path.is_symlink():
            content = os.readlink(path)
        elif path.is_file():
            content = hashlib.sha256(path.read_bytes()).hexdigest()
        else:
            content = None
        entries[str(path.relative_to(directory))] = (
            st.st_mode, st.st_ino, st.st_size, st.st_mtime_ns, content,
        )
    return entries


class OutputModeTests(unittest.TestCase):
    def assert_rejected_without_changes(self, directory, prefix, compressed):
        before = snapshot(directory)
        result = run_pair(prefix, compressed)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn("opposite compression mode", result.stderr)
        self.assertIn("new output prefix or directory", result.stderr)
        self.assertEqual(snapshot(directory), before)

    def test_completed_run_rejects_opposite_mode(self):
        for initial_compressed in (False, True):
            with self.subTest(initial_compressed=initial_compressed):
                with tempfile.TemporaryDirectory() as tmp:
                    directory = Path(tmp)
                    prefix = directory / "run ; literal"
                    first = run_pair(prefix, initial_compressed)
                    self.assertEqual(first.returncode, 0, first.stderr)
                    self.assertEqual(len(list(directory.iterdir())), 18)
                    self.assert_rejected_without_changes(
                        directory, prefix, not initial_compressed)

    def test_each_partial_table_rejects_before_writing(self):
        for compressed in (False, True):
            for table in TABLES:
                with self.subTest(compressed=compressed, table=table):
                    with tempfile.TemporaryDirectory() as tmp:
                        directory = Path(tmp)
                        prefix = directory / "run"
                        opposite = ".tsv" if compressed else ".tsv.gz"
                        current = ".tsv.gz" if compressed else ".tsv"
                        Path(str(prefix) + "." + table + opposite).write_bytes(
                            b"existing partial table\n")
                        Path(str(prefix) + ".evidence" + current).write_bytes(
                            b"existing current-mode table\n")
                        Path(str(prefix) + ".run.json").write_bytes(
                            b'{"previous_run": true}\n')
                        self.assert_rejected_without_changes(
                            directory, prefix, compressed)

    def test_opposite_directory_and_dangling_symlink_are_conflicts(self):
        for compressed in (False, True):
            for kind in ("directory", "dangling_symlink"):
                with self.subTest(compressed=compressed, kind=kind):
                    with tempfile.TemporaryDirectory() as tmp:
                        directory = Path(tmp)
                        prefix = directory / "run"
                        opposite = ".tsv" if compressed else ".tsv.gz"
                        path = Path(str(prefix) + ".te_contexts" + opposite)
                        if kind == "directory":
                            path.mkdir()
                            (path / "keep").write_bytes(b"keep\n")
                        else:
                            path.symlink_to("missing-target")
                        self.assert_rejected_without_changes(
                            directory, prefix, compressed)

    def test_same_mode_rerun_remains_supported(self):
        for compressed in (False, True):
            with self.subTest(compressed=compressed):
                with tempfile.TemporaryDirectory() as tmp:
                    directory = Path(tmp)
                    prefix = directory / "run"
                    for _ in range(2):
                        result = run_pair(prefix, compressed)
                        self.assertEqual(result.returncode, 0, result.stderr)
                    marker = json.loads(Path(str(prefix) + ".run.json").read_text())
                    extension = ".tsv.gz" if compressed else ".tsv"
                    self.assertEqual(set(marker["outputs"]),
                                     {table + extension for table in TABLES})
                    self.assertEqual(marker["output_compression"],
                                     "gzip" if compressed else "none")
                    self.assertEqual(len(list(directory.iterdir())), 18)


if __name__ == "__main__":
    unittest.main()
