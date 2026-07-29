import importlib.util
import json
import pickle
import tempfile
import unittest
from pathlib import Path
from unittest import mock


SCRIPT = Path("/home/beoungle/docker_work/phenogpt2/phenogpt2_to_hpo.py")


spec = importlib.util.spec_from_file_location("phenogpt2_to_hpo", SCRIPT)
phenogpt2_to_hpo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(phenogpt2_to_hpo)


class PhenoGPT2ToHpoTests(unittest.TestCase):
    def setUp(self):
        patcher = mock.patch.object(phenogpt2_to_hpo, "_probe_model_loader")
        self.loader_probe = patcher.start()
        self.addCleanup(patcher.stop)

    def make_model(self, root, weight_name="weights.custom"):
        (root / "runtime-metadata.custom").write_text("metadata")
        (root / weight_name).write_bytes(b"weights")

    def make_embedding_model(self, root):
        self.make_model(root, weight_name="embedding-weights.custom")

    def make_upstream_source(self, root, repo_id="vendor/embedding-model"):
        root.mkdir()
        (root / "helpers.py").write_text(
            "from sentence_transformers import SentenceTransformer\n"
            f"MODEL_ID = {repo_id!r}\n"
            "def load_embedding():\n"
            "    return SentenceTransformer(MODEL_ID)\n"
        )
        return root

    def test_make_input_uses_json_null_image(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            note = tmp / "note.txt"
            out_json = tmp / "input.json"
            note.write_text("Patient has seizures.")

            phenogpt2_to_hpo.make_input(note, "sample1", out_json)

            payload = json.loads(out_json.read_text())
            self.assertEqual(payload["sample1"]["clinical_note"], "Patient has seizures.")
            self.assertIsNone(payload["sample1"]["image"])
            self.assertEqual(payload["sample1"]["pid"], "sample1")

    def test_convert_prefers_filtered_when_requested(self):
        result = {
            "sample1": {
                "text": {
                    "complete": True,
                    "phenotypes": {
                        "seizure": {"HPO_ID": "HP:0001250"},
                        "ataxia": {"HPO_ID": "HP:0001251"},
                    },
                    "filtered_phenotypes": {
                        "seizure": {"HPO_ID": "HP:0001250"},
                    },
                    "negation_analysis": {
                        "demographics": {},
                        "phenotypes": {"seizure": {"correct": "True"}},
                    },
                },
                "image": {},
            }
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result_path = tmp / "phenogpt2_rep0.pkl"
            out_hpo = tmp / "hpo.txt"
            with result_path.open("wb") as handle:
                pickle.dump(result, handle)

            phenogpt2_to_hpo.convert(result_path, out_hpo, prefer_filtered=True)

            self.assertEqual(out_hpo.read_text(), "HP:0001250\n")

    def test_convert_treats_empty_filtered_phenotypes_as_authoritative(self):
        result = {
            "sample1": {
                "text": {
                    "complete": True,
                    "phenotypes": {
                        "seizure": {"HPO_ID": "HP:0001250"},
                        "ataxia": "HP:0001251",
                    },
                    "filtered_phenotypes": {},
                    "negation_analysis": {
                        "demographics": {},
                        "phenotypes": {
                            "seizure": {"correct": "False"},
                            "ataxia": {"correct": "False"},
                        },
                    },
                }
            }
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result_path = tmp / "phenogpt2_rep0.json"
            out_hpo = tmp / "hpo.txt"
            result_path.write_text(json.dumps(result))

            phenogpt2_to_hpo.convert(result_path, out_hpo, prefer_filtered=True)

            self.assertEqual(out_hpo.read_text(), "")

    def test_convert_rejects_incomplete_negation_output(self):
        result = {
            "sample1": {
                "text": {
                    "complete": False,
                    "phenotypes": {"seizure": {"HPO_ID": "HP:0001250"}},
                    "filtered_phenotypes": {},
                    "negation_analysis": {"demographics": {}, "phenotypes": {}},
                }
            }
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result_path = tmp / "phenogpt2_rep0.json"
            out_hpo = tmp / "hpo.txt"
            result_path.write_text(json.dumps(result))

            with self.assertRaisesRegex(ValueError, "incomplete"):
                phenogpt2_to_hpo.convert(result_path, out_hpo, prefer_filtered=True)

    def test_convert_rejects_missing_negation_audit(self):
        result = {
            "sample1": {
                "text": {
                    "complete": True,
                    "phenotypes": {"seizure": {"HPO_ID": "HP:0001250"}},
                    "filtered_phenotypes": {},
                }
            }
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result_path = tmp / "phenogpt2_rep0.json"
            out_hpo = tmp / "hpo.txt"
            result_path.write_text(json.dumps(result))

            with self.assertRaisesRegex(ValueError, "valid negation_analysis"):
                phenogpt2_to_hpo.convert(result_path, out_hpo, prefer_filtered=True)

    def test_convert_accepts_empty_negation_candidates(self):
        result = {
            "sample1": {
                "text": {
                    "complete": True,
                    "phenotypes": {},
                    "filtered_phenotypes": {},
                    "negation_analysis": {"demographics": {}, "phenotypes": {}},
                }
            }
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result_path = tmp / "phenogpt2_rep0.json"
            out_hpo = tmp / "hpo.txt"
            result_path.write_text(json.dumps(result))

            phenogpt2_to_hpo.convert(result_path, out_hpo, prefer_filtered=True)

            self.assertEqual(out_hpo.read_text(), "")

    def test_convert_rejects_missing_filtered_phenotypes(self):
        result = {
            "sample1": {
                "text": {
                    "complete": True,
                    "phenotypes": {"seizure": {"HPO_ID": "HP:0001250"}},
                    "negation_analysis": {"demographics": {}, "phenotypes": {}},
                }
            }
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result_path = tmp / "phenogpt2_rep0.json"
            out_hpo = tmp / "hpo.txt"
            result_path.write_text(json.dumps(result))

            with self.assertRaisesRegex(ValueError, "missing filtered_phenotypes"):
                phenogpt2_to_hpo.convert(result_path, out_hpo, prefer_filtered=True)

    def test_convert_deduplicates_recursive_hpo_ids(self):
        result = {
            "sample1": {
                "text": {
                    "phenotypes": {
                        "term": {
                            "HPO_ID": "HP:0001250",
                            "evidence": "HP:0001250 and HP:0004322",
                        }
                    }
                }
            }
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result_path = tmp / "phenogpt2_rep0.pkl"
            out_hpo = tmp / "hpo.txt"
            with result_path.open("wb") as handle:
                pickle.dump(result, handle)

            phenogpt2_to_hpo.convert(result_path, out_hpo, prefer_filtered=False)

            self.assertEqual(out_hpo.read_text(), "HP:0001250\nHP:0004322\n")

    def test_validate_model_accepts_complete_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)

            fingerprint = phenogpt2_to_hpo.validate_model(model)

            self.assertRegex(fingerprint, r"^[0-9a-f]{64}$")

    def test_validate_negation_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "negation"
            model.mkdir()
            self.make_model(model)

            fingerprint = phenogpt2_to_hpo.validate_model(model, kind="negation")

            self.assertRegex(fingerprint, r"^[0-9a-f]{64}$")

    def test_validate_embedding_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "embedding"
            model.mkdir()
            self.make_embedding_model(model)

            fingerprint = phenogpt2_to_hpo.validate_model(model, kind="embedding")

            self.assertRegex(fingerprint, r"^[0-9a-f]{64}$")

    def test_prepare_embedding_cache_exposes_standalone_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            model = tmp / "embedding"
            cache = tmp / "cache"
            model.mkdir()
            cache.mkdir()
            self.make_embedding_model(model)
            source = self.make_upstream_source(tmp / "app")
            fingerprint = phenogpt2_to_hpo.validate_model(model, kind="embedding")

            snapshot = phenogpt2_to_hpo.prepare_embedding_cache(
                model, cache, fingerprint, source
            )

            self.assertTrue(snapshot.is_symlink())
            self.assertEqual(snapshot.resolve(), model.resolve())
            self.assertEqual(
                (
                    cache
                    / "models--vendor--embedding-model"
                    / "refs"
                    / "main"
                ).read_text(),
                f"{fingerprint[:40]}\n",
            )

    def test_prepare_embedding_cache_rejects_conflicting_snapshot(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            model = tmp / "embedding"
            cache = tmp / "cache"
            model.mkdir()
            cache.mkdir()
            self.make_embedding_model(model)
            source = self.make_upstream_source(tmp / "app")
            fingerprint = phenogpt2_to_hpo.validate_model(model, kind="embedding")
            snapshot = (
                cache
                / "models--vendor--embedding-model"
                / "snapshots"
                / fingerprint[:40]
            )
            snapshot.mkdir(parents=True)

            with self.assertRaisesRegex(ValueError, "conflicts"):
                phenogpt2_to_hpo.prepare_embedding_cache(
                    model, cache, fingerprint, source
                )

    def test_prepare_embedding_cache_rejects_conflicting_ref(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            model = tmp / "embedding"
            cache = tmp / "cache"
            model.mkdir()
            cache.mkdir()
            self.make_embedding_model(model)
            source = self.make_upstream_source(tmp / "app")
            fingerprint = phenogpt2_to_hpo.validate_model(model, kind="embedding")
            refs = cache / "models--vendor--embedding-model" / "refs"
            refs.mkdir(parents=True)
            (refs / "main").write_text("0" * 40 + "\n")

            with self.assertRaisesRegex(ValueError, "different model revision"):
                phenogpt2_to_hpo.prepare_embedding_cache(
                    model, cache, fingerprint, source
                )

    def test_validate_model_accepts_checkpoint_with_incidental_empty_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            (model / "download.lock").touch()

            fingerprint = phenogpt2_to_hpo.validate_model(model)

            self.assertRegex(fingerprint, r"^[0-9a-f]{64}$")
            self.loader_probe.assert_called_once_with(model.resolve(), "main")

    def test_discovers_embedding_repository_from_upstream_source(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = self.make_upstream_source(Path(tmpdir) / "app", "org/model-v2")

            repo_id = phenogpt2_to_hpo.discover_embedding_repo_id(source)

            self.assertEqual(repo_id, "org/model-v2")

    def test_rejects_ambiguous_embedding_repositories(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "app"
            self.make_upstream_source(source, "org/model-a")
            (source / "other.py").write_text(
                "from sentence_transformers import SentenceTransformer\n"
                "SentenceTransformer('org/model-b')\n"
            )

            with self.assertRaisesRegex(ValueError, "multiple"):
                phenogpt2_to_hpo.discover_embedding_repo_id(source)

    def test_validate_model_rejects_empty_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()

            with self.assertRaisesRegex(ValueError, "has no files"):
                phenogpt2_to_hpo.validate_model(model)

    def test_validate_model_defers_layout_to_loader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model, weight_name="unexpected-layout.weights")
            self.loader_probe.side_effect = ValueError("loader rejected checkpoint")

            with self.assertRaisesRegex(ValueError, "loader rejected checkpoint"):
                phenogpt2_to_hpo.validate_model(model)

    def test_validate_model_fingerprint_changes_with_content(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            before = phenogpt2_to_hpo.validate_model(model)
            (model / "weights.custom").write_bytes(b"new-weights")

            after = phenogpt2_to_hpo.validate_model(model)

            self.assertNotEqual(before, after)

    def test_resolve_model_directory_supports_huggingface_cache(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = Path(tmpdir) / "models--vendor--model"
            revision = "a" * 40
            snapshot = cache / "snapshots" / revision
            blob = cache / "blobs" / "content"
            snapshot.mkdir(parents=True)
            blob.parent.mkdir()
            blob.write_bytes(b"weights")
            (cache / "refs").mkdir()
            (cache / "refs" / "main").write_text(f"{revision}\n")
            (snapshot / "weights-any-name").symlink_to("../../blobs/content")

            resolved = phenogpt2_to_hpo.resolve_model_directory(cache)
            fingerprint = phenogpt2_to_hpo.validate_model(cache)

            self.assertEqual(resolved, snapshot.resolve())
            self.assertRegex(fingerprint, r"^[0-9a-f]{64}$")
            self.loader_probe.assert_called_once_with(snapshot.resolve(), "main")

    def test_validate_model_rejects_incomplete_huggingface_cache(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = Path(tmpdir) / "models--vendor--model"
            revision = "b" * 40
            (cache / "snapshots" / revision).mkdir(parents=True)
            (cache / "refs").mkdir()
            (cache / "refs" / "main").write_text(f"{revision}\n")

            with self.assertRaisesRegex(ValueError, "has no files"):
                phenogpt2_to_hpo.validate_model(cache)

    def test_validate_model_rejects_escaping_symlink(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            model = tmp / "model-v1"
            model.mkdir()
            self.make_model(model)
            outside = tmp / "outside.safetensors"
            outside.write_bytes(b"weights")
            weight = model / "weights.custom"
            weight.unlink()
            weight.symlink_to(outside)

            with self.assertRaisesRegex(ValueError, "outside its mount"):
                phenogpt2_to_hpo.validate_model(model)


if __name__ == "__main__":
    unittest.main()
