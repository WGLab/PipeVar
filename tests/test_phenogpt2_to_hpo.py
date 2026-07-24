import importlib.util
import json
import pickle
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path("/home/beoungle/docker_work/phenogpt2/phenogpt2_to_hpo.py")


spec = importlib.util.spec_from_file_location("phenogpt2_to_hpo", SCRIPT)
phenogpt2_to_hpo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(phenogpt2_to_hpo)


class PhenoGPT2ToHpoTests(unittest.TestCase):
    def make_model(self, root, shard="model-00001-of-00001.safetensors"):
        for name in phenogpt2_to_hpo.REQUIRED_MODEL_FILES[:-1]:
            (root / name).write_text("metadata")
        (root / shard).write_bytes(b"weights")
        (root / "model.safetensors.index.json").write_text(
            json.dumps({"weight_map": {"model.weight": shard}})
        )

    def make_embedding_model(self, root):
        for name in phenogpt2_to_hpo.MODEL_REQUIRED_FILES["embedding"]:
            if name == "model.safetensors":
                (root / name).write_bytes(b"weights")
            else:
                (root / name).write_text("metadata")

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
            phenogpt2_to_hpo.write_manifest(model)

            fingerprint = phenogpt2_to_hpo.validate_model(model)

            self.assertRegex(fingerprint, r"^[0-9a-f]{64}$")

    def test_validate_negation_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "negation"
            model.mkdir()
            self.make_model(model)
            phenogpt2_to_hpo.write_manifest(model)

            fingerprint = phenogpt2_to_hpo.validate_model(model, kind="negation")

            self.assertRegex(fingerprint, r"^[0-9a-f]{64}$")

    def test_validate_embedding_checkpoint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "embedding"
            model.mkdir()
            self.make_embedding_model(model)
            phenogpt2_to_hpo.write_manifest(model)

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
            phenogpt2_to_hpo.write_manifest(model)
            fingerprint = phenogpt2_to_hpo.validate_model(model, kind="embedding")

            snapshot = phenogpt2_to_hpo.prepare_embedding_cache(
                model, cache, fingerprint
            )

            self.assertTrue(snapshot.is_symlink())
            self.assertEqual(snapshot.resolve(), model.resolve())
            self.assertEqual(
                (
                    cache
                    / phenogpt2_to_hpo.EMBEDDING_CACHE_REPO
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
            phenogpt2_to_hpo.write_manifest(model)
            fingerprint = phenogpt2_to_hpo.validate_model(model, kind="embedding")
            snapshot = (
                cache
                / phenogpt2_to_hpo.EMBEDDING_CACHE_REPO
                / "snapshots"
                / fingerprint[:40]
            )
            snapshot.mkdir(parents=True)

            with self.assertRaisesRegex(ValueError, "conflicts"):
                phenogpt2_to_hpo.prepare_embedding_cache(
                    model, cache, fingerprint
                )

    def test_prepare_embedding_cache_rejects_conflicting_ref(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            model = tmp / "embedding"
            cache = tmp / "cache"
            model.mkdir()
            cache.mkdir()
            self.make_embedding_model(model)
            phenogpt2_to_hpo.write_manifest(model)
            fingerprint = phenogpt2_to_hpo.validate_model(model, kind="embedding")
            refs = cache / phenogpt2_to_hpo.EMBEDDING_CACHE_REPO / "refs"
            refs.mkdir(parents=True)
            (refs / "main").write_text("0" * 40 + "\n")

            with self.assertRaisesRegex(ValueError, "different model revision"):
                phenogpt2_to_hpo.prepare_embedding_cache(
                    model, cache, fingerprint
                )

    def test_validate_model_rejects_missing_shard(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            (model / "model-00001-of-00001.safetensors").unlink()
            phenogpt2_to_hpo.write_manifest(model)

            with self.assertRaisesRegex(ValueError, "required model file is missing"):
                phenogpt2_to_hpo.validate_model(model)

    def test_validate_model_rejects_malformed_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            (model / "model.safetensors.index.json").write_text("not-json")
            phenogpt2_to_hpo.write_manifest(model)

            with self.assertRaisesRegex(ValueError, "invalid model index"):
                phenogpt2_to_hpo.validate_model(model)

    def test_validate_model_fingerprint_changes_with_shard(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            phenogpt2_to_hpo.write_manifest(model)
            before = phenogpt2_to_hpo.validate_model(model)
            (model / "model-00001-of-00001.safetensors").write_bytes(b"new-weights")
            phenogpt2_to_hpo.write_manifest(model)

            after = phenogpt2_to_hpo.validate_model(model)

            self.assertNotEqual(before, after)

    def test_validate_model_rejects_stale_manifest(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            phenogpt2_to_hpo.write_manifest(model)
            (model / "model-00001-of-00001.safetensors").write_bytes(b"changed")

            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                phenogpt2_to_hpo.validate_model(model)

    def test_validate_model_rejects_unmanifested_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            phenogpt2_to_hpo.write_manifest(model)
            (model / "generation_config.json").write_text("{}")

            with self.assertRaisesRegex(ValueError, "unmanifested"):
                phenogpt2_to_hpo.validate_model(model)

    def test_write_manifest_is_atomic_and_leaves_no_temporary_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)

            manifest = phenogpt2_to_hpo.write_manifest(model)

            self.assertTrue(manifest.is_file())
            self.assertEqual(list(model.glob(f".{phenogpt2_to_hpo.MANIFEST_NAME}.*")), [])

    def test_validate_model_rejects_traversal_in_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            (model / "model.safetensors.index.json").write_text(
                json.dumps({"weight_map": {"model.weight": "../weights.safetensors"}})
            )
            phenogpt2_to_hpo.write_manifest(model)

            with self.assertRaisesRegex(ValueError, "unsafe model path"):
                phenogpt2_to_hpo.validate_model(model)

    def test_validate_model_rejects_escaping_symlink(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            model = tmp / "model-v1"
            model.mkdir()
            self.make_model(model)
            outside = tmp / "outside.safetensors"
            outside.write_bytes(b"weights")
            shard = model / "model-00001-of-00001.safetensors"
            shard.unlink()
            shard.symlink_to(outside)

            with self.assertRaisesRegex(ValueError, "materialized files"):
                phenogpt2_to_hpo.write_manifest(model)


if __name__ == "__main__":
    unittest.main()
