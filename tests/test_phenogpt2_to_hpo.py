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
                    "phenotypes": {
                        "seizure": {"HPO_ID": "HP:0001250"},
                        "ataxia": {"HPO_ID": "HP:0001251"},
                    },
                    "filtered_phenotypes": {
                        "seizure": {"HPO_ID": "HP:0001250"},
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

    def test_convert_falls_back_to_phenotypes(self):
        result = {
            "sample1": {
                "text": {
                    "phenotypes": {
                        "seizure": {"HPO_ID": "HP:0001250"},
                        "ataxia": "HP:0001251",
                    },
                    "filtered_phenotypes": {},
                }
            }
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            result_path = tmp / "phenogpt2_rep0.json"
            out_hpo = tmp / "hpo.txt"
            result_path.write_text(json.dumps(result))

            phenogpt2_to_hpo.convert(result_path, out_hpo, prefer_filtered=True)

            self.assertEqual(out_hpo.read_text(), "HP:0001250\nHP:0001251\n")

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

    def test_validate_model_rejects_missing_shard(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            (model / "model-00001-of-00001.safetensors").unlink()

            with self.assertRaisesRegex(ValueError, "required model file is missing"):
                phenogpt2_to_hpo.validate_model(model)

    def test_validate_model_rejects_malformed_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            (model / "model.safetensors.index.json").write_text("not-json")

            with self.assertRaisesRegex(ValueError, "invalid model index"):
                phenogpt2_to_hpo.validate_model(model)

    def test_validate_model_fingerprint_changes_with_shard(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            before = phenogpt2_to_hpo.validate_model(model)
            (model / "model-00001-of-00001.safetensors").write_bytes(b"new-weights")

            after = phenogpt2_to_hpo.validate_model(model)

            self.assertNotEqual(before, after)

    def test_validate_model_rejects_traversal_in_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model = Path(tmpdir) / "model-v1"
            model.mkdir()
            self.make_model(model)
            (model / "model.safetensors.index.json").write_text(
                json.dumps({"weight_map": {"model.weight": "../weights.safetensors"}})
            )

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

            with self.assertRaisesRegex(ValueError, "resolves outside checkpoint root"):
                phenogpt2_to_hpo.validate_model(model)


if __name__ == "__main__":
    unittest.main()
