import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from prepare_inputs import resolve_gene


class MappingTests(unittest.TestCase):
    def test_unique_old_mapping(self):
        self.assertEqual(resolve_gene("oldA", {"oldA": {"Glyma.01G000100"}}),
                         ("Glyma.01G000100", "mapped_v1_to_v2"))

    def test_ambiguity_is_not_first_match(self):
        self.assertEqual(resolve_gene("oldA", {"oldA": {"Glyma.01G000100", "Glyma.01G000200"}}),
                         ("", "excluded_ambiguous_mapping"))

    def test_blank_unknown_and_unplaced(self):
        self.assertEqual(resolve_gene(" ", {}), ("", "excluded_blank_id"))
        self.assertEqual(resolve_gene("oldA", {}), ("", "excluded_unmapped_id"))
        self.assertEqual(resolve_gene("Glyma.u000100", {}), ("Glyma.u000100", "v2_id"))

    def test_capitalization_is_not_version_mapping(self):
        self.assertEqual(resolve_gene("Glyma.01g000100", {}), ("Glyma.01G000100", "v2_id"))
        self.assertEqual(resolve_gene("Glyma01g00100", {}), ("", "excluded_unmapped_id"))


if __name__ == "__main__":
    unittest.main()
