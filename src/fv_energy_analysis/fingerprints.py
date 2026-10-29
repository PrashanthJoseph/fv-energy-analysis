"""Command-line entry point for saved interaction fingerprint reports."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .fingerprint_reports import generate


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-results", required=True)
    parser.add_argument("--cif-folder", type=Path, help="Original ModelCIF folder for H97–H113 pLDDT")
    parser.add_argument("--parent-cif", type=Path, help="Optional original parent ModelCIF for reference pLDDT")
    parser.add_argument("--scores-csv", type=Path, help="Merged Boltz-2/GNINA score CSV to append to saved reports")
    parser.add_argument("--export-xlsx", action="store_true",
                        help="Also export Excel (requires the reporting optional dependency)")
    parser.add_argument("--include-empty-metadata-sections", action="store_true",
                        help="Keep pLDDT/Boltz-2/GNINA sections blank when their original sources are unavailable")
    for analysis in ("discovery-studio", "e2-architecture", "hcdr3-plddt", "boltz2-scores", "gnina-scores"):
        parser.add_argument("--skip-" + analysis, action="store_true")
    args = parser.parse_args(argv)
    batch = Path(args.batch_results).resolve()
    summary = generate(batch, args.cif_folder, args.parent_cif, args.scores_csv,
                       include_empty_metadata_sections=args.include_empty_metadata_sections,
                       analysis_options={key: False for key in ("discovery_studio", "e2_architecture", "hcdr3_plddt", "boltz2_scores", "gnina_scores")
                                         if getattr(args, "skip_" + key)})
    print(f"Grouped fingerprint summary: {summary}")
    if args.export_xlsx:
        from .report import write_json
        from .spreadsheet import export_workbook

        payload = json.loads((batch / "batch_summary_grouped_fingerprints.json").read_text(encoding="utf-8"))
        output = export_workbook(payload, batch / "batch_summary_grouped_fingerprints.xlsx")
        manifest_path = batch / "batch_manifest.json"
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest.setdefault("interaction_fingerprints", {})["grouped_summary_xlsx"] = str(output)
            write_json(manifest_path, manifest)
        print(f"Grouped Excel summary: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
