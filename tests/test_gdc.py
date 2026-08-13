"""GDC client tests (mocked HTTP, golden rows from real COAD files)."""

import json

import pytest

from sigmutselcovs.gdc import (
    build_files_filter,
    query_gdc_files,
    write_gdc_manifest,
    write_gdc_sample_sheet,
)

# Two hits transcribed from the real TCGA-COAD artifacts
# (gdc_manifest.2025-08-13.134154.txt / gdc_sample_sheet.2025-08-18.tsv).
HIT_1 = {
    "file_id": "2c73d8ac-84b9-4094-9f5e-fbd6c2e8eabf",
    "file_name": (
        "c8c7600d-2c13-4df1-b058-2d02432bd49e"
        ".rna_seq.augmented_star_gene_counts.tsv"
    ),
    "md5sum": "4f9639a7ecab8b0dab8db9b17f145a2c",
    "file_size": 4216302,
    "state": "released",
    "data_category": "Transcriptome Profiling",
    "data_type": "Gene Expression Quantification",
    "cases": [
        {
            "project": {"project_id": "TCGA-COAD"},
            "submitter_id": "TCGA-AA-3688",
            "samples": [
                {
                    "submitter_id": "TCGA-AA-3688-01A",
                    "tissue_type": "Tumor",
                    "tumor_descriptor": "Primary",
                    "specimen_type": None,
                    "preservation_method": None,
                }
            ],
        }
    ],
}
HIT_2 = {
    "file_id": "7389dfaa-fa36-431f-a3b9-960b4fb9d836",
    "file_name": (
        "20ed86b0-2e3d-44ec-88cd-eac5486fe15c"
        ".rna_seq.augmented_star_gene_counts.tsv"
    ),
    "md5sum": "1f28a2c99a32dcef16c7779242c41d7d",
    "file_size": 4222773,
    "state": "released",
    "data_category": "Transcriptome Profiling",
    "data_type": "Gene Expression Quantification",
    "cases": [
        {
            "project": {"project_id": "TCGA-COAD"},
            "submitter_id": "TCGA-G4-6298",
            "samples": [
                {
                    "submitter_id": "TCGA-G4-6298-01A",
                    "tissue_type": "Tumor",
                    "tumor_descriptor": "Primary",
                    "specimen_type": "Solid Tissue",
                    "preservation_method": None,
                }
            ],
        }
    ],
}

MANIFEST_LINE_1 = (
    "2c73d8ac-84b9-4094-9f5e-fbd6c2e8eabf\t"
    "c8c7600d-2c13-4df1-b058-2d02432bd49e"
    ".rna_seq.augmented_star_gene_counts.tsv\t"
    "4f9639a7ecab8b0dab8db9b17f145a2c\t4216302\treleased"
)
SHEET_LINE_1 = (
    "2c73d8ac-84b9-4094-9f5e-fbd6c2e8eabf\t"
    "c8c7600d-2c13-4df1-b058-2d02432bd49e"
    ".rna_seq.augmented_star_gene_counts.tsv\t"
    "Transcriptome Profiling\t"
    "Gene Expression Quantification\tTCGA-COAD\t"
    "TCGA-AA-3688\tTCGA-AA-3688-01A\tTumor\tPrimary\t"
    "Unknown\tUnknown"
)
SHEET_LINE_2 = (
    "7389dfaa-fa36-431f-a3b9-960b4fb9d836\t"
    "20ed86b0-2e3d-44ec-88cd-eac5486fe15c"
    ".rna_seq.augmented_star_gene_counts.tsv\t"
    "Transcriptome Profiling\t"
    "Gene Expression Quantification\tTCGA-COAD\t"
    "TCGA-G4-6298\tTCGA-G4-6298-01A\tTumor\tPrimary\t"
    "Solid Tissue\tUnknown"
)


def test_build_files_filter_golden():
    got = build_files_filter(
        "TCGA-COAD",
        data_type="Gene Expression Quantification",
        workflow_type="STAR - Counts",
    )
    expected = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": ["TCGA-COAD"],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "data_type",
                    "value": ["Gene Expression Quantification"],
                },
            },
            {
                "op": "in",
                "content": {"field": "access", "value": ["open"]},
            },
            {
                "op": "in",
                "content": {
                    "field": "analysis.workflow_type",
                    "value": ["STAR - Counts"],
                },
            },
        ],
    }
    assert got == expected


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload


class _FakeSession:
    """Serves the two hits over two pages."""

    def __init__(self):
        self.calls = []

    def post(self, url, json=None, timeout=None):
        self.calls.append(json)
        start = json["from"]
        hits = [HIT_2, HIT_1]  # unsorted on purpose
        page = hits[start : start + 1]
        return _FakeResponse(
            {
                "data": {
                    "hits": page,
                    "pagination": {
                        "count": len(page),
                        "total": len(hits),
                    },
                }
            }
        )


def test_query_paginates_and_sorts():
    session = _FakeSession()
    hits = query_gdc_files("TCGA-COAD", session=session, page_size=1)
    assert len(session.calls) == 2
    assert [h["file_id"] for h in hits] == [
        HIT_1["file_id"],
        HIT_2["file_id"],
    ]  # sorted
    sent_filter = session.calls[0]["filters"]
    assert json.dumps(sent_filter) == json.dumps(
        build_files_filter(
            "TCGA-COAD",
            data_type="Gene Expression Quantification",
            workflow_type="STAR - Counts",
        )
    )


def test_manifest_matches_gdc_client_format(tmp_path):
    path = write_gdc_manifest(
        [HIT_1, HIT_2], tmp_path / "manifest.txt"
    )
    lines = path.read_text().splitlines()
    assert lines[0] == "id\tfilename\tmd5\tsize\tstate"
    assert lines[1] == MANIFEST_LINE_1


def test_sample_sheet_golden_rows(tmp_path):
    path = write_gdc_sample_sheet(
        [HIT_1, HIT_2], tmp_path / "sheet.tsv"
    )
    lines = path.read_text().splitlines()
    assert lines[0] == (
        "File ID\tFile Name\tData Category\tData Type\t"
        "Project ID\tCase ID\tSample ID\tTissue Type\t"
        "Tumor Descriptor\tSpecimen Type\t"
        "Preservation Method"
    )
    assert lines[1] == SHEET_LINE_1
    assert lines[2] == SHEET_LINE_2


def test_sample_sheet_dated_name_and_loader_roundtrip(tmp_path):
    path = write_gdc_sample_sheet([HIT_1], directory=tmp_path)
    assert path.name.startswith("gdc_sample_sheet.")
    assert path.name.endswith(".tsv")

    import pandas as pd

    sheet = pd.read_csv(path, sep="\t")
    # the columns import_tcga_gene_expression relies on
    assert {
        "File ID",
        "File Name",
        "Sample ID",
        "Tissue Type",
    } <= set(sheet.columns)
    assert sheet.loc[0, "Tissue Type"] == "Tumor"


def test_sample_sheet_requires_destination():
    with pytest.raises(ValueError, match="path or directory"):
        write_gdc_sample_sheet([HIT_1])
