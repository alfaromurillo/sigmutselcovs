"""CLI smoke tests (no network)."""

import pytest

from sigmutselcovs.cli import build_parser, main


def test_projects_lists_registry(capsys):
    assert main(["projects"]) == 0
    out = capsys.readouterr().out
    assert "COAD" in out and "BRCA" in out
    assert "repliseq" in out


def test_download_dry_run(tmp_path, capsys):
    assert (
        main(
            [
                "download",
                "COAD",
                "--data-dir",
                str(tmp_path),
                "--dry-run",
            ]
        )
        == 0
    )
    out = capsys.readouterr().out
    assert "would-download" in out


def test_download_subset_choice_enforced(tmp_path):
    with pytest.raises(SystemExit):
        build_parser().parse_args(
            [
                "download",
                "COAD",
                "--data-dir",
                str(tmp_path),
                "--which",
                "bogus",
            ]
        )


def test_download_which_accepts_encode_chromatin(tmp_path):
    """Regression: --which's choices list is separate from
    download.py's _WHICH tuple and had drifted out of sync when
    encode_chromatin was added there."""
    args = build_parser().parse_args(
        [
            "download",
            "ACC",
            "--data-dir",
            str(tmp_path),
            "--which",
            "encode_chromatin",
        ]
    )
    assert args.which == ["encode_chromatin"]


def test_build_include_exclude_accept_encode_chromatin(tmp_path):
    for flag in ("--include", "--exclude"):
        args = build_parser().parse_args(
            [
                "build",
                "ACC",
                "--data-dir",
                str(tmp_path),
                flag,
                "encode_chromatin",
            ]
        )
        assert (
            args.include == ["encode_chromatin"]
            if flag == "--include"
            else args.exclude == ["encode_chromatin"]
        )


def test_fetch_unpublished_raises(tmp_path):
    from sigmutselcovs.fetch import CovariateArtifactsUnavailable

    # An unregistered project code is never in the packaged
    # zenodo.json's `records`, so this raises before any network
    # call -- matching this file's own "no network" contract. Using
    # a real TCGA code here would depend on whether that project
    # happens to be published on Zenodo at test-run time, which this
    # file must not assume either way.
    with pytest.raises(CovariateArtifactsUnavailable):
        main(["fetch", "ZZZZ", "--dest", str(tmp_path)])


def test_command_required():
    with pytest.raises(SystemExit):
        build_parser().parse_args([])


def test_fetch_columns_and_blocks_imply_full(tmp_path, monkeypatch):
    import pandas as pd

    from sigmutselcovs import fetch

    calls = {}

    def fake(project, which, **kwargs):
        calls.update(project=project, which=which, **kwargs)
        return pd.DataFrame(
            {"atac_1": [0.5]},
            index=pd.Index(["ENSG1"], name="ensembl_gene_id"),
        )

    monkeypatch.setattr(fetch, "fetch_covariate_matrix", fake)
    out = tmp_path / "sel.csv"
    code = main(
        [
            "fetch",
            "COAD",
            "--blocks",
            "atac",
            "chromatin",
            "--columns",
            "mrt",
            "--out",
            str(out),
        ]
    )
    assert code == 0
    assert calls["which"] == "full"
    assert calls["blocks"] == ["atac", "chromatin"]
    assert calls["columns"] == ["mrt"]
    assert pd.read_csv(out, index_col=0).loc["ENSG1", "atac_1"] == 0.5


def test_pca_artifact_needs_generic_dir_unless_own_only(tmp_path):
    args = build_parser().parse_args(
        ["pca-artifact", "COAD", "--data-dir", str(tmp_path)]
    )
    assert args.generic_dir is None and not args.own_only
    with pytest.raises(SystemExit):
        main(["pca-artifact", "COAD", "--data-dir", str(tmp_path)])
    args = build_parser().parse_args(
        [
            "pca-artifact",
            "COAD",
            "--data-dir",
            str(tmp_path),
            "--generic-dir",
            str(tmp_path / "GENERIC"),
            "--variance-threshold",
            "0.95",
        ]
    )
    assert args.generic_dir.endswith("GENERIC")
    assert args.variance_threshold == 0.95
