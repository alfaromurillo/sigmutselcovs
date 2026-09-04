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
