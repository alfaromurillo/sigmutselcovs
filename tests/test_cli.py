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


def test_fetch_unpublished_raises(tmp_path):
    from sigmutselcovs.fetch import CovariateArtifactsUnavailable

    with pytest.raises(CovariateArtifactsUnavailable):
        main(["fetch", "BRCA", "--dest", str(tmp_path)])


def test_command_required():
    with pytest.raises(SystemExit):
        build_parser().parse_args([])
