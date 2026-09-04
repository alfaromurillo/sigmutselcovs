"""Console entry point: the ``sigmutselcovs`` command.

Subcommands:

- ``projects`` — list registered cancer types
- ``download <PROJECT> --data-dir DIR [--which ...]`` — fetch source
  data
- ``build <PROJECT> --data-dir DIR`` — build and cache the matrices
- ``validate <PROJECT> --data-dir DIR`` — run the sanity and
  plausibility checks
- ``fetch <PROJECT> [--which pca]`` — download pre-built matrices
  from Zenodo
- ``check-updates`` — compare covariate sources to their recorded
  state
- ``download-gtex`` — fetch the GTEx GCT into the user cache
"""

import argparse
import logging
import sys


def _add_project_arg(parser):
    parser.add_argument(
        "project", help="TCGA study code (e.g. COAD, BRCA)"
    )


def _add_data_dir_arg(parser):
    parser.add_argument(
        "--data-dir", required=True, help="project data directory"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sigmutselcovs",
        description="Tissue-specific covariate matrices for "
        "sigmutsel",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="debug logging"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("projects", help="list registered cancer types")

    p = sub.add_parser("download", help="download source data")
    _add_project_arg(p)
    _add_data_dir_arg(p)
    p.add_argument(
        "--which",
        nargs="+",
        choices=[
            "gexp",
            "repliseq",
            "roadmap",
            "atac",
            "encode_chromatin",
        ],
        default=[
            "gexp",
            "repliseq",
            "roadmap",
            "atac",
            "encode_chromatin",
        ],
    )
    p.add_argument("--force", action="store_true")
    p.add_argument("--dry-run", action="store_true")

    p = sub.add_parser("build", help="build the covariate matrices")
    _add_project_arg(p)
    _add_data_dir_arg(p)
    p.add_argument(
        "--force-generation",
        action="store_true",
        help="recompute the per-gene caches",
    )
    p.add_argument(
        "--include",
        nargs="+",
        default=None,
        choices=[
            "gtex",
            "gexp",
            "repliseq",
            "roadmap",
            "atac",
            "encode_chromatin",
        ],
    )
    p.add_argument(
        "--exclude",
        nargs="+",
        default=(),
        choices=[
            "gtex",
            "gexp",
            "repliseq",
            "roadmap",
            "atac",
            "encode_chromatin",
        ],
    )
    p.add_argument(
        "--no-cache",
        action="store_true",
        help="do not write covariate_matrices/",
    )
    p.add_argument(
        "--strict",
        action="store_true",
        help="fail on validation errors",
    )

    p = sub.add_parser("validate", help="run covariate checks")
    _add_project_arg(p)
    _add_data_dir_arg(p)
    p.add_argument("--strict", action="store_true")

    p = sub.add_parser(
        "fetch", help="download pre-built matrices from Zenodo"
    )
    _add_project_arg(p)
    p.add_argument(
        "--which",
        default="pca",
        choices=["pca", "full", "simple", "tcga"],
    )
    p.add_argument("--dest", default=None)
    p.add_argument("--force", action="store_true")

    p = sub.add_parser(
        "check-updates", help="compare sources to recorded state"
    )
    p.add_argument(
        "--update-file",
        action="store_true",
        help="rewrite known states in sources.json",
    )
    p.add_argument("--sources", nargs="+", default=None)

    p = sub.add_parser(
        "download-gtex", help="fetch the GTEx GCT into the user cache"
    )
    p.add_argument("--force", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )

    if args.command == "projects":
        from .registry import load_registry

        for code, spec in sorted(load_registry().items()):
            sources = [
                s
                for s in (
                    "gtex",
                    "gexp",
                    "repliseq",
                    "roadmap",
                    "atac",
                )
                if getattr(spec, s) is not None
            ]
            print(
                f"{code}: {spec.description} "
                f"[{', '.join(sources)}]"
            )
        return 0

    if args.command == "download":
        from .download import download_covariates

        report = download_covariates(
            args.project,
            args.data_dir,
            which=tuple(args.which),
            force=args.force,
            dry_run=args.dry_run,
        )
        print(report.summary())
        failed = any(
            info["status"] == "failed"
            for info in report.sources.values()
        )
        return 1 if failed else 0

    if args.command == "build":
        from .builder import build_covariate_matrix

        matrices = build_covariate_matrix(
            args.project,
            args.data_dir,
            force_generation=args.force_generation,
            include=args.include,
            exclude=tuple(args.exclude),
            cache_matrices=not args.no_cache,
            strict_validation=args.strict,
        )
        for name in ("full", "simple", "tcga"):
            shape = getattr(matrices, name).shape
            print(
                f"cov_matrix_{name}: {shape[0]} genes x "
                f"{shape[1]} columns"
            )
        return 0

    if args.command == "validate":
        from .validate import (
            print_validation_report,
            validate_covariates,
        )

        frame = validate_covariates(
            args.project, args.data_dir, strict=False
        )
        print_validation_report(frame)
        failed = (frame["status"] == "fail").any()
        if args.strict and failed:
            return 1
        return 0

    if args.command == "fetch":
        from .fetch import fetch_covariate_matrix

        frame = fetch_covariate_matrix(
            args.project, args.which, dest=args.dest, force=args.force
        )
        print(
            f"cov_matrix_{args.which}: {frame.shape[0]} genes x "
            f"{frame.shape[1]} columns"
        )
        return 0

    if args.command == "check-updates":
        from .updates import check_updates, print_update_report

        frame = check_updates(
            sources=args.sources, update_file=args.update_file
        )
        print_update_report(frame)
        return 1 if (frame["status"] == "changed").any() else 0

    if args.command == "download-gtex":
        from .download import ensure_gtex_gct

        path = ensure_gtex_gct(force=args.force)
        print(path)
        return 0

    raise AssertionError(f"unhandled command {args.command}")


if __name__ == "__main__":
    sys.exit(main())
