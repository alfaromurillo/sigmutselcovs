"""sigmutselcovs — covariate matrix builders for sigmutsel."""

from .covariates_checks import (
    check_all,
    check_collinearity,
    check_missingness,
    check_skewness,
    check_variance,
    fix_all,
    fix_skewness,
    fix_variance,
)

try:
    from importlib.metadata import version

    __version__ = version("sigmutselcovs")
except Exception:  # noqa: BLE001 - version is informational only
    __version__ = "unknown"

# Lazy exports (PEP 562): keep `import sigmutselcovs` light — the
# builder pulls in pyBigWig and (via sigmutsel defaults) pymc.
_LAZY = {
    "build_covariate_matrix": "builder",
    "CovariateMatrices": "builder",
    "load_registry": "registry",
    "available_projects": "registry",
    "get_project": "registry",
    "ProjectSpec": "registry",
    "project_paths": "paths",
    "ProjectPaths": "paths",
    "bigwig_files": "paths",
    "download_covariates": "download",
    "DownloadReport": "download",
    "validate_covariates": "validate",
    "print_validation_report": "validate",
    "fetch_covariate_matrix": "fetch",
    "fetch_covariate_artifacts": "fetch",
    "CovariateArtifactsUnavailable": "fetch",
    "check_updates": "updates",
    "print_update_report": "updates",
}

__all__ = [
    "check_all",
    "check_collinearity",
    "check_missingness",
    "check_skewness",
    "check_variance",
    "fix_all",
    "fix_skewness",
    "fix_variance",
    *sorted(_LAZY),
]


def __getattr__(name):
    if name in _LAZY:
        import importlib

        module = importlib.import_module(f".{_LAZY[name]}", __name__)
        return getattr(module, name)
    raise AttributeError(
        f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(_LAZY))
