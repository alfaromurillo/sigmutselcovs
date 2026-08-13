# Contributing to sigmutselcovs

Thank you for considering contributing to sigmutselcovs!

## Development Setup

1. Clone the repository:
```bash
git clone https://github.com/alfaromurillo/sigmutselcovs.git
cd sigmutselcovs
```

2. Create a virtual environment:
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

3. Install sigmutsel (not yet on PyPI) and this package in
   development mode:
```bash
pip install "sigmutsel @ git+https://github.com/alfaromurillo/sigmutsel.git"
pip install -e ".[dev]"
```

See `DEVELOPMENT.md` for the module map, the generalized
download/build/validate workflow, and how to register a new cancer
type.

## Running Tests

```bash
pytest
```

Network-marked tests are deselected by default
(`-m 'not network'`, see `pyproject.toml`); slow tests are tagged
`slow`. Most of the suite runs against synthetic fixtures with no
external calls.

## Code Style

This project uses:
- **black** for code formatting (70 character line length, matching
  sigmutsel)
- **ruff** for linting

Format your code before committing:
```bash
black src/ tests/
ruff check src/ tests/ --fix
```

A `pre-commit` hook enforces both automatically (`black --check` +
`ruff check`, blocking the commit on failure). Install it once per
clone:
```bash
pip install pre-commit
pre-commit install
```

## Pull Request Process

1. Create a new branch for your feature/fix
2. Make your changes
3. Add tests for new functionality
4. Ensure all tests pass
5. Format your code
6. Update `CHANGELOG.md` under `[Unreleased]`
7. Submit a pull request

## Reporting Issues

Please use the GitHub issue tracker to report bugs or suggest
features.
