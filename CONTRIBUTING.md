# Contributing to PyGaia

Contributions are welcome: bug reports, documentation fixes, new features, tests, or ideas.

## Reporting issues

Before opening an issue, check the [issue tracker](https://github.com/agabrown/PyGaia/issues) to see whether it has already been reported. Use the bug report or feature request template, and include a minimal example that reproduces the problem.

## Development setup

```bash
git clone https://github.com/agabrown/PyGaia.git
cd PyGaia
python -m pip install -e ".[test]"
```

The version number is derived from the git tags by [hatch-vcs](https://github.com/ofek/hatch-vcs), so make sure your clone has the tags (`git fetch --tags`).

## Tests

Tests live in `src/pygaia/tests` and are run with [pytest](https://docs.pytest.org/):

```bash
python -m pytest
```

Please add tests for any new functionality or bug fix.

## Pull requests

1. Fork the repository and create a branch for your change.
2. Make sure `python -m pytest` passes.
3. Add an entry to [HISTORY.md](HISTORY.md) describing the change.
4. Open a pull request and fill in the template.

## Releases

Releases are made by the maintainer by creating a GitHub release with a tag of the form `vX.Y.Z`. The tag sets the package version, and the publishing workflow uploads the release to [PyPI](https://pypi.org/project/PyGaia/). Remember to update `version` and `date-released` in [CITATION.cff](CITATION.cff) before tagging.
