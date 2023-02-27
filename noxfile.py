import difflib
import os
import pathlib
import shutil
import tempfile

import nox

PROJECT = "barrier3d"
ROOT = pathlib.Path(__file__).parent


nox.options.sessions = ["lint", "test", "test-bmi", "test-cli", "test-notebooks"]


@nox.session
def test(session: nox.Session) -> None:
    """Run the tests."""
    session.install("-r", "requirements-testing.txt")
    session.install(".")

    args = [
        "-n",
        "auto",
        "--cov",
        PROJECT,
        "-vvv",
    ] + session.posargs

    if "CI" in os.environ:
        args.append(f"--cov-report=xml:{ROOT.absolute()!s}/coverage.xml")
    session.run("pytest", *args)

    if "CI" not in os.environ:
        session.run("coverage", "report", "--ignore-errors", "--show-missing")


@nox.session(name="test-cli")
def test_cli(session: nox.Session):
    session.install(".")

    session.run("b3d", "--help")
    session.run("b3d", "--version")
    session.run("b3d", "plot", "--help")
    session.run("b3d", "run", "--help")
    session.run("b3d", "setup", "--help")
    session.run("b3d", "show", "--help")

    for infile in ["dunes", "elevations", "growthparam", "parameters", "storms"]:
        out = session.run("b3d", "show", infile, external=True, silent=True)
        if not out.strip():
            session.error(f"The command `b3d show {infile}` produced no output")

    expected = {
        "barrier3d-default-dunes.csv",
        "barrier3d-default-elevations.csv",
        "barrier3d-default-growthparam.csv",
        "barrier3d-default-parameters.yaml",
        "barrier3d-default-storms.csv",
    }

    with tempfile.TemporaryDirectory(ignore_cleanup_errors=False) as tmpdirname:
        with session.chdir(tmpdirname):
            session.run("b3d", "setup")
            actual = {str(p) for p in pathlib.Path(".").iterdir()}

    if actual != expected:
        session.log(
            os.linesep.join(
                ["diff actual expected"]
                + list(
                    difflib.unified_diff(
                        sorted(actual),
                        sorted(expected),
                        fromfile="actual",
                        tofile="expected",
                    )
                )
            )
        )
        session.error("The command `b3d setup` generated an unexpected set of files")


@nox.session(name="test-bmi", venv_backend="mamba")
def test_bmi(session: nox.Session) -> None:
    """Run the tests."""
    session.conda_install("bmi-tester")
    session.conda_install("--file", "requirements.txt")
    session.install(".", "--no-deps")

    session.run(
        "bmi-test",
        "--config-file=tests/test_bmi/barrier3d.yaml",
        "--root-dir=tests/test_bmi",
        "-vvv",
        "barrier3d.bmi:Barrier3dBmi",
    )


@nox.session(name="test-notebooks")
def test_notebooks(session: nox.Session) -> None:
    """Run the notebooks."""
    args = [
        "pytest",
        "notebooks",
        "--nbmake",
        "--nbmake-kernel=python3",
        "--nbmake-timeout=3000",
        "-n",
        "auto",
        "-vvv",
    ] + session.posargs

    session.install("-r", "requirements-testing.txt")
    session.install("nbmake")
    session.install("-r", "notebooks/requirements.txt")
    session.install(".")

    session.run(*args)


@nox.session
def lint(session: nox.Session) -> None:
    """Look for lint."""
    session.install("pre-commit")
    session.run("pre-commit", "run", "--all-files")


@nox.session
def build(session: nox.Session) -> None:
    """Build sdist and wheel dists."""
    session.install("pip")
    session.install("build")
    session.run("python", "--version")
    session.run("pip", "--version")
    session.run("python", "-m", "build", "--outdir", "./build/wheelhouse")


@nox.session
def release(session):
    """Tag, build and publish a new release to PyPI."""
    session.install("zest.releaser[recommended]")
    session.run("fullrelease")


@nox.session(name="publish-testpypi")
def publish_testpypi(session):
    """Publish wheelhouse/* to TestPyPI."""
    session.run("twine", "check", "build/wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "--repository-url",
        "https://test.pypi.org/legacy/",
        "build/wheelhouse/*.tar.gz",
    )


@nox.session(name="publish-pypi")
def publish_pypi(session):
    """Publish wheelhouse/* to PyPI."""
    session.run("twine", "check", "build/wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "build/wheelhouse/*.tar.gz",
    )


@nox.session(python=False)
def clean(session):
    """Remove all .venv's, build files and caches in the directory."""
    for folder in _args_to_folders(session.posargs):
        with session.chdir(folder):
            shutil.rmtree("build", ignore_errors=True)
            shutil.rmtree("build/wheelhouse", ignore_errors=True)
            shutil.rmtree(f"{PROJECT}.egg-info", ignore_errors=True)
            shutil.rmtree(".pytest_cache", ignore_errors=True)
            shutil.rmtree(".venv", ignore_errors=True)

            for pattern in ["*.py[co]", "__pycache__"]:
                _clean_rglob(pattern)


@nox.session(python=False, name="clean-checkpoints")
def clean_checkpoints(session):
    """Remove jupyter notebook checkpoint files."""
    for folder in _args_to_folders(session.posargs):
        with session.chdir(folder):
            _clean_rglob("*-checkpoint.ipynb")
            _clean_rglob(".ipynb_checkpoints")


def _args_to_folders(args):
    return [ROOT] if not args else [pathlib.Path(f) for f in args]


def _clean_rglob(pattern):
    nox_dir = pathlib.Path(".nox")

    for p in pathlib.Path(".").rglob(pattern):
        if nox_dir in p.parents:
            continue
        if p.is_dir():
            p.rmdir()
        else:
            p.unlink()
