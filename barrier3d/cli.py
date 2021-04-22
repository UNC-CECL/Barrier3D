#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import pathlib
import re
import sys
from collections import OrderedDict
from functools import partial

import click
import numpy as np
import pandas as pd
import pkg_resources
import yaml

from .barrier3d import Barrier3d, Barrier3dError
from .bmi import Barrier3dBmi
from .configuration import Barrier3dConfiguration

__version__ = "0.1"


out = partial(click.secho, bold=True, err=True)
err = partial(click.secho, fg="red", err=True)


def csv_contents(infile):
    datadir = pathlib.Path(pkg_resources.resource_filename("barrier3d", "data"))
    with open(datadir / f"barrier3d-{infile}.csv") as fp:
        return fp.read()


INFILES = {
    "parameters": Barrier3dConfiguration().to_yaml(),
    "elevations": csv_contents("elevations"),
    "dunes": csv_contents("dunes"),
    "growthparam": csv_contents("growthparam"),
    "storms": csv_contents("storms"),
}


class Barrier3dOutputWriter:
    def __init__(self, bmi, output_interval=1, filepath="output.csv"):
        self._bmi = bmi
        self._output_interval = output_interval
        self._steps = 0

        with open(filepath, mode="w") as fp:
            print("# version: {0}".format(__version__), file=fp)
            print(
                ",".join(
                    [
                        "# Time step [-]",
                        "Overwash flux [m^3/m]",
                        "Shoreface flux [m^3/m]",
                        "Change in shoreface toe position [dam]",
                        "Change in shoreline position [dam]",
                        "Back-barrier shoreline position [dam]",
                    ]
                ),
                file=fp,
            )

    def _save_output(self):
        data = pd.DataFrame(
            {
                "time_stamp": self._steps,
                "q_ow": self._bmi.QowTS[-1],
                "q_sf": self._bmi.QsfTS[-1],
                "dx_toe_dt": np.diff(self._bmi.x_t_TS[-2:]),
                "dx_sf_dt": np.diff(self._bmi.x_s_TS[-2:]),
                "x_b": self._bmi.x_b_TS[-1],
            },
            index=(0,),
        )
        data.to_csv("output.csv", mode="a", index=False, sep=",", header=False)

    def update(self, n_steps):
        self._steps += n_steps
        if self._steps % self._output_interval == 0:
            self._save_output()


@click.group(chain=True)
@click.version_option()
@click.option(
    "--cd",
    default=".",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="chage to directory, then execute",
)
def barrier3d(cd) -> None:
    """Barrier3D: A spatially explicit exploratory model of barrier island evolution in three dimensions

    \b
    Examples:

      Create a folder with example input files,

        $ mkdir example && b3d setup

      Run a simulation using the examples input files,

        $ cd example && b3d run

      Commands can also be chained together. The following will setup
      a simulation, run it, and then plot elevations.

        $ mkdir example && b3d setup run plot elevation
    """
    os.chdir(cd)


@barrier3d.command()
@click.option("-v", "--verbose", is_flag=True, help="Emit status messages to stderr.")
@click.option("--dry-run", is_flag=True, help="Do not actually run the model")
def run(dry_run: bool, verbose: bool) -> None:
    """Run a simulation."""
    run_dir = pathlib.Path.cwd()
    # config_file = run_dir / "barrier3d.yaml"

    message = []
    # if not config_file.is_file():
    #     message.append("missing Barrier3D configuration file: {0}".format(config_file))
    if (run_dir / "output.csv").exists():
        message.append(
            "Barrier3D output file already exists: {0}".format(run_dir / "output.csv")
        )
    if message:
        err(os.linesep.join(message))
        raise click.Abort(os.linesep.join(message))

    barrier3d = Barrier3d.from_path(run_dir, fmt="yaml")
    output = Barrier3dOutputWriter(barrier3d, filepath=run_dir / "output.csv")

    n_steps = barrier3d._TMAX - 1

    if dry_run:
        out(str(run_dir))
    elif n_steps == 0:
        out("Nothing to do (years == 0). 😴")
    else:
        with click.progressbar(
            range(n_steps),
            label=" ".join(["🚀", str(run_dir)]),
            item_show_func=lambda step: "step {0} of {1}".format(
                int(0 if step is None else step), n_steps
            ),
        ) as bar:
            for step in bar:
                barrier3d.update()
                output.update(1)

        out("💥 Finished! 💥")


@barrier3d.command()
@click.argument("infile", type=click.Choice(sorted(INFILES)))
def show(infile: str) -> None:
    """Show example input files."""
    print(_contents_of_input_file(infile))


@barrier3d.command()
def setup() -> None:
    """Setup a folder of input files for a simulation."""
    files = {
        "parameters": pathlib.Path("barrier3d-parameters.yaml"),
        "elevations": pathlib.Path("barrier3d-elevations.csv"),
        "dunes": pathlib.Path("barrier3d-dunes.csv"),
        "growthparam": pathlib.Path("barrier3d-growthparam.csv"),
        "storms": pathlib.Path("barrier3d-storms.csv"),
    }

    existing_files = [name for name in files.values() if name.exists()]
    if existing_files:
        for name in existing_files:
            err(
                f"{name}: File exists. Either remove and then rerun or setup in a different folder",
            )
    else:
        for infile, fname in files.items():
            with open(fname, "w") as fp:
                print(_contents_of_input_file(infile), file=fp)

    if existing_files:
        raise click.Abort()


def _contents_of_input_file(infile: str) -> str:
    if infile not in INFILES:
        raise ValueError(
            "unknown input file type ({infile} not one of {infiles})".format(
                infile=infile, infiles=", ".join(sorted(infiles))
            )
        )

    return INFILES[infile]


@barrier3d.command()
@click.argument(
    "value",
    type=click.Choice(["time_step", "q_ow", "q_sf", "dx_toe_dt", "dx_sf_dt", "x_bb"]),
)
def plot(value: str) -> None:
    """Plot output from a simulation."""
    import matplotlib.pyplot as plt

    filepath = pathlib.Path("output.csv")
    out(str(filepath.resolve()))

    data = pd.read_csv(
        filepath,
        names=("time_step", "q_ow", "q_sf", "dx_toe_dt", "dx_sf_dt", "x_bb"),
        comment="#",
    )
    plt.plot(data["time_step"], data[value])

    plt.show()
