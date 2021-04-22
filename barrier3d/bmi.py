#! /usr/bin/env python
"""Basic Model Interface implementation for River Module"""
import pathlib

import numpy as np
from bmipy import Bmi

from .barrier3d import Barrier3d


class Barrier3dBmi(Bmi):

    _name = "Barrier3D"
    _input_var_names = ()
    _output_var_names = (
        "overwash_flux",
        "shoreface_flux",
        "shoreface_toe_position",
        "shoreline_position",
        "back_barrier_shoreline_position",
        "height_of_barrier",
    )

    def __init__(self):
        """Create a Barrier3dBmi module that is ready for initialization."""
        self._model = None
        self._values = {}
        self._var_units = {}

    def initialize(self, config_file):
        filepath = pathlib.Path(config_file)

        if filepath.name != "barrier3d-parameters.yaml":
            raise ValueError(
                "barrier3d parameter file must be named barrier3d-parameters.yaml"
            )

        self._model = Barrier3d.from_path(filepath.parent, fmt="yaml")

        self._values = {
            "overwash_flux": lambda: np.array(self._model.QowTS[:]),
            "shoreface_flux": lambda: np.array(self._model.QsfTS[:]),
            "shoreface_toe_position": lambda: np.array(self._model.x_t_TS[:]),
            "shoreline_position": lambda: np.array(self._model.x_s_TS[:]),
            "back_barrier_shoreline_position": lambda: np.array(self._model.x_b_TS[:]),
            "height_of_barrier": lambda: np.array(self._model.h_b_TS[:]),
        }

        self._var_units = {
            "overwash_flux": "m^3 / m",
            "shoreface_flux": "m^3 / m",
            "shoreface_toe_position": "dam",
            "shoreline_position": "dam",
            "back_barrier_shoreline_position": "dam",
            "height_of_barrier": "dam",
        }

        self._var_type = {}
        for name in self._input_var_names + self._output_var_names:
            self._var_type[name] = str(np.dtype(float))

    def update(self):
        self._model.update()
        self._model.update_dune_domain()

    def update_frac(self, time_frac):
        """Update model by a fraction of a time step."""
        raise NotImplementedError("update_frac")

    def update_until(self, then):
        """Update model until a particular time."""
        raise NotImplementedError("update_until")

    def finalize(self):
        pass

    def get_var_type(self, name):
        return self._var_type[name]

    def get_var_units(self, name):
        return self._var_units[name]

    def get_var_nbytes(self, name):
        return self._values[name]().nbytes

    def get_var_itemsize(self, name):
        return np.dtype(self.get_var_type(name)).itemsize

    def get_var_size(self, name):
        return self._values[name]().size

    def get_var_location(self, name):
        return "none"

    def get_var_grid(self, name):
        raise NotImplementedError("get_var_grid")

    def get_grid_rank(self, grid):
        raise NotImplementedError("get_grid_rank")

    def get_grid_size(self, grid):
        raise NotImplementedError("get_grid_size")

    def get_grid_shape(self, grid, shape):
        raise NotImplementedError("get_grid_shape")

    def get_grid_spacing(self, grid, spacing):
        raise NotImplementedError("get_grid_spacing")

    def get_grid_origin(self, grid, origin):
        raise NotImplementedError("get_grid_origin")

    def get_grid_type(self, grid):
        raise NotImplementedError("get_grid_type")

    def get_value(self, name, dest):
        dest[:] = self._values[name]()
        return dest

    def set_value(self, var_name, new_vals):
        """Set model values."""
        raise NotImplementedError("set_value")

    def get_component_name(self):
        return self._name

    def get_input_var_names(self):
        return self._input_var_names

    def get_output_var_names(self):
        return self._output_var_names

    def get_input_item_count(self):
        return len(self._input_var_names)

    def get_output_item_count(self):
        return len(self._output_var_names)

    def get_start_time(self):
        return 0.0

    def get_end_time(self):
        return float(self._model._TMAX)

    def get_current_time(self):
        return float(self._model.time_step)

    def get_time_step(self):
        return 1.0

    def get_time_units(self):
        return "yr"

    def get_grid_node_count(self, grid):
        raise NotImplementedError("get_grid_node_count")

    def get_grid_edge_count(self, grid):
        raise NotImplementedError("get_grid_edge_count")

    def get_grid_face_count(self, grid):
        raise NotImplementedError("get_grid_face_count")

    def get_grid_edge_nodes(self, grid, edge_nodes):
        raise NotImplementedError("get_grid_edge_nodes")

    def get_grid_face_edges(self, grid, face_edges):
        raise NotImplementedError("get_grid_face_edges")

    def get_grid_face_nodes(self, grid, face_nodes):
        raise NotImplementedError("get_grid_edge_nodes")

    def get_grid_nodes_per_face(self, grid, nodes_per_face):
        raise NotImplementedError("get_grid_nodes_per_face")

    def get_grid_x(self, grid, x):
        raise NotImplementedError("get_grid_x")

    def get_grid_y(self, grid, y):
        raise NotImplementedError("get_grid_y")

    def get_grid_z(self, grid, z):
        raise NotImplementedError("get_grid_z")

    def get_value_at_indices(self, name, dest, inds):
        raise NotImplementedError("get_value_at_indices")

    def get_value_ptr(self, name):
        raise NotImplementedError("get_value_ptr")

    def set_value_at_indices(self, name, ids, src):
        raise NotImplementedError("set_value_at_indices")
