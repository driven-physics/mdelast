"""Time-series containers and HDF5 serialization for MD observables."""

import numpy as np

from averages import *

import logging

if len(logging.root.handlers) == 0:
    logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Series")


class ThermalResults:
    """Namespace populated with post-thermalization aggregate properties."""

    pass


class TimeSeries:
    """Scalar samples with full-run and cumulative statistics."""

    def __init__(self):
        self.data = []
        self.mean = None
        self.std = None
        self.std_mean = None
        self.mean_cum = None
        self.std_cum = None

        self.data_gauss = None
        self.gauss = None

    def clean(self):
        self.data = []
        self.mean_cum = []
        self.std_cum = []

    def append(self, x):
        self.data.append(x)

    def set_averages(self, cycles_thermalization=0, cycles_thermalization_cum=0):
        """Compute statistics after the requested sample offsets."""

        self.cycles_thermalization = cycles_thermalization
        data_cut = self.data[cycles_thermalization:]
        self.mean, self.std = get_mean_std(data_cut)
        self.std_mean = self.std / abs(self.mean)
        self.mean_cum, self.std_cum = get_mean_std_cum(self.data[cycles_thermalization_cum:])

        self.data_gauss = np.linspace(self.mean - 4 * self.std, self.mean + 4 * self.std, 1000)
        self.gauss = gaussian(self.data_gauss, self.mean, self.std)

    def write_h5(self, file_h5, group_name):
        """Write available data and statistics below ``group_name``."""

        logger.info(f"writing TimeSeries to h5: {group_name!r}...")
        for key in [
            "data",
            "mean",
            "std",
            "std_mean",
            "mean_cum",
            "std_cum",
            "data_gauss",
            "gauss",
            "cycles_thermalization",
        ]:
            if hasattr(self, key):
                value = getattr(self, key)
                if not (value is None):
                    logger.info(f"writing {group_name!r}: {key!r}")
                    file_h5[f"{group_name}/{key}"] = value
        logger.info(f"Finished writing TimeSeries to h5: {group_name!r}")


class TimeSeriesVector:
    """Fixed-length vector whose components are independent ``TimeSeries``."""

    def __init__(self, N):
        self.data = [TimeSeries() for i in range(N)]
        self.N = N

    def clean(self):
        for i in range(self.N):
            self.data[i].clean()

    def append(self, vec):
        assert len(vec) == self.N, len(vec)
        for i, v in enumerate(vec):
            self.data[i].append(v)

    def get_time_array(self):
        return np.copy(np.array([self.data[i].data for i in range(self.N)]))

    def set_averages(self, cycles_thermalization=0, cycles_thermalization_cum=0):
        for i in range(self.N):
            self.data[i].set_averages(cycles_thermalization, cycles_thermalization_cum)

    def write_h5(self, file_h5, group_name):
        logger.info(f"writing TimeSeriesVector to h5: {group_name!r}")
        for i in range(self.N):
            self.data[i].write_h5(file_h5=file_h5, group_name=f"{group_name}_{i+1}")


class TimeSeriesMatrix:
    """Square matrix whose components are independent ``TimeSeries``."""

    def __init__(self, N):
        self.data = [[TimeSeries() for i in range(N)] for i in range(N)]
        self.N = N

    def clean(self):
        for i in range(self.N):
            for j in range(self.N):
                self.data[i][j].clean()

    def append(self, mat):
        assert mat.shape == (self.N, self.N), mat.shape
        for i in range(self.N):
            for j in range(self.N):
                self.data[i][j].append(mat[i, j])

    def get_time_array(self):
        N00 = len(self.data[0][0].data)
        assert all(len(self.data[i][j].data) == N00 for j in range(self.N) for i in range(self.N))
        return np.copy(
            np.array(
                [[np.array(self.data[i][j].data) for j in range(self.N)] for i in range(self.N)]
            )
        )

    def set_averages(self, cycles_thermalization=0, cycles_thermalization_cum=0):
        for i in range(self.N):
            for j in range(self.N):
                self.data[i][j].set_averages(cycles_thermalization, cycles_thermalization_cum)

    def write_h5(self, file_h5, group_name):
        logger.info(f"writing TimeSeriesMatrix to h5: {group_name!r}")
        for i in range(self.N):
            for j in range(self.N):
                self.data[i][j].write_h5(file_h5=file_h5, group_name=f"{group_name}_{i+1}_{j+1}")


class Series:
    """Collection of thermodynamic, stress, energy, and elasticity series."""

    def __init__(self):
        self.time = []
        self.temperature_kin = TimeSeries()
        self.temperature_conf_inverse = TimeSeries()
        self.total_laplacian = TimeSeries()
        self.total_force2 = TimeSeries()
        self.pressure = TimeSeries()
        self.stress_tensor = TimeSeriesVector(6)
        self.etot_atom = TimeSeries()
        self.epot_atom = TimeSeries()
        self.epot_eam_atom = TimeSeries()
        self.ekin_atom = TimeSeries()
        self.etot_cons = TimeSeries()
        self.etot_cons_err = TimeSeries()
        self.born_stress_tensor = TimeSeriesVector(6)
        self.born_elasticity_tensor = TimeSeriesMatrix(6)
        self.V_V0 = TimeSeries()
        self.L_L0 = TimeSeries()

        self.thermal_results = ThermalResults()

    def clean(self):
        logger.info("cleaning series...")
        for key, value in self.__dict__.items():
            try:
                value.clean()
                logger.info(f"cleaned '{key}' data")
            except:
                pass
