"""Population statistics used by the MD time-series containers."""

import numpy as np

gaussian = (
    lambda x, mu, sigma: 1.0
    / (2.0 * np.pi * sigma**2.0) ** 0.5
    * np.exp(-((x - mu) ** 2.0) / (2.0 * sigma**2))
)


def get_mean_std(values):
    """Return the arithmetic mean and population standard deviation."""

    mean = np.average(values)
    variance = np.average((values - mean) ** 2)
    return mean, np.sqrt(variance)


def get_mean_cum(values):
    """Return the cumulative arithmetic mean."""

    cum = np.cumsum(values)
    return cum / np.array(list(range(1, len(values) + 1)))


def get_mean_std_cum(values):
    """Return cumulative population means and standard deviations."""

    values_array = np.array(values)
    mean_cum = get_mean_cum(values_array)
    std_cum = np.sqrt(get_mean_cum(values_array**2.0) - mean_cum**2)
    return mean_cum, std_cum


def test_mean_std_cum(N):
    """Display an interactive demonstration for ``N`` random samples."""

    from matplotlib import pyplot as plt

    values = np.random.rand(N)
    mean_cum, std_cum = get_mean_std_cum(values)
    mean, std = get_mean_std(values)

    plt.plot(list(range(1, N + 1)), mean_cum, c="b")
    plt.axhline(y=mean, c="r", ls="--")
    plt.grid()
    plt.show()

    plt.plot(list(range(1, N + 1)), std_cum)
    plt.axhline(y=std, ls="--")
    plt.grid()
    plt.show()
