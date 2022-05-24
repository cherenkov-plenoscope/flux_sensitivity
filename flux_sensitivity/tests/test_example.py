import flux_sensitivity
import pkg_resources
import glob
import os
import pandas


RESOURCE_DIR = pkg_resources.resource_filename(
    "flux_sensitivity",
    os.path.join(
        'tests',
        'resources',
    )
)


def read_resources(resource_dir):
    resources = {}
    for p in glob.glob(os.path.join(resource_dir, "*.csv")):
        fname = os.path.basename(p)
        fname = os.path.splitext(fname)[0]

        a = pandas.read_csv(p).to_numpy()
        if len(a.shape) == 2:
            if a.shape[1] == 1:
                a = a.reshape((a.shape[0],))
        resources[fname] = a
    return resources


def test_resources():
    RES = read_resources(resource_dir=RESOURCE_DIR)

    print(RES)