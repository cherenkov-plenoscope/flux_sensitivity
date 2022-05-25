import flux_sensitivity
import numpy as np


def test_bell_spectrum_mask_shape():

    for n in range(3, 19):
        Prgt = np.eye(n)

        mask = flux_sensitivity.differential.scenarios.bell_spectrum.init_B(
            probability_reco_given_true=Prgt, containment=0.68,
        )

        assert mask.shape == Prgt.shape
