import flux_sensitivity
import numpy as np


def test_bell_spectrum_mask_shape():

    for n in range(3, 19):
        Prgt = np.eye(n)

        mask = flux_sensitivity.differential_sensitivity.make_mask_for_energy_confusion_matrix_for_bell_spectrum(
            probability_reco_given_true=Prgt, containment=0.68,
        )

        assert mask.shape == Prgt.shape
