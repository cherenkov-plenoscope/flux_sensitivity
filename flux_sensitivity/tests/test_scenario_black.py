import flux_sensitivity
import numpy as np


def test_scenario_black_mask_shape():

    for n in range(3, 19):
        Prgt = np.eye(n)

        mask = flux_sensitivity.differential.scenarios.black.init_matrix_B(
            probability_reco_given_true=Prgt, containment=0.68,
        )

        assert mask.shape == Prgt.shape
