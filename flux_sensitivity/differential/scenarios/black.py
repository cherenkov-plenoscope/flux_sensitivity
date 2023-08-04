import numpy as np


def init_matrix_B(probability_reco_given_true, containment=0.68):
    # ax0 -> true
    # ax1 -> reco
    num_bins = probability_reco_given_true.shape[0]
    M = probability_reco_given_true
    mask = np.zeros(shape=(num_bins, num_bins))

    # estimate containment regions:
    for etrue in range(num_bins):
        if np.sum(M[etrue, :]) > 0.0:
            assert 0.99 < np.sum(M[etrue, :]) < 1.01

            accumulated_containment = 0.0
            reco_best = np.argmax(M[etrue, :])

            accumulated_containment, weight = next_containment_and_weight(
                accumulated_containment=accumulated_containment,
                bin_containment=M[etrue, reco_best],
                target_containment=containment,
            )

            mask[etrue, reco_best] = weight
            start = reco_best - 1
            stop = reco_best + 1
            i = 0
            while accumulated_containment < containment:
                if start > 0:
                    accumulated_containment, w = next_containment_and_weight(
                        accumulated_containment=accumulated_containment,
                        bin_containment=M[etrue, start],
                        target_containment=containment,
                    )
                    mask[etrue, start] = w
                    start -= 1
                if accumulated_containment == containment:
                    break

                if stop + 1 < num_bins:
                    accumulated_containment, w = next_containment_and_weight(
                        accumulated_containment=accumulated_containment,
                        bin_containment=M[etrue, stop],
                        target_containment=containment,
                    )
                    mask[etrue, stop] = w
                    stop += 1
                if accumulated_containment == containment:
                    break

                if start == 0 and stop + 1 == num_bins:
                    break

                i += 1
                assert i < 2 * num_bins
    return mask


def next_containment_and_weight(
    accumulated_containment, bin_containment, target_containment,
):
    assert 0 <= accumulated_containment <= 1
    assert 0 <= bin_containment <= 1
    assert 0 < target_containment <= 1

    missing_containment = target_containment - accumulated_containment
    assert missing_containment > 0

    if bin_containment > 0:
        weight = np.min([missing_containment / bin_containment, 1])
    else:
        weight = 0

    if weight == 1:
        return accumulated_containment + bin_containment, 1
    else:
        return target_containment, weight
