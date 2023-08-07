###########################################
Estimate the Sensitivity of your Instrument
###########################################
|TestStatus| |PyPiStatus| |BlackStyle| 

Made for astronomy with gamma-rays in the atmospheric Cherenkov-method to compute the integral or differential sensitivity of your instrument with respect to the gamma-rays energy.


************************
Differential Sensitivity
************************
For demonstration, we show how to compute the differential sensitivity for the Cherenkov-Telescope-Array (CTA) based on a public estimate for its response-function [``Prod5-South-20deg-AverageAz-14MSTs37SSTs.1800s-v0.1.fits``].
From the response-function we need the following three estimates:

Instrument response function
============================
To be precise: This response-function is what I (as a non CTA-member) extracted from the provided documentation. To extract these quantities I interpolate and/or average over e.g. solid angle.

Rate of background
    |img_irf_background_rate_onregion|
    
    The rate of background (the sum of all contributions) w.r.t. to reconstructed gamma-ray-energy.

Area for signal
    |img_irf_signal_area|

    The effective area to collect gamma-rays w.r.t. the true gamma-ray-energy.

Probability to confuse energy
    |img_irf_probability_reco_given_true|

    The probabilty to confuse the true with the reconstructed gamma-ray-energy.
    This is also called `energy migration'.

Scenarios
=========

This package offeres multiple scenarios how to handle the instrument's non perfec reconstruction in energy.
The scenarios are named using colors to stress their differences while not implying a hierachy.

Blue
----
The ``blue`` scenario handles the instrument's non perfec reconstruction in energy by simply ignoring it. This scenario assumes that the reconstructed gamma-ray-energy is equal to the true gamma-ray-energy.

+-------------------------------------+--------------------------------------------+
| Matrix ``G``                        | Matrix ``B``                               |
+=====================================+============================================+
| |img_G_matrix_in_scenario_blue|     | |img_B_matrix_in_scenario_blue|            |
+-------------------------------------+--------------------------------------------+

+-------------------------------------+--------------------------------------------+
| Area for signal                     | Rate of background                         |
+=====================================+============================================+
| |img_signal_area_in_scenario_blue|  | |img_background_rate_in_scenario_blue|     |
+-------------------------------------+--------------------------------------------+

Yellow
------
The ``yellow`` scenario...

+-------------------------------------+--------------------------------------------+
| Matrix ``G``                        | Matrix ``B``                               |
+=====================================+============================================+
| |img_G_matrix_in_scenario_yellow|   | |img_B_matrix_in_scenario_yellow|          |
+-------------------------------------+--------------------------------------------+

Green
-----
The ``green`` scenario...

+-------------------------------------+--------------------------------------------+
| Matrix ``G``                        | Matrix ``B``                               |
+=====================================+============================================+
| |img_G_matrix_in_scenario_green|    | |img_B_matrix_in_scenario_green|           |
+-------------------------------------+--------------------------------------------+


Black
-----
The ``black`` scenario...

+-------------------------------------+--------------------------------------------+
| Matrix ``G``                        | Matrix ``B``                               |
+=====================================+============================================+
| |img_G_matrix_in_scenario_black|    | |img_B_matrix_in_scenario_black|           |
+-------------------------------------+--------------------------------------------+


.. |BlackStyle| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

.. |TestStatus| image:: https://github.com/cherenkov-plenoscope/flux_sensitivity/actions/workflows/test.yml/badge.svg?branch=main
   :target: https://github.com/cherenkov-plenoscope/flux_sensitivity/actions/workflows/test.yml

.. |PyPiStatus| image:: https://img.shields.io/pypi/v/flux-sensitivity-sebastian-achim-mueller
   :target: https://pypi.org/project/flux-sensitivity-sebastian-achim-mueller/

.. |img_irf_background_rate_onregion| image:: flux_sensitivity/tests/resources/cta/plot/irf_background_rate_onregion.jpg

.. |img_irf_signal_area| image:: flux_sensitivity/tests/resources/cta/plot/irf_signal_area.jpg

.. |img_irf_probability_reco_given_true| image:: flux_sensitivity/tests/resources/cta/plot/irf_probability_reco_given_true.jpg


.. |img_G_matrix_in_scenario_blue| image:: flux_sensitivity/tests/resources/cta/plot/G_matrix_in_scenario_blue.jpg

.. |img_G_matrix_in_scenario_yellow| image:: flux_sensitivity/tests/resources/cta/plot/G_matrix_in_scenario_yellow.jpg

.. |img_G_matrix_in_scenario_green| image:: flux_sensitivity/tests/resources/cta/plot/G_matrix_in_scenario_green.jpg

.. |img_G_matrix_in_scenario_black| image:: flux_sensitivity/tests/resources/cta/plot/G_matrix_in_scenario_black.jpg


.. |img_B_matrix_in_scenario_blue| image:: flux_sensitivity/tests/resources/cta/plot/B_matrix_in_scenario_blue.jpg

.. |img_B_matrix_in_scenario_yellow| image:: flux_sensitivity/tests/resources/cta/plot/B_matrix_in_scenario_yellow.jpg

.. |img_B_matrix_in_scenario_green| image:: flux_sensitivity/tests/resources/cta/plot/B_matrix_in_scenario_green.jpg

.. |img_B_matrix_in_scenario_black| image:: flux_sensitivity/tests/resources/cta/plot/B_matrix_in_scenario_black.jpg


.. |img_signal_area_in_scenario_blue| image:: flux_sensitivity/tests/resources/cta/plot/signal_area_in_scenario_blue.jpg

.. |img_signal_area_in_scenario_yellow| image:: flux_sensitivity/tests/resources/cta/plot/signal_area_in_scenario_yellow.jpg

.. |img_signal_area_in_scenario_green| image:: flux_sensitivity/tests/resources/cta/plot/signal_area_in_scenario_green.jpg

.. |img_signal_area_in_scenario_black| image:: flux_sensitivity/tests/resources/cta/plot/signal_area_in_scenario_black.jpg


.. |img_background_rate_in_scenario_blue| image:: flux_sensitivity/tests/resources/cta/plot/background_rate_in_scenario_blue.jpg

.. |img_background_rate_in_scenario_yellow| image:: flux_sensitivity/tests/resources/cta/plot/background_rate_in_scenario_yellow.jpg

.. |img_background_rate_in_scenario_green| image:: flux_sensitivity/tests/resources/cta/plot/background_rate_in_scenario_green.jpg

.. |img_background_rate_in_scenario_black| image:: flux_sensitivity/tests/resources/cta/plot/background_rate_in_scenario_black.jpg
