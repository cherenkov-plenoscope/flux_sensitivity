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

Instrument Response Function
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
Each scenario results in a different estimate for the instrument's differential sensitivity.
The scenarios are named using colors to point out their differences while not implying a hierachy.
We compiled this list of scenarios from what we found in the wild.

The figure below shows the settle differences of the scenarios.

+-----------------------+-----------------------+-----------------------+-----------------------+
| blue                  | yellow                | green                 | black                 |
+=======================+=======================+=======================+=======================+
| |img_diff_sens_blue|  | |img_diff_sens_yellow|| |img_diff_sens_green| | |img_diff_sens_black| |
+-----------------------+-----------------------+-----------------------+-----------------------+
| The differential sensitivity computed by this package is in black.                            |
| The differential sensitivity provided by CTA is in blue                                       |
| (The blue bars are canted because CTA only provides these when multiplying the                |
| flux-axis with the energy to some power).                                                     |
| For reference, the differential sensitivity of Fermi-LAT (10 years) is shown in orange,       |   
| and the flux of the Crab-nebula (``1e0``, ``1e-1``, ``1e-2``, ``1e-3``) is in dashed lines.   |
+-----------------------------------------------------------------------------------------------+

Each scenario is represented by two matrices ``G`` and ``B``.
Matrix ``G`` defines how a scenario takes the effective area for the signal into account,
and matrix ``B`` defines how a scenario takes the rate of background into account.

Blue
----
The ``blue`` scenario handles the instrument's non perfec reconstruction in energy by simply ignoring it.
Thus its matrix ``G`` is the unit-matrix with ones on the diagonal, and zeors everywhere else.
As a result, this scenarios area for the signal is just the area for the signal estimated in the instrument's response-function.
The matrix ``B`` is also the unit-matrix so that this scenario's rate of background is just the one estimated in instrument's response-function.
The obvious problem with this is of course, that the computed differential sensitivity will be poor when the instrumen's confusion in energy is significant.
The big advantage of this scenario is, that its energy-axis actually is the true gamma-ray-energy.

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

|img_diff_sens_blue|

Yellow
------
The ``yellow`` scenario handles the non-perfect reconstruction in energy by allowing true gamma-rays which where falsely
reconstructed to have energies outside of a certain energy-bin to still contribute to the signal found in this energy-bin
based on the probabilities found in the instruments confusion in energy.


+-------------------------------------+--------------------------------------------+
| Matrix ``G``                        | Matrix ``B``                               |
+=====================================+============================================+
| |img_G_matrix_in_scenario_yellow|   | |img_B_matrix_in_scenario_yellow|          |
+-------------------------------------+--------------------------------------------+

+-------------------------------------+--------------------------------------------+
| Area for signal                     | Rate of background                         |
+=====================================+============================================+
| |img_signal_area_in_scenario_yellow|||img_background_rate_in_scenario_yellow|    |
+-------------------------------------+--------------------------------------------+

|img_diff_sens_yellow|

Green
-----
The ``green`` scenario...

+-------------------------------------+--------------------------------------------+
| Matrix ``G``                        | Matrix ``B``                               |
+=====================================+============================================+
| |img_G_matrix_in_scenario_green|    | |img_B_matrix_in_scenario_green|           |
+-------------------------------------+--------------------------------------------+

+-------------------------------------+--------------------------------------------+
| Area for signal                     | Rate of background                         |
+=====================================+============================================+
| |img_signal_area_in_scenario_green| | |img_background_rate_in_scenario_green|    |
+-------------------------------------+--------------------------------------------+

|img_diff_sens_green|


Black
-----
The ``black`` scenario was proposed by Werner Hofmann and takes a different approach.
Instead of altering the area of the signal, this scenario alters the rate of the background.
The black scenario widens the energy-range in the background to the range required to collect one sigma (68%) of the signal.
This means that matrix ``B`` now collects contributions from multiple bins in reconstructed gamma-ray-energy.
The wider range in energy is estimated using the instruments confusion in energy by estimating the range in reconstructed gamma-ray-energy which contains 68% of the gamma-rays.
To represent the containment of 68% in the signal, this scenarios area in signal uses a matrix ``G`` with the elements on its diagonal being ``0.68``.
The advantage here is, that matrix ``G`` has only zeros off its diagonal and thus the black scenario can show the true gamma-rays-energy on its energy-axis. 

+-------------------------------------+--------------------------------------------------------------------------+
| Matrix ``G``                        | Matrix ``B``                                                             |
+=====================================+==========================================================================+
| |img_G_matrix_in_scenario_black|    | |img_B_matrix_in_scenario_black|                                         |
+-------------------------------------+--------------------------------------------------------------------------+
| Elements on diagonal are ``0.68``.  | At low energies, the range in energy is wider to collect enough signal.  |
+-------------------------------------+--------------------------------------------------------------------------+

+-------------------------------------+--------------------------------------------+
| Area for signal                     | Rate of background                         |
+=====================================+============================================+
| |img_signal_area_in_scenario_black| | |img_background_rate_in_scenario_black|    |
+-------------------------------------+--------------------------------------------+

|img_diff_sens_black|

Critical Rate
=============
Independent of the scenarios we listed, one additional degree of freedom when computing a differential sensitivity is how one computes the critical rate which is required in order to claim a detection.


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


.. |img_diff_sens_blue| image:: flux_sensitivity/tests/resources/cta/plot/sed_style_portal/differential_sensitivity_blue.jpg

.. |img_diff_sens_yellow| image:: flux_sensitivity/tests/resources/cta/plot/sed_style_portal/differential_sensitivity_yellow.jpg

.. |img_diff_sens_green| image:: flux_sensitivity/tests/resources/cta/plot/sed_style_portal/differential_sensitivity_green.jpg

.. |img_diff_sens_black| image:: flux_sensitivity/tests/resources/cta/plot/sed_style_portal/differential_sensitivity_black.jpg
