Estimate the Sensitivity of your Instrument
===========================================
|TestStatus| |PyPiStatus| |BlackStyle| 

Made for astronomy with gamma-rays in the atmospheric Cherenkov-method.

This tool estimates either the integral or differential sensitivity of your instrument.
This is integral or differential with respect to the energy of the gamma-rays.

This tool puts special emphasis on instruments which have a poor reconstruction of the gamma-rays energies.
When the instrument's reconstructed energies of the gamma-rays differ from the true energies, it turns out there are multiple correct, but nervertheless different ways to express the instrument's differential sensitivity.
This tool offers, and discusses these different ways.

+-------------------------------------+--------------------------------------------+-----------+
| blue                                | yellow                                     | Header 3  |
+=====================================+============================================+===========+
| |img_signal_area_in_scenario_blue|  | |img_signal_area_in_scenario_yellow|       | column 3  |
+-------------------------------------+--------------------------------------------+-----------+

.. |BlackStyle| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

.. |TestStatus| image:: https://github.com/cherenkov-plenoscope/flux_sensitivity/actions/workflows/test.yml/badge.svg?branch=main
   :target: https://github.com/cherenkov-plenoscope/flux_sensitivity/actions/workflows/test.yml

.. |PyPiStatus| image:: https://img.shields.io/pypi/v/flux-sensitivity-sebastian-achim-mueller
   :target: https://pypi.org/project/flux-sensitivity-sebastian-achim-mueller/

.. |img_signal_area_in_scenario_blue| image:: flux_sensitivity/tests/resources/cta/plot/signal_area_in_scenario_blue.jpg

.. |img_signal_area_in_scenario_yellow| image:: flux_sensitivity/tests/resources/cta/plot/signal_area_in_scenario_yellow.jpg

.. |img_signal_area_in_scenario_green| image:: flux_sensitivity/tests/resources/cta/plot/signal_area_in_scenario_green.jpg

.. |img_signal_area_in_scenario_black| image:: flux_sensitivity/tests/resources/cta/plot/signal_area_in_scenario_black.jpg
