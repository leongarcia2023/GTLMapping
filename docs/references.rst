References
==========

Core MIREX method
-----------------

Butler, M. J. & Tan, J. C. (2009), *Mid-Infrared Extinction Mapping
of Infrared Dark Clouds: Probing the Initial Conditions for Massive
Stars and Star Clusters*, ApJ 696, 484.
`doi:10.1088/0004-637X/696/1/484
<https://doi.org/10.1088/0004-637X/696/1/484>`__

Butler, M. J. & Tan, J. C. (2012), *Mid-Infrared Extinction Mapping
of Infrared Dark Clouds. II. The Structure of Massive Starless Cores
and Clumps*, ApJ 754, 5.
`doi:10.1088/0004-637X/754/1/5
<https://doi.org/10.1088/0004-637X/754/1/5>`__

Butler, M. J. & Tan, J. C. (2013), erratum to the 2012 paper, ApJ
766, 66.
`doi:10.1088/0004-637X/766/1/66
<https://doi.org/10.1088/0004-637X/766/1/66>`__

Butler, M. J., Tan, J. C., & Kainulainen, J. (2014), *The Darkest
Shadows: Deep Mid-Infrared Extinction Mapping of a Massive
Protocluster*, ApJL 782, L30. This paper applies local saturation
searches to deeper Spitzer imaging and reduces the independence scale
from 8 to 4 arcsec in that dataset.
`doi:10.1088/2041-8205/782/2/L30
<https://doi.org/10.1088/2041-8205/782/2/L30>`__

Cloud geometry
--------------

Simon, R., Jackson, J. M., Rathborne, J. M., & Chambers, E. T. (2006),
*A Catalog of MSX Infrared Dark Cloud Candidates*, ApJ 639, 227.
`doi:10.1086/499342 <https://doi.org/10.1086/499342>`__

Dust opacity
------------

Ossenkopf, V. & Henning, T. (1994), *Dust Opacities for Protostellar
Cores*, A&A 291, 943.
`opacity tables <https://hera.ph1.uni-koeln.de/~ossk/Jena/tables.html>`__

Rodrigo, C., Solano, E., Bayo, A., et al. (2024), Spanish Virtual
Observatory Filter Profile Service. The supplied Sgr C analysis uses
these throughput profiles for filter convolution.
`SVO Filter Profile Service <https://svo2.cab.inta-csic.es/theory/fps/>`__

Related high-dynamic-range mapping
----------------------------------

Kainulainen, J. & Tan, J. C. (2013), *High-dynamic-range extinction
mapping of infrared dark clouds*, A&A 549, A53.
This work combines NIR and MIR extinction, adopts
:math:`\tau_8=0.29\tau_K`, and discusses the roughly 30-percent
absolute MIR-opacity uncertainty.
`doi:10.1051/0004-6361/201219526
<https://doi.org/10.1051/0004-6361/201219526>`__

Lim, W. & Tan, J. C. (2014), *Far-Infrared Extinction Mapping of
Infrared Dark Clouds*, ApJL 780, L29. This paper extends the saturation
and foreground logic to 24 and 70 microns and emphasizes
wavelength-dependent noise, beam size, and opacity.
`doi:10.1088/2041-8205/780/2/L29
<https://doi.org/10.1088/2041-8205/780/2/L29>`__

Lim, W., Tan, J. C., Kainulainen, J., Ma, B., & Butler, M. J. (2016),
*The Distribution of Mass Surface Densities in a High-Mass
Protocluster*, ApJL 829, L19. The foreground/background comparison
shows that diffuse-emission subtraction can materially change derived
surface-density distributions.

JWST-era applications
----------------------

Fedriani, R., Tan, J. C., Law, C.-Y., Crowe, S., et al. (in preparation),
*The JWST-NIRCam View of Sagittarius C. IV. Mid-Infrared Extinction
(MIREX) Mapping*. The study applies the BT12 saturation prescription to
JWST F480M data with an ERR-derived noise level and adjacent background
regions. GTLMapping uses its filter-convolved opacity table. The F480M values
9.76 and 15.23 cm\ :sup:`2` g\ :sup:`-1` correspond to gas-to-dust ratios 156
and 100.

André, P., Mattern, M., Arzoumanian, D., Shimajiri, Y., et al. (2025),
*Structure and Fragmentation Scale of a Massive Star-Forming Filament
in NGC 6334: High-Resolution Mid-Infrared Absorption Imaging with
JWST*, ApJL 984, L59. It compares a large-scale median-background
method with a method calibrated to an independent submillimeter column
density map and explicitly discusses saturation.

Complementary extinction-law mapping
------------------------------------

Fahrion, K. & De Marchi, G. (2023), *Extending the Extinction Law in
30 Doradus to the Infrared with JWST*, A&A 671, L14. This is a
red-clump-star extinction-law and extinction-map analysis rather than
a diffuse-background MIREX algorithm; it informs future variable-law
work but is not implemented by GTLMapping.
`doi:10.1051/0004-6361/202346240
<https://doi.org/10.1051/0004-6361/202346240>`__
