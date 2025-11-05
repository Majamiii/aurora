# Aurora Borealis

AuroraSim — Modified Version with Electron Scattering

This is a modified version of AuroraSim
, originally developed by Kyle Mills. [https://github.com/aurora-sim/Aurora-Sim](https://github.com/millskyle/AuroraSim)
I added electron scattering to the simulation, which noticeably improved the accuracy and realism of the auroral emission results.
The original paper’s test data had some issues, so I redid the runs to get more reliable results.

Project done at Petnica Science Seminar, September 2024


✨ What's New

Electron Scattering Implemented
Added scattering effects for energetic electrons, improving how the simulation handles particle interactions in the upper atmosphere.
This leads to more realistic brightness distributions and emission profiles.

Tests Re-run with Correct Data
The original test setup had inconsistencies, so all simulations were repeated using corrected parameters and data files.
Results are now more consistent with observed auroral behavior.

Improved Output Quality
Smoother intensity transitions and better agreement between model predictions and real auroral measurements.


The papers for the conferences and texts about results can be found in the /papers directory.
The code is in src.
