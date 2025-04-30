# u-track 2.5

![Alt Text](img/utrack.gif?raw=true)

u-track is a multiple-particle tracking MATLAB software that is designed to (1) track dense particle fields, (2) close gaps in particle trajectories resulting from detection failure, and (3) capture particle merging and splitting events resulting from occlusion or genuine aggregation and dissociation events. Its core is based on formulating correspondence problems as linear assignment problems and searching for a globally optimal solution.

- Version 2.5 adds support for running on Mac with Apple silicon chips (tested on MacBook M4 Pro using MATLAB 2024b and GSL 2.8).
- Version 2.4 "Movie selection" graphic user interface has been renamed to "u-quantify", user can use the new Matlab command `u_quantify` to launch the GUI and start to run u-track.
- Version 2.3 runs much faster on large/long movies (improvements made by Carmen Klein Herenbrink and Brian Devree; reference coming soon).
- Version 2.2 adds parallel processing functionality for multi-movie datasets when using the GUI.
- Version 2.1 enables the analysis of movies stored on an OMERO server.
- Version 2.0 includes two new tracking applications: microtubule plus-end tracking (previously distributed as plusTipTracker) and nuclei tracking.
- A third optional processing step has been added to the analysis workflow, track analysis, with two methods: motion analysis and microtubule plus-end classification.

For more information, please see [**Robust single-particle tracking in live-cell time-lapse sequences**](http://www.nature.com/nmeth/journal/v5/n8/full/nmeth.1237.html), *Nature Methods*, 2008, 5, 695–702, written by Khuloud Jaqaman, Dinah Loerke, Marcel Mettlen, Hirotaka Kuwata, Sergio Grinstein, Sandra L Schmid, [Gaudenz Danuser](https://www.danuserlab-utsw.org/).

Besides basic particle tracking, the software supports the features for tracking microtubule plus end markers described in [**plusTipTracker: Quantitative image analysis software for the measurement of microtubule dynamics**](https://www.ncbi.nlm.nih.gov/pubmed/21821130), *J. Struct. Biol.*, 2011, 176(2):168-84, written by Kathryn T Applegate, Sebastien Besson, Alexandre Matov, Maria H Bagonis, Khuloud Jaqaman, [Gaudenz Danuser](https://www.danuserlab-utsw.org/); and for tracking fluorescently-labeled cell nuclei described in [**Substrate stiffness regulates cadherin-dependent collective migration through myosin-II contractility**](https://www.ncbi.nlm.nih.gov/pubmed/23091067), *J. Cell Biol.*, 2012, 199(3):545-63, written by Mei Rosa Ng, Achim Besser, [Gaudenz Danuser](https://www.danuserlab-utsw.org/), Joan S Brugge.

A containerized u-track software is available to [download](https://hub.docker.com/repository/docker/jennyzouutsw/utrack2.3) and can be run in Docker following [this instructions](https://github.com/JennyZouUTSW/auto-docker/blob/develop/utrack2018b/README.md).

[Danuser Lab Website](https://www.danuserlab-utsw.org/)

[Jaqaman Lab Website](https://www.utsouthwestern.edu/labs/jaqaman/)

[Software Links](https://github.com/DanuserLab/)
