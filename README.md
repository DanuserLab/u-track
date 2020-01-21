# u-track 2.0

![Alt Text](img/utrack.gif?raw=true)

u-track is a multiple-particle tracking MATLAB software that is designed to (1) track dense particle fields, (2) close gaps in particle trajectories resulting from detection failure, and (3) capture particle merging and splitting events resulting from occlusion or genuine aggregation and dissociation events. Its core is based on formulating correspondence problems as linear assignment problems and searching for a globally optimal solution.

- Version 2.3 runs much faster on large/long movies (improvements made by Carmen Klein Herenbrink and Brian Devree; reference coming soon).
- Version 2.2 adds parallel processing functionality for multi-movie datasets when using the GUI.
- Version 2.1 enables the analysis of movies stored on an OMERO server.
- Version 2.0 includes two new tracking applications: microtubule plus-end tracking (previously distributed as plusTipTracker) and nuclei tracking.
- A third optional processing step has been added to the analysis workflow, track analysis, with two methods: motion analysis and microtubule plus-end classification.

For more information, please see [Jaqaman et al., Nature Methods 5, pp. 695-702 (2008)](http://www.nature.com/nmeth/journal/v5/n8/full/nmeth.1237.html). Besides basic particle tracking, the software supports the features described in [Applegate et al. J. Struct. Biol. 176(2):168-84. 2011](https://www.ncbi.nlm.nih.gov/pubmed/21821130) for tracking microtubule plus end markers; and in [Ng et al. J. Cell Biol. 199(3):545-63. 2012](https://www.ncbi.nlm.nih.gov/pubmed/23091067) for tracking fluorescently-labeled cell nuclei.


[Danuser Lab Website](https://www.utsouthwestern.edu/labs/danuser/)

[Jaqaman Lab Website](https://www.utsouthwestern.edu/labs/jaqaman/)

[Software Links](https://www.utsouthwestern.edu/labs/danuser/software/)
