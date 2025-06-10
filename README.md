# TreeSort

## Overview

The idea behind TreeSort is the observation that *if there is no reassortment, then the evolutionary histories of different segments should be identical*. TreeSort then uses a phylogenetic tree for one segment (e.g., the HA influenza A virus segment) as an evolutionary hypothesis for another segment (e.g., the NA segment). We will refer to the first segment as the *reference* and the second segment as the *challenge*. By trying to fit the sequence alignment of the challenge segment to the reference tree, TreeSort identifies points on that tree, where this evolutionary hypothesis breaks. The "breaking" manifests in the mismatch between the divergence time on the reference tree (e.g., 1 year divergence between sister clades) and an unlikely high number of substitutions in the challenge segment that are required to explain the reference tree topology under the null hypothesis of no reassortment.

TreeSort has demonstrated very high accuracy in reassortment inference in simulations (manuscript in preparation). TreeSort can process datasets with tens of thousands of virus strains in just a few minutes and can scale to very large datasets with hundreds of thousands of strains. (This overview is from [https://github.com/flu-crew/TreeSort/blob/main/README.md](https://github.com/flu-crew/TreeSort/blob/main/README.md))

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

There is one application service specification defined here:

1. [TreeSort](app_specs/TreeSort.md): The TreeSort tool infers both recent and ancestral reassortment events along the branches of a phylogenetic tree of a fixed genomic segment. It uses a statistical hypothesis testing framework to identify branches where reassortment with other segments has occurred and reports these events.

The code in this module provides the BV-BRC application service wrapper scripts for the TreeSort service as well
as some backend utilities:

| Script name | Purpose |
| ----------- | ------- |
| [App-TreeSort.pl](service-scripts/App-TreeSort.pl) | App script for the [TreeSort service](https://www.bv-brc.org/docs/quick_references/services/treesort_service.html) |

## See also

* [TreeSort Service](https://www.bv-brc.org/app/TreeSort)
* [Quick Reference](https://www.bv-brc.org/docs/quick_references/services/treesort_service.html)
* [TreeSort Service Tutorial](https://www.bv-brc.org/docs/tutorial/treesort/treesort.html)

## References


1. Alexey Markin, Catherine A. Macken, Amy L. Baker, Tavis K. Anderson, "Revealing reassortment in influenza A viruses with TreeSort"
[bioRxiv 2024.11.15.623781](https://www.biorxiv.org/content/10.1101/2024.11.15.623781v1); doi: https://doi.org/10.1101/2024.11.15.623781

2. GitHub: [https://github.com/flu-crew/TreeSort](https://github.com/flu-crew/TreeSort)




