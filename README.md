# Permutation-based pathway enrichment analysis

Python tools to perform a permutation-based pathway enrichment analysis. Currently supporting KEGG pathways.

# Usage

```python
from pathwayenrichment.representation import ClusterPermutator
from pathwayenrichment.databaseparser import KEGGPathwayParser
from pathwayenrichment.utils import randomPartition
```

First, let's download the KEGG database for Dokdonia, a marine bacterium. To this end, we employ KEGG's entry code for Dokdonia (dok). We will then parse the database to obtain a list of genes and associated cellular pathways and systems.


```python
KEGGparser = KEGGPathwayParser.fromKEGGidentifier('dok', only_curated_pathways=True)
gene_pathways, gene_systems = KEGGparser.getGenePathways()
system_pathways = KEGGparser.getSystemPathways()
gene_info = KEGGparser.getGeneInfoFromKEGGorthology()
gene_list = list(gene_pathways.keys())
print(f'There are a total of {len(gene_list)} genes')
```

    There are a total of 786 genes


Now, we simulate a set of gene clusters to perform a pathway enrichment analysis on them. To this end, we will randomly partition the set of genes into clusters.


```python
genes_under_study = gene_list[:300]
clusters = dict(zip(
    ['A', 'B', 'C', 'D'],
    randomPartition(gene_list, bin_sizes=[75, 25, 150, 50])
))
```

Now we are ready to instantiate a ClusterPermutator to run the enrichment analysis. We will permute the total set of genes to form new random clusters 10000 times, our sample size to compute the sample p-value.


```python
permutator = ClusterPermutator(clusters, gene_pathways, system_pathways)
res = permutator.sampleClusterPermutationSpace(sample_size=10000, n_processes=4)
```

    Finished permutation sampling



```python
# Here are the first 10 pathways with lowest sample p-value
{k:v for k,v in list(res['pathway']['A'].items())[:10]}
```




    {'03018 RNA degradation [PATH:dok03018]': (0.2777777777777778, 0.0484),
     '00020 Citrate cycle (TCA cycle) [PATH:dok00020]': (0.18181818181818182,
      0.0691),
     '02020 Two-component system [PATH:dok02020]': (0.2, 0.1527),
     '00541 O-Antigen nucleotide sugar biosynthesis [PATH:dok00541]': (0.19047619047619047,
      0.1641),
     '03060 Protein export [PATH:dok03060]': (0.2, 0.1683),
     '02024 Quorum sensing [PATH:dok02024]': (0.14814814814814814, 0.218),
     '00520 Amino sugar and nucleotide sugar metabolism [PATH:dok00520]': (0.14285714285714285,
      0.2211),
     '02010 ABC transporters [PATH:dok02010]': (0.15, 0.2422),
     '00040 Pentose and glucuronate interconversions [PATH:dok00040]': (0.3333333333333333,
      0.25),
     '00053 Ascorbate and aldarate metabolism [PATH:dok00053]': (0.2, 0.25)}



Here, we see the 10 pathways with lowest sample p-values within cluster _A_.
