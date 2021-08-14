
import numpy as np
import random
import multiprocessing as mp
import warnings


class PathwayRepresentation:
    """
    Methods to compute pathway and system representation in gene clusters.
    A system is a collection of pathways and each pathway a collection of genes.
    """
    def __init__(self, clusters: dict, gene_pathways: dict,
                 system_pathways: dict = None) -> None:
        self._gene_pathways = gene_pathways
        self._pathway_sizes = self.getCounts(
            [p for pathways in self._gene_pathways.values() for p in pathways]
        )
        self._pathway_ids = list(self._pathway_sizes.keys())
        self._system_pathways = system_pathways
        if system_pathways is not None:
            self.checkSystemPathways()
        self._clusters = clusters
        self._cluster_ids = list(self._clusters.keys())

    @staticmethod
    def getCounts(array_like: list, sort_by_value=True) -> list:
        """
        Get dict of array item counts. Dict sorted
        by keys unless sort_by_value set to True.
        """
        unique_cond, counts_cond = np.unique(array_like, return_counts=True)
        counts = dict(zip(unique_cond, counts_cond))
        if sort_by_value:
            return dict(sorted(
                counts.items(), key=lambda item: item[1], reverse=True))
        else:
            return counts

    def getPathwaysFromGeneList(self, gene_list: list) -> list:
        """
        Get list of pathways (not unique) represented in gene list
        """
        genes_with_assigned_pathway = [g for g in gene_list if g in self._gene_pathways]
        return [p for g in genes_with_assigned_pathway for p in self._gene_pathways[g]]

    def initializeClusterRepresentationDict(self) -> dict:
        """
        Initialize cluster pathway representation dict with keys equal to pathways
        and  values equal to 0.
        """
        return {pathway: 0 for pathway in self._pathway_ids}

    def initializeClustersRepresentationDict(self) -> dict:
        """
        Initialize a cluster representation dictionary with 0s. Dict of dicts
        with keys equal to cluster ids and pathway names within.
        """
        return {
                cluster_id: self.initializeClusterRepresentationDict()
                for cluster_id in self._cluster_ids
        }

    def computeClusterPathwayRepresentation(self, cluster_pathway_sizes: dict) -> dict:
        pathway_representation = self.initializeClusterRepresentationDict()
        for pathway, total_counts in self._pathway_sizes.items():
            if pathway in cluster_pathway_sizes.keys():
                pathway_representation[pathway] = cluster_pathway_sizes[pathway] / total_counts
        return pathway_representation

    def getClustersPathwayRepresentation(self, clusters: dict) -> dict:
        """
        Compute fraction of genes of each pathway represented in each cluster
        gene_pathways: dict with keys being gene ids and values are lists of pathways (str)
        to which each gene belongs.
        """
        clusters_pathway_representation = {}
        for cluster_id, cluster in clusters.items():
            cluster_pathway_sizes = self.getCounts(
                self.getPathwaysFromGeneList(cluster)
            )
            clusters_pathway_representation[cluster_id] = self.computeClusterPathwayRepresentation(
                cluster_pathway_sizes)
        return clusters_pathway_representation

    def checkSystemPathways(self):
        """
        Check if all pathways defined within systems are included in the
        pathway set gathered for current gene set. Used to check if input
        systems_pathway_dict is consistent with pathway data.
        """
        if self._system_pathways is None:
            raise ValueError('Undefined system_pathways dictionary')
        system_included_pathways = {}
        ill_defined_systems = []
        for system, pathways in self._system_pathways.items():
            included_pathways = [
                pathway for pathway in pathways if pathway in self._pathway_ids
            ]
            system_included_pathways[system] = included_pathways
            if len(included_pathways) < len(pathways):
                ill_defined_systems.append(system)
        if ill_defined_systems:
            warnings.warn('Some systems contained pathways without assigned genes. Removing these pathways from dictionary')
        self._system_pathways = system_included_pathways

    def getClustersSystemRepresentation(self,
                                        clusters_pathway_representation: dict) -> dict:
        """
        Compute system representation out of pathway
        representation values (a system is a set of pathways or subsystems)

        R_sys = (1/N_sys) * Sum_path (R_path * N_path), for all pathways in system,
        where: N_sys, number of genes in system
               R_sys, system representation in cluster (genes in cluster/total genes in system)
               N_path, number of genes in pathway
               R_path, pathway representation in cluster (genes in cluster/total genes in pathway)
        """
        if self._system_pathways is None:
            raise ValueError('Undefined system_pathways dictionary')
        clusters_system_representation = {}
        for cluster_id, cluster_rep in clusters_pathway_representation.items():
            clusters_system_representation[cluster_id] = {}
            for system_id, pathways in self._system_pathways.items():
                system_size = sum([self._pathway_sizes[pathway] for pathway in pathways])
                system_rep = (1 / system_size) * sum([cluster_rep[pathway] * self._pathway_sizes[pathway] for pathway in pathways])
                clusters_system_representation[cluster_id][system_id] = system_rep
        return clusters_system_representation


class ClusterPermutator(PathwayRepresentation):
    """
    Methods to run a permutation analysis on gene clusters
    pathway_representation: class instance of PathwayRepresentation
    NOTE: The inits of these two classes are redundant. Improve it.
    """
    def __init__(self, clusters: dict, gene_pathways: dict,
                 system_pathways: dict = None) -> None:
        super().__init__(clusters, gene_pathways, system_pathways)
        self._cluster_pathway_representation = self.getClustersPathwayRepresentation(
            self._clusters)
        if system_pathways is not None:
            self._cluster_system_representation = self.getClustersSystemRepresentation(
                self._cluster_pathway_representation
            )
        self._gene_list = [g for c in self._clusters.values() for g in c]
        self._cluster_sizes = [len(v) for v in self._clusters.values()]
        self._cluster_representation_p_values = self.initializeClustersRepresentationDict()
        self._sample_size = None

    @staticmethod
    def randomPartition(elems: list, bin_sizes: list) -> list:
        """
        Randomly partition list elements into
        bins of given sizes
        """
        def shuffleList(elems):
            random.shuffle(elems)
            return elems

        elems = shuffleList(elems)
        partition = []
        start, end = 0, 0
        for bin_size in bin_sizes:
            end += bin_size
            partition.append(elems[start:end])
            start += bin_size
        return partition

    @staticmethod
    def copyProxyDictionary(proxydict) -> dict:
        """
        Convert ProxyDict, as returned by multiprocessing.Manager to regular dict.
        """
        return {k: v.copy() for k, v in proxydict.items()}

    def getPermutedClusterPathwayRepresentation(self, i=None):
        """
        Compute pathway representation in permuted clusters
        (Ignore 'i' argument, required by Multiprocessing.map)
        """
        permuted_clusters = dict(
            zip(self._cluster_ids, self.randomPartition(self._gene_list, self._cluster_sizes))
        )
        return self.getClustersPathwayRepresentation(permuted_clusters)

    def _aggregatePermutationPathwayResults(self, samples: list) -> dict:
        aggregated_samples = {
            cluster_id: {pathway: []
                         for pathway in self._pathway_ids}
            for cluster_id in self._cluster_ids
        }
        for sample in samples:
            for cluster_id in self._cluster_ids:
                for pathway in self._pathway_ids:
                    aggregated_samples[cluster_id][pathway].append(sample[cluster_id][pathway])
        return aggregated_samples

    def _aggregatePermutationSystemResults(self, samples: list) -> dict:
        aggregated_samples = {
            cluster_id: {system: []
                         for system in self._system_pathways.keys()}
            for cluster_id in self._cluster_ids
        }
        for sample in samples:
            for cluster_id in self._cluster_ids:
                for system in self._system_pathways.keys():
                    aggregated_samples[cluster_id][system].append(sample[cluster_id][system])
        return aggregated_samples

    def aggregatePermutationSampleResults(self, samples: list, system_type='pathway') -> dict:
        """
        Aggregate permutation sample results (generated via multiprocessing)
        into a single object
        """
        if 'pathway' in system_type.lower():
            return self._aggregatePermutationPathwayResults(samples)
        else:
            return self._aggregatePermutationSystemResults(samples)

    def _sortPathwaysByPvalue(self, cluster_pathway_pvalue: dict) -> dict:
        """
        Sort pathways within a cluster based on their sample pvalue.
        cluster_pathway_pvalue, dictionary where keys are pathway names
        and values are tuples: (representation, p-value)
        """
        pathways = cluster_pathway_pvalue
        sorted_keys = np.array(
                               list(pathways.keys())
        )[np.argsort([pvalue for rep, pvalue in pathways.values()])]
        return {k: pathways[k] for k in sorted_keys}

    @staticmethod
    def _computeSamplePvalue(distribution: list, statistic_value: float) -> float:
        """
        Compute sample p-value given null distribution and statistic value
        """
        N_equal_more_extreme = len(
            np.where(np.array(distribution) >= statistic_value)[0]
            )
        N_total = len(distribution)
        return N_equal_more_extreme / N_total

    def computePermutationSamplePvalue(self, aggregated_samples: dict,
                                       pathway_representation: dict) -> dict:
        """
        Compute sample p-value for each pathway and cluster
        """
        rep_pvalues = {cluster_id: {} for cluster_id in pathway_representation.keys()}
        for cluster_id, rep_values in pathway_representation.items():
            pathways_rep = {}
            for pathway_id, rep_value in rep_values.items():
                rep_distribution = aggregated_samples[cluster_id][pathway_id]
                p_value = self._computeSamplePvalue(rep_distribution, rep_value)
                pathways_rep[pathway_id] = (rep_value, p_value)
            rep_pvalues[cluster_id] = self._sortPathwaysByPvalue(pathways_rep)
        return rep_pvalues

    def sampleClusterPermutationSpace(self, sample_size=1000, n_processes=None) -> list:
        """
        Permute genes in clusters n times.
        NOTE: Postprocessing time can be optimized, e.g., by choosing list of lists
        instead of dictionary to store pathway representation sample. Currently it takes
        over a minute to process samples (when both pathway and system representation
        are computed)
        """
        self._sample_size = sample_size
        if n_processes is None:
            n_processes = mp.cpu_count() - 1
        pool = mp.Pool(n_processes)
        samples = pool.map(
            self.getPermutedClusterPathwayRepresentation,
            [() for _ in range(sample_size)]
        )
        pool.close()
        print('Finished permutation sampling')
        pathway_rep_result = self.copyProxyDictionary(
            self.computePermutationSamplePvalue(self.aggregatePermutationSampleResults(
                samples, system_type='pathway'), self._cluster_pathway_representation)
            )

        if self._system_pathways is None:
            return pathway_rep_result
        else:
            system_samples = [
                self.getClustersSystemRepresentation(sample) for sample in samples
                ]
            system_rep_result = self.copyProxyDictionary(self.computePermutationSamplePvalue(
                self.aggregatePermutationSampleResults(system_samples, system_type='system'),
                self._cluster_system_representation)
                )
            return {'pathway': pathway_rep_result, 'system': system_rep_result}

    # Shared memory in multiprocessing (not very efficient)

    # def getPermutedClusterPathwayRepresentationPvalue(self, shared_res_dict):
    #     """
    #     Compute pathway representation in permuted clusters and compare with
    #     actual value
    #     """
    #     permuted_representation = self.getPermutedClusterPathwayRepresentation()
    #     for cluster_id, cluster in self._cluster_pathway_representation.items():
    #         for pathway, representation in cluster.items():
    #             if permuted_representation[cluster_id][pathway] >= representation:
    #                 self._cluster_representation_p_values[cluster_id][pathway] += 1 / self._sample_size
    #
    # def sampleClusterPermutationSpaceSharedMemory(self,
    #                                               sample_size=1000,
    #                                               n_processes=None) -> list:
    #     """
    #     Permute genes in clusters n times
    #     NOTE: Shared memory doesn't seem to work properly. When using a shared dict
    #     one process heavely uses CPU and the others do not. Executing time is far greater
    #     than multiprocessing with no shared memory.
    #     """
    #     self._sample_size = sample_size
    #     manager = mp.Manager()
    #     shared_res_dict = {
    #                        cluster_id: manager.dict(self.initializeClusterRepresentationDict())
    #                        for cluster_id in self._cluster_ids
    #                        }
    #     self._cluster_representation_p_values = shared_res_dict
    #
    #     if n_processes is None:
    #         n_processes = mp.cpu_count() - 1
    #
    #     pool = mp.Pool(n_processes)
    #     pool.map(
    #         self.getPermutedClusterPathwayRepresentationPvalue,
    #         [_ for _ in range(sample_size)]
    #     )
    #     pool.close()
    #     return shared_res_dict
