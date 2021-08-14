import re
import io
import numpy as np
import pandas as pd
import json
import requests
import warnings


def downloadKEGGpathwaysForID(KEGG_entry_ID: str) -> dict:
    """
    KEGG_estry_ID: genome entry ID or organism code within KEGG's database.
    """
    url = f'http://rest.kegg.jp/list/pathway/{KEGG_entry_ID}'
    text = requests.get(url).text
    columns = ['dok', 'pathway name']
    df = pd.read_csv(io.StringIO(text), sep='\t', header=None, names=columns)
    df = df.applymap(lambda s: s.replace('path:', ''))
    pathway_ids, pathway_names = df.values[:, 0], df.values[:, 1]
    return dict(zip(pathway_ids, pathway_names))


class KEGGPathwayParser:
    """
    Methods to parse KEGG pathway orthology
    """
    def __init__(self, KEGG_dict):
        """
        kid: KEGG identifier for desired organism
        """
        self._KEGG_dict = KEGG_dict
        self._gene_kos = None

    @classmethod
    def fromKEGGidentifier(cls, kid, only_curated_pathways=True):
        """
        Download KEGG functional orthology for specified KEGG identifier (kid).
        Pathways can be restricted to those that have been curated, as retrieved by the
        KEGG API.
        """
        kegg_pathways = cls._downloadKEGGorthology(kid)
        KEGG_dict = cls._getKEGGdict(kegg_pathways)
        if only_curated_pathways:
            return cls(cls.filterKEGGorthology(KEGG_dict, kid))
        else:
            return cls(KEGG_dict)

    @classmethod
    def fromJSON(self, path_to_json):
        """
        """

    @staticmethod
    def filterKEGGorthology(KEGG_dict, kid):
        """
        Filter organism-specific KEGG dict by curated pathways
        """
        filtered_KEGG_dict = {}
        curated_pathways = downloadKEGGpathwaysForID(kid)
        for supersystem_id, supersystem in KEGG_dict.items():
            filtered_KEGG_dict[supersystem_id] = {}
            for system_id, system in supersystem.items():
                pathways_to_keep = []
                for pathway_id, pathway in system.items():
                    if any([ko in pathway_id for ko in curated_pathways.keys()]):
                        pathways_to_keep.append(pathway_id)
                if pathways_to_keep:
                    filtered_KEGG_dict[supersystem_id][system_id] = {
                        pathway_id: pathway for pathway_id, pathway in system.items()
                        if pathway_id in pathways_to_keep
                    }
        return filtered_KEGG_dict

    @staticmethod
    def _downloadKEGGorthology(kid: str, path_to_save_json: str = None):
        """
        Download KEGG BRITE functional hierarchy
        Hierarchy can be downloaded for a specific organism, as long as it is
        listed in the KEGG database
        """
        if kid is None:
            kid = 'ko'
        url = f'http://rest.kegg.jp/get/br:{kid}00001/json'
        file = requests.get(url, allow_redirects=True)
        if file.status_code == 404:
            warnings.warn(f'Organism not found in KEGG. Downloading generic orthology')
            url = f'http://rest.kegg.jp/get/br:ko00001/json'
            file = requests.get(url, allow_redirects=True)
        kegg_pathways = json.loads(file.content)['children']
        if path_to_save_json is not None:
            with open(path_to_save_json, 'w') as outfile:
                json.dump(kegg_pathways, outfile)
        return kegg_pathways

    def _downloadKEGGgenePathwaysForID(self, KEGG_entry_ID: str) -> dict:
        """
        KEGG_entry_ID: genome entry ID or organism code within KEGG's database.
        """
        url = f'http://rest.kegg.jp/link/{KEGG_entry_ID}/pathway'
        text = requests.get(url).text
        columns = ['ko', 'gene_id']
        df = pd.read_csv(io.StringIO(text), sep='\t', header=None, names=columns)
        kid = re.match('^[^:]*', df.iloc[0, 1]).group()
        df = df.applymap(lambda s: re.sub('^[^:]*:', '', s))
        df['ko'] = df['ko'].map(lambda s: s.replace(kid, 'ko'))
        df = df.reindex(columns=['gene_id', 'ko'])
        genes, kos = df.values[:, 0], df.values[:, 1]
        return dict(zip(genes, kos))

    def getKEGGgenePathwaysOLD(self, KEGG_entry_ID: str) -> dict:
        """
        KEGG_estry_ID: genome entry ID or organism code within KEGG's database.
        """
        gene_pathways, gene_systems = {}, {}
        if self._gene_kos is None:
            self._gene_kos = self._downloadKEGGgenePathwaysForID(KEGG_entry_ID)

        for gene_id, ko in self._gene_kos.items():
            try:
                ko_name = self.getPathwayOrthologyFromKO(ko)
            except Exception:
                ko_name = {'subsystem': '', 'system': ''}
            gene_pathways[gene_id] = ko_name['subsystem']
            gene_systems[gene_id] = ko_name['system']
        return (gene_pathways, gene_systems)

    @staticmethod
    def extractKoID(KEGG_pathway_name: str) -> str:
        """
        Extract KEGG's ko identifier from pathway name string.
        """
        try:
            return 'ko' + re.search('ko(.*?)\]', KEGG_pathway_name).group(1)
        except Exception:
            return ''

    @staticmethod
    def simplifyKEGGpathwayName(KEGG_pathway_name: str) -> str:
        """
        Remove ko identifier from KEGG pathway name
        """
        return re.sub('\d{5}', '', re.sub('\[.*?\]', '', KEGG_pathway_name)).strip()

    def parseKEGGPathwayJSON(self, path_to_json: str, simplified_names=False,
                             supersystems_to_keep: list=None) -> dict:
        """
        Parse KEGG pathway JSON into a nested dictionary
        """
        with open(path_to_json) as json_file:
            kegg_pathways = json.load(json_file)['children']
        if supersystems_to_keep is not None:
            kegg_pathways = [supersystem for supersystem in kegg_pathways
                             if supersystem['name'] in supersystems_to_keep]

        KEGG_dict = self._getKEGGdict(kegg_pathways)
        ko_orthology_dict = self._getKoOrthologyDict(kegg_pathways)
        return (KEGG_dict, ko_orthology_dict)

    @staticmethod
    def _getKEGGdict(kegg_pathways: dict) -> dict:
        """
        Obtain a simplified python dictionry from the direct JSON load object.
        """
        KEGG_dict = {}
        for supersystem in kegg_pathways:
            KEGG_dict[supersystem['name']] = {}
            for system in supersystem['children']:
                KEGG_dict[supersystem['name']][system['name']] = {}
                for pathway in system['children']:
                    if 'children' in pathway.keys():
                        KEGG_dict[supersystem['name']][
                            system['name']][pathway['name']] = []
                        for gene in pathway['children']:
                            KEGG_dict[supersystem['name']][
                                system['name']][pathway['name']].append(gene['name'])
        return KEGG_dict

    def _getGeneListFromKEGGorthology(self) -> list:
        """
        Retrieve list of included genes in (organism-specific) KEGG orthology
        """
        genes = []
        for supersystem_id, supersystem in self._KEGG_dict.items():
            for system_id, system in supersystem.items():
                for pathway_id, pathway in system.items():
                    for gene_str in pathway:
                        gene_id = gene_str.split(' ')[0]
                        genes.append(gene_id)
        return np.unique(genes).tolist()

    def getGeneInfoFromKEGGorthology(self) -> list:
        """
        Retrieve gene product info of included genes in (organism-specific) KEGG orthology
        """
        gene_info = {}
        for supersystem_id, supersystem in self._KEGG_dict.items():
            for system_id, system in supersystem.items():
                for pathway_id, pathway in system.items():
                    for gene_str in pathway:
                        gene_id = gene_str.split(' ')[0]
                        if gene_id not in gene_info.keys():
                            gene_meta = gene_str.split('\t')[1]
                            gene_info[gene_id] = gene_meta
        return gene_info

    def getGenePathways(self):
        """
        Extract assigned KEGG pathways and systems for each gene in organism-specific
        KEGG function orthology.
        """
        genes = self._getGeneListFromKEGGorthology()
        gene_pathways = {gene_id: [] for gene_id in genes}
        gene_systems = {gene_id: [] for gene_id in genes}

        for supersystem_id, supersystem in self._KEGG_dict.items():
            for system_id, system in supersystem.items():
                for pathway_id, pathway in system.items():
                    for gene_str in pathway:
                        gene_id = gene_str.split(' ')[0]
                        gene_pathways[gene_id].append(pathway_id)
                        gene_systems[gene_id].append(system_id)
        return (gene_pathways, gene_systems)

    def getSystemPathways(self):
        """
        Retrieve pathways (subsystems) within each KEGG system
        """
        system_pathways = {}
        for supersystem_id, supersystem in self._KEGG_dict.items():
            for system_id, system in supersystem.items():
                pathways = [path_id for path_id in system.keys()]
                if pathways:
                    system_pathways[system_id] = pathways
        return system_pathways

    def _getKoOrthologyDict(self, kegg_pathways: dict) -> dict:
        """
        Obtain dictionary containing KEGG pathway classification for each
        KEGG pathway id (koXXXXX)
        simplified_names: self.simplifyKEGGpathwayName
        """
        ko_orthology = {}

        for supersystem in kegg_pathways:
            for system in supersystem['children']:
                for subsystem in system['children']:
                    try:
                        ko_id = self.extractKoID(subsystem['name'])
                        ko_path_name = subsystem['name']
                    except Exception:
                        pass
                    ko_orthology[ko_id] = {
                        'subsystem': ko_path_name,
                        'system': system['name'],
                        'supersystem': supersystem['name']
                    }
        return ko_orthology

    def getPathwayOrthologyFromKO(self, ko_identifier: str) -> dict:
        """
        Extract KEGG orthology for given ko identifier
        subsystem (pathway) < system < supersystem)
        """
        try:
            return self._ko_orthology[ko_identifier]
        except Exception:
            raise ValueError('KO not in local database')

    def getGenePathwaysDictionaryFromEGGNOG(self, path_to_eggnog_data: str) -> dict:
        """
        Make dictionary where keys are gene IDs and values lists of pathways.
        """
