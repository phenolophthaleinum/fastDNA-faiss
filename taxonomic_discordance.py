#!/usr/bin/env python3
"""
author: Jakub Barylski
Created on 10.12.2021 (v0.0.1)
"""

import json
import pickle
from pathlib import Path
from typing import Dict, List, Tuple, Any, Union
from tqdm import tqdm
from colorama import Fore, init
from timeit import default_timer as timer
import joblib
import argparse

import numpy as np

PathLike = Union[Path, str]

# ful taxon rank spectrum
# TAXONOMIC_RANKS = ('species', 'subgenus', 'genus', 'subfamily', 'family', 'superfamily', 'suborder', 'order', 'superorder', 'subclass', 'class', 'superclass', 'subphylum', 'phylum', 'superphylum', 'kingdom', 'superkingdom', 'subroot')

# selected only consistent ranks annotated in all prokaryotic lineages
TAXONOMIC_RANKS = ('species',
                   'genus',
                   'family',
                   'order',
                   'class',
                   'phylum',
                   'kingdom')


class Species:

    def __init__(self,
                 lineage: tuple,
                 ranks: tuple):
        self.lineage = lineage
        self.ranks = ranks
        self.id = self.lineage[0]
        self.genomes = []

    def taxonomic_distance(self,
                           species: 'Species') -> int:
        """
        Score taxonomic relation between two species.
        If the both are identical return 0.
        Otherwise, return number of taxonomic levels
        one need to go up to find a common taxon.
        E.g. species = 0, genus = 1, family = 2 etc.
        @return: taxonomic similarity score
        """

        compared_ids = {i for i in species.lineage}
        for rank, taxon in zip(self.ranks, self.lineage):
            if taxon in compared_ids and rank in TAXONOMIC_RANKS:
                return TAXONOMIC_RANKS.index(rank)
        return len(TAXONOMIC_RANKS) + 1


class DistanceMatrix:
    """
    Dictionary-like object with species -> genome distances
    :param master_host_dict: dictionary of host metadata stored in 'host.json'
    :return: {'species_id_0': {'genome0': distance0, (...) 'genomeN': distanceN}, species_id_1 ...}
    """

    def __init__(self, master_host_dict: Dict[str, Dict[str, Any]]):

        tmp_dict = {}

        for genome_id, genome_metadata in master_host_dict.items():
            lineage = tuple(reversed(genome_metadata['lineage_names']))
            ranks = tuple(reversed(genome_metadata['lineage_ranks']))
            assert ranks[0] == 'species'
            if lineage[0] not in tmp_dict:
                tmp_dict[lineage[0]] = Species(lineage, ranks)
            tmp_dict[lineage[0]].genomes.append(genome_id)

        remaining_species = list(tmp_dict.values())
        self.species = {}
        self.genomes = {}
        n_species = len(remaining_species)

        self.matrix = np.empty([n_species, n_species], dtype=np.int16)
        np.fill_diagonal(self.matrix, 0)

        query_index = 0
        bar = tqdm(total=n_species)
        bar.set_description('creating distance matrix')
        while remaining_species:
            query_s = remaining_species.pop(0)
            for target_index, target_s in enumerate(remaining_species, query_index + 1):
                self.matrix[query_index][target_index] = self.matrix[target_index][query_index] = query_s.taxonomic_distance(target_s)
            for genome_id in query_s.genomes:
                self.genomes[genome_id] = query_index
            self.species[query_s.id] = self.matrix[query_index]
            query_index += 1
            bar.update()
        bar.close()

    def get_distance(self, genome_id, species_id):
        """
        Retrieve a taxonomic distance between a genome and a species
        :param genome_id: e.g. NC_017548
        :param species_id:
        :return:
        """
        return self.species[species_id][self.genomes[genome_id]]

    def save(self, path: PathLike):
        """
        Serialize the matrix for further use
        :param path:
        :return:
        """
        path = Path(path)
        with path.open('wb') as out:
            pickle.dump(self, out)

    def save_parallel(self, path: PathLike, compress_param: Union[int, bool, Tuple[str, int]]):
        path = Path(path)
        with path.open('wb') as out:
            joblib.dump(self, out, compress=compress_param)


def load_matrix(path: PathLike) -> DistanceMatrix:
    path = Path(path)
    with path.open('rb') as out:
        matrix = pickle.load(out)
    return matrix


def load_matrix_parallel(path: PathLike) -> DistanceMatrix:
    path = Path(path)
    with path.open('rb') as out:
        matrix = joblib.load(out)
    return matrix


def taxonomic_discordance(distances: DistanceMatrix,
                          sorted_predictions: Dict[str, List[Tuple[str, float]]],
                          master_virus_dict: Dict[str, Dict[str, Any]],
                          top_n: int = 3):
    """
    Calculate a "taxonomic discordance" for provided predictions.
    This is metric that is a ration of:
    taxonomic distances between true and predicted hosts
    to
    all possible on the hierarchical taxon tree (worst of possible predictions)
    The calculated for "top_n" predictions extracted from probability-sorted prediction dict
    Weight of each is prediction adjusted inversely to its rank.
    :param distances: species -> genome distances form TaxonomicDiscordance.get_distances()
    :param sorted_predictions: {virus_id_0: [('host_id_BEST', p_BEST), (...), ('host_id_WORST', p_WORST)], virus_id_1 ...}
    :param master_virus_dict: dictionary of virus metadata stored in 'virus.json'
    :param top_n: number of top prediction to assess
    :return: 0-1 (0 is perfect prediction 1 is total confusion)
    """
    prediction_count = 0
    observed_shift = 0
    overall = 0
    skipped = 0
    bar = tqdm(total=len(sorted_predictions))
    bar.set_description('getting prediction-standard distances')
    for virus_id, predictions in sorted_predictions.items():
        true_host = master_virus_dict[virus_id]['host']['organism_name']  # TODO Edwards notation is VERY CONFUSING (should be master_virus_dict[virus_id]['host']['species/taxid']) !!!
        if len(predictions) > top_n:
            predictions = predictions[:top_n]
        for rank, (host_id, confidence) in enumerate(predictions, 1):
            prediction_count += 1 / rank
            overall += 1
            if host_id in distances.genomes and true_host in distances.species:
                observed_shift += distances.get_distance(host_id, true_host) / rank
            else:
                # print(f'{host_id} OR {true_host} not in the matrix (please check original host.json)')
                skipped += 1
        bar.update()

    bar.close()

    print(f'{skipped} of {overall} predictions skipped due to missing matrix values')

    max_shift = prediction_count * (len(TAXONOMIC_RANKS) + 1)

    return observed_shift / max_shift


def taxonomic_accordance(distance_dict: DistanceMatrix,
                         sorted_predictions: Dict[str, List[Tuple[str, float]]],
                         master_virus_dict: Dict[str, Dict[str, Any]],
                         top_n: int = 3):
    """
    Reverse of taxonomic discordance (1 - taxonomic_discordance)
    This is 0-1 metric that is getting closer to maximum
    when predictions align with known host taxonomy
    :param distance_dict: species -> genome distances form TaxonomicDiscordance.get_distances()
    :param sorted_predictions: {virus_id_0: [('host_id_BEST', p_BEST), (...), ('host_id_WORST', p_WORST)], virus_id_1 ...}
    :param master_virus_dict: dictionary of virus metadata stored in 'virus.json'
    :param top_n: number of top prediction to assess
    :return: 0-1 (1 is perfect prediction 0 is total confusion)
    """
    return 1 - taxonomic_discordance(distance_dict,
                                     sorted_predictions,
                                     master_virus_dict,
                                     top_n)


# uncomment for test
# if __name__ == '__main__':
#
#     # exemplar use cse
#
#     # colorama
#     init()
#
#     # global timer
#     total_start = timer()
#
#     host_json = Path('host.json')
#     with host_json.open() as hj:
#         host_dict = json.load(hj)
#
#     virus_json = Path('virus.json')
#     with virus_json.open() as hj:
#         virus_dict = json.load(hj)
#
#     # dists = DistanceMatrix(host_dict)
#     # # dists.save('tax_matrix.pkl')
#     # dists.save_parallel('tax_matrix_p2.lzma', ('lzma', 3))
#
#     o = load_matrix_parallel('tax_matrix_p2.lzma')
#
#     print(o.matrix)
#
#     with open("D:/edwards2016/runs/run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim/rank/rank-test.json", 'r') as preds:
#         preds_d = json.load(preds)
#     pred_json = Path('predictions/blastn.json')
#
#     with pred_json.open() as pj:
#         preds = json.load(pj)
#
#     print(taxonomic_accordance(o,
#                                preds,
#                                virus_dict))
#
#     total_end = timer()
#     total_runtime = total_end - total_start
#     print(f"{Fore.GREEN} Total elapsed time: {total_runtime:.6f} seconds")
