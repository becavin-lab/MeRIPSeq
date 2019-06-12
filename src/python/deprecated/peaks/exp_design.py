import csv
import sys
from m6a.utils import m6a_utils
from collections import defaultdict


#
#   Read exp design organization
#
def read_exp_design(exp_design_name):
    print(exp_design_name)
    file_name = m6a_utils.PATH + 'ExpDesign/target_' + exp_design_name + '.txt'
    with open(file_name, "rU") as exp_design_file:
        exp_design = csv.DictReader(exp_design_file, delimiter='\t')
        listdata_to_biocond = dict()
        for row in exp_design:
            listdata_to_biocond[row.get('label')] = row.get('Cond')
        return listdata_to_biocond


def get_biocond_to_dataset(exp_design_name):
    print(exp_design_name)
    file_name = m6a_utils.PATH + 'ExpDesign/target_' + exp_design_name + '.txt'
    with open(file_name, "rU") as exp_design_file:
        exp_design = csv.DictReader(exp_design_file, delimiter='\t')
        biocond_to_dataset = defaultdict(list)
        for row in exp_design:
            biocond_to_dataset[row.get('Cond')].append(row.get('label'))
        return biocond_to_dataset


def get_data_list(exp_design_name):
    print(exp_design_name)
    file_name = m6a_utils.PATH + 'ExpDesign/target_' + exp_design_name + '.txt'
    with open(file_name, "rU") as exp_design_file:
        exp_design = csv.DictReader(exp_design_file, delimiter='\t')
        listdata = []
        for row in exp_design:
            listdata.append(row.get('label'))
        return listdata


def get_biocond_list(exp_design_name):
    print(exp_design_name)
    bioconds = set()
    file_name = m6a_utils.PATH + 'ExpDesign/target_' + exp_design_name + '.txt'
    with open(file_name, "rU") as exp_design_file:
        exp_design = csv.DictReader(exp_design_file, delimiter='\t')
        for row in exp_design:
            bioconds.add(row.get('Cond'))
        return bioconds