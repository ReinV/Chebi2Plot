#!/usr/bin/python

import networkx
import obonet
import csv
import os.path
import urllib
import json
import time

def return_latest_ontology():
    '''
    This function imports the latest updated version of the ChEBI ontology, and returns the version number and ontology.
    '''

    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo'
    # file = open('files/chebi.obo', encoding = 'utf8')
    print('Importing latest ontology...')
    graph = obonet.read_obo(url)
    # file.close()
    version = graph.graph['data-version']
    return version, graph

def return_current_version():
    '''
    This function opens the ChEBI files with id's and names, and returns the version number used to update this file.
    '''
    file = open('files/ontology_version.txt', 'r')
    version = file.read()
    return version

def return_archived_ontology(version):
    '''
    This function returns an archived ontology based on the version number.
    '''
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel' + version + '/ontology/chebi.obo'
    graph = obonet.read_obo(url)
    # id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    return graph

def show_updates(graph_new, graph_old):
    '''
    This function compares two ontologies and returnes the difference in nodes and edges (chemicals and relations).
    '''
    total_nodes = len(graph_new)
    total_edges = graph_new.number_of_edges()
    difference_nodes = total_nodes - len(graph_old)
    difference_edges = total_edges - graph_old.number_of_edges()
    message = 'Latest ChEBI ontology contains %d new chemicals of %d total chemicals...\nand %d new relations of %d total relations' % (difference_nodes, total_nodes, difference_edges, total_edges)
    return message

def get_mass(node, graph):
    '''
    This function retrieves the mass of a molecule from the ontology.
    '''
    mass = "-"
    try:
        for string in graph.node[node]['property_value']:
            if 'mass' in string and 'monoisotopicmass' not in string:
                mass = string.split('\"')[1] # wat te doen met een 0 ?
    except:
        pass
    return mass

def get_smile(node, graph):
    '''
    This function retrieves Smiles from the ontology.
    '''
    smile = ''
    try:
        for value in graph.node[node]['property_value']:
            if 'smile' in value:
                smile = value.split('\"')[1]
    except:
        pass
    return smile

def get_relations(nodes, graph):
    '''
    This function recieves a node in the graph, and for that node, it retrieves every parent with 'is_a' relationships.
    It returns a dictionary with those parents.
    '''
    parent_to_key = dict()

    for node in nodes:
        for child, parent, key in graph.out_edges(node, keys=True):
            if key == 'is_a':
                try:
                    parent_to_key[parent]
                except:
                    parent_to_key[parent] = key

    return parent_to_key

def get_superterms(id, graph):
    '''
    This function recieves a ChEBI ID and the latest ChEBI ontology.
    It returns a list of parent ChEBI ids with 'is_a' relationship, up to the root.
    '''
    list_relations = []
    nodes = [id]
    end = False

    while end == False:
        parent_to_key = get_relations(nodes, graph) # get the 'is a' relationships
        if len(parent_to_key) == 0: #if there are no 'is a' relationships, end the search
            end = True
        else:
            nodes = []
            for parent in parent_to_key.keys():
                nodes.append(parent)
                new_id = parent.split(":")[1]
                list_relations.append(new_id)

    return list_relations

def update_version_number(number):
    '''
    This function updates the ontology version text file with the version number of the ontology by which the files have been updated.
    '''
    file = open('files/ontology_version.txt', 'w')
    file.write(number)
    return

def get_new_names(graph, file):
    '''
    This function recieves the latest ChEBI ontolog, and the file containing the ChEBI names.
    It returns a dictionary with names of new molecules not yet stored in the file.
    '''
    id_to_info = read_file(file)
    id_to_name = dict()
    for key, data in graph.nodes(data=True):
        id = key.split(":")[1]
        try:
            id_to_info[id]
        except:
            name = data.get('name')
            id_to_name[id] = name
    # Mapping from term ID to name
    # id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    return id_to_name

def get_new_smiles(graph, file):
    '''
    This function recieves the latest ChEBI ontology, and the file containing the smiles.
    It returns a dictionary with smiles of new molecules not yet stored in the file.
    '''
    id_to_info= read_file(file)
    id_to_smile = dict()
    for key in graph.nodes():
        id = key.split(":")[1]
        try:
            id_to_info[id]
        except:
            smile = get_smile(key, graph)
            if smile != '':
                id_to_smile[id] = smile

    return id_to_smile

def get_new_predictions(id_to_smile):
    '''
    This function recieves a dictionary with new smiles. It returns two dictionaries with predicted logP and logS values using the smiles.
    Smiles are put in a string, and added to the url. A task is created which is applying the model AlogPS3.0.
    After the model is done, the output is downloaded using the task_id.
    '''
    print('Getting new predictions from OCHEM...')
    MODEL_ID = 4 # AlogPS3.0 model id
    timeout = 100
    sleep_time = 60

    string_of_smiles = ''
    for smile in id_to_smile.values():
        string_of_smiles += str(smile) + '$$$$'

    url = 'https://ochem.eu/modelservice/postModel.do?'
    values = {'modelId': MODEL_ID, 'mol': string_of_smiles}

    # Start task, apply model
    connection = False
    while connection == False:
        try:
            data = urllib.parse.urlencode(values)
            data = data.encode('ascii')
            response = urllib.request.urlopen(url, data, timeout=timeout)
            json_data = json.loads(response.read())
            connection = True
        except:
            print('connection failed, trying again in %d seconds' % sleep_time)
            time.sleep(sleep_time) # in seconds

    print('OCHEM taskId: %d' % json_data['taskId'])

    # After model is done (while loop), results will be donwloaded
    task_id = json_data["taskId"]
    status = 'pending'
    request_url = 'http://ochem.eu/modelservice/fetchModel.do?taskId='+str(task_id)
    while status == 'pending':
        time.sleep(sleep_time)
        try:
            response = urllib.request.urlopen(request_url, timeout=timeout)
            task_output = json.loads(response.read())
            status = task_output['status']
        except:
            print('applying model, trying again in %d seconds' % sleep_time)

        print('task status: %s' % status)

    count = 0
    exceptions = 0
    id_to_logP = dict()
    id_to_logS = dict()

    # Predicted values will be extracted from the json data and put in dictionaries
    for output_list in task_output['predictions']:
        try:
            predictions = output_list['predictions']
            logP = predictions[0]['value']
            logS = predictions[1]['value']
        except:
            exceptions += 1
            logP = '-'
            logS = '-'
        id = list(id_to_smile.keys())[count]
        id_to_logP[id] = logP
        id_to_logS[id] = logS
        count += 1
    print('Got errors for %d of %d predictions in total' % (exceptions, count))

    return id_to_logP, id_to_logS

def read_file(file):
    '''
    This function reads a file and returns a dictionary of the CHEBI ID's.
    If the file does not exists, the file is made and an empty dictionary is returned.
    '''
    id_to_info = dict()

    if os.path.exists(file):
        f = open(file, 'r')
        lines = f.readlines()

        for line in lines:
            line_to_list = line.split('\t')
            id = line_to_list[0]
            info = line_to_list[1].strip() # nodig?
            id_to_info[id] = info
    else:
        f = open(file, 'w') # make file

    return id_to_info

def rewrite_file(graph, file):
    '''
    This function recieves the latest ChEBI ontology, and the file that will be overwritten.
    Depending on the file argument, it will get the corresponding information (mass, superterms), and writes this to the file.
    '''

    with open(file, 'w', newline='', encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter = '\t')
        for key in graph.nodes():
            id = key.split(":")[1]
            if file == 'files/ChEBI2Mass.tsv':
                info = get_mass(key, graph)
            elif file == 'files/ChEBI2Superterms.tsv':
                info = get_superterms(key, graph)
            writer.writerow([id, info])
    print('%s updated' % file)

    return

def update_file(id_to_info, file):
    '''
    This function recieves a dictionary and a file.
    Information in the dictionary will be added to that file.
    '''
    with open(file, 'a', newline='', encoding="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter = '\t')
        for id in id_to_info:
            info = id_to_info[id]
            writer.writerow([id, info])
    print('%s updated' % file)
    return
    
def main():
    files = ['files/ChEBI2Names.tsv','files/ChEBI2Smiles.tsv', 'files/ChEBI2Superterms.tsv', 'files/ChEBI2Mass.tsv', 'files/ChEBI2logP.tsv', 'files/ChEBI2logS.tsv']
    current_version = return_current_version()
    latest_version, graph = return_latest_ontology() # graph = ontology
    print('ChEBI version used to update files: %s' % current_version)
    print('ChEBI latest version: %s' % latest_version)

    if current_version == latest_version:
        print('Files are up-to-date')
    else:
        print('Files need updating ...')
        graph_old = return_archived_ontology(current_version)
        updates = show_updates(graph, graph_old)
        print(updates)

        id_to_name = get_new_names(graph, file='files/ChEBI2Names.tsv')
        id_to_smile = get_new_smiles(graph, file='files/ChEBI2Smiles.tsv')
        id_to_logP, id_to_logS = get_new_predictions(id_to_smile)

        update_file(id_to_name, file='files/ChEBI2Names.tsv')
        update_file(id_to_smile, file='files/ChEBI2Smiles.tsv')
        update_file(id_to_logP, file='files/ChEBI2logP.tsv')
        update_file(id_to_logS, file='files/ChEBI2logS.tsv')

        rewrite_file(graph, file='files/ChEBI2Mass.tsv')
        rewrite_file(graph, file='files/ChEBI2Superterms.tsv')

        update_version_number(latest_version)

if __name__ == '__main__':
    main()
