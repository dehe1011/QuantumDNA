import pandas as pd
import json
import yaml
import os
import sys

frontend = __file__[:__file__.rfind('Quantum_DNA_1.0')]+ 'Quantum_DNA_1.0'
if frontend not in sys.path:
    del sys.path[0]
    sys.path.insert(0, frontend)

# ------------------------------------------------- save and load data --------------------------------------------

def my_save(data, metadata, filename, directory = 'stored_data', save_excel = False):
    ''' Saves data and metadata to a json file. The data (without metadata) can additionally be saved as excel file.
    Note: 
        Usually the filepath starts with 'stored_data/...'
    '''
    full_directory = os.path.join(frontend, directory)
    filename += '.json'
    os.makedirs(full_directory, exist_ok=True)
    filepath = os.path.join(full_directory, filename)
    data = (data, metadata)
    with open(filepath, 'w') as f:
        json.dump(data, f)
    if save_excel: 
        df = pd.DataFrame(data)
        df.to_excel('data.xlsx', index=False)

def my_load(filename, load_metadata = False, directory = 'stored_data'):
    ''' Loads the data and the metadata (if wanted) from a json file. 
    Note: 
        Usually the filepath starts with 'stored_data/...'
    '''
    full_directory = os.path.join(frontend, directory)
    filename += '.json'
    os.makedirs(full_directory, exist_ok=True)
    filepath = os.path.join(full_directory, filename)
    with open(filepath, 'r') as f: 
        data = json.load(f)
    if not load_metadata: 
        return data[0]
    else: 
        return data[0], data[1]

# ----------------------------------------- load configuration file with defualt values --------------------------

def get_config(filename = 'config'):
    ''' 
    Load data stored in a .yaml configuration file (e.g., global variables). 
    '''
    filename += '.yaml'
    filepath = os.path.join(frontend, filename)
    with open(filepath, 'r') as file:
        config = yaml.safe_load(file)
    return config

# ---------------------------------------------------- save figures ------------------------------------

def save_fig(fig, filename: str, directory: str = 'stored_data\plots', format: str = '.svg'):
    """
    Save Figures.
    """
    full_directory = os.path.join(frontend, directory)
    os.makedirs(full_directory, exist_ok=True)
    full_filepath = os.path.join(full_directory, filename) + format
    filepath = os.path.join(directory, filename) + format
    if os.path.exists(full_filepath):
        print(f"Figure {filepath} already exists!")
        filename += '_new'
        full_filepath = os.path.join(full_directory, filename) +format
        filepath = os.path.join(directory, filename) + format
    fig.savefig(full_filepath)
    print(f"Figure saved as {filepath}")
