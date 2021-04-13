"""Process an SDF file and a csv file containing data to load into database
"""
import logging
import os
import gzip
import traceback
import csv
import sys
import uuid

from typing import Dict

from rdkit import Chem, RDLogger
from standardize_molecule import standardize_to_noniso_smiles
from utils.sdf_utils import sdf_get_next_record, sdf_write_record, sdf_add_property, is_valid_uuid

# The columns *every* standard file is expected to contain.
# Use UPPER_CASE.
# All standard files must start with these columns.
_OUTPUT_COLUMNS = ['osmiles', 'smiles', 'inchis', 'inchik', 'hac', 'molecule-uuid']

# Two loggers - one for basic logging, one for events.
basic_logger = logging.getLogger('basic')
basic_logger.setLevel(logging.INFO)
basic_handler = logging.StreamHandler()
basic_formatter = logging.Formatter('%(asctime)s # %(levelname)s %(message)s')
basic_handler.setFormatter(basic_formatter)
basic_logger.addHandler(basic_handler)

event_logger = logging.getLogger('event')
event_logger.setLevel(logging.INFO)
event_handler = logging.StreamHandler()
event_formatter = logging.Formatter('%(asctime)s # %(levelname)s -EVENT- %(message)s')
event_handler.setFormatter(event_formatter)
event_logger.addHandler(event_handler)

# Get and display the environment material
# (guaranteed to be provided)
# using the basic (non-event) logger
dataset_filename = os.getenv('DT_DATASET_FILENAME')
dataset_input_path = os.getenv('DT_DATASET_INPUT_PATH')
dataset_output_path = os.getenv('DT_DATASET_OUTPUT_PATH')
dataset_extra_variables = os.getenv('DT_DATASET_EXTRA_VARIABLES')
processing_vars: Dict = {}

def get_processing_variables():
    """Validate and return extra variables provided by user in API call
    The assumption is that this will be a text block of format name=value
    separated by line feeds.

    :returns: process_vars - dictionary of processing variables.
    """
    process_vars = {}

    # Set defaults
    _valid_params = ['generate_uuid', 'existing_title_fieldname']
    process_vars['generate_uuid'] = True
    process_vars['existing_title_fieldname'] = ''

    # Split into list of pairs
    params = dataset_extra_variables.split('&')
    for row in params:
        param = row.split('=')
        if param[0].lower() in _valid_params:
            process_vars[param[0]] = param[1]

    if isinstance(process_vars['generate_uuid'], str) and \
            process_vars['generate_uuid'].lower() in ['false']:
        process_vars['generate_uuid'] = False

    if process_vars['existing_title_fieldname'] and not process_vars['generate_uuid']:
        event_logger.info('Invalid combination of parameters')
        sys.exit(1)

    return process_vars


def write_output_sdf(output_sdf_file, molecule_block, molecule_name, properties,
                     rec_nr):
    """Write the given record to the output file

    Some of this may end up in std_utils..

    :param processing_vars: Dict
    :param output_sdf_file: File Object
    :param molecule_block:
    :param molecule_name:
    :param properties:
    :param rec_nr: to be inserted into file
    :returns:
    """

    # User has named a property name to receive the existing molecule name that will be
    # overridden by the UUID.
    if processing_vars['existing_title_fieldname'] and molecule_name:
        properties = sdf_add_property(properties, processing_vars['existing_title_fieldname'],
                                      [molecule_name])

    molecule_uuid = str(uuid.uuid4())

    sdf_write_record(output_sdf_file, molecule_block, molecule_uuid, properties, rec_nr)

    return molecule_uuid


def process_file(output_writer, input_sdf_file, output_sdf_file):
    """Process the given dataset and process the molecule
    information, writing it as csv-separated fields to the output.

    As we load the molecule we 'standardise' the SMILES and
    add inchi information.

    :param output_writer: DictWriter instance of csv output file
    :param input_sdf_file: SDF File to process
    :param output_sdf_file: SDF File to write to if required.
    :returns: The number of items processed and the number of failures
    """

    num_processed = 0
    num_failed = 0
    num_mols = 0

    while True:

        molecule_block, molecule_name, properties = sdf_get_next_record(input_sdf_file)

        if not molecule_block:
            if processing_vars['generate_uuid']:
                # End and Close the SDF file if there are no more molecules.
                output_sdf_file.close()
            break

        # If something exists in the molecule_block, then a record has been found,
        # otherwise end of file has been reached.
        num_processed += 1
        mol = Chem.MolFromMolBlock(molecule_block)

        if not mol:
            # RDKit could not handle the record
            num_failed += 1
            event_logger.info('RDkit cannot handle record %s - empty mol', num_processed)
            continue

        try:
            osmiles = Chem.MolToSmiles(mol)
            # Standardise the smiles and return a standard molecule.
            noniso = standardize_to_noniso_smiles(osmiles)

            if not noniso[1]:
                num_failed += 1
                event_logger.info('Record %s failed to standardize in RDKit', num_processed)
                continue

            num_mols += 1
            if processing_vars['generate_uuid']:
                # If we are generating a UUID for the molecules then we need to rewrite
                # the record to a new output SDF file.
                molecule_uuid = write_output_sdf(output_sdf_file,
                                                 molecule_block, molecule_name, properties,
                                                 num_mols)
            else:
                # if we are not generating a UUID then the molecule name must already contain
                # a UUID.
                if is_valid_uuid(molecule_name):
                    molecule_uuid = molecule_name
                else:
                    num_failed += 1
                    event_logger.info('Record %s did not contain a uuid', num_processed)
                    continue

            inchis = Chem.inchi.MolToInchi(noniso[1], '')
            inchik = Chem.inchi.InchiToInchiKey(inchis)

            # Write the standardised data to the csv file
            output_writer.writerow({'osmiles': osmiles, 'smiles': noniso[0],
                             'inchis': inchis, 'inchik': inchik,
                             'hac': noniso[1].GetNumHeavyAtoms(), 'molecule-uuid': molecule_uuid})

        except: # pylint: disable=bare-except
            num_failed += 1
            traceback.print_exc()
            event_logger.info('%s Caused a failure in RDKit', osmiles)
            sys.exit(1)

    return num_processed, num_failed, num_mols


if __name__ == '__main__':
    # Say Hello
    basic_logger.info('sdf-format-support')

    # Display environment variables
    basic_logger.info('DT_DATASET_FILENAME=%s', dataset_filename)
    basic_logger.info('DT_DATASET_INPUT_PATH=%s', dataset_input_path)
    basic_logger.info('DT_DATASET_OUTPUT_PATH=%s', dataset_output_path)
    basic_logger.info('DT_DATASET_EXTRA_VARIABLES=%s', dataset_extra_variables)

    processing_vars = get_processing_variables()
    basic_logger.info('generate_uuid=%s', processing_vars['generate_uuid'])
    basic_logger.info('existing_title_fieldname=%s', processing_vars['existing_title_fieldname'])

    basic_logger.info('SDF Data Loader')

    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Open the file we'll write the standardised data set to.
    loader_filename = os.path.join(dataset_output_path, 'tmploaderfile.csv')
    basic_logger.info('Writing to %s...', loader_filename)

    with open(loader_filename, 'wt') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=_OUTPUT_COLUMNS)
        writer.writeheader()

        input_filename = os.path.join(dataset_input_path, dataset_filename)
        output_filename = os.path.join(dataset_output_path, dataset_filename)

        event_logger.info('Processing %s...', input_filename)
        output_sdf: object = None

        if dataset_filename.endswith('.gz'):
            input_sdf = gzip.open(input_filename, 'rt')
            if processing_vars['generate_uuid']:
                output_sdf = gzip.open(output_filename, 'wt')
        else:
            input_sdf = open(input_filename, 'rt')
            if processing_vars['generate_uuid']:
                output_sdf = open(output_filename, 'wt')

        processed, failed, mols =\
            process_file(writer, input_sdf, output_sdf)

    # Summary
    event_logger.info('{:,} processed molecules'.format(processed))
    basic_logger.info('{:,} molecule failures'.format(failed))
    basic_logger.info('{:,} molecule added'.format(mols))
    basic_logger.info('Process complete')
