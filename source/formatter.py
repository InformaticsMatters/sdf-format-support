"""Process an SDF file into a csv file containing data to load into database"""
import logging
import os
import traceback
import csv
import sys
import uuid
import datetime
import json
import gzip
import shutil

from typing import Dict

from standardize_molecule import standardize_to_noniso_smiles
from data_manager_metadata.metadata import Metadata, FieldsDescriptorAnnotation
from data_manager_metadata.data_tier_api import get_metadata_filenames
from data_manager_metadata.annotation_utils import est_schema_field_type

from rdkit import Chem, RDLogger

from utils.sdf_utils import (sdf_get_next_record,
                             sdf_write_record,
                             sdf_add_property,
                             is_valid_uuid)

# The columns *every* standard file is expected to contain.
# All standard files must start with these columns.
_OUTPUT_COLUMNS = ['smiles', 'inchis', 'inchik', 'hac', 'molecule-uuid',
                   'rec_number']

# Base FieldsDescriptor fields to create an SDF annotation with.
_BASE_FIELD_NAMES = {
    'molecule': {'type': 'object', 'description': 'SDF molecule block',
               'required': True, 'active': True},
    'uuid': {'type': 'string', 'description': 'Unique Identifier',
             'required': True, 'active': True},
    }

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
event_formatter = logging.Formatter\
    ('%(asctime)s # %(levelname)s -EVENT- %(message)s')
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
    separated by '&'.

    :returns: process_vars - dictionary of processing variables.
    """
    process_vars = {}

    # Set defaults
    _valid_params = ['generate_uuid', 'existing_title_fieldname']
    process_vars['generate_uuid'] = True
    process_vars['existing_title_fieldname'] = ''

    if not dataset_extra_variables:
        return process_vars

    # Split into list of pairs
    try:
        params = dataset_extra_variables.split('&')
        for row in params:
            param = row.split('=')
            if param[0].lower() in _valid_params:
                process_vars[param[0]] = param[1]

        if isinstance(process_vars['generate_uuid'], str) and \
                process_vars['generate_uuid'].lower() in ['false']:
            process_vars['generate_uuid'] = False

        if process_vars['existing_title_fieldname'] and \
                not process_vars['generate_uuid']:
            event_logger.error('Invalid combination of parameters')
            sys.exit(1)

        return process_vars

    except:  # pylint: disable=bare-except
        event_logger.error('Problem decoding parameters - please check format')
        sys.exit(1)


# Gzip Supporting functions
def uncompress_file():
    """"
    Uncompress the file into the same location
    Remove .gz from filename
    """
    uncompressed_filename = os.path.splitext(dataset_filename)[0]
    uncompressed_path = os.path.join(dataset_input_path, uncompressed_filename)
    with gzip.open(input_filename, 'rt') as f_in:
        with open(uncompressed_path, 'w') as f_out:
            shutil.copyfileobj(f_in, f_out)

    event_logger.info('Input file uncompressed')

    return uncompressed_path, uncompressed_filename


def compress_file():
    """"
    Recompress the file into the same location
    Remove the uncompressed file
    """
    compressed_path = \
        os.path.join(dataset_output_path, dataset_filename)
    with open(output_filename, 'rb') as f_in:
        with gzip.open(compressed_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    event_logger.info('Output file compressed')
    # Tidy up
    os.remove(output_filename)
    os.remove(input_filename)
    event_logger.info('Temporary files removed')


def check_name_in_fields(properties, fields) -> dict:
    """ check the name in the properties. If the name does not exist
    then add the name and type to the fields dictionary.
    """

    for name in properties:
        if name not in fields:
            fields[name] = est_schema_field_type(properties[name])

    return fields


def _log_progress(num_processed):

    if not num_processed % 50000:
        event_logger.info('%s records processed', num_processed)


def write_output_sdf_error(output_sdf_file, molecule_block, molecule_name,
                           properties, rec_nr):
    """Write the given record to the output file in the case of an invalid
    molecule or smiles

    :param output_sdf_file: File Object
    :param molecule_block:
    :param molecule_name:
    :param properties:
    :param rec_nr: to be inserted into file
    """

    # User has named a property name to receive the existing molecule name
    # that will be overridden by the UUID.
    if processing_vars['existing_title_fieldname'] and molecule_name:
        properties = sdf_add_property\
            (properties, processing_vars['existing_title_fieldname'].lower(),
             [molecule_name])

    molecule_uuid = 'Error'
    sdf_write_record(output_sdf_file, molecule_block, molecule_uuid,
                     properties, rec_nr)


def write_output_sdf(output_sdf_file, molecule_block, molecule_name,
                     properties,
                     rec_nr):
    """Write the given record to the output file

    :param output_sdf_file: File Object
    :param molecule_block:
    :param molecule_name:
    :param properties:
    :param rec_nr: to be inserted into file
    :returns: uuid csv file
    """

    # User has named a property name to receive the existing molecule name
    # that will be overridden by the UUID.
    if processing_vars['existing_title_fieldname'] and molecule_name:
        properties = sdf_add_property\
            (properties, processing_vars['existing_title_fieldname'].lower(),
             [molecule_name])

    molecule_uuid = str(uuid.uuid4())

    sdf_write_record(output_sdf_file, molecule_block,
                     molecule_uuid, properties, rec_nr)

    return molecule_uuid


def process_file(output_writer, input_sdf_file, output_sdf_file):
    """Process the given dataset and process the molecule
    information, writing it as csv-separated fields to the output.

    As we load the molecule we 'standardise' the SMILES and
    add inchi information.

    :param output_writer: DictWriter instance of csv output file
    :param input_sdf_file: SDF File to process
    :param output_sdf_file: SDF File to write to if required.
    :returns: The number of items processed, fails, molecules and field types
    """

    num_processed = 0
    num_failed = 0
    num_mols = 0
    osmiles = ''

    fields = {'molecule': 'object', 'uuid': 'string'}

    while True:

        molecule_block, molecule_name, properties = \
            sdf_get_next_record(input_sdf_file)

        if not molecule_block:
            if processing_vars['generate_uuid']:
                # End and Close the SDF file if there are no more molecules.
                output_sdf_file.close()
            break

        # If something exists in the molecule_block, then a record has been
        # found, otherwise end of file has been reached.
        num_processed += 1
        _log_progress (num_processed)

        mol = Chem.MolFromMolBlock(molecule_block)

        if not mol:
            # RDKit could not handle the record
            num_failed += 1
            if processing_vars['generate_uuid']:
                write_output_sdf_error(output_sdf_file,
                                       molecule_block,
                                       molecule_name,
                                       properties,
                                       num_processed)
            event_logger.info('RDkit cannot handle record %s - empty mol',
                              num_processed)
            continue

        try:
            osmiles = Chem.MolToSmiles(mol)
            # Standardise the smiles and return a standard molecule.
            noniso = standardize_to_noniso_smiles(osmiles)

            if not noniso[1]:
                num_failed += 1
                if processing_vars['generate_uuid']:
                    write_output_sdf_error(output_sdf_file, molecule_block,
                                           molecule_name,
                                           properties, num_processed)
                event_logger.info('Record %s failed to standardize in RDKit',
                                  num_processed)
                continue

            num_mols += 1
            if processing_vars['generate_uuid']:
                # If we are generating a UUID for the molecules then we need
                # to rewrite the record to a new output SDF file.
                molecule_uuid = write_output_sdf(output_sdf_file,
                                                 molecule_block,
                                                 molecule_name,
                                                 properties,
                                                 num_processed)
            else:
                # if we are not generating a UUID then the molecule name must
                # already contain a UUID.
                if is_valid_uuid(molecule_name):
                    molecule_uuid = molecule_name
                else:
                    num_failed += 1
                    event_logger.info('Record %s did not contain a uuid',
                                      num_processed)
                    continue

            inchis = Chem.inchi.MolToInchi(noniso[1], '')

            # Save any new fields found in list to create Fields descriptor
            fields = check_name_in_fields(properties, fields)

            # Write the standardised data to the csv file
            output_writer.writerow({'smiles': noniso[0],
                                    'inchis': inchis,
                                    'inchik':
                                        Chem.inchi.InchiToInchiKey(inchis),
                                    'hac': noniso[1].GetNumHeavyAtoms(),
                                    'molecule-uuid': molecule_uuid,
                                    'rec_number': num_processed})

        except:  # pylint: disable=bare-except
            num_failed += 1
            traceback.print_exc()
            event_logger.error('%s Caused a failure in RDKit', osmiles)
            sys.exit(1)

    return num_processed, num_failed, num_mols, fields


def process_fields_descriptor(fields):
    """ Create the FieldsDescriptor annotation that will be uploaded with the
    dataset.
    If a schema has been included, then the FieldsDescriptor is initialised
    from this and then upserted with the fields present in the file.
    Any fields in the FD that are not in the file are set to inactive
    """
    origin =  'Automatically created from ' + dataset_filename + ' on ' \
              + str(datetime.datetime.utcnow())
    anno_in_desc = ''
    anno_in_fields = {}
    # If a FieldsDescriptor has been generated from an existing file
    # (say it's a new version of an existing file or derived from an
    # existing file), then prime the fields list
    if os.path.isfile(meta_in_filename):
        with open(meta_in_filename, 'rt') as meta_in_file:
            metadata = Metadata(json.load(meta_in_file))
            f_desc = metadata.get_compiled_fields()
            fd_in_desc = f_desc['description']
            fd_in_fields = f_desc['fields']

    if fd_in_fields:
        event_logger.info('Gernerating annotations from existing '
                          'FieldsDescriptor')
    else:
        fd_in_fields=_BASE_FIELD_NAMES
        event_logger.info('Gernerating new FieldsDescriptor')

    fd_new = FieldsDescriptorAnnotation(origin=origin,
                                        description=fd_in_desc,
                                        fields=fd_in_fields)

    # Match old and new fields
    # If field exists in fields and fd_new then ignore
    # If field in fields but not in fd_new then add
    # If field exists in fd_new but not in fields then make inactive.
    existing_fields = fd_new.get_fields()
    for field in fields:
        if field not in existing_fields:
            fd_new.add_field(field, True, fields[field])

    for field in existing_fields:
        if field not in fields:
            fd_new.add_field(field, False)

    # Recreate output and write the list of annotations to it.
    with open(meta_out_filename, "w") as meta_file:
        metadata.add_annotation(fd_new)
        json.dump(metadata.to_dict(), meta_file)
    event_logger.info('FieldsDescriptor generated')


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
    basic_logger.info('existing_title_fieldname=%s',
                      processing_vars['existing_title_fieldname'])

    basic_logger.info('SDF Data Loader')

    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Non-invasive way of allowing gzip files to be sent.
    # If the input file is a gzip, then uncompress for processing
    # Recompress on exit
    input_filename = os.path.join(dataset_input_path, dataset_filename)
    process_filename = dataset_filename
    compress: bool = False
    if dataset_filename.endswith('.gz'):
        compress: bool = True
        input_filename, process_filename = \
            uncompress_file()

    loader_filename = os.path.join(dataset_output_path, 'tmploaderfile.csv')
    meta_filename, dummy = get_metadata_filenames(process_filename)
    meta_in_filename = os.path.join(
        dataset_input_path,
        meta_filename)
    meta_out_filename = os.path.join(
        dataset_output_path,
        meta_filename)

    basic_logger.info('Writing data to %s...', loader_filename)
    basic_logger.info('Looking for current annotation in %s...',
                      meta_in_filename)
    basic_logger.info('Writing annotations to %s...', meta_out_filename)

    # Open the file we'll write the standardised data set to.
    with open(loader_filename, 'wt') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=_OUTPUT_COLUMNS)
        writer.writeheader()

        output_filename = os.path.join(dataset_output_path, process_filename)

        event_logger.info('Processing %s...', input_filename)
        output_sdf: object = None

        input_sdf = open(input_filename, 'rt')
        if processing_vars['generate_uuid']:
            output_sdf = open(output_filename, 'wt')

        processed, failed, mols, file_fields =\
            process_file(writer, input_sdf, output_sdf)

    if compress:
        compress_file()

    process_fields_descriptor(file_fields)

    # Summary
    event_logger.info('{:,} processed molecules'.format(processed))
    basic_logger.info('{:,} molecule failures'.format(failed))
    basic_logger.info('{:,} molecule added'.format(mols))
    basic_logger.info('Process complete')
