#!/usr/bin/env python
"""Convert from SDF file to file with a selected mime-type
"""
import json
import uuid
import datetime
import logging
import os
import gzip
import sys
from typing import Dict
from utils.sdf_utils import sdf_get_next_record, is_valid_uuid

_MIME_TYPE_MAP = {'chemical/x-mdl-sdfile': 'sdf',
                  'application/x-squonk-dataset-molecule-v2+json': 'json',
                  'application/schema+json': 'json_schema'}

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
dataset_output_format = os.getenv('DT_DATASET_OUTPUT_FORMAT')
dataset_output_filename = os.getenv('DT_DATASET_OUTPUT_FILENAME')


# Supporting function for json_schema
def is_number(value: str) -> int:
    """"
    Determines the number type for the property value
    :return 0 for text, 1 for float, 2 for integer
    """
    try:
        float_number = float(value)
    except ValueError:
        return 0
    else:
        if float_number.is_integer():
            return 2
        return 1


def check_name_in_properties(properties, prop_types):
    """ check the name in the properties and set type.
    """

    for name in properties:
        if name not in prop_types:
            prop_types[name] = 2
        prop_name_type = is_number(properties[name])
        if prop_name_type < prop_types[name]:
            prop_types[name] = prop_name_type

    return prop_types


class ConvertFile:
    """Class ConvertFile

    Purpose: Converts from an input mime-type to an output mime-type.

    """
    errors: int = 0
    lines: int = 0
    records: int = 0

    def convert(self, from_mime_type: str, to_mime_type: str, infile: str, outfile: str) -> bool:
        """Dispatch method"""
        from_type: str = _MIME_TYPE_MAP.get(from_mime_type)
        to_type: str = _MIME_TYPE_MAP.get(to_mime_type)

        method_name = 'convert_' + str(from_type) + '_to_' + str(to_type)
        # Get the method from 'self'.
        # method = getattr(self, method_name, lambda: "Invalid conversion")
        try:
            method = getattr(self, method_name)
        except AttributeError:
            event_logger.info('Method to support %s not found', to_mime_type)
            return False

        # Call the method as we return it
        return method(infile, outfile)  # pylint: disable=too-many-function-args

    def process_molecules_json(self, infile, outfile):
        """ process molecules in SDF file
        """

        # Loop through molecules
        while True:
            molecule: Dict = {}
            molecule_block, molecule_name, properties = sdf_get_next_record(infile)

            if not molecule_block:
                break

            if self.records:
                outfile.write(',')
            self.records += 1

            # title line
            if molecule_name:
                molecule['name'] = molecule_name

            if is_valid_uuid(molecule_name):
                molecule_uuid = molecule_name
            else:
                molecule_uuid = str(uuid.uuid4())

            molecule['molblock'] = molecule_block

            record = {'uuid': molecule_uuid, 'molecule': molecule,
                      'values': properties}
            json_str = json.dumps(record)
            outfile.write(json_str)

        outfile.write(']')

    def convert_sdf_to_json(self, infile, outfile) -> bool:
        """converts the given SDF file into a Squonk json file.
           Returns True if file successfully converted.
        """

        if infile.endswith('.gz'):
            infile_handle = gzip.open(infile, 'rt')
        else:
            infile_handle = open(infile, 'rt')

        if outfile.endswith('.gz'):
            outfile_handle = gzip.open(outfile, 'wt')
        else:
            outfile_handle = open(outfile, 'wt')

        outfile_handle.write('[')

        try:
            self.process_molecules_json(infile_handle, outfile_handle)
        finally:
            outfile_handle.close()

        if self.errors > 0:
            return False
        return True

    def process_properties_json(self, infile):
        """ process properties in SDF file
        """

        prop_types = {}
        # Loop through molecules
        while True:
            molecule_block, dummy, properties = sdf_get_next_record(infile)

            if not molecule_block:
                break

            self.records += 1
            prop_types = check_name_in_properties(properties, prop_types)

        return prop_types

    def convert_sdf_to_json_schema(self, infile, outfile) -> bool:
        """converts the given SDF file into a Squonk json schema.
           Returns True if file successfully converted.
        """

        # sdf properties for internal json format - used in the json schema on internal conversion.
        _json_properties_sdf = {
            'uuid': {'type': 'string', 'description': 'Unique UUID'},
            'molecule': {'type': 'object', 'properties': {
                'name': {'type': 'string', 'description': 'Molecule name'},
                'molblock': {'type': 'string',
                             'description': 'Molfile format 2D or 3D representation'},
                'isosmiles': {'type': 'string',
                              'description': 'Isomeric SMILES representation'},
                'stdsmiles': {'type': 'string',
                              'description': 'Standardised non-isomeric SMILES representation'}
            }}}
        _lookups = {0: 'string', 1: 'number', 2: 'integer'}

        # Extract the properties from the file
        if infile.endswith('.gz'):
            infile_handle = gzip.open(infile, 'rt')
        else:
            infile_handle = open(infile, 'rt')

        try:
            prop_types = self.process_properties_json(infile_handle)
        finally:
            infile_handle.close()

        # Add the properties from the file to the schema
        values = {}
        _json_properties_sdf['values'] = {'type': 'object', 'properties': values}
        for prop_name in prop_types:
            prop_type = prop_types[prop_name]
            values[prop_name] = {'type': _lookups[prop_type]}

        schema_sdf = {'$schema': 'http://json-schema.org/draft/2019-09/schema#',
                      'description': 'Automatically created from ' + infile + ' on '
                                     + str(datetime.datetime.now()),
                      'properties': _json_properties_sdf}

        json_str: str = json.dumps(schema_sdf)

        if outfile.endswith('.gz'):
            outfile_handle = gzip.open(outfile, 'w')
        else:
            outfile_handle = open(outfile, 'w')

        outfile_handle.write(json_str)

        if self.errors > 0:
            return False
        return True


if __name__ == "__main__":
    # Say Hello
    basic_logger.info('sdf-format-support')

    # Display environment variables
    basic_logger.info('DT_DATASET_FILENAME=%s', dataset_filename)
    basic_logger.info('DT_DATASET_INPUT_PATH=%s', dataset_input_path)
    basic_logger.info('DT_DATASET_OUTPUT_PATH=%s', dataset_output_path)
    basic_logger.info('DT_DATASET_OUTPUT_FILENAME=%s', dataset_output_filename)
    basic_logger.info('DT_DATASET_OUTPUT_FORMAT=%s', dataset_output_format)

    # format input file path and output file path.
    input_file: str = os.path.join(dataset_input_path, dataset_filename)
    output_file: str = os.path.join(dataset_output_path, dataset_output_filename)

    # if no input file then raise error and exit
    if not os.path.isfile(input_file):
        event_logger.info('File %s is not present', input_file)
        sys.exit(1)

    basic_logger.info('SDF Converter')
    event_logger.info('Processing %s...', input_file)

    converter = ConvertFile()

    processed: bool = converter.convert('chemical/x-mdl-sdfile', dataset_output_format, input_file,
                                        output_file)

    if processed:
        basic_logger.info('SDF Converter finished successfully')
        basic_logger.info('records processed=%s', converter.records)
        basic_logger.info('errors=%s', converter.errors)
    else:
        basic_logger.info('SDF Converter failed')
        sys.exit(1)
