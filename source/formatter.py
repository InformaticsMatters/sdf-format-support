import logging
import os
import gzip
import traceback
import csv
from rdkit import Chem, RDLogger
from standardize_molecule import standardize_to_noniso_smiles

# The columns *every* standard file is expected to contain.
# Use UPPER_CASE.
# All standard files must start with these columns.
_OUTPUT_COLUMNS = ['smiles', 'inchis', 'inchik', 'hac']

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


# Now enter the formatting logic...
def process_file(writer, dataset_file):
    """Process the given dataset and process the molecule
    information, writing it as csv-separated fields to the output.

    As we load the molecule we 'standardise' the SMILES and
    add inchi information.

    :param writer: DictWriter instance of csv output file
    :param dataset_file: The (compressed) file to process
    :returns: The number of items processed and the number of failures
    """

    num_processed =0
    num_failed = 0
    num_mols = 0

    event_logger.info('Processing %s...', dataset_file)

    if dataset_file.endswith('.gz'):
        gz = gzip.open(dataset_file)
        supplr = Chem.ForwardSDMolSupplier(gz)
    else:
        supplr = Chem.ForwardSDMolSupplier(dataset_file)

    for mol in supplr:

        if not mol:
            # RDKit could not handle the record
            num_failed += 1
            continue
        try:
            osmiles = Chem.MolToSmiles(mol)
            # Standardise the smiles and return a standard molecule.
            noniso = standardize_to_noniso_smiles(osmiles)

            if not noniso[1]:
                num_failed += 1
                continue

            num_mols += 1
            smiles = noniso[0]
            inchis = Chem.inchi.MolToInchi(noniso[1], '')
            inchik = Chem.inchi.InchiToInchiKey(inchis)
            hac = noniso[1].GetNumHeavyAtoms()

            # Write the standardised data to the csv file
            writer.writerow({'smiles': smiles, 'inchis': inchis, 'inchik': inchik, 'hac': hac})

        except:
            num_failed += 1
            traceback.print_exc()

        # Enough?
        num_processed += 1

    return num_processed, num_failed, num_mols


if __name__ == '__main__':

    # Say Hello
    basic_logger.info('sdf-format-support')

    # Display environment variables
    basic_logger.info('DT_DATASET_FILENAME=%s', dataset_filename)
    basic_logger.info('DT_DATASET_INPUT_PATH=%s', dataset_input_path)
    basic_logger.info('DT_DATASET_OUTPUT_PATH=%s', dataset_output_path)

    basic_logger.info('SDF Data Loader')

    # Suppress basic RDKit logging...
    RDLogger.logger().setLevel(RDLogger.ERROR)

    # Open the file we'll write the standardised data set to.
    # A text, tab-separated file.
    output_filename = os.path.join(dataset_output_path, 'molfile.csv')
    basic_logger.info('Writing to %s...', output_filename)

    num_processed = 0
    with open(output_filename, 'wt') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=_OUTPUT_COLUMNS)
        writer.writeheader()
        num_processed, number_failed, number_mols =\
            process_file(writer, os.path.join(dataset_input_path,
                                              dataset_filename))

    # Summary
    event_logger.info('{:,} processed molecules'.format(num_processed))
    basic_logger.info('{:,} molecule failures'.format(number_failed))
    basic_logger.info('{:,} molecule added'.format(number_mols))
    basic_logger.info('Process complete'.format(number_failed))
