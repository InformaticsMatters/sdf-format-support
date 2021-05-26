# Tests

## success/1
A copy of the file in the min-apps-data-tier repo

## Example 1.1: SDF import - UUID added

Generate file with UUID (default).

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=test.sdf
    $ export DATASET_OUTPUT_FORMAT=
    $ export DATASET_EXTRA_VARIABLES=
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up

## Example 1.2: SDF import - No UUID 

Generate file with UUID not overidden. 
This will produce error if the input file has mol blocks that do not already have a UUID in the first line.  

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=test.sdf
    $ export DATASET_OUTPUT_FORMAT=
    $ export DATASET_EXTRA_VARIABLES='generate_uuid=False'
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up

## Example 1.3: SDF import - UUID added and Name rewritten

Generate file with UUID and any existing molecule name written to a new parameter.

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=poses.sdf
    $ export DATASET_OUTPUT_FORMAT=
    $ export DATASET_EXTRA_VARIABLES='existing_title_fieldname=Smiles'
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up

## Example 1.4: SDF import - UUID added ad Name not rewritten

Generate file no UUID and any existing molecule name written to a new parameter - invalid combination.

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=poses.sdf
    $ export DATASET_OUTPUT_FORMAT=
    $ export DATASET_EXTRA_VARIABLES='existing_title_fieldname=Smiles&generate_uuid=False'
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up


## Example 2.1: SDF convert to json

Here we covert our test/success/**1** dataset to a JSON file...

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=test.sdf
    $ export DATASET_OUTPUT_FILENAME=test.json
    $ export DATASET_OUTPUT_FORMAT=application/x-squonk-dataset-molecule-v2+json
    $ export DATASET_EXTRA_VARIABLES=
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up

## Example 2.2: SDF convert to json schema

Here we covert our test/success/**1** dataset to a JSON schema file...

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=test.sdf
    $ export DATASET_OUTPUT_FILENAME=test.schema.json
    $ export DATASET_OUTPUT_FORMAT=application/schema+json
    $ export DATASET_EXTRA_VARIABLES=
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up
