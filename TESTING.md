# Testing
You must one subdirectory in `test/success`, and a dataset file in it's `input`
directory that your formatter can process successfully. The template's
`dummy.txt` file in `test/success/1` illustrates the expected file structure.

Remember, to allow for possible future automated testing, each test dataset
should be in a separate `input` directory, and you must have at least one
successful test and (ideally) some failure tests: -

-   `test/success/1/input/<dataset file>`
-   `test/failure/1/input/<dataset file>`

Following this structure means your test data can be easily applied to your
container for functional testing using the `docker-compose.yaml` file provided.

You'll need to have installed `docker-compose` - ideally from within a
Python virtual environment: -

    $ python -m venv ~/.venv/format-support
    $ source ~/.venv/format-support/bin/activate

    $ pip install --upgrade pip
    $ pip install -r build-requirements.txt

...before building your formatter image using `docker-compose`: -

    $ IMAGE_NAME=${PWD##*/} docker-compose build

...and then, simulating the provision of the output directory that would
normally be created by the DataTier Manager, we can run a specific test...

## Example 1.1: SDF import - UUID added

Generate file with UUID (default).

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=test.sdf.gz
    $ export DATASET_OUTPUT_FORMAT=
    $ export DATASET_EXTRA_VARIABLES=
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up

## Example 1.2: SDF import - No UUID 

Generate file with UUID not overidden. 
This will produce error if the input file has mol blocks that do not already have a UUID in the first line.  

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=test.sdf.gz
    $ export DATASET_OUTPUT_FORMAT=
    $ export DATASET_EXTRA_VARIABLES='generate_uuid=False'
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up

## Example 1.3: SDF import - UUID added and Name rewritten

Generate file with UUID and any existing molecule name written to a new parameter.

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=poses.sdf.gz
    $ export DATASET_OUTPUT_FORMAT=
    $ export DATASET_EXTRA_VARIABLES='existing_title_fieldname=Smiles'
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up

## Example 1.4: SDF import - UUID added ad Name not rewtitten

Generate file no UUID and any existing molecule name written to a new parameter - invalid combination.

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=poses.sdf.gz
    $ export DATASET_OUTPUT_FORMAT=
    $ export DATASET_EXTRA_VARIABLES='existing_title_fieldname=Smiles&generate_uuid=False'
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up


## Example 2.1: SDF convert to json

Here we covert our test/success/**1** dataset to a JSON file...

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=test.sdf.gz
    $ export DATASET_OUTPUT_FILENAME=test.json
    $ export DATASET_OUTPUT_FORMAT=application/x-squonk-dataset-molecule-v2+json
    $ export DATASET_EXTRA_VARIABLES=
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up

## Example 2.2: SDF convert to json schema

Here we covert our test/success/**1** dataset to a JSON schema file...

    $ export DATASET_NAME=1
    $ export DATASET_FILENAME=test.sdf.gz
    $ export DATASET_OUTPUT_FILENAME=test.schema.json
    $ export DATASET_OUTPUT_FORMAT=application/schema+json
    $ export DATASET_EXTRA_VARIABLES=
    $ mkdir -p test/success/${DATASET_NAME}/output
    $ IMAGE_NAME=${PWD##*/} docker-compose up
