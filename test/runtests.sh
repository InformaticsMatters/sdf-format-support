#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# Run all the Tests.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Success 1 - Formatter.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# test 1.1: SDF import - UUID added
# -----------------------------------------------------------------------------

export TEST_TYPE=success
export TEST_DIR=1
export DATASET_FILENAME=test.sdf
export DATASET_OUTPUT_FORMAT=
export DATASET_EXTRA_VARIABLES=
rm -rf -f test/${TEST_TYPE}/${TEST_DIR}/output
mkdir -p test/${TEST_TYPE}/${TEST_DIR}/output
IMAGE_NAME=${PWD##*/} docker-compose up

# -----------------------------------------------------------------------------
# test 1.2 - SDF import - No UUID
# -----------------------------------------------------------------------------

export TEST_TYPE=success
export TEST_DIR=1
export DATASET_FILENAME=test_uuid.sdf
export DATASET_OUTPUT_FORMAT=
export DATASET_EXTRA_VARIABLES='generate_uuid=False'
mkdir -p test/${TEST_TYPE}/${TEST_DIR}/output
IMAGE_NAME=${PWD##*/} docker-compose up

# -----------------------------------------------------------------------------
# test 1.3 - SDF import - UUID added and Name rewritten
# -----------------------------------------------------------------------------

export TEST_TYPE=success
export TEST_DIR=1
export DATASET_FILENAME=poses.sdf
export DATASET_OUTPUT_FORMAT=
export DATASET_EXTRA_VARIABLES='existing_title_fieldname=Smiles'
mkdir -p test/${TEST_TYPE}/${TEST_DIR}/output
IMAGE_NAME=${PWD##*/} docker-compose up

# -----------------------------------------------------------------------------
# test 1.4 - SDF import - Invalid UUID in header
# -----------------------------------------------------------------------------

export TEST_TYPE=success
export TEST_DIR=1
export DATASET_FILENAME=test_uuid_invalid.sdf
export DATASET_OUTPUT_FORMAT=
export DATASET_EXTRA_VARIABLES='generate_uuid=False'
mkdir -p test/${TEST_TYPE}/${TEST_DIR}/output
IMAGE_NAME=${PWD##*/} docker-compose up

# -----------------------------------------------------------------------------
# test 1.5: SDF convert to json
# -----------------------------------------------------------------------------

export TEST_TYPE=success
export TEST_DIR=1
export DATASET_FILENAME=test.sdf
export DATASET_OUTPUT_FILENAME=test.json
export DATASET_OUTPUT_FORMAT=application/x-squonk-dataset-molecule-v2+json
export DATASET_EXTRA_VARIABLES=
mkdir -p test/${TEST_TYPE}/${TEST_DIR}/output
IMAGE_NAME=${PWD##*/} docker-compose up

# -----------------------------------------------------------------------------
# test 1.6: SDF convert to json schema
# -----------------------------------------------------------------------------

export TEST_TYPE=success
export TEST_DIR=1
export DATASET_FILENAME=test.sdf
export DATASET_OUTPUT_FILENAME=test.schema.json
export DATASET_OUTPUT_FORMAT=application/schema+json
export DATASET_EXTRA_VARIABLES=
mkdir -p test/${TEST_TYPE}/${TEST_DIR}/output
IMAGE_NAME=${PWD##*/} docker-compose up

# -----------------------------------------------------------------------------
# Failure 1 - Formatter.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Generate file no UUID and any existing molecule name written to a new parameter - invalid combination.
# -----------------------------------------------------------------------------

export TEST_TYPE=failure
export TEST_DIR=1
export DATASET_FILENAME=poses.sdf
export DATASET_OUTPUT_FORMAT=
export DATASET_EXTRA_VARIABLES='existing_title_fieldname=Smiles&generate_uuid=False'
mkdir -p test/${TEST_TYPE}/${TEST_DIR}/output
IMAGE_NAME=${PWD##*/} docker-compose up
