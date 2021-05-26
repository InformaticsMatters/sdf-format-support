#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# Run all the Tests.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# test 1.1: SDF import - UUID added
# -----------------------------------------------------------------------------

export DATASET_NAME=1
export DATASET_FILENAME=test.sdf
export DATASET_OUTPUT_FORMAT=
export DATASET_EXTRA_VARIABLES=
mkdir -p test/success/${DATASET_NAME}/output
IMAGE_NAME=${PWD##*/} docker-compose up

# -----------------------------------------------------------------------------
# test 1.2 - SDF import - No UUID
# -----------------------------------------------------------------------------

export DATASET_NAME=1
export DATASET_FILENAME=test.sdf
export DATASET_OUTPUT_FORMAT=
export DATASET_EXTRA_VARIABLES='generate_uuid=False'
mkdir -p test/success/${DATASET_NAME}/output
IMAGE_NAME=${PWD##*/} docker-compose up

# -----------------------------------------------------------------------------
# test 1.3 - SDF import - UUID added and Name rewritten
# -----------------------------------------------------------------------------

export DATASET_NAME=1
export DATASET_FILENAME=poses.sdf
export DATASET_OUTPUT_FORMAT=
export DATASET_EXTRA_VARIABLES='existing_title_fieldname=Smiles'
mkdir -p test/success/${DATASET_NAME}/output
IMAGE_NAME=${PWD##*/} docker-compose up
