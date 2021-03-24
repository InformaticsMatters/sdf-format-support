#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# The DataTier Dataset 'format support' container entrypoint.
# -----------------------------------------------------------------------------
# Refer to RULES.md for comprehensive format support rules,
# summarised below...
#
# 1. The environment variable `DT_DATASET_FILENAME` will be set to
#    a string value representing the filename of the dataset
# 2. The environment variable `DT_DATASET_INPUT_PATH` will be
#    a directory you can find the file defined in part 1
# 3. The environment variable `DT_DATASET_OUTPUT_PATH` will be
#    a directory you can write to
# 4. The environment variable `DT_DATASET_OUTPUT_FORMAT` will be
#    set to a MIME type if a file format conversion,
#    rather than data processing, is to be performed
# -----------------------------------------------------------------------------

python source/formatter.py
