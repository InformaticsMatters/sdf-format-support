import logging
import os
import time

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

# Say Hello
basic_logger.info('sdf-format-support')

# Get and display the environment material
# (guaranteed to be provided)
# using the basic (non-event) logger
dataset_name = os.getenv('DT_DATASET_NAME')
dataset_file = os.getenv('DT_DATASET_FILE')
dataset_output_path = os.getenv('DT_DATASET_OUTPUT_PATH')

basic_logger.error('DT_DATASET_NAME=%s', dataset_name)
basic_logger.info('DT_DATASET_FILE=%s', dataset_file)
basic_logger.info('DT_DATASET_OUTPUT_PATH=%s', dataset_output_path)

# Now enter the formatting logic...
# Here we use the event logger.
event_logger.info('Progress %d%%', 0)
time.sleep(4.0)
event_logger.info('Progress %d%%', 50)
time.sleep(4.0)
event_logger.info('Progress %d%%', 100)
