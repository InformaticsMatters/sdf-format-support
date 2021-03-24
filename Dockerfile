# Base your contaier on something of use and,
# by convention, launch your run-time logic
# from a 'docker-entrypoint.sh' shell-script.
#
# Note: Regardeless of the user you have avalable at build-time,
#       the container's actual user ID and group ID will be assigned by the
#       Data Tier Manager ... so plan accordingly.

# Base image RDKit
ARG from_image=informaticsmatters/rdkit-python3-debian:Release_2019_09
FROM ${from_image}

# All formatter images MUST place their
# implementations (and expect to start) in /home/format-support
WORKDIR /home/format-support

# Add python libraries to environment
COPY requirements.txt ./
USER root
RUN python -m pip install --upgrade pip && \
    pip install -r requirements.txt

# Copy source code in
COPY source/ ./source/

# Inject the entrypoint,
# making sure anyone can read and execute it.
COPY docker-entrypoint.sh ./
RUN chmod a+rx ./*.sh

# Start the formattter
CMD ./docker-entrypoint.sh
