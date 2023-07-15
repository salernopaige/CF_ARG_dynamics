#!/bin/bash

# Clean out old databases
# Perform all in the CARD directory that you want
rgi clean --local
echo "Old CARD databases cleaned out"

# Download the most up to date CARD Database
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
echo "Updated CARD database downloaded"

# Load the newest database to the local directory
rgi load --card_json card.json --local
echo "CARD loaded to working directory"

# Annotate CARD database and add to working directory -- make sure to adjust version
rgi card_annotation -i card.json > card_annotation.log 2>&1
rgi load -i card.json --card_annotation card_database_v3.2.7.fasta --local
echo "Annotated db added to working directory"

# Download the latest wildcard database
wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
mkdir -p wildcard
tar -xjf wildcard_data.tar.bz2 -C wildcard
gunzip wildcard/*.gz
echo "Wildcard database downloaded"

# Annotate the wildcard database and load it -- make sure to adjust version
rgi wildcard_annotation -i wildcard --card_json card.json -v 3.2.7 > wildcard_annotation.log 2>&1
rgi load --wildcard_annotation wildcard_database_v3.2.7.fasta --card_json card.json --wildcard_index wildcard/index-for-model-sequences.txt --card_annotation card_database_v3.2.7.fasta --local
 echo "Wildcard database annotated and loaded"
