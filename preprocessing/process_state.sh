#!/bin/bash

# Usage: ./process_state.sh {State Abbrev.} {2-digit State ID} {Full Name (capitalized, underscores for spaces)} {District Type} {Geo. Level}
# ex. ./process_state.sh RI CD tract ./data

set -e

# Check if the number of arguments is correct.
if [ $# -ne 4 ]; then
	echo "Usage: $0 <STATE_ABBR> <DISTRICT_TYPE> <GEO_LEVEL> <DATA_DIR>"
	exit 1
fi

GREEN='\033[0;32m'
NOCOLOR='\033[0m'

STATE_ABBR=${1}
STATE_ID=$(python3 -c "from us import states; print(states.${STATE_ABBR}.fips)")
STATE_NAME=$(python3 -c "from us import states; print('_'.join(states.${STATE_ABBR}.name.split()))")

# CD = Congressional Districts
# SS = State Senate
# SH = State House
# BG = Block Group
DISTRICT=${2}

# tabblock = block level data
# bg = block group level data
# tract = tract level data
LEVEL=${3}

DATA_DIR=${4}
mkdir -p ${DATA_DIR}

# Specify if all racial groups should be included in the processed shapefile
DETAILED=false

COUNTIES=false

DISTRICTING_METHOD="districts"

VINTAGE="20210608"
# VINTAGE="20200527"
# VINTAGE="20210428_12-2"
DATE_STR="${VINTAGE:0:4}-${VINTAGE:4:2}-${VINTAGE:6:2}"

LEVEL_NAME=$LEVEL
if [ "$LEVEL" = "tabblock" ]
then
	DP_URL="https://assets.nhgis.org/differential-privacy/v${VINTAGE}/nhgis_ppdd_${VINTAGE}_block_${STATE_ABBR}.zip"
	DP_FILE="nhgis_ppdd_${VINTAGE}_block_${STATE_ABBR}.csv"
	LEVEL_NAME="block"
elif [ "$LEVEL" = "bg" ]
then
	DP_URL="https://assets.nhgis.org/differential-privacy/v${VINTAGE}/nhgis_ppdd_${VINTAGE}_blck_grp.zip"
	DP_FILE="nhgis_ppdd_${VINTAGE}_blck_grp.csv"
elif [ "$LEVEL" = "tract" ]
then
	DP_URL="https://assets.nhgis.org/differential-privacy/v${VINTAGE}/nhgis_ppdd_${VINTAGE}_tract.zip"
	DP_FILE="nhgis_ppdd_${VINTAGE}_tract.csv"
fi

TIGER_DIR=${DATA_DIR}/2010-SF1/states/${STATE_ABBR}/TIGER
SAVE_DIR=${DATA_DIR}/merged_shp_${VINTAGE}/${STATE_ABBR}

if [ "$LEVEL" = "tabblock" ]
then
	DP_DIR=${DATA_DIR}/${DATE_STR}/state_data/${STATE_ABBR}
else 
	DP_DIR=${DATA_DIR}/${DATE_STR}/${LEVEL}
fi

# Download district outlines
echo -e "${GREEN}Downloading district data${NOCOLOR}"
if [ "$DISTRICT" = "CD" ]
then
	CD_URL="https://www2.census.gov/geo/tiger/TIGER2016/CD/tl_2016_us_cd115.zip"
	CD_DIR=${DATA_DIR}/2010-SF1/districts
	CD_FILE=${DATA_DIR}/2010-SF1/districts/tl_2016_us_cd115.zip
	BEF=${DATA_DIR}/2010-SF1/${STATE_ID}_${STATE_ABBR}_CD113.txt
elif [ "$DISTRICT" = "SH" ]
then
	CD_URL="https://www2.census.gov/geo/tiger/TIGER2014/SLDL/tl_2014_${STATE_ID}_sldl.zip"
	CD_DIR=${DATA_DIR}/2010-SF1/districts/state_house
	CD_FILE=${DATA_DIR}/2010-SF1/districts/state_house/tl_2014_${STATE_ID}_sldl.zip
	BEF=${DATA_DIR}/2010-SF1/state_house/${STATE_ID}_${STATE_ABBR}_SLDL.txt
elif [ "$DISTRICT" = "SS" ]
then
	CD_URL="https://www2.census.gov/geo/tiger/TIGER2014/SLDU/tl_2014_${STATE_ID}_sldu.zip"
	CD_DIR=${DATA_DIR}/2010-SF1/districts/state_sen
	CD_FILE=${DATA_DIR}/2010-SF1/districts/state_sen/tl_2014_${STATE_ID}_sldu.zip
	BEF=${DATA_DIR}/2010-SF1/state_sen/${STATE_ID}_${STATE_ABBR}_SLDU.txt
elif [ "$DISTRICT" = "BG" ]
then
	CD_URL="https://www2.census.gov/geo/pvs/tiger2010st/${STATE_ID}_${STATE_NAME}/${STATE_ID}/tl_2010_${STATE_ID}_bg10.zip"
	CD_DIR=${DATA_DIR}/2010-SF1/districts/block_group
	CD_FILE=${DATA_DIR}/2010-SF1/districts/block_group/tl_2010_${STATE_ID}_bg10.zip
	BEF=""
else
	echo "Invalid District Type"
	exit 1
fi
# echo "wget -nc $CD_URL -P $CD_DIR"
wget -nc $CD_URL -P $CD_DIR -nv

mkdir -p $TIGER_DIR
mkdir -p $DP_DIR
mkdir -p ${SAVE_DIR}/${LEVEL_NAME}

# Download TIGER shapefile for state
echo -e "${GREEN}Downloading TIGER data${NOCOLOR}"
wget -nc https://www2.census.gov/geo/pvs/tiger2010st/${STATE_ID}_${STATE_NAME}/${STATE_ID}/tl_2010_${STATE_ID}_${LEVEL}10.zip \
	-P $TIGER_DIR -nv
unzip -n ${TIGER_DIR}/tl_2010_${STATE_ID}_${LEVEL}10.zip -d ${TIGER_DIR}

wget -nc https://www2.census.gov/geo/pvs/tiger2010st/${STATE_ID}_${STATE_NAME}/${STATE_ID}/tl_2010_${STATE_ID}_cousub10.zip \
	-P $TIGER_DIR -nv
unzip -n ${TIGER_DIR}/tl_2010_${STATE_ID}_cousub10.zip -d ${TIGER_DIR}

# Download DP data for state
echo -e "${GREEN}Downloading differential privacy data${NOCOLOR}"
# echo "wget -nc $DP_URL -P $DP_DIR"
wget -nc $DP_URL -P $DP_DIR -nv
unzip -n ${DP_DIR}/*.zip -d ${DP_DIR}


# Process Shapefile
if [ "$DETAILED" = true ]
then
	DETAILED_FLAG="--detailed"
	DETAILED_STRING="_detailed"
	PLACE_FLAG="--place ${TIGER_DIR}/tl_2010_${STATE_ID}_cousub10.shp"
else
	DETAILED_FLAG=""
	DETAILED_STRING=""
	PLACE_FLAG=""
fi

echo -e "${GREEN}Processing shapefile${NOCOLOR}"
python3 /home/crc/proj/topdown_redistricting/Gingles1/preprocessing/make_jl_shapefile.py ${TIGER_DIR}/tl_2010_${STATE_ID}_${LEVEL}10.shp ${CD_FILE} ${DP_DIR}/${DP_FILE} ${SAVE_DIR}/${LEVEL_NAME}/${STATE_ABBR}_${LEVEL_NAME}_${DISTRICT}${DETAILED_STRING}.shp ${LEVEL_NAME} $PLACE_FLAG $DETAILED_FLAG $COUNTIES_FLAG $NOISE_STR --districting_method ${DISTRICTING_METHOD} --bef ${BEF} --force