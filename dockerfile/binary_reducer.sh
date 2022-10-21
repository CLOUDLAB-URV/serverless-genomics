#!/bin/bash

###
# USAGE: ./binary_reducer.sh {wholw path to the merger script} {number of threads {} {files*}}
###

COMMAND="$1"
THREADS="$2"
INPUTS=( "${@:3}" )

exit 0
if [[ "$COMMAND" == "" || "$THREADS" == "" || ! ( "$THREADS" =~ ^[1-9][0-9]*$ ) || "${#INPUTS[@]}" == 0 ]]; then
        echo "Usage: binary_reducer <command> <num_threads> <input_1> ..." > /dev/stderr
          exit 1
fi

set -e

#echo "$COMMAND"
#echo "$THREADS"
#echo "${#INPUTS[@]}"

OUT_DIR=$(mktemp -d -p .)
if [[ "$OUT_DIR" == "" ]]; then
        echo "ERROR: Could not create temporary directory" > /dev/stderr
          exit 1
fi

LEVEL=0
while (( ${#INPUTS[@]} > 1 )); do
        LEVEL=$(( $LEVEL + 1 ))
        OUTPUTS=()
        #echo ===Level $LEVEL=== > /dev/stderr
        #echo " INPUTS: ${INPUTS[@]}" > /dev/stderr
        while (( ${#INPUTS[@]} > 1 )); do
                LEVEL_THREADS=$(( ${#INPUTS[@]} / 2 ))
                if (( $LEVEL_THREADS > $THREADS )); then
                        LEVEL_THREADS=$THREADS
                fi
                for ((I=1;I<=LEVEL_THREADS;++I)); do
                        OUTPUT="${OUT_DIR}/${LEVEL}_${#OUTPUTS[@]}"
                        $COMMAND "${INPUTS[0]}" "${INPUTS[1]}" > "$OUTPUT" &
                        OUTPUTS+=( "$OUTPUT" )
                        INPUTS=( "${INPUTS[@]:2}" )
                done;
                wait
        done
        if (( ${#INPUTS[@]} == 1 )); then
                # Odd number of inputs
                OUTPUT="${OUT_DIR}/${LEVEL}_${#OUTPUTS[@]}"
                cp "${INPUTS[0]}" "$OUTPUT"
                OUTPUTS+=( "$OUTPUT" )
        fi
        #echo " OUTPUTS: ${OUTPUTS[@]}" > /dev/stderr
        INPUTS=( "${OUTPUTS[@]}" )
done

cat "${INPUTS[0]}"
rm -rf "$OUT_DIR"