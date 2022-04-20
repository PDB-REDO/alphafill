#!/bin/bash

idfile="$1"

while IFS= read -r line; do
        read afid pdbid <<< "$line"
        validate-fill --quiet --config ../alphafill.conf $afid $pdbid
done < "$idfile"
