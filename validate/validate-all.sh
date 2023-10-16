#!/bin/bash

set -e

idfile="$1"
outfile="${2:-result.tsv}"

while IFS= read -r line; do
        read afid pdbid <<< "$line"
		echo $afid $pdbid
        ../build/alphafill validate --quiet --config ../alphafill.conf --af-id $afid --pdb-id $pdbid >> $outfile
done < "$idfile"
