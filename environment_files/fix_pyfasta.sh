#!/bin/bash

PYFASTA_DIR=$(conda env list | awk '/utils/ {print $NF "/lib/python3.13/site-packages/pyfasta/"}')

sed -i -e 's/from fasta/from pyfasta.fasta/' -e 's/from records/from pyfasta.records/' -e 's/from split_fasta/from pyfasta.split_fasta/' $PYFASTA_DIR/__init__.py &&
sed -i 's/from collections/from collections.abc/' $PYFASTA_DIR/fasta.py &&
sed -i 's/from records/from pyfasta.records/' $PYFASTA_DIR/fasta.py &&
sed -i 's/import cPickle/import _pickle as cPickle/;s/tostring/tobytes/' $PYFASTA_DIR/records.py &&
sed -i 's/from cStringIO import StringIO/from io import StringIO/' $PYFASTA_DIR/split_fasta.py
