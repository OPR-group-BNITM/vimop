#!/usr/bin/env python

# This file includes or is derived from code developed by Oxford Nanopore Technologies Plc.
# Copyright (c) 2020-, Oxford Nanopore Technologies Plc. All rights reserved. 
# This file is subject to the terms described in the LICENSE file in the root of this repository.

"""Map a basecalling model to a Medaka model.

An unknown basecalling model or basecalling model without an appropriate
Medaka model will explode with a large stderr message and exit non-zero.
A happy basecalling model will print a Medaka model to stdout and exit 0.
"""

# Delegating this to a Python script seems overkill but allows us to
# expand to more complex logic trivially in future.
# Plus I don't want to write this in Groovy right now.


import csv
import os
import sys
import textwrap

from .util import wf_parser  # noqa: ABS101


def exit_obvious_error(header, error_msg, advice, args, width=80):
    """Write an obvious looking error to stderr and quit."""
    line = ("-" * width) + '\n'
    msg = (
        f"The input basecaller configuration '{args.basecaller_cfg}' does not "
        "have a suitable Medaka model "
    )
    sys.stderr.write(line)
    sys.stderr.write(f"[CRITICAL ERROR] {header}\n\n")
    sys.stderr.write(textwrap.fill(msg + error_msg, width))
    sys.stderr.write('\n\n')
    sys.stderr.write(textwrap.fill(advice, width))
    sys.stderr.write('\n')
    sys.stderr.write(line)
    sys.exit(os.EX_DATAERR)


def main(args):
    """Entry point."""
    with open(args.model_tsv) as tsv:
        for row in csv.DictReader(tsv, delimiter='\t'):
            if row["basecall_model_name"] == args.basecaller_cfg:
                model_type = str(args.model_type)
                model = row[model_type]
                reason = row["medaka_nomodel_reason"]
                if model == "-" or model == "":
                    # Basecalling model valid but no Medaka model
                    exit_obvious_error(
                        header="No appropriate Medaka model.",
                        error_msg=f"because {reason}.",
                        advice=(
                            "It is not possible to run Medaka"
                            "with this data.\n"),
                        args=args
                        )
                    break  # exit before here but this keeps my intention obvious
                else:
                    # Good model found
                    sys.stdout.write(model)
                    break
        else:
            # No model found (loop not broken)
            exit_obvious_error(
                header="Unknown basecaller configuration.",
                error_msg=(
                    "because the basecaller configuration has not been recognised."
                ),
                advice=(
                    "Check your --basecaller_cfg has been provided correctly. "
                ),
            )


def argparser():
    """Argument parser for entry point."""
    parser = wf_parser("resolve_medaka_model")
    parser.add_argument("model_tsv")
    parser.add_argument("basecaller_cfg")
    parser.add_argument("model_type")
    return parser
