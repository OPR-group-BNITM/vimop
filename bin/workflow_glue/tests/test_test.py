# This file includes or is derived from code developed by Oxford Nanopore Technologies Plc.
# Copyright (c) 2020-, Oxford Nanopore Technologies Plc. All rights reserved. 
# This file is subject to the terms described in the LICENSE file in the root of this repository.

"""A dummy test."""

import argparse
from workflow_glue import report


def test():
    """Just showing that we can import using the workflow-glue."""
    assert isinstance(report.argparser(), argparse.ArgumentParser)
