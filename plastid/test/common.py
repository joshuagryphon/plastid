#!/usr/bin/env python
from numpy.testing import suppress_warnings
from plastid.util.services.exceptions import DataWarning

#===============================================================================
# Warnings suppression
#
# Use within bodies of test functions as e.g. `with sup_data: foo`
# Do NOT use as function decorators, as these will hose `yield` tests
#===============================================================================

sup_data = suppress_warnings()
sup_data.filter(category=DataWarning)
