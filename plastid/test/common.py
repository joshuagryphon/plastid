#!/usr/bin/env python
from numpy.testing import suppress_warnings
from plastid.util.services.exceptions import DataWarning

sup_data = suppress_warnings()
sup_data.filter(category=DataWarning)
