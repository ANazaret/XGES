import logging

import sys
from time import time

sys.path.append('..')

from evaluation.benchmarks import simulation, run_xges
from xges.search import XGES

data, dag = simulation(500, 10000, 3)
# setup logger

run_xges(data, extended_search=False, verbose=1)

model = XGES()
logging.basicConfig(level=logging.INFO)
start_time = time()
model.fit(data, extended_search=False)
print("Execution time: ", time() - start_time)
