"""One J() evaluation in its own process, started by calibrate_athabasca_asyncdds_mpi.py.
Usage: python eval_one.py <result_file> <names_json> <values_json>
SUMMA_MASTER and SUMMA_CONFIG come from the calling worker's environment; the working
directory is the worker's rank folder. Writes the KGE to <result_file>; a fatal model error
ends this process with a non-zero exit code instead."""
import json, sys
from pathlib import Path

res, names, values = Path(sys.argv[1]), json.loads(sys.argv[2]), json.loads(sys.argv[3])
basin = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(basin.parent / "calibration"))
from summa_python_interface import J, PARAM_NAMES

assert list(PARAM_NAMES) == names, "athabasca_params.csv order must match PARAM_NAMES"
kge = float(J(values))
res.write_text(repr(kge))
