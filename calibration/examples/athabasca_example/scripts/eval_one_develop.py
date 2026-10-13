"""One J() evaluation with the NCAR/summa develop build, started by calibrate_athabasca_asyncdds_mpi_develop.py.
Same as eval_one.py, but loads develop's libsumma.so through the SUMMA_LIB-aware interface in develop_port/.
Usage: python eval_one_develop.py <result_file> <names_json> <values_json>"""
import json, os, sys
from pathlib import Path

res, names, values = Path(sys.argv[1]), json.loads(sys.argv[2]), json.loads(sys.argv[3])
basin = Path(__file__).resolve().parent.parent
os.environ["SUMMA_LIB"] = str(basin.parent.parent / "summa-develop" / "build" / "cmake_build" / "libsumma.so")
sys.path.insert(0, str(basin / "scripts" / "develop_port"))
from summa_python_interface import J, PARAM_NAMES

assert list(PARAM_NAMES) == names, "athabasca_params.csv order must match PARAM_NAMES"
kge = float(J(values))
res.write_text(repr(kge))
