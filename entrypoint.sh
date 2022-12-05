#!/bin/bash --login
set +euo pipefail
conda activate lfmtools-docker

set -euo pipefail
exec streamlit run /app/lfmtools/Home.py --server.port=8501 --server.address=0.0.0.0