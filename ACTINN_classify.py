# ACTINN Script
# !git clone https://github.com/mafeiyang/ACTINN.git

import pandas as pd
import numpy as np
import gzip
import shutil

# Check the format of ACTINN label input
with gzip.open("/content/ACTINN/test_data/train_label.txt.gz", "rb") as f_in:
    with open("/content/sample_label.txt", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
