#!/usr/bin/env python3

import GEOparse
import os

# Directory where the SOFT file will be downloaded
dest_dir = "./data"
os.makedirs(dest_dir, exist_ok=True)

# Download GSE55763
gse = GEOparse.get_GEO(
    geo="GSE55763",
    destdir=dest_dir,
    silent=False
)
