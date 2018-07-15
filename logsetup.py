

import logging
import sys

def setup_log():

    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    df = logging.FileHandler("debug.log")
    df.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    df.setFormatter(formatter)
    root.addHandler(df)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    root.addHandler(ch)