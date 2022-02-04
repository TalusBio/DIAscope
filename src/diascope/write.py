"src/diascope/write.py"
import pandas as pd

from src.diascope.helpers import create_connection


def write_sql(out_path: str, metadata: pd.DataFrame, precursor: pd.DataFrame, ranges: pd.DataFrame, spectra: pd.DataFrame):
    dia_connection = create_connection(out_path)
    metadata.to_sql(name="metadata", con=dia_connection, index=False)
    precursor.to_sql(name="precursor", con=dia_connection, index=False)
    ranges.to_sql(name="ranges", con=dia_connection, index=False)
    spectra.to_sql(name="spectra", con=dia_connection, index=False)