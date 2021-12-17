"""test/test_mzml_dia_equal.py"""
from pathlib import Path
import sqlite3

import pandas as pd
from pyteomics.mzml import MzML

from src.diascope.parse import parse_metadata, parse_spectra, parse_precursor, parse_ranges

DATA_DIR = Path(__file__).resolve().parent.joinpath("data")

FILE_NAME = "small"
MZML_PATH = DATA_DIR.joinpath(f"{FILE_NAME}.mzML")
DIA_PATH = DATA_DIR.joinpath(f"{FILE_NAME}.dia")


def test_mzml_dia_metadata_equal():
    dia_connection = sqlite3.connect(DIA_PATH)

    dia_metadata = pd.read_sql_query(sql="SELECT * FROM metadata;", con=dia_connection)
    mzml_metadata = parse_metadata(mzml_path=str(MZML_PATH))

    pd.testing.assert_frame_equal(dia_metadata, mzml_metadata)


def test_mzml_dia_precursor_equal():
    mzml_file = MzML(str(MZML_PATH))
    dia_connection = sqlite3.connect(DIA_PATH)

    dia_precursor = pd.read_sql_query(sql="SELECT * FROM precursor;", con=dia_connection)
    mzml_precursor = parse_precursor(mzml_file)

    dia_precursor = dia_precursor.drop(["MassArray", "IntensityArray"], axis=1)
    mzml_precursor = mzml_precursor.drop(["MassArray", "IntensityArray"], axis=1)
    
    pd.testing.assert_frame_equal(dia_precursor, mzml_precursor)


def test_mzml_dia_spectra_equal():
    mzml_file = MzML(str(MZML_PATH))
    dia_connection = sqlite3.connect(DIA_PATH)

    dia_spectra = pd.read_sql_query(sql="SELECT * FROM spectra;", con=dia_connection)
    mzml_spectra = parse_spectra(mzml_file)
    
    dia_spectra = dia_spectra.drop(["MassArray", "IntensityArray"], axis=1)
    mzml_spectra = mzml_spectra.drop(["MassArray", "IntensityArray"], axis=1)

    pd.testing.assert_frame_equal(dia_spectra, mzml_spectra)


def test_mzml_dia_ranges_equal():
    mzml_file = MzML(str(MZML_PATH))
    dia_connection = sqlite3.connect(DIA_PATH)

    dia_ranges = pd.read_sql_query(sql="SELECT * FROM ranges;", con=dia_connection)
    mzml_ranges = parse_ranges(mzml_file)

    dia_ranges = dia_ranges.sort_values(by="Start").reset_index(drop=True)
    mzml_ranges = mzml_ranges.sort_values(by="Start").reset_index(drop=True)

    pd.testing.assert_frame_equal(dia_ranges, mzml_ranges)