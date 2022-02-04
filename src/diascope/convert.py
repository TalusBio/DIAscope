"""src/diascope/main.py"""
import click

from src.diascope.parse import parse_metadata, parse_spectra, parse_precursor, parse_ranges
from src.diascope.write import write_sql


@click.command()
@click.option("-i", "--input_file", help="Input file path")
@click.option("-o", "--output_file", default=None, help="Output file path")
@click.option("-f", "--format", default="mzML", help="Input file format")
def convert(input_file, output_file, format):
    """Convert file to .dia format."""
    if format == "mzML":
        if not output_file:
            output_file = input_file.replace("mzML", "dia")
        mzml_metadata = parse_metadata(mzml_path=input_file)
        mzml_precursor = parse_precursor(input_file)
        mzml_spectra = parse_spectra(input_file)
        mzml_ranges = parse_ranges(input_file)
        write_sql(
            out_path=output_file, 
            metadata=mzml_metadata, 
            precursor=mzml_precursor, 
            spectra=mzml_spectra, 
            ranges=mzml_ranges
        )
    else:
        raise Exception("File format not supported.")


if __name__ == "__main__":
    convert()