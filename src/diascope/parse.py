"src/diascope/parse.py"
import os
import zlib
from statistics import mean

from lxml import etree
import pandas as pd
from pyteomics.mzml import MzML

from .helpers import etree_to_dict

MULTI_INSTRUMENT_DELIMITER = "[MULTI-INSTRUMENT-DELIMITER]"
INSTRUMENT_ID_COMPONENT_DELIMITER = "[INSTRUMENT-ID-COMPONENT-DELIMITER]"

mzml_file = MzML("/Users/ricomeinl/Desktop/talus/DIAscope/test/data/small.mzML")


def parse_metadata(mzml_path):
    tree = etree.parse(mzml_path)
    metadata_list = []
    # TODO: replace with the real
    metadata_list.append(("totalPrecursorTIC", 0.0))
    software_nodes = etree_to_dict([elem for elem in tree.getroot()[0] if "software" in str(elem)][0])
    metadata_list += [(f"SoftwareVersion_{node['cvParam']['accession']}", node["version"]) for node in software_nodes["softwareList"]["software"]]
    instrument_nodes = etree_to_dict([elem for elem in tree.getroot()[0] if "instrument" in str(elem)][0])
    instrument_str = ""
    for instrument in instrument_nodes["instrumentConfigurationList"]["instrumentConfiguration"]:
        instrument_str += f"configurationId:{instrument.get('id', 'null')},"
        instrument_str += f"accession:{instrument.get('accession', 'null')},"
        instrument_str += f"name:{instrument.get('name', 'null')}"
        instrument_str += INSTRUMENT_ID_COMPONENT_DELIMITER
        for type, type_obj in instrument["componentList"].items():
            if type not in ["source", "analyzer", "detector"]:
                continue
            name_objs = type_obj["cvParam"] if isinstance(type_obj["cvParam"], list) else [type_obj["cvParam"]]
            for name_obj in name_objs:
                instrument_str += f"order:{type_obj.get('order', 'null')},"
                instrument_str += f"cvRef:{name_obj.get('cvRef', 'null')},"
                instrument_str += f"accessionId:{name_obj.get('accession', 'null')},"
                instrument_str += f"name:{name_obj.get('name', 'null')},"
                instrument_str += f"type:{type};"
        instrument_str += MULTI_INSTRUMENT_DELIMITER
    metadata_list.append(("InstrumentConfigurations", instrument_str))
    metadata_list.append(("filelocation", mzml_path))
    metadata_list.append(("filename", os.path.basename(mzml_path)))
    metadata_list.append(("sourcename", tree.getroot()[0].attrib["id"]))
    metadata_list.append(("version", "0.0.0"))
    return pd.DataFrame(metadata_list, columns=["Key", "Value"])


def parse_spectra(mzml_file):
    spectra_list = []
    for spectrum in mzml_file:
        if spectrum["ms level"] > 1:
            # TODO: fix mass and intensity array encoding
            isolation_window_center = spectrum["precursorList"]["precursor"][0]["isolationWindow"]["isolation window target m/z"]
            spectrum_dict = dict(
                Fraction=0,
                SpectrumName=spectrum["id"],
                PrecursorName=spectrum["precursorList"]["precursor"][0]["spectrumRef"],
                SpectrumIndex=spectrum["index"],
                ScanStartTime=spectrum["scanList"]["scan"][0]["scan start time"]*60,
                IonInjectionTime=-1.0,
                IsolationWindowLower=isolation_window_center - spectrum["precursorList"]["precursor"][0]["isolationWindow"]["isolation window lower offset"],
                IsolationWindowCenter=isolation_window_center,
                IsolationWindowUpper=isolation_window_center + spectrum["precursorList"]["precursor"][0]["isolationWindow"]["isolation window upper offset"],
                MassEncodedLength=len(bytearray(spectrum["m/z array"])),
                MassArray=zlib.compress(bytearray(spectrum["m/z array"])),
                IntensityEncodedLength=len(bytearray(spectrum["intensity array"])),
                IntensityArray=zlib.compress(bytearray(spectrum["intensity array"])),
            )
            spectra_list.append(spectrum_dict)
    return pd.DataFrame(spectra_list)


def parse_precursor(mzml_file):
    precursor_list = []
    for spectrum in mzml_file:
        if spectrum["ms level"] == 1:
            # TODO: fix mass and intensity array encoding; decode their version and reencode it
            precursor_dict = dict(
                Fraction=0,
                SpectrumName=spectrum["id"],
                SpectrumIndex=spectrum["index"],
                ScanStartTime=spectrum["scanList"]["scan"][0]["scan start time"]*60,
                IonInjectionTime=-1.0,
                IsolationWindowLower=spectrum["scanList"]["scan"][0]["scanWindowList"]["scanWindow"][0]["scan window lower limit"],
                IsolationWindowUpper=spectrum["scanList"]["scan"][0]["scanWindowList"]["scanWindow"][0]["scan window upper limit"],
                MassEncodedLength=len(bytearray(spectrum["m/z array"])),
                MassArray=zlib.compress(bytearray(spectrum["m/z array"])),
                IntensityEncodedLength=len(bytearray(spectrum["intensity array"])),
                IntensityArray=zlib.compress(bytearray(spectrum["intensity array"])),
                TIC=spectrum["total ion current"],
            )
            precursor_list.append(precursor_dict)
    return pd.DataFrame(precursor_list)


def parse_ranges(mzml_file):
    def get_deltas(times): return [times[i]-times[i-1] for i in range(1, len(times))] or [0]

    ranges_dict = {}
    for spectrum in mzml_file:
        if spectrum["ms level"] > 1:
            isolation_window_center = spectrum["precursorList"]["precursor"][0]["isolationWindow"]["isolation window target m/z"]
            range_start = isolation_window_center - spectrum["precursorList"]["precursor"][0]["isolationWindow"]["isolation window lower offset"]
            range_stop = isolation_window_center + spectrum["precursorList"]["precursor"][0]["isolationWindow"]["isolation window upper offset"]
            range_ = (range_start, range_stop)
            ranges_dict[range_] = ranges_dict.get(range_, []) + [spectrum["scanList"]["scan"][0]["scan start time"]*60]

    ranges_list = []
    for (range_start, range_stop), rts in ranges_dict.items():
        deltas = get_deltas(rts)
        avg_duty_cycle = mean(deltas)
        range_dict = dict(
            Start=range_start,
            Stop=range_stop,
            DutyCycle=avg_duty_cycle,
            NumWindows=len(rts),
        )
        ranges_list.append(range_dict)
    return pd.DataFrame(ranges_list)
