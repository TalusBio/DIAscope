"src/diascope/parse.py"
import os
import zlib
from statistics import mean

from lxml import etree
from numpy import array
import pandas as pd
from pyteomics.mzml import MzML

from .helpers import etree_to_dict, array_pack

MULTI_INSTRUMENT_DELIMITER = "[MULTI-INSTRUMENT-DELIMITER]"
INSTRUMENT_ID_COMPONENT_DELIMITER = "[INSTRUMENT-ID-COMPONENT-DELIMITER]"

mzml_file = MzML("/Users/ricomeinl/Desktop/talus/DIAscope/test/data/small.mzML")


def get_list_from_obj(list_obj, target_name):
    count = list_obj.get("count", 0)
    return [list_obj[target_name][i] for i in range(count)]


def get_scan_times(scan_list, time_type, unit_multiplier=1.0, default_val=-1.0):
    scan_list_objs = get_list_from_obj(scan_list, target_name="scan")
    scan_times = [scan_obj.get(time_type) for scan_obj in scan_list_objs]
    return [scan_t * unit_multiplier if scan_t else default_val for scan_t in scan_times]


def get_isolation_windows(precursor_list):
    def extract_window(window_obj):
        center = window_obj.get("isolation window target m/z", 0.0)
        lower = center - window_obj.get("isolation window lower offset", 0.0)
        upper = center + window_obj.get("isolation window upper offset", 0.0)
        return (center, lower, upper)
    precursor_objs = get_list_from_obj(precursor_list, target_name="precursor")
    return [extract_window(precursor_obj.get("isolationWindow")) for precursor_obj in precursor_objs]


def get_scan_windows(scan_window_list):
    def extract_window(window_obj):
        lower = window_obj.get("scan window lower limit", 0.0)
        upper = window_obj.get("scan window upper limit", 0.0)
        return (lower, upper)
    scan_list_objs = [scan_list_obj.get("scanWindowList") for scan_list_obj in get_list_from_obj(scan_window_list, target_name="scan")]
    scan_window_objs = [get_list_from_obj(scan_list, target_name="scanWindow") for scan_list in scan_list_objs]
    return [[extract_window(scan_window) for scan_window in scan_window_obj_list] for scan_window_obj_list in scan_window_objs]


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
        if spectrum.get("ms level", 0) > 1:
            packed_mz_arr = array_pack(spectrum.get("m/z array"), format_character='d')
            packed_intensity_arr = array_pack(spectrum.get("intensity array"), format_character='f')
            scan_start_time = get_scan_times(scan_list=spectrum.get("scanList"), time_type="scan start time", unit_multiplier=60.0)[0]
            ion_injection_time = get_scan_times(scan_list=spectrum.get("scanList"), time_type="ion injection time", unit_multiplier=0.001)[0]
            iso_window_center, iso_window_lower, iso_window_upper = get_isolation_windows(precursor_list=spectrum.get("precursorList"))[0]
            precursor_name = [precursor_obj.get("spectrumRef") for precursor_obj in get_list_from_obj(spectrum.get("precursorList"), target_name="precursor")][0]
            
            spectrum_dict = dict(
                Fraction=0,
                SpectrumName=spectrum.get("id"),
                PrecursorName=precursor_name,
                SpectrumIndex=spectrum.get("index"),
                ScanStartTime=scan_start_time,
                IonInjectionTime=ion_injection_time,
                IsolationWindowLower=iso_window_lower,
                IsolationWindowCenter=iso_window_center,
                IsolationWindowUpper=iso_window_upper,
                MassEncodedLength=len(packed_mz_arr),
                MassArray=zlib.compress(packed_mz_arr),
                IntensityEncodedLength=len(packed_intensity_arr),
                IntensityArray=zlib.compress(packed_intensity_arr),
            )
            spectra_list.append(spectrum_dict)
    return pd.DataFrame(spectra_list)


def parse_precursor(mzml_file):
    precursor_list = []
    for spectrum in mzml_file:
        if spectrum.get("ms level", 0) == 1:
            packed_mz_arr = array_pack(spectrum.get("m/z array"), format_character='d')
            packed_intensity_arr = array_pack(spectrum.get("intensity array"), format_character='f')
            scan_start_time = get_scan_times(scan_list=spectrum.get("scanList"), time_type="scan start time", unit_multiplier=60.0)[0]
            ion_injection_time = get_scan_times(scan_list=spectrum.get("scanList"), time_type="ion injection time", unit_multiplier=0.001)[0]
            scan_window_lower, scan_window_upper = get_scan_windows(scan_window_list=spectrum.get("scanList"))[0][0]

            precursor_dict = dict(
                Fraction=0,
                SpectrumName=spectrum.get("id"),
                SpectrumIndex=spectrum.get("index"),
                ScanStartTime=scan_start_time,
                IonInjectionTime=ion_injection_time,
                IsolationWindowLower=scan_window_lower,
                IsolationWindowUpper=scan_window_upper,
                MassEncodedLength=len(packed_mz_arr),
                MassArray=zlib.compress(packed_mz_arr),
                IntensityEncodedLength=len(packed_intensity_arr),
                IntensityArray=zlib.compress(packed_intensity_arr),
                TIC=spectrum.get("total ion current"),
            )
            precursor_list.append(precursor_dict)
    return pd.DataFrame(precursor_list)


def parse_ranges(mzml_file):
    def get_deltas(times): return [times[i]-times[i-1] for i in range(1, len(times))] or [0]

    ranges_dict = {}
    for spectrum in mzml_file:
        if spectrum.get("ms level", 0) > 1:
            _, iso_window_lower, iso_window_upper = get_isolation_windows(precursor_list=spectrum.get("precursorList"))[0]
            range_ = (iso_window_lower, iso_window_upper)
            scan_start_time = get_scan_times(scan_list=spectrum.get("scanList"), time_type="scan start time", unit_multiplier=60.0)[0]
            ranges_dict[range_] = ranges_dict.get(range_, []) + [scan_start_time]

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
