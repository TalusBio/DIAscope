"src/diascope/helpers.py"
import re
import struct
from collections import defaultdict


def etree_to_dict(t):
    pattern = r"{.*[^{}]}"
    parsed_tag = re.sub(pattern, "", t.tag)
    d = {parsed_tag : {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {parsed_tag: {k: v[0] if len(v) == 1 else v
                          for k, v in dd.items()}}
    if t.attrib:
        d[parsed_tag].update((k, v)
                        for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children or t.attrib:
            if text:
                d[parsed_tag]['#text'] = text
        else:
            d[parsed_tag] = text
    return d


def array_pack(arr, byte_order='>', format_character='f'):
    data_format = byte_order + format_character
    return b''.join(struct.pack(data_format, i) for i in arr)