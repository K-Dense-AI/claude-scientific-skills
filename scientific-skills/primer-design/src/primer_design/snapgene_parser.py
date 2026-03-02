"""SnapGene .dna binary file parser."""

import xml.etree.ElementTree as ET


def parse_snapgene(filepath):
    """Parse SnapGene .dna binary file.

    Returns
    -------
    tuple : (sequence: str, is_circular: bool, features: list[dict])
        features: [{"name", "type", "start", "end", "strand"}, ...]
        strand: 1 = forward, 2 = reverse complement
    """
    with open(filepath, "rb") as fh:
        data = fh.read()

    sequence = None
    is_circular = False
    features = []

    i = 0
    while i < len(data) - 5:
        btype = data[i]
        blen = int.from_bytes(data[i + 1:i + 5], "big")
        if i + 5 + blen > len(data):
            break

        if btype == 0:
            flags = data[i + 5]
            is_circular = bool(flags & 0x01)
            sequence = data[i + 6:i + 5 + blen].decode("ascii")

        elif btype == 10:
            xml_str = data[i + 5:i + 5 + blen].decode("utf-8", errors="replace")
            root = ET.fromstring(xml_str)
            for feat in root.findall(".//Feature"):
                name = feat.get("name", "")
                ftype = feat.get("type", "")
                direc = int(feat.get("directionality", "1"))
                for seg in feat.findall("Segment"):
                    rng = seg.get("range", "")
                    if "-" in rng:
                        s, e = rng.split("-", 1)
                        features.append({
                            "name": name,
                            "type": ftype,
                            "start": int(s),
                            "end": int(e),
                            "strand": direc,
                        })

        i += 5 + blen

    return sequence, is_circular, features
