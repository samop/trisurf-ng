#!/usr/bin/python3
import xml.etree.ElementTree as ET
import base64
import zlib

tree = ET.parse('../src/timestep_000000.vtu')
root = tree.getroot()
trisurf=root.find('trisurf')
version=root.find('trisurfversion')
print(version.text);
print(trisurf.items())

#xml=zlib.decompress(base64.b64decode(trisurf.text))
#xml=trisurf.text
#print(xml)

#tree2=ET.ElementTree(ET.fromstring(str(xml))
#tree2=ET.ElementTree(ET.fromstring("<root>"+str(xml)[0:-1]+"</root>"))
