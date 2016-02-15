#!/usr/bin/python3
import xml.etree.ElementTree as ET


tree = ET.parse('../src/timestep_000000.vtu')
root = tree.getroot()
trisurf=root.find('trisurf')
print(trisurf.items())


