import os
import pdb
import xml.etree.ElementTree as ET

def xml_list(folder):
    """Return a list of all *xml file in folder"""
    list_of_xml_file = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(".xml"):
                 list_of_xml_file.append(os.path.join(root, file))

    return list_of_xml_file

def read_single_xml(filepath):
    tree = ET.parse(filepath)
    elements = tree.getroot().getchildren()
    # TODO: how to raise exception to not evaludate this function
    # and all functions after that one
    obs = elements[0].text
    lgn = elements[1].text

    return obs, lgn

if __name__ == '__main__':
    xml_list = xml_list(XML_MAPPING_FOLDER)

    for xml_file in xml_list:
        read_single_xml(xml_file)
    print xml_list