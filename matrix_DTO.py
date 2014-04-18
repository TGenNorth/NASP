#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "0.9.3"
__email__ = "dlemmer@tgen.org"

'''
Created on April 3, 2014

@author: dlemmer
'''
from xml.etree import ElementTree

def _write_parameters( node, data ):
    for k,v in data.items():
        subnode = ElementTree.SubElement(node, k)
        subnode.text = v
    return node

def _add_input_file( node, filetype, attributes, file ):
    subnode = ElementTree.SubElement(node, filetype, attributes)
    subnode.text = file
    return subnode

def _parse_parameters( node ):
    parms = {}
    for element in node.iter():
        parms[element.tag] = element.text
    return parms

def _parse_files( node ):
    input_files = []
    for frankenfasta in node.iter("frankenfasta"):
        input_files.append("%s,%s,::%s" % ("frankenfasta", frankenfasta.get("aligner"), frankenfasta.text))
    for vcf in node.iter("vcf"):
        input_files.append("%s,%s,%s::%s" % ("vcf", vcf.get("aligner"), vcf.get("snpcaller"), vcf.text))

def write_dto( matrix_parms, franken_fastas, vcf_files, xml_file ):
    root = ElementTree.Element("matrix_data")
    parm_node = ElementTree.SubElement(root, "parameters")
    _write_parameters(parm_node, matrix_parms)
    files_node = ElementTree.SubElement(root, "files")
    for (name, aligner, file) in franken_fastas:
        attributes = {'name':name, 'aligner':aligner}
        _add_input_file(files_node, "frankenfasta", attributes, file)
    for (name, aligner, snpcaller, file) in vcf_files:
        attributes = {'name':name, 'aligner':aligner, 'snpcaller':snpcaller}
        _add_input_file(files_node, "vcf", attributes, file)
    ElementTree.ElementTree(root).write(xml_file)
    return xml_file
    
def parse_dto( xml_file ):
    xmltree = ElementTree.parse( xml_file )
    root = xmltree.getroot()
    matrix_parms = _parse_parameters(root.find("parameters"))
    input_files = _parse_files(root.find("files"))
    return (matrix_parms, input_files)

def main():
    pass

if __name__ == "__main__": main()
