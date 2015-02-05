#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "0.9.8"
__email__ = "dlemmer@tgen.org"

'''
Created on April 3, 2014

@author: dlemmer
'''
import logging
from xml.etree import ElementTree
from concurrent.futures import ProcessPoolExecutor
from collections import namedtuple
from nasp2.parse import Vcf, Fasta

MatrixParameters = namedtuple('MatrixParameters', ['reference_fasta', 'reference_dups', 'minimum_coverage', 'minimum_proportion', 'matrix_folder', 'stats_folder', 'filter_matrix_format'])
MatrixParameters.__new__.__defaults__ = ("", "", "", "", "", "", "")


def _write_parameters(node, data):
    for k, v in data.items():
        subnode = ElementTree.SubElement(node, k)
        subnode.text = v
    return node


def _add_input_file(node, filetype, attributes, file):
    subnode = ElementTree.SubElement(node, filetype, attributes)
    subnode.text = file
    return subnode


def _parse_files(node):
    """
    Returns a tuple of Vcf and Fasta objects indexed in parallel.

    Args:
        node:

    Return:
        tuple of Sample:
    """
    with ProcessPoolExecutor() as executor:
        futures = []
        for frankenfasta in node.iter("frankenfasta"):
            futures.append(executor.submit(Fasta, frankenfasta.text, frankenfasta.get("name"), frankenfasta.get("aligner")))
            # input_files.append("%s,%s,::%s" % ("frankenfasta", frankenfasta.get("aligner"), frankenfasta.text))
        for vcf in node.iter("vcf"):
            futures.append(executor.submit(Vcf, vcf.text, vcf.get("name"), vcf.get("aligner"), vcf.get("snpcaller")))
            # input_files.append("%s,%s,%s,::%s" % ("vcf", vcf.get("aligner"), vcf.get("snpcaller"), vcf.text))
        return tuple(future.result() for future in futures)


def write_dto(matrix_parms, franken_fastas, vcf_files, xml_file):
    from xml.dom import minidom

    root = ElementTree.Element("matrix_data")
    parm_node = ElementTree.SubElement(root, "parameters")
    _write_parameters(parm_node, matrix_parms)
    files_node = ElementTree.SubElement(root, "files")
    for (name, aligner, file) in franken_fastas:
        attributes = {'name': name, 'aligner': aligner}
        _add_input_file(files_node, "frankenfasta", attributes, file)
    for (name, aligner, snpcaller, file) in vcf_files:
        attributes = {'name': name, 'aligner': aligner, 'snpcaller': snpcaller}
        _add_input_file(files_node, "vcf", attributes, file)
    dom = minidom.parseString(ElementTree.tostring(root, 'utf-8'))
    output = open(xml_file, 'w')
    output.write(dom.toprettyxml(indent="    "))
    output.close()
    return xml_file


def parse_dto(xml_file):
    xmltree = ElementTree.parse(xml_file)
    root = xmltree.getroot()
    matrix_parms = MatrixParameters({element.tag: element.text for element in root.find("parameters").iter()})
    input_files = _parse_files(root.find("files"))
    return matrix_parms, input_files


def main():
    pass


if __name__ == "__main__":
    main()
