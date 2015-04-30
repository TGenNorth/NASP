#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "0.9.8"
__email__ = "dlemmer@tgen.org"

'''
Created on April 3, 2014

@author: dlemmer
'''
from xml.etree import ElementTree
from collections import namedtuple

# TODO: Remove
# MatrixParameters = namedtuple('MatrixParameters', [
#     'reference_fasta',
#     'reference_dups',
#     'minimum_coverage',
#     'minimum_proportion',
#     'matrix_folder',
#     'stats_folder',
#     'parameters',
#     'filter_matrix_format'
# ])
# MatrixParameters.__new__.__defaults__ = ("", "", "", "", "", "", "")

NaspFile = namedtuple('NaspFile', ['path', 'name', 'aligner', 'snpcaller'])


def _write_parameters(node, data):
    for k, v in data.items():
        subnode = ElementTree.SubElement(node, k)
        subnode.text = v
    return node


def _add_input_file(node, filetype, attributes, filepath):
    subnode = ElementTree.SubElement(node, filetype, attributes)
    subnode.text = filepath
    return subnode


def write_dto(matrix_parms, franken_fastas, vcf_files, xml_file):
    from xml.dom import minidom
    import re

    root = ElementTree.Element("matrix_data")
    parm_node = ElementTree.SubElement(root, "parameters")
    _write_parameters(parm_node, matrix_parms)
    files_node = ElementTree.SubElement(root, "files")
    for (name, aligner, filepath) in franken_fastas:
        attributes = {'name': name, 'aligner': aligner}
        _add_input_file(files_node, "frankenfasta", attributes, filepath)
    for (name, aligner, snpcaller, filepath) in vcf_files:
        pattern = re.compile('^(.*?)(?:-(pre-aligned|bwa(mem)?|novo|bowtie2|snap]))?(?:-(pre-called|gatk|solsnp|varscan|samtools]))?$')
        match = pattern.match(name)
        if match:
            name = match.group(1)
        attributes = {'name': name, 'aligner': aligner, 'snpcaller': snpcaller}
        _add_input_file(files_node, "vcf", attributes, filepath)
    dom = minidom.parseString(ElementTree.tostring(root, 'utf-8'))
    with open(xml_file, 'w') as output:
        output.write(dom.toprettyxml(indent=" " * 5))
    return xml_file


def parse_dto(xml_file):
    """
    Args:
        xml_file:

    Returns:
        dict: The xml file as a dictionary.
    """
    xmltree = ElementTree.parse(xml_file)
    root = xmltree.getroot()

    matrix_params = {element.tag.replace('-', '_'): element.text for element in root.find("parameters").iter()}
    matrix_params['frankenfasta'] = tuple(
        NaspFile(
            path=frankenfasta.text,
            name=frankenfasta.get("name"),
            aligner=frankenfasta.get("aligner"),
            snpcaller=None
        ) for frankenfasta in root.find("files").iter("frankenfasta")
    )
    matrix_params['vcf'] = tuple(
        NaspFile(
            path=vcf.text,
            name=vcf.get("name"),
            aligner=vcf.get("aligner"),
            snpcaller=vcf.get("snpcaller")
        ) for vcf in root.find("files").iter("vcf")
    )
    return matrix_params


def main():
    pass


if __name__ == "__main__":
    main()
