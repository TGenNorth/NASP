#!/usr/bin/env python3

__author__ = "Darrin Lemmer"
__version__ = "0.9.8"
__email__ = "dlemmer@tgen.org"

'''
Created on April 3, 2014

@author: dlemmer
'''
from xml.etree import ElementTree
from concurrent.futures import ProcessPoolExecutor
from collections import namedtuple
import itertools

from nasp2.parse import Vcf, Fasta


MatrixParameters = namedtuple('MatrixParameters', ['reference_fasta', 'reference_dups', 'minimum_coverage', 'minimum_proportion', 'matrix_folder', 'stats_folder', 'parameters', 'filter_matrix_format'])
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


def write_dto(matrix_parms, franken_fastas, vcf_files, xml_file):
    from xml.dom import minidom

    root = ElementTree.Element("matrix_data")
    parm_node = ElementTree.SubElement(root, "parameters")
    _write_parameters(parm_node, matrix_parms)
    files_node = ElementTree.SubElement(root, "files")
    # TODO: Can 'name' be redefined to be just the sample name instead of the filename without the extension?
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
    """
    Args:
        xml_file:

    Returns:
        MatrixParameters, tuple of Analyses
    """
    xmltree = ElementTree.parse(xml_file)
    root = xmltree.getroot()

    futures = []
    with ProcessPoolExecutor() as executor:
        # Parse matrix parameters indexing the reference and duplicates file in parallel.
        matrix_params = {element.tag.replace('-', '_'): element.text for element in root.find("parameters").iter()}
        matrix_params['reference_fasta'] = executor.submit(Fasta, matrix_params['reference_fasta'], 'reference', None, True)
        # TODO: handle undefined reference_dups
        matrix_params['reference_dups'] = executor.submit(Fasta, matrix_params['reference_dups'], 'reference', None)
        matrix_params['minimum_coverage'] = float(matrix_params['minimum_coverage'])
        matrix_params['minimum_proportion'] = float(matrix_params['minimum_proportion'])

        # Index Vcf and Frankenfastas in parallel.
        for frankenfasta in root.find("files").iter("frankenfasta"):
            futures.append(executor.submit(Fasta, frankenfasta.text, frankenfasta.get("name"), frankenfasta.get("aligner")))
        for vcf in root.find("files").iter("vcf"):
            futures.append(executor.submit(Vcf, vcf.text, vcf.get("name"), vcf.get("aligner"), vcf.get("snpcaller")))

        # Return when the indexing processes complete.
        matrix_params['reference_fasta'] = matrix_params['reference_fasta'].result()
        matrix_params['reference_dups'] = matrix_params['reference_dups'].result()
        # TODO: try/catch futures exception to return a failed genome object and write the error to the parse log.
        # Group the analyses by sample name in order to collect sample-level statistics.
        # The SampleAnalyses are sorted before grouping because groups are determined by when the key changes.
        # See analyse.sample_positions() for more details regarding the structure of sample_groups
        sample_analyses = tuple(tuple(v) for _, v in itertools.groupby(sorted(future.result() for future in futures), lambda x: x.name))
        return MatrixParameters(**matrix_params), sample_analyses


def main():
    pass


if __name__ == "__main__":
    main()
