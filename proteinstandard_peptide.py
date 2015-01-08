#!/usr/bin/env python
"""
Generates peptide data for standards

@author: Nick Ruggero
@date: Created 12/01/2014
"""

import argparse
import os
import csv

from standard_curve_library import IGNORE_CHARACTERS
from standard_curve_library import removeIgnoreCharacters, averageList

def main(data_infile, data_outfile, PSM):
	# Open file and load data into lists
	print "Starting protein standard peptide analysis"
	print "Opening file {}".format(data_infile)
	protein_id_full_list	= []
	peptide_full_list		= []
	intensity_full_list		= []
	area_full_list			= []
	with open(data_infile, 'rU') as csv_infile:
		dict_reader = csv.DictReader(csv_infile)
		for row in dict_reader:
			protein_id_full_list.append(row['Reference'])
			peptide_full_list.append(row['Peptide'])
			intensity_full_list.append(float(row['Precursor Intensity']))
			area_full_list.append(float(row['Area']))

	# Analyze data
	print "Performing analysis..."

	## Creating unique list of protein id's
	## Index of protein id is used to index future lists of peptides, etc.
	protein_id_unique = list(set(protein_id_full_list))

	## Sorting peptides and associated area and intensity data based on protein
	peptide_data_by_protein = [{} for x in range(len(protein_id_unique))]

	for idx in range(len(peptide_full_list)):
		protein_dict = peptide_data_by_protein[protein_id_unique.index(protein_id_full_list[idx])]
		peptide = peptide_full_list[idx]

		if not PSM:
			peptide = removeIgnoreCharacters(peptide)

		if protein_dict.has_key(peptide):
			protein_dict[peptide]['area'].append(area_full_list[idx])
			protein_dict[peptide]['intensity'].append(intensity_full_list[idx])
		else:
			protein_dict[peptide] = {'area' : [area_full_list[idx]], 'intensity' : [intensity_full_list[idx]]}

	## Save output computing averages/intensities on the fly
	print "Writing to file {}".format(data_outfile)
	with open(data_outfile, 'wb') as csv_outfile:
		csv_writer = csv.writer(csv_outfile)

		# Output area data
		csv_writer.writerow(['protein id', 'peptide', 'mean area'])

		for protein_idx in range(len(protein_id_unique)):
			firstRow = True
			for peptide in peptide_data_by_protein[protein_idx].iterkeys():
				if firstRow:
					protein_identifier = protein_id_unique[protein_idx]
					firstRow = False
				else:
					protein_identifier = ''

				csv_writer.writerow([
					protein_identifier,
					peptide,
					averageList(peptide_data_by_protein[protein_idx][peptide]['area'])
					])

		# Write blank row
		csv_writer.writerow([])
		csv_writer.writerow(['protein id', 'average of top 3 intensities', 'peptides used'])

		# Output intensity data
		for protein_idx in range(len(protein_id_unique)):
			peptide_maxes_for_protein = [(max(peptide_data_by_protein[protein_idx][peptide]['intensity']), peptide) for peptide in peptide_data_by_protein[protein_idx].iterkeys()]
			peptide_maxes_for_protein.sort()
			top_3 = peptide_maxes_for_protein[-3:]
			average_top_3 = averageList([x[0] for x in top_3])

			csv_writer.writerow([protein_id_unique[protein_idx], average_top_3, top_3])

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("data_infile", help = "CSV file containing analyzed output", type = str)
	parser.add_argument("data_outfile", help = "CSV file containing raw input", type = str)
	parser.add_argument("PSM", help = "Boolean toggle for ignoring special characters (1 = peptide spectral match)", type = int)

	args = parser.parse_args().__dict__

	main(args["data_infile"], args["data_outfile"], args["PSM"])