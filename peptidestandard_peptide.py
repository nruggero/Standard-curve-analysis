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
	print "Starting peptide standard peptide analysis"
	print "Opening file {}".format(data_infile)
	peptide_reference_full_list	= []
	peptide_full_list		= []
	intensity_full_list		= []
	area_full_list			= []
	with open(data_infile, 'rU') as csv_infile:
		dict_reader = csv.DictReader(csv_infile)
		for row in dict_reader:
			peptide_reference_full_list.append(row['Reference'])
			peptide_full_list.append(row['Peptide'])
			intensity_full_list.append(float(row['Precursor Intensity']))
			area_full_list.append(float(row['Area']))

	# Analyze data
	print "Performing analysis..."

	## Scrubbing out all peptide references that do not start with "gi"
	peptide_reference_idx = [i for i,x in enumerate(peptide_reference_full_list) if x[:2] == 'gi']

	peptide_reference_list 	= [peptide_reference_full_list[idx] for idx in peptide_reference_idx]
	peptide_list 			= [peptide_full_list[idx] for idx in peptide_reference_idx]
	intensity_list 			= [intensity_full_list[idx] for idx in peptide_reference_idx]
	area_list 				= [area_full_list[idx] for idx in peptide_reference_idx]

	## Sorting peptides and associated area and intensity data based on protein
	peptide_dict = {}

	for idx in range(len(peptide_list)):
		peptide = peptide_list[idx]

		if not PSM:
			peptide = removeIgnoreCharacters(peptide)

		if peptide_dict.has_key(peptide_list[idx]):
			peptide_dict[peptide]['area'].append(area_list[idx])
			peptide_dict[peptide]['intensity'].append(intensity_list[idx])
		else:
			peptide_dict[peptide] = {'area' : [area_list[idx]], 'intensity' : [intensity_list[idx]], 'reference' : [peptide_reference_list[idx]]}

	## Save output computing averages/intensities on the fly
	print "Writing to file {}".format(data_outfile)
	with open(data_outfile, 'wb') as csv_outfile:
		csv_writer = csv.writer(csv_outfile)

		# Output area data
		csv_writer.writerow(['peptide', 'mean area', 'peptide id'])

		for peptide in peptide_dict.iterkeys():
			csv_writer.writerow([
				peptide,
				averageList(peptide_dict[peptide]['area']),
				peptide_dict[peptide]['reference'][0]
				])
			if len(peptide_dict[peptide]['reference']) > 1:
				raise Exception, "Peptide {} has more than one reference ID!\n".format(peptide)

		# # Write blank row
		# csv_writer.writerow([])
		# csv_writer.writerow(['average of top 3 intensities', 'peptides used'])

		# # Output intensity data
		# peptide_maxes_for_protein = [(max(peptide_dict[peptide]['intensity']), peptide) for peptide in peptide_dict.iterkeys()]
		# peptide_maxes_for_protein.sort()
		# top_3 = peptide_maxes_for_protein[-3:]
		# average_top_3 = averageList([x[0] for x in top_3])

		# csv_writer.writerow([average_top_3, top_3])

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("data_infile", help = "CSV file containing analyzed output", type = str)
	parser.add_argument("data_outfile", help = "CSV file containing raw input", type = str)
	parser.add_argument("PSM", help = "Boolean toggle for ignoring special characters (1 = peptide spectral match)", type = int)

	args = parser.parse_args().__dict__

	main(args["data_infile"], args["data_outfile"], args["PSM"])