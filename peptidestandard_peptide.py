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
from standard_curve_library import removeIgnoreCharacters, averageList, findIndexClosestToNumber, median

def main(data_infile, data_outfile, PSM):
	# Open file and load data into lists
	print "Starting peptide standard peptide analysis"
	print "Opening file {}".format(data_infile)
	peptide_reference_full_list	= []
	peptide_full_list		= []
	intensity_full_list		= []
	elutiontime_full_list	= []
	maxtime_full_list 		= []
	area_full_list			= []
	with open(data_infile, 'rU') as csv_infile:
		dict_reader = csv.DictReader(csv_infile)
		for row in dict_reader:
			peptide_reference_full_list.append(row['Reference'])
			peptide_full_list.append(row['Peptide'])
			elutiontime_full_list.append(float(row['Time']))
			maxtime_full_list.append(float(row['Max Time']))
			intensity_full_list.append(float(row['Precursor Intensity']))
			area_full_list.append(float(row['Area']))

	# Analyze data
	print "Performing analysis..."

	## Scrubbing out all peptide references that do not start with "gi"
	peptide_reference_idx = [i for i,x in enumerate(peptide_reference_full_list) if x[:2] == 'gi']

	peptide_reference_list 	= [peptide_reference_full_list[idx] for idx in peptide_reference_idx]
	peptide_list 			= [peptide_full_list[idx] for idx in peptide_reference_idx]
	elutiontime_list 		= [elutiontime_full_list[idx] for idx in peptide_reference_idx]
	maxtime_list 			= [maxtime_full_list[idx] for idx in peptide_reference_idx]
	intensity_list 			= [intensity_full_list[idx] for idx in peptide_reference_idx]
	area_list 				= [area_full_list[idx] for idx in peptide_reference_idx]

	## Sorting peptides and associated area and intensity data based on protein
	peptide_dict = {}

	for idx in range(len(peptide_list)):
		peptide = peptide_list[idx]

		if not PSM:
			peptide = removeIgnoreCharacters(peptide)

		if peptide_dict.has_key(peptide):
			peptide_dict[peptide]['elutiontime'].append(elutiontime_list[idx])
			peptide_dict[peptide]['maxtime'].append(maxtime_list[idx])
			peptide_dict[peptide]['area'].append(area_list[idx])
			peptide_dict[peptide]['intensity'].append(intensity_list[idx])
		else:
			peptide_dict[peptide] = {
									'elutiontime' : [elutiontime_list[idx]],
									'maxtime' : [maxtime_list[idx]],
									'area' : [area_list[idx]],
									'intensity' : [intensity_list[idx]],
									'reference' : [peptide_reference_list[idx]],
									'closest_to_maxtime_idx' : None
									}


	## Find elution time closest to maxtime
	for peptide in peptide_dict.iterkeys():
		if len(set(peptide_dict[peptide]['maxtime'])) > 1:
			print '{} has more than one maxtime. Picking largest area {}!'.format(peptide, max(peptide_dict[peptide]['area']))
			largest_area = max(peptide_dict[peptide]['area'])
			largest_area_index = peptide_dict[peptide]['area'].index(largest_area)
			maxtime = peptide_dict[peptide]['maxtime'][largest_area_index]
		else:
			maxtime = peptide_dict[peptide]['maxtime'][0]
		peptide_dict[peptide]['closest_to_maxtime_idx'] = findIndexClosestToNumber(
															maxtime, peptide_dict[peptide]['elutiontime']
															)

	## Save output computing averages/intensities on the fly
	print "Writing to file {}".format(data_outfile)
	with open(data_outfile, 'wb') as csv_outfile:
		csv_writer = csv.writer(csv_outfile)

		# Output area data
		csv_writer.writerow(['peptide', 'area with elution time closest to maxtime', 'sum of area', 'median of area', 'peptide id'])

		allAreas = []
		for x in peptide_dict.itervalues():
			allAreas.extend(x['area'])
		sumOfRunArea = sum(allAreas)
		medianOfRunArea = median(allAreas)

		for peptide in peptide_dict.iterkeys():
			csv_writer.writerow([
				peptide,
				peptide_dict[peptide]['area'][peptide_dict[peptide]['closest_to_maxtime_idx']],
				sumOfRunArea, #sum(peptide_dict[peptide]['area']),
				medianOfRunArea, #median(peptide_dict[peptide]['area']),
				peptide_dict[peptide]['reference'][0]
				])
			if len(peptide_dict[peptide]['reference']) > 1:
				raise Exception, "Peptide {} has more than one reference ID!\n".format(peptide)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("data_infile", help = "CSV file containing analyzed output", type = str)
	parser.add_argument("data_outfile", help = "CSV file containing raw input", type = str)
	parser.add_argument("PSM", help = "Boolean toggle for ignoring special characters (1 = peptide spectral match)", type = int)

	args = parser.parse_args().__dict__

	main(args["data_infile"], args["data_outfile"], args["PSM"])