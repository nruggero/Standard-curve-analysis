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
	print "Starting protein standard peptide analysis"
	print "Opening file {}".format(data_infile)
	protein_id_full_list		= []
	peptide_full_list			= []
	intensity_full_list			= []
	elutiontime_full_list		= []
	maxtime_full_list 			= []
	area_full_list				= []
	with open(data_infile, 'rU') as csv_infile:
		dict_reader = csv.DictReader(csv_infile)
		for row in dict_reader:
			protein_id_full_list.append(row['Reference'])
			peptide_full_list.append(row['Peptide'])
			elutiontime_full_list.append(float(row['Time']))
			maxtime_full_list.append(float(row['Max Time']))
			intensity_full_list.append(float(row['Precursor Intensity']))
			area_full_list.append(float(row['Area']))

	# Analyze data
	print "Performing analysis..."

	## Creating unique list of protein id's
	## Index of protein id is used to index future lists of peptides, etc.
	protein_id_unique = list(set(protein_id_full_list))

	## Sorting peptides and associated area and intensity data based on protein
	peptide_data_by_protein = dict([(x, {}) for x in protein_id_unique])
	
	for idx in range(len(peptide_full_list)):
		protein_dict = peptide_data_by_protein[protein_id_full_list[idx]]
		peptide = peptide_full_list[idx]

		if not PSM:
			peptide = removeIgnoreCharacters(peptide)

		if protein_dict.has_key(peptide):
			protein_dict[peptide]['elutiontime'].append(elutiontime_full_list[idx])
			protein_dict[peptide]['maxtime'].append(maxtime_full_list[idx])
			protein_dict[peptide]['area'].append(area_full_list[idx])
			protein_dict[peptide]['intensity'].append(intensity_full_list[idx])
		else:
			protein_dict[peptide] = {
									'elutiontime' : [elutiontime_full_list[idx]],
									'maxtime' : [maxtime_full_list[idx]],
									'area' : [area_full_list[idx]],
									'intensity' : [intensity_full_list[idx]],
									'closest_to_maxtime_idx' : None
									}
									
	## Find elution time closest to maxtime
	for protein in peptide_data_by_protein.iterkeys():
		peptide_dict = peptide_data_by_protein[protein]
		for peptide in peptide_dict.iterkeys():
			if len(set(peptide_dict[peptide]['maxtime'])) > 1:
				print '{} in {} has more than one maxtime. Picking largest area {}!'.format(peptide, protein, max(peptide_dict[peptide]['area']))
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
		csv_writer.writerow(['protein id', 'peptide', 'area with elution time closest to maxtime', 'sum of area', 'median of area'])

		for protein in protein_id_unique:
			firstRow = True
			peptide_dict  = peptide_data_by_protein[protein]
			if firstRow:
				protein_identifier = protein
				firstRow = False
			else:
				protein_identifier = ''
			
			allAreas = []
			for x in peptide_dict.itervalues():
				allAreas.extend(x['area'])
			sumOfRunArea = sum(allAreas)
			medianOfRunArea = median(allAreas)
			
			for peptide in peptide_dict.iterkeys():
				csv_writer.writerow([
					protein_identifier,
					peptide,
					peptide_dict[peptide]['area'][peptide_dict[peptide]['closest_to_maxtime_idx']],
					sumOfRunArea, #sum(peptide_dict[peptide]['area']),
					medianOfRunArea, #median(peptide_dict[peptide]['area'])
					])

		# Write blank row
		csv_writer.writerow([])
		csv_writer.writerow(['protein id', 'average of top 3 intensities', 'peptides used'])

		# Output intensity data
		for protein in protein_id_unique:
			peptide_maxes_for_protein = [(max(peptide_data_by_protein[protein][peptide]['intensity']), peptide) for peptide in peptide_data_by_protein[protein].iterkeys()]
			peptide_maxes_for_protein.sort()
			top_3 = peptide_maxes_for_protein[-3:]
			average_top_3 = averageList([x[0] for x in top_3])

			csv_writer.writerow([protein, average_top_3, top_3])

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("data_infile", help = "CSV file containing analyzed output", type = str)
	parser.add_argument("data_outfile", help = "CSV file containing raw input", type = str)
	parser.add_argument("PSM", help = "Boolean toggle for ignoring special characters (1 = peptide spectral match)", type = int)

	args = parser.parse_args().__dict__

	main(args["data_infile"], args["data_outfile"], args["PSM"])