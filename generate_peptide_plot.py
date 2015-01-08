import subprocess
import os
import argparse
import csv
from collections import OrderedDict

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np

NCOL = 3

def main(input_directory, output_directory):
	# Get all files in input directory
	input_files = os.listdir(input_directory)
	input_files = [x for x in input_files if x[0] != '.']

	concentrations = [float(x[:x.find('f')]) for x in input_files]
	if len([x for x in input_files if x.find('f') == -1]):
		raise Exception, 'Could not use f to find fmol in this file!\n'

	peptide_dicts_at_conc = []
	for input_file in input_files:
		print 'loading file {}'.format(input_file)
		with open(os.path.join(input_directory, input_file), 'rU') as csvfile:
			dict_reader = csv.DictReader(csvfile)
			peptide_dict = {}
			for row in dict_reader:
				peptide = row["peptide"]
				area = float(row["mean area"])
				peptide_id = row["peptide id"]

				if peptide_dict.has_key(peptide):
					raise Exception, "Peptide already in dict!\n"
				else:
					peptide_dict[peptide] = {'area' : area, 'id' : peptide_id}
			peptide_dicts_at_conc.append(peptide_dict)

	peptide_set_list = []
	for peptide_dict in peptide_dicts_at_conc:
		peptide_set_list.append(set(peptide_dict.keys()))
	peptide_set = list(set.intersection(*peptide_set_list))

	output_dict = dict(
		[(peptide, {'area' : [], 'id' : []}) for peptide in peptide_set]
			)
	for peptide_dict_at_conc in peptide_dicts_at_conc:
		for peptide in peptide_set:
			output_dict[peptide]['area'].append(peptide_dict_at_conc[peptide]['area'])

			if not len(output_dict[peptide]['id']):
				output_dict[peptide]['id'].append(peptide_dict_at_conc[peptide]['id'])
	# Sort output_dict so it always prints in same order (by peptide ID not peptide sequence)
	output_dict = OrderedDict(sorted(output_dict.iteritems(), key=lambda (k,v): (v['id'],k)))

	# Write plot file
	plt.figure(figsize = (8.5, 11))

	number_rows = np.ceil(len(output_dict.keys()) / float(NCOL))

	for i, peptide in enumerate(output_dict.keys()):
		area = output_dict[peptide]['area']

		# Plot data
		axis = plt.subplot(number_rows, NCOL, i+1)
		axis.plot(concentrations, area, 'x')

		# Calculate fit
		coeffs, res, sing_val, rank, shape = np.polyfit(np.array(concentrations), np.array(area), deg = 1, full = True)
		fit = np.poly1d(coeffs)
		m = coeffs[0]
		b = coeffs[1]

		# Calculate R^2 of fit
		yhat = fit(concentrations)
		ybar = np.sum(area)/len(area)
		ssreg = np.sum((yhat-ybar)**2)
		sstot = np.sum((area - ybar)**2)
		r_squared = ssreg / sstot

	    # Plot fit
		fit_area = fit(concentrations)
		axis.plot(concentrations, fit_area, '-')
		axis.text(max(concentrations) - 0.05 * max(concentrations), max(area) - 0.05 * max(area), 'm={0}\nb={1}\nR^2={2:f}'.format(m,b,r_squared), fontsize=8, horizontalalignment='right', verticalalignment='top')
		# Format axis labels
		locs, labels = plt.xticks()
		plt.setp(labels, rotation=90)
		plt.xlabel(output_dict[peptide]['id'][0] + '\nConcentration (fmol)')
		plt.ylabel('Area')

	plt.subplots_adjust(hspace = 0.7, wspace = 0.7)
	plt.savefig(os.path.join(output_directory, 'standard_curve.pdf'))

	# Write CSV file
	with open(os.path.join(output_directory, 'complied_output.csv'), 'wb') as csv_outfile:
		csv_writer = csv.writer(csv_outfile)

		headers = [x for x in concentrations]
		headers.append('id')
		headers.append('peptide')
		csv_writer.writerow(headers)
		units = ['fmol']*len(concentrations)
		csv_writer.writerow(units)

		for peptide in output_dict.iterkeys():
			row = output_dict[peptide]['area']
			row.append(output_dict[peptide]['id'][0])
			row.append(peptide)
			csv_writer.writerow(row)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("input_directory", type = str)
	parser.add_argument("output_directory", type = str)

	args = parser.parse_args().__dict__

	main(args["input_directory"], args["output_directory"])