import subprocess
import os
import argparse

def main(input_directory, output_directory, PSM):
	# Get all files in input directory
	input_files = os.listdir(input_directory)
	input_files = [x for x in input_files if x[0] != '.']
	import ipdb; ipdb.set_trace()
	for input_file in input_files:
		subprocess.call(
			["python",
			"peptidestandard_peptide.py",
			os.path.join(input_directory,input_file),
			os.path.join(output_directory, input_file[:-4] + "_output.csv"),
			str(PSM)])

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("input_directory", type = str)
	parser.add_argument("output_directory", type = str)
	parser.add_argument("PSM", help = "Boolean toggle for ignoring special characters (1 = peptide spectral match)", type = int)

	args = parser.parse_args().__dict__

	main(args["input_directory"], args["output_directory"], args["PSM"])