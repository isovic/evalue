#! /usr/bin/python

import sys;
import math;

# This is the description taken from BLAST's source code:
# /** Supported substitution and gap costs with corresponding quality values
#  * for nucleotide sequence comparisons.
#  * NB: the values 0 and 0 for the gap costs are treated as the defaults used for
#  * the greedy gapped extension, i.e. 
#  * gap opening = 0, 
#  * gap extension = 1/2 match - mismatch.
#  * 
#  * The fields are:
#  * 
#  * 1. Gap opening cost,
#  * 2. Gap extension cost,
#  * 3. Lambda,
#  * 4. K,
#  * 5. H,
#  * 6. Alpha,
#  * 7. Beta,
#  * 8. Theta
#  */
class EvalueParams:
	def __init__(self, values=None):
		if (values == None):
			self.match = 0;
			self.mismatch = 0;
			self.gap_open = 0;
			self.gap_extend = 0;
			self.Lambda = 0;
			self.K = 0;
			self.H = 0;
			self.Alpha = 0;
			self.Beta = 0;
			self.Theta = 0;
		else:
			self.match = values[0];
			self.mismatch = values[1];
			self.gap_open = values[2];
			self.gap_extend = values[3];
			self.Lambda = values[4];
			self.K = values[5];
			self.H = values[6];
			self.Alpha = values[7];
			self.Beta = values[8];
			self.Theta = values[9];
			self.values = values;



# Returns a dict of pre-calculated E-value parameters, parsed from BLAST's source code.
# The keys of the dict are tuples of 4 values: (match, mismatch, gap open, gap extend).
# The value for each key consists of 
def load_blast_params(blast_stat_path):
	try:
		fp = open(blast_stat_path);
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % blast_stat_path);
		return;

	lines = fp.readlines();
	fp.close();

	evalue_params_dict = {};

	i = 0;
	while (i < len(lines)):
		line = lines[i].strip();
		if ('static const array_of_8 blastn_values_' in line):
			# The line is formatted as: static const array_of_8 blastn_values_1_2[] = {
			# Last two numbers are match and mismatch values.
			split_line = line.split('static const array_of_8 blastn_values_')[1].split('_');
			match = float(split_line[0]);
			mismatch = float(split_line[1].split('[')[0].strip());

			j = i + 1;
			while (('};' in lines[j]) == False):
				line_j = lines[j].strip();
				values_string = line_j.split('}')[0].split('{')[1].strip();
				values = [match, mismatch] + [float(value) for value in values_string.split(',')];
				current_evalue_params = EvalueParams(values);
				key = tuple([int(value) for value in values[0:4]]);
				evalue_params_dict[key] = current_evalue_params;
				j += 1;
		i += 1;

	# print evalue_params_dict[(1, 4, 1, 2)].values;
	# print '\n'.join([str(key) for key in sorted(evalue_params_dict.keys())]);
	return evalue_params_dict;

def calc_evalue(evalue_params_dict, match, mismatch, gap_open, gap_extend, query_length, target_length, alignment_score):
	try:
		evalue_params = evalue_params_dict[(match, mismatch, gap_open, gap_extend)];
	except:
		sys.stderr.write('ERROR: Could not find parameters for combination: match = %d, mismatch = %d, gap_open = %d, gap_close = %d!\n' % (match, mismatch, gap_open, gap_close));
		return 1000000000.0;

	return evalue_params.K * float(query_length) * float(target_length) * math.exp(-evalue_params.Lambda * float(alignment_score) );

def main():
	evalue_params_dict = load_blast_params('blast_stat.c');
	[match, mismatch, gap_open, gap_extend] = [5, 4, 8, 6];
	query_length = 20000;
	target_length = 500000000;
	alignment_score = 1400;

	sys.stdout.write('Running example for E-value calculation of nucleotide alignments.\n');
	sys.stdout.write('Using parameters:\n');
	sys.stdout.write('alignment_score = %d\n' % alignment_score);
	sys.stdout.write('query_length = %d\n' % query_length);
	sys.stdout.write('target_length = %d\n' % target_length);
	sys.stdout.write('match = %d, mismatch = %d, gap_open = %d, gap_extend = %d\n' % (match, mismatch, gap_open, gap_extend));
	sys.stdout.write('\n');

	evalue = calc_evalue(evalue_params_dict, match, mismatch, gap_open, gap_extend, query_length, target_length, alignment_score);
	sys.stdout.write('evalue = %s\n' % (str(evalue)) );

if __name__ == "__main__":
	main();