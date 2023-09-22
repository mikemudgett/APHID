import primer3
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
import numpy as np
import regex as re
from Bio.Seq import Seq
from io import StringIO

class Gene:
	def __init__(self):
		self.coding_exons = []
		self.five_prime_UTR = []
		self.three_prime_UTR = []
		
	def set_polarity(self, strand):
		if strand == '+':
			self.polarity = 0
		else:
			self.polarity = -1
			
	def set_coords(self, start, end):
		self.mRNA_start = start
		self.mRNA_end = end
		
	def add_feature(self, line):
		fields = line.split('\t')
		match fields[2]:
			case 'mRNA':
				self.set_polarity(fields[6])
				self.set_coords(int(fields[3]), int(fields[4]))
			case 'exon':
				pass
			case 'CDS':
				#Start, End, Length
				#Assumes the exons are all in order
				self.coding_exons.append((int(fields[3]), int(fields[4]), int(fields[4]) - int(fields[3])+1))
			case 'five_prime_UTR':
				self.five_prime_UTR.append((int(fields[3]), int(fields[4])))
			case 'three_prime_UTR':
				self.three_prime_UTR.append((int(fields[3]), int(fields[4])))
			case _:
				pass

#DOES NOT give the reverse complement if you have '-' direction.
def grow_primer_to_tm(seq, start_pos, direction, min_Tm=59, min_length=15):
	melt_temp = 0
	increment = min_length
	start = start_pos
	while melt_temp < min_Tm:
		if direction == '+':
			primer = seq[start_pos:start+increment+1]
			melt_temp = primer3.calc_tm(str(primer))
		else:
			primer = seq[start+1-increment:start+1]
			melt_temp = primer3.calc_tm(str(primer))

		increment += 1
	#This returns a "Seq" not a string. A little sloppy.
	return primer


#BLAST each candidate we have at this point against the whole genome database.
#Maybe just select the candidate with the least number of significant hits?
#Make the final decision based on the longest homology arm.
def get_blast_hits(seq):

#for candidate in candidate_primer_list.copy():
	temp_fa = open('candidate.fa', 'w')
	temp_fa.write('>candidate' + '\n' + seq)
	temp_fa.close()
	blast_cline = NcbiblastnCommandline(query='candidate.fa',db='athal', evalue=10, task='blastn-short', outfmt='10 nident qstart qend')
	blast_result = pd.read_csv(StringIO(blast_cline()[0]))
	#print(len(blast_result))
	#return len(blast_result)
	#return blast_result.sum()[1]
	return blast_result
	

#Figure out the strandedness yet. I don't know if there is an issue here.
#Faster now.
#Everything might kind of be off by one. Might not really matter?
def choose_primer(template_seq, lower_bound, upper_bound, strand, opposite_primer='', min_Tm=55, max_Tm=62, min_length=20, max_length=27):

	primer_candidate_list = pd.DataFrame(columns = ['sequence', 'local_position', 'length', 'tm', 'hairpin_tm', 'homodimer_tm', 'heterodimer_tm', 'thermo_tm_sum', 'blast'])

	min_pos = lower_bound
	max_pos = upper_bound

	for length in range(min_length, max_length+1):
		if strand == '+':
			start = min_pos
			stop = max_pos
		elif strand == '-':
			start = min_pos - length
			stop = max_pos - length
		#Get BLAST hits for the entire range of possible primers of this length.
		blast_result = get_blast_hits(str(template_seq[start:stop + length]))
		for x in range(start, stop):
			cand_seq = str(template_seq[x:x+length])
			if strand == '-':
				cand_seq = str(Seq(cand_seq).reverse_complement())
			cand_tm = primer3.calc_tm(cand_seq)
			if (cand_tm > min_Tm) and (cand_tm < max_Tm):
				cand_hairpin_tm = primer3.calc_hairpin(cand_seq, dna_conc=200).tm
				cand_homodimer_tm = primer3.calc_homodimer(cand_seq, dna_conc=200).tm
				if opposite_primer == '':
					cand_heterodimer_tm = 0
				else:
					cand_heterodimer_tm = primer3.calc_heterodimer(cand_seq, opposite_primer, dna_conc=200).tm
				#Don't need to bother with meeting a threshold if we aren't blasting every single one of them.
				#if (cand_hairpin_tm < thermo_thresh) and (cand_homodimer_tm < thermo_thresh) and (cand_heterodimer_tm < thermo_thresh):
					#cand_blast = get_blast_hits(cand_seq)


				#Now we just have to figure out how to get the number of blast hits that match a given primer
				cand_blast = 0
				for entry in blast_result.index:
					#print(str(x) + '  ' + str(blast_result.iloc[entry,1]))

					if ((x+1-start)<=blast_result.iloc[entry,1]) and ((x+length-start)>=blast_result.iloc[entry,2]):
						cand_blast += 1


				primer_candidate_list = pd.concat([primer_candidate_list, pd.Series({'sequence': cand_seq, 
																					 'local_position': int(x), 
																					 'length': int(length), 
																					 'tm': cand_tm, 
																					 'hairpin_tm': cand_hairpin_tm, 
																					 'homodimer_tm': cand_homodimer_tm, 
																					 'heterodimer_tm': cand_heterodimer_tm,
																					 'thermo_tm_sum': cand_hairpin_tm + cand_homodimer_tm + cand_heterodimer_tm,
																					 'blast': cand_blast}).to_frame().T], ignore_index = True,)
					
	sort_list = ['blast', 'thermo_tm_sum', 'length', 'local_position']
	if strand == '+':
		sort_list_order = [True, True, True, True]
	elif strand == '-':
		sort_list_order = [True, True, True, False]

	primer_candidate_list = primer_candidate_list.sort_values(sort_list, ascending=sort_list_order)

	return primer_candidate_list.iloc[0]
	#this is a pandas series with all the info generated about the chosen primer


def pull_guides(chr_seq, pam_seq, guide_length):
	IUPAC_dict = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'U':'T',
				  'R':'[AG]', 'Y':'[CT]', 'S':'[GC]', 'W':'[AT]',
				  'K':'[GT]', 'M':'[AC]', 'B':'[CGT]', 'D':'[AGT]',
				  'H':'[ACT]', 'V':'[ACG]', 'N':'[ACGT]', 'X':'[ACGT]'}

	guide_list = []

	regex_token = '.{' + str(guide_length) + '}'
	for x in pam_seq:
		regex_token += IUPAC_dict[x]


	rc_pam_seq = Seq(pam_seq).reverse_complement()
	rc_regex_token = ''
	for x in rc_pam_seq:
		rc_regex_token += IUPAC_dict[x]
	rc_regex_token += '.{' + str(guide_length) + '}'

	sites = re.findall(regex_token, chr_seq, overlapped=True)
	for guide in sites:
		guide_list.append(guide[0:guide_length])
		
	rc_sites = re.findall(rc_regex_token, chr_seq, overlapped=True)
	for rc_guide in rc_sites:
		guide_list.append(str(Seq(rc_guide).reverse_complement())[0:guide_length])

	return guide_list