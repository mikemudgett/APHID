import sys, os
from Bio import SeqIO
import numpy as np
from Bio.Blast.Applications import NcbimakeblastdbCommandline

from aphid_utils import pull_guides


#This script prepares the necessary files to run APHID. It comprises of three steps:
#1. Download the genome files if they are not already in the current directory.
#2. Pull a list of all the guides in the genome, specific to your PAM sequences and guide length.
#   This will make it easier to search for possible off-targets.
#3. Create a BLAST database of the Arabidopsis genome. This will make finding specific primers easier.


args = sys.argv
if len(args) != 3:
	print('Input Error.')
	print('Usage: python aphid_prep.py {PAM Seq} {Guide Length}')
	print('Example: python aphid_prep.py NGG 20')
	sys.exit()

try:
	pam_seq = str(args[1])
	guide_length = int(args[2])
except:
	print('Input Error.')
	print('Usage: python aphid_prep.py {PAM Seq} {Guide Length}')
	print('Example: python aphid_prep.py NGG 20')
	sys.exit()

#pam_seq = 'NGG'
#guide_length = 20

#PART 1: Genome Download

#Download individual chromosome files if they are not already in the directory.
#We are using Arabidopsis thaliana release 56 from ensembl plants. This could be updated with new releases as long as the naming conventions stay the same.
print('Downloading genome and annotation files...')
for num in ['1', '2', '3', '4', '5', 'Mt', 'Pt']:
	fn = 'Arabidopsis_thaliana.TAIR10.dna.chromosome.' + num + '.fa'
	if not os.path.exists(fn):
		os.system('wget -nc https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/fasta/arabidopsis_thaliana/dna/' + fn + '.gz')
		os.system('gunzip ' + fn)
	fn = 'Arabidopsis_thaliana.TAIR10.56.chromosome.' + num + '.gff3'
	if not os.path.exists(fn):
		os.system('wget -nc https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-56/plants/gff3/arabidopsis_thaliana/' + fn + '.gz')
		os.system('gunzip ' + fn)
print('Genome files are downloaded.')
print('')

#PART 2: Guide Search

#this is pretty slow, takes ~7 minutes for NGG and 20bp guides
print('Creating list of all nuclear guide sites with PAM ' + pam_seq + ' of length ' + str(guide_length) + '...')

out_fn = 'all_guides_' + pam_seq + '_' + str(guide_length) + '.txt'

if os.path.exists(out_fn):
	print('Guide list already exists in directory. Moving on to creating BLAST database.')
else:
	out_file = open(out_fn, 'w')
	print('Scanning genome for guides.')
	for chr_num in (1,2,3,4,5):
		with open('Arabidopsis_thaliana.TAIR10.dna.chromosome.' + str(chr_num) + '.fa') as fa:
			for record in SeqIO.parse(fa, "fasta"):
				chr_seq = str(record.seq)
				chr_guides = pull_guides(chr_seq, pam_seq, guide_length)
				for guide in chr_guides:
					out_file.write(str(guide) + '\n')
		print('Completed scanning chromosome ' + str(chr_num) + '.')
		
	out_file.close()
	print('Guide list complete.')
print('')




#PART 3: BLAST Database
print('Building BLAST database...')


if os.path.exists('Arabidopsis_thaliana.TAIR10.dna.full.fa'):
	print('Concatenated genome already in directory.')
else:
	print('Concatenating genome.')
	os.system('cat Arabidopsis_thaliana.TAIR10.dna.chromosome.*.fa > Arabidopsis_thaliana.TAIR10.dna.full.fa')

if os.path.exists('athal.nhr'):
	print('BLAST database already exists in directory.')
else:
	print('Generating database.')
	makedb_cline = NcbimakeblastdbCommandline(dbtype='nucl', input_file='Arabidopsis_thaliana.TAIR10.dna.full.fa', out='athal')
	makedb_cline()

print('')
print('Setup complete. You may move on to using APHID.')

