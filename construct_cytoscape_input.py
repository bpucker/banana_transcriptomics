### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python construct_cytoscape_input.py
					--genes <GENE_FILE>
					--coexp <COEXP_DATA_FOLDER>
					--out <OUTPUT_FILE>
					
					optional:
					--tf <TF_CLASSIFICATION_FILE>
					--r <MINIMAL_CORRELATION>[0.7]
					--p <MAXIMAL_P_VALUE>[0.05]
					"""

import os, sys, glob

# --- end of imports --- #

def load_genes( gene_file ):
	"""! @brief load genes from given file """
	
	genes = {}
	with open( gene_file, "r" ) as f:
		line = f.readline()
		while line:
			gene = line.strip()
			if "\t" in gene:
				parts = gene.split('\t')
				genes.update( { parts[0]: parts[1] } )
			else:
				genes.update( { gene: gene } )
			line = f.readline()
	return genes


def load_coexp_data( coexp_file_folder, genes, TFs_genes, r_cutoff, p_cutoff ):
	"""! @brief load coexpression relationships between genes of interest """
	
	coexp = []
	black = {}
	filenames = glob.glob( coexp_file_folder + "*.txt" )
	for filename in filenames:
		with open( filename, "r" ) as f:
			f.readline()	#remove header
			line = f.readline()
			while line:
				parts = line.strip().split('\t')
				if float( parts[2] ) > r_cutoff and float( parts[3] ) < p_cutoff:
					status = 0
					try:
						genes[ parts[0] ]
						status += 1
					except KeyError:
						pass
					try:
						TFs_genes[ parts[1] ]
						status += 1
					except KeyError:
						pass
					if status > 1:
						try:
							black[ parts[1] + "_%_" + parts[0] ]
						except KeyError:
							coexp.append( [ parts[0], parts[1], parts[2] ] )
							black.update( { parts[0] + "_%_" + parts[1]: None } )
				line = f.readline()
	return coexp


def load_TF_classification( TF_classification_file ):
	"""! @brief load TF classification """
	
	TF_classification = {}
	with open( TF_classification_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if "bHLH" in parts[1]:
				TF_classification.update( { parts[0]: "bHLH" } )
			elif "WRKY" in parts[1]:
				TF_classification.update( { parts[0]: "WRKY" } )
			elif "NAC" in parts[1]:
				TF_classification.update( { parts[0]: "NAC" } )
			elif "MYB" in parts[1]:
				TF_classification.update( { parts[0]: "MYB" } )
			elif "MADS" in parts[1]:
				TF_classification.update( { parts[0]: "MADS" } )
			elif "Hsf" in parts[1]:
				TF_classification.update( { parts[0]: "Hsf" } )
			elif "bZIP" in parts[1]:
				TF_classification.update( { parts[0]: "bZIP" } )
			elif "ARF" in parts[1]:
				TF_classification.update( { parts[0]: "ARF" } )
			line = f.readline()
	return TF_classification


def main( arguments ):
	"""! @brief run everything """
	
	gene_file = arguments[ arguments.index('--genes')+1 ]
	coexp_file_folder = arguments[ arguments.index('--coexp')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	
	if '--tf' in arguments:
		TF_classification_file = arguments[ arguments.index('--tf')+1 ]
		TF_classification = load_TF_classification( TF_classification_file )
	else:
		TF_classification = {}
	
	if '--r' in arguments:
		r_cutoff = float( arguments[ arguments.index('--r')+1 ] )
	else:
		r_cutoff = 0.7
	
	if '--p' in arguments:
		p_cutoff = float( arguments[ arguments.index('--p')+1 ] )
	else:
		p_cutoff = 0.05


	genes = load_genes( gene_file )
	print "number of genes of interest: " + str( len( genes.keys() ) )

	coexp_data = load_coexp_data( coexp_file_folder, genes, TF_classification, r_cutoff, p_cutoff )

	with open( output_file, "w" ) as out:
		for entry in coexp_data:
			new_line = []
			try:
				new_line.append( genes[ entry[0] ] )
			except KeyError:
				new_line.append( entry[0] )
			try:
				new_line.append( genes[ entry[1] ] )	#TF should be the second ID in pairs
			except KeyError:
				new_line.append( entry[1] )
			new_line.append( entry[2] )
			try:
				new_line.append( TF_classification[ entry[1] ] )
			except KeyError:
				new_line.append( "other" )
			out.write( "\t".join( new_line ) + "\n" )


if '--genes' in sys.argv and '--coexp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
