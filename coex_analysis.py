### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python coex_analysis.py
					--exp <EXPRESSION_FILE>
					--anno <ANNOTATION_FILE>
					--out <OUTPUT_FOLDER>
					--candidates <CANDIDATE_GENE_FILE>
					
					optional:
					--minexp <MINIMAL_EXPRESSION_PER_GENE>[10]
					--r <MIN_CORRELATION_COEFFICIENT>[0.7]
					--p <MAX_P_VALUE>[0.05]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

from operator import itemgetter
import numpy as np
import re, math, os, sys
from scipy import stats

# ---- end of imports --- #

def load_expression_values( filename ):
	"""! @brief load all expression values """
	
	expression_data = {}
	with open( filename, "r" ) as f:
		tissues = f.readline().strip().split('\t')
		line = f.readline()
		if len( line.strip().split('\t') ) > len( tissues ):	#no header of the first column present
			tissues = tissues[1:]
		while line:
			parts = line.strip().split('\t')
			expression = {}
			for idx, each in enumerate( parts[1:] ):
				expression.update( { tissues[  idx ] : float( parts[ idx+1 ] ) } )
			line = f.readline()
			expression_data.update( { parts[0]: expression } )
	return expression_data


def compare_candidates_against_all( candidate, gene_expression, r_cut, p_cut, min_exp ):
	"""! @brief compare candidate gene expression against all genes to find co-expressed genes """
	
	tissues = sorted( gene_expression[ gene_expression.keys()[0] ].keys() )
	coexpressed_genes = []
	for i, gene2 in enumerate( gene_expression.keys() ):
		if candidate != gene2:
			values = []
			total_expression = 0
			for tissue in tissues:
				try:
					x = gene_expression[ candidate ][ tissue ]
					y = gene_expression[ gene2 ][ tissue ]
					total_expression += y
					if not math.isnan( x ) and not math.isnan( y ) :
						values.append( [ x, y ] )
				except KeyError:
					print tissue
			r, p = stats.spearmanr( values )
			if not math.isnan( r ) and total_expression > min_exp:
				if abs( r ) > r_cut and p < p_cut:
					coexpressed_genes.append( { 'id': gene2, 'correlation': r, 'p_value': p } )
		
	return coexpressed_genes


def load_annotation( annotation_file ):
	"""! @brief load annotation mapping table """
	
	annotation_mapping_table = {}
	
	with open( annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			annotation_mapping_table.update( { parts[0]: parts[1] } )
			line = f.readline()
	return annotation_mapping_table


def main( arguments ):
	"""! @brief run everything """
	
	expression_file = arguments[ arguments.index('--exp')+1 ]
	annotation_file = arguments[ arguments.index('--anno')+1 ]
	
	output_prefix = arguments[ arguments.index('--out')+1 ]
	candidate_gene_file = arguments[ arguments.index('--candidates')+1 ]
	
	if '--r' in arguments:
		r_cut = float( arguments[ arguments.index('--r')+1 ] )
	else:
		r_cut = 0.7
	
	if '--p' in arguments:
		p_cut = float( arguments[ arguments.index('--p')+1 ] )
	else:
		p_cut = 0.05
	
	if '--minexp' in arguments:
		min_exp = float( arguments[ arguments.index('--minexp')+1 ] )
	else:
		min_exp=10.0
	
	if not os.path.exists( output_prefix ):
		os.makedirs( output_prefix )
	else:
		print output_prefix
		
	annotation_mapping_table = load_annotation(  annotation_file)
	gene_expression = load_expression_values( expression_file )
	
	high_impact_candidates = [ ]
	
	with open( candidate_gene_file, "r" ) as f:
		line = f.readline()
		while line:
			gene = line.strip()
			high_impact_candidates.append( gene )
			line = f.readline()
	
	high_impact_candidates = list( set( high_impact_candidates ) )
	
	number_of_genes = float( len( gene_expression.keys() ) )
	for x, candidate in enumerate( high_impact_candidates ):
		coexpression_output_file = output_prefix + candidate + ".txt"
		if not os.path.isfile( coexpression_output_file ):
			with open( coexpression_output_file, "w", 0 ) as out:
				print str( x+1 ) + "/" + str( len( high_impact_candidates ) )
				coexp_values = compare_candidates_against_all( candidate, gene_expression, r_cut, p_cut, min_exp )
				pos_cor_vals = []
				neg_cor_vals = []
				for val in coexp_values:
					if val['correlation'] > 0:
						pos_cor_vals.append( val )
					else:
						neg_cor_vals.append( val )
				coexpressed_genes = sorted( pos_cor_vals, key=itemgetter( 'correlation' ) )[::-1] + sorted( neg_cor_vals, key=itemgetter( 'correlation' ) )[::-1]
				
				out.write( 'CandidateGene\tGeneID\tSpearmanCorrelation\tadjusted_p-value\tFunctionalAnnotation\n' )
				for entry in coexpressed_genes:
					try:
						out.write( "\t".join( map( str, [ candidate, entry['id'], entry['correlation'], entry['p_value'] * number_of_genes, annotation_mapping_table[ entry['id'] ] ] ) ) + '\n' )
					except KeyError:
						out.write( "\t".join( map( str, [ candidate, entry['id'], entry['correlation'], entry['p_value'] * number_of_genes, "N/A" ] ) ) + '\n' )


if "--exp" in sys.argv and '--anno' in sys.argv and '--out' in sys.argv and '--candidates' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
