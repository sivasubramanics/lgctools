"""
Few default variables for the lgctools
"""

MISSING_CALLS = ('N:N', '?:?', "?", "Uncallable",
                 "Unused", "missing", "Empty", "NTC")
BLANK_SAMPLES = ('NTC', 'Empty', "")
CSV = ','
TAB = '\t'
NEWLINE = '\n'
SPACE = ' '
ALLELE_SEP = ':'
CONSENSUSCUTOFF = 75
MISSING_ALLELES = ['?', 'Uncallable', 'Unused']
COLORS = ['#ff3333', '#4d4dff', '#00cc66', '#cc66ff', '#336699']
HMPHEAD = ("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID", "panelLSID",
           "QCcode")
HMPNA = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA']
MULTIPROCESSING = True
IUPAC = {'A':'A:A', 
         'T':'T:T', 
         'C':'C:C', 
         'G':'G:G', 
         'N':'N:N', 
         'R':'A:G', 
         'S':'C:G', 
         'M':'A:C', 
         'Y':'C:T', 
         'W':'A:T', 
         'K':'G:T'}
SAMPLE_MISSING_CUTOFF = 80
MARKER_MISSING_CUTOFF = 100
MARKER_MAF_CUTOFF = 0
