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
IUPAC = {'A': 'A:A',
         'T': 'T:T',
         'C': 'C:C',
         'G': 'G:G',
         'N': 'N:N',
         'R': 'A:G',
         'S': 'C:G',
         'M': 'A:C',
         'Y': 'C:T',
         'W': 'A:T',
         'K': 'G:T'}
SAMPLE_MISSING_MAXIMUM = 80
MARKER_MISSING_MAXIMUM = 100
MARKER_PIC_MINIMUM = 0
MARKER_MAF_MINIMUM = 0
CUTOFF_FONE = 60
CUTOFF_CONSENSUS = 75
TASKS = ['convert', 'reheader', 'makeplots', 'summary',
         'merge', 'differences', 'performance', 'bestmarkers', ]

FORMATS = ['lgc', 'grid', 'hmp', 'fjk']
USAGE = """
Program: lgctools (Tools for processing LGC files)
Version: 0.1

Usage:   lgctools <task> [options]

Tasks:
  -- File Operations
     covert              convert different file formats
     reheader            change the sample name/id from the given genotype file
     makeplots           generate genotype call plots
     summary             generate sample and marker summary for the data from genotype file
     merge               merge multiple genotype data files

  -- Performance
     differences         find polymorphic markers between genoytpes
     performance         check performance of the provided marker set
     bestmarkers         find best markers from the given marker set

  -- Parental Purity
     checkparents        returns the list of parents for which the data is available in
     consensus           call consensus for the replicated genotypes (parents)
     filter              filter or remove parent replicates from the genotype files

  -- Pedigree Verification
     pedver              analyze genotype data for F1 Verification
     
  -- Forward Breeding
     fwdbreed            analyze genotype data for favorable allele propotions
     
"""
