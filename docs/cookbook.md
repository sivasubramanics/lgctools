# lgctools Cookbook

From discussion: Siva, Victor

## Installation (Siva will verify prerequisites for install)

1. Clone the repository:  
`$ git clone https://github.com/sivasubramanics/lgctools`  
`$ cd lgctools`  
2. Run setup:  
`$ python3 setup.py build`  
3. Invoke lgctools:  
`$ lgctools`  
4. Display lgctools options:  
`$ lgctools --help`  
`usage: lgctools [-h] [--summary] [--rename] [--out-grid] [--out-flapjack] [--out-hapmap]`   
`                [--differences] [--performance] [--make-plots] [--bestmarkers] [--markercloud]`  
`                [--filter] [--pedver] [--consensus] [--fwdbreed] [--lgc-file <FILE>]`  
`                [--lgc-files <FILES> [<FILES> ...]] [--grid-file <FILE>] [--grid-files <FILES>]`  
`                [--hmp-file <FILE>] [--samplemap-file <FILE>] [--ped-file <FILE>]`  
`                [--designation-file <FILE>] [--meta-data <FILE>] [--marker-list <FILE>]`  
`                [--sample-list <FILE>] [--male-parents-list <FILE>] [--female-parents-list <FILE>]`  
`                [--qtl-file <FILE>] [--out <STR>] [--f1het-cutOff <INT>]`  
`                [--consensus-cutOff <INT>] [--max-missing-site <INT>] [--min-pic-site <FLOAT>]`  
`                [--min-maf-site <FLOAT>] [--max-missing-sample <INT>]`  
`QC pipeline: Processes the LGC file for tasks involved in purity check`   
`optional arguments:`  
`  -h, --help            show this help message and exit`  
`  --summary             Generate sample and marker summary for the data from genotype file`  
`  --rename              Rename the samples provided in the genotype file`  
`  --out-grid            Write genotype data in grid file format`  
`  --out-flapjack        Write genotype data in flapjack file format`  
`  --out-hapmap          Write genotype data in hapmap file format`  
`  --differences         Find polymorphic markers between genotypes`  
`  --performance         Check performance of the provided marker set`  
`  --make-plots          Create plots based on the LGC data`  
`  --bestmarkers         Find best markers from the given marker set`  
`  --markercloud         Based on best marker summary make word cloud`  
`  --filter              Filter genotype data`  
`  --pedver              Analyze genotype data for F1 Verification`  
`  --consensus           Call consensus and make purity reports`  
`  --fwdbreed            Analyze genotype data for favorable allele propotions`  
`  --lgc-file <FILE>     LGC raw data File`  
`  --lgc-files <FILES> [<FILES> ...]`  
`                        LGC raw data File`  
`  --grid-file <FILE>    LGC Grid Matrix File`  
`  --grid-files <FILES>  Comma seperated LGC Grid Matrix Files`  
`  --hmp-file <FILE>     Hapmap genotype File`  
`  --samplemap-file <FILE>`  
`                        Tab seperated sample map file. <SAMPLE_ID> <SAMPLE_NAME>`  
`  --ped-file <FILE>     Tab seperated Pedigree File. <F_ONE> <PAR_A> <PAR_B>`  
`  --designation-file <FILE>`  
`                        Tab seperated Designations parent File. <SAMPLE_NAME> <DESIGNATION>`  
`  --meta-data <FILE>    Parent information file`  
`  --marker-list <FILE>  File with list of snps to analyze`  
`  --sample-list <FILE>  File with list of samples to analyze`  
`  --male-parents-list <FILE>`  
`                        File with list of male parents to consider`  
`  --female-parents-list <FILE>`  
`                        File with list of female parents to consider`  
`  --qtl-file <FILE>     QTL file in GOBii format`  
`  --out <STR>           Output filename prefix`  
`  --f1het-cutOff <INT>  Percentage expected heterozygosity for F1 verification`  
`  --consensus-cutOff <INT>`  
`                        Percentage propotion to be considered to call consensus`  
`  --max-missing-site <INT>`  
`                        Percentage propotion to be considered for maximum missing per marker`  
`  --min-pic-site <FLOAT>`  
`                        Percentage propotion to be considered for PIC per marker`  
`  --min-maf-site <FLOAT>`  
`                        Percentage propotion to be considered for MAF per marker`  
`  --max-missing-sample <INT>`  
`                        Percentage propotion to be considered for maximum missing per sample`  

## Get summary information of an LGC file.

1. __Input:__  LGC file in csv format
2. __Command:__  `$ lgctools --summary --lgc-file <filename>`
3. __Ouput(s):__  `tsv files: out_marker_summary.txt, out_sample_summary.txt`

__Note:__  By default, the script automatically creates files with filenames prefixed with `out_`. You can specify a custom prefix to identify the output. The command now looks like and will create filenames with `my_output` prefix:  
`$ lgctools --summary --lgc-file <filename> --out my_output` 

## Rename sample names in genotype file

1. __Input:__ LGC file; Sample map file (see note below)  
2. __Command:__ `$ lgctools –-rename --samplemap-file <sample map file> –-lgc-file <lgc file> –-out my_out`  
3. __Output:__ grid file  
  
__Note:__ Sample map file is a two-column file where the first column is the LGC sample ID and the second column is the BMS sample name (unique ID).  
  
`LGC_SampleID    BMS_SampleName`  
`8T3QSA1xhsFY0   ICCV 181114:9`  
`8T3QSY1HUsiNz   ICCV 181114:10`  
`8T3QSC3xbel7z   ICCV 181114:11`  
`8T3QSNYuhxw1N   ICCV 181114:12`  
`8T3QSTpBTNk5n   ICCV 181114:13`  
`8T3QSABgN751a   ICCV 181114:14`  
`8T3QSgixO4FwT   ICCV 181114:15`  
`8T3QSk0EeAAzQ   ICCV 181114:16`  
`8T3QS0euXBPTF   ICCV 181117:9`  
`8T3QS83sNaLkE   ICCV 181117:10`  
`8T3QSZ6MthiVs   ICCV 181117:11`  
`8T3QSRzHV8meN   ICCV 181117:12`  
`8T3QS0ZrSVhZB   ICCV 181117:13`  
`8T3QS7MH6iQvz   ICCV 181117:14`  
`8T3QSHcS8awvu   ICCV 181117:15`  

## Convert an lgc file to hapmap (to do, create example for all compatible file formats)

1. __Input:__ LGC file  
2. __Command:__ `$ lgctools –-lgc-file <lgc file> –-out my_out –-out-hapmap`  
3. __Output:__ `my_out.hmp.txt`  
  
__Note:__ You can also convert an LGC file into grid and flapjack formats by replacing `--out-hapmap` with `--out-flapjack` or `--out-grid`. Similarly, you can convert a grid file into hapmap by replacing `--lgc-file` with `--grid-file` and corresponding output file: `--out-hapmap`.  

## Get all possible polymorphic markers in all pairwise genotypic combinations. (to do, include --prefix, create examples for other file formats hapmap, grid).

1. __Input:__ LGC or hapmap or grid file  
2. __Command:__ `$ lgctools –-lgc-file <lgc-file> –-differences –-out diff`  
3. __Output:__ `diff.txt`  
  
__Note:__ You can specify a hapmap or grid file by replacing `--lgc-file` with corresponding input type (`--hapmap-file` or `--grid-file`). Output is a table of markers for all pairwise genotypes combinations.  

__Discussion of output (Siva). Show screenshot of output and discuss the result.__  
