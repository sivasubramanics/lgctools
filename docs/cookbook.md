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
`                [--differences] [--performance] [--make-plots] [--bestmarkers] [--markercloud]`  
`                [--filter] [--pedver] [--consensus] [--fwdbreed] [--lgc-file <FILE>]`  
`                [--lgc-files <FILES> [<FILES> ...]] [--grid-file <FILE>] [--grid-files <FILES>]`  
`                [--hmp-file <FILE>] [--samplemap-file <FILE>] [--ped-file <FILE>]`  
`                [--designation-file <FILE>] [--meta-data <FILE>] [--marker-list <FILE>]`  
`                [--sample-list <FILE>] [--male-parents-list <FILE>] [--female-parents-list <FILE>]`  
`                [--qtl-file <FILE>] [--out <STR>] [--f1het-cutOff <INT>]`  
`                [--consensus-cutOff <INT>] [--max-missing-site <INT>] [--min-pic-site <FLOAT>]`  
`                [--min-maf-site <FLOAT>] [--max-missing-sample <INT>]`  <br />
`QC pipeline: Processes the LGC file for tasks involved in purity check`   
`optional arguments:`  
`  -h, --help            show this help message and exit`  
`  --summary             Generate sample and marker summary for the data from genotype file`  
`  --rename              Rename the samples provided in the genotype file`  
`  --out-grid            Write genotype data in grid file format`  
`  --out-flapjack        Write genotype data in flapjack file format`  
`  --out-hapmap          Write genotype data in hapmap file format`  
`  --differences         Find polymorphic markers between genotypes`  
`  --performance         Check performance of the provided marker set`  
`  --make-plots          Create plots based on the LGC data`  
`  --bestmarkers         Find best markers from the given marker set`  
`  --markercloud         Based on best marker summary make word cloud`  
`  --filter              Filter genotype data`  
`  --pedver              Analyze genotype data for F1 Verification`  
`  --consensus           Call consensus and make purity reports`  
`  --fwdbreed            Analyze genotype data for favorable allele propotions`  
`  --lgc-file <FILE>     LGC raw data File`  
`  --lgc-files <FILES> [<FILES> ...]`  
`                        LGC raw data File`  
`  --grid-file <FILE>    LGC Grid Matrix File`  
`  --grid-files <FILES>  Comma seperated LGC Grid Matrix Files`  
`  --hmp-file <FILE>     Hapmap genotype File`  
`  --samplemap-file <FILE>`  
`                        Tab seperated sample map file. <SAMPLE_ID> <SAMPLE_NAME>`  
`  --ped-file <FILE>     Tab seperated Pedigree File. <F_ONE> <PAR_A> <PAR_B>`  
`  --designation-file <FILE>`  
`                        Tab seperated Designations parent File. <SAMPLE_NAME> <DESIGNATION>`  
`  --meta-data <FILE>    Parent information file`  
`  --marker-list <FILE>  File with list of snps to analyze`  
`  --sample-list <FILE>  File with list of samples to analyze`  
`  --male-parents-list <FILE>`  
`                        File with list of male parents to consider`  
`  --female-parents-list <FILE>`  
`                        File with list of female parents to consider`  
`  --qtl-file <FILE>     QTL file in GOBii format`  
`  --out <STR>           Output filename prefix`  
`  --f1het-cutOff <INT>  Percentage expected heterozygosity for F1 verification`  
`  --consensus-cutOff <INT>`  
`                        Percentage propotion to be considered to call consensus`  
`  --max-missing-site <INT>`  
`                        Percentage propotion to be considered for maximum missing per marker`  
`  --min-pic-site <FLOAT>`  
`                        Percentage propotion to be considered for PIC per marker`  
`  --min-maf-site <FLOAT>`  
`                        Percentage propotion to be considered for MAF per marker`  
`  --max-missing-sample <INT>`  
`                        Percentage propotion to be considered for maximum missing per sample`  

