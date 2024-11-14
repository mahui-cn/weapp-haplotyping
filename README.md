# weapp-haplotyping

Y-DNA and mtDNA Haplotyping for WeGene, It features DFS algorithm and specific rules to search the most appropriate terminal haplotype based on haplogroup tree which should contain SNP information with position and derived mutation.

The scripts should be used by weapp of WeGene (https://github.com/wegene-llc/weapp-developer-guide), which will let users to authorize to share their personal genome data. Once authorized, user's data will be sent to the script as input and final output will be used by WeGene to show in web page.

run: `python main.py < data/data.json`

online version: https://www.wegene.com/crowdsourcing/details/1265