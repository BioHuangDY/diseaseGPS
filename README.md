# diseaseGPS

A bioinformatics software tool for diagnosis of genetic diseases based on genes and phenotypes.

## SYNOPSIS

diseasegps.py [options]

## WHAT DOES IT DO

diseaseGPS is a python script for diagnosis of genetic diseases based on genes and phenotypes.

## FOLDER structure

.
├── annovardb
├── intervardb
├── humandb
├── test
├── example
├── gpsdb
├─── config_gps.ini
├─── README.md
└─── diseasegps.py

## PREREQUISITE

1. You need install Python >=3.7.10.

2. You need install Perl >=5.22.1

3. You need install gcc version >=5.4.0

4. You need install [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) version >= 2018-04-16.  Put the three files(annotate_variation.pl, convert2annovar.pl, table_annovar.pl) into the folder annovardb

5. You need install [INTERVAR](https://github.com/WGLab/InterVar#intervar) version >= 2019-03-27.  Put all the files into the folder intervardb

6. You need download dataset for annotation from [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) and put all the dataset into the folder humandb. You can use the code by following:

   ```
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene  ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar esp6500siv2_all  ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar 1000g2015aug_all   ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar exac03    ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar gnomad211_genome    ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar dbnsfp42a ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar dbscsnv11 ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar avsnp150 ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar clinvar_20210501       ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar dbnsfp31a_interpro      ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar ensGene   ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar dbnsfp31a_interpro      ./humandb/
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar knownGene ./humandb/
   ```

   You need download the datasets and their idx as follows:

   - refGene 
   - esp6500siv2_all 
   - 1000g2015aug_all 
   - exac03 
   - gnomad211_genome 
   - dbnsfp42a 
   - dbscsnv11 
   - avsnp150 
   - clinvar_20210501 
   - dbnsfp31a_interpro 
   - ensGene 
   - knownGene

   except the dataset rmsk, you can use the code by following:

   ```
   ./annovardb/annotate_variation.pl -downdb -buildver hg19 rmsk ./humandb/
   ```

7. You need download other files such as hp.obo and gene_to_phenotype.txt from [HPO](https://hpo.jax.org/app/download/)

8. Please use the updated files(should be generated: >= 2021-09-13) from all dataset, outdated files will bring problems of InterVar.

## OPTIONS

#### Basic options

- **-h,  --help**            
  - show this help message and exit

* **--version**             
  * show program's version number and exit
* **-c config_gps.ini,  --config=config_gps.ini**
  * The config file of all options. it is for your own configure file.You can edit all the options in the
    configure.

if you use this options,you can ignore all the other options bellow.

* **-b hg19,  --buildver=hg19**

  * The genomic build version, it can be hg19 or hg38 , will support other version later 

* **-i test/testNjob1.vcf,  --input=test/testNjob1.vcf**

  * The input file contains your variants

* **-p test/testNjob1_hpo.txt,  --phenotype=test/testNjob1_hpo.txt**

  * The input file contains your hpoterms separated by "\n" , which means one hpo term per line.

* **-k test/testNjob1_hpo.json,  --phenopackets=test/testNjob1_hpo.json**

  * The phenopackets formatted phenotype input file. The file can only contain information for one patient and cannot be used simultaneously with the "-p" option.

* **--input_type=VCF**      

  * The input file type, it can be  AVinput(Annovar's format), VCF(VCF with single sample)

* **-o test/testNjob1,  --output=test/testNjob1**

  * The prefix of output file which contains the results, the file of results will be as *[prefix].gps*

  * In the meantime, you can get the intermediate results as *[prefix]\_gene.list*   and   *[prefix]\_hpo.list* 

    

#### Database location Options

* **-t gpsdb, --database_gps=gpsdb**
  * The  database location/dir for the gps dataset files

* **-d humandb, --database_locat=humandb**
  * The  database location/dir for the annotation datasets

* **-r intervardb, --database_intervar=intervardb**
  * The  database location/dir for the intervar dataset files
* **-a annovardb, --database_annovar=annovardb**
  * The  database location/dir for the annovar dataset files



#### Annovar Options

* **--table_annovar=./table_annovar.pl**
  * The Annovar perl script of table_annovar.pl

* **--convert2annovar=./convert2annovar.pl**
  * The Annovar perl script of convert2annovar.pl

* **--annotate_variation=./annotate_variation.pl**
  * The Annovar perl script of annotate_variation.pl

* **--skip_annovar**      
  * Skip the Annovar annotation, this can be true only after you  already got the annovar's annotation results



## EXAMPLE

```shell
./diseasegps.py -c config_gps.ini  # Run the examples in config_gps.ini
./diseasegps.py  -b hg19 -i your_vcf_input  --input_type=VCF -p your_phenotype_input  -o your_output
```

## HOW DOES IT WORK

DiseaseGPS takes vcf files containing gene mutation information and hpo files containing phenotype information as input, and takes the sorted result of OMIM disease diagnosis as output. The execution of DiseaseGPS mainly consists of three major steps: 1) Annotate mutation files to get the gene scores;  2) Analyze the phenotype information to obtain the phenotype score; and 3) Integrate phenotype score and gene score to get the final OMIM score ranking.  DiseaseGPS annotate mutations using the ANNOVAR tool, including refGene, esp6500siv2_all, 1000g2015aug_all, exac03, gnomad211_genome, dbnsfp42a, dbscsnv11, avsnp150, clinvar_20210501, dbnsfp31a_interpro, ensGene, knownGene and rmsk. Before performing pathogenicity analysis, we need to pre-screen the frequency of mutations in various population databases.  Users can choose to use any combination of 1000G, ExAC, ESP, gnomAD to pre-screen mutations. The threshold can be determined by the user, and the default is 0.1. According to the annotation file, DiseaseGPS use the InterVar tool to implement the ACMG2015 guidelines. Bases on the evidence, we can automatically recognize 18 indicators as following: PVS1, PS1, PS4, PM1, PM2, PM4, PM5, PP2, PP3, PP5, BA1, BS1, BS2, BP1, BP3, BP4,  BP6, BP7 and calculate a Bayesian posterior probability.
$$
OddsPath=O_{PVSt}^{(\frac{N_{PSu}}{8}+\frac{N_{PM}}{4}+\frac{N_{PSt}}{2}+\frac{N_{PVSt}}{1}-\frac{N_{BSu}}{8}-\frac{N_{BSt}}{1})}
$$

$$
Post\_P=\frac{OddsPath×Prior\_P}{((OddsPath-1)×Prior\_P+1)}
$$

Among them, N<sub>PSu</sub>, N<sub>PM</sub>, N<sub>PSt</sub>, N<sub>PVSt</sub>, N<sub>BSu</sub>, and N<sub>BSt</sub> are the number of supporting pathogenic evidence PP, the number of moderate pathogenic evidence PM, the number of strong pathogenic evidence PS, the number of very strong pathogenic evidence PVS, the number of supporting benign evidence BP, and the strong benign evidence BS, respectively. OPVSt is the pathogenic strength of strong evidence of pathogenicity, and the default is 350. OddsPath is the conditional pathogenicity probability of gene mutation site. Prior_P is the prior probability of the pathogenicity of the gene mutation site, and the default is 0.1. Post_P is the posterior probability of the pathogenicity of the gene mutation site, which is obtained by calculation. Post_P can be used to evaluate the pathogenicity of gene mutation sites, which is used as gene score.

The phenotype set entered by the user is marked as v<sub>1</sub>, and the annotated phenotype set of any OMIM disease is marked as v<sub>2</sub>.  Based on the hp.obo file, we can build a tree structure diagram of HPO. The phenotypic score is composed of TopSimilarity and MeicSimilarity.  TopSimilarity can be caculated by following:
$$
𝑇opSimilarity(𝑣_1,𝑣_2 )=\frac{𝑇𝑜𝑝(𝑣_1 \bigcap v_2)}{𝑇𝑜𝑝(𝑣_1)}
$$
Where Top(v) refers to the ancestors number of set v.

MeicSimilarity can be caculated by following:
$$
MeicSimilarity (𝑣_1,𝑣_2 )=\max\limits_{𝑝∈𝑣_1 \bigcap v_2}⁡𝐼𝐶_p ∗𝑎𝑛𝑐𝑒𝑠𝑡𝑜𝑟𝑠(𝑝)
$$

$$
𝐼𝐶_p=\ln⁡(\frac{𝑁}{𝑛_𝑝})
$$

Where N is the number of all OMIM diseases, n<sub>p</sub> is the phenotype annotated by phenotype p.

TopSimilarity and MeicSimilarity can be added by a cutoff 0.9 to get phenotype score.   So far we have got the gene score and phenotype score. The final diseasegps score can be obtained by adding the two in a certain ratio, which is preset to 0.5. Users can reset the value.

We also developed a web server of DiseaseGPS, which can be accessed at https://diseasegps.sjtu.edu.cn/. Users can directly input their vcf file containing mutation information. The phenotypic part can be filtered in the HPO tree structure. Users can choose not to register to use, or register to manage their own input and results every time. The web version of DiseaseGPS runs very fast and can bring a good experience to users.

## Web server

diseaseGPS: https://diseasegps.sjtu.edu.cn/

## LICENSE

diseaseGPS is free for non-commercial use without warranty. Users need to obtain licenses such as HPO,ANNOVAR and INTERVAR by themselves. Please contact the authors for commercial use.

## REFERENCE

Jianping Jiang. VCF-Server: A web-based visualization tool for high-throughput variant data mining and management. Mol Genet Genomic Med. 2019;7(7):e00641. https://pubmed.ncbi.nlm.nih.gov/31127704/

Xiaofeng Gong. A new method to measure the semantic similarity from query phenotypic abnormalities to diseases based on the human phenotype ontology. BMC Bioinformatics. 2018;19(Suppl 4):162. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2064-y