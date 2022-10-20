#!/usr/bin/env python
#########################################################################
# Author: Daoyi Huang (huangdaoyi@sjtu.edu.cn)
# Created Time: 2021-3-27 Tuesday
# File Name: diseasegps.py File type: python
# Last Change:.
# Description: python script for  diagnosing genetic diseases
#########################################################################
import os, re, time, sys, platform, optparse, glob

prog = "diseasegps"

version = """%prog 1.0.1 20210913
Written by Daoyi Huang, huangdaoyi@sjtu.edu.cn. 
diseasegps is free for non-commercial use without warranty.
Please contact the authors for commercial use.
Copyright (C) 2021 Hui Lv Lab
============================================================================
"""

usage = """Usage: %prog [OPTION] -i  INPUT -o  OUTPUT ...
       %prog  --config=config_gps.ini ...
"""

description = """=============================================================================
diseasegps
Gene Phenotype Search                                                                       
Diagnose genetic diseases based on genetic data and phenotypic data using python scripts of diseasegps.
=============================================================================
"""
end = """=============================================================================
Thanks for using diseasegps!
Report bugs to huangdaoyi@sjtu.edu.cn;
diseasegps homepage: <http://diseasegps.sjtu.edu.cn/>
=============================================================================
"""

if platform.python_version() < '3.0.0':
    import ConfigParser
else:
    import configparser

paras = {}
triplet_dict = {}


def ConfigSectionMap(config, section):
    global paras
    options = config.options(section)
    for option in options:
        try:
            paras[option] = config.get(section, option)
            if paras[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            paras[option] = None
    return


user_evidence_dict = {}


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False


def flip_freq(feature1=0, feature2=0, feature3=0, feature4=0):
    FREQ_THRESHOLD = float(paras['freq_threshold'])
    value1 = 0
    value2 = 0
    value3 = 0
    value4 = 0
    if is_number(feature1):
        value1 = float(feature1)
    if is_number(feature2):
        value2 = float(feature2)
    if is_number(feature3):
        value3 = float(feature3)
    if is_number(feature4):
        value4 = float(feature4)

    if (value1 > FREQ_THRESHOLD) or (value2 > FREQ_THRESHOLD) or (
            value3 > FREQ_THRESHOLD) or (value4 > FREQ_THRESHOLD):
        return False
    else:
        return True


def flip_genotype(genotype, genes, func, exonicfunc, compound_dict):
    flag = False
    for gene in genes.split(';'):
        if triplet_dict.__contains__((gene, genotype)):
            flag = True
        else:
            if genotype == "het" and compound_dict[gene] >= 2:
                genotype2 = "compound_het"
                if triplet_dict.__contains__(
                    (gene, genotype2)
                ) and func != 'intronic' and exonicfunc != 'synonymous SNV':
                    flag = True
    return flag


########### 检查数据库中的所有数据是否完整 ###########
def check_downdb():
    path = paras['database_locat']
    path = path.strip()
    path = path.rstrip("\/")
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        print("Notice: the folder of %s is created!" % path)
    else:
        print("Warning: the folder of %s is already created!" % path)
    ds = paras['database_names']
    ds.expandtabs(1)
    # database_names = refGene esp6500siv2_all 1000g2015aug_all exac03 gnomad211_genome dbnsfp42a dbscsnv11 avsnp150 clinvar_20210501 dbnsfp31a_interpro rmsk ensGene knownGene
    if not os.path.isfile(paras['annotate_variation']):
        print(
            "Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
            % paras['annotate_variation'])
        if paras['skip_annovar'] != True:
            sys.exit()

    for dbs in ds.split():
        # os.path.isfile(options.table_annovar)
        file_name = dbs
        #if dbs=="1000g2014oct":
        #    file_name="ALL.sites.2014_10"
        if dbs == "1000g2015aug_all":
            file_name = "ALL.sites.2015_08"  # hg19_ALL.sites.2015_08.txt

        dataset_file = paras['database_locat'] + "/" + paras[
            'buildver'] + "_" + file_name + ".txt"
        if dbs != 'rmsk':
            cmd = "perl " + paras['annotate_variation'] + " -buildver " + paras[
                'buildver'] + " -downdb -webfrom annovar " + file_name + " " + paras[
                    'database_locat']
        if dbs == 'rmsk':
            cmd = "perl " + paras['annotate_variation'] + " -buildver " + paras[
                'buildver'] + " -downdb " + file_name + " " + paras[
                    'database_locat']
        if not os.path.isfile(dataset_file):
            if dbs == "1000g2015aug_all":
                file_name = "1000g2015aug"
                dataset_file = paras['database_locat'] + "/" + paras[
                    'buildver'] + "_" + file_name + ".txt"
                cmd = "perl " + paras['annotate_variation'] + " -buildver " + paras[
                    'buildver'] + " -downdb -webfrom annovar " + file_name + " " + paras[
                        'database_locat']
            if paras['skip_annovar'] != True:
                print(
                    "Warning: The Annovar dataset file of %s is not in %s,begin to download this %s ..."
                    % (dbs, paras['database_locat'], dataset_file))

            if paras['skip_annovar'] != True:
                print("%s" % cmd)
                os.system(cmd)


########### 检查输入的基因vcf文件，并生成avinput文件，用于后续的annovar分析 ###########
def check_input():
    inputft = paras['inputfile_type']
    if inputft.lower() == 'avinput':
        return
    if inputft.lower() == 'vcf':
        if os.path.isfile(paras['convert2annovar']):
            #convert2annovar.pl -format vcf4 variantfile > variant.avinput
            cmd = "perl " + paras[
                'convert2annovar'] + " -format vcf4 -withzyg -includeinfo  " + paras[
                    'inputfile'] + "> " + paras['inputfile'] + ".avinput"
            print(
                "Warning: Begin to convert your vcf file of %s to AVinput of Annovar ..."
                % paras['inputfile'])
            print("%s" % cmd)
            os.system(cmd)
        else:
            print(
                "Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                % paras['convert2annovar'])
            if paras['skip_annovar'] != True:
                sys.exit()

    return


########### 调用annovar，生成annovar分析文件 ###########
def check_annovar_result():
    # table_annovar.pl example/ex1.avinput humandb/ -buildver hg19 -out myanno -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp147,ljb26_all,CLINSIG,gnomad_genome   -operation  g,f,f,f,f,f,f   -nastring . -csvout
    inputft = paras['inputfile_type']
    annovar_options = " "
    if re.findall('true', paras['otherinfo'], flags=re.IGNORECASE):
        annovar_options = annovar_options + "--otherinfo "
    if re.findall('true', paras['onetranscript'], flags=re.IGNORECASE):
        annovar_options = annovar_options + "--onetranscript "

    if not os.path.isfile(paras['table_annovar']):
        print(
            "Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
            % paras['table_annovar'])
        if paras['skip_annovar'] != True:
            sys.exit()
    if inputft.lower() == 'avinput':
        cmd = "perl " + paras['table_annovar'] + " " + paras[
            'inputfile'] + " " + paras['database_locat'] + " -buildver " + paras[
                'buildver'] + " -remove -out " + paras[
                    'outfile'] + " -protocol refGene,esp6500siv2_all,1000g2015aug_all,exac03,gnomad211_genome,dbnsfp42a,dbscsnv11,avsnp150,clinvar_20210501,dbnsfp31a_interpro,rmsk,ensGene,knownGene  -operation  g,f,f,f,f,f,f,f,f,f,r,g,g   -nastring ." + annovar_options
        print("%s" % cmd)
        os.system(cmd)
    if inputft.lower() == 'vcf':
        cmd = "perl " + paras['table_annovar'] + " " + paras[
            'inputfile'] + ".avinput " + paras[
                'database_locat'] + " -buildver " + paras[
                    'buildver'] + " -remove -out " + paras[
                        'outfile'] + " -protocol refGene,esp6500siv2_all,1000g2015aug_all,exac03,gnomad211_genome,dbnsfp42a,dbscsnv11,avsnp150,clinvar_20210501,dbnsfp31a_interpro,rmsk,ensGene,knownGene  -operation  g,f,f,f,f,f,f,f,f,f,r,g,g    -nastring ." + annovar_options
        print("%s" % cmd)
        os.system(cmd)

    return


########### 根据种群频率筛选annovar注释后的突变 ###########
def prescreen_frequence(anvfile):
    os.rename(anvfile, anvfile + ".prescreen")
    newoutfile = anvfile
    anvfile = anvfile + ".prescreen"
    dfreqs = paras['freq_prscreen']
    dfreqs.expandtabs(1)
    Freqs_flgs = {}
    for dbfreq in dfreqs.split():
        Freqs_flgs[dbfreq] = 0
    Gene_flgs = {
        'Otherinfo': 0,
        'Gene.refGene': 0,
        'ExonicFunc.refGene': 0,
        'Func.refGene': 0
    }

    # Freqs_flgs = {
    #     'esp6500siv2_all': 0,
    #     '1000g2015aug_all': 0,
    #     'ExAC_ALL': 0,
    #     'AF_popmax': 0
    # }
    try:
        fh = open(anvfile, "r")
        fw = open(newoutfile, "w")
        strs = fh.read()
        # 新增的compound统计部分
        line_sum = 0
        compound_dict = {}
        for line in strs.split('\n'):
            cls = line.split('\t')
            if len(cls) < 2: break
            if line_sum == 0:
                search_key_index(line, Gene_flgs)
            else:
                for gene in cls[Gene_flgs['Gene.refGene']].split(';'):
                    if cls[Gene_flgs['Func.refGene']] != 'intronic' and cls[
                            Gene_flgs[
                                'ExonicFunc.refGene']] != 'synonymous SNV':
                        if not compound_dict.__contains__(gene):
                            compound_dict[gene] = 1
                        else:
                            compound_dict[gene] += 1
                    else:
                        if not compound_dict.__contains__(gene):
                            compound_dict[gene] = 0
            line_sum += 1

        line_sum = 0
        for line in strs.split('\n'):
            cls = line.split('\t')
            if len(cls) < 2: break
            if line_sum == 0:
                search_key_index(line, Freqs_flgs)
                line = line.replace('AF\t', 'gnomAD_genome_ALL\t')
                line = line.replace('AF_afr', 'gnomAD_genome_AFR')
                line = line.replace('AF_amr', 'gnomAD_genome_AMR')
                line = line.replace('AF_eas', 'gnomAD_genome_EAS')
                line = line.replace('AF_fin', 'gnomAD_genome_FIN')
                line = line.replace('AF_nfe', 'gnomAD_genome_NFE')
                line = line.replace('AF_oth', 'gnomAD_genome_OTH')
                line = line.replace('AF_asj', 'gnomAD_genome_ASJ')
                line = line.replace('CLNDN', 'CLNDBN')
                line = line.replace('CLNACC', 'CLNALLELEID')
                line = line.replace('CLNDSDB', 'CLNDISDB')
                line = line.replace('phyloP100way_vertebrate',
                                    'phyloP46way_placental')
                line = line.replace('avsnp150', 'avsnp147')
                fw.write(line + '\n')
                line_sum += 1

            else:
                freq_dic = {
                    'esp6500siv2_all': 0,
                    '1000g2015aug_all': 0,
                    'ExAC_ALL': 0,
                    'AF_popmax': 0
                }
                for k in Freqs_flgs.keys():
                    freq_dic[k] = cls[Freqs_flgs[k]]

                if flip_freq(freq_dic['esp6500siv2_all'],
                             freq_dic['1000g2015aug_all'],
                             freq_dic['ExAC_ALL'],
                             freq_dic['AF_popmax']) and flip_genotype(
                                 cls[Gene_flgs['Otherinfo']],
                                 cls[Gene_flgs['Gene.refGene']],
                                 cls[Gene_flgs['Func.refGene']],
                                 cls[Gene_flgs['ExonicFunc.refGene']],
                                 compound_dict):
                    fw.write(line + '\n')

                # if flip_freq(cls[Freqs_flgs['esp6500siv2_all']],
                #              cls[Freqs_flgs['1000g2015aug_all']],
                #              cls[Freqs_flgs['ExAC_ALL']],
                #              cls[Freqs_flgs['AF_popmax']]):
                #     fw.write(line + '\n')
                line_sum += 1

    except IOError:
        print("Error: can\'t read/write the annovar output file %s %s" %
              (anvfile, newoutfile))
        sys.exit()
        return
    except IndexError:
        print(line)
        print(cls)
        sys.exit()
    else:
        pass
        fh.close()
        fw.close()

    return (sum)


def prescreen_annovar_result():
    for annovar_outfile in glob.iglob(
            paras['outfile'] + "*." + paras['buildver'] +
            "_multianno.txt"):  # 在目录下搜索指定匹配模式文件，*代表任意0到多字符
        prescreen_frequence(annovar_outfile)


def read_datasets():
    # triplet_OmimGeneGenotype
    try:
        fh = open(paras['triplet'], "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            cls2 = line2.split('\t')
            if len(cls2) > 1:
                if not triplet_dict.__contains__((cls2[0], cls2[3])):
                    triplet_dict[(cls2[0], cls2[3])] = set()
                    triplet_dict[(cls2[0], cls2[3])].add(cls2[2])
                else:
                    triplet_dict[(cls2[0], cls2[3])].add(cls2[2])

                if not triplet_dict.__contains__((cls2[0], '.')):
                    triplet_dict[(cls2[0], '.')] = set()
                    triplet_dict[(cls2[0], '.')].add(cls2[2])
                else:
                    triplet_dict[(cls2[0], '.')].add(cls2[2])

    except IOError:
        print("Error: can\'t read the triplet_OmimGeneGenotype file %s" %
              paras['triplet'])
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()


########### 读取输入的phenotype文件中的表型HPO ###########
def check_hpoterms():
    if not os.path.isfile(paras['phenotype']):
        print(
            "Warning: The hpoterms file [ %s ] is not here, please input the correct hpoterms file"
            % paras['phenotype'])
        sys.exit()
    sample = set()
    with open(paras['phenotype'], 'r') as f:
        hpoterms = f.readlines()
        for hp_tmp in hpoterms:
            hp_tmp = hp_tmp.strip()
            m1 = re.match('HP:\d{7}', hp_tmp)
            if m1 is not None:
                sample.add(hp_tmp)

    if len(sample) == 0:
        print(
            "Warning: The hpoterms file [ %s ] don't meet specifications, please enter one HPOterm per line"
            % paras['phenotype'])
        sys.exit()
    else:
        return sample


########### 调用intervar生成intervar结果文件 ###########
def check_intervar_result():
    inputft = paras['inputfile_type']
    if not os.path.isfile(paras['intervarpy']):
        print(
            "Warning: The InterVar file [ %s ] is not here,please download InterVar firstly: https://github.com/WGLab/InterVar"
            % paras['intervarpy'])

    if inputft.lower() == 'vcf':
        cmd = 'cd ' + paras['database_intervar'] + " && python3 " + paras[
            'intervar_file'] + " -b " + paras[
                'buildver'] + " --skip_annovar " + " -o " + paras['outfile']
        print("%s" % cmd)
        os.system(cmd)

    if inputft.lower() == 'avinput':
        cmd = 'cd ' + paras['database_intervar'] + " && python3 " + paras[
            'intervar_file'] + " -b " + paras[
                'buildver'] + " --skip_annovar " + " -o " + paras['outfile']
        print("%s" % cmd)
        os.system(cmd)

    print('--------intervar over------------')

    return


def search_key_index(line, dict):
    cls = line.split('\t')
    for key in dict.keys():
        for i in range(1, len(cls)):
            ii = i - 1
            if key == cls[ii]:
                dict[key] = ii
                break
    return


def my_disease_gps(hpoterms):

    hpo_list = list(hpoterms)
    newoutfile = paras['outfile'] + ".ACMG"
    newoutfile2 = paras['outfile'] + "_gene.list"
    newoutfile3 = paras['outfile'] + "_hpo.list"
    newoutfile4 = paras['outfile'] + ".gps"

    # call function
    try:
        for intervar_outfile in glob.iglob(
                paras['outfile'] + "*." + paras['buildver'] +
                "_multianno.txt.intervar"):  # 在目录下搜索指定匹配模式文件，*代表任意0到多字符
            os.system('perl {} {}  {}  {}'.format(paras['intervar2vcf'],
                                                  intervar_outfile, newoutfile,
                                                  newoutfile2))

            os.system('{} -omim {} -obo {} -term {} -o {}'.format(
                paras['hpoc'], paras['all_source_omim'], paras['obo'],
                ','.join(hpo_list), newoutfile3))
    except IOError:
        print(
            "Error: can\'t readi/write the annovar output files %s\t%s\t%s\t%s"
            % (newoutfile, newoutfile2, newoutfile3, newoutfile4))
        sys.exit()
        return

    # combine
    try:
        result_gene = {}
        with open(newoutfile2, 'r') as fh:
            fh_readlines = fh.readlines()
            for line in fh_readlines:
                line = line.strip()
                line_split = line.split('\t')
                if re.match('ORPHA', line_split[0]):
                    continue
                if not result_gene.__contains__(line_split[0]):
                    result_gene[line_split[0]] = float(line_split[1])

        # rank_result_gene = sorted(result_gene.items(),
        #                         key=lambda result_gene: result_gene[1],
        #                         reverse=True)

        # f_result_gene = open(
        #     paras['outfile'] + "_rank_gene.list", 'w')
        # current_count = 0
        # while current_count < len(rank_result_gene):
        #     min_count = current_count
        #     current_list = []
        #     current_list.append(rank_result_gene[current_count])

        #     while current_count + 1 < len(rank_result_gene):
        #         if rank_result_gene[current_count][1] == rank_result_gene[
        #                 current_count + 1][1]:
        #             current_list.append(rank_result_gene[current_count + 1])
        #             current_count += 1
        #         else:
        #             break

        #     max_count = current_count
        #     final_count = int((min_count + max_count) / 2) + 1
        #     current_count += 1

        #     for j in current_list:
        #         f_result_gene.write('{}\t{}\t{}\n'.format(final_count, j[0], j[1]))

        result_hpo = {}
        with open(newoutfile3, 'r') as fh:
            fh_readlines = fh.readlines()
            for line in fh_readlines:
                line = line.strip()
                line_split = line.split('\t')
                if re.match('ORPHA', line_split[0]):
                    continue
                if not result_hpo.__contains__(line_split[0]):
                    result_hpo[line_split[0]] = float(line_split[1])
                else:
                    result_hpo[line_split[0]] = max(float(line_split[1]),
                                                    result_hpo[line_split[0]])

        # rank_result_hpo = sorted(result_hpo.items(),
        #                         key=lambda result_hpo: result_hpo[1],
        #                         reverse=True)

        # f_result_hpo = open(
        #     paras['outfile'] + "_rank_hpo.list", 'w')
        # current_count = 0
        # while current_count < len(rank_result_hpo):
        #     min_count = current_count
        #     current_list = []
        #     current_list.append(rank_result_hpo[current_count])

        #     while current_count + 1 < len(rank_result_hpo):
        #         if rank_result_hpo[current_count][1] == rank_result_hpo[
        #                 current_count + 1][1]:
        #             current_list.append(rank_result_hpo[current_count + 1])
        #             current_count += 1
        #         else:
        #             break

        #     max_count = current_count
        #     final_count = int((min_count + max_count) / 2) + 1
        #     current_count += 1

        #     for j in current_list:
        #         f_result_hpo.write('{}\t{}\t{}\n'.format(final_count, j[0], j[1]))

        f_result_combine = open(newoutfile4, 'w')
        result_combine = {}
        for i in set(result_gene.keys()) | set(result_hpo.keys()):
            result_combine[i] = result_gene.get(
                i, 0.0) * float(paras['gene_hpo_cutoff']) + result_hpo.get(
                    i, 0.0) * (1 - float(paras['gene_hpo_cutoff']))

        rank_result_combine = sorted(
            result_combine.items(),
            key=lambda result_combine: result_combine[1],
            reverse=True)

        current_count = 0
        while current_count < len(rank_result_combine):
            min_count = current_count
            current_list = []
            current_list.append(rank_result_combine[current_count])

            while current_count + 1 < len(rank_result_combine):
                if rank_result_combine[current_count][
                        1] == rank_result_combine[current_count + 1][1]:
                    current_list.append(rank_result_combine[current_count + 1])
                    current_count += 1
                else:
                    break

            max_count = current_count
            final_count = int((min_count + max_count) / 2) + 1
            current_count += 1

            for j in current_list:
                f_result_combine.write('{}\t{}\t{}\n'.format(
                    final_count, j[0], j[1]))

    except IOError:
        print(
            "Error: can\'t readi/write the annovar output files %s\t%s\t%s\t%s"
            % (newoutfile, newoutfile2, newoutfile3, newoutfile4))
        sys.exit()
        return
    else:
        # f_result_gene.close()
        # f_result_hpo.close()
        f_result_combine.close()


########### 主程序 ###########
def main():

    if platform.python_version() < '3.0.0':
        config = ConfigParser.ConfigParser()
    else:
        config = configparser.ConfigParser()

    parser = optparse.OptionParser(usage=usage,
                                   version=version,
                                   description=description)

    parser.add_option("-?",
                      action="help",
                      help=optparse.SUPPRESS_HELP,
                      dest="help")
    parser.add_option("-v",
                      action="version",
                      help=optparse.SUPPRESS_HELP,
                      dest="version")

    parser.add_option(
        "-c",
        "--config",
        dest="config",
        action="store",
        help=
        "The config file of all options. it is for your own configure file.You can edit all the options in the configure and if you use this options,you can ignore all the other options bellow",
        metavar="config_gps.ini")

    parser.add_option(
        "-b",
        "--buildver",
        dest="buildver",
        action="store",
        help=
        "The genomic build version, it can be hg19 or hg38 , will support other version later",
        metavar="hg19")

    parser.add_option("-i",
                      "--input",
                      dest="input",
                      action="store",
                      help="The input file contains your variants",
                      metavar="test/testNjob1.vcf")

    parser.add_option(
        "-p",
        "--phenotype",
        dest="phenotype",
        action="store",
        help="The input file contains your hpoterms separated by commas",
        metavar="test/testNjob1_hpo.txt")

    parser.add_option(
        "--input_type",
        dest="input_type",
        action="store",
        help=
        "The input file type, it can be  AVinput(Annovar's format),VCF(VCF with single sample)",
        metavar="VCF")

    parser.add_option(
        "-o",
        "--output",
        dest="output",
        action="store",
        help=
        "The prefix of output file which contains the results, the file of results will be as [$$prefix].gps ",
        metavar="test/testNjob1")

    group = optparse.OptionGroup(
        parser, "Database location Options",
        "Caution: check these options from manual of INTERVAR. The INTERVAR version should be >=  2019-03-27, older verions of INTERVAR will bring problems."
    )
    group.add_option(
        "-t",
        "--database_gps",
        dest="database_gps",
        action="store",
        help="The  database location/dir for the gps dataset files",
        metavar="gpsdb")
    group.add_option(
        "-d",
        "--database_locat",
        dest="database_locat",
        action="store",
        help="The  database location/dir for the annotation datasets",
        metavar="humandb")
    group.add_option(
        "-r",
        "--database_intervar",
        dest="database_intervar",
        action="store",
        help="The  database location/dir for the intervar dataset files",
        metavar="intervardb")
    group.add_option(
        "-a",
        "--database_annovar",
        dest="database_annovar",
        action="store",
        help="The  database location/dir for the annovar dataset files",
        metavar="annovardb")
    parser.add_option_group(group)

    group = optparse.OptionGroup(
        parser, "Annovar Options",
        "Caution: check these options from manual of Annovar. The ANNOVAR version should be >=  2018-04-16, older verions of ANNOVAR will bring problems."
    )
    group.add_option("--table_annovar",
                     action="store",
                     help="The Annovar perl script of table_annovar.pl",
                     metavar="./table_annovar.pl",
                     dest="table_annovar")
    group.add_option("--convert2annovar",
                     action="store",
                     help="The Annovar perl script of convert2annovar.pl",
                     metavar="./convert2annovar.pl",
                     dest="convert2annovar")
    group.add_option("--annotate_variation",
                     action="store",
                     help="The Annovar perl script of annotate_variation.pl",
                     metavar="./annotate_variation.pl",
                     dest="annotate_variation")

    group.add_option(
        "--skip_annovar",
        action="store_true",
        help=
        "Skip the Annovar annotation, this can be true only after you  already got the annovar's annotation results",
        dest="skip_annovar")

    parser.add_option_group(group)
    group = optparse.OptionGroup(
        parser, "Examples",
        """./diseasegps.py -c config_gps.ini  # Run the examples in config_gps.ini
        ./diseasegps.py  -b hg19 -i your_vcf_input  --input_type=VCF -p your_phenotype_input  -o your_output """
    )
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    #(options,args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    print("%s" % description)
    print("%s" % version)
    print("Notice: Your command of diseasegps is %s" % sys.argv[:])

    config_file = os.path.join(os.path.dirname(__file__), "config_gps.ini")
    if os.path.isfile(config_file):
        config.read(config_file)
        sections = config.sections()
        for section in sections:
            ConfigSectionMap(config, section)
    else:
        print(
            "Error: The default configure file of [ config_gps.ini ] is not here, exit! Please redownload the diseasegps."
        )
        sys.exit()


#begin to process user's options:
    if options.config is not None:
        if os.path.isfile(options.config):
            config.read(options.config)
            sections = config.sections()
            for section in sections:
                ConfigSectionMap(config, section)
        else:
            print(
                "Error: The config file [ %s ] is not here,please check the path of your config file."
                % options.config)
            sys.exit()

    if options.buildver != None:
        paras['buildver'] = options.buildver
    if options.database_locat != None:
        paras['database_locat'] = options.database_locat
    if options.input != None:
        paras['inputfile'] = options.input
    if options.input_type != None:
        paras['inputfile_type'] = options.input_type
    if options.phenotype != None:
        paras['phenotype'] = options.phenotype
    if options.output != None:
        paras['outfile'] = options.output
    if options.database_gps != None:
        paras['database_gps'] = options.database_gps
    if options.database_intervar != None:
        paras['database_intervar'] = options.database_intervar
    if options.database_annovar != None:
        paras['database_annovar'] = options.database_annovar

    paras['skip_annovar'] = False

    if options.skip_annovar == True:
        paras['skip_annovar'] = True

    if options.table_annovar != None and options.skip_annovar != True:
        if os.path.isfile(options.table_annovar):
            paras['table_annovar'] = options.table_annovar
        else:
            print(
                "Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                % options.table_annovar)
            if options.skip_annovar != True:
                sys.exit()
    if options.convert2annovar != None and options.skip_annovar != True:
        if os.path.isfile(options.convert2annovar):
            paras['convert2annovar'] = options.convert2annovar
        else:
            print(
                "Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                % options.convert2annovar)
            if options.skip_annovar != True:
                sys.exit()
    if options.annotate_variation != None and options.skip_annovar != True:
        if os.path.isfile(options.annotate_variation):
            paras['annotate_variation'] = options.annotate_variation
        else:
            print(
                "Warning: The Annovar file [ %s ] is not here,please download ANNOVAR firstly: http://www.openbioinformatics.org/annovar"
                % options.annotate_variation)
            if options.skip_annovar != True:
                sys.exit()

    if not os.path.isfile(paras['inputfile']):
        print(
            "Error: Your input file [ %s ] is not here,please check the path of your input file."
            % paras['inputfile'])
        sys.exit()
    if not os.path.isfile(paras['phenotype']):
        print(
            "Error: Your phenotype file [ %s ] is not here,please check the path of your input file."
            % paras['phenotype'])
        sys.exit()

    print("INFO: The options are %s " % paras)
    check_downdb()
    check_input()
    hpoterms = check_hpoterms()

    if options.skip_annovar != True:
        check_annovar_result()  #  to obtain myanno.hg19_multianno.csv
    else:
        print(
            "Warning: You activated the option of --skip_annovar, the Annovar will not run!"
        )
        print(
            "Warning: The diseasegps will interpret the variants based on your old annotation information!"
        )

    prescreen_annovar_result()
    check_intervar_result()  # to obtain .txt.intervar
    my_disease_gps(hpoterms)

    print('---------------annovar, prescreen and intervar over--------------')

if __name__ == "__main__":
    main()