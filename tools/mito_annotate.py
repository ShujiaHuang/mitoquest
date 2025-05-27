"""
This script is modified according to the following script:
https://github.com/leklab/mitochondrial_constraint/blob/main/build_model/annotate_mutations.py
"""

import argparse
import sys
import gzip
import csv
import datetime
import json
from itertools import zip_longest


def rcrs_pos_to_ref(anno_file_path):
    """Generate dictionary linking each position to its reference nucleotide in the rCRS.

    :return: dictionary where the key is the position in rCRS, and the value is its reference nucleotide
    """
    dictionary = {}
    with gzip.open(anno_file_path+'/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf.gz', 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            dictionary[row["POS"]] = row["REF"]
        
    return dictionary


def rcrs_pos_to_trinucleotide(anno_file_path):
    """Generate dictionary linking each position to its reference trinucleotide in the rCRS.

    :return: dictionary where the key is the position in rCRS, and the value is its reference trinucleotide
    """
    # first, generate dictionary to convert coordinate to reference nucleotide
    rcrs_pos2ref = rcrs_pos_to_ref(anno_file_path)
    # now, generate dictionary of coordinate to trinucleotide
    dict = {}
    with gzip.open(anno_file_path+'/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf.gz', 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pos = int(row["POS"])
            ref = row["REF"]

            if pos == 16569:
                # dealing with circular genome
                trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(1)]
            elif pos == 1:
                # dealing with circular genome
                trinucleotide = rcrs_pos2ref[str(16569)] + ref + rcrs_pos2ref[str(pos + 1)]
            elif ref == "N":
                continue  # ie skip, to handle the 'N' spacer expected at m.3107
            elif rcrs_pos2ref[str(pos + 1)] == "N":
                trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(pos + 2)]
            elif rcrs_pos2ref[str(pos - 1)] == "N":
                trinucleotide = rcrs_pos2ref[str(pos - 2)] + ref + rcrs_pos2ref[str(pos + 1)]
            else:
                trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(pos + 1)]
            dict[str(pos)] = trinucleotide
            
    return dict


def dbSNP_annotate(anno_file_path):
    """Generate dictionary with the rsid for each variant in dbSNP.
    """
    dict = {}
    with gzip.open(anno_file_path+"/databases/dbSNP.chrM.vcf.gz", "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            chrom, pos, rsid, refs, alts, *_skip = line.strip().split('\t')
            
            # 这一步是必须的，因为有些变异的REF和ALT由于多突变的原因后缀是相同的
            REF_ALT_list = [remove_common_suffix(REF, ALT) for REF, ALT in list(
                zip_longest(refs.split(','), alts.split(','), fillvalue=refs.split(',')[0]))]
            
            # work for multiallelic variants 
            for r, a in REF_ALT_list:
                # handle the case where there are multiple alts, e.g. m.1234A>G,T
                if a == '.': continue
                dict[(chrom, pos, r, a)] = rsid
            
    return dict


def gnomad_annotate(anno_file_path):
    """Generate dictionary with the maximum observed heteroplasmy (max_hl) and other annotations for each variant in gnomAD.
    soure:https://gnomad.broadinstitute.org/downloads
    check updated 20250403
    :return: a dictionary, where the key is a tuple of the position, ref and alt, and the values are annotations
    """
    dict = {}
    with gzip.open(anno_file_path+'/databases/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv.gz', 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row["filters"] == "PASS":
                dict[(row["ref"], row["position"], row["alt"])] = (row["max_observed_heteroplasmy"], 
                                                                   row["AF_hom"], 
                                                                   row["AF_het"], 
                                                                   row["AC_hom"],
                                                                   row["AC_het"])
    return dict


def phylop_annotate(anno_file_path):
    """Generate dictionary with the phyloP scores for conservation, from 100 vertebrates.

    :return: a dictionary, where the key is the position, and the value is the phyloP conservation score
    """
    dict = {}
    pos = 0  # to handle a header row
    for row in gzip.open(anno_file_path+'/insilicos/chrM.phyloP100way.wigFix.gz', 'rt'):
        dict[pos] = row.replace('\n', '')
        pos += 1
    return dict


def vep_annotate(anno_file_path, anno_column_info=None):
    """Create a dictionary of the VEP annotations for every possible single nucleotide variant in the mtDNA.

    :return: dictionary where tuple of the variant and value identifier is key, and value is list of annotations
    """
    if anno_column_info is None:
        return {}
    
    vep = {}
    # use vcf where variants in two genes are split across two rows, for easy parsing
    with gzip.open(anno_file_path+"/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            csq = {}  # to store the VEP annotations for each variant in dictionary
            for label in anno_column_info:
                if label in row:
                    if label == "Allele":
                        # The `Allele` column is the first column in the VEP output, and the format is: 'CSQ=Allele' in the input file
                        csq[label] = row[label].split("=")[-1]
                    else:
                        csq[label] = row[label]
                else:
                    csq[label] = ""
            
            vep[(row["REF"], row["POS"], row["ALT"])] = csq  # if variant in two genes, will only keep the last one
            
    return vep


def tRNA_positions(anno_file_path):
    """Create a dictionary of the tRNA position numbers (ranging from 1-73) encoded by mtDNA.

    :return: dictionary where position in the mtDNA is key and the value is the tRNA position
    """
    dict = {}
    for row in csv.DictReader(open(anno_file_path+'/other_annotations/tRNA_positions.txt'), delimiter='\t'):
        if row["m.pos"] in dict:  # if the mtDNA position is in two tRNA genes
            dict[row["m.pos"]].append(row["tRNA_position"])
        else:
            dict[row["m.pos"]] = [row["tRNA_position"]]
    return dict


def RNA_domains_mods(anno_file_path):
    """Create a dictionary of modified RNA bases and tRNA domains encoded by mtDNA.

    :return: dictionary where position in the mtDNA is key and values are whether the base is modified and tRNA domain
    """
    dict = {}
    with gzip.open(anno_file_path+"/other_annotations/mito_RNA_modifications_and_domains.txt.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # handle the fact few bases are in two RNA genes
            if row["MODIFICATION"] != "#N/A":
                if "," in row["GENE"]:  # bases in two RNA genes coded together
                    dict[row["POS"], "modified"] = [row["MODIFICATION"], row["MODIFICATION"]]
                else:
                    dict[row["POS"], "modified"] = [row["MODIFICATION"]]

            if row["DOMAIN"]:
                if "," in row["GENE"]:  # bases in two RNA genes coded together
                    dict[row["POS"], "domain"] = [row["DOMAIN"], row["DOMAIN"]]
                else:
                    dict[row["POS"], "domain"] = [row["DOMAIN"]]
                    
    # Note, there are 4 discrepancies from these domain annotations and the pair type
    # at pair m.8307+8311 and m.10039+10046, in mamit-tRNAdb they are drawn as pairs, but in the pdf are in loop domain
    return dict


def RNA_base_type(anno_file_path):
    """Create a dictionary of RNA base types encoded by mtDNA (Watson-Crick (WC), non-WC, or loop/other).

    :return: dictionary where position in the mtDNA is key and the value is the base type
    """
    # first, create a dictionary for bases in pairs so can look-up if WC or not
    RNA_dict = {}
    with gzip.open(anno_file_path+"/other_annotations/all_RNA_bases.tsv.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # type b excludes base m.3107N, filter to bases in pairs
            if (row["Type"] == "b") and row["Pair_coordinate"]:
                # note there are four bases in two tRNAs, hence using gene as second key
                RNA_dict[(row["Genomic_coordinate"], row["file"])] = row["RNA.base"]
            
    # now iterate through RNA bases
    dict = {}
    with gzip.open(anno_file_path+"/other_annotations/all_RNA_bases.tsv.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["Type"] == "b":  # only excludes base m.3107N
                # first determine base type
                if row["Pair_coordinate"]:  # if in pair
                    base1 = RNA_dict[(row["Genomic_coordinate"], row["file"])]
                    # pairing base
                    base2 = RNA_dict[(row["Pair_coordinate"], row["file"])]
                    if (base1 == 'A' and base2 == 'T') or (base1 == 'T' and base2 == 'A'):
                        base_type = "WC"
                    elif (base1 == 'C' and base2 == 'G') or (base1 == 'G' and base2 == 'C'):
                        base_type = "WC"
                    else:
                        base_type = "non-WC"
                else:
                    base_type = "loop-or-other"
                    
                # now build dictionary, appending to handle bases in two genes
                if row["Genomic_coordinate"] not in dict:
                    dict[row["Genomic_coordinate"]] = [base_type]
                else:
                    dict[row["Genomic_coordinate"]].append(base_type)
                    
    return dict


def uniprot_annotations(anno_file_path):
    """Create a dictionary of the UniProt annotations for mtDNA-encoded proteins.

    :return: dictionary where position in the mtDNA is key and the value is the UniProt annotation
    """
    # on the UniProt website, they display annotations for 'sites' and 'topology' for each protein, if available
    # but no family or domains for mtDNA-encoded proteins in UniProt
    dict = {}
    with gzip.open(anno_file_path+"/other_annotations/uniprot_beds_chrMT_combined.20250403.txt.gz", "rt") as f:
        for line in f:
            row = line.strip().split('\t')
            # restrict to annotations of interest
            if any(x in row[14] for x in ["binding", "metal"]):  # metal binding or binding site
                annotation = "site:" + row[14].split("UP000005640_9606_")[1].split(".bed")[0] + "-" + row[13].replace(";", "|")
                # per UniProt: start_coord = row[1] and end_coord = row[2], but there can be intervals/blocks between these
                # number of blocks representing the annotation is row[9]
                # row[10] are the block sizes, a comma separated list
                # row[11] are block starts, a comma separated list of block offsets relative to the annotation start
                nblock = list(range(1, int(row[9]) + 1))
                for block in nblock:
                    start = int(row[1]) + int(row[11].split(",")[(block - 1)]) + 1  # need the plus 1 offset
                    # need the minus 1 offset
                    end = start + int(row[10].split(",")[(block - 1)]) - 1
                    if end > int(row[2]):
                        sys.exit('[ERROR]: the predicted end coordinate is greater than the provided')
                    for pos in list(range(1, 16570)):
                        if (pos >= start) and (pos <= end):
                            if pos in dict:
                                if annotation not in dict[pos]:
                                    dict[pos].append(annotation)
                            else:
                                dict[pos] = [annotation]
    return dict


def curated_func_sites(anno_file_path):
    """Residues involved in complex I proton transfer, curated from PMID:32972993.

    :return: dictionary where position in the mtDNA is key and the value is the annotation label
    """
    # file with CI proton residues are labelled by their protein position/residue number
    # so first create a lookup to go from protein and protein position to mtDNA coordinate
    res_to_pos = {}
    with gzip.open(anno_file_path+"/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (row["SYMBOL"], row["Protein_position"]) not in res_to_pos:
                res_to_pos[(row["SYMBOL"], row["Protein_position"])] = [
                    int(row["POS"])]
            else:
                res_to_pos[(row["SYMBOL"], row["Protein_position"])
                        ].append(int(row["POS"]))

    dict = {}
    for row in csv.DictReader(open(anno_file_path+'/other_annotations/CI_proton_residues_PMID32972993.txt'), delimiter='\t'):
        for pos in list(range(1, 16570)):
            if pos in res_to_pos[(row["locus"], row["residue"])]:
                dict[pos] = "proton-transfer-PMID32972993"
                
    return dict


def apogee(anno_file_path):
    """Create a dictionary of the APOGEE scores for missense.

    :return: dictionary where position in the mtDNA is key and the value is the APOGEE score
    """
    dict = {}
    with gzip.open(anno_file_path+"/insilicos/MitImpact_db_3.1.3.txt.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # note, for variants in two genes this won't distangle which gene has which score
            if (row["Ref"], row["Start"], row["Alt"]) in dict:  # i.e. variant is in two genes
                dict[(row["Ref"], row["Start"], row["Alt"])].append(row["APOGEE1"])
            else:
                dict[(row["Ref"], row["Start"], row["Alt"])] = [row["APOGEE1"]]

    return dict


def mitotip(anno_file_path):
    """Create a dictionary of the MitoTip in silico scores for tRNA variants.

    :return: dictionary where position in the mtDNA is key and the value is the MitoTip classification
    """
    # note variants in two genes have same score
    dict = {}
    with gzip.open(anno_file_path+"/insilicos/mitotip_scores.txt.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            prediction = ''
            if row["Quartile"] == "Q1":
                prediction = "confirmed_pathogenic" if row["Mitomap_Status"] == "Confirmed" else "likely_pathogenic"
            elif row["Quartile"] == "Q2":
                prediction = "possibly_pathogenic"
            elif row["Quartile"] == "Q3":
                prediction = "possibly_benign"
            elif row["Quartile"] == "Q4":
                prediction = "likely_benign"
                
            dict[(row["rCRS"], row["Position"], row["Alt"])] = prediction

    return dict


def hmtvar(anno_file_path):
    """Create a dictionary of the HmtVar in silico scores for tRNA variants.

    :return: dictionary where the key is a tuple of the ref, position and alt, and the value is the HmtVar classification
    """
    dict = {}
    with gzip.open(anno_file_path+"/insilicos/hmtvar_annotations.txt.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            insilico = ""
            # extract the in silico prediction from the annotation
            if len(row["HmtVar"]) > 3:
                annotation = json.loads(row["HmtVar"])
                insilico = str(annotation["pathogenicity"])
            dict[(row["REF"], row["POS"], row["ALT"])] = insilico

    return dict


def in_helix(anno_file_path):
    """Generate dictionary with the maximum observed heteroplasmy (max_hl) for each variant in HelixMTdb.
    source:https://www.helix.com/mitochondrial-variant-database
    check updated 20250403
    :return: a dictionary, where the key is a tuple of the ref, position and alt, and the value is the max_hl
    """
    dict = {}
    with gzip.open(anno_file_path+"/databases/HelixMTdb_20200327.tsv.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["alleles"].count(",") == 1:  # skip multiallelic variant calls
                pos = row["locus"].split("chrM:")[1]
                ref = row["alleles"].split("\"")[1]
                alt = row["alleles"].split("\"")[3]

                # calculate maximum heteroplasmy, as only provided for variants observed at heteroplasmy
                if float(row["AF_hom"]) == 0:
                    max_het = float(row["max_ARF"])
                else:  # then homoplasmic are observed
                    max_het = 1
                    
                dict[(ref, pos, alt)] = (max_het, row["AF_hom"], row["AF_het"]) 
                
    return dict


def mitomap_locus(anno_file_path):
    """Generate a list of the genome locus in MITOMAP: 
    https://www.mitomap.org/foswiki/bin/view/MITOMAP/GenomeLoci
    
    Data format:
        "Map Locus","Starting","Ending","Shorthand","Description","Reference"
        "MT-3H","384","391","CR:mt3H","mt3 H-strand control element","1"
        "MT-3L","16499","16506","CR:mt3L","L-strand control element [on complement]","1"

    Args:
        anno_file_path (_type_): path to the annotation files

    Returns:
        _type_: Array
    """
    locus = []
    for row in csv.DictReader(open(anno_file_path+'/databases/mitomap_genome_loci.csv')):
        start, end = int(row["Starting"]), int(row["Ending"])
        if end < start: 
            continue
        
        locus.append([start, end, row["Map Locus"]])
    
    locus.sort(key=lambda x: x[0])
    return locus


def mitomap(anno_file_path):
    """Generate dictionary with the disease association statys for each variant in MITOMAP.

    :return: a dictionary, where the key is a tuple of the ref, position and alt, and the value is the status
    """
    dict1, dict2 = {}, {}
    with gzip.open(anno_file_path+"/databases/MITOMAP_disease_20250403.txt.gz", "rt", encoding="ISO-8859-1") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if ("Cfrm" in str(row["status"])) or ("Reported" in str(row["status"])):
                if (len(row["ref"]) == 1) and (len(row["alt"]) == 1) and row["alt"].isalpha() and (row["ref"] != row["alt"]):  # if SNVs
                    dict1[(row["ref"], row["pos"], row["alt"])] = (row["status"].replace(";", "|"),  # replace ';' with '|' to prevent error in VCF INFO 
                                                                   row["homoplasmy"].replace(";", "|"), 
                                                                   row["heteroplasmy"].replace(";", "|"), 
                                                                   row["disease"].replace(";", "|").strip().replace(
                                                                       "/ ", "/").replace(" /", "/").replace(" ", "_"))

    with gzip.open(anno_file_path+"/databases/MITOMAP_polymorphisms_20250403.txt.gz", "rt", encoding="ISO-8859-1") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if (len(row["ref"]) == 1) and (len(row["alt"]) == 1) and row["alt"].isalpha() and (row["ref"] != row["alt"]):  # if SNVs
                # 61,883 is total number of gb sequences to convert to allele freq
                # 2025-01-28: [GenBank Frequency Information. The current GenBank frequency data in our variant tables is derived from 61,883 human mitochondrial DNA](https://www.mitomap.org/MITOMAP/GBFreqInfo)
                dict2[(row["ref"], row["pos"], row["alt"])] = (int(row["gbcnt"]), (int(row["gbcnt"]) / 61883))

    return dict1, dict2


def clinvar(anno_file_path):
    """Generate dictionary with the clinical significance interpretation for each variant in ClinVar.

    :return: a dictionary, where the key is a tuple of the ref, position and alt, and the value is the interpretation
    """
    dict = {}
    with gzip.open(anno_file_path+"/databases/clinvar_result_chrMT_SNVonly_05252022.txt.gz", "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pos = row["GRCh38Location"]
            alt = row["Canonical SPDI"].split(':')[-1]
            ref = row["Canonical SPDI"].split(':')[2]
            interp = row["Clinical significance (Last reviewed)"].split('(')[0]
            
            if (len(ref) == 1) and (len(alt) == 1) and (ref != alt):  # if SNVs    
                # if "no assertion criteria" not in row["Review status"]:
                # exclude those only listed for cancer, some are cancer and mito diseases keep those
                # there shouldn't be any however that are not 'no assertion criteria provided for release used
                # identified through manual inspection of conditions, those only annotated for cancer excluded
                if interp != "Conflicting interpretations of pathogenicity":
                    if (row["Condition(s)"] != "Familial colorectal cancer" and
                        row["Condition(s)"] != "Familial cancer of breast" and
                        row["Condition(s)"] != "Acute megakaryoblastic leukemia|Mediastinal germ cell tumor" and
                        row["Condition(s)"] != "Neoplasm of ovary"):
                        
                        dict[(ref, pos, alt)] = interp
    return dict


def chimp_ref_lookup(anno_file_path):
    """Parse an alignment of the human and chimpanzee reference mtDNA sequences and determine the ancestral chimp allele.
    Note this uses a shifted version of the human reference sequence.

    :return: a dictionary, where the key is a tuple of the ref and position, and the value is the chimp allele
    """
    # read in alignment file generated by the msa package in R
    # create dictionary to parse results
    dict = {}
    with open(anno_file_path+'/other_annotations/human-shifted_chimp_mt_aln.txt') as file:
        # mark position in the alignment file
        h_aln_pos = 1
        c_aln_pos = 1
        while True:
            line = file.readline()
            if not line:
                break
            if line.startswith('[1]'):  # these are the human sequence in the alignment
                sequence = line.strip().split(' ')[1]
                for base in sequence:
                    dict[('human', h_aln_pos)] = base
                    h_aln_pos += 1
            if line.startswith('[2]'):  # these are the chimp sequence in the alignment
                sequence = line.strip().split(' ')[1]
                for base in sequence:
                    dict[('chimp', c_aln_pos)] = base
                    c_aln_pos += 1

    # now parse dict to return chimp reference allele at each position
    # this is shifted human mtDNA, so starts at m.577 - this is to match the start position of the chimpanzee reference
    new_dict = {}
    pos = 577
    for key in dict:
        if (key[0] == 'human') and (dict[key] != '-'):  # these are gaps in human alignment
            aln_pos = key[1]
            human_ref = dict[key]
            chimp_ref = dict[('chimp', aln_pos)]
            new_dict[(human_ref, pos)] = chimp_ref
            pos += 1
            if pos == 16570:  # will happen to 16569
                pos = 1  # to renumber
                
    return new_dict


def remove_common_suffix(s1, s2):
    """
    Remove the common suffix of two strings.
    """
    i = 0
    while i < min(len(s1), len(s2)) and s1[-1 - i] == s2[-1 - i]:
        i += 1

    new_s1 = s1[:-i] if i > 0 else s1
    new_s2 = s2[:-i] if i > 0 else s2
    new_s1 = new_s1 if new_s1 else s1
    new_s2 = new_s2 if new_s2 else s2
    return (new_s1, new_s2)


def annotate(input_file, annotated_txt, annotated_vcf, anno_file_path):
    """Annotate the file with all possible mitochondrial mutations and their likelihood scores.

    :param input_file: the file with mutation likelihood scores, output of composite_likelihood_mito.py
    """
    # generate required dictionaries
    vep_anno_list = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", 
                     "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", 
                     "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", 
                     "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "HGVS_OFFSET"]

    rcrs_pos2trinuc = rcrs_pos_to_trinucleotide(anno_file_path)
    dbsnp = dbSNP_annotate(anno_file_path)
    gnomad = gnomad_annotate(anno_file_path)
    phylop = phylop_annotate(anno_file_path)
    vep = vep_annotate(anno_file_path, vep_anno_list)
    tRNA_position = tRNA_positions(anno_file_path)
    RNA_dom_mod = RNA_domains_mods(anno_file_path)
    RNA_type = RNA_base_type(anno_file_path)
    uniprot = uniprot_annotations(anno_file_path)
    other_prot = curated_func_sites(anno_file_path)
    apogee_scores = apogee(anno_file_path)
    mitotip_scores = mitotip(anno_file_path)
    hmtvar_scores = hmtvar(anno_file_path)
    helix = in_helix(anno_file_path)
    mitomap_genome_loci = mitomap_locus(anno_file_path)
    mitomap_vars1, mitomap_vars2 = mitomap(anno_file_path)
    clinvar_vars = clinvar(anno_file_path)
    chimp_dict = chimp_ref_lookup(anno_file_path)

    f = gzip.open(annotated_txt, "wt") if annotated_txt.endswith(".gz") else open(annotated_txt, "w")
    output_vcf = gzip.open(annotated_vcf, "wt") if annotated_vcf.endswith(".gz") else open(annotated_vcf, "w")
    header_list = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER', 'trinucleotide', 'AF', 'AC', 'AN', 'TOTAL_N', 'HOM_PF',
                   'HET_PF', 'SUM_PF', 'PT', 'mitomap_locus', 'symbol', 'transcript', 'feature', 'biotype', 
                   'consequence', 'impact', 'HGVSc', 'HGVSp', 'amino_acids', 'protein_position', 'codon_change', 
                   'gnomad_max_hl', 'gnomad_af_hom', 'gnomad_af_het', 'gnomad_ac_hom', 'gnomad_ac_het', 'in_phylotree', 
                   'phyloP_score', 'tRNA_position', 'tRNA_domain', 'RNA_base_type', 'RNA_modified', 'rRNA_bridge_base', 
                   'uniprot_annotation', 'other_prot_annotation', 'apogee_class', 'mitotip_class', 'hmtvar_class', 
                   'helix_max_hl', 'helix_af_hom', 'helix_af_het', 'mitomap_gbcnt', 'mitomap_af', 'mitomap_status', 
                   'mitomap_plasmy', 'mitomap_disease', 'clinvar_interp', 'chimp_ref']
    header = '\t'.join(header_list)
    f.write(header + '\n')

    with gzip.open(input_file, "rt") if input_file.endswith(".gz") else open(input_file, "r") as IN:
        ovlp_start_idx = 0
        for line in IN:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    output_vcf.write("##INFO=<ID=trinucleotide,Number=A,Type=String,Description=\"Position to its reference trinucleotide\">\n")
                    output_vcf.write("##INFO=<ID=symbol,Number=A,Type=String,Description=\"Gene symbol\">\n")
                    output_vcf.write("##INFO=<ID=mitomap_locus,Number=A,Type=String,Description=\"Mitochondrial DNA Function Locations in MITOMAP: https://www.mitomap.org/foswiki/bin/view/MITOMAP/GenomeLoci\">\n")
                    output_vcf.write(f"##INFO=<ID=VEP_CSQ,Number=.,Type=String,Description=\"Consequence annotations of Ensembl VEP. Format: {'|'.join(vep_anno_list)}\">\n")
                    output_vcf.write("##INFO=<ID=gnomad_max_hl,Number=A,Type=String,Description=\"Maximum observed mtDNA heteroplasmy (max_hl) in gnomAD\">\n")
                    output_vcf.write("##INFO=<ID=gnomad_af_hom,Number=A,Type=String,Description=\"Observed homoplasmy allele frequence in gnomAD\">\n")
                    output_vcf.write("##INFO=<ID=gnomad_af_het,Number=A,Type=String,Description=\"Observed heteroplasmy allele frequence in gnomAD\">\n")
                    output_vcf.write("##INFO=<ID=gnomad_ac_hom,Number=A,Type=String,Description=\"Observed homoplasmy allele count in gnomAD\">\n")
                    output_vcf.write("##INFO=<ID=gnomad_ac_het,Number=A,Type=String,Description=\"Observed heteroplasmy allele count in gnomAD\">\n")
                    output_vcf.write("##INFO=<ID=in_phylotree,Number=A,Type=String,Description=\"phylotree variants\">\n")
                    output_vcf.write("##INFO=<ID=phyloP_score,Number=A,Type=String,Description=\"phyloP scores for conservation\">\n")
                    output_vcf.write("##INFO=<ID=tRNA_position,Number=A,Type=String,Description=\"The tRNA position numbers (ranging from 1-73) encoded by mtDNA\">\n")
                    output_vcf.write("##INFO=<ID=tRNA_domain,Number=A,Type=String,Description=\"Modified tRNA domains encoded by mtDNA\">\n")
                    output_vcf.write("##INFO=<ID=RNA_base_type,Number=A,Type=String,Description=\"RNA base types encoded by mtDNA (Watson-Crick (WC), non-WC, or loop/other)\">\n")
                    output_vcf.write("##INFO=<ID=RNA_modified,Number=A,Type=String,Description=\"Modified RNA bases encoded by mtDNA\">\n")
                    output_vcf.write("##INFO=<ID=rRNA_bridge_base,Number=A,Type=String,Description=\"It's rRNA bridge base or not\">\n")
                    output_vcf.write("##INFO=<ID=uniprot_annotation,Number=A,Type=String,Description=\"The UniProt annotations for mtDNA-encoded proteins\">\n")
                    output_vcf.write("##INFO=<ID=other_prot_annotation,Number=A,Type=String,Description=\"Residues involved in complex I proton transfer, curated from PMID:32972993\">\n")
                    output_vcf.write("##INFO=<ID=apogee_class,Number=A,Type=String,Description=\"The APOGEE scores for missense\">\n")
                    output_vcf.write("##INFO=<ID=mitotip_class,Number=A,Type=String,Description=\"The MitoTip in silico scores for tRNA variants\">\n")
                    output_vcf.write("##INFO=<ID=hmtvar_class,Number=A,Type=String,Description=\"The HmtVar in silico scores for tRNA variants\">\n")
                    output_vcf.write("##INFO=<ID=helix_max_hl,Number=A,Type=String,Description=\"Maximum observed heteroplasmy (max_hl) in HelixMTdb\">\n")
                    output_vcf.write("##INFO=<ID=helix_af_hom,Number=A,Type=String,Description=\"Observed homoplasmy allele frequence in HelixMTdb\">\n")
                    output_vcf.write("##INFO=<ID=helix_af_het,Number=A,Type=String,Description=\"Observed heteroplasmy allele frequence in HelixMTdb\">\n")
                    output_vcf.write("##INFO=<ID=mitomap_gbcnt,Number=A,Type=String,Description=\"Total number of gb sequences to convert to allele freq in MITOMAP\">\n")
                    output_vcf.write("##INFO=<ID=mitomap_af,Number=A,Type=String,Description=\"AF in MITOMAP\">\n")
                    output_vcf.write("##INFO=<ID=mitomap_status,Number=A,Type=String,Description=\"Disease association status for variant in MITOMAP\">\n")
                    output_vcf.write("##INFO=<ID=mitomap_plasmy,Number=A,Type=String,Description=\"Plasmy in MITOMAP\">\n")
                    output_vcf.write("##INFO=<ID=mitomap_disease,Number=A,Type=String,Description=\"Disease in MITOMAP\">\n")
                    output_vcf.write("##INFO=<ID=clinvar_interp,Number=A,Type=String,Description=\"The clinical significance interpretation for variant in ClinVar\">\n")
                    output_vcf.write("##INFO=<ID=chimp_ref,Number=A,Type=String,Description=\"The ancestral chimpanzee allele for thsi position\">\n")
                    output_vcf.write(f"##annotate_command=python {' '.join(sys.argv)}\n")

                output_vcf.write(line)
                continue
            
            CHROM, POS, ID, REFs, ALTs, QUAL, FILTER, INFO, FORMAT, *SAMPLES = line.strip().split('\t')
            if POS not in rcrs_pos2trinuc:  # Skip the position of ref=='N'
                continue
            
            # 这一步是必须的，因为有些变异的REF和ALT由于多突变的原因后缀是相同的
            REF_ALT_list = [remove_common_suffix(REF, ALT) for REF, ALT in list(
                zip_longest(REFs.split(','), ALTs.split(','), fillvalue=REFs.split(',')[0]))]

            variant_list   = [REF + POS + ALT for REF, ALT in REF_ALT_list]
            var_tuple_list = [(REF, POS, ALT) for REF, ALT in REF_ALT_list]
            
            IDs = []
            for REF, ALT in REF_ALT_list:
                # check if the variant is in dbSNP
                if (CHROM, POS, REF, ALT) in dbsnp:
                    IDs.append(dbsnp[(CHROM, POS, REF, ALT)])
            ID = ','.join(IDs) if IDs else '.'
                
            # POS in mitomap_genome_loci or not
            mitomap_locus_id = ''
            for i in range(ovlp_start_idx, len(mitomap_genome_loci)):
                start, end, locus = mitomap_genome_loci[i]
                if int(POS) > end: continue
                if int(POS) < start: break
                
                ovlp_start_idx   = i
                mitomap_locus_id = locus
                break
            
            in_phylo_list = [1 if "\n" + variant + "\n" in open(anno_file_path+'/databases/phylotree_variants.txt').read() else 0 for variant in variant_list]
            max_hl_list   = [gnomad[var_tuple][0] if var_tuple in gnomad else 0 for var_tuple in var_tuple_list]
            
            vep_csq_list     = []
            vep_symbol_list  = []
            vep_gene_list    = []
            vep_feature_list = []
            vep_biotype_list = []
            vep_conseq_list  = []
            vep_impact_list  = []
            vep_hgvsc_list   = []
            vep_hgvsp_list   = []
            vep_aa_list      = []
            vep_pp_list      = []
            vep_codon_change_list = []
            
            for REF, POS, alt in var_tuple_list:
                if (REF, POS, alt) in vep:
                    vep_csq_list.append("|".join([vep[(REF, POS, alt)][label] for label in vep_anno_list]))

                    vep_symbol_list.append(vep[(REF, POS, alt)]["SYMBOL"])
                    vep_gene_list.append(vep[(REF, POS, alt)]["Gene"])
                    vep_feature_list.append(vep[(REF, POS, alt)]["Feature"])
                    vep_biotype_list.append(vep[(REF, POS, alt)]["BIOTYPE"])
                    vep_conseq_list.append(vep[(REF, POS, alt)]["Consequence"])
                    vep_impact_list.append(vep[(REF, POS, alt)]["IMPACT"])
                    vep_hgvsc_list.append(vep[(REF, POS, alt)]["HGVSc"].split(':')[-1])  # ENST00000361390.2:c.271A>C
                    vep_hgvsp_list.append(vep[(REF, POS, alt)]["HGVSp"].split(':')[-1])  # ENSP00000354189.1:p.Thr91Pro
                    vep_aa_list.append(vep[(REF, POS, alt)]["Amino_acids"])
                    vep_pp_list.append(vep[(REF, POS, alt)]["Protein_position"])
                    vep_codon_change_list.append(vep[(REF, POS, alt)]["Codons"])
            
            gnomad_af_hom_list = [gnomad[var_tuple][1] if var_tuple in gnomad else 0 for var_tuple in var_tuple_list]
            gnomad_af_het_list = [gnomad[var_tuple][2] if var_tuple in gnomad else 0 for var_tuple in var_tuple_list]
            gnomad_ac_hom_list = [gnomad[var_tuple][3] if var_tuple in gnomad else 0 for var_tuple in var_tuple_list]
            gnomad_ac_het_list = [gnomad[var_tuple][4] if var_tuple in gnomad else 0 for var_tuple in var_tuple_list]
            
            tRNA_pos = tRNA_position[POS] if POS in tRNA_position else ''
            tRNA_dom = RNA_dom_mod[(POS, "domain")] if (POS, "domain") in RNA_dom_mod else ''
            RNA_mod  = RNA_dom_mod[(POS, "modified")] if (POS, "modified") in RNA_dom_mod else ''
            RNA_base = str(RNA_type[POS]).strip('[]').replace("'", "").replace(" ", "") if POS in RNA_type else ''
            RNA_bridge = "Yes" if ("\n" + POS + "\n") in open(anno_file_path+'/other_annotations/rRNA_bridge_bases.txt').read() else "No"

            uniprot_annot = str(uniprot[int(POS)]).strip('[]').replace("'", "").replace(" ", "") if int(POS) in uniprot else ''
            other_prot_annot = str(other_prot[int(POS)]).strip('[]').replace("'", "").replace(" ", "") if int(POS) in other_prot else ''
            apogee_score_list = [str(apogee_scores[var_tuple]).strip('[]').replace("'", "").replace(" ", "") if var_tuple in apogee_scores else '' for var_tuple in var_tuple_list]
            mitotip_score_list = [mitotip_scores[var_tuple] if var_tuple in mitotip_scores else '' for var_tuple in var_tuple_list]
            
            helix_max_hl_list = [helix[var_tuple][0] if var_tuple in helix else 0 for var_tuple in var_tuple_list]
            helix_af_hom_list = [helix[var_tuple][1] if var_tuple in helix else 0 for var_tuple in var_tuple_list]
            helix_af_het_list = [helix[var_tuple][2] if var_tuple in helix else 0 for var_tuple in var_tuple_list]
            
            mitomap_ac_list     = [mitomap_vars2[var_tuple][0] if var_tuple in mitomap_vars2 else 0 for var_tuple in var_tuple_list]
            mitomap_af_list     = [mitomap_vars2[var_tuple][1] if var_tuple in mitomap_vars2 else 0 for var_tuple in var_tuple_list]
            mitomap_status_list = [mitomap_vars1[var_tuple][0] if var_tuple in mitomap_vars1 else '' for var_tuple in var_tuple_list]
            mitomap_plasmy_list = [(mitomap_vars1[var_tuple][1] + '/' + mitomap_vars1[var_tuple][2]) if var_tuple in mitomap_vars1 else '' for var_tuple in var_tuple_list]
            mitomap_dz_list     = [mitomap_vars1[var_tuple][3] if var_tuple in mitomap_vars1 else '' for var_tuple in var_tuple_list]
            clinvar_int_list    = [clinvar_vars[var_tuple] if var_tuple in clinvar_vars else '' for var_tuple in var_tuple_list]
            hmtvar_scores_list  = [str(hmtvar_scores[(REF, POS, alt)]).strip('[]').replace("'", "").replace(" ", "") 
                                   if (REF, POS, alt) in hmtvar_scores else '' for REF, POS, alt in var_tuple_list]

            tRNA_pos_str  = str(tRNA_pos).strip('[]').replace("'", "").replace(" ", "")
            tRNA_dom_str  = str(tRNA_dom).strip('[]').replace("'", "").replace(" ", "")
            RNA_mod_str   = str(RNA_mod).strip('[]').replace("'", "").replace(" ", "")
            
            info_dict = {}
            for tab in INFO.split(";"):
                v = tab.split("=")
                info_dict[v[0]] = v[1] if len(v) > 1 else None

            for i, (ref, pos, alt) in enumerate(var_tuple_list):
                chimp_ref_str = "".join([chimp_dict[(ref[c], int(pos)+c)] for c in range(len(ref))])
                rsid = dbsnp[(CHROM, pos, ref, alt)] if (CHROM, pos, ref, alt) in dbsnp else '.'
                af = info_dict['AF'].split(',')[i]
                ac = info_dict['AC'].split(',')[i]
                an = info_dict['AN']
                f.write('\t'.join([CHROM, pos, rsid, ref, alt, FILTER, rcrs_pos2trinuc[pos], af, ac, an]) + '\t' + 
                        '\t'.join([info_dict['Total_N'], info_dict['HOM_PF'], info_dict['HET_PF']])  + '\t' +
                        '\t'.join([info_dict['SUM_PF'], info_dict['PT'], mitomap_locus_id])          + '\t' +
                        (vep[(ref, pos, alt)]["SYMBOL"]           if (ref, pos, alt) in vep else "") + '\t' +
                        (vep[(ref, pos, alt)]["Gene"]             if (ref, pos, alt) in vep else "") + '\t' +
                        (vep[(ref, pos, alt)]["Feature"]          if (ref, pos, alt) in vep else "") + '\t' +
                        (vep[(ref, pos, alt)]["BIOTYPE"]          if (ref, pos, alt) in vep else "") + '\t' +
                        (vep[(ref, pos, alt)]["Consequence"]      if (ref, pos, alt) in vep else "") + '\t' + 
                        (vep[(ref, pos, alt)]["IMPACT"]           if (ref, pos, alt) in vep else "") + '\t' +
                        (vep[(ref, pos, alt)]["HGVSc"]            if (ref, pos, alt) in vep else "") + '\t' +
                        (vep[(ref, pos, alt)]["HGVSp"]            if (ref, pos, alt) in vep else "") + '\t' +
                        (vep[(ref, pos, alt)]["Amino_acids"]      if (ref, pos, alt) in vep else "") + '\t' +
                        (vep[(ref, pos, alt)]["Protein_position"] if (ref, pos, alt) in vep else "") + '\t' +
                        (vep[(ref, pos, alt)]["Codons"]           if (ref, pos, alt) in vep else "") + '\t' +
                        str(max_hl_list[i])                                                          + '\t' +
                        str(gnomad_af_hom_list[i])                                                   + '\t' +
                        str(gnomad_af_het_list[i])                                                   + '\t' +
                        str(gnomad_ac_hom_list[i])                                                   + '\t' +
                        str(gnomad_ac_het_list[i])                                                   + '\t' +
                        str(in_phylo_list[i])                                                        + '\t' +
                        phylop[int(POS)]                                                             + '\t' +
                        tRNA_pos_str                                                                 + '\t' +
                        tRNA_dom_str                                                                 + '\t' +
                        RNA_base                                                                     + '\t' +
                        RNA_mod_str                                                                  + '\t' +
                        RNA_bridge                                                                   + '\t' +
                        uniprot_annot                                                                + '\t' +
                        other_prot_annot                                                             + '\t' +
                        apogee_score_list[i]                                                         + '\t' +
                        mitotip_score_list[i]                                                        + '\t' +
                        hmtvar_scores_list[i]                                                        + '\t' +
                        str(helix_max_hl_list[i])                                                    + '\t' +
                        str(helix_af_hom_list[i])                                                    + '\t' +
                        str(helix_af_het_list[i])                                                    + '\t' +
                        str(mitomap_ac_list[i])                                                      + '\t' +
                        str(mitomap_af_list[i])                                                      + '\t' +
                        mitomap_status_list[i]                                                       + '\t' +
                        mitomap_plasmy_list[i]                                                       + '\t' +
                        mitomap_dz_list[i]                                                           + '\t' +
                        clinvar_int_list[i]                                                          + '\t' +
                        chimp_ref_str                                                                + '\n')

            # f.write(CHROM + '\t' + POS  + '\t' + ID        + '\t' + 
            #         REFs  + '\t' + ALTs + '\t' + INFO      + '\t' +
            #         rcrs_pos2trinuc[POS]                   + '\t' +
            #         ','.join(vep_symbol_list)              + '\t' + 
            #         mitomap_locus_id                       + '\t' +
            #         ','.join(vep_conseq_list)              + '\t' +
            #         ','.join(vep_aa_list)                  + '\t' +
            #         ','.join(vep_pp_list)                  + '\t' +
            #         ','.join(vep_codon_change_list)        + '\t' +
            #         ','.join(map(str, max_hl_list))        + '\t' +
            #         ','.join(map(str, gnomad_af_hom_list)) + '\t' +
            #         ','.join(map(str, gnomad_af_het_list)) + '\t' +
            #         ','.join(map(str, gnomad_ac_hom_list)) + '\t' +
            #         ','.join(map(str, gnomad_ac_het_list)) + '\t' +
            #         ','.join(map(str, in_phylo_list))      + '\t' +
            #         phylop[int(POS)]                       + '\t' +
            #         tRNA_pos_str                           + '\t' +
            #         tRNA_dom_str                           + '\t' +
            #         RNA_base                               + '\t' +
            #         RNA_mod_str                            + '\t' +
            #         str(RNA_bridge)                        + '\t' +
            #         uniprot_annot                          + '\t' +
            #         other_prot_annot                       + '\t' +
            #         ','.join(apogee_score_list)            + '\t' +
            #         ','.join(map(str, mitotip_score_list)) + '\t' +
            #         ','.join(map(str, hmtvar_scores_list)) + '\t' +
            #         ','.join(map(str, helix_max_hl_list))  + '\t' +
            #         ','.join(map(str, helix_af_hom_list))  + '\t' +
            #         ','.join(map(str, helix_af_het_list))  + '\t' +
            #         ','.join(map(str, mitomap_ac_list))    + '\t' +
            #         ','.join(map(str, mitomap_af_list))    + '\t' +
            #         ','.join(mitomap_status_list)          + '\t' +
            #         ','.join(mitomap_plasmy_list)          + '\t' +
            #         ','.join(mitomap_dz_list)              + '\t' +
            #         ','.join(clinvar_int_list)             + '\t' +
            #         chimp_ref_str                          + '\n')

            chimp_ref_str = "".join([chimp_dict[(REFs[i], int(POS)+i)] for i in range(len(REFs))])
            anno_info = "trinucleotide="         + rcrs_pos2trinuc[POS]                                + ';' + \
                        "symbol="                + ','.join(vep_symbol_list)                           + ';' + \
                        "mitomap_locus="         + mitomap_locus_id                                    + ';' + \
                        "VEP_CSQ="               + ','.join(vep_csq_list)                              + ';' + \
                        "gnomad_max_hl="         + ','.join(map(str, max_hl_list))                     + ';' + \
                        "gnomad_af_hom="         + ','.join(map(str, gnomad_af_hom_list))              + ';' + \
                        "gnomad_af_het="         + ','.join(map(str, gnomad_af_het_list))              + ';' + \
                        "gnomad_ac_hom="         + ','.join(map(str, gnomad_ac_hom_list))              + ';' + \
                        "gnomad_ac_het="         + ','.join(map(str, gnomad_ac_het_list))              + ';' + \
                        "in_phylotree="          + ','.join(map(str, in_phylo_list))                   + ';' + \
                        "phyloP_score="          + phylop[int(POS)]                                    + ';' + \
                        "tRNA_position="         + tRNA_pos_str                                        + ';' + \
                        "tRNA_domain="           + tRNA_dom_str                                        + ';' + \
                        "RNA_base_type="         + RNA_base                                            + ';' + \
                        "RNA_modified="          + RNA_mod_str                                         + ';' + \
                        "rRNA_bridge_base="      + RNA_bridge                                          + ';' + \
                        "uniprot_annotation="    + uniprot_annot                                       + ';' + \
                        "other_prot_annotation=" + other_prot_annot                                    + ';' + \
                        "apogee_class="          + ','.join([r for r in apogee_score_list if r!=''])   + ';' + \
                        "mitotip_class="         + ','.join([r for r in mitotip_score_list if r!=''])  + ';' + \
                        "hmtvar_class="          + ','.join([r for r in hmtvar_scores_list if r!=''])  + ';' + \
                        "helix_max_hl="          + ','.join(map(str, helix_max_hl_list))               + ';' + \
                        "helix_af_hom="          + ','.join(map(str, helix_af_hom_list))               + ';' + \
                        "helix_af_het="          + ','.join(map(str, helix_af_het_list))               + ';' + \
                        "mitomap_gbcnt="         + ','.join(map(str, mitomap_ac_list))                 + ';' + \
                        "mitomap_af="            + ','.join(map(str, mitomap_af_list))                 + ';' + \
                        "mitomap_status="        + ','.join([r for r in mitomap_status_list if r!='']) + ';' + \
                        "mitomap_plasmy="        + ','.join([r for r in mitomap_plasmy_list if r!='']) + ';' + \
                        "mitomap_disease="       + ','.join([r for r in mitomap_dz_list if r!=''])     + ';' + \
                        "clinvar_interp="        + ','.join([r for r in clinvar_int_list if r!=''])    + ';' + \
                        "chimp_ref="             + chimp_ref_str
                        
            INFO = INFO.strip() + ';' + anno_info
            line = '\t'.join([CHROM, POS, ID, REFs, ALTs, QUAL, FILTER, INFO, FORMAT, *SAMPLES])
            output_vcf.write(line + '\n')
    
    f.close()
    output_vcf.close()


def main():
    parser = argparse.ArgumentParser(description='To annotate the VCF file using the annotated files')
    parser.add_argument('-i', "--input", type=str, help="The input VCF file")
    parser.add_argument('-d', "--directory", type=str, help="The directory containing the annotated files")
    parser.add_argument('-o', "--annotated_txt", type=str, help="The output annotated file in txt format")
    parser.add_argument('-v', "--annotated_vcf", type=str, help="The output annotated file in vcf format containing "
                                                                "the information about all samples")
    args = parser.parse_args()
    annotate(input_file=args.input,
             annotated_txt=args.annotated_txt,
             annotated_vcf=args.annotated_vcf,
             anno_file_path=args.directory)


if __name__ == "__main__":
    print(datetime.datetime.now(), "Annotating the mitochondrial variants.")
    main()
    print(datetime.datetime.now(), "All done!" + '\n')
