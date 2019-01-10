

from Bio import SeqIO
import os
import sys

filename = sys.argv[1]
out_directory = sys.argv[2]
filebase = os.path.splitext(filename)[0]

handle = open(filename, 'rU')


def extract_cds_genes():
    for record in SeqIO.parse(handle, 'genbank') :
        seq = str(record.seq)
        for feature in record.features:
            print feature.type
            print feature.location

            if feature.type == 'CDS' or feature.type == 'rRNA':

                if 'gene' in feature.qualifiers:
                    geneName = feature.qualifiers['gene'][0]
                elif 'product' in feature.qualifiers:
                    geneName = feature.qualifiers['product'][0]
                else:
                    print 'ERROR when parsing feature:'
                    print feature.qualifiers
                    quit()

                geneName = geneName.replace(' ', '-');
                print geneName

                geneFile = open(geneName + '.fa', 'w+')
                geneFile.write('>')
                geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                geneFile.write(seq[feature.location.start.position:feature.location.end.position])
                geneFile.write("\n\n")
            print '----------------------------'







def extract_individual_gene():
    gene_array = ["fusA", "gyrB", "leuS", "pyrG", "rplB", "rpoB"]
    genome_name = os.path.basename(filebase)
    print genome_name
    absent_genomes = []
    for record in SeqIO.parse(handle, 'genbank') :
            seq = str(record.seq)
            for feature in record.features:
                if feature.type == 'CDS':
                    #print "CDS"
                    #print feature
                    if 'gene' not in feature.qualifiers:
                        if filename not in absent_genomes:
                            absent_genomes.append(filename)
                        continue
                        #print "No gene"
                    if feature.qualifiers['gene'][0] == "fusA_1":
                        #print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".fa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(seq[feature.location.start.position:feature.location.end.position])
                        geneFile.write("\n\n")
		    if feature.qualifiers['gene'][0] == "fusA_2":
                        #print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".fa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(seq[feature.location.start.position:feature.location.end.position])
                        geneFile.write("\n\n")
                    if feature.qualifiers['gene'][0] == "pyrG_1":
                        #print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".fa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(seq[feature.location.start.position:feature.location.end.position])
                        geneFile.write("\n\n")
                    if feature.qualifiers['gene'][0] == "pyrG_2":
                        #print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".fa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(seq[feature.location.start.position:feature.location.end.position])
                        geneFile.write("\n\n")

    print absent_genomes

#extract_individual_gene()



def extract_individual_protein():
    gene_array = ["fusA", "gyrB", "leuS", "pyrG", "rplB", "rpoB"]
    genome_name = os.path.basename(filebase)
    print genome_name
    absent_genomes = []
    for record in SeqIO.parse(handle, 'genbank') :
            seq = str(record.seq)
            for feature in record.features:
                if feature.type == 'CDS':
                    #print "CDS"
                    #print feature
                    if 'gene' not in feature.qualifiers:
                        if filename not in absent_genomes:
                            absent_genomes.append(filename)
                        continue
                        #print "No gene"
                    if feature.qualifiers['gene'][0] == "dnaA":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")
		    if feature.qualifiers['gene'][0] == "fusA_1":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")
                    if feature.qualifiers['gene'][0] == "gyrB":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")
		    if feature.qualifiers['gene'][0] == "leuS":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")
		    if feature.qualifiers['gene'][0] == "pyrG":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")
                    if feature.qualifiers['gene'][0] == "rplB":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")
                    if feature.qualifiers['gene'][0] == "rpoB":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")
		    if feature.qualifiers['gene'][0] == "pyrG_2":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")
                    if feature.qualifiers['gene'][0] == "fusA_2":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
                        geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")
    print absent_genomes


#extract_individual_protein()





def extract_individual_protein_2():
    gene_array = ["fusA", "gyrB", "leuS", "pyrG", "rplB", "rpoB"]
    genome_name = os.path.basename(filebase)
    print genome_name
    absent_genomes = []
    for record in SeqIO.parse(handle, 'genbank') :
            seq = str(record.seq)
            for feature in record.features:
                if feature.type == 'CDS':
                    #print "CDS"
                    #print feature
                    if 'gene' not in feature.qualifiers:
                        if filename not in absent_genomes:
                            absent_genomes.append(filename)
                        continue
                        #print "No gene"
                    if feature.qualifiers['gene'][0] == "pyrG_1":
                        print feature
                        geneName = feature.qualifiers['gene'][0]
                        geneName = geneName.replace(' ', '-')
                        translation_seq = feature.qualifiers['translation'][0]
			geneFilename = out_directory + "/" + genome_name +  "_" + geneName + ".faa"
                        geneFile = open(geneFilename, 'w+')
                        geneFile.write('>')
                        geneFile.write(os.path.basename(filebase) + '_' + geneName + "\n")
                        geneFile.write(translation_seq)
                        geneFile.write("\n\n")


extract_individual_protein_2()
