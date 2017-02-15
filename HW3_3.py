# From 2)
def fasta_read(genome_file):
	names_line_num = []
	names = []
	genome = []
	for i, line in enumerate(genome_file):
		line = line.strip()
		genome.append(line)
		if line.startswith(">"):
			names_line_num.append(i)
			chromosome = line[1:]
			chromosome = chromosome.replace("_"," ")
			chromosome = " ".join((chromosome, "+"))
			names.append(chromosome)

	names_line_num.append(len(genome))
	sequences = []
	for i in range(len(names_line_num)-1):
		seq = [''.join(genome[names_line_num[i]+1:names_line_num[i+1]])]
		sequences.extend(seq)

	dictionary = dict(zip(names, sequences))
	return dictionary 

def reverse_complement(dictionary):
	reverse_complement_seq = []
	complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
	for chromosome in dictionary:
		reverse_complement = [''.join(complement.get(base,base) for base in reversed(dictionary[chromosome]))]
		reverse_complement_seq.extend(reverse_complement)
	dictionary_rev = dict(zip(dictionary.keys(),reverse_complement_seq))
	dictionary_rev = {key.replace("+","-"): dictionary_rev[key] for key in dictionary_rev.keys()}
	return dictionary_rev

def read_fastq(reads_file):
	name_list = []
	seq_list = []
	name2_list = []
	qual_list = []
	while True:
		name = reads_file.readline()
		if name == "":
			break
		name_list.append(name.strip())
		seq = reads_file.readline()
		seq_list.append(seq.strip())
		name2 = reads_file.readline()
		name2_list.append(name2.strip())
		qual = reads_file.readline()
		qual_list.append(qual.strip())
	dictionary_reads = dict(zip(name_list, seq_list))
	return dictionary_reads

# 3-a)
def preprocessing(genome_dic, reads_dic):
	genome_hashes = []
	chromo_position_hashes = []
	for chromosome, seq in genome_dic.items():
		m = len(reads_dic["@SOMEREAD_0"])
		for n in range(len(seq)-m+1):
			genome_hashes.append(seq[n:(n+m)])
			chromo_position_hashes.append("%s_%s" %(chromosome, n))
	dictionary_hashes = dict(zip(genome_hashes, chromo_position_hashes))
	return dictionary_hashes

x = open("/mnt/c/Users/SMG/Desktop/my_genome.fa")
genome_dic_reg = fasta_read(x)
genome_dic_rev = reverse_complement(genome_dic_reg)

y = open("/mnt/c/Users/SMG/Desktop/my_reads.fastq")
reads_dic = read_fastq(y)

a = preprocessing(genome_dic_reg, reads_dic)
b = preprocessing(genome_dic_rev, reads_dic)
preprocessed_genome = a.copy()
preprocessed_genome.update(b)

# 3-b)
myOutput = open("Reads_map2.txt","w")

for name, read in reads_dic.items():
	if preprocessed_genome.has_key(read):
		chromosome, sep, position = preprocessed_genome[read].partition("_")
		myOutput.write("Read %s mapped to %s at position %s.\n" % (name, chromosome, position))

# 3-c)
print "The script took 9 seconds to run."

# 3-d)
print "The big O time complexity is O(n), where n is the length of the dictionary of reads (reads_dic)."

# 3-e)
print "The big O space complexity is O(n*m), where n is the length of the dictionary of reads (reads_dic) and m is the length of the dictionary of hashed genome (preprocessed_genome)" 

x.close()
y.close()
myOutput.close()	
