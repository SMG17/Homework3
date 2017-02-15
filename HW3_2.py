# 2-b)
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
			chromosome = " ".join((chromosome, "+"))
			names.append(chromosome)

	names_line_num.append(len(genome))
	sequences = []
	for i in range(len(names_line_num)-1):
		seq = [''.join(genome[names_line_num[i]+1:names_line_num[i+1]])]
		sequences.extend(seq)

	dictionary = dict(zip(names, sequences))
	return dictionary 

x = open("/mnt/c/Users/SMG/Desktop/my_genome.fa")
genome_dic_reg = fasta_read(x)

# 2-c)
def reverse_complement(dictionary):
	reverse_complement_seq = []
	complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
	for chromosome in dictionary:
		reverse_complement = [''.join(complement.get(base,base) for base in reversed(dictionary[chromosome]))]
		reverse_complement_seq.extend(reverse_complement)
	dictionary_rev = dict(zip(dictionary.keys(),reverse_complement_seq))
	dictionary_rev = {key.replace("+","-"): dictionary_rev[key] for key in dictionary_rev.keys()}
	return dictionary_rev

genome_dic_rev = reverse_complement(genome_dic_reg)

# 2-d)
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

def read_search(genome, read):
	if genome.find(read) == -1:
		return None
	else:
		return genome.find(read)

y = open("/mnt/c/Users/SMG/Desktop/my_reads.fastq")
reads_dic = read_fastq(y)

# 2-e)
myOutput = open("Reads_map.txt","w")
whole_genome = genome_dic_reg.copy()
whole_genome.update(genome_dic_rev)

for name, read in reads_dic.items():
	for chromosome, seq in whole_genome.items():
		position = read_search(seq, read)
		if position != None:
			myOutput.write("Read %s mapped to %s at position %s.\n" % (name, chromosome, position))

# 2-f)
print "The script took 1m20s to run."

x.close()
y.close()
myOutput.close()	
