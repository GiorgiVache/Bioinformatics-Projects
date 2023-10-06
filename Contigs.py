ASCII_to_qualityscore_dict = {!}

def readfastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences[:100]

def contig_generation(reads, y, x):
    # y = minimum number of same nucleotides
    # x = maximum number of same nucleotides
    with open('Contigs.txt', 'w') as c:
        try:
            for read1 in reads:
                for read2 in reads:
                    if read1 == read2:
                        continue
                    elif read1[0:x] == read2[-x:]:
                        contig = read2[:-x] + read1
                        print(contig)
                        reads.remove(read1)
                        reads.remove(read2)
                        reads.append(contig)
                        break
            if x > y:
                contig_generation(reads, y, (x - 1))
            for contig in reads:
                c.write(contig + '\n')
        except KeyboardInterrupt:
            c.close()
    return reads

contig_generation((readfastq('SRR396636.sra_1.fastq')), 3, 6)
