"""The program reads a fastq file and from the reads it generates contigs and stores them in an sql database"""
import sqlite3


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
    return sequences[:50]


conn = sqlite3.connect("contigs.db")
cursor = conn.cursor()
cursor.execute('''DROP TABLE IF EXISTS reads''')
cursor.execute('''CREATE TABLE IF NOT EXISTS reads (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    read TEXT,
                    length TEXT
                 )''')
for sequence in readfastq('SRR396636.sra_1.fastq'):
    cursor.execute("INSERT INTO reads (read, length) VALUES (?, ?)", (sequence, len(sequence)))
conn.commit()

cursor.execute('''DROP TABLE IF EXISTS contigs''')
cursor.execute('''CREATE TABLE IF NOT EXISTS contigs (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    ID_from_reads_1 TEXT,
                    ID_from_reads_2 TEXT,
                    contig TEXT,
                    length TEXT
                 )''')


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
                        cursor.execute('SELECT id FROM reads WHERE read = ? ', (read1,))
                        read1_id = cursor.fetchone()[0]
                        cursor.execute('SELECT id FROM reads WHERE read = ? ', (read2,))
                        read2_id = cursor.fetchone()[0]
                        cursor.execute("INSERT INTO contigs (ID_from_reads_1, ID_from_reads_2, contig, length) VALUES "
                                       "(?, ?, ?, ?)",
                                       (read1_id, read2_id, contig, len(contig)))
                        conn.commit()
                        reads.remove(read1)
                        reads.remove(read2)
                        reads.append(contig)
                        cursor.execute("INSERT INTO reads (read, length) VALUES (?, ?)", (contig, len(contig)))
                        break
            if x > y:
                contig_generation(reads, y, (x - 1))
            for contig in reads:
                c.write(contig + '\n')
        except KeyboardInterrupt:
            c.close()
    conn.close()
    return reads




contig_generation(readfastq('SRR396636.sra_1.fastq'), 3, 5)
