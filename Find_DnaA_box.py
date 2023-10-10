def minimum_skew(genome):
    skew_list = [0]
    for i in range(len(genome)):
        nucleotide = genome[i]
        if nucleotide == 'C':
            skew_list.append(skew_list[-1] - 1)
        elif nucleotide == 'G':
            skew_list.append(skew_list[-1] + 1)
        else:
            skew_list.append(skew_list[-1])
    min_value = min(skew_list)
    index_list = [i for i, x in enumerate(skew_list) if x == min_value]
    index_list_string = ''
    for _ in index_list:
        index_list_string += str(_) + ' '
    return index_list_string.rstrip()

f = open('Salmonella_enterica.txt', 'r')
genome_text = f.read().replace('\n', '')
window_start = minimum_skew(genome_text).split()[0]
window = genome_text[int(window_start):((int(window_start)) + 500)]

def maxmap(worddict):
    frequentpatterns = []
    max_count = max(worddict.values())
    for key, value in worddict.items():
        if worddict[key] == max_count:
            frequentpatterns.append(key)
    return frequentpatterns

def hamming_distance(p, q):
    list1 = list(p)
    list2 = list(q)
    count = 0
    for i in range(len(list1)):
        if list1[i] != list2[i]:
            count += 1
    return count

def suffix(pattern):
    suffix_pattern = pattern[1:]
    return suffix_pattern

def first_symbol(pattern):
    first_symbol = pattern[0]
    return first_symbol

def neighbors(pattern, d):
    nucleotides = ['A', 'T', 'G', 'C']
    if d == 0:
        neighborhood =[]
        neighborhood.append(pattern)
        return neighborhood
    elif len(pattern) == 1:
        return nucleotides
    neighborhood = []
    suffix_neighbors = neighbors(suffix(pattern), d)
    for string_text in suffix_neighbors:
        if hamming_distance(suffix(pattern), string_text) < d:
            for nucleotide in nucleotides:
                neighbor = nucleotide + string_text
                neighborhood.append(neighbor)
        else:
            neighbor = first_symbol(pattern) + string_text
            neighborhood.append(neighbor)
    return neighborhood


def reversecomplement(forward_pattern):
    large_to_small = {'A':'a',
               'T':'t',
               'G':'g',
               'C':'c'
               }
    for key, value in large_to_small.items():
        if key in forward_pattern:
            forward_pattern = forward_pattern.replace(key, value)
    complementary = {'a':'T',
               't':'A',
               'g':'C',
               'c':'G'
               }
    for key, value in complementary.items():
        if key in forward_pattern:
            forward_pattern = forward_pattern.replace(key, value)
    reversecomplementary = forward_pattern[::-1].strip()
    return reversecomplementary

def frequent_words_with_mismatches(text, k, d):
    patterns = []
    freqmap = {}
    n = len(text)
    for i in range(n-k+1):
        pattern = text[i:(i+k)]
        reverse_pattern = reversecomplement(pattern)
        neighborhood = neighbors(pattern, d)
        reverse_neighborhood = neighbors(reverse_pattern, d)
        for j in range(len(neighborhood)):
            neighbor = neighborhood[j]
            if neighbor not in freqmap:
                freqmap[neighbor] = 1
            else:
                freqmap[neighbor] += 1
        for j in range(len(reverse_neighborhood)):
            reverse_neighbor = reverse_neighborhood[j]
            if reverse_neighbor not in freqmap:
                freqmap[reverse_neighbor] = 1
            else:
                freqmap[reverse_neighbor] += 1
    m = maxmap(freqmap)
    for pattern, count in freqmap.items():
        if pattern in m:
            patterns.append(pattern)
    return patterns



x = frequent_words_with_mismatches(window, 9, 1)
print(x)

new = open('DnaA box frequent kmers.txt', 'w')
for element in x:
    new.write(element + ' ')
new.close()
