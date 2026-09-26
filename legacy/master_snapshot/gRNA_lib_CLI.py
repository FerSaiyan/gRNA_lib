def find_gRNA(gene_sequence, PAM):
    PAM_idx_list = []
    PAM_list = []
    gRNA_list = []
    #print("Searching for PAM:", PAM)
    for i in range(len(gene_sequence)):
        #aux = gene_sequence[i+1:i+len(PAM)]
        #print(aux)
        
        if gene_sequence[i+1:i+len(PAM)] == PAM[1:]:
            PAM_idx_list.append(i)
            #print("PAM found at index", i)
    for j in PAM_idx_list:
        PAM_list.append([j, gene_sequence[j:j+3]])
    
    for k,l in PAM_list:
        if k < 20:
            continue
        gRNA_list.append([gene_sequence[k-20:k], k, l])
    return gRNA_list

def complimentary_sequence(gene_sequence):
    complimentary_sequence = []
    for i in range(len(gene_sequence)):
        if gene_sequence[i] == 'A':
            complimentary_sequence.append('T')
        elif gene_sequence[i] == 'T':
            complimentary_sequence.append('A')
        elif gene_sequence[i] == 'C':
            complimentary_sequence.append('G')
        elif gene_sequence[i] == 'G':
            complimentary_sequence.append('C')
    return ''.join(complimentary_sequence[::-1])

def gRNA_ranking(gRNA_list, genome_sequence, Sp_cas9=True):
    ranking = {}
    off_target = {}
    for n_gRNA in range(len(gRNA_list)):
        gRNA = gRNA_list[n_gRNA]
        off_target[gRNA] = genome_sequence.count(gRNA)
        ranking[gRNA] = 0

        if gRNA[len(gRNA)-1] != 'G':
            gRNA_list[n_gRNA] += 'G'
            gRNA = gRNA_list[n_gRNA]
        
        if gRNA.find("TTTT") == -1:
            ranking[gRNA] = -1
            continue
        
        GC_content = (gRNA.count('G') + gRNA.count('C'))/len(gRNA)
        if GC_content < 0.3 or GC_content >= 0.8:
            ranking[gRNA] += -30
        else:
            ranking[gRNA] += 30

        if gRNA[-4] == 'A' or gRNA[-4] == 'T' and Sp_cas9 == True:
            ranking[gRNA] += 20
        
        if gRNA.count('TTT') > 1 and gRNA.count('AAA') > 1:
            ranking[gRNA] += -10
        
            

    str1_print = "A 'G' may have been added to the end of some gRNA sequences to improve the efficiency"
    
            

    



def main():
    #sample_seq = tataaatcgtccaatggtacctttaacaggtggtccatccGGGGaaaaaaaatttatatatggttattgtcggctaaggcctacctggactccggta
    gene_sequence_plus = input("Enter the gene sequence from 5' to 3': ").upper()
    gene_sequence_minus = complimentary_sequence(gene_sequence_plus)
    print("The complimentary bases (5' to 3') are: ", gene_sequence_minus)

    PAM = input("Enter the PAM sequence: ").upper()
    print("On the given sequence (gRNA, cut index, PAM):")
    print(find_gRNA(gene_sequence_plus, PAM))
    print("On the complimentary sequence (gRNA, cut index, PAM):")
    print(find_gRNA(gene_sequence_minus, PAM))


if __name__ == "__main__":
    main()