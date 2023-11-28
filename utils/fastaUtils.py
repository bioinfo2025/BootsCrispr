# 功能：读取解压后的*.fasta文件，将内容按键值对形式存为字典
# 输入：
## fasta_name：解压后的文件名
# 输出：
## fa_dict：包含描述和序列的字典
def fasta2dict(fasta_name):
    fa_dict = {}
    global seq_name
    with open(fasta_name) as fa:
        for line in fa:
            # 去除末尾换行符
            line = line.replace('\n','')
            if line.startswith('>'):
                # 去除 > 号
                seq_name = line[1:] if line[1:] is None or  line[1:]=='' else "00"
                fa_dict[seq_name] = ''
            else:
                # 去除末尾换行符并连接多行序列
                fa_dict[seq_name] += line.replace('\n','')
    return fa_dict


def DNA_rever_complement(sequence):
    trantab = str.maketrans('ACGTacgt', 'TGCAtgca')     # trantab = str.maketrans(intab, outtab)   # 制作翻译表
    string = sequence.translate(trantab)     # str.translate(trantab)  # 转换字符
    string_tmp = list(string)[::-1]
    string_fin = ''.join(string_tmp)

    return string_fin


# 功能： 获取fasta文件反向互补序列、计算fasta文件总碱基个数、N50/N90、GC含量、返回指定位置指定长度的碱基、将所有碱基以300连续“N”作为间隔得到一整条序列
# 输入：
## fasta_dict：值为序列信息的字典
## N_percent: N50时为50, N90时为90, 也可以为0~100内的其他整数
## *cut_local：可变长度参数，指定需要提取的碱基位置信息，数据为偶数且为整数。注：此参数可以缺省，表示不进行”返回指定位置指定长度的碱基“的计算
# 输出：总碱基长度，N50或N90以及任意Nn，GC含量，截取指定长度碱基存入字典，将所有碱基以300连续“N”作为间隔得到一整条序列
## rever_comple_seq：获取fasta文件反向互补序列，输出为字典
## all_dna_len: fasta文件中总碱基数
## N_len：N50或N90长度, 输入参数（即，N_percent）需要对应
## GC_percent：GC含量
## fasta_dict_cut：获得指定位置，指定长度的碱基，输出为字典
## all_dna_fin：将所有碱基以300连续“N”作为间隔得到一整条序列，输出为列表
def return_Nn(fasta_dict, N_percent: int, *cut_local: int):
    import sys

    rever_comple_seq = {}
    fasta_dict_counts = {}
    fasta_dict_cut = {}
    all_dna = []
    all_dna_len = 0
    len_tmp = 0
    nGC = 0

    # 将原来字典中碱基序列换成序列长度
    ## key为序列描述信息，value为碱基信息
    for key, value in fasta_dict.items():
        # GC碱基数
        nGC_tmp = value.count("g") + value.count("G") + value.count("c") + value.count("C")
        nGC += nGC_tmp
        # 在每条序列后面接300个连续“N”
        all_dna_tmp = value + 'N' * 300
        all_dna += all_dna_tmp

        # 调用上述自定义的“DNA_rever_complement”函数， 获取反向互补序列，并存为字典
        rever_comple_seq[key] = DNA_rever_complement(value)
        fasta_dict_counts[key] = len(value)

        # 获取指定长度，指定位置碱基序列
        ## 从第三个参数开始为指定碱基的位置参数（即，cut_local参数），可选参数
        ## 每两个为一对，该参数个数必须为偶数
        cut_dna = []
        # 技巧：以步长为2取cut_local参数信息
        for i in range(len(cut_local))[::2]:
            if i + 1 < len(cut_local):
                cut_dna.append(value[cut_local[i]:cut_local[i + 1]])
            else:
                print("Error:the lengths of cut_local parameters must be even!!!")
                sys.exit()
        fasta_dict_cut[key] = cut_dna

        # *.fasta文件总碱基数
        all_dna_len += len(value)
    # 所有碱基整合为一条序列
    all_dna_fin = ''.join(all_dna)
    # GC含量
    GC_percent = nGC / all_dna_len
    # print("GC_percent:%f" %GC_percent)
    # N50或N90
    cutoff = int(all_dna_len * (N_percent / 100))
    # 将得到的序列长度字典按序列长度降序排序，用于计算N50,N90
    fasta_dict_counts_order = dict(sorted(fasta_dict_counts.items(), key=lambda k: k[1], reverse=True))  # 按碱基长度排序得到新字典
    # 注：lambda为匿名函数
    for key, value in fasta_dict_counts_order.items():
        len_tmp += value
        if len_tmp >= cutoff:
            # print("N%d为：%d" %(N_percent, value))
            N_len = value
            break

    return rever_comple_seq, all_dna_len, N_len, GC_percent, fasta_dict_cut, all_dna_fin



