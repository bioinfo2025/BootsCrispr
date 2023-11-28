from decimal import Decimal
import pandas as pd
import numpy as np
import requests
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def accGcper(rnaSrq):
    gCount = rnaSrq.count("G")
    cCount = rnaSrq.count("C")
    gcCount = gCount + cCount
    gcperc= gcCount/len(rnaSrq)
    return Decimal(gcperc).quantize(Decimal('0.00')) * 100
#
import re
def getGGAll(genSeq):
    allLoc = []
    for m in re.finditer("GG", genSeq):
        allLoc.append(m.start())
    return allLoc

def sortByGene(geneName,preWays):
    #hg19 CNN > RNN >Transformer
    sorrPreWays = []
    if geneName is "hg19":
     for preway in preWays:
         if "CNN" == preway:
             sorrPreWays.append(preway)
     for preway in preWays:
        if "RNN" == preway:
            sorrPreWays.append(preway)
     for preway in preWays:
        if "Transformer" == preway:
            sorrPreWays.append(preway)
     #针对其他基因序列权重 后续继续添加
    else:
        for preway in preWays:
            if "CNN" == preway:
                sorrPreWays.append(preway)
        for preway in preWays:
            if "RNN" == preway:
                sorrPreWays.append(preway)
        for preway in preWays:
            if "Transformer" == preway:
                sorrPreWays.append(preway)

    return sorrPreWays

def httpPost():
    url = "http://127.0.0.1:8000/transCrispr/predict/"
    # 请求Header(字典形式储存)
    header = {"content-type": "application/json"}
    # 请求Body(字典形式储存)
    datas = []

    datas.append(["CACCCTACAAATCCTCCTCCGGT","0.343251947"])
    datas.append(["CACCCTACAAATCCTCCTCCGGT","0.343251947"])
    post_data = {"inputData": datas}
    # 发送POST请求
    response = requests.post(url, data=post_data, headers=header)










