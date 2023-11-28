import functools

from django.contrib.auth.decorators import login_required
from django.core import serializers
from django.http import JsonResponse
from django.shortcuts import render,HttpResponse
import json
import numpy as np
import tensorflow as tf
from django_redis import get_redis_connection

import deepCrispr.deepcrispr as dc
import os
from sgRnaDesigner.beans.beans import RNABean
import pandas as pd
import CFD.CFD.cfdScore as cfd
import re
import datetime
import uuid
import re
import sgRnaDesigner.utils.fastaUtils as faUtils
import  sgRnaDesigner.utils.utils as utils
from django.core.cache import cache
from django_redis import get_redis_connection
import sgRnaDesigner.utils.initSeqScore as initscore
import RNNCrispr.RNNCrispr as RNNPredict
import requests








# Create your views here.
from sgRnaDesigner.beans.beans import Pamdef
import sgRnaDesigner.utils.utils as sgUtils
preWays = []



def index(request):
    return render(request,"index.html")

def welcome(request):
    return render(request,"analysis/sgHome.html")


def home(request):
    return render(request,"sgRnaIndex.html")

def result(request):
    return render(request,"sgRNAResult.html")

def contact(request):
    return render(request,"sgRnaContact.html")

def about(request):
    return render(request,"sgRnaAbout.html")


'''
def analysis(request):
    if request.method == "POST":
        geneSeq = request.POST.get("geneSeq", None)
        print(geneSeq)


    return render(request,"modelanalysis.html")
    '''



def design(request):

    pam1 = Pamdef('1','SpCas9','NGG','Native Streptococcus pyogenes Cas9',20)
    pams = []
    pams.append(pam1)
    #pam2 = Pamdef('VRER SpCas9*', 'NGCG', 'D1135V.G1218R.R1335E.T1337R', 20)
    return render(request, "sgRNADesign.html",{"pams": pams})





def sgRNAFigure(request):
    # 获取request中数据
    userToken = request.GET.get("userToken", None)
    dataIndex = request.GET.get("dataIndex", None)
    genName = request.GET.get("genName", None)
    return render(request, "sgRNAFigure.html", {"userToken":userToken,"dataIndex":dataIndex,"genName":genName})


def sgRNAFigureHeatMap(request):
    return render(request, "sgRNAFigureHeatMap.html")


def login_required_manual(func):
    def wrapper(*args, **kwargs):
        request = args[1]
        if not request.session.get("username"):
            return JsonResponse({"code": -1, "msg": "not login"}, status=401)
        return func(*args, **kwargs)
    return wrapper






def sgRNAFigureCircle(request):

    return render(request, "sgRNAFigureCircle.html")

def sgRNAFigureLife(request):

    return render(request, "sgRNAFigureLife.html")

def sgRNASun(request):

    return render(request, "sgRNAsun.html")






def liftJson(request):
    file = open('/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/sgRnaDesigner/static/json/life-expectancy.json', 'r')
    json_data = file.read()
    print(json_data)
    file.close()
    return HttpResponse(json_data, content_type="application/json", status=200)



def circleJson(request):
    file = open('/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/sgRnaDesigner/static/json/les-miserables.json', 'r')
    json_data = file.read()
    print(json_data)
    file.close()
    return HttpResponse(json_data, content_type="application/json", status=200)





def analysis(request):

    #获取页面数据

    if request.method == "POST":

        pamtype = request.POST.get("pamtype", None)
        seqLen = int(request.POST.get("seqLen", "0"))
        pamCode = request.POST.get("pamCode", None)
        seqtype = request.POST.get("seqtype", None)
        preWayTemp = request.POST.get("preWay", None)
        preWays = preWayTemp.split(",")
        genName = request.POST.get("genName", None)
        # 获取基因名字 针对不同的基因进行预测权重分配 获取权重分配比
        preWays = utils.sortByGene(genName, preWays)

        userEmail = request.POST.get("userEmail", None)

       # print("current predict ways has：%d",len(preWays))

        dataList = []
        sgRNAS = []
        sgRNAS_RNN = []
        sgRNAS_Transformer = []

        if seqtype == '1':
            genSeq = request.POST.get("genSeq", None).replace('\n', '').replace('\r', '')
            genSeq = genSeq.strip()
        if seqtype == '2':

            fastaFile = request.FILES.get("fastaFile")
            # 保存文件到服务器
            now = datetime.datetime.now()
            uname = '/Users/FryTsui/anna/study/fasta/' + str(uuid.uuid1()) + fastaFile.name
            with open(uname, 'wb+') as destination:
                for chunk in fastaFile.chunks():
                    destination.write(chunk)
            # 读取fasta文件中的基因序列
            destination.close()

            fa_dict = faUtils.fasta2dict(uname)
            for key, value in fa_dict.items():
                genSeq = value

        # 得到反向互补序列
        fxgenSeq = faUtils.DNA_rever_complement(genSeq)
        # 获取所有GG序列
        zxLocs = utils.getGGAll(genSeq)
        nxLocs = utils.getGGAll(fxgenSeq)
        # 找到含有GG序列的字符串 向前移动21位 同时选择GC含量在40~60%之间
        for index in range(len(zxLocs)):
            if zxLocs[index] > 21:
                # print("==========="+str(zxLocs[index]))
                # print("===========" + str(genSeq[zxLocs[index] - seqLen - 1:zxLocs[index]]))
                bakSeq = genSeq[zxLocs[index] - seqLen - 1:zxLocs[index]] + pamCode[1:]
                bakSeqGC = sgUtils.accGcper(bakSeq)
                if bakSeqGC > 40 and bakSeqGC < 60:
                    indexEnd = int(index) + int(seqLen)
                    # CNN序列格式封装
                    sgRNAS.append(["", float(index), float(indexEnd), "+", bakSeq, float(0)])
                    # RNN序列组添加
                    sgRNAS_RNN.append(["", float(index), float(indexEnd), "+", bakSeq, "AAAAAAAAAAAAAAAAAAAAAAA",
                                       "NNNNNNNNNNNNNNNNNNNNNNN", "AAAAAAAAAAAAAAAAAAAAAAA", "NNNNNNNNNNNNNNNNNNNNNNN",
                                       float(0)])
                    #Transformer格式封装


        for index in range(len(nxLocs)):
            if nxLocs[index] > 21:
                bakSeq = fxgenSeq[nxLocs[index] - seqLen - 1:nxLocs[index]] + pamCode[1:]
                bakSeqGC = sgUtils.accGcper(bakSeq)
                if bakSeqGC > 40 and bakSeqGC < 60:
                    indexEnd = int(index) + int(seqLen)
                    # CNN序列格式封装
                    sgRNAS.append(["", float(index), float(indexEnd), "-", bakSeq, float(0)])
                    # RNN序列组添加
                    sgRNAS_RNN.append(["", float(index), float(indexEnd), "-", bakSeq, "AAAAAAAAAAAAAAAAAAAAAAA",
                                       "NNNNNNNNNNNNNNNNNNNNNNN", "AAAAAAAAAAAAAAAAAAAAAAA", "NNNNNNNNNNNNNNNNNNNNNNN",
                                       float(0)])
                    # Transformer格式封装

        #根据excellent先级最高的序列获得前20个序列 再根据后续预测方式excellent化排名
        firstPreWay = preWays[0]
        preWays.remove(firstPreWay)

        if firstPreWay == "CNN":
            #excellent先得到CNN的预测前20序列
            cnnSortRNAs =  CNN_Predict(sgRNAS)
            for listzis in range(0, 20):
                #初步封装结果
                rnaBean = {}
                mismatchSeq = []
                rnaBean["genName"] = ""
                rnaBean["startLoc"] = str(cnnSortRNAs[listzis][1])
                rnaBean["endLoc"] = str(cnnSortRNAs[listzis][2])
                rnaBean["RNAValue"] = str(cnnSortRNAs[listzis][4])
                rnaBean["gcPerc"] = str(sgUtils.accGcper(cnnSortRNAs[listzis][4]))
                rnaBean["CFDscore"] = "0"
                rnaBean["cnnScore"] = str(cnnSortRNAs[listzis][5])
                rnaBean["rnnScore"] = "0"
                rnaBean["transScore"] = ''
                rnaBean["strand"] = str(cnnSortRNAs[listzis][3])
                rnaBean["onemismatch"] = int(0)
                rnaBean["twomismatch"] = int(0)
                rnaBean["threemismatch"] = int(0)
                rnaBean["fourmismatch"] = int(0)
                rnaBean["misMatch"] = ""
                rnaBean["misRNA"] = mismatchSeq
                rnaBean["averageScore"] = "0"
                dataList.append(rnaBean)

            #重新封装RNN 预测格式
            sgRNAS_RNN = []
            sgRNAS_Transformer = []
            for listzis in range(0, 20):
                # 重新封装RNN Transformer 预测格式
                sgRNAS_RNN.append(
                    ["", cnnSortRNAs[listzis][1], cnnSortRNAs[listzis][2], cnnSortRNAs[listzis][3], cnnSortRNAs[listzis][4],
                     "AAAAAAAAAAAAAAAAAAAAAAA",
                     "NNNNNNNNNNNNNNNNNNNNNNN", "AAAAAAAAAAAAAAAAAAAAAAA",
                     "NNNNNNNNNNNNNNNNNNNNNNN",
                     float(0)])
                #重新封装Transformer 预测格式
                #需要补充
        if firstPreWay == "RNN":
            rnnSortRNAs = RNN_Predict(sgRNAS_RNN)
            #初步封装结果
            for listzis in range(0, 20):
                #初步封装结果
                rnaBean = {}
                mismatchSeq = []
                rnaBean["genName"] = ""
                rnaBean["startLoc"] = str(rnnSortRNAs[listzis][1])
                rnaBean["endLoc"] = str(rnnSortRNAs[listzis][2])
                rnaBean["RNAValue"] = str(rnnSortRNAs[listzis][4])
                rnaBean["gcPerc"] = str(sgUtils.accGcper(rnnSortRNAs[listzis][4]))
                rnaBean["CFDscore"] = "0"
                rnaBean["cnnScore"] = "0"
                rnaBean["rnnScore"] = str(rnnSortRNAs[listzis][9])
                rnaBean["transScore"] = ''
                rnaBean["offTargetScore"] = "0"
                rnaBean["strand"] = str(cnnSortRNAs[listzis][3])
                rnaBean["misMatch"] = ""
                rnaBean["onemismatch"] = int(0)
                rnaBean["twomismatch"] = int(0)
                rnaBean["threemismatch"] = int(0)
                rnaBean["fourmismatch"] = int(0)
                rnaBean["misRNA"] = mismatchSeq
                #misRNA 为对象 包括4个属性
                rnaBean["averageScore"] = "0"
                dataList.append(rnaBean)

            #重新封装CNN Transfomer的预测格式
            sgRNAS = []
            sgRNAS_Transformer = []
            # 按顺序重新封装CNN
            for listzis in range(0, 20):
                sgRNAS.append(
                    ["", rnnSortRNAs[listzis][1], rnnSortRNAs[listzis][2], rnnSortRNAs[listzis][3],
                     rnnSortRNAs[listzis][4],
                     float(0)]
                )
                # 重新封装Transformer 预测格式
                # 需要补充


        if firstPreWay == "Transformer":
            #重新封装CNN序列
            sgRNAS = []
            #重新封装RNN序列
            sgRNAS_Transformer = []

        #根据后续预测方式excellent化序列排名
        #进行后续预测
        cnnNextPreRNAS = []
        rnnNextPreRNAS = []
        transNextPreRNAS = []
        for preWay in preWays:
            if preWay == "CNN":
                cnnNextPreRNAS = CNN_Predict(sgRNAS,False)
                for i in range(len(dataList)):
                    dataList[i]["cnnScore"] = str(cnnNextPreRNAS[i][5])
            if preWay == "RNN":
                rnnNextPreRNAS = RNN_Predict(sgRNAS_RNN,False)
                for i in range(len(dataList)):
                    dataList[i]["rnnScore"] = str(rnnNextPreRNAS[i][9])
            if preWay == "Transformer":
                transNextPreRNAS = Transformer_Predict(sgRNAS_Transformer)
                for i in range(len(dataList)):
                    dataList[i]["rnnScore"] = str(transNextPreRNAS[i][9])

        #预测脱靶位置 首先获得不匹配字符位置 对位置采用CNN脱靶方式预测是否会脱靶，如果会再进行计算CFD_Score ？？？本次疑问 是否excellent先进行CFD计算
        searchOffTargetByGeneMap(genName, dataList)

        #计算CFD得分
        cal_CFDScore(dataList)

        calOffTargetScore(dataList,len(pamCode)+seqLen )


        #进行最终excellent化排序
        #dataList.sort(key=functools.cmp_to_key(compare_prescore))
        #将结果在服务器上进行保存 生成用户的token 方便后续使用
        userToken = str(uuid.uuid1())

        #将数据保存近服务器
        userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/"+userToken+".txt"
        json_list = json.dumps(dataList)
        with open(userTokenFileName, 'w') as destination:
            destination.write(json_list + "\n")
        destination.close()

    return render(request, "sgRNAResult.html", {"data":dataList,"userToken": userToken,"genName":genName})


def compare_prescore(x,y):
    x_cnn_score = x["cnnScore"]
    y_cnn_score = y["cnnScore"]

    x_rnn_score = x["rnnScore"]
    y_rnn_score = y["rnnScore"]

    x_trans_score = x["transScore"]
    y_trans_score = y["transScore"]

    x_off_score = x["offTargetScore"]
    y_off_score = y["offTargetScore"]

    for preway in preWays:
        if preway == "CNN":
            return x_cnn_score - y_cnn_score;
        if preway == "RNN":
            return x_rnn_score - y_rnn_score
        if preway == "Transformer":
            return x_trans_score - y_trans_score




def calOffTargetScore(dataList,seqLen):
    sess = tf.InteractiveSession()
    off_target_model_dir = '/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/deepCrispr/trained_models/offtar_pt_cnn_reg'
    dcmodel = dc.DCModelOfftar(sess, off_target_model_dir, is_reg=True)
    for data in dataList:
        deepOffScore = float(0)
        targetRNAOffScore = float(0)
        #得到不匹配序列
        misRNAs = data["misRNA"]
        RNASeq = data["RNAValue"]

        for targetRNAs in misRNAs:
            offRNAs = []
            #得到当前misRNA的序列之
            targetRNA = targetRNAs["misSeq"]
            if targetRNA.find("-") == -1 and len(targetRNA) == seqLen:

                offRNA = ['', RNASeq, 'AAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAA',
                          'NNNNNNNNNNNNNNNNNNNNNNN', targetRNA, 'AAAAAAAAAAAAAAAAAAAAAAA']
                offRNA.append('AAAAAAAAAAAAAAAAAAAAAAA')
                offRNA.append('AAAAAAAAAAAAAAAAAAAAAAA')
                offRNA.append('NNNNNNNNNNNNNNNNNNNNNNN')
                offRNA.append(0)
                offRNAs.append(offRNA)
                tf.reset_default_graph()

                #得到在每一个潜在位置的脱靶分数
                input_data = dc.Epiotrt(fpath="", sgRNAs=offRNAs, num_epi_features=4, with_y=True)
                (x_on, x_off), y = input_data.get_dataset()
                x_on = np.expand_dims(x_on, axis=2)  # shape(x) = [100, 8, 1, 23]
                x_off = np.expand_dims(x_off, axis=2)  # shape(x) = [100, 8, 1, 23]

                targetRNAOffScore = dcmodel.offtar_predict(x_on, x_off)
                offTargetScoreList = targetRNAOffScore.tolist()
                offScore = offTargetScoreList[0]
                targetRNAs["misOfftargetScore"] = str(offScore)
                deepOffScore += offScore

        data["offTargetScore"] = str(deepOffScore)
    sess.close()


def cal_CFDScore(dataList):
    for data in dataList:
        data["misMatch"] = "1-"+str(data["onemismatch"])+"-2-"+str(data["twomismatch"])+"-3-"+str(data["threemismatch"])+"-4-"+str(data["fourmismatch"])
        cfd_score = float(0)
        #取data下的misRNA
        misRNAs = data["misRNA"]
        RNASeq = data["RNAValue"]
        for tbRNABak in misRNAs:
            mm_scores, pam_scores = cfd.get_mm_pam_scores()
            wt = RNASeq  # args.wt
            off = tbRNABak["misSeq"]  # args.off
            m_wt = re.search('[^ATCG]', wt)
            m_off = re.search('[^ATCG]', off)
            if (m_wt is None) and (m_off is None):
                pam = off[-2:]
                sg = off[:-3]
                print("wt is "+wt)
                print("sg is " + sg)
                print("pam is " + pam)
                cfd_score += cfd.calc_cfd(wt, sg, pam)
                print("CFD score: %s", str(cfd_score))
        data["CFDscore"] = str(cfd_score)



def compare_score(x,y):
    x_score = x[5]
    y_score = y[5]
    return  y_score -x_score

def compare_RNNscore(x,y):
    x_score = x[9]
    y_score = y[9]
    return   y_score-x_score

def CNN_OffPredict(sgRNA,allOffRNA,doSort=True):
    # deepcrispr 脱靶预测得分
    deepOffScore = float(0)
    offRNAs = []
    if len(allOffRNA) > 0:
        for targetRNAs in allOffRNA:
            if len(targetRNAs) > 0:
             for  targetRNA in targetRNAs:
                offRNA = ['', sgRNA, 'AAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAA',
                          'NNNNNNNNNNNNNNNNNNNNNNN', targetRNA, 'AAAAAAAAAAAAAAAAAAAAAAA']
                offRNA.append('AAAAAAAAAAAAAAAAAAAAAAA')
                offRNA.append('AAAAAAAAAAAAAAAAAAAAAAA')
                offRNA.append('NNNNNNNNNNNNNNNNNNNNNNN')
                offRNA.append(0)
                offRNAs.append(offRNA)
                tf.reset_default_graph()
                input_data = dc.Epiotrt(fpath="", sgRNAs=offRNAs, num_epi_features=4, with_y=True)
                (x_on, x_off), y = input_data.get_dataset()
                x_on = np.expand_dims(x_on, axis=2)  # shape(x) = [100, 8, 1, 23]
                x_off = np.expand_dims(x_off, axis=2)  # shape(x) = [100, 8, 1, 23]
                sess = tf.InteractiveSession()
                off_target_model_dir = '/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/deepCrispr/trained_models/offtar_pt_cnn_reg'
                dcmodel = dc.DCModelOfftar(sess, off_target_model_dir, is_reg=True)
                deepOffScore += dcmodel.offtar_predict(x_on, x_off)
    return deepOffScore

def CNN_Predict(cnn_sgRNAS,doSort=True):
    # 打靶预测得分
    # 使用CNN进行预测
    input_data = dc.Sgt(sgRNAS=cnn_sgRNAS, fpath="", with_y=True)
    x, y = input_data.get_dataset()
    x = np.expand_dims(x, axis=2)  # shape(x) = [10, 4, 1, 23]
    sess = tf.InteractiveSession()
    on_target_model_dir = '/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/deepCrispr/trained_models/ontar_cnn_reg_seq'
    dcmodel = dc.DCModelOntar(sess, on_target_model_dir, is_reg=True, seq_feature_only=True)
    deepOnPre = dcmodel.ontar_predict(x)
    #
    print("predicted_on_CNN_target is %s:", deepOnPre)
    sgRNASON = []

    for sgRNAIndex in range(0, len(cnn_sgRNAS)):
        cnn_sgRNAS[sgRNAIndex][5] = deepOnPre[sgRNAIndex]
    # 按打靶得分排序取前20位
    if doSort:
            # 排序取前20位
            cnn_sgRNAS.sort(key=functools.cmp_to_key(compare_score))
            for listzis in range(0, 20):
                sgRNASON.append(cnn_sgRNAS[listzis])
    else:
        sgRNASON = cnn_sgRNAS
    return sgRNASON

def RNN_Predict(rnn_sgRNAS,doSort=True):
    RNNPredicts = RNNPredict.RNNPredict(rnn_sgRNAS)
    print("RNN Predicts  is %s:", RNNPredicts)
    for sgRNAIndex in range(0, len(rnn_sgRNAS)):
        rnn_sgRNAS[sgRNAIndex][9] = RNNPredicts[sgRNAIndex][0]
    sgRNASRNN = []
    # 排序取前20位
    if doSort:
        rnn_sgRNAS.sort(key=functools.cmp_to_key(compare_RNNscore))
        # 取预测前20个序列
        for listzis in range(0, 20):
            sgRNASRNN.append(rnn_sgRNAS[listzis])
    else:
        sgRNASRNN = rnn_sgRNAS
    return sgRNASRNN


def Transformer_Predict(tranformer_sgRNAS):
    url = "http://127.0.0.1:8000/preTransCrispr/"
    # Transformer 输入参数
    tranformersortRNAS = requests.post(url, data=tranformer_sgRNAS)
    return tranformersortRNAS


#寻找潜在的脱靶位置 利用cnn预测是否会脱靶 并返回CFD_Score
def searchOffTargetByGeneMap(genaName,dataList):


    #生成1个不匹配
    #本地基因组存储位置
    seneLoc = "/Users/FryTsui/anna/study/cas-finder/build/" + genaName
    firstLine = 'NNNNNNNNNNNNNNNNNNNNNRG 2 1'


    firstInputTxt ="/Users/FryTsui/anna/study/cas-finder/mismatch/input"+str(uuid.uuid1()) +".txt"
    firstOnputTxt = "/Users/FryTsui/anna/study/cas-finder/mismatch/output" + str(uuid.uuid1()) + ".txt"

    '''

    twoInputTxt = '/Users/FryTsui/anna/study/cas-finder/mismatch2/' + str(uuid.uuid1()) + ".txt"
    twoOnputTxt = '/Users/FryTsui/anna/study/cas-finder/mismatch2/' + str(uuid.uuid1()) + ".txt"

    threeInputTxt = '/Users/FryTsui/anna/study/cas-finder/mismatch3/' + str(uuid.uuid1()) + ".txt"
    threeOnputTxt = '/Users/FryTsui/anna/study/cas-finder/mismatch3/' + str(uuid.uuid1()) + ".txt"

    fourInputTxt = '/Users/FryTsui/anna/study/cas-finder/mismatch4/' + str(uuid.uuid1()) + ".txt"
    fourOnputTxt = '/Users/FryTsui/anna/study/cas-finder/mismatch4/' + str(uuid.uuid1()) + ".txt"
    '''

    with open(firstInputTxt, 'w') as destination:
        destination.writelines(seneLoc+"\n")
        destination.writelines(firstLine + "\n")
        for i in range(len(dataList)):
            # 寻找潜在不配序列
            sgRNASON = dataList[i]["RNAValue"]

            destination.write(sgRNASON + ' 1' + ' Seq' + str(i) +"\n")
            destination.write(sgRNASON + ' 2' + ' Seq' + str(i) + "\n")
            destination.write(sgRNASON + ' 3' + ' Seq' + str(i) + "\n")
            destination.write(sgRNASON + ' 4' + ' Seq' + str(i) + "\n")
    destination.close()




    #调用cas-offinder指令 测试时暂时屏蔽
    #os.system("/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/Cas-offinder/bin/cas-offinder "+firstInputTxt+" C  "+firstOnputTxt)
    firstOnputTxt = "/Users/FryTsui/anna/study/cas-finder/mismatch1/output.txt"
    '''
    os.system(
        "/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/Cas-offinder/bin/cas-offinder " + twoInputTxt + " C  " + twoOnputTxt)
    os.system(
        "/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/Cas-offinder/bin/cas-offinder " + threeInputTxt + " C  " + threeOnputTxt)
    os.system(
        "/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/Cas-offinder/bin/cas-offinder " + fourInputTxt + " C  " + fourOnputTxt)
'''
    #遍历输出的错配文件
    with open(firstOnputTxt,'r') as destination:
        for mismatchLine in destination:
            if not mismatchLine.startswith("#"):
                mismatchArr = mismatchLine.split()#默认以空格为分割符
                seqTemp = mismatchArr[0] #代表是第几个
                misRNA = {}
                #对应datalist的下标
                misNum = int(seqTemp[3:])
                misMatchNum = int(mismatchArr[7])
                misSeq = mismatchArr[3]
                misChrome = mismatchArr[4]
                misLocation = mismatchArr[5]
                misDirection  = mismatchArr[6]
                misRNA["misSeq"] = misSeq
                misRNA["misChrome"] = misChrome
                misRNA["misLocation"] = misLocation
                misRNA["misDirection"] = misDirection
                misRNA["misOfftargetScore"] = "0"
                dataList[misNum]["misRNA"].append(misRNA)

                if misMatchNum == 1:
                    dataList[misNum]["onemismatch"] = int(dataList[misNum]["onemismatch"]) + 1
                if misMatchNum == 2:
                    dataList[misNum]["twomismatch"] = int(dataList[misNum]["twomismatch"]) + 1
                if misMatchNum == 3:
                    dataList[misNum]["threemismatch"] = int(dataList[misNum]["threemismatch"]) + 1
                if misMatchNum == 4:
                    dataList[misNum]["fourmismatch"] = int(dataList[misNum]["fourmismatch"]) + 1

    destination.close()






















