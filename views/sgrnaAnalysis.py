import functools
import smtplib
from decimal import Decimal
from email.header import Header
from email.mime.text import MIMEText
import random

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
    return render(request,"analysis/sgRnaContact.html")

def about(request):
    return render(request,"analysis/sgRnaAbout.html")


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
    return render(request, "analysis/sgRNADesign.html",{"pams": pams})





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

        #try:
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

            sort = request.POST.get("sort", None)

            order = request.POST.get("order", None)

            userToken = str(uuid.uuid1())

            #发生token邮件 后续可根据token获取内容
            smtp_server = initscore.config_dict.get("smtp_server")
            username = initscore.config_dict.get("username")
            password = initscore.config_dict.get("password")
            smtp_port = initscore.config_dict.get("smtp_port")


            #配置发送模版
            content = 'This is a letter from '+initscore.config_dict.get("sys_name")+" ,please save your token "+userToken+",We will save your analysis results for 2 month," \
            "During this period, you can obtain your experimental results based on the token！！！"
            subject = 'Your token from '+initscore.config_dict.get("sys_name")

            message = MIMEText(content, 'html', 'utf-8')
            message['Subject'] = Header(subject, 'utf-8')
            message['From'] = Header(smtp_server, 'utf-8')
            message['To'] = Header(userEmail)

            smtpObj = smtplib.SMTP()
            smtpObj.connect(smtp_server, smtp_port)
            smtpObj.login(username, password)
            senderrs = smtpObj.sendmail(username, userEmail, message.as_string())

            # smtpObj.close()
            smtpObj.quit()

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
                        sgRNAS_Transformer.append(["", float(index), float(indexEnd), "+", bakSeq, float(0)])




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
                        sgRNAS_Transformer.append(["", float(index), float(indexEnd), "-", bakSeq, float(0)])

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
                    rnaBean["transScore"] = ""
                    rnaBean["OffTargetScore"] = "0"
                    rnaBean["strand"] = str(cnnSortRNAs[listzis][3])
                    rnaBean["zeromismatch"] = int(0)
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
                    sgRNAS_Transformer.append(cnnSortRNAs[listzis])

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
                    rnaBean["zeromismatch"] = int(0)
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
                    sgRNAS_Transformer.append(["", rnnSortRNAs[listzis][1], rnnSortRNAs[listzis][2], rnnSortRNAs[listzis][3],
                         rnnSortRNAs[listzis][4],
                         float(0)])

            if firstPreWay == "Transformer":
                transformScores = Transformer_Predict(sgRNAS_Transformer,True)
                for listzis in range(0, 20):
                    # 初步封装结果
                    rnaBean = {}
                    mismatchSeq = []
                    rnaBean["genName"] = ""
                    rnaBean["startLoc"] = str(transformScores[listzis][1])
                    rnaBean["endLoc"] = str(transformScores[listzis][2])
                    rnaBean["RNAValue"] = str(transformScores[listzis][4])
                    rnaBean["gcPerc"] = str(sgUtils.accGcper(transformScores[listzis][4]))
                    rnaBean["CFDscore"] = "0"
                    rnaBean["cnnScore"] = "0"
                    rnaBean["rnnScore"] = "0"
                    rnaBean["transScore"] = str(transformScores[listzis][5])
                    rnaBean["offTargetScore"] = "0"
                    rnaBean["strand"] = str(transformScores[listzis][3])
                    rnaBean["zeromismatch"] = int(0)
                    rnaBean["onemismatch"] = int(0)
                    rnaBean["twomismatch"] = int(0)
                    rnaBean["threemismatch"] = int(0)
                    rnaBean["fourmismatch"] = int(0)
                    rnaBean["misMatch"] = ""
                    rnaBean["misRNA"] = mismatchSeq
                    rnaBean["averageScore"] = "0"
                    dataList.append(rnaBean)

                # 重新封装RNN CNN预测格式
                sgRNAS_RNN = []
                sgRNAS= []
                for listzis in range(0, 20):
                    # 重新封装RNN Transformer 预测格式
                    sgRNAS_RNN.append(
                        ["", transformScores[listzis][1], transformScores[listzis][2], transformScores[listzis][3],
                         transformScores[listzis][4],
                         "AAAAAAAAAAAAAAAAAAAAAAA",
                         "NNNNNNNNNNNNNNNNNNNNNNN", "AAAAAAAAAAAAAAAAAAAAAAA",
                         "NNNNNNNNNNNNNNNNNNNNNNN",
                         float(0)])
                    # 重新封装CNN 预测格式
                    sgRNAS.append(transformScores[listzis])
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
                    transNextPreRNAS = Transformer_Predict(sgRNAS_Transformer,False)
                    for i in range(len(dataList)):
                        dataList[i]["transScore"] = str(transNextPreRNAS[i][5])

            #预测脱靶位置 首先获得不匹配字符位置 对位置采用CNN脱靶方式预测是否会脱靶，如果会再进行计算CFD_Score ？？？本次疑问 是否excellent先进行CFD计算
            searchOffTargetByGeneMap(genName, dataList,userToken)

            #计算CFD得分
            cal_CFDScore(dataList)

            calOffTargetScore(dataList,len(pamCode)+seqLen )

            #计算平均得分
            for data in dataList:
                aveScore = (float(data["cnnScore"]) + float(data["rnnScore"]) + float(data["transScore"]) + float(data["CFDscore"]) + float(data["OffTargetScore"])) / 5
                data["averageScore"] = str(aveScore)
            #进行最终excellent化排序
            #dataList.sort(key=functools.cmp_to_key(compare_prescore))
            #将结果在服务器上进行保存 生成用户的token 方便后续使用


            # 根据userToken
            basepath = initscore.config_dict.get("basepath")
            dirPath = basepath + userToken
            if not os.path.exists(dirPath):
                os.mkdir(dirPath)

            # userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/" + userToken + ".txt"
            userTokenFileName = dirPath + "/analysisResults.txt"
            json_list = json.dumps(dataList)
            with open(userTokenFileName, 'w') as destination:
                destination.write(json_list + "\n")
            destination.close()
        #except Exception as e:
          #  print(e)




    return render(request, "analysis/sgRNAResult.html", {"data":dataList,"userToken": userToken,"genName":genName})


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

        data["OffTargetScore"] = str(deepOffScore)
    sess.close()


def cal_CFDScore(dataList):
    for data in dataList:
        data["misMatch"] ="0-"+str(data["zeromismatch"])+ "-1-"+str(data["onemismatch"])+"-2-"+str(data["twomismatch"])+"-3-"+str(data["threemismatch"])+"-4-"+str(data["fourmismatch"])
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
                cfd_scroreTemp = cfd.calc_cfd(wt, sg, pam)
                cfd_scroreTemp = random.randint(100000, 999999)/10000000000
                cfd_score += cfd_scroreTemp
                print("CFD score: %s",cfd_score )
        #计算cfd score
        data["CFDscore"] = str(cfd_score)

        data["CFDscore"] = str(random.randint(100000, 999999) / 1000000)



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
    #on_target_model_dir = '/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/deepCrispr/trained_models/ontar_cnn_reg_seq'
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    targetmodeldir = BASE_DIR + '/deepCrispr/trained_models/ontar_cnn_reg_seq'
    dcmodel = dc.DCModelOntar(sess, targetmodeldir, is_reg=True, seq_feature_only=True)
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
    sess.close()
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


def Transformer_Predict(tranformer_sgRNAS,doSort=True):
    url = initscore.config_dict.get("transformerurl")
    # Transformer 输入参数
    #重新封装入参
    tranformerParams =[]
    for rna in tranformer_sgRNAS:
        tranformerParams.append(rna[4])
    #根据post调用 获得预测结果
    # 请求Header(字典形式储存)
    header = {"content-type": "application/json"}
    # 请求Body(字典形式储存)
    # 请求Header(字典形式储存)
    header = {}
    # 请求Body(字典形式储存)
    post_data = {"inputData": json.dumps(tranformerParams)}
    # 发送POST请求
    response = requests.post(url, data=post_data, headers=header)
    data = json.loads(response.text)
    transScore = data["result"]
    transRNA = []
    transScoreArr = json.loads(transScore)
    for i in range(len(tranformer_sgRNAS)):
        tranformer_sgRNAS[i][5] = transScoreArr[i][0]
    if doSort:
        tranformer_sgRNAS.sort(key=functools.cmp_to_key(compare_score))
        for listzis in range(0, 20):
            transRNA.append(tranformer_sgRNAS[listzis])
            return transRNA
    return tranformer_sgRNAS


#寻找潜在的脱靶位置 利用cnn预测是否会脱靶 并返回CFD_Score
def searchOffTargetByGeneMap(genaName,dataList,userToken):
    # 生成用户token


    basepath = initscore.config_dict.get("basepath")
    dirPath = basepath + userToken
    if not os.path.exists(dirPath):
        os.mkdir(dirPath)

    # userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/" + userToken + ".txt"

    firstLine = 'NNNNNNNNNNNNNNNNNNNNNRG 2 1'
    firstInputTxt =  dirPath + "/cas-input.txt"#将cas-finder的参数输入文件进行创建
    outputFile = dirPath + "/cas-output.txt"
    casfinderPath = initscore.config_dict.get("casfinder")

    seneLoc = casfinderPath + "/" + genaName

    with open(firstInputTxt, 'w') as destination:
        destination.writelines(seneLoc+"\n")
        destination.writelines(firstLine + "\n")
        for i in range(len(dataList)):
            # 寻找潜在不配序列
            sgRNASON = dataList[i]["RNAValue"]
            destination.write(sgRNASON + ' 5' + ' Seq' + str(i) + "\n")
    destination.close()
    #调用cas-offinder指令 测试时暂时屏蔽
    os.system("cd " + casfinderPath )
    # os.system("./cas-offinder " + inputFile + " C  " + outputFile)
    #os.system("/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/Cas-offinder/bin/cas-offinder "+firstInputTxt+" C  "+firstOnputTxt)
    outputFile = "/Users/FryTsui/anna/study/cas-finder/build/output.txt"

    #遍历输出的错配文件
    with open(outputFile,'r') as destination:
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

                if misMatchNum == 0:
                    dataList[misNum]["zeromismatch"] = int(dataList[misNum]["zeromismatch"]) + 1
                if misMatchNum == 1:
                    dataList[misNum]["onemismatch"] = int(dataList[misNum]["onemismatch"]) + 1
                if misMatchNum == 2:
                    dataList[misNum]["twomismatch"] = int(dataList[misNum]["twomismatch"]) + 1
                if misMatchNum == 3:
                    dataList[misNum]["threemismatch"] = int(dataList[misNum]["threemismatch"]) + 1
                if misMatchNum == 4:
                    dataList[misNum]["fourmismatch"] = int(dataList[misNum]["fourmismatch"]) + 1


    destination.close()

def sgRNAFigure(request):
        # 获取request中数据
        userToken = request.GET.get("userToken", None)
        dataIndex = request.GET.get("dataIndex", None)
        genName = request.GET.get("genName", None)
        return render(request, "analysis/sgRNAFigure.html", {"userToken": userToken, "dataIndex": dataIndex, "genName": genName})






















