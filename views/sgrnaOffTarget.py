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

from django.shortcuts import render


def design(request):


    return render(request, "design/sgRNAOfftargetMap.html",{})


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
                print("targetRNA is :"+targetRNA)
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

def design(request):

    #获取参数
    if request.method == "POST":

        pamtype = request.POST.get("pamtype", None)
        pamCode = request.POST.get("pamCode", None)

        genName = request.POST.get("genName", None)

        genSeq = request.POST.get("genSeq", None)

        misMatchnumber = request.POST.get("misMatchnumber", None)

        dnabulgeSize = request.POST.get("dnabulgeSize", None)

        rnabulgeSize = request.POST.get("rnabulgeSize", None)

        genSeqs = genSeq.split()

        #拼接cas-offider参数
        #生成用户token
        userToken = str(uuid.uuid1())
        inputFile = "/Users/FryTsui/anna/study/cas-finder/mismatch/input" + userToken+ ".txt"

        seneLoc = "/Users/FryTsui/anna/study/cas-finder/build/" + genName
        firstLine = "NNNNNNNNNNNNNNNNNNNN"+pamCode+"  " +dnabulgeSize+"  "+rnabulgeSize

        with open(inputFile, 'w') as destination:
            destination.writelines(seneLoc + "\n")
            destination.writelines(firstLine + "\n")
            for i in range(len(genSeqs)):
                # 寻找潜在不配序列
                sgRNASON = genSeqs[i]
                destination.write(sgRNASON + " "+str(misMatchnumber)+ " Seq" + str(i) +"\n")
        destination.close()

        outputFile  = "/Users/FryTsui/anna/study/cas-finder/mismatch/output" + userToken + ".txt"

        #执行cas-finder命令
        os.system(
            "/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/Cas-offinder/bin/cas-offinder " + inputFile + " C  " + outputFile)

        #分析输出文件 将misMatch数据进行存储
        #分别将不同序列的错配值存储进datalist
        dataList = []
        for i in range(len(genSeqs)):
            misRNA = []
            dataList.append(misRNA)
        sess = tf.InteractiveSession()
        off_target_model_dir = '/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/deepCrispr/trained_models/offtar_pt_cnn_reg'
        dcmodel = dc.DCModelOfftar(sess, off_target_model_dir, is_reg=True)

        with open(outputFile, 'r') as destination:
            for mismatchLine in destination:
                if not mismatchLine.startswith("#"):
                    mismatchArr = mismatchLine.split()  # 默认以空格为分割符
                    seqTemp = mismatchArr[0]  # 代表是第几个
                    misRNA = {}
                    # 对应datalist的下标
                    misNum = int(seqTemp[3:])
                    misMatchNum = int(mismatchArr[7])
                    misBuType = mismatchArr[1]
                    misSeq = mismatchArr[3]
                    misChrome = mismatchArr[4]
                    misLocation = mismatchArr[5]
                    misDirection = mismatchArr[6]
                    misRNA["misBuType"] = misBuType
                    misRNA["misSeq"] = misSeq
                    misRNA["misChrome"] = misChrome
                    misRNA["misLocation"] = misLocation
                    misRNA["misDirection"] = misDirection
                    misRNA["misOfftargetScore"] = "0"
                    #预测脱靶得分
                    offRNAs = []
                    if misSeq.find("-") == -1 and len(misSeq) ==len(genSeqs[misNum]) :
                        offRNA = ['',genSeqs[misNum] , 'AAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAA',
                                  'AAAAAAAAAAAAAAAAAAAAAAA',
                                  'NNNNNNNNNNNNNNNNNNNNNNN', misSeq, 'AAAAAAAAAAAAAAAAAAAAAAA']
                        offRNA.append('AAAAAAAAAAAAAAAAAAAAAAA')
                        offRNA.append('AAAAAAAAAAAAAAAAAAAAAAA')
                        offRNA.append('NNNNNNNNNNNNNNNNNNNNNNN')
                        offRNA.append(0)
                        offRNAs.append(offRNA)
                        tf.reset_default_graph()
                        # 得到在每一个潜在位置的脱靶分数
                        input_data = dc.Epiotrt(fpath="", sgRNAs=offRNAs, num_epi_features=4, with_y=True)
                        (x_on, x_off), y = input_data.get_dataset()
                        x_on = np.expand_dims(x_on, axis=2)  # shape(x) = [100, 8, 1, 23]
                        x_off = np.expand_dims(x_off, axis=2)  # shape(x) = [100, 8, 1, 23]

                        targetRNAOffScore = dcmodel.offtar_predict(x_on, x_off)
                        offTargetScoreList = targetRNAOffScore.tolist()
                        offScore = offTargetScoreList[0]
                        misRNA["misOfftargetScore"] = str(offScore)
                    #计算CFDScore

                    dataList[misNum].append(misRNA)

        sess.close()
        destination.close()


    return render(request, "design/sgRNAOfftargetResult.html",{"datalist":dataList,"userToken":userToken,"genName":genName})
