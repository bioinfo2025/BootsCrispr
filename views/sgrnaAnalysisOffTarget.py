import functools
import random
import smtplib
from email.header import Header
from email.mime.text import MIMEText

from django.contrib.auth.decorators import login_required
from django.core import serializers
from django.http import JsonResponse
from django.shortcuts import render,HttpResponse
import json
import numpy as np
import tensorflow as tf
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
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
from django.core.paginator import Paginator, Page, EmptyPage, PageNotAnInteger


from django.shortcuts import render


def design(request):


    return render(request, "analysis/sgRNAOfftargetMap.html",{})


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

def analysis(request):

    #获取参数
    if request.method == "POST":

        pamtype = request.POST.get("pamtype", None)
        pamCode = request.POST.get("pamCode", None)

        genName = request.POST.get("genName", None)

        genSeq = request.POST.get("genSeq", None)

        misMatchnumber = request.POST.get("misMatchnumber", None)

        dnabulgeSize = request.POST.get("dnabulgeSize", None)

        rnabulgeSize = request.POST.get("rnabulgeSize", None)

        userEmail = request.POST.get("userEmail", None)

        userToken = str(uuid.uuid1())

        genSeqs = genSeq.split()
        genList = []
        for genSeq in genSeqs:
            genList.append(genSeq)

            # 发生token邮件 后续可根据token获取内容
            smtp_server = initscore.config_dict.get("smtp_server")
            username = initscore.config_dict.get("username")
            password = initscore.config_dict.get("password")
            smtp_port = initscore.config_dict.get("smtp_port")

            # 配置发送模版
            content = 'This is a letter from ' + initscore.config_dict.get(
                "sys_name") + " ,please save your token \" "+ userToken + "\".We will save your results for 2 months." \
                                                                        "During this period, you can obtain your experimental results based on your token through \"Contact\" section！！！"
            subject = 'Your token from ' + initscore.config_dict.get("sys_name")

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



        #拼接cas-offider参数
        #生成用户token
        basepath = initscore.config_dict.get("basepath")

        #根据userToken
        dirPath = basepath + userToken
        #判断目录是否存在

        if not os.path.exists(dirPath):
            os.mkdir(dirPath)

        #将cas-finder的参数输入文件进行创建
        inputFile = dirPath + "/cas-input.txt"

        casfinderPath = initscore.config_dict.get("casfinder")

        seneLoc = casfinderPath +"/" +genName
        firstLine = "NNNNNNNNNNNNNNNNNNNN"+pamCode+"  " +dnabulgeSize+"  "+rnabulgeSize

        with open(inputFile, 'w') as destination:
            destination.writelines(seneLoc + "\n")
            destination.writelines(firstLine + "\n")
            for i in range(len(genSeqs)):
                # 寻找潜在不配序列
                sgRNASON = genSeqs[i]
                destination.write(sgRNASON + " "+str(misMatchnumber)+ " Seq" + str(i) +"\n")
        destination.close()

        outputFile  =  dirPath + "/cas-output.txt"

        #执行cas-finder命令
        #


        os.system("cd "+casfinderPath)
        outputFile = dirPath + "/cas-output.txt"
        file = open(outputFile, "w")
        file.close()
        #os.system("cas-offinder " + inputFile + " C  " + outputFile)
        outputFile = "/Users/FryTsui/anna/study/cas-finder/build/output.txt"


        #分析输出文件 将misMatch数据进行存储
        #分别将不同序列的错配值存储进datalist
        dataList = []
        sess = tf.InteractiveSession()
        #off_target 文件路径


        off_target_model_dir =os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) +"/deepCrispr/trained_models/offtar_pt_cnn_reg"
        #off_target_model_dir = '/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/deepCrispr/trained_models/offtar_pt_cnn_reg'
        dcmodel = dc.DCModelOfftar(sess, off_target_model_dir, is_reg=True)
        try:
            for i in range(len(genSeqs)):
                misRNA = []
                dataList.append(misRNA)

            offtargetNum = 0
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
                        misRNA["seq"] = mismatchArr[2]
                        misRNA["misBuType"] = misBuType
                        misRNA["misSeq"] = misSeq
                        misRNA["misChrome"] = misChrome
                        misRNA["misLocation"] = misLocation
                        misRNA["misDirection"] = misDirection
                        #记得恢复
                        # misRNA["misOfftargetScore"] = "0"
                        misRNA["misOfftargetScore"] = str(random.randint(10000, 999999)/1000000)
                        misRNA["misCFDScore"] = str(random.randint(10000, 999999)/1000000)
                        #记得恢复
                        #预测脱靶得分
                        offRNAs = []
                        ll = len(genSeqs[misNum-1]) - len(pamCode)+1
                        tarSeq = genSeqs[misNum-1][:ll]+pamCode[1:]


                        if misSeq.find("-") == -1 and len(misSeq) == len(tarSeq) and tarSeq.find("N") == -1  :
                            offRNA = ['',tarSeq , 'AAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAA',
                                      'AAAAAAAAAAAAAAAAAAAAAAA',
                                      'NNNNNNNNNNNNNNNNNNNNNNN', misSeq.upper(), 'AAAAAAAAAAAAAAAAAAAAAAA']
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
                            if offScore > 0:
                                offtargetNum = offtargetNum  + 1
                            misRNA["misOfftargetScore"] = str(offScore)
                        #计算CFDScore

                        dataList[misNum-1].append(misRNA)
        except Exception as e:
            print(str(e))
            sess.close()
            destination.close()
        finally:
            sess.close()
            destination.close()

      #将datalist写入results文件
        # 将数据保存进usertoken服务器的results

        #userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/" + userToken + ".txt"

        userTokenFileName = dirPath + "/OfftargetResults.txt"
        json_list = json.dumps(dataList)
        with open(userTokenFileName, 'w') as destination:
            destination.write(json_list + "\n")
        destination.close()
    return render(request, "analysis/sgRNAofftartgetScroll.html",{"userToken":userToken,"genName":genName,"seqs":json.dumps(genList)})


def figureoverview(request):

  #获取request中数据
  userToken = request.GET.get("userToken", None)

  genName = request.GET.get("genName", None)

  return  render(request, "analysis/sgFigureOverview.html",{"userToken":userToken,"genName":genName})


def pierose(request):

  #获取request中数据
  userToken = request.GET.get("userToken", None)

  genName = request.GET.get("genName", None)
  chroms = initscore.chromoMap.get(genName)


  #从usertoken中分别读多个数据
  json_str = ""
  # 读取用户数据
  basepath = initscore.config_dict.get("basepath")
  dirPath = basepath + userToken
  userTokenFileName = dirPath + "/OfftargetResults.txt"
  #userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/" + userToken + ".txt"
  # 遍历输出的错配文件
  with open(userTokenFileName, 'r') as destination:
      for line in destination:
          json_str += line
  # 读取文件
  datalist = json.loads(json_str)
  datas = []#多个搜索序列
  data = []#每个搜索序列的list 里面append { value: 40, name: 'chr1' },
  for data in datalist:
      data = []
      for i in range(chroms):
          #data.append({"value":int(0),"name":"chr"+str(i+1)})
          #记得恢复
          data.append({"value": random.randint(100,200), "name": "chr" + str(i + 1)})
      for misrna in data:
          #取misChrome
          misChrome = misrna["misChrome"]

  return  render(request, "analysis/sgpieroseType.html",{"userToken":userToken,"genName":json.dumps()})



def sgRNAofftartgetScroll(request):
    return render(request, "analysis/sgRNAofftartgetScroll.html", {})

def sgRNAOfftargetResult(request):
    dataIndex = request.GET.get("dataIndex", "1")

    userToken = request.GET.get("userToken", None)

    genName = request.GET.get("genName", None)

    #rows（每页多少行）、page（当前显示第几页）

    rows = int(request.GET.get("rows", 20))
    page = int(request.GET.get("page", 1)) -1

    #读取文件
    # 从usertoken中分别读多个数据
    json_str = ""
    # 读取用户数据

    # 读取用户数据
    basepath = initscore.config_dict.get("basepath")
    dirPath = basepath + userToken
    userTokenFileName = dirPath + "/OfftargetResults.txt"

    # 遍历输出的错配文件
    with open(userTokenFileName, 'r') as destination:
        for line in destination:
            json_str += line
    # 读取文件
    datalist = json.loads(json_str)
    pageStart = page*rows
    pageEnd = pageStart +rows

    total = len(datalist[int(dataIndex)-1])

  #获取染色体数据
    chromsDatas = []
    chromsData = {}
    chromeLen = initscore.chromoMap.get(genName)
    for i in range(chromeLen):
        #拼接数据
        #记得恢复
        #chromsData = {"value":0,"name":"chr"+str(i)}
        chromsData = {"value": random.randint(10,300), "name": "chr" + str(i)}
        chromsDatas.append(chromsData)
    #取对应的datalist数据
    rnaDatas = datalist[int(dataIndex) -1]
    for rna in rnaDatas:
        misChrome = rna["misChrome"]
        if misChrome[3:4].isdigit():
            misChrNum = int(misChrome[3:4]) - 1  # 得到是在哪个染色体上的错配 得到后相应的染色体错配数加1
            chromsDatas[misChrNum]["value"] = int(chromsDatas[misChrNum]["value"]) + 1

    #return render(request, "analysis/sgRNAOfftargetResult.html", {"paper":paper,"userToken":userToken,"genName":genName,"pageSize":10,"page":pageNumber,"total":len(datalist[int(dataIndex)-1])})
    return render(request, "analysis/sgRNAOfftargetResult.html",
                  {"piedata":json.dumps(chromsDatas),"total":total,"page":page+1, "userToken":userToken,"genName":genName,"dataIndex":dataIndex})




@require_http_methods(["POST"])
@csrf_exempt
def sgRNAOfftargetResultJson(request):
    dataIndex = request.POST.get("dataIndex", "1")

    userToken = request.POST.get("userToken", None)

    genName = request.POST.get("genName", None)

    #rows（每页多少行）、page（当前显示第几页）

    rows = int(request.POST.get("rows", 20))
    page = int(request.POST.get("page", 1)) -1

    #读取文件
    # 从usertoken中分别读多个数据
    json_str = ""
    # 读取用户数据
   # userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/" + userToken + ".txt"
    basepath = initscore.config_dict.get("basepath")
    dirPath = basepath + userToken
    userTokenFileName = dirPath + "/OfftargetResults.txt"
    # 遍历输出的错配文件
    with open(userTokenFileName, 'r') as destination:
        for line in destination:
            json_str += line
    # 读取文件
    datalist = json.loads(json_str)
    pageStart = page*rows
    pageEnd = pageStart +rows
    datalistIndex = int(dataIndex)-1
    datas = datalist[datalistIndex][pageStart:pageEnd]
    total = len(datalist[datalistIndex])

    #return render(request, "analysis/sgRNAOfftargetResult.html", {"paper":paper,"userToken":userToken,"genName":genName,"pageSize":10,"page":pageNumber,"total":len(datalist[int(dataIndex)-1])})
    return JsonResponse({"data":json.dumps(datas)}, status=200, content_type='application/json; charset=utf-8')





