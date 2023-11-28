import functools
import random
from decimal import Decimal

from django.contrib.auth.decorators import login_required
from django.core import serializers
from django.http import JsonResponse
from django.shortcuts import render,HttpResponse

import sgRnaDesigner.utils.initSeqScore as initscore

import json

def figureoverview(request):

  #获取request中数据
  userToken = request.GET.get("userToken", None)

  genName = request.GET.get("genName", None)

  return  render(request, "figures/sgFigureOverview.html",{"userToken":userToken,"genName":genName})


def lineStack(request):

  # var legenddatalist = ['CNN-Score','RNN-Score','Transformer-Score','CFD-Score','Offtarget-Score','Average-Score']
  # 获取request中数据
  userToken = request.GET.get("userToken", None)

  json_str = ""
  # 读取用户数据


  #userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/" + userToken + ".txt"
  # 根据userToken
  basepath = initscore.config_dict.get("basepath")
  dirPath = basepath + userToken

  # userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/" + userToken + ".txt"
  userTokenFileName = dirPath + "/analysisResults.txt"
  # 遍历输出的错配文件
  with open(userTokenFileName, 'r') as destination:
    for line in destination:
      json_str += line
  destination.close()
  # 读取文件
  datalist = json.loads(json_str)

 #拼接CNNscore
  cnnScores = []
  rnnScores = []
  transformScores= []
  cfdScores = []
  offtargetScores = []
  averageScores = []
  scores = []

  for data in datalist:
      cnnScore = "0"
      if data["cnnScore"] != "":
        cnnScore = str(Decimal(100*float(data["cnnScore"])).quantize(Decimal("00")) )
      cnnScores.append(cnnScore)

      rnnScore = "0"
      if data["rnnScore"] != "":
        #rnnScore = str(round(float(data["rnnScore"]),4) * 100)
        rnnScore = str(Decimal(100 * float(data["rnnScore"])).quantize(Decimal("00")))
      rnnScores.append(rnnScore)

      transScore = "0"
      if data["transScore"] != "":
        #transScore = str(round(float(data["transScore"]), 4) * 100)
        transScore = str(Decimal(100 * float(data["transScore"])).quantize(Decimal("00")))
      transformScores.append(transScore)

      CFDscore = "0"
      if data["CFDscore"] != "":
       # CFDscore = str(round(float(data["CFDscore"]), 4) * 100)
       CFDscore = str(Decimal(100 * float(data["CFDscore"])).quantize(Decimal("00")))
      cfdScores.append(CFDscore)

      offTargetScore = "0"
      if data["OffTargetScore"] != "":
        #offTargetScore = str(round(float(data["offTargetScore"]), 4) * 100)
        offTargetScore = str(Decimal(10000 * float(data["OffTargetScore"])).quantize(Decimal("00")))
      offtargetScores.append(offTargetScore)

      aveScore = (float(cnnScore) + float(rnnScore) +float(transScore) + float(CFDscore) +float(offTargetScore))/5
      averageScores.append(str(Decimal(float(aveScore)).quantize(Decimal("00"))))


  scores.append(cnnScores)
  scores.append(rnnScores)
  scores.append(transformScores)
  scores.append(cfdScores)
  scores.append(offtargetScores)
  scores.append(averageScores)

  scores_list = json.dumps(scores)
  return  render(request, "figures/lineStack.html",{"scores":scores_list})


def sgRNASunburst(request):

    return render(request, "figures/sgRNAsunburst.html")


def sgRNAPienet(request):
    # 封装饼图数据
    userToken = request.GET.get("userToken", None)
    dataIndex = int(request.GET.get("dataIndex", '1')) - 1
    genName = request.GET.get("genName", None)

    #获取染色体条数
    chromeLen = initscore.chromoMap.get(genName)

    #封装legenddata 包括4个不匹配 和染色体条数
    legenddata = []
    legenddata.append("0-misMatch")
    legenddata.append("1-misMatch")
    legenddata.append("2-misMatch")
    legenddata.append("3-misMatch")
    legenddata.append("4-misMatch")


    for i in range(chromeLen):
        legenddata.append("chr"+str(i+1))

    #读取用户结果文件
    json_str = ""
    # 读取用户数据
    # 根据userToken
    basepath = initscore.config_dict.get("basepath")
    dirPath = basepath + userToken

    #userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/" + userToken + ".txt"
    userTokenFileName = dirPath + "/analysisResults.txt"
    # 遍历输出的错配文件
    with open(userTokenFileName, 'r') as destination:
        for line in destination:
            json_str += line
    # 读取文件
    datalist = json.loads(json_str)
    data = datalist[dataIndex]
    seriDatas1 = []
    seriData = {"value": int(data["zeromismatch"]), "name": "0-misMatch"}
    seriDatas1.append(seriData)
    seriData = {"value":int(data["onemismatch"]),"name":"1-misMatch"}
    seriDatas1.append(seriData)
    seriData = {"value": int(data["twomismatch"]), "name": "2-misMatch"}
    seriDatas1.append(seriData)
    seriData = {"value": int(data["threemismatch"]), "name": "3-misMatch"}
    seriDatas1.append(seriData)
    seriData = {"value": int(data["fourmismatch"]), "name": "4-misMatch"}
    seriDatas1.append(seriData)

    seriDatas2 = []
    for i in range(chromeLen):
        seriData = {"value": 0, "name":"chr"+str(i+1) }
        seriDatas2.append(seriData)

    #封装在染色体上的不匹配数
    misMatchs = data["misRNA"]
    for misRNA in misMatchs:
        misChrome = misRNA["misChrome"]
        try:
            misChrNum = int(misChrome[3:4]) - 1 #得到是在哪个染色体上的错配 得到后相应的染色体错配数加1
            seriDatas2[misChrNum]["value"] = int(seriDatas2[misChrNum]["value"]) + 1
        except Exception as e:
            print(e)
    return render(request, "figures/sgpienest.html",{"legenddata":json.dumps(legenddata),"seriDatas1":json.dumps(seriDatas1),"seriDatas2":json.dumps(seriDatas2)})


def sgRNAPie(request):
    #封装饼图数据
    userToken = request.GET.get("userToken", None)
    dataIndex = int(request.GET.get("dataIndex", '1')) -1

    '''
    json_str = ""
    # 读取用户数据
    userTokenFileName = "/Users/FryTsui/anna/study/cas-finder/results/" + userToken + ".txt"
    # 遍历输出的错配文件
    with open(userTokenFileName, 'r') as destination:
        for line in destination:
            json_str += line
    # 读取文件
    datalist = json.loads(json_str)
    data = datalist[dataIndex]'''

    #封装

    return render(request, "figures/sgRNAPie.html")


def sgRNAsingle(request):


    #封装各个染色体上脱靶的位置和脱靶率

    return render(request, "figures/sgscattersingle.html")

