import json
import random

from django.http import HttpResponse
from django.shortcuts import render



def sgScattermatrix(request):

    #生成数据
    rawData = []

    #循环生成数据

    for i in range(50):
        data = []
        #生成cnn 数据
        #生成150条数据 每个模型50条 6个物种 human mouse monkey rice wheat corn
        #[55, 9, 56, 0.46, 18, 6, 'good', 'CNN'],

        #生成human 预测值
        humanScore = random.randint(40, 100)
        mouseScore = random.randint(50, 70)
        monkeyScore = random.randint(60, 80)
        riceScore = random.randint(0, 20)
        cornScore = random.randint(10, 50)
        wheatScore = random.randint(0, 30)

        data.append(humanScore)
        data.append(mouseScore)
        data.append(monkeyScore)
        data.append(riceScore)
        data.append(cornScore)
        data.append(wheatScore)
        averageScore = 0
        averageScore = (humanScore + mouseScore + monkeyScore + riceScore + cornScore + wheatScore) / 6
        if averageScore > 70:
            data.append('poor')
        if averageScore > 50 and averageScore < 70:
            data.append('good')
        if averageScore < 50:
            data.append('excellent')#poor
        data.append('CNN')

        rawData.append(data)

    for i in range(50):
        data = []
        # 生成rnn 数据
        # 生成150条数据 每个模型50条 6个物种 human mouse monkey rice wheat corn
        # [55, 9, 56, 0.46, 18, 6, 'good', 'CNN'],
        humanScore = random.randint(20, 80)
        mouseScore = random.randint(30, 80)
        monkeyScore = random.randint(40, 100)
        riceScore = random.randint(10, 30)
        cornScore = random.randint(20, 40)
        wheatScore = random.randint(10, 30)
        data.append(humanScore)
        data.append(mouseScore)
        data.append(monkeyScore)
        data.append(riceScore)
        data.append(cornScore)
        data.append(wheatScore)
        averageScore = 0
        averageScore = (humanScore + mouseScore + monkeyScore + riceScore + cornScore + wheatScore )/ 6
        if averageScore > 70:
            data.append('poor')
        if averageScore > 50 and averageScore < 70:
            data.append('good')
        if averageScore < 50:
            data.append('excellent')
        data.append('RNN')
        rawData.append(data)

    for i in range(50):
        data = []
        # 生成transformer 数据
        # 生成150条数据 每个模型50条 6个物种 human mouse monkey rice wheat corn
        # [55, 9, 56, 0.46, 18, 6, 'good', 'CNN'],
        humanScore = random.randint(30, 70)
        mouseScore = random.randint(20, 80)
        monkeyScore = random.randint(40, 100)
        riceScore = random.randint(40, 90)
        cornScore = random.randint(40, 100)
        wheatScore = random.randint(20, 100)
        data.append(humanScore)
        data.append(mouseScore)
        data.append(monkeyScore)
        data.append(riceScore)
        data.append(cornScore)
        data.append(wheatScore)
        averageScore = 0
        averageScore = (humanScore + mouseScore + monkeyScore + riceScore + cornScore + wheatScore )/ 6
        if averageScore > 70:
            data.append('poor')
        if averageScore > 50 and averageScore < 70:
            data.append('good')
        if averageScore < 50:
            data.append('excellent')
        data.append('Transformer')
        rawData.append(data)

    datas = json.dumps(rawData)
    return render(request, "pictures/scatter-matrix.html", {"datas":datas})




def datasetEncode(request):
    rawData = []
    datas = json.dumps(rawData)
    return render(request, "pictures/dataset-encode.html", {"datas": datas})


def datasetEncodeJson(request):
    file = open(
        '/Users/FryTsui/anna/study/pmworkspace/crispr-designer/crisprDesigner/sgRnaDesigner/static/json/pis.json',
        'r')
    json_data = file.read()

    file.close()
    return HttpResponse(json_data, content_type="application/json", status=200)



def heatMap(request):
    rawData = []
    datas = json.dumps(rawData)
    return render(request, "pictures/heatmap-cartesian.html", {"datas": datas})



def cluster(request):
    rawData = []
    datas = json.dumps(rawData)
    return render(request, "pictures/scatter-clustering.html", {"datas": datas})







