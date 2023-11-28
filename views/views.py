from django.core import serializers
from django.http import JsonResponse
from django.shortcuts import render,HttpResponse


def page_not_found(request, exception=None):
    """
    404 页面
    :param request:
    :return:
    """
    return render(request, "404.html", status=404)


def server_error(request, exception=None):
    """
    500 页面
    :param request:
    :return:
    """
    return render(request, "404.html", status=500)








