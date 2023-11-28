import functools

from django.contrib.auth.decorators import login_required
from django.core import serializers
from django.http import JsonResponse
from django.shortcuts import render,HttpResponse

def login(request):

  return  render(request, "sglogin.html")

def register(request):
    return render(request, "register.html")