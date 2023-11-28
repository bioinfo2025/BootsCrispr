import json
import smtplib
import datetime
import uuid
from email.header import Header
from email.mime.text import MIMEText
from random import random

from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods

import sgRnaDesigner.utils.initSeqScore as initscore





@require_http_methods(["POST"])
@csrf_exempt
def contact(request):
    if request.method == "POST":
        userEmail = request.POST.get("userEmail", None)
        emailConent = request.POST.get("emailConent", None)


        smtp_server = initscore.config_dict.get("smtp_server")
        username = initscore.config_dict.get("username")
        password = initscore.config_dict.get("password")
        smtp_port = initscore.config_dict.get("smtp_port")

        content = 'Thank you for your letter,'+\
        'We will provide feedback on your issue in the future!!!'
        subject = 'Letter of thanks from' + initscore.config_dict.get("sys_name")

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

    
        contactLetterPath = initscore.config_dict.get("contactLetterPath")

      

        now = datetime.datetime.now()

        fileName =str(now.year)+"-"+str(now.month)+"-"+str(now.day)+str(uuid.uuid1())+".txt"

     
        with open(fileName, 'w') as destination:
            destination.write(userEmail + "\n")
            destination.write(str(emailConent))
        destination.close()

    return JsonResponse({"code": 200}, status=200, content_type='application/json; charset=utf-8')
