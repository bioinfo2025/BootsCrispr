from django.apps import AppConfig
import sgRnaDesigner.utils.initSeqScore as initscore


class SgrnadesignerConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "sgRnaDesigner"

    def ready(self):
        # 初始化database
        #指定地方
        #initscore.make_database('/Users/FryTsui/anna/study/fasta/chr1_173.txt', '/Users/FryTsui/anna/study/search/', 23)
        #initscore.make_database('/Users/FryTsui/anna/study/fasta/chr1.fa', '/Users/FryTsui/anna/study/fasta/', 23)
        initscore.init_chrome()
        initscore.init_params()
        print("sgRnaDesigner is ready!")
