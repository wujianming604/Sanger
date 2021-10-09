#!/usr/bin/env python
# coding: utf-8

# In[27]:


from email.header import decode_header
from html.parser import HTMLParser
from email.utils import parseaddr
import imaplib, email, os, poplib
from email.parser import Parser
from urllib.parse import quote
from pathlib import Path
import urllib.request
import pandas as pd
import logging
import random
import zipfile
import shutil
import time
import glob
import re
import os
os.chdir("/share_data/clin_result/clin_epilepsy_result/Sanger")
import bio_function
import sys

try:
    import xlrd
except Exception as import_error:
    os.system('pip install xlrd')
    import xlrd
try:
    import pandas
except Exception as import_error:
    os.system('pip install pandas')
    import pandas as pd
try:
    import pymysql
except Exception as import_error:
    os.system('pip install pymysql')
    import pymysql

##Logging配置
logFormat = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(filename='run.log', level=logging.DEBUG, format=logFormat)
##基本配置
email = "wujianming@nhwa-group.com"
password = "benmeiZHANYUE604"
pop3_server = "mail.nhwa-group.com"

server = poplib.POP3(pop3_server)

server.user(email)
server.pass_(password)

class Mysql:
    def __init__(self, sql_user='NA', sql_pswd='NA'):
        self.user = 'wujm' if sql_user == 'NA' else sql_user
        self.pswd = 'wjM123456++' if sql_pswd == 'NA' else sql_pswd

    def connect(self, sql_db='clinepilepsy_pipeline'):
        try:
            sql_con = pymysql.connect(host='192.168.99.7', user=self.user, passwd=self.pswd, db=sql_db)
        except Exception as connect_error:
            logging.error('连接数据库-%s错误！' % sql_db)
            logging.error(connect_error)
            sys.exit()
        return sql_con

    # 创建游标
    def cursor(self, sql_con):
        sql_cur = sql_con.cursor()
        return sql_cur

    # 获取
    def run_cmd(self, sql_cur, sql_cmd):
        try:
            sql_cur.execute(sql_cmd)
        except Exception as run_error:
            logging.error('运行sql命令问题:\n%s' % sql_cmd)
            logging.error(run_error)
            sys.exit()
            

    
#检查表
def checkXLSX(xlsxFile):
    xls_data = xlrd.open_workbook(xlsxFile)
    sheetName = xls_data.sheet_names() #暂考虑只有一张Sheet1
    data = pd.read_excel(xlsxFile,sheet_name = sheetName[0],index_col = None,na_values= ['9999'])

    #检查标题是否符合要求
    targetTitle = "样本编号 样本名 基因 检测区域 检测位点 参考碱基 检测结果 突变类型 测序方向 峰图位置 峰图文件"
    dataTitle = ' '.join(i for i in data.columns.values)
    if dataTitle != targetTitle:
        logging.critical("表格的标题与以前不一致，程序退出！！")
        sys.exit()
    else:
        ##检查是否有多余的空行， 为了防止map随便写，现制定"检测位点"和"峰图文件"有空行的列，就删除
        data.dropna(subset=['测序方向'],inplace=True)
        data.dropna(subset=['峰图文件'],inplace=True)

        ##检查检测区域,已发现的一个问题是：除了区域外，还有其他字符串，以空格分隔
        data['检测区域'] = data['检测区域'].str.split(' ',expand=True)[0]

        ##再检查峰图文件，有的时候，没有ab1后缀
        for index, row in data.iterrows():
            if "ab1" not in row['峰图文件']:
                logging.warning("%s-样本的ab1文件没有后缀！"% row['样本编号'])
                data.loc[index, '峰图文件'] = row['峰图文件'] + '.ab1'
            if "del" in row['检测区域'] and  "ins" in row['检测区域']:
                data.loc[index,'检测区域'] = row['检测区域'].split("_")[0]
    return data





#分析数据
def analysisSanger(rawSampleName, zipFile):

    handleResultPath = "/share_data/clin_result/clin_epilepsy_result/Sanger/otherSangerResult/"
    os.chdir(handleResultPath)
    os.system("unzip -o -q %s" % zipFile)
    
    ##解压之后，进入到目录里面
    os.chdir(rawSampleName)
    
    #获取xlsx名称
    xlsxFiles = glob.glob('*.xlsx')
    if len(xlsxFiles) == 0:
        logging.error("该目录内没有xlsx结果表, 退出！")
        sys.exit()
        
    logging.info("开始打开{}-{}进行分析".format(zipFile,xlsxFiles))
    for eachFile in xlsxFiles:
        data = checkXLSX(eachFile)
        #接下来开始进行正常分析了。
        data['sampleId'] = data['样本名']
        data['chrom'] = data['检测区域'].str.split(":",0).str[0]
        data['pos'] = data['检测区域'].str.split(":",0).str[1]
        data['pos'].replace('[A-Z]|[a-z]|>','',regex=True,inplace=True)
        data['geneType'] = data['突变类型']
        data['ab1Name'] = data['峰图文件']
        data['gene'] = data['检测位点']
        data['fangx'] = data['测序方向']
        sampleInfos = data[["sampleId","chrom","pos","geneType","ab1Name","gene","fangx"]]
        for index, row in sampleInfos.iterrows():
            sampleId = row['sampleId']
            chromId = row['chrom'].replace("chr","")
            posId = row['pos']
            if row['fangx'] == "正向":
                fangx = "F"
            else:
                fangx = "R"
            upStart = int(posId) - 30  #往前截取30bp
            upEnd = int(posId) - 1
            downStart = int(posId) + 1
            downEnd = int(posId) + 30  #往后截取30bp
            #获取snp位点前后30bp碱基序列
            upSeq = bio_function.extract_fa(str(chromId),str(upStart),str(upEnd)).upper()
            downSeq = bio_function.extract_fa(str(chromId),str(downStart),str(downEnd)).upper()
            sampleAb1 = row['ab1Name']
            #将ab1文件名前缀作为图片的name
            #rawName = os.path.splitext(sampleAb1)[0].split("-")
            pngName = sampleId + '_' + row['gene'] + '_' + fangx + '.png'
            #定位sample_id 对应的ab1文件
            logging.info("'ab1':'{}'--'png':'{}'--'upSeq':'{}'--'downSeq':{}.".format(sampleAb1,pngName,upSeq,downSeq))
            #运行程序
            logging.info("开始Sanger截图...")
            os.system("/usr/bin/Rscript /share_data/wujm/Config/script/clinepilepsy/plot_by_sangerseqR-copy.R %s %s %s %s" % (sampleAb1,pngName,upSeq,downSeq))
            #程序运行结束之后，将sampleId相关的文件全放到同一个sampleId目录下


os.chdir("/share_data/clin_result/clin_epilepsy_result/Sanger")
localZips = glob.glob('*.zip')

for localSample in localZips:
    logging.warning("在本地，有需要分析的zip文件----{}".format(localSample))
    logging.warning("在本地，有需要分析的zip文件----{}".format(localSample))
    fn = time.strftime("%m%d%H")
    rawSampleName = localSample.replace(".zip","")
    #将文件移动至downLoad
    shutil.move(localSample, "/share_data/clin_result/clin_epilepsy_result/Sanger/otherSangerResult/"+localSample)
    #解析完了之后，将压缩包信息上传至sangerzipfiles
    logging.info("开始分析本地文件--{}".format(localSample))
    analysisSanger(rawSampleName, localSample)
    logging.info("本地文件分析完毕--{}".format(localSample))