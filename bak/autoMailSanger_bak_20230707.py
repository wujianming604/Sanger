#!/usr/bin/env python
# coding: utf-8
# update: 20210628 
# content: 如果样本是研发样本或者臻智选样本，直接分析，不用再凑家系
# update: 20210731
# content: 替换压缩包目录
# update: 20220701
# content: pathlib搜索目录里面的ab1文件
# update: 20230206
# content: 修改chrom和pos列 为基因和检测区域
# update: 20230217
# content: 修改png图片名称
# update: 20230516
# content: 添加sanger返回时间 , 以png图片为准（同时，只有先证者的数据有，才会更新）
# update: 20230525
# content: 兼容[检测区域]里面，有- or _的情况

from email.header import decode_header
from html.parser import HTMLParser
from email.utils import parseaddr
import imaplib, email, os, poplib
from email.parser import Parser
from urllib.parse import quote
from pathlib import Path
from icecream import ic
import urllib.request
import pandas as pd
import pysnooper
import pymssql
import logging
import random
import zipfile
import shutil
import time
import glob
import re
import os
#os.chdir("/share_data/wujm/project/Django/Project/uploadSanger")
os.chdir("/share_data/clin_result/clin_epilepsy_result/Sanger/downLoad")
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
nowTime = time.strftime("%Y-%m-%d", time.localtime())

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

#获取系统里贝安臻和臻智选的所有先证者样本
def getAllSample():
    lis = []
    conn = pymssql.connect('47.101.134.95','whn','nyuen2018', 'Gene',charset='utf8')
    cursor = conn.cursor()
    cursor.execute("SELECT InfoId FROM dbo.J_SampleRegister WHERE  ProductType = 243 or ProductType = 249 ")
    allSampleNames = cursor.fetchall()
    for i in allSampleNames:
        lis.append(i[0].strip())
    conn.close()

    return lis

#更新Sanger已返回 状态
# 215：计算分析完成；
# 216：报告生成中
# 217：报告生成完成
# 309 数据重分析完成
# 323 Sanger已返回
def update_bioAnaState(sampleName):
    conn = pymssql.connect('47.101.134.95','whn','nyuen2018', 'Gene',charset='utf8')
    cursor = conn.cursor()
    cursor.execute("SELECT inner_state FROM sample_product_receive WHERE sample_num='{}'".format(sampleName))
    state = cursor.fetchone()[0]
    #如果报告生成完成，就不需要更新状态20230605
    if state != 217:
        cursor.execute("UPDATE sample_product_receive set inner_state=323 WHERE sample_num='{}' and ( product_type = 243 or product_type = 249 ) and product != '{}'".format(sampleName,'抗癫痫药物分析'))
    #设置线粒体
        conn.commit()
    conn.close()

#更新Sanger已返回的时间
def updateBAZjindubiao(sampleName):
    conn2 = pymysql.connect('192.168.99.7', 'wujm','wjM123456++', 'clinepilepsy_pipeline',charset='utf8')
    cursor = conn2.cursor()
    cursor.execute("update Genetic_Counseling_processanalysis set sanger_complete_time='{}' where sample_name='{}';".format(nowTime, sampleName))
    conn2.commit()
    conn2.close()

def updateZZXjindubiao(sampleName):
    conn2 = pymysql.connect('192.168.99.7', 'wujm','wjM123456++', 'clinepilepsy_pipeline',charset='utf8')
    cursor = conn2.cursor()
    cursor.execute("update Genetic_Counseling_zzxprocessanalysis set sanger_complete_time='{}' where sample_name='{}';".format(nowTime, sampleName))
    conn2.commit()
    conn2.close()

#获取亲属及患者的样本编号
def get_dict():
    sampleDict = {}# 从300例之后开始
    logging.info("开始连接PM数据库，获取亲属样本信息。")
    sql_con = pymysql.connect(host='192.168.99.7', user='wujm', passwd='wjM123456++', db='clinepilepsy_pipeline')
    sql_cur = sql_con.cursor()
    cmd = "select type_id, user_name, sample_name, kinsfolk_sample from PM_pm;"
    sql_cur.execute(cmd)
    infos = sql_cur.fetchall()
   
    ##先判断有没有亲属送样
    for each in infos[300:]:
        if each[-1] == None or each[-1] == "Empty" or each[-1].strip() == '':
            logging.warning("%s-这个家系没有父母样本" % each[2])
            sampleDict[each[2]] = each[0]+'_'+each[1]+'_'+each[2]
        else:
            #print(each)
            sampleDict[each[2]] = each[0]+'_'+each[1]+'_'+each[2]
            for jias in re.split("[,，]",each[-1].strip()):
                jias = re.split("[:：]",jias.strip())
                sampleDict[jias[1]] = each[0]+'_'+each[1]+'_'+each[2]
  
    return sampleDict


##查看所有的压缩包上传列表，看看是否已存在

@pysnooper.snoop(output=nowTime+".getSangerZipList.log",prefix="getSangerZipList:")
def getSangerZipList(rawSampleName, newSampleName):
    #randomNumber = random.randint(1,11)
    randomNumber = 0
    nowTime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    
    my_obj = Mysql()
    my_con = my_obj.connect()
    my_cur = my_obj.cursor(my_con)
    get_sample_cmd = "select fileName from Sanger_sangerzipfiles;"
    my_obj.run_cmd(my_cur, get_sample_cmd)
    infos = my_cur.fetchall()
    sampleList = [i[0] for i in infos]
    #剩余结果全部当作新样本来分析
    if "剩余结果" in rawSampleName:
        insertCMD = "insert into Sanger_sangerzipfiles set fileName = '%s', sampleNumbers = %d, status = 0,createdTime = '%s', analysisStatus = 1, updateTime = '%s';"  % (newSampleName, randomNumber, nowTime, nowTime)
        my_obj.run_cmd(my_cur, insertCMD)
        logging.info("插入%s至数据库" % rawSampleName)
        my_con.commit()
    #如果压缩包在库里存在，跳过不分析
    
    status = 0
    for sample in sampleList:
        if rawSampleName in sample:
            logging.warning("%s-该压缩包里的样本已经分析过了。" % rawSampleName)
            status = 1
    
    if status == 0:
        insertCMD = "insert into Sanger_sangerzipfiles set fileName = '%s', sampleNumbers = %d, status = 0,createdTime = '%s', analysisStatus = 1, updateTime = '%s';"  % (newSampleName, randomNumber, nowTime, nowTime)
        my_obj.run_cmd(my_cur, insertCMD)
        logging.info("插入%s至数据库" % rawSampleName)
        my_con.commit()
    my_con.close()
    
    
#检查表
@pysnooper.snoop(output=nowTime+".checkXLSX.log",watch="data")
def checkXLSX(xlsxFile):
    xls_data = xlrd.open_workbook(xlsxFile)
    sheetName = xls_data.sheet_names() #暂考虑只有一张Sheet1
    data = pd.read_excel(xlsxFile,sheet_name = sheetName[0],index_col = None,na_values= ['9999'])

    #检查标题是否符合要求
    targetTitle = "样本编号 样本名 基因 检测区域 检测位点 参考碱基 检测结果 突变类型 测序方向 峰图位置 峰图文件"
    dataTitle = ' '.join(i for i in data.columns.values)
    print(targetTitle)
    print(dataTitle)
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



def updateSangersamples(List):
    nowTime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    my_obj = Mysql()
    my_con = my_obj.connect()
    my_cur = my_obj.cursor(my_con)
    no_ana_cmd = "select sampleName, chrom, pos, genoType from Sanger_sangersamples;"
    my_obj.run_cmd(my_cur, no_ana_cmd)
    allSampleNames = my_cur.fetchall()
    for each in List:
        sampleName, chrom, pos, GT, ab1File, ab1Png = each.split(',')
        #判断sampleName，chrom, pos，GT是否存在，如果存在，更新
        if (sampleName, chrom, pos, GT) in allSampleNames:
            #update
            update_cmd = "update Sanger_sangersamples set abiFileUrl='{abiFile}', abiPngUrl='{abiPng}', update_time='{updateTime}'  where sampleName='{name}' and chrom='{chr}' and pos={pos} and genoType='{gt}';".format(abiFile=ab1File,abiPng=ab1Png,updateTime=nowTime, name=sampleName, chr=chrom, pos=pos, gt=GT)
            my_obj.run_cmd(my_cur, update_cmd)
            my_con.commit()
        else:
            #insert
            print("insert into Sanger_sangersamples (sampleName, chrom, pos, genoType, abiFileUrl, abiPngUrl, created_time, update_time) values ('{sample}','{chr}',{pos},'{gt}','{abiFile}','{abiPng}','{createTime}','{updateTime}');".format(sample=sampleName,chr=chrom,pos=pos,gt=GT,abiFile=ab1File,abiPng=ab1Png,createTime=nowTime, updateTime=nowTime))
            insert_cmd = "insert into Sanger_sangersamples (sampleName, chrom, pos, genoType, abiFileUrl, abiPngUrl, created_time, update_time) values ('{sample}','{chr}','{pos}','{gt}','{abiFile}','{abiPng}','{createTime}','{updateTime}');".format(sample=sampleName,chr=chrom,pos=pos,gt=GT,abiFile=ab1File,abiPng=ab1Png,createTime=nowTime, updateTime=nowTime)
            my_obj.run_cmd(my_cur, insert_cmd)
            my_con.commit()
    my_con.close()
    


#分析数据

@pysnooper.snoop(output=nowTime+".analysisSanger.log",watch="resultList")
def analysisSanger(rawSampleName, zipFile):
    status = 0
    sampleDict = get_dict()
    allSampleList = getAllSample()
    #不管是Django上传，压缩包存放在下面目录
    #defaultPath = "/share_data/wujm/project/Django/Project/uploadSanger"
    #程序自动识别邮箱并下载，压缩包存放在下面目录
    defaultPath = "/share_data/clin_result/clin_epilepsy_result/Sanger/downLoad"
    #allResultPath = "/share_data/wujm/project/Django/Project/uploadSanger/sangerResult/"
    allResultPath = "/share_data/clin_result/clin_epilepsy_result/Sanger/sangerResult/"
    os.chdir(defaultPath)   #切换上传文件所在目录里面
    
    timePath = time.strftime("%Y%m%d%H%M%S", time.localtime())
    if not os.path.exists(timePath):
        os.mkdir(timePath)  #生成日期目录，
    else:
        pass
    os.chdir(timePath)
    
    zipFileAbosultePath = defaultPath + '/' + zipFile   #压缩文件的绝对路径
    newPath = defaultPath +'/'+timePath   #新建的日期目录，上传的压缩包要 复制 至该目录下进行分析。
    
    ##默认上传或者下载的  都是zip文件格式
    shutil.move(zipFileAbosultePath, zipFile)
    os.system("unzip -o -q %s" % zipFile)
    
    ##解压之后，进入到目录里面
    
    #os.chdir(rawSampleName)
    p = Path(".")
    for f in p.iterdir():
        if f.is_dir():
            newDirName = str(f).replace(" ","")
            os.rename(str(f),newDirName)  #为了去除解压缩之后的空格

            #切换目录，开始分析
            os.chdir(newDirName)
            #获取xlsx名称
            xlsxFiles = glob.glob('*.xlsx')
            if len(xlsxFiles) == 0:
                logging.error("该目录内没有xlsx结果表, 退出！")
                status = 1
                
            logging.info("开始打开{}-{}进行分析".format(zipFile,xlsxFiles))
            for eachFile in xlsxFiles:
                data = checkXLSX(eachFile)
                #接下来开始进行正常分析了。
                resultList = []
                data['sampleId'] = data['样本名']
                data['chrom'] = data['检测区域'].str.split(":",0).str[0]
                data['pos'] = data['检测区域'].str.split(":",0).str[1]
                data['pos'].replace('[A-Z]|[a-z]|>','',regex=True,inplace=True)
                data['geneType'] = data['突变类型']
                data['ab1Name'] = data['峰图文件']
                #data['gene'] = data['检测位点']
                data['gene'] = data['基因']
                data['fangx'] = data['测序方向']
                data['area'] = data['检测区域']
                sampleInfos = data[["sampleId","chrom","pos","geneType","ab1Name","gene","fangx","area"]]
                for index, row in sampleInfos.iterrows():
                    sampleId = row['sampleId']
                    #update20220701 ,用p.rglob找到ab1文件,复制到当前目录下
                    for i in p.rglob(row['ab1Name']):
                        os.system(f"cp {i} .")
                    chromId = row['chrom'].replace("chr","")
                    #考虑有-的情况
                    if '-' in row['pos']:
                        posId = int(row['pos'].split('-')[0])
                    elif '_' in row['pos']:
                        posId = int(row['pos'].split('_')[0])
                    else:
                        posId = row['pos']
                    if row['fangx'] == "正向":
                        fangx = "F"
                    else:
                        fangx = "R"
                    upStart = int(posId) - 20  #往前截取30bp
                    upEnd = int(posId) - 1
                    downStart = int(posId) + 1
                    downEnd = int(posId) + 20  #往后截取30bp
                    #获取snp位点前后30bp碱基序列
                    upSeq = bio_function.extract_fa(str(chromId),str(upStart),str(upEnd)).upper()
                    print("upSeq:",upSeq,"长度",len(upSeq),flush=True)
                    downSeq = bio_function.extract_fa(str(chromId),str(downStart),str(downEnd)).upper()
                    print("downSeq:",downSeq,"长度",len(downSeq),flush=True)
                    sampleAb1 = row['ab1Name']
                    if row['gene'] == "":
                        geneName = "NA"
                    else:
                        geneName = row['gene']
                    ic("正在分析{}-{}。".format(sampleId, row['geneType']))
                    #将ab1文件名前缀作为图片的name
                    #rawName = os.path.splitext(sampleAb1)[0].split("-")
                    #pngName = sampleId + '_' + row['gene'] + '_' + fangx + '.png'
                    #更新图片的命名方式20230217
                    pngName = sampleId + '_' + geneName + '_' + row['area'] + '.png'
                    pngName = pngName.replace(":","_").replace(">","_")
                    #定位sample_id 对应的ab1文件
                    logging.info("'ab1':'{}'--'png':'{}'--'upSeq':'{}'--'downSeq':{}.".format(sampleAb1,pngName,upSeq,downSeq))
                    #运行程序
                    logging.info("开始Sanger截图...")
                    os.system("/usr/bin/Rscript /share_data/wujm/Config/script/clinepilepsy/plot_by_sangerseqR-copy.R %s %s %s %s" % (sampleAb1,pngName,upSeq,downSeq))
                    #print("/usr/bin/Rscript /share_data/wujm/Config/script/clinepilepsy/plot_by_sangerseqR-copy.R %s %s %s %s" % (sampleAb1,pngName,upSeq,downSeq))
                    #程序运行结束之后，将sampleId相关的文件全放到同一个sampleId目录下,非贝安臻样本直接新建样本名目录
                    if sampleId not in sampleDict:
                        logging.warning("{} 不在数据库里面，可能是臻智选样本还是研发样本。".format(sampleId))
                        sampleIdPath = allResultPath + sampleId
                        if not os.path.exists(sampleIdPath):
                            os.mkdir(sampleIdPath)
                        else:
                            pass
                        #模糊匹配
                        localSampleList = glob.glob('*%s*[ab1,png]' % sampleId)
                        logging.info("Sanger截图完毕，开始截图结果拷贝至sangerResult目录...")
                        for i in localSampleList:
                            os.system("cp %s %s" % (i, sampleIdPath))
                        ab1AboPath = allResultPath+sampleId+'/'+sampleAb1
                        newPngAboPath = allResultPath+sampleId+'/'+pngName
                        if sampleId in allSampleList:
                            updateZZXjindubiao(sampleId)  #更新臻智选sanger返回时间
                            update_bioAnaState(sampleId)  #更新sanger返回状态
                    else:
                        logging.info("'{}'-的家系信息为{}".format(sampleId,sampleDict[sampleId]))
                        sampleIdPath = allResultPath + sampleDict[sampleId]
                        if not os.path.exists(sampleIdPath):
                            os.mkdir(sampleIdPath);#生成日期目录，
                        else:
                            pass
                        #模糊匹配
                        localSampleList = glob.glob('*%s*[ab1,png]' % sampleId)
                        logging.info("Sanger截图完毕，开始截图结果拷贝至sangerResult目录...")
                        for i in localSampleList:
                            os.system("cp %s %s" % (i, sampleIdPath))
                        ab1AboPath = allResultPath+sampleDict[sampleId]+'/'+sampleAb1
                        newPngAboPath = allResultPath+sampleDict[sampleId]+'/'+pngName
                        if sampleId in allSampleList:
                            updateBAZjindubiao(sampleId)  #更新贝安臻sanger返回时间
                            update_bioAnaState(sampleId)  #更新sanger返回状态
                    #resultList.append(str(sampleId)+','+chromId+','+str(posId)+','+row['geneType']+','+ab1AboPath+','+newPngAboPath)
                    print(str(sampleId)+','+geneName+','+row['area']+','+row['geneType']+','+ab1AboPath+','+newPngAboPath)
                    
                    resultList.append(str(sampleId)+','+geneName+','+row['area']+','+row['geneType']+','+ab1AboPath+','+newPngAboPath)
                    status = 1
                os.system("chmod 777 %s" % allResultPath)
                os.system("chmod -R 777  %s/*" % allResultPath)
            
            #进行数据库更新
            updateSangersamples(resultList)
            return status


    
    
class MyHTMLParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.data = []   # 定义data数组用来存储html中的数据
        self.links = [] 
            
    def handle_starttag(self, tag, attrs):
        print('<%s>' % tag)
        if tag == "a":
            if len(attrs) == 0: pass
            else:
                for (variable, value)  in attrs:
                    if variable == "href":
                        self.links.append(value)
     
    def handle_endtag(self, tag):
        print('</%s>' % tag)
 
    def handle_startendtag(self, tag, attrs):
        print('<%s/>' % tag)
 
    def handle_data(self, data):
        print('data===>', data)
 
    def handle_comment(self, data):
        print('<!--', data, '-->')
 
    def handle_entityref(self, name):
        print('&%s;' % name)
 
    def handle_charref(self, name):
        print('&#%s;' % name)

def guess_charset(msg):
    charset = msg.get_charset()
    if charset is None:
        content_type = msg.get('Content-Type', '').lower()
        pos = content_type.find('charset=')
        if pos >= 0:
            charset = content_type[pos + 8:].strip()
    return charset

def decode_str(s):
    value, charset = decode_header(s)[0]
    if charset:
        value = value.decode(charset)
    return value

def print_info(msg, indent=0):
    if indent == 0:
        for header in ['From', 'To', 'Subject']:
            value = msg.get(header, '')
            if value:
                if header == 'Subject':
                    value = decode_str(value)
                elif header == 'From':
                    value = decode_str(value)
                    if value != "迈浦测序部":
                        continue
                else:
                    hdr, addr = parseaddr(value)
                    name = decode_str(hdr)
                    value = u'%s <%s>' % (name, addr)

    if msg.is_multipart():
        parts = msg.get_payload()
        for n, part in enumerate(parts):
            print_info(part, indent + 1)
    else:
        content_type = msg.get_content_type()
        if content_type == 'text/plain' or content_type == 'text/html':
            content = msg.get_payload(decode=True)
            charset = guess_charset(msg)
            if charset:
                content = content.decode(charset)
                parserStr(content)
        else:
            logging.error(content_type)
            
def getUrl(urlName, newSampleName):
    urllib.request.urlretrieve(urlName, newSampleName)
            
            
def parserStr(content):
    adjunct = re.findall(r'<li><a href=\"(.*)\">2023', content)
    if adjunct:
        for each in adjunct:
            #替换里面的字符串
            each = each.replace("amp;","")
            #print(each)
            #后缀为月 日 时
            fn = time.strftime("%m%d%H")
            #原始样本名称
            rawSampleName = each.split("=")[-1].replace(".zip","")
            #加了日期后缀的的名称
            newSampleName =  rawSampleName+"_"+fn+".zip"
            #url中包含ASCII以外的字符 使用 quote(‘汉字’）解决
            urlName = "=".join(i for i in each.split("=")[:-1])+"="+quote(each.split("=")[-1])
            #下载文件，另存在sampleName
            logging.info("开始下载文件--{}".format(urlName))
            getUrl(urlName, newSampleName)
            time.sleep(5)
            logging.info("下载文件完成--{}".format(newSampleName))
            
            #解析完了之后，将压缩包信息上传至sangerzipfiles
            logging.info("开始分析--{}".format(newSampleName))
            getSangerZipList(rawSampleName, newSampleName)
            logging.info("分析完毕--{}".format(newSampleName))
            #上传完成之后开始分析
    
# In[26]:

@pysnooper.snoop(output=nowTime+".run.log",prefix="run:")
def run():
    #判断Sanger_sangerzipfiles表中的analysisStatus状态是否为1，如果为1，表示未分析，需要进行run(FILE)分析； 如果是0 ，表示已分析，pass。
    my_obj = Mysql()
    my_con = my_obj.connect()
    my_cur = my_obj.cursor(my_con)
    no_ana_cmd = "select id, fileName from Sanger_sangerzipfiles where analysisStatus = 1"
    my_obj.run_cmd(my_cur, no_ana_cmd)
    info = my_cur.fetchall()
    if len(info) == 0:
        pass
    else:
        for index in range(len(info)):
            no_ana_id = info[index][0]
            no_ana_filename = info[index][1]
            rawName = "_".join(i for i in no_ana_filename.split("_")[:-1])
            logging.info("开始分析{}".format(no_ana_filename))
            
            status = analysisSanger(rawName, no_ana_filename)
            #运行完成之后，更新analysisStatus为0
            if status == 1:
                logging.info("{}--analysisStatus更新为0".format(no_ana_id))
                update_cmd = "update Sanger_sangerzipfiles set analysisStatus = 0 where id={};".format(no_ana_id)
                my_obj.run_cmd(my_cur, update_cmd)
                my_con.commit()
    my_con.close() 
    



# In[28]:


#update time 20210419
#由于邮箱里NY开头的样本检索不到，又或者有额外的sanger截图需求，所以添加：
#检索当前目录是否有zip结尾的样本，如果有，则加入分析。（因为分析过的样本都会移动至相应目录里面，所以在当前目录下有zip文件，代表是需要分析的sanger）
#先分析本地的
os.chdir("/share_data/clin_result/clin_epilepsy_result/Sanger")
localZips = glob.glob('*.zip')
print(localZips)

for localSample in localZips:
    logging.warning("在本地，有需要分析的zip文件----{}".format(localSample))
    fn = time.strftime("%m%d%H")
    rawSampleName = localSample.replace(".zip","")
    #加了日期后缀的的名称
    newSampleName =  rawSampleName+"_"+fn+".zip"
    #将文件移动至downLoad
    shutil.move(localSample, "/share_data/clin_result/clin_epilepsy_result/Sanger/downLoad/"+newSampleName)
    #解析完了之后，将压缩包信息上传至sangerzipfiles
    logging.info("开始分析本地文件--{}".format(newSampleName))
    
    getSangerZipList(rawSampleName, newSampleName)
    logging.info("本地文件分析完毕--{}".format(newSampleName))

os.chdir("downLoad")

number = len(server.list()[1])
#选择前5个邮箱
if sys.argv[1]:    
   selectNum = int(sys.argv[1]) -1 
else:
   selectNum = 4
#messages = [server.retr(i) for i in range(number-4, number)]  5个
messages = [server.retr(i) for i in range(number-selectNum, number+1)]
messages = [b'\r\n'.join(mssg[1]).decode('gbk') for mssg in messages]
messages = [Parser().parsestr(mssg) for mssg in messages]
messages = messages[::-1]
for i in messages:
    print_info(i)

run() 
