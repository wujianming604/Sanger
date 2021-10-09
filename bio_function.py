#!/usr/bin/env python3
#coding:utf-8

from __future__ import division
import argparse
import os
import re
import sys
import requests
import time
import json


# 获取随机数
class Random:
    def __init__(self):
        self.web_dict = {
            'url': 'https://api.random.org/json-rpc/2/invoke',
            'api_key': [
                {'key': '783f9fb8-49d0-484f-bf80-071119c43355', 'id': 5208},
            ]
        }
        '''
        self.web_dict = {
            'url': 'https://www.random.org/integers',
            'proxy': {
                'http': '222.218.122.5:9999',
                # 'http': 'http://thomas:qwer159753++@185.201.226.26:9090'
            },
        }
        '''
        self.web_header = {
            'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3',
            'accept-encoding': 'utf-8',  # 注意
            'accept-language': 'zh-CN,zh;q=0.9,zh-TW;q=0.8,en;q=0.7,ja;q=0.6',
            'cache-control': 'max-age=0',
            'sec-fetch-mode': 'navigate',
            'sec-fetch-site': 'none',
            'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_0) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/77.0.3865.120 Safari/537.36',
        }

    def get_integer(self, total_numbers=10, inter_min=1, inter_max=10, replacement='False'):
        if re.search('true', replacement, re.IGNORECASE):
            replacement = True
        else:
            replacement = False
        # web req
        request_dict = {
            'jsonrpc': '2.0',
            'method': 'generateIntegers',
            'params': {
                'apiKey': self.web_dict['api_key'][0]['key'],
                'n': total_numbers,
                'min': inter_min,
                'max': inter_max,
                'replacement': replacement,
                'base': 10,
            },
            'id': self.web_dict['api_key'][0]['id'],
        }
        req_obj = requests.post(self.web_dict['url'], json=request_dict, headers=self.web_header)
        return json.loads(req_obj.text)['result']['random']['data']
        '''
        request_dict = {
            'num': total_numbers,
            'min': inter_min, 'max': inter_max,
            'replacement': replacement,
            'col': 1, 'base': 10,
            'format': 'plain',
        }
        # print(json.dumps(request_dict))
        req_obj = requests.get(
            self.web_dict['url'], params=request_dict,
            headers=self.web_header, proxies=self.web_dict['proxy'],
        )
        if req_obj.status_code == 200:
            return list(map(lambda x: int(x), list(set(re.split('\s+', req_obj.text.strip())))))
        else:
            print(req_obj.text)
            print('无法获取random随机数')
            sys.exit()
        '''


# 格式化染色体
def format_chr(chrom):
    if re.match('ENST', chrom):
        chrom = re.split('\.', re.split('\|', chrom)[0])[0]
    elif re.match('[0-9]|M', chrom):
        chrom = 'chr'+chrom
    #if re.search('23|24|X|Y', chrom):
    #    chrom = 'chrX' if re.search('23|X', chrom) else 'chrY'
    #elif re.match('chrMT|chr25|25|MT', chrom):
    #    chrom = 'chrM'
    if re.search('X|chrX', chrom):
        chrom = 'chrX'
    elif re.search('Y|chrY', chrom):
        chrom = 'chrY'
    elif re.match('chrMT|MT', chrom):
        chrom = 'chrM'
    return chrom



# 提取染色体序列以及rnaxulie
def extract_fa(chrom, pos_start, pos_end, genome_file='NA'):
    # 检测文件是否存在
    if genome_file == 'NA':
        genome_file = '/share_data/wujm/Config/reference/hg19/ucsc.hg19.fasta'
    fai_file = genome_file+'.fai'
    if os.path.exists(genome_file+'.fai') is False:
        #print('基因组fai文件不存在-%s！' % genome_file+'.fai')
        sys.exit()
    # 染色体标准化
    chrom = format_chr(chrom)
    # 检测输入物理位置是否包含非数字字符
    if re.search('\D', pos_start) or re.search('\D', pos_end):
        #print('输入的物理位置包含非数字:%s-%s!' % (pos_start, pos_end))
        sys.exit()
    if int(pos_end) < int(pos_start):
        #print('输入开始位置-%s大于结束位置-%s!' % (pos_end, pos_start))
        sys.exit()
    # 读取fai文件，获取字典
    # head_tag 首序列开始字节 一行字节数目 一行字节加上\n数目
    fai_dict = {}
    with open(fai_file, 'r') as fai_h:
        for line in fai_h:
            infos = re.split('\s+', line.strip())
            infos[0] = format_chr(infos[0])
            fai_dict[infos[0]] = infos[1:]
    fai_h.close()
    # 检测染色体是否在fai内
    if chrom not in fai_dict:
        #print('%s-fai文件没有发现对应染色体！' % chrom)
        sys.exit()
    elif int(pos_end) > int(fai_dict[chrom][0]):
        pos_end = fai_dict[chrom][0]
    elif int(pos_start) <= 0:
        pos_start = '1'
    # 读取fata文件
    with open(genome_file, 'r') as genome_h:
        file_pos_start = int(int(fai_dict[chrom][1])+int(pos_start)-1+float(pos_start)//50)
        genome_h.seek(file_pos_start)
        getLeftInt = int(pos_start) // 50
        getRightInt = int(pos_end) // 50
        if  (getRightInt-getLeftInt) >= 1:
            #print("序列跨50或者100,甚至更多了啊！")
            seq_str = re.sub('\s+', '', genome_h.read(int(pos_end)-int(pos_start)+1 + (getRightInt-getLeftInt)))
        elif (getRightInt-getLeftInt) == 0:
            #print("序列没有跨50或者100！")
            seq_str = re.sub('\s+', '', genome_h.read(int(pos_end)-int(pos_start)+1))
        # print(file_pos_start); print(seq_str); sys.exit()
        # print(seq_str)
        return seq_str


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom","-chr",help="input chrom number, no need chr string.",required=True,type=str)
    parser.add_argument("--start","-s",help="start position of region.",required=True,type=str)
    parser.add_argument("--end","-e",help="end position of region.",required=True,type=str)
    parser.add_argument("--reference","-r",help="human reference.",required=False,default='NA')
    args = parser.parse_args()
    extract_fa(args.chrom,args.start,args.end,args.reference)
    #extract_fa('11', '88241166', '88241227')
    # random_obj = Random()
    # print(random_obj.get_integer(total_numbers=200, inter_max=396, inter_min=0))
