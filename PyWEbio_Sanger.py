#!/bin/env python
#coding=utf-8
import pywebio
import pywebio.output as output
import pywebio.input as input
import pywebio.pin as pin
from pywebio.session import hold
from pywebio import start_server
import os
import glob
import sys
os.chdir("/share_data/clin_result/clin_epilepsy_result/Sanger/PyWebIO_Sanger")


def get_ab1_list():
    ab1List = glob.glob("*ab1")
    return ab1List



def get_sanger():
    ab1List = get_ab1_list()
    print(ab1List)
    pin.put_radio("ab1",label="请选择ab1文件：",options=ab1List)
    #pin.put_input("pos",label="请输入想要截取的图片位置",type=input.NUMBER)
    pin.put_input("pos",label="请输入想要截取的图片位置")
    #pin.put_input("outName",label="请输入结果图片的名字")

    def cut():
        ab1Name = pin.pin["ab1"]
        pos = pin.pin["pos"]
        aList = ab1Name.split("-")
        outName = aList[1]+"_"+aList[2]+"-"+aList[4]+'_'+aList[5].split("_")[0]+'.png'
        os.system("/usr/bin/Rscript sangerLocalPosPlot_linux_v2_pywebio.R {} {} {}".format(ab1Name, pos, outName))
        output.put_image(open(outName,'rb').read())
        
    output.put_buttons(["分析"], lambda _:cut())
    hold()


   
if __name__ == '__main__':
    start_server(get_sanger,port=8083)
    #get_sanger()
