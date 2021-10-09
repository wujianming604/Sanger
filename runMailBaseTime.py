#!/usr/bin/env python3
import os
import re
import sys
import datetime
import time


if __name__ == '__main__':
    while True:
        time_now = time.strftime("%H:%M:%S", time.localtime()) # 刷新
        if time_now == "09:00:00" or time_now == "09:00:01" or time_now == "13:00:00" or time_now == "13:00:01" or time_now == "17:00:00" or time_now == "17:00:01": #此处设置每天定时的时间
            run_cmd = ' '.join(['python', 'autoMailSanger.py', '5'])
            os.system(run_cmd)
            time.sleep(2)
