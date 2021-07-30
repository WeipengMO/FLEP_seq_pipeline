#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-05-08 09:25:59
@LastEditTime : 2020-05-08 10:19:55
@Description  : 
'''


import shutil
import click
from glob import glob
import os
from math import ceil
import logging


@click.command()
@click.option('--dir_path', required=True)
@click.option('--suffix', required=True)
@click.option('--split_num', required=False, default=8)
def main(dir_path, suffix, split_num=8):
    logging.basicConfig(
        level=logging.DEBUG,  
        format='%(asctime)s %(filename)s: %(message)s',  
        datefmt='%m/%d/%Y %I:%M:%S %p',  
        filename=f'{dir_path}/split_files_into_dirs.log',  
        filemode='w'
        )
        
    files = glob(f'{dir_path}/*.{suffix}')
    files_count = len(files)
    split_size = ceil(len(files) / split_num)
    for n, i in enumerate(range(0, files_count, split_size)):
        target_file_path = os.path.join(dir_path, str(n))
        os.makedirs(target_file_path)

        with open(target_file_path+'.sm', 'w'):
            pass  # 生成空文件
        
        step = i+split_size if i+split_size < files_count else files_count
        for origin_file_path in files[i: step]:
            shutil.move(origin_file_path, target_file_path)
            logging.info(f'move {origin_file_path} to {target_file_path}')


if __name__ == "__main__":
    main()
