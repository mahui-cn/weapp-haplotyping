# -*- coding: utf-8 -*-
import sys
import json
import warnings
from os import cpu_count
from concurrent.futures import ThreadPoolExecutor, as_completed

# wegene_utils 库会包含在每个应用的环境当中，无需自行打包上传
# 这里提供源代码以供应用完整运行
import wegene_utils

from haplotyping import *

"""
当输入是部分位点时, 基因位点数据以 json 形式输入:
    {"inputs": {"RS671": "AA", "RS12203592": "CA", "format": "wegene_affy_2"}}
当输入全部位点时，全部位点对应的字符串序列会被 gzip 压缩并以 base64 转码，
转码前的数据在 data 域中
    {"inputs": {"data": "xfgakljdflkja...", format: "wegene_affy_2"}}
你需要解码并解压该数据，解压后的字符串如下:
    AACCTACCCCCC...
进一步，你需要利用相应格式的索引文件对每个位点进行解析
我们提供了整个流程的示例代码以及索引文件
"""

warnings.filterwarnings("ignore")

# 从 stdin 读取输入数据
body = sys.stdin.read()

try:
    # 如果输入的数据是全部位点数据，可以使用 wegene_utils的方法进行解析后使用
    inputs = json.loads(body)["inputs"]
    # 使用 wegene_utils 解析以后，数据会被解析成下面这样的 json 格式:
    #   {'rs123': {'genotype': 'AA', 'chromosome': '1', position: '1236'}, ...}
    user_genome = wegene_utils.process_raw_genome_data(inputs)

    # 筛选出用户的Y和mt-SNP的pos和genotype
    user_y_dict = {}
    user_mt_dict = {}
    for values in user_genome.values():
        if (
            "chromosome" in values
            and "position" in values
            and "genotype" in values
            and values["chromosome"] in {"Y", "MT"}
            and values["genotype"][0]
            in {
                "A",
                "T",
                "G",
                "C",
            }
        ):
            if values["chromosome"] == "Y":
                user_y_dict[values["position"]] = values["genotype"]
            elif values["chromosome"] == "MT":
                user_mt_dict[values["position"]] = values["genotype"]

    # 使用23Mofang单倍群树
    source = "mf"

    # 单倍群可信度阈值
    haplo_tol: float = 0.5

    # 单倍群分型对象
    yHaplo: Haplotyping = None
    mtHaplo: Haplotyping = None

    # 并发Y和mt单倍群分型，注：女性没有Y
    cpu_count = cpu_count()
    with ThreadPoolExecutor(
        max_workers=cpu_count if cpu_count != None and cpu_count > 0 else 1
    ) as t:
        task_list = []
        if len(user_y_dict) > 0:
            yHaplo = Haplotyping(
                "haplotree/{}_y_snp_tree.json".format(source.lower()), source, "y"
            )
            task_list.append(t.submit(yHaplo.analyse, user_y_dict))
        if len(user_mt_dict) > 0:
            mtHaplo = Haplotyping(
                "haplotree/{}_mt_snp_tree.json".format(source.lower()), source, "mt"
            )
            task_list.append(t.submit(mtHaplo.analyse, user_mt_dict))

        for task in as_completed(task_list):
            task.result()

    # 输出HTML
    result = []
    if (yHaplo != None and len(yHaplo.HaplogroupList) > 0) or (
        mtHaplo != None and len(mtHaplo.HaplogroupList) > 0
    ):
        if yHaplo != None and len(yHaplo.HaplogroupList) > 0:
            # 显示Y单倍群列表
            result.append(str(yHaplo))

            # 获取用户的Y家族信息
            user_y_family_dict = {}
            with open(
                "haplotree/{}_y_dict.json".format(source.lower()),
                "r",
                encoding="utf-8-sig",
            ) as y_dict_file:
                y_dict = json.load(y_dict_file)
                y_dict = y_dict if "dict" not in y_dict else y_dict["dict"]
                for y_haplo_dict in yHaplo.HaplogroupList[0]["haplo_path"]:
                    if "hf" in y_dict[y_haplo_dict["haplo"]]:
                        user_y_family_dict[y_haplo_dict["haplo"]] = y_dict[
                            y_haplo_dict["haplo"]
                        ]

        if mtHaplo != None and len(mtHaplo.HaplogroupList) > 0:
            # 显示mt单倍群列表
            result.append(str(mtHaplo))

        result.append(
            "<div class='alert alert-info' style='margin-top: 10px;' role='alert'><ul>"
        )

        if yHaplo != None and len(yHaplo.HaplogroupList) > 0:
            y_family_str = ""
            if len(user_y_family_dict) > 0:
                y_family_str += "上下游关联家族：<ul>"
                for y in user_y_family_dict.keys():
                    y_family_str += "<li>{} (共祖{}年)".format(
                        y, user_y_family_dict[y]["a"]
                    )
                    for family_dict in user_y_family_dict[y]["hf"]:
                        y_family_str += "<a href='https://www.23mofang.com/ancestry/family/{}' target='_blank' title='查看家族详情'><span class='badge margin-left-5 margin-right-5'>{}</span>{}</a>".format(
                            family_dict["fi"],
                            family_dict["sn"],
                            family_dict["ft"],
                        )
                    y_family_str += "</li>"
                y_family_str += "</ul>"
            result.append(
                "<li>您的Y父系单倍群最有可能是<span style='color: red;'>{}</span>，{}。{}</li>".format(
                    yHaplo.HaplogroupList[0]["haplo"],
                    (
                        "可信度较高"
                        if yHaplo.HaplogroupList[0]["haplo_score"] >= haplo_tol
                        else "可信度较低，可能是由于您的微基因数据有误，或不适用于此分型计算器"
                    ),
                    y_family_str,
                )
            )

        if mtHaplo != None and len(mtHaplo.HaplogroupList) > 0:
            result.append(
                "<li>您的mt母系单倍群最有可能是<span style='color: red;'>{}</span>，{}。</li>".format(
                    mtHaplo.HaplogroupList[0]["haplo"],
                    (
                        "可信度较高"
                        if mtHaplo.HaplogroupList[0]["haplo_score"] >= haplo_tol
                        else "可信度较低，可能是由于您的微基因数据有误，或不适用于此分型计算器"
                    ),
                )
            )

        result.append(
            "<li>此单倍群分型计算器基于{}的{}{}{}。采用“均衡型策略”处理SNP假阳的情况。</li>".format(
                source.upper(),
                (
                    "父系单倍群树（{:,}个单倍群）".format(yHaplo.HaploCount)
                    if yHaplo != None
                    else ""
                ),
                "，" if yHaplo != None and mtHaplo != None else "",
                (
                    "母系单倍群树（{:,}个单倍群）".format(mtHaplo.HaploCount)
                    if mtHaplo != None
                    else ""
                ),
            )
        )
        result.append(
            "<li>此次分析使用您的微基因数据中{:,}个Y-SNP和{:,}个mt-SNP，不包含nocall和indel位点。</li>".format(
                len(user_y_dict), len(user_mt_dict)
            )
        )
        result.append(
            "<li>非WeGene用户，可在<a href='http://geneu.xyz' target='_blank'>基因助手GeneU.xyz</a>上传样本使用此计算器。</li>"
        )
        result.append("</ul></div>")
    else:
        result.append(
            "<div class='alert alert-warning' role='alert>您的微基因数据没有可用于单倍群分型的信息，请更换样本重试尝试。</div>"
        )

    # 输出给用户的结果只需要通过 print 输出即可，print只可调用一次
    print("".join(result))

except Exception as e:
    # 错误信息需要被从 stderr 中输出，否则会作为正常结果输出
    for msg in e.args:
        sys.stderr.write(msg)
    exit(2)
