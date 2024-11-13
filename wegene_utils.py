# -*- coding: utf-8 -*-

__all__ = [
    "process_raw_genome_data",
    "is_genotype_exist",
    "is_wegene_format",
    "get_genome_from_tsv",
    "get_genome_from_json",
]

import sys
import gzip
import base64
from io import BytesIO
import json


def sort_genotype(genotyope):
    return "".join(sorted(genotyope))


"""
Reads the genome string anmd format and parse into a dict of
    {'rs1234': {'genotype': 'AA', 'chromosome': '1', position: '123456'}, ...}
"""


def parse_genome_string(genome_str, genome_format):
    try:
        genome_dict = {}
        # Index files for all posible formats will be provided automatically
        # Do not change the default path below if you wish to use those
        with open("./indexes/index_" + genome_format + ".idx", "r") as idx_f:
            # 统计基因数据模板的行数，
            tmpl_line_count = 0
            for line in idx_f:
                if not line.startswith("NA") and line.strip() != "":
                    tmpl_line_count += 1
                    fields = line.strip().split("\t")
                    index_pos = fields[0].strip()
                    rsid = fields[1].strip()
                    chromosome = fields[2].strip()
                    position = fields[3].strip()
                    start_pos = int(index_pos) * 2
                    genome_dict[rsid] = {
                        "genotype": sort_genotype(
                            genome_str[start_pos : start_pos + 2]
                        ),
                        "chromosome": chromosome,
                        "position": position,
                    }

            # 校验用户分型的位点数和基因索引是否一致
            if len(genome_str) // 2 != tmpl_line_count:
                raise Exception(
                    "您的基因样本位点数：{:,}，基因模板（{}）位点数：{:,}，两者不一致无法解析！".format(
                        len(genome_str) // 2, genome_format, tmpl_line_count
                    )
                )

        return genome_dict
    except Exception as e:
        raise e


def process_raw_genome_data(raw_inputs):
    try:
        genome = str(
            gzip.GzipFile(fileobj=BytesIO(base64.b64decode(raw_inputs["data"]))).read(),
            encoding="utf8",
        )
        genome_format = raw_inputs["format"]
        return parse_genome_string(genome, genome_format)
    except Exception as e:
        raise e


def is_genotype_exist(input, rsid):
    return rsid in input and input[rsid] != "--" and input[rsid] != "__"


def is_wegene_format(format_str):
    return "wegene_" in format_str


# 从TSV格式文件加载基因数据
def get_genome_from_tsv(tsvFileName):
    user_genome = {}
    with open(tsvFileName, "r", encoding="utf-8") as tsvFile:
        for tsvLine in tsvFile.readlines():
            if len(tsvLine) > 0 and not tsvLine.startswith(("#", "\n", "\t", '"')):
                tsvLineArray = tsvLine.split("\t")
                if len(tsvLineArray) == 4 and tsvLineArray[3][0] in [
                    "A",
                    "T",
                    "G",
                    "C",
                ]:
                    user_genome[tsvLineArray[0]] = {
                        "chromosome": tsvLineArray[1].strip(),
                        "position": tsvLineArray[2].strip(),
                        "genotype": tsvLineArray[3].strip(),
                    }

    return user_genome


# 从JSON格式文件加载基因数据
def get_genome_from_json(jsonFileName):
    user_genome = {}
    with open(jsonFileName, "r", encoding="utf-8") as jsonFile:
        jFile = jsonFile.read()
        if len(jFile) > 0:
            user_genome = json.loads(jFile)

    return user_genome
